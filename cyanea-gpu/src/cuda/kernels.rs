//! CUDA C kernel source for runtime compilation via NVRTC.

/// Combined CUDA kernel source for all GPU operations.
pub const KERNEL_SOURCE: &str = r#"
extern "C" {

// ── Reductions ─────────────────────────────────────────────────────
// 256-thread blocks, tree reduction in shared memory, two-pass.

__global__ void reduce_sum(
    const double* input,
    double* output,
    unsigned int count
) {
    __shared__ double sdata[256];
    unsigned int tid = threadIdx.x;
    unsigned int gid = blockIdx.x * blockDim.x + tid;

    sdata[tid] = (gid < count) ? input[gid] : 0.0;
    __syncthreads();

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) output[blockIdx.x] = sdata[0];
}

__global__ void reduce_min(
    const double* input,
    double* output,
    unsigned int count
) {
    __shared__ double sdata[256];
    unsigned int tid = threadIdx.x;
    unsigned int gid = blockIdx.x * blockDim.x + tid;

    sdata[tid] = (gid < count) ? input[gid] : __longlong_as_double(0x7FF0000000000000LL);
    __syncthreads();

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (sdata[tid + s] < sdata[tid]) sdata[tid] = sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) output[blockIdx.x] = sdata[0];
}

__global__ void reduce_max(
    const double* input,
    double* output,
    unsigned int count
) {
    __shared__ double sdata[256];
    unsigned int tid = threadIdx.x;
    unsigned int gid = blockIdx.x * blockDim.x + tid;

    sdata[tid] = (gid < count) ? input[gid] : __longlong_as_double(0xFFF0000000000000LL);
    __syncthreads();

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (sdata[tid + s] > sdata[tid]) sdata[tid] = sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) output[blockIdx.x] = sdata[0];
}

// ── Pairwise distance ──────────────────────────────────────────────
// 2D grid (n, n), upper triangle + mirror.

__global__ void pairwise_euclidean(
    const double* data,
    double* result,
    unsigned int n,
    unsigned int dim
) {
    unsigned int i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n || j >= n) return;

    if (i == j) { result[i * n + j] = 0.0; return; }
    if (i > j) return;

    double sum_sq = 0.0;
    for (unsigned int d = 0; d < dim; d++) {
        double diff = data[i * dim + d] - data[j * dim + d];
        sum_sq += diff * diff;
    }
    double dist = sqrt(sum_sq);
    result[i * n + j] = dist;
    result[j * n + i] = dist;
}

__global__ void pairwise_manhattan(
    const double* data,
    double* result,
    unsigned int n,
    unsigned int dim
) {
    unsigned int i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n || j >= n) return;

    if (i == j) { result[i * n + j] = 0.0; return; }
    if (i > j) return;

    double sum = 0.0;
    for (unsigned int d = 0; d < dim; d++) {
        double diff = data[i * dim + d] - data[j * dim + d];
        sum += (diff < 0.0) ? -diff : diff;
    }
    result[i * n + j] = sum;
    result[j * n + i] = sum;
}

__global__ void pairwise_cosine(
    const double* data,
    double* result,
    unsigned int n,
    unsigned int dim
) {
    unsigned int i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n || j >= n) return;

    if (i == j) { result[i * n + j] = 0.0; return; }
    if (i > j) return;

    double dot_val = 0.0, norm_a = 0.0, norm_b = 0.0;
    for (unsigned int d = 0; d < dim; d++) {
        double a = data[i * dim + d];
        double b = data[j * dim + d];
        dot_val += a * b;
        norm_a += a * a;
        norm_b += b * b;
    }
    double denom = sqrt(norm_a) * sqrt(norm_b);
    double dist = (denom == 0.0) ? 0.0 : (1.0 - dot_val / denom);
    result[i * n + j] = dist;
    result[j * n + i] = dist;
}

// ── Matrix multiply ────────────────────────────────────────────────
// Tiled with 16x16 shared memory. C[m,n] = A[m,k] * B[k,n], row-major.

#define TILE_SIZE 16

__global__ void matmul(
    const double* A,
    const double* B,
    double* C,
    unsigned int M,
    unsigned int K,
    unsigned int N
) {
    __shared__ double tileA[TILE_SIZE][TILE_SIZE];
    __shared__ double tileB[TILE_SIZE][TILE_SIZE];

    unsigned int row = blockIdx.y * TILE_SIZE + threadIdx.y;
    unsigned int col = blockIdx.x * TILE_SIZE + threadIdx.x;
    unsigned int ty = threadIdx.y;
    unsigned int tx = threadIdx.x;

    double sum = 0.0;
    unsigned int num_tiles = (K + TILE_SIZE - 1) / TILE_SIZE;

    for (unsigned int t = 0; t < num_tiles; t++) {
        unsigned int a_col = t * TILE_SIZE + tx;
        tileA[ty][tx] = (row < M && a_col < K) ? A[row * K + a_col] : 0.0;

        unsigned int b_row = t * TILE_SIZE + ty;
        tileB[ty][tx] = (b_row < K && col < N) ? B[b_row * N + col] : 0.0;

        __syncthreads();

        for (unsigned int p = 0; p < TILE_SIZE; p++) {
            sum += tileA[ty][p] * tileB[p][tx];
        }

        __syncthreads();
    }

    if (row < M && col < N) {
        C[row * N + col] = sum;
    }
}

} // extern "C"
"#;
