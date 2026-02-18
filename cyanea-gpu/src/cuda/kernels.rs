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

/// CUDA kernel source for GPU k-mer counting.
pub const KMER_KERNEL_SOURCE: &str = r#"
extern "C" {

__global__ void kmer_count(
    const unsigned char* seq,
    unsigned int seq_len,
    unsigned int k,
    unsigned int* counts
) {
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid + k > seq_len) return;

    unsigned int hash = 0;
    for (unsigned int i = 0; i < k; i++) {
        unsigned char base = seq[gid + i];
        unsigned int code;
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else return; // skip k-mers with N
        hash = (hash << 2) | code;
    }

    atomicAdd(&counts[hash], 1u);
}

} // extern "C"
"#;

/// CUDA kernel source for GPU Smith-Waterman protein alignment.
pub const SW_KERNEL_SOURCE: &str = r#"
extern "C" {

// Amino acid to index lookup table
__constant__ int AA_MAP[128] = {
    22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,
    22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,
    22,22,22,22,22,22,22,22,22,22,23,22,22,22,22,22,
    22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,
    22, 0,20, 4, 3, 6,13, 7, 8, 9,22,11,10,12, 2,22,
    14, 5, 1,15,16,22,19,17,22,18,21,22,22,22,22,22,
    22, 0,20, 4, 3, 6,13, 7, 8, 9,22,11,10,12, 2,22,
    14, 5, 1,15,16,22,19,17,22,18,21,22,22,22,22,22
};

__global__ void smith_waterman_batch(
    const unsigned char* queries,
    const unsigned char* targets,
    const unsigned int* query_lengths,
    const unsigned int* target_lengths,
    const unsigned int* query_offsets,
    const unsigned int* target_offsets,
    const int* submat,
    int gap_open, int gap_extend,
    int* scores,
    unsigned int* query_ends,
    unsigned int* target_ends,
    unsigned int n_pairs
) {
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n_pairs) return;

    unsigned int qlen = query_lengths[gid];
    unsigned int tlen = target_lengths[gid];
    unsigned int qoff = query_offsets[gid];
    unsigned int toff = target_offsets[gid];

    int best = 0;
    unsigned int best_qi = 0, best_ti = 0;

    // Allocate DP arrays on stack (limited to ~512 residues)
    int prev_h[513];
    int curr_e[513];

    for (unsigned int j = 0; j <= tlen; j++) {
        prev_h[j] = 0;
        curr_e[j] = 0;
    }

    for (unsigned int i = 1; i <= qlen; i++) {
        int qi = AA_MAP[queries[qoff + i - 1] & 0x7F];
        int diag = 0;
        int f = -1000000;

        for (unsigned int j = 1; j <= tlen; j++) {
            int ti = AA_MAP[targets[toff + j - 1] & 0x7F];
            int sc = submat[qi * 24 + ti];

            int eo = prev_h[j] + gap_open + gap_extend;
            int ee = curr_e[j] + gap_extend;
            curr_e[j] = (eo > ee) ? eo : ee;

            int fo = diag + gap_open + gap_extend;
            int fe = f + gap_extend;
            f = (fo > fe) ? fo : fe;

            int old_h = prev_h[j];
            int h = diag + sc;
            if (curr_e[j] > h) h = curr_e[j];
            if (f > h) h = f;
            if (h < 0) h = 0;

            diag = old_h;
            prev_h[j] = h;

            if (h > best) {
                best = h;
                best_qi = i - 1;
                best_ti = j - 1;
            }
        }
    }

    scores[gid] = best;
    query_ends[gid] = best_qi;
    target_ends[gid] = best_ti;
}

} // extern "C"
"#;

/// CUDA kernel source for GPU MinHash k-mer hashing.
pub const MINHASH_KERNEL_SOURCE: &str = r#"
extern "C" {

__device__ unsigned long long murmur3_mix(unsigned long long h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

__global__ void minhash_hash_kmers(
    const unsigned char* seq,
    unsigned int seq_len,
    unsigned int k,
    unsigned long long* hashes
) {
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid + k > seq_len) return;

    unsigned long long kmer_val = 0;
    for (unsigned int i = 0; i < k; i++) {
        unsigned char base = seq[gid + i];
        unsigned int code;
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else {
            hashes[gid] = 0xFFFFFFFFFFFFFFFFULL;
            return;
        }
        kmer_val = (kmer_val << 2) | (unsigned long long)code;
    }

    hashes[gid] = murmur3_mix(kmer_val ^ (unsigned long long)k);
}

} // extern "C"
"#;
