// Tiled matrix multiply for WebGPU backend.
// C[m,n] = A[m,k] * B[k,n], row-major.
// 16×16 tiles in workgroup shared memory.

@group(0) @binding(0) var<storage, read> A: array<f32>;
@group(0) @binding(1) var<storage, read> B: array<f32>;
@group(0) @binding(2) var<storage, read_write> C: array<f32>;
@group(0) @binding(3) var<uniform> params: MatmulParams;

struct MatmulParams {
    M: u32,
    K: u32,
    N: u32,
}

const TILE: u32 = 16;

var<workgroup> tileA: array<array<f32, 16>, 16>;
var<workgroup> tileB: array<array<f32, 16>, 16>;

@compute @workgroup_size(16, 16)
fn matmul(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) tid: vec3<u32>,
) {
    let row = gid.y;
    let col = gid.x;
    let ty = tid.y;
    let tx = tid.x;

    var sum: f32 = 0.0;
    let num_tiles = (params.K + TILE - 1) / TILE;

    for (var t = 0u; t < num_tiles; t++) {
        let a_col = t * TILE + tx;
        if (row < params.M && a_col < params.K) {
            tileA[ty][tx] = A[row * params.K + a_col];
        } else {
            tileA[ty][tx] = 0.0;
        }

        let b_row = t * TILE + ty;
        if (b_row < params.K && col < params.N) {
            tileB[ty][tx] = B[b_row * params.N + col];
        } else {
            tileB[ty][tx] = 0.0;
        }

        workgroupBarrier();

        for (var p = 0u; p < TILE; p++) {
            sum += tileA[ty][p] * tileB[p][tx];
        }

        workgroupBarrier();
    }

    if (row < params.M && col < params.N) {
        C[row * params.N + col] = sum;
    }
}
