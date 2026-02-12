#include <metal_stdlib>
using namespace metal;

// Tiled matrix multiply: C[m,n] = A[m,k] * B[k,n], row-major.
// 16x16 tiles in threadgroup shared memory.

constant uint TILE_SIZE = 16;

kernel void matmul(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float* C [[buffer(2)]],
    device const uint& M [[buffer(3)]],
    device const uint& K [[buffer(4)]],
    device const uint& N [[buffer(5)]],
    uint2 gid [[thread_position_in_grid]],
    uint2 tid [[thread_position_in_threadgroup]]
) {
    uint row = gid.y;
    uint col = gid.x;

    threadgroup float tileA[TILE_SIZE][TILE_SIZE];
    threadgroup float tileB[TILE_SIZE][TILE_SIZE];

    float sum = 0.0f;
    uint num_tiles = (K + TILE_SIZE - 1) / TILE_SIZE;

    for (uint t = 0; t < num_tiles; t++) {
        // Load tile from A
        uint a_col = t * TILE_SIZE + tid.x;
        if (row < M && a_col < K) {
            tileA[tid.y][tid.x] = A[row * K + a_col];
        } else {
            tileA[tid.y][tid.x] = 0.0f;
        }

        // Load tile from B
        uint b_row = t * TILE_SIZE + tid.y;
        if (b_row < K && col < N) {
            tileB[tid.y][tid.x] = B[b_row * N + col];
        } else {
            tileB[tid.y][tid.x] = 0.0f;
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        // Compute partial product
        for (uint p = 0; p < TILE_SIZE; p++) {
            sum += tileA[tid.y][p] * tileB[p][tid.x];
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (row < M && col < N) {
        C[row * N + col] = sum;
    }
}
