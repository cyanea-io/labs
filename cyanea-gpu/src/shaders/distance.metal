#include <metal_stdlib>
using namespace metal;

// Pairwise distance matrix: 2D grid (n, n), each thread computes one cell.
// Upper-triangle optimization: compute once, mirror to lower triangle.

kernel void pairwise_euclidean(
    device const float* data [[buffer(0)]],
    device float* result [[buffer(1)]],
    device const uint& n [[buffer(2)]],
    device const uint& dim [[buffer(3)]],
    uint2 gid [[thread_position_in_grid]]
) {
    uint i = gid.y;
    uint j = gid.x;
    if (i >= n || j >= n) return;

    if (i == j) {
        result[i * n + j] = 0.0f;
        return;
    }
    // Only compute upper triangle
    if (i > j) return;

    float sum_sq = 0.0f;
    for (uint d = 0; d < dim; d++) {
        float diff = data[i * dim + d] - data[j * dim + d];
        sum_sq += diff * diff;
    }
    float dist = sqrt(sum_sq);
    result[i * n + j] = dist;
    result[j * n + i] = dist;
}

kernel void pairwise_manhattan(
    device const float* data [[buffer(0)]],
    device float* result [[buffer(1)]],
    device const uint& n [[buffer(2)]],
    device const uint& dim [[buffer(3)]],
    uint2 gid [[thread_position_in_grid]]
) {
    uint i = gid.y;
    uint j = gid.x;
    if (i >= n || j >= n) return;

    if (i == j) {
        result[i * n + j] = 0.0f;
        return;
    }
    if (i > j) return;

    float sum = 0.0f;
    for (uint d = 0; d < dim; d++) {
        sum += abs(data[i * dim + d] - data[j * dim + d]);
    }
    result[i * n + j] = sum;
    result[j * n + i] = sum;
}

kernel void pairwise_cosine(
    device const float* data [[buffer(0)]],
    device float* result [[buffer(1)]],
    device const uint& n [[buffer(2)]],
    device const uint& dim [[buffer(3)]],
    uint2 gid [[thread_position_in_grid]]
) {
    uint i = gid.y;
    uint j = gid.x;
    if (i >= n || j >= n) return;

    if (i == j) {
        result[i * n + j] = 0.0f;
        return;
    }
    if (i > j) return;

    float dot_val = 0.0f;
    float norm_a = 0.0f;
    float norm_b = 0.0f;
    for (uint d = 0; d < dim; d++) {
        float a = data[i * dim + d];
        float b = data[j * dim + d];
        dot_val += a * b;
        norm_a += a * a;
        norm_b += b * b;
    }
    float denom = sqrt(norm_a) * sqrt(norm_b);
    float dist = (denom == 0.0f) ? 0.0f : (1.0f - dot_val / denom);
    result[i * n + j] = dist;
    result[j * n + i] = dist;
}
