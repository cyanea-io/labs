// Pairwise distance matrix shaders for WebGPU backend.
// 2D dispatch (n, n), each thread computes one cell.
// Upper-triangle optimization: compute once, mirror to lower triangle.

@group(0) @binding(0) var<storage, read> data: array<f32>;
@group(0) @binding(1) var<storage, read_write> result: array<f32>;
@group(0) @binding(2) var<uniform> params: DistanceParams;

struct DistanceParams {
    n: u32,
    dim: u32,
}

@compute @workgroup_size(16, 16)
fn pairwise_euclidean(
    @builtin(global_invocation_id) gid: vec3<u32>,
) {
    let i = gid.y;
    let j = gid.x;
    if (i >= params.n || j >= params.n) { return; }

    if (i == j) {
        result[i * params.n + j] = 0.0;
        return;
    }
    if (i > j) { return; }

    var sum_sq: f32 = 0.0;
    for (var d = 0u; d < params.dim; d++) {
        let diff = data[i * params.dim + d] - data[j * params.dim + d];
        sum_sq += diff * diff;
    }
    let dist = sqrt(sum_sq);
    result[i * params.n + j] = dist;
    result[j * params.n + i] = dist;
}

@compute @workgroup_size(16, 16)
fn pairwise_manhattan(
    @builtin(global_invocation_id) gid: vec3<u32>,
) {
    let i = gid.y;
    let j = gid.x;
    if (i >= params.n || j >= params.n) { return; }

    if (i == j) {
        result[i * params.n + j] = 0.0;
        return;
    }
    if (i > j) { return; }

    var sum: f32 = 0.0;
    for (var d = 0u; d < params.dim; d++) {
        sum += abs(data[i * params.dim + d] - data[j * params.dim + d]);
    }
    result[i * params.n + j] = sum;
    result[j * params.n + i] = sum;
}

@compute @workgroup_size(16, 16)
fn pairwise_cosine(
    @builtin(global_invocation_id) gid: vec3<u32>,
) {
    let i = gid.y;
    let j = gid.x;
    if (i >= params.n || j >= params.n) { return; }

    if (i == j) {
        result[i * params.n + j] = 0.0;
        return;
    }
    if (i > j) { return; }

    var dot_val: f32 = 0.0;
    var norm_a: f32 = 0.0;
    var norm_b: f32 = 0.0;
    for (var d = 0u; d < params.dim; d++) {
        let a = data[i * params.dim + d];
        let b = data[j * params.dim + d];
        dot_val += a * b;
        norm_a += a * a;
        norm_b += b * b;
    }
    let denom = sqrt(norm_a) * sqrt(norm_b);
    var dist: f32;
    if (denom == 0.0) {
        dist = 0.0;
    } else {
        dist = 1.0 - dot_val / denom;
    }
    result[i * params.n + j] = dist;
    result[j * params.n + i] = dist;
}
