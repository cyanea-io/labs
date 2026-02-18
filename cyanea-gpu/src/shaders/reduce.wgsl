// Reduction shaders for WebGPU backend (sum, min, max).
// 256 threads per workgroup, tree reduction in workgroup shared memory.
// Two-phase: GPU reduces within workgroups, host sums partial results.

@group(0) @binding(0) var<storage, read> input: array<f32>;
@group(0) @binding(1) var<storage, read_write> output: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;

struct Params {
    count: u32,
}

var<workgroup> shared_data: array<f32, 256>;

@compute @workgroup_size(256)
fn reduce_sum(
    @builtin(local_invocation_id) tid: vec3<u32>,
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(workgroup_id) wg_id: vec3<u32>,
) {
    let local_id = tid.x;
    let global_id = gid.x;

    if (global_id < params.count) {
        shared_data[local_id] = input[global_id];
    } else {
        shared_data[local_id] = 0.0;
    }
    workgroupBarrier();

    for (var stride = 128u; stride > 0u; stride >>= 1u) {
        if (local_id < stride) {
            shared_data[local_id] += shared_data[local_id + stride];
        }
        workgroupBarrier();
    }

    if (local_id == 0u) {
        output[wg_id.x] = shared_data[0];
    }
}

@compute @workgroup_size(256)
fn reduce_min(
    @builtin(local_invocation_id) tid: vec3<u32>,
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(workgroup_id) wg_id: vec3<u32>,
) {
    let local_id = tid.x;
    let global_id = gid.x;

    if (global_id < params.count) {
        shared_data[local_id] = input[global_id];
    } else {
        shared_data[local_id] = 3.402823466e+38; // f32 max as infinity substitute
    }
    workgroupBarrier();

    for (var stride = 128u; stride > 0u; stride >>= 1u) {
        if (local_id < stride) {
            shared_data[local_id] = min(shared_data[local_id], shared_data[local_id + stride]);
        }
        workgroupBarrier();
    }

    if (local_id == 0u) {
        output[wg_id.x] = shared_data[0];
    }
}

@compute @workgroup_size(256)
fn reduce_max(
    @builtin(local_invocation_id) tid: vec3<u32>,
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(workgroup_id) wg_id: vec3<u32>,
) {
    let local_id = tid.x;
    let global_id = gid.x;

    if (global_id < params.count) {
        shared_data[local_id] = input[global_id];
    } else {
        shared_data[local_id] = -3.402823466e+38; // f32 min as -infinity substitute
    }
    workgroupBarrier();

    for (var stride = 128u; stride > 0u; stride >>= 1u) {
        if (local_id < stride) {
            shared_data[local_id] = max(shared_data[local_id], shared_data[local_id + stride]);
        }
        workgroupBarrier();
    }

    if (local_id == 0u) {
        output[wg_id.x] = shared_data[0];
    }
}
