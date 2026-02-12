#include <metal_stdlib>
using namespace metal;

// Threadgroup-shared memory for tree reduction, 256 threads per group.
// Two-phase: GPU reduces within threadgroups, host sums partial results.

kernel void reduce_sum(
    device const float* input [[buffer(0)]],
    device float* output [[buffer(1)]],
    device const uint& count [[buffer(2)]],
    uint tid [[thread_position_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint group_id [[threadgroup_position_in_grid]],
    uint group_size [[threads_per_threadgroup]]
) {
    threadgroup float shared_data[256];
    shared_data[tid] = (gid < count) ? input[gid] : 0.0f;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = group_size / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            shared_data[tid] += shared_data[tid + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (tid == 0) {
        output[group_id] = shared_data[0];
    }
}

kernel void reduce_min(
    device const float* input [[buffer(0)]],
    device float* output [[buffer(1)]],
    device const uint& count [[buffer(2)]],
    uint tid [[thread_position_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint group_id [[threadgroup_position_in_grid]],
    uint group_size [[threads_per_threadgroup]]
) {
    threadgroup float shared_data[256];
    shared_data[tid] = (gid < count) ? input[gid] : INFINITY;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = group_size / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            shared_data[tid] = min(shared_data[tid], shared_data[tid + stride]);
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (tid == 0) {
        output[group_id] = shared_data[0];
    }
}

kernel void reduce_max(
    device const float* input [[buffer(0)]],
    device float* output [[buffer(1)]],
    device const uint& count [[buffer(2)]],
    uint tid [[thread_position_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint group_id [[threadgroup_position_in_grid]],
    uint group_size [[threads_per_threadgroup]]
) {
    threadgroup float shared_data[256];
    shared_data[tid] = (gid < count) ? input[gid] : -INFINITY;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = group_size / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            shared_data[tid] = max(shared_data[tid], shared_data[tid + stride]);
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (tid == 0) {
        output[group_id] = shared_data[0];
    }
}
