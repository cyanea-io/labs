#include <metal_stdlib>
using namespace metal;

// GPU k-mer counting: each thread processes one starting position in the sequence.
// Encodes k bases into a 2-bit hash and atomically increments the count table.
// Base encoding: A=0, C=1, G=2, T=3. Positions with N are skipped.

kernel void kmer_count(
    device const uchar* seq [[buffer(0)]],
    device const uint& seq_len [[buffer(1)]],
    device const uint& k [[buffer(2)]],
    device atomic_uint* counts [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid + k > seq_len) return;

    uint hash = 0;
    for (uint i = 0; i < k; i++) {
        uchar base = seq[gid + i];
        uint code;
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else return; // skip k-mers containing N or other ambiguous bases
        hash = (hash << 2) | code;
    }

    atomic_fetch_add_explicit(&counts[hash], 1u, memory_order_relaxed);
}
