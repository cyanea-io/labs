#include <metal_stdlib>
using namespace metal;

// GPU MinHash: parallel k-mer hashing with MurmurHash3 (64-bit).
// Each thread hashes one k-mer position and writes the hash value.
// Host-side selects the bottom-k smallest hashes.

// MurmurHash3 64-bit finalizer
inline ulong murmur3_mix(ulong h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdUL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53UL;
    h ^= h >> 33;
    return h;
}

kernel void minhash_hash_kmers(
    device const uchar* seq [[buffer(0)]],
    device const uint& seq_len [[buffer(1)]],
    device const uint& k [[buffer(2)]],
    device ulong* hashes [[buffer(3)]],
    device atomic_uint* valid_count [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid + k > seq_len) return;

    // Encode k-mer as a 64-bit value using 2-bit encoding
    // For k > 32, we just hash the bytes directly
    ulong kmer_hash = 0;
    for (uint i = 0; i < k; i++) {
        uchar base = seq[gid + i];
        uint code;
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else {
            hashes[gid] = 0xFFFFFFFFFFFFFFFFUL; // sentinel for invalid
            return;
        }
        kmer_hash = (kmer_hash << 2) | ulong(code);
    }

    // Apply MurmurHash3 finalizer
    kmer_hash = murmur3_mix(kmer_hash ^ ulong(k));
    hashes[gid] = kmer_hash;
    atomic_fetch_add_explicit(valid_count, 1u, memory_order_relaxed);
}
