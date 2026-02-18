//! GPU-accelerated MinHash sketch computation.
//!
//! Provides parallel MinHash sketching for DNA sequences using GPU
//! (Metal/CUDA) with automatic CPU fallback. K-mers are hashed with
//! MurmurHash3 (64-bit) and the `sketch_size` smallest hash values
//! are retained (bottom-sketch).

use cyanea_core::{CyaneaError, Result};

/// A MinHash sketch for a single sequence.
#[derive(Debug, Clone)]
pub struct MinHashSketch {
    /// Sorted smallest hash values (bottom-sketch).
    pub hashes: Vec<u64>,
    /// The k-mer length used.
    pub k: usize,
}

/// Computes MinHash sketches for the given DNA sequences using the best
/// available backend.
///
/// Each sequence produces one `MinHashSketch` containing the `sketch_size`
/// smallest 64-bit hash values of its k-mers.
///
/// # Errors
///
/// Returns an error if `k == 0`, `sketch_size == 0`, or `sequences` is empty.
pub fn gpu_minhash(
    sequences: &[&[u8]],
    k: usize,
    sketch_size: usize,
) -> Result<Vec<MinHashSketch>> {
    validate_inputs(sequences, k, sketch_size)?;

    #[cfg(feature = "metal")]
    {
        if let Ok(result) = gpu_minhash_metal(sequences, k, sketch_size) {
            return Ok(result);
        }
    }

    gpu_minhash_cpu(sequences, k, sketch_size)
}

/// CPU-only MinHash computation.
pub fn gpu_minhash_cpu(
    sequences: &[&[u8]],
    k: usize,
    sketch_size: usize,
) -> Result<Vec<MinHashSketch>> {
    validate_inputs(sequences, k, sketch_size)?;

    let mut sketches = Vec::with_capacity(sequences.len());
    for seq in sequences {
        let sketch = compute_sketch_cpu(seq, k, sketch_size);
        sketches.push(sketch);
    }
    Ok(sketches)
}

/// Estimates Jaccard similarity from two MinHash sketches.
///
/// Both sketches must have the same k-mer length. Returns the fraction
/// of shared hash values.
///
/// # Errors
///
/// Returns an error if sketches have different k values or are empty.
pub fn gpu_minhash_jaccard(a: &MinHashSketch, b: &MinHashSketch) -> Result<f64> {
    if a.k != b.k {
        return Err(CyaneaError::InvalidInput(format!(
            "minhash_jaccard: k mismatch ({} vs {})",
            a.k, b.k
        )));
    }
    if a.hashes.is_empty() || b.hashes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "minhash_jaccard: sketches must not be empty".to_string(),
        ));
    }

    // Merge-based intersection count on sorted sketches
    let mut i = 0;
    let mut j = 0;
    let mut shared = 0usize;
    let mut total = 0usize;

    while i < a.hashes.len() && j < b.hashes.len() {
        if a.hashes[i] == b.hashes[j] {
            shared += 1;
            total += 1;
            i += 1;
            j += 1;
        } else if a.hashes[i] < b.hashes[j] {
            total += 1;
            i += 1;
        } else {
            total += 1;
            j += 1;
        }
    }
    total += (a.hashes.len() - i) + (b.hashes.len() - j);

    if total == 0 {
        return Ok(0.0);
    }
    Ok(shared as f64 / total as f64)
}

/// MurmurHash3 64-bit finalizer.
fn murmur3_mix(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}

/// Encodes a DNA base to 2-bit. Returns `None` for non-ACGT.
fn base_to_2bit(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn compute_sketch_cpu(seq: &[u8], k: usize, sketch_size: usize) -> MinHashSketch {
    use std::collections::BinaryHeap;

    if seq.len() < k {
        return MinHashSketch {
            hashes: vec![],
            k,
        };
    }

    // Use a max-heap of size sketch_size to keep track of the bottom-k hashes
    let mut heap = BinaryHeap::with_capacity(sketch_size + 1);

    for pos in 0..=(seq.len() - k) {
        let kmer = &seq[pos..pos + k];
        // Encode k-mer
        let mut kmer_val = 0u64;
        let mut valid = true;
        for &base in kmer {
            if let Some(code) = base_to_2bit(base) {
                kmer_val = (kmer_val << 2) | code;
            } else {
                valid = false;
                break;
            }
        }
        if !valid {
            continue;
        }

        let hash = murmur3_mix(kmer_val ^ (k as u64));

        if heap.len() < sketch_size {
            heap.push(hash);
        } else if let Some(&max) = heap.peek() {
            if hash < max {
                heap.pop();
                heap.push(hash);
            }
        }
    }

    let mut hashes: Vec<u64> = heap.into_vec();
    hashes.sort_unstable();

    MinHashSketch { hashes, k }
}

fn validate_inputs(sequences: &[&[u8]], k: usize, sketch_size: usize) -> Result<()> {
    if k == 0 {
        return Err(CyaneaError::InvalidInput(
            "minhash: k must be > 0".to_string(),
        ));
    }
    if sketch_size == 0 {
        return Err(CyaneaError::InvalidInput(
            "minhash: sketch_size must be > 0".to_string(),
        ));
    }
    if sequences.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "minhash: sequences must not be empty".to_string(),
        ));
    }
    Ok(())
}

// ── Metal GPU dispatch ────────────────────────────────────────────

#[cfg(feature = "metal")]
fn gpu_minhash_metal(
    sequences: &[&[u8]],
    k: usize,
    sketch_size: usize,
) -> Result<Vec<MinHashSketch>> {
    use metal_rs::{CompileOptions, Device, MTLResourceOptions, MTLSize};
    use std::ffi::c_void;

    let device = Device::system_default()
        .ok_or_else(|| CyaneaError::Other("no Metal device available".into()))?;
    let queue = device.new_command_queue();

    let shader_src = include_str!("shaders/minhash.metal");
    let opts = CompileOptions::new();
    let lib = device
        .new_library_with_source(shader_src, &opts)
        .map_err(|e| CyaneaError::Other(format!("Metal minhash shader compile: {e}")))?;
    let func = lib
        .get_function("minhash_hash_kmers", None)
        .map_err(|e| CyaneaError::Other(format!("Metal minhash function: {e}")))?;
    let pipeline = device
        .new_compute_pipeline_state_with_function(&func)
        .map_err(|e| CyaneaError::Other(format!("Metal minhash pipeline: {e}")))?;

    let k_u32 = k as u32;
    let k_buf = device.new_buffer_with_data(
        &k_u32 as *const u32 as *const c_void,
        4,
        MTLResourceOptions::StorageModeShared,
    );

    let mut sketches = Vec::with_capacity(sequences.len());

    for seq in sequences {
        if seq.len() < k {
            sketches.push(MinHashSketch {
                hashes: vec![],
                k,
            });
            continue;
        }

        let n_kmers = seq.len() - k + 1;
        let seq_len = seq.len() as u32;

        let seq_buf = device.new_buffer_with_data(
            seq.as_ptr() as *const c_void,
            seq.len() as u64,
            MTLResourceOptions::StorageModeShared,
        );
        let len_buf = device.new_buffer_with_data(
            &seq_len as *const u32 as *const c_void,
            4,
            MTLResourceOptions::StorageModeShared,
        );
        let hashes_buf = device.new_buffer(
            (n_kmers * 8) as u64,
            MTLResourceOptions::StorageModeShared,
        );
        let valid_buf = device.new_buffer(4, MTLResourceOptions::StorageModeShared);
        unsafe {
            std::ptr::write_bytes(valid_buf.contents() as *mut u8, 0, 4);
        }

        let cmd = queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&pipeline);
        enc.set_buffer(0, Some(&seq_buf), 0);
        enc.set_buffer(1, Some(&len_buf), 0);
        enc.set_buffer(2, Some(&k_buf), 0);
        enc.set_buffer(3, Some(&hashes_buf), 0);
        enc.set_buffer(4, Some(&valid_buf), 0);

        let grid = MTLSize::new(n_kmers as u64, 1, 1);
        let tg = MTLSize::new(256.min(n_kmers as u64).max(1), 1, 1);
        enc.dispatch_threads(grid, tg);
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();

        // Read hashes and select bottom-k
        let ptr = hashes_buf.contents() as *const u64;
        let all_hashes = unsafe { std::slice::from_raw_parts(ptr, n_kmers) };

        // Filter out sentinel values and select bottom-k
        let mut valid_hashes: Vec<u64> = all_hashes
            .iter()
            .copied()
            .filter(|&h| h != u64::MAX)
            .collect();
        valid_hashes.sort_unstable();
        valid_hashes.truncate(sketch_size);

        sketches.push(MinHashSketch {
            hashes: valid_hashes,
            k,
        });
    }

    Ok(sketches)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minhash_deterministic() {
        let seq = b"ACGTACGTACGTACGT";
        let s1 = gpu_minhash_cpu(&[seq.as_slice()], 4, 10).unwrap();
        let s2 = gpu_minhash_cpu(&[seq.as_slice()], 4, 10).unwrap();
        assert_eq!(s1[0].hashes, s2[0].hashes);
    }

    #[test]
    fn minhash_identical_jaccard_one() {
        let seq = b"ACGTACGTACGTACGT";
        let sketches = gpu_minhash_cpu(&[seq.as_slice(), seq.as_slice()], 4, 100).unwrap();
        let j = gpu_minhash_jaccard(&sketches[0], &sketches[1]).unwrap();
        assert!((j - 1.0).abs() < 1e-10, "Jaccard should be 1.0, got {j}");
    }

    #[test]
    fn minhash_disjoint_jaccard_zero() {
        // Two sequences with completely different k-mer sets
        let seq1 = b"AAAAAAAAAA"; // only AA, AAA, AAAA k-mers
        let seq2 = b"CCCCCCCCCC"; // only CC, CCC, CCCC k-mers
        let sketches =
            gpu_minhash_cpu(&[seq1.as_slice(), seq2.as_slice()], 4, 100).unwrap();
        let j = gpu_minhash_jaccard(&sketches[0], &sketches[1]).unwrap();
        assert!(j < 0.01, "Jaccard should be ~0.0, got {j}");
    }

    #[test]
    fn minhash_sketch_size_respected() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let sketch_size = 5;
        let sketches = gpu_minhash_cpu(&[seq.as_slice()], 3, sketch_size).unwrap();
        assert!(sketches[0].hashes.len() <= sketch_size);
        assert!(!sketches[0].hashes.is_empty());
        // Verify sorted
        for w in sketches[0].hashes.windows(2) {
            assert!(w[0] <= w[1]);
        }
    }

    #[test]
    fn minhash_empty_error() {
        assert!(gpu_minhash_cpu(&[], 4, 10).is_err());
        let seq = b"ACGT";
        assert!(gpu_minhash_cpu(&[seq.as_slice()], 0, 10).is_err());
        assert!(gpu_minhash_cpu(&[seq.as_slice()], 4, 0).is_err());
    }
}
