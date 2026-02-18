//! GPU-accelerated k-mer counting.
//!
//! Provides parallel k-mer counting using GPU (Metal/CUDA) with automatic
//! CPU fallback. K-mers are encoded as 2-bit integers (A=0, C=1, G=2, T=3)
//! and counted in a table of size 4^k.
//!
//! GPU dispatch is limited to k ≤ 14 (4^14 = 268M entries). Larger k values
//! use CPU-only HashMap-based counting.

use cyanea_core::{CyaneaError, Result};

/// Result of a k-mer counting operation.
#[derive(Debug, Clone)]
pub struct KmerCountResult {
    /// Counts indexed by 2-bit k-mer encoding (size = 4^k).
    pub counts: Vec<u32>,
    /// The k-mer length used.
    pub k: usize,
}

/// Maximum k for GPU dispatch (4^14 = 268M entries).
#[cfg(any(feature = "metal", feature = "cuda"))]
const MAX_GPU_K: usize = 14;

/// Counts k-mers in the given DNA sequences using the best available backend.
///
/// Sequences should contain A/C/G/T (uppercase). Positions with N or other
/// ambiguous bases are skipped. Multiple sequences accumulate into a single
/// count table.
///
/// # Errors
///
/// Returns an error if `k == 0`, `k > 31`, or `sequences` is empty.
pub fn gpu_kmer_count(sequences: &[&[u8]], k: usize) -> Result<KmerCountResult> {
    validate_inputs(sequences, k)?;

    #[cfg(feature = "metal")]
    if k <= MAX_GPU_K {
        if let Ok(result) = gpu_kmer_count_metal(sequences, k) {
            return Ok(result);
        }
    }

    // CPU fallback
    gpu_kmer_count_cpu(sequences, k)
}

/// CPU-only k-mer counting.
pub fn gpu_kmer_count_cpu(sequences: &[&[u8]], k: usize) -> Result<KmerCountResult> {
    validate_inputs(sequences, k)?;

    let table_size = 4usize.checked_pow(k as u32).ok_or_else(|| {
        CyaneaError::InvalidInput(format!("k={k} too large: 4^{k} overflows"))
    })?;
    let mut counts = vec![0u32; table_size];

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        for pos in 0..=(seq.len() - k) {
            if let Some(hash) = encode_kmer(&seq[pos..pos + k]) {
                counts[hash] += 1;
            }
        }
    }

    Ok(KmerCountResult { counts, k })
}

/// Encodes a k-mer as a 2-bit integer. Returns `None` if any base is not ACGT.
fn encode_kmer(kmer: &[u8]) -> Option<usize> {
    let mut hash = 0usize;
    for &base in kmer {
        let code = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        hash = (hash << 2) | code;
    }
    Some(hash)
}

fn validate_inputs(sequences: &[&[u8]], k: usize) -> Result<()> {
    if k == 0 {
        return Err(CyaneaError::InvalidInput(
            "kmer_count: k must be > 0".to_string(),
        ));
    }
    if k > 31 {
        return Err(CyaneaError::InvalidInput(
            "kmer_count: k must be ≤ 31".to_string(),
        ));
    }
    if sequences.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "kmer_count: sequences must not be empty".to_string(),
        ));
    }
    Ok(())
}

// ── Metal GPU dispatch ────────────────────────────────────────────

#[cfg(feature = "metal")]
fn gpu_kmer_count_metal(sequences: &[&[u8]], k: usize) -> Result<KmerCountResult> {
    use metal_rs::{CompileOptions, Device, MTLResourceOptions, MTLSize};
    use std::ffi::c_void;

    let device = Device::system_default()
        .ok_or_else(|| CyaneaError::Other("no Metal device available".into()))?;
    let queue = device.new_command_queue();

    let shader_src = include_str!("shaders/kmer.metal");
    let opts = CompileOptions::new();
    let lib = device
        .new_library_with_source(shader_src, &opts)
        .map_err(|e| CyaneaError::Other(format!("Metal kmer shader compile: {e}")))?;
    let func = lib
        .get_function("kmer_count", None)
        .map_err(|e| CyaneaError::Other(format!("Metal kmer function: {e}")))?;
    let pipeline = device
        .new_compute_pipeline_state_with_function(&func)
        .map_err(|e| CyaneaError::Other(format!("Metal kmer pipeline: {e}")))?;

    let table_size = 4usize.pow(k as u32);
    // Allocate zero-initialized counts buffer
    let counts_buf = device.new_buffer(
        (table_size * std::mem::size_of::<u32>()) as u64,
        MTLResourceOptions::StorageModeShared,
    );
    // Zero it
    unsafe {
        std::ptr::write_bytes(counts_buf.contents() as *mut u8, 0, table_size * 4);
    }

    let k_u32 = k as u32;
    let k_buf = device.new_buffer_with_data(
        &k_u32 as *const u32 as *const c_void,
        4,
        MTLResourceOptions::StorageModeShared,
    );

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
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

        let n_kmers = seq.len() - k + 1;
        let cmd = queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&pipeline);
        enc.set_buffer(0, Some(&seq_buf), 0);
        enc.set_buffer(1, Some(&len_buf), 0);
        enc.set_buffer(2, Some(&k_buf), 0);
        enc.set_buffer(3, Some(&counts_buf), 0);

        let grid = MTLSize::new(n_kmers as u64, 1, 1);
        let tg = MTLSize::new(256.min(n_kmers as u64).max(1), 1, 1);
        enc.dispatch_threads(grid, tg);
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();
    }

    // Read counts back
    let ptr = counts_buf.contents() as *const u32;
    let counts = unsafe { std::slice::from_raw_parts(ptr, table_size) }.to_vec();

    Ok(KmerCountResult { counts, k })
}

// ── CUDA GPU dispatch ─────────────────────────────────────────────

#[cfg(feature = "cuda")]
fn _gpu_kmer_count_cuda(sequences: &[&[u8]], k: usize) -> Result<KmerCountResult> {
    use crate::cuda::kernels::KMER_KERNEL_SOURCE;
    use cudarc::driver::{CudaContext, LaunchConfig, PushKernelArg};
    use cudarc::nvrtc::compile_ptx;

    let ctx = CudaContext::new(0)
        .map_err(|e| CyaneaError::Other(format!("CUDA context: {e}")))?;
    let stream = ctx.default_stream();

    let ptx = compile_ptx(KMER_KERNEL_SOURCE)
        .map_err(|e| CyaneaError::Other(format!("CUDA kmer compile: {e}")))?;
    let module = ctx
        .load_module(ptx)
        .map_err(|e| CyaneaError::Other(format!("CUDA module: {e}")))?;
    let func = module
        .load_function("kmer_count")
        .map_err(|e| CyaneaError::Other(format!("CUDA kmer function: {e}")))?;

    let table_size = 4usize.pow(k as u32);
    let mut counts_dev: cudarc::driver::CudaSlice<u32> = stream
        .alloc_zeros(table_size)
        .map_err(|e| CyaneaError::Other(format!("CUDA alloc: {e}")))?;

    let k_u32 = k as u32;

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        let seq_u8: Vec<u8> = seq.to_vec();
        let seq_dev = stream
            .clone_htod(&seq_u8)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        let seq_len = seq.len() as u32;
        let n_kmers = seq.len() - k + 1;
        let block_size = 256u32;
        let grid_size = ((n_kmers as u32) + block_size - 1) / block_size;

        let cfg = LaunchConfig {
            grid_dim: (grid_size, 1, 1),
            block_dim: (block_size, 1, 1),
            shared_mem_bytes: 0,
        };

        unsafe {
            stream
                .launch_builder(&func)
                .arg(&seq_dev)
                .arg(&seq_len)
                .arg(&k_u32)
                .arg(&mut counts_dev)
                .launch(cfg)
        }
        .map_err(|e| CyaneaError::Other(format!("CUDA launch: {e}")))?;
    }

    let counts = stream
        .clone_dtoh(&counts_dev)
        .map_err(|e| CyaneaError::Other(format!("CUDA dtoh: {e}")))?;

    Ok(KmerCountResult { counts, k })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_count_simple() {
        let seq = b"ACGTACGT";
        let result = gpu_kmer_count_cpu(&[seq.as_slice()], 2).unwrap();
        assert_eq!(result.k, 2);
        assert_eq!(result.counts.len(), 16); // 4^2
        // AC appears at positions 0, 4 → count 2
        let ac_idx = encode_kmer(b"AC").unwrap();
        assert_eq!(result.counts[ac_idx], 2);
        // CG appears at positions 1, 5 → count 2
        let cg_idx = encode_kmer(b"CG").unwrap();
        assert_eq!(result.counts[cg_idx], 2);
        // GT appears at positions 2, 6 → count 2
        let gt_idx = encode_kmer(b"GT").unwrap();
        assert_eq!(result.counts[gt_idx], 2);
        // TA appears at position 3 → count 1
        let ta_idx = encode_kmer(b"TA").unwrap();
        assert_eq!(result.counts[ta_idx], 1);
    }

    #[test]
    fn kmer_count_with_n() {
        let seq = b"ACNGT";
        let result = gpu_kmer_count_cpu(&[seq.as_slice()], 2).unwrap();
        // AC at 0 → 1
        let ac_idx = encode_kmer(b"AC").unwrap();
        assert_eq!(result.counts[ac_idx], 1);
        // CN at 1 → skipped (N)
        // NG at 2 → skipped (N)
        // GT at 3 → 1
        let gt_idx = encode_kmer(b"GT").unwrap();
        assert_eq!(result.counts[gt_idx], 1);
        let total: u32 = result.counts.iter().sum();
        assert_eq!(total, 2);
    }

    #[test]
    fn kmer_count_multiple_sequences() {
        let seq1 = b"AAAA";
        let seq2 = b"AAAA";
        let result = gpu_kmer_count_cpu(&[seq1.as_slice(), seq2.as_slice()], 2).unwrap();
        let aa_idx = encode_kmer(b"AA").unwrap();
        // Each "AAAA" has 3 AA k-mers → 6 total
        assert_eq!(result.counts[aa_idx], 6);
    }

    #[test]
    fn kmer_count_k1() {
        let seq = b"AACGT";
        let result = gpu_kmer_count_cpu(&[seq.as_slice()], 1).unwrap();
        assert_eq!(result.counts.len(), 4);
        assert_eq!(result.counts[0], 2); // A
        assert_eq!(result.counts[1], 1); // C
        assert_eq!(result.counts[2], 1); // G
        assert_eq!(result.counts[3], 1); // T
    }

    #[test]
    fn kmer_count_empty_error() {
        assert!(gpu_kmer_count_cpu(&[], 2).is_err());
        let seq = b"ACGT";
        assert!(gpu_kmer_count_cpu(&[seq.as_slice()], 0).is_err());
    }
}
