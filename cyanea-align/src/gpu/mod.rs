//! GPU-accelerated batch alignment.
//!
//! Provides banded affine-gap alignment on Metal (macOS) and CUDA (NVIDIA) GPUs.
//! Pairs that exceed the bandwidth or when no GPU is available fall back to CPU.

pub mod common;

#[cfg(feature = "metal")]
pub mod metal_align;

#[cfg(feature = "cuda")]
pub mod cuda_align;

pub use common::GpuAlignConfig;

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentMode, AlignmentResult};
use cyanea_core::{CyaneaError, Result};

/// Available GPU backend for alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpuBackend {
    /// Apple Metal (macOS / Apple Silicon).
    Metal,
    /// NVIDIA CUDA.
    Cuda,
}

/// Detect which GPU backends are available at runtime.
pub fn available_backends() -> Vec<GpuBackend> {
    #[allow(unused_mut)]
    let mut backends = Vec::new();

    #[cfg(feature = "metal")]
    {
        if metal_align::MetalAligner::new().is_ok() {
            backends.push(GpuBackend::Metal);
        }
    }

    #[cfg(feature = "cuda")]
    {
        if cuda_align::CudaAligner::new().is_ok() {
            backends.push(GpuBackend::Cuda);
        }
    }

    backends
}

/// Batch-align sequence pairs on the GPU with default configuration.
///
/// Automatically selects the best available GPU backend (Metal on macOS,
/// CUDA on NVIDIA). Falls back to CPU alignment if no GPU is available
/// or the batch is too small.
///
/// # Errors
///
/// Returns an error if any alignment fails.
pub fn align_batch_gpu(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
) -> Result<Vec<AlignmentResult>> {
    align_batch_gpu_with_config(pairs, mode, scoring, &GpuAlignConfig::default())
}

/// Batch-align sequence pairs on the GPU with explicit configuration.
///
/// # Errors
///
/// Returns an error if any alignment fails.
pub fn align_batch_gpu_with_config(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
    config: &GpuAlignConfig,
) -> Result<Vec<AlignmentResult>> {
    if pairs.is_empty() {
        return Ok(Vec::new());
    }

    // Below threshold: CPU is faster
    if pairs.len() < config.min_pairs_for_gpu {
        return cpu_fallback(pairs, mode, scoring);
    }

    // Try GPU backends in preference order
    #[cfg(feature = "metal")]
    {
        match metal_align::MetalAligner::new() {
            Ok(aligner) => return aligner.align_batch(pairs, mode, scoring, config),
            Err(_) => {} // fall through
        }
    }

    #[cfg(feature = "cuda")]
    {
        match cuda_align::CudaAligner::new() {
            Ok(aligner) => return aligner.align_batch(pairs, mode, scoring, config),
            Err(_) => {} // fall through
        }
    }

    // No GPU available
    cpu_fallback(pairs, mode, scoring)
}

/// Batch-align on a specific GPU backend.
///
/// # Errors
///
/// Returns an error if the requested backend is not available or any alignment fails.
pub fn align_batch_on(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
    backend: GpuBackend,
    config: &GpuAlignConfig,
) -> Result<Vec<AlignmentResult>> {
    if pairs.is_empty() {
        return Ok(Vec::new());
    }

    // Suppress unused-variable warnings when no GPU features are enabled.
    let _ = (mode, scoring, config);

    match backend {
        GpuBackend::Metal => {
            #[cfg(feature = "metal")]
            {
                let aligner = metal_align::MetalAligner::new()?;
                return aligner.align_batch(pairs, mode, scoring, config);
            }
            #[cfg(not(feature = "metal"))]
            Err(CyaneaError::Other("Metal feature not enabled".into()))
        }
        GpuBackend::Cuda => {
            #[cfg(feature = "cuda")]
            {
                let aligner = cuda_align::CudaAligner::new()?;
                return aligner.align_batch(pairs, mode, scoring, config);
            }
            #[cfg(not(feature = "cuda"))]
            Err(CyaneaError::Other("CUDA feature not enabled".into()))
        }
    }
}

/// CPU fallback: align each pair using the standard CPU implementation.
fn cpu_fallback(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
) -> Result<Vec<AlignmentResult>> {
    pairs
        .iter()
        .map(|(q, t)| crate::align(q, t, mode, scoring))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringMatrix;

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn gpu_fallback_empty() {
        let pairs: Vec<(&[u8], &[u8])> = vec![];
        let results = align_batch_gpu(&pairs, AlignmentMode::Global, &dna_scheme()).unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn gpu_fallback_small_batch() {
        // Below min_pairs_for_gpu threshold → CPU fallback
        let pairs: Vec<(&[u8], &[u8])> = vec![(b"ACGT", b"ACGT"), (b"AAA", b"TTT")];
        let results = align_batch_gpu(&pairs, AlignmentMode::Global, &dna_scheme()).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].score, 8);
    }

    #[test]
    fn available_backends_returns_list() {
        let backends = available_backends();
        // Regardless of hardware, this should not panic
        assert!(backends.len() <= 2);
    }

    #[test]
    fn config_default() {
        let cfg = GpuAlignConfig::default();
        assert_eq!(cfg.max_bandwidth, 128);
        assert_eq!(cfg.min_pairs_for_gpu, 64);
    }

    #[test]
    fn cpu_fallback_correctness() {
        let pairs: Vec<(&[u8], &[u8])> =
            vec![(b"ACGT", b"ACGT"), (b"AAAA", b"TTTT"), (b"ACGT", b"ACT")];
        let results = cpu_fallback(&pairs, AlignmentMode::Global, &dna_scheme()).unwrap();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].score, 8);
    }
}
