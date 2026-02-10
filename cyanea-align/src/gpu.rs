//! GPU-accelerated batch alignment backends.
//!
//! Provides feature-gated CUDA and Metal backends for high-throughput pairwise
//! alignment on GPUs. When neither GPU feature is enabled, the public
//! [`align_batch_gpu`] function falls back to CPU-based batch alignment.
//!
//! # Architecture
//!
//! ```text
//! gpu::align_batch_gpu()
//!   ├── cuda::CudaAligner   (feature = "cuda")
//!   ├── metal::MetalAligner  (feature = "metal")
//!   └── CPU fallback          (default)
//! ```

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentMode, AlignmentResult};

/// Which GPU backend to prefer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpuBackend {
    /// Automatically select the best available backend.
    Auto,
    /// Force CUDA (requires `cuda` feature).
    Cuda,
    /// Force Metal (requires `metal` feature).
    Metal,
    /// Force CPU fallback.
    Cpu,
}

/// Batch-align sequence pairs, dispatching to the best available backend.
///
/// Falls back to CPU if no GPU backend is compiled in or if `backend` is
/// [`GpuBackend::Cpu`].
///
/// # Errors
///
/// Returns an error if any sequence pair fails to align or if a requested
/// GPU backend is not available.
pub fn align_batch_gpu(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
    backend: GpuBackend,
) -> Result<Vec<AlignmentResult>> {
    match backend {
        #[cfg(feature = "cuda")]
        GpuBackend::Cuda => {
            Err(CyaneaError::Other("CUDA alignment backend is not yet implemented".into()))
        }
        #[cfg(not(feature = "cuda"))]
        GpuBackend::Cuda => {
            Err(CyaneaError::Other("CUDA feature is not enabled".into()))
        }
        #[cfg(feature = "metal")]
        GpuBackend::Metal => {
            Err(CyaneaError::Other("Metal alignment backend is not yet implemented".into()))
        }
        #[cfg(not(feature = "metal"))]
        GpuBackend::Metal => {
            Err(CyaneaError::Other("Metal feature is not enabled".into()))
        }
        GpuBackend::Auto | GpuBackend::Cpu => {
            cpu_fallback(pairs, mode, scoring)
        }
    }
}

/// Returns which GPU backends are compiled in.
pub fn available_backends() -> Vec<GpuBackend> {
    #[allow(unused_mut)]
    let mut backends = vec![GpuBackend::Cpu];
    #[cfg(feature = "cuda")]
    backends.push(GpuBackend::Cuda);
    #[cfg(feature = "metal")]
    backends.push(GpuBackend::Metal);
    backends
}

fn cpu_fallback(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
) -> Result<Vec<AlignmentResult>> {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        pairs
            .par_iter()
            .map(|(q, t)| crate::align(q, t, mode, scoring))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    pairs
        .iter()
        .map(|(q, t)| crate::align(q, t, mode, scoring))
        .collect()
}

/// CUDA alignment backend.
#[cfg(feature = "cuda")]
pub mod cuda {
    // Placeholder — requires CUDA toolkit and `cudarc` or raw driver bindings.
}

/// Metal alignment backend.
#[cfg(feature = "metal")]
pub mod metal {
    // Placeholder — requires `metal-rs` and Apple GPU access.
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoringMatrix, ScoringScheme};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn batch_gpu_cpu_fallback() {
        let pairs: Vec<(&[u8], &[u8])> = vec![
            (b"ACGT", b"ACGT"),
            (b"ACGT", b"ACGA"),
        ];
        let results = align_batch_gpu(&pairs, AlignmentMode::Global, &dna_scheme(), GpuBackend::Cpu).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].score, 8);
    }

    #[test]
    fn batch_gpu_auto() {
        let pairs: Vec<(&[u8], &[u8])> = vec![(b"ACGT", b"ACGT")];
        let results = align_batch_gpu(&pairs, AlignmentMode::Global, &dna_scheme(), GpuBackend::Auto).unwrap();
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn available_backends_includes_cpu() {
        let backends = available_backends();
        assert!(backends.contains(&GpuBackend::Cpu));
    }

    #[test]
    fn cuda_not_available_without_feature() {
        #[cfg(not(feature = "cuda"))]
        {
            let pairs: Vec<(&[u8], &[u8])> = vec![(b"ACGT", b"ACGT")];
            let result = align_batch_gpu(&pairs, AlignmentMode::Global, &dna_scheme(), GpuBackend::Cuda);
            assert!(result.is_err());
        }
    }

    #[test]
    fn metal_not_available_without_feature() {
        #[cfg(not(feature = "metal"))]
        {
            let pairs: Vec<(&[u8], &[u8])> = vec![(b"ACGT", b"ACGT")];
            let result = align_batch_gpu(&pairs, AlignmentMode::Global, &dna_scheme(), GpuBackend::Metal);
            assert!(result.is_err());
        }
    }
}
