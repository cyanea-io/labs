//! GPU compute abstraction for the Cyanea bioinformatics ecosystem.
//!
//! Provides a [`Backend`] trait that abstracts over CPU, CUDA, Metal,
//! and WebGPU compute backends. The [`CpuBackend`] is always available as
//! a reference implementation and fallback. GPU backends are feature-gated
//! behind `metal` (Apple Silicon), `cuda` (NVIDIA), and `wgpu` (WebGPU).
//!
//! # Quick start
//!
//! ```
//! use cyanea_gpu::{CpuBackend, Backend};
//!
//! let backend = CpuBackend::new();
//! let buf = backend.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
//! let sum = backend.reduce_sum(&buf).unwrap();
//! assert!((sum - 10.0).abs() < 1e-12);
//! ```
//!
//! # Auto-selecting a backend
//!
//! ```
//! use cyanea_gpu::auto_backend;
//!
//! let backend = auto_backend();
//! let info = backend.device_info();
//! println!("Using: {}", info.name);
//! ```

pub mod backend;
pub mod buffer;
pub mod cpu;
pub mod ops;

// Domain-specific GPU modules (standalone, not part of Backend trait)
pub mod kmer;
pub mod minhash_gpu;
pub mod smith_waterman;

#[cfg(any(feature = "metal", feature = "wgpu"))]
pub mod shaders;

#[cfg(feature = "cuda")]
pub mod cuda;

#[cfg(feature = "metal")]
pub mod metal;

#[cfg(feature = "wgpu")]
pub mod wgpu_backend;

pub use backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
pub use buffer::Buffer;
pub use cpu::CpuBackend;
pub use kmer::{gpu_kmer_count, gpu_kmer_count_cpu, KmerCountResult};
pub use minhash_gpu::{gpu_minhash, gpu_minhash_cpu, gpu_minhash_jaccard, MinHashSketch};
pub use ops::tiled_pairwise_distance;
pub use smith_waterman::{gpu_smith_waterman_batch, gpu_sw_cpu, SwResult, BLOSUM62_24};

#[cfg(feature = "cuda")]
pub use cuda::CudaBackend;

#[cfg(feature = "metal")]
pub use metal::MetalBackend;

#[cfg(feature = "wgpu")]
pub use wgpu_backend::WgpuBackend;

/// Returns the best available backend for the current platform.
///
/// Selection order (compile-time feature gates, runtime availability):
/// 1. `metal` — Apple Metal (if feature enabled and GPU available)
/// 2. `cuda` — NVIDIA CUDA (if feature enabled and GPU available)
/// 3. `wgpu` — WebGPU (if feature enabled and adapter available)
/// 4. CPU fallback (always available)
pub fn auto_backend() -> Box<dyn Backend> {
    #[cfg(feature = "metal")]
    if let Ok(m) = MetalBackend::new() {
        return Box::new(m);
    }
    #[cfg(feature = "cuda")]
    if let Ok(c) = CudaBackend::new() {
        return Box::new(c);
    }
    #[cfg(feature = "wgpu")]
    if let Ok(w) = WgpuBackend::new() {
        return Box::new(w);
    }
    Box::new(CpuBackend::new())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_backend_returns_cpu() {
        // Without GPU features, always returns CPU.
        #[cfg(not(any(feature = "metal", feature = "cuda", feature = "wgpu")))]
        {
            let b = auto_backend();
            assert_eq!(b.device_info().kind, BackendKind::Cpu);
        }
    }

    #[test]
    fn auto_backend_usable_through_trait() {
        let b = auto_backend();
        let buf = b.buffer_from_slice(&[10.0, 20.0, 30.0]).unwrap();
        let sum = b.reduce_sum(&buf).unwrap();
        assert!((sum - 60.0).abs() < 0.1);
    }
}
