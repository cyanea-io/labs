//! GPU compute abstraction for the Cyanea bioinformatics ecosystem.
//!
//! Provides a [`Backend`] trait that abstracts over CPU, CUDA, and Metal
//! compute backends. The [`CpuBackend`] is always available as a reference
//! implementation and fallback. GPU backends are feature-gated stubs that
//! will be replaced with real implementations in the future.
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

#[cfg(feature = "cuda")]
pub mod cuda;

#[cfg(feature = "metal")]
pub mod metal;

pub use backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
pub use buffer::Buffer;
pub use cpu::CpuBackend;

#[cfg(feature = "cuda")]
pub use cuda::CudaBackend;

#[cfg(feature = "metal")]
pub use metal::MetalBackend;

/// Returns the best available backend for the current platform.
///
/// Selection order (compile-time feature gates):
/// 1. `metal` — Apple Metal (if feature enabled)
/// 2. `cuda` — NVIDIA CUDA (if feature enabled)
/// 3. CPU fallback (always available)
///
/// Note: Since GPU backends are currently stubs that panic on construction,
/// this always returns the CPU backend in practice.
pub fn auto_backend() -> Box<dyn Backend> {
    Box::new(CpuBackend::new())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_backend_returns_cpu() {
        let b = auto_backend();
        assert_eq!(b.device_info().kind, BackendKind::Cpu);
    }

    #[test]
    fn auto_backend_usable_through_trait() {
        let b = auto_backend();
        let buf = b.buffer_from_slice(&[10.0, 20.0, 30.0]).unwrap();
        let sum = b.reduce_sum(&buf).unwrap();
        assert!((sum - 60.0).abs() < 1e-12);
    }
}
