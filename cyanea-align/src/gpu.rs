//! GPU-accelerated batch alignment backends (not yet implemented).
//!
//! This module will provide CUDA and Metal backends for high-throughput pairwise
//! alignment of large sequence batches on the GPU.
//!
//! # Planned features
//!
//! - **CUDA backend** (`cuda` feature) — Launch grid of thread blocks where each
//!   block aligns one sequence pair using shared-memory tiled DP. Target throughput:
//!   millions of short-read alignments per second on modern NVIDIA GPUs.
//!
//! - **Metal backend** (`metal` feature) — Equivalent compute-shader implementation
//!   for Apple Silicon, using Metal Performance Shaders and threadgroup memory.
//!
//! - **Batch alignment** — Accept `Vec<(query, target)>` pairs, upload to device
//!   memory, execute in parallel, and download results. The host-side API mirrors
//!   [`crate::batch::align_batch`] for drop-in replacement.
//!
//! - **Streaming** — For datasets larger than device memory, automatically partition
//!   into chunks and pipeline upload/compute/download across multiple command buffers
//!   or CUDA streams.
//!
//! - **Mixed precision** — Support both i32 and i16 score accumulators, with
//!   automatic overflow detection and fallback. 16-bit mode doubles throughput on
//!   architectures with wide SIMD ALUs (e.g. Tensor Cores).
//!
//! # Architecture
//!
//! ```text
//! gpu::align_batch_gpu()
//!   ├── cuda::CudaAligner   (feature = "cuda")
//!   └── metal::MetalAligner  (feature = "metal")
//! ```

/// CUDA alignment backend.
#[cfg(feature = "cuda")]
pub mod cuda {
    // Placeholder — will contain CudaAligner, kernel launch wrappers, etc.
}

/// Metal alignment backend.
#[cfg(feature = "metal")]
pub mod metal {
    // Placeholder — will contain MetalAligner, compute pipeline setup, etc.
}
