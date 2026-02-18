//! Embedded GPU shader sources for the Metal and WebGPU backends.

/// Metal Shading Language source for reduction kernels (sum, min, max).
pub const REDUCE_MSL: &str = include_str!("reduce.metal");

/// Metal Shading Language source for pairwise distance kernels.
pub const DISTANCE_MSL: &str = include_str!("distance.metal");

/// Metal Shading Language source for tiled matrix multiplication.
pub const MATMUL_MSL: &str = include_str!("matmul.metal");

/// Metal Shading Language source for GPU k-mer counting.
pub const KMER_MSL: &str = include_str!("kmer.metal");

/// Metal Shading Language source for GPU Smith-Waterman protein alignment.
pub const SW_MSL: &str = include_str!("sw.metal");

/// Metal Shading Language source for GPU MinHash k-mer hashing.
pub const MINHASH_MSL: &str = include_str!("minhash.metal");

// ── WGSL shaders (WebGPU) ─────────────────────────────────────────

/// WGSL source for reduction kernels (sum, min, max).
#[cfg(feature = "wgpu")]
pub const REDUCE_WGSL: &str = include_str!("reduce.wgsl");

/// WGSL source for pairwise distance kernels.
#[cfg(feature = "wgpu")]
pub const DISTANCE_WGSL: &str = include_str!("distance.wgsl");

/// WGSL source for tiled matrix multiplication.
#[cfg(feature = "wgpu")]
pub const MATMUL_WGSL: &str = include_str!("matmul.wgsl");
