//! Embedded GPU shader sources for the Metal backend.

/// Metal Shading Language source for reduction kernels (sum, min, max).
pub const REDUCE_MSL: &str = include_str!("reduce.metal");

/// Metal Shading Language source for pairwise distance kernels.
pub const DISTANCE_MSL: &str = include_str!("distance.metal");

/// Metal Shading Language source for tiled matrix multiplication.
pub const MATMUL_MSL: &str = include_str!("matmul.metal");
