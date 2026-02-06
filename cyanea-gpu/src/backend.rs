//! Backend trait and core types for GPU compute abstraction.

use core::fmt;

use cyanea_core::{Result, Summarizable};

use crate::buffer::Buffer;

/// Identifies which compute backend is in use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BackendKind {
    /// CPU-only backend using plain iterators.
    Cpu,
    /// NVIDIA CUDA backend.
    Cuda,
    /// Apple Metal backend.
    Metal,
}

impl fmt::Display for BackendKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cpu => write!(f, "CPU"),
            Self::Cuda => write!(f, "CUDA"),
            Self::Metal => write!(f, "Metal"),
        }
    }
}

/// Information about a compute device.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DeviceInfo {
    /// Human-readable device name.
    pub name: String,
    /// Backend kind.
    pub kind: BackendKind,
    /// Total device memory in bytes (0 if unknown).
    pub total_memory: u64,
    /// Maximum parallelism (threads, warps, etc.).
    pub max_parallelism: usize,
}

impl Summarizable for DeviceInfo {
    fn summary(&self) -> String {
        format!(
            "{} ({}, {} bytes, {} parallel)",
            self.name, self.kind, self.total_memory, self.max_parallelism
        )
    }
}

/// Distance metric for GPU pairwise computations.
///
/// Separate from `cyanea_ml::DistanceMetric` to avoid circular dependencies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DistanceMetricGpu {
    /// Euclidean (L2) distance.
    Euclidean,
    /// Manhattan (L1) distance.
    Manhattan,
    /// Cosine distance (1 - cosine similarity).
    Cosine,
}

impl fmt::Display for DistanceMetricGpu {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Euclidean => write!(f, "Euclidean"),
            Self::Manhattan => write!(f, "Manhattan"),
            Self::Cosine => write!(f, "Cosine"),
        }
    }
}

/// Trait for GPU compute backends.
///
/// All backends must be `Send + Sync` for safe use across threads.
/// Methods use `&dyn Fn` closures to maintain object safety.
pub trait Backend: Send + Sync {
    // ── Resource info ──────────────────────────────────────────────

    /// Returns information about the underlying compute device.
    fn device_info(&self) -> DeviceInfo;

    // ── Buffer management ──────────────────────────────────────────

    /// Creates a buffer from a host slice, copying data to the device.
    fn buffer_from_slice(&self, data: &[f64]) -> Result<Buffer>;

    /// Creates a zero-initialized buffer of the given length.
    fn buffer_zeros(&self, len: usize) -> Result<Buffer>;

    /// Reads device buffer contents back to host memory.
    fn read_buffer(&self, buf: &Buffer) -> Result<Vec<f64>>;

    /// Overwrites device buffer contents from a host slice.
    ///
    /// Returns an error if `data.len() != buf.len()`.
    fn write_buffer(&self, buf: &mut Buffer, data: &[f64]) -> Result<()>;

    /// Returns the number of `f64` elements in the buffer.
    fn buffer_len(&self, buf: &Buffer) -> usize;

    // ── Parallel primitives ────────────────────────────────────────

    /// Applies `f` element-wise to `input` and writes results to `output`.
    ///
    /// `output` must have the same length as `input`.
    fn elementwise_map(
        &self,
        input: &Buffer,
        output: &mut Buffer,
        f: &dyn Fn(f64) -> f64,
    ) -> Result<()>;

    /// Reduces a buffer to a single sum.
    fn reduce_sum(&self, buf: &Buffer) -> Result<f64>;

    /// Reduces a buffer to its minimum value.
    fn reduce_min(&self, buf: &Buffer) -> Result<f64>;

    /// Reduces a buffer to its maximum value.
    fn reduce_max(&self, buf: &Buffer) -> Result<f64>;

    // ── Bioinformatics operations ──────────────────────────────────

    /// Computes an n×n pairwise distance matrix.
    ///
    /// `data` is a flat row-major matrix of shape `(n, dim)`.
    /// Returns a flat n×n distance matrix.
    fn pairwise_distance_matrix(
        &self,
        data: &Buffer,
        n: usize,
        dim: usize,
        metric: DistanceMetricGpu,
    ) -> Result<Buffer>;

    /// Multiplies two matrices: C = A × B.
    ///
    /// `a` is row-major (m, k), `b` is row-major (k, n).
    /// Returns a row-major (m, n) buffer.
    fn matrix_multiply(
        &self,
        a: &Buffer,
        b: &Buffer,
        m: usize,
        k: usize,
        n: usize,
    ) -> Result<Buffer>;

    /// Computes all-pairs results for `n` items each of length `item_len`.
    ///
    /// `items` is a flat row-major matrix of shape `(n, item_len)`.
    /// `f` is called with two item slices and returns a scalar.
    /// Returns a flat n×n matrix.
    fn batch_pairwise(
        &self,
        items: &Buffer,
        n: usize,
        item_len: usize,
        f: &dyn Fn(&[f64], &[f64]) -> f64,
    ) -> Result<Buffer>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backend_kind_display() {
        assert_eq!(BackendKind::Cpu.to_string(), "CPU");
        assert_eq!(BackendKind::Cuda.to_string(), "CUDA");
        assert_eq!(BackendKind::Metal.to_string(), "Metal");
    }

    #[test]
    fn backend_kind_debug_clone() {
        let kind = BackendKind::Cpu;
        let cloned = kind;
        assert_eq!(kind, cloned);
        assert_eq!(format!("{:?}", kind), "Cpu");
    }

    #[test]
    fn device_info_summary() {
        let info = DeviceInfo {
            name: "Test Device".to_string(),
            kind: BackendKind::Cpu,
            total_memory: 1024,
            max_parallelism: 8,
        };
        assert_eq!(info.summary(), "Test Device (CPU, 1024 bytes, 8 parallel)");
    }

    #[test]
    fn distance_metric_gpu_equality() {
        assert_eq!(DistanceMetricGpu::Euclidean, DistanceMetricGpu::Euclidean);
        assert_ne!(DistanceMetricGpu::Euclidean, DistanceMetricGpu::Manhattan);
        assert_ne!(DistanceMetricGpu::Manhattan, DistanceMetricGpu::Cosine);
    }

    #[test]
    fn distance_metric_gpu_display() {
        assert_eq!(DistanceMetricGpu::Euclidean.to_string(), "Euclidean");
        assert_eq!(DistanceMetricGpu::Manhattan.to_string(), "Manhattan");
        assert_eq!(DistanceMetricGpu::Cosine.to_string(), "Cosine");
    }
}
