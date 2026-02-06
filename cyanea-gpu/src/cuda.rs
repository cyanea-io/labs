//! CUDA backend stub.
//!
//! This module provides a placeholder `CudaBackend` type that will be replaced
//! with a real CUDA implementation in the future. The planned architecture:
//!
//! - Use `cudarc` or raw CUDA driver API bindings
//! - Async host↔device memory transfers via CUDA streams
//! - Custom kernels for pairwise distance, matrix multiply, reductions
//! - Shared memory tiling for matrix operations
//! - Warp-level primitives for reductions

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
use crate::buffer::Buffer;

/// Stub CUDA backend. All methods return errors.
pub struct CudaBackend {
    _private: (),
}

impl CudaBackend {
    /// Attempts to create a CUDA backend.
    ///
    /// # Panics
    ///
    /// Always panics — CUDA support is not yet implemented.
    pub fn new() -> Self {
        panic!("CUDA backend is not yet implemented")
    }
}

fn not_impl() -> CyaneaError {
    CyaneaError::Other("CUDA backend is not yet implemented".to_string())
}

impl Backend for CudaBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: "CUDA (stub)".to_string(),
            kind: BackendKind::Cuda,
            total_memory: 0,
            max_parallelism: 0,
        }
    }

    fn buffer_from_slice(&self, _data: &[f64]) -> Result<Buffer> {
        Err(not_impl())
    }

    fn buffer_zeros(&self, _len: usize) -> Result<Buffer> {
        Err(not_impl())
    }

    fn read_buffer(&self, _buf: &Buffer) -> Result<Vec<f64>> {
        Err(not_impl())
    }

    fn write_buffer(&self, _buf: &mut Buffer, _data: &[f64]) -> Result<()> {
        Err(not_impl())
    }

    fn buffer_len(&self, buf: &Buffer) -> usize {
        buf.len()
    }

    fn elementwise_map(
        &self,
        _input: &Buffer,
        _output: &mut Buffer,
        _f: &dyn Fn(f64) -> f64,
    ) -> Result<()> {
        Err(not_impl())
    }

    fn reduce_sum(&self, _buf: &Buffer) -> Result<f64> {
        Err(not_impl())
    }

    fn reduce_min(&self, _buf: &Buffer) -> Result<f64> {
        Err(not_impl())
    }

    fn reduce_max(&self, _buf: &Buffer) -> Result<f64> {
        Err(not_impl())
    }

    fn pairwise_distance_matrix(
        &self,
        _data: &Buffer,
        _n: usize,
        _dim: usize,
        _metric: DistanceMetricGpu,
    ) -> Result<Buffer> {
        Err(not_impl())
    }

    fn matrix_multiply(
        &self,
        _a: &Buffer,
        _b: &Buffer,
        _m: usize,
        _k: usize,
        _n: usize,
    ) -> Result<Buffer> {
        Err(not_impl())
    }

    fn batch_pairwise(
        &self,
        _items: &Buffer,
        _n: usize,
        _item_len: usize,
        _f: &dyn Fn(&[f64], &[f64]) -> f64,
    ) -> Result<Buffer> {
        Err(not_impl())
    }
}
