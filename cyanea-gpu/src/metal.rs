//! Metal backend stub.
//!
//! This module provides a placeholder `MetalBackend` type that will be replaced
//! with a real Metal implementation in the future. The planned architecture:
//!
//! - Use `metal-rs` crate for Apple GPU access
//! - Metal Performance Shaders for matrix operations
//! - Compute command encoders for custom kernels
//! - Shared memory between CPU and GPU via unified memory architecture

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
use crate::buffer::Buffer;

/// Stub Metal backend. All methods return errors.
pub struct MetalBackend {
    _private: (),
}

impl MetalBackend {
    /// Attempts to create a Metal backend.
    ///
    /// # Panics
    ///
    /// Always panics — Metal support is not yet implemented.
    pub fn new() -> Self {
        panic!("Metal backend is not yet implemented")
    }
}

fn not_impl() -> CyaneaError {
    CyaneaError::Other("Metal backend is not yet implemented".to_string())
}

impl Backend for MetalBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: "Metal (stub)".to_string(),
            kind: BackendKind::Metal,
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
