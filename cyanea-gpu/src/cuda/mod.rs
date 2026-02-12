//! CUDA backend — GPU compute via NVIDIA CUDA with runtime compilation.
//!
//! Uses `cudarc` for safe CUDA driver API access and NVRTC for runtime
//! compilation of CUDA C kernels to PTX. All operations use `double` (f64)
//! since NVIDIA GPUs have native f64 support.

pub mod kernels;

use std::sync::Arc;

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::compile_ptx;

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
use crate::buffer::Buffer;

/// NVIDIA CUDA compute backend.
///
/// Compiles all CUDA kernels at construction via NVRTC and caches
/// function handles. Uses the cudarc crate which dynamically loads
/// libcuda — no compile-time CUDA SDK dependency.
pub struct CudaBackend {
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    reduce_sum_fn: CudaFunction,
    reduce_min_fn: CudaFunction,
    reduce_max_fn: CudaFunction,
    pairwise_euclidean_fn: CudaFunction,
    pairwise_manhattan_fn: CudaFunction,
    pairwise_cosine_fn: CudaFunction,
    matmul_fn: CudaFunction,
}

impl CudaBackend {
    /// Creates a CUDA backend, compiling all kernels via NVRTC.
    ///
    /// # Errors
    ///
    /// Returns an error if no CUDA device is available, the driver cannot
    /// be loaded, or kernel compilation fails.
    pub fn new() -> Result<Self> {
        let ctx = CudaContext::new(0)
            .map_err(|e| CyaneaError::Other(format!("CUDA context init: {e}")))?;
        let stream = ctx.default_stream();

        let ptx = compile_ptx(kernels::KERNEL_SOURCE)
            .map_err(|e| CyaneaError::Other(format!("CUDA kernel compile: {e}")))?;

        let module = ctx
            .load_module(ptx)
            .map_err(|e| CyaneaError::Other(format!("CUDA module load: {e}")))?;

        let load = |name: &str| -> Result<CudaFunction> {
            module
                .load_function(name)
                .map_err(|e| CyaneaError::Other(format!("CUDA function '{name}': {e}")))
        };

        Ok(Self {
            ctx,
            stream,
            reduce_sum_fn: load("reduce_sum")?,
            reduce_min_fn: load("reduce_min")?,
            reduce_max_fn: load("reduce_max")?,
            pairwise_euclidean_fn: load("pairwise_euclidean")?,
            pairwise_manhattan_fn: load("pairwise_manhattan")?,
            pairwise_cosine_fn: load("pairwise_cosine")?,
            matmul_fn: load("matmul")?,
        })
    }

    /// Expose the CUDA context for reuse by other crates (e.g. cyanea-align).
    pub fn context(&self) -> &Arc<CudaContext> {
        &self.ctx
    }

    /// Runs a two-phase reduction.
    fn run_reduce(
        &self,
        buf: &Buffer,
        func: &CudaFunction,
        identity: f64,
        host_combine: fn(f64, f64) -> f64,
    ) -> Result<f64> {
        let host_data = buf
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".into()))?;
        if host_data.is_empty() {
            return Err(CyaneaError::InvalidInput("buffer is empty".into()));
        }

        let count = host_data.len() as u32;
        let block_size = 256u32;
        let num_blocks = (count + block_size - 1) / block_size;

        let input_dev: CudaSlice<f64> = self
            .stream
            .clone_htod(host_data)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;

        let output_dev: CudaSlice<f64> = self
            .stream
            .alloc_zeros(num_blocks as usize)
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc: {e}")))?;

        let cfg = LaunchConfig {
            grid_dim: (num_blocks, 1, 1),
            block_dim: (block_size, 1, 1),
            shared_mem_bytes: 0,
        };

        unsafe {
            self.stream
                .launch_builder(func)
                .arg(&input_dev)
                .arg(&output_dev)
                .arg(&count)
                .launch(cfg)
        }
        .map_err(|e| CyaneaError::Other(format!("CUDA launch: {e}")))?;

        let partials = self
            .stream
            .clone_dtoh(&output_dev)
            .map_err(|e| CyaneaError::Other(format!("CUDA dtoh: {e}")))?;

        let mut result = identity;
        for &p in &partials {
            result = host_combine(result, p);
        }
        Ok(result)
    }
}

impl Backend for CudaBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: "CUDA GPU".to_string(),
            kind: BackendKind::Cuda,
            total_memory: 0,
            max_parallelism: 1024,
        }
    }

    fn buffer_from_slice(&self, data: &[f64]) -> Result<Buffer> {
        let cuda_slice: CudaSlice<f64> = self
            .stream
            .clone_htod(data)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        Ok(Buffer::from_cuda(
            Some(data.to_vec()),
            cuda_slice,
            data.len(),
        ))
    }

    fn buffer_zeros(&self, len: usize) -> Result<Buffer> {
        let cuda_slice: CudaSlice<f64> = self
            .stream
            .alloc_zeros(len)
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc: {e}")))?;
        Ok(Buffer::from_cuda(Some(vec![0.0; len]), cuda_slice, len))
    }

    fn read_buffer(&self, buf: &Buffer) -> Result<Vec<f64>> {
        if let Some(ref data) = buf.host_data {
            return Ok(data.clone());
        }
        if let Some(ref cs) = buf.cuda_slice {
            return self
                .stream
                .clone_dtoh(cs)
                .map_err(|e| CyaneaError::Other(format!("CUDA dtoh: {e}")));
        }
        Err(CyaneaError::Other("buffer has no data".into()))
    }

    fn write_buffer(&self, buf: &mut Buffer, data: &[f64]) -> Result<()> {
        if data.len() != buf.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "write_buffer length mismatch: buffer has {} elements, got {}",
                buf.len(),
                data.len()
            )));
        }
        let cuda_slice: CudaSlice<f64> = self
            .stream
            .clone_htod(data)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        buf.host_data = Some(data.to_vec());
        buf.cuda_slice = Some(cuda_slice);
        Ok(())
    }

    fn buffer_len(&self, buf: &Buffer) -> usize {
        buf.len()
    }

    fn elementwise_map(
        &self,
        input: &Buffer,
        output: &mut Buffer,
        f: &dyn Fn(f64) -> f64,
    ) -> Result<()> {
        // CPU fallback — closures can't be sent to GPU.
        let src = input
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("input buffer has no host data".into()))?;
        if output.len() != input.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "elementwise_map: output length {} != input length {}",
                output.len(),
                input.len()
            )));
        }
        let mapped: Vec<f64> = src.iter().map(|&x| f(x)).collect();
        let cuda_slice: CudaSlice<f64> = self
            .stream
            .clone_htod(&mapped)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        output.host_data = Some(mapped);
        output.cuda_slice = Some(cuda_slice);
        Ok(())
    }

    fn reduce_sum(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.reduce_sum_fn, 0.0, |a, b| a + b)
    }

    fn reduce_min(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.reduce_min_fn, f64::INFINITY, f64::min)
    }

    fn reduce_max(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.reduce_max_fn, f64::NEG_INFINITY, f64::max)
    }

    fn pairwise_distance_matrix(
        &self,
        data: &Buffer,
        n: usize,
        dim: usize,
        metric: DistanceMetricGpu,
    ) -> Result<Buffer> {
        let flat = data
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".into()))?;
        if flat.len() != n * dim {
            return Err(CyaneaError::InvalidInput(format!(
                "pairwise_distance_matrix: expected {} elements ({}x{}), got {}",
                n * dim,
                n,
                dim,
                flat.len()
            )));
        }

        let func = match metric {
            DistanceMetricGpu::Euclidean => &self.pairwise_euclidean_fn,
            DistanceMetricGpu::Manhattan => &self.pairwise_manhattan_fn,
            DistanceMetricGpu::Cosine => &self.pairwise_cosine_fn,
        };

        let input_dev: CudaSlice<f64> = self
            .stream
            .clone_htod(flat)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        let result_dev: CudaSlice<f64> = self
            .stream
            .alloc_zeros(n * n)
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc: {e}")))?;

        let block_dim = 16u32;
        let grid_x = ((n as u32) + block_dim - 1) / block_dim;
        let grid_y = grid_x;
        let n_u32 = n as u32;
        let dim_u32 = dim as u32;

        let cfg = LaunchConfig {
            grid_dim: (grid_x, grid_y, 1),
            block_dim: (block_dim, block_dim, 1),
            shared_mem_bytes: 0,
        };

        unsafe {
            self.stream
                .launch_builder(func)
                .arg(&input_dev)
                .arg(&result_dev)
                .arg(&n_u32)
                .arg(&dim_u32)
                .launch(cfg)
        }
        .map_err(|e| CyaneaError::Other(format!("CUDA launch: {e}")))?;

        let result_data = self
            .stream
            .clone_dtoh(&result_dev)
            .map_err(|e| CyaneaError::Other(format!("CUDA dtoh: {e}")))?;

        Ok(Buffer::from_cuda(Some(result_data), result_dev, n * n))
    }

    fn matrix_multiply(
        &self,
        a: &Buffer,
        b: &Buffer,
        m: usize,
        k: usize,
        n: usize,
    ) -> Result<Buffer> {
        let a_data = a
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer a has no host data".into()))?;
        let b_data = b
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer b has no host data".into()))?;
        if a_data.len() != m * k {
            return Err(CyaneaError::InvalidInput(format!(
                "matrix_multiply: a has {} elements, expected {}x{}={}",
                a_data.len(),
                m,
                k,
                m * k
            )));
        }
        if b_data.len() != k * n {
            return Err(CyaneaError::InvalidInput(format!(
                "matrix_multiply: b has {} elements, expected {}x{}={}",
                b_data.len(),
                k,
                n,
                k * n
            )));
        }

        let a_dev: CudaSlice<f64> = self
            .stream
            .clone_htod(a_data)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        let b_dev: CudaSlice<f64> = self
            .stream
            .clone_htod(b_data)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        let c_dev: CudaSlice<f64> = self
            .stream
            .alloc_zeros(m * n)
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc: {e}")))?;

        let tile = 16u32;
        let grid_x = ((n as u32) + tile - 1) / tile;
        let grid_y = ((m as u32) + tile - 1) / tile;
        let m_u32 = m as u32;
        let k_u32 = k as u32;
        let n_u32 = n as u32;

        let cfg = LaunchConfig {
            grid_dim: (grid_x, grid_y, 1),
            block_dim: (tile, tile, 1),
            shared_mem_bytes: 0,
        };

        unsafe {
            self.stream
                .launch_builder(&self.matmul_fn)
                .arg(&a_dev)
                .arg(&b_dev)
                .arg(&c_dev)
                .arg(&m_u32)
                .arg(&k_u32)
                .arg(&n_u32)
                .launch(cfg)
        }
        .map_err(|e| CyaneaError::Other(format!("CUDA launch: {e}")))?;

        let result_data = self
            .stream
            .clone_dtoh(&c_dev)
            .map_err(|e| CyaneaError::Other(format!("CUDA dtoh: {e}")))?;

        Ok(Buffer::from_cuda(Some(result_data), c_dev, m * n))
    }

    fn batch_pairwise(
        &self,
        items: &Buffer,
        n: usize,
        item_len: usize,
        f: &dyn Fn(&[f64], &[f64]) -> f64,
    ) -> Result<Buffer> {
        // CPU fallback — closures can't be sent to GPU.
        let flat = items
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".into()))?;
        if flat.len() != n * item_len {
            return Err(CyaneaError::InvalidInput(format!(
                "batch_pairwise: expected {} elements ({}x{}), got {}",
                n * item_len,
                n,
                item_len,
                flat.len()
            )));
        }
        let mut result = vec![0.0_f64; n * n];
        for i in 0..n {
            let row_i = &flat[i * item_len..(i + 1) * item_len];
            for j in (i + 1)..n {
                let row_j = &flat[j * item_len..(j + 1) * item_len];
                let val = f(row_i, row_j);
                result[i * n + j] = val;
                result[j * n + i] = val;
            }
        }
        let cuda_slice: CudaSlice<f64> = self
            .stream
            .clone_htod(&result)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod: {e}")))?;
        Ok(Buffer::from_cuda(Some(result), cuda_slice, n * n))
    }
}

#[cfg(all(test, feature = "cuda"))]
mod tests {
    use super::*;

    fn backend() -> CudaBackend {
        CudaBackend::new().expect("CUDA should be available")
    }

    #[test]
    fn device_info_kind() {
        let info = backend().device_info();
        assert_eq!(info.kind, BackendKind::Cuda);
    }

    #[test]
    fn buffer_round_trip() {
        let b = backend();
        let data = vec![1.0, 2.0, 3.0];
        let buf = b.buffer_from_slice(&data).unwrap();
        assert_eq!(b.buffer_len(&buf), 3);
        assert_eq!(b.read_buffer(&buf).unwrap(), data);
    }

    #[test]
    fn reduce_sum_ok() {
        let b = backend();
        let buf = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let sum = b.reduce_sum(&buf).unwrap();
        assert!((sum - 10.0).abs() < 1e-10);
    }

    #[test]
    fn pairwise_euclidean() {
        let b = backend();
        let data = b
            .buffer_from_slice(&[0.0, 0.0, 3.0, 0.0, 0.0, 4.0])
            .unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 3, 2, DistanceMetricGpu::Euclidean)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 3 + 1] - 3.0).abs() < 1e-10);
        assert!((mat[0 * 3 + 2] - 4.0).abs() < 1e-10);
        assert!((mat[1 * 3 + 2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn matrix_multiply_known() {
        let b = backend();
        let a = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let mb = b.buffer_from_slice(&[5.0, 6.0, 7.0, 8.0]).unwrap();
        let result = b.matrix_multiply(&a, &mb, 2, 2, 2).unwrap();
        let data = b.read_buffer(&result).unwrap();
        assert!((data[0] - 19.0).abs() < 1e-10);
        assert!((data[1] - 22.0).abs() < 1e-10);
        assert!((data[2] - 43.0).abs() < 1e-10);
        assert!((data[3] - 50.0).abs() < 1e-10);
    }
}
