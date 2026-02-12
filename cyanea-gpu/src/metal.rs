//! Metal backend — GPU compute via Apple Metal on macOS/Apple Silicon.
//!
//! Uses Metal Shading Language (MSL) compute shaders for reductions,
//! pairwise distance matrices, and matrix multiplication. All GPU
//! computation uses `float` (f32); data is converted to/from `f64` at
//! the host-device boundary since Apple Silicon GPUs lack native f64.

use std::ffi::c_void;

use metal_rs::{
    Buffer as MtlBuffer, CommandQueue, CompileOptions, ComputePipelineState, Device, MTLSize,
    MTLResourceOptions,
};

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
use crate::buffer::Buffer;
use crate::shaders;

/// Thread group size for reduction kernels.
const REDUCE_THREADS: u64 = 256;

/// Pre-compiled Metal pipeline states for all kernels.
struct MetalPipelines {
    reduce_sum: ComputePipelineState,
    reduce_min: ComputePipelineState,
    reduce_max: ComputePipelineState,
    pairwise_euclidean: ComputePipelineState,
    pairwise_manhattan: ComputePipelineState,
    pairwise_cosine: ComputePipelineState,
    matmul: ComputePipelineState,
}

/// Apple Metal compute backend.
///
/// Compiles all MSL shaders at construction time and caches pipeline
/// states. Buffers use `StorageModeShared` for unified memory access
/// on Apple Silicon.
pub struct MetalBackend {
    device: Device,
    queue: CommandQueue,
    pipelines: MetalPipelines,
}

impl MetalBackend {
    /// Creates a Metal backend, compiling all shaders.
    ///
    /// # Errors
    ///
    /// Returns an error if no Metal device is available or shader compilation fails.
    pub fn new() -> Result<Self> {
        let device = Device::system_default()
            .ok_or_else(|| CyaneaError::Other("no Metal device available".into()))?;
        let queue = device.new_command_queue();

        let opts = CompileOptions::new();

        let reduce_lib = device
            .new_library_with_source(shaders::REDUCE_MSL, &opts)
            .map_err(|e| CyaneaError::Other(format!("Metal reduce shader compile: {e}")))?;
        let distance_lib = device
            .new_library_with_source(shaders::DISTANCE_MSL, &opts)
            .map_err(|e| CyaneaError::Other(format!("Metal distance shader compile: {e}")))?;
        let matmul_lib = device
            .new_library_with_source(shaders::MATMUL_MSL, &opts)
            .map_err(|e| CyaneaError::Other(format!("Metal matmul shader compile: {e}")))?;

        let make_pipeline =
            |lib: &metal_rs::Library, name: &str| -> Result<ComputePipelineState> {
                let func = lib
                    .get_function(name, None)
                    .map_err(|e| CyaneaError::Other(format!("Metal function '{name}': {e}")))?;
                device
                    .new_compute_pipeline_state_with_function(&func)
                    .map_err(|e| CyaneaError::Other(format!("Metal pipeline '{name}': {e}")))
            };

        let pipelines = MetalPipelines {
            reduce_sum: make_pipeline(&reduce_lib, "reduce_sum")?,
            reduce_min: make_pipeline(&reduce_lib, "reduce_min")?,
            reduce_max: make_pipeline(&reduce_lib, "reduce_max")?,
            pairwise_euclidean: make_pipeline(&distance_lib, "pairwise_euclidean")?,
            pairwise_manhattan: make_pipeline(&distance_lib, "pairwise_manhattan")?,
            pairwise_cosine: make_pipeline(&distance_lib, "pairwise_cosine")?,
            matmul: make_pipeline(&matmul_lib, "matmul")?,
        };

        Ok(Self {
            device,
            queue,
            pipelines,
        })
    }

    /// Uploads f64 host data to a Metal buffer as f32.
    fn upload_f32(&self, data: &[f64]) -> MtlBuffer {
        let f32_data: Vec<f32> = data.iter().map(|&x| x as f32).collect();
        let byte_len = (f32_data.len() * std::mem::size_of::<f32>()) as u64;
        if byte_len == 0 {
            return self
                .device
                .new_buffer(4, MTLResourceOptions::StorageModeShared);
        }
        self.device.new_buffer_with_data(
            f32_data.as_ptr() as *const c_void,
            byte_len,
            MTLResourceOptions::StorageModeShared,
        )
    }

    /// Allocates a zero-initialized Metal buffer of `count` f32 elements.
    fn alloc_f32(&self, count: usize) -> MtlBuffer {
        let byte_len = ((count * std::mem::size_of::<f32>()) as u64).max(4);
        self.device
            .new_buffer(byte_len, MTLResourceOptions::StorageModeShared)
    }

    /// Creates a Metal buffer containing a single u32 value.
    fn scalar_u32(&self, val: u32) -> MtlBuffer {
        self.device.new_buffer_with_data(
            &val as *const u32 as *const c_void,
            std::mem::size_of::<u32>() as u64,
            MTLResourceOptions::StorageModeShared,
        )
    }

    /// Reads f32 data from a Metal buffer and converts to f64.
    fn read_f32(buf: &MtlBuffer, count: usize) -> Vec<f64> {
        let ptr = buf.contents() as *const f32;
        let slice = unsafe { std::slice::from_raw_parts(ptr, count) };
        slice.iter().map(|&x| x as f64).collect()
    }

    /// Reads f32 data from a Metal buffer as f32.
    fn read_f32_raw(buf: &MtlBuffer, count: usize) -> Vec<f32> {
        let ptr = buf.contents() as *const f32;
        unsafe { std::slice::from_raw_parts(ptr, count) }.to_vec()
    }

    /// Runs a two-phase reduction (GPU threadgroup reduce + CPU final sum).
    fn run_reduce(
        &self,
        buf: &Buffer,
        pipeline: &ComputePipelineState,
        identity: f32,
        host_combine: fn(f32, f32) -> f32,
    ) -> Result<f64> {
        let host_data = buf
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".into()))?;
        if host_data.is_empty() {
            return Err(CyaneaError::InvalidInput("buffer is empty".into()));
        }

        let count = host_data.len() as u32;
        let num_groups = ((count as u64) + REDUCE_THREADS - 1) / REDUCE_THREADS;

        let input_mtl = match &buf.metal_buffer {
            Some(mb) => mb.clone(),
            None => self.upload_f32(host_data),
        };

        let output_mtl = self.alloc_f32(num_groups as usize);
        let count_buf = self.scalar_u32(count);

        let cmd = self.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(pipeline);
        enc.set_buffer(0, Some(&input_mtl), 0);
        enc.set_buffer(1, Some(&output_mtl), 0);
        enc.set_buffer(2, Some(&count_buf), 0);

        let grid = MTLSize::new(num_groups * REDUCE_THREADS, 1, 1);
        let tg = MTLSize::new(REDUCE_THREADS, 1, 1);
        enc.dispatch_threads(grid, tg);
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();

        // CPU final reduction over partial results
        let partials = Self::read_f32_raw(&output_mtl, num_groups as usize);
        let mut result = identity;
        for &p in &partials {
            result = host_combine(result, p);
        }
        Ok(result as f64)
    }
}

// Safety: Metal device/queue/pipeline objects are thread-safe Obj-C wrappers.
unsafe impl Send for MetalBackend {}
unsafe impl Sync for MetalBackend {}

impl Backend for MetalBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: self.device.name().to_string(),
            kind: BackendKind::Metal,
            total_memory: self.device.recommended_max_working_set_size(),
            max_parallelism: self.pipelines.reduce_sum.max_total_threads_per_threadgroup()
                as usize,
        }
    }

    fn buffer_from_slice(&self, data: &[f64]) -> Result<Buffer> {
        let mtl_buf = self.upload_f32(data);
        Ok(Buffer::from_metal(Some(data.to_vec()), mtl_buf, data.len()))
    }

    fn buffer_zeros(&self, len: usize) -> Result<Buffer> {
        let mtl_buf = self.alloc_f32(len);
        Ok(Buffer::from_metal(Some(vec![0.0; len]), mtl_buf, len))
    }

    fn read_buffer(&self, buf: &Buffer) -> Result<Vec<f64>> {
        if let Some(ref data) = buf.host_data {
            return Ok(data.clone());
        }
        if let Some(ref mtl_buf) = buf.metal_buffer {
            return Ok(Self::read_f32(mtl_buf, buf.len));
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
        buf.host_data = Some(data.to_vec());
        buf.metal_buffer = Some(self.upload_f32(data));
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
        // CPU fallback — closures can't be sent to GPU shaders.
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
        output.metal_buffer = Some(self.upload_f32(&mapped));
        output.host_data = Some(mapped);
        Ok(())
    }

    fn reduce_sum(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.pipelines.reduce_sum, 0.0, |a, b| a + b)
    }

    fn reduce_min(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.pipelines.reduce_min, f32::INFINITY, f32::min)
    }

    fn reduce_max(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.pipelines.reduce_max, f32::NEG_INFINITY, f32::max)
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

        let pipeline = match metric {
            DistanceMetricGpu::Euclidean => &self.pipelines.pairwise_euclidean,
            DistanceMetricGpu::Manhattan => &self.pipelines.pairwise_manhattan,
            DistanceMetricGpu::Cosine => &self.pipelines.pairwise_cosine,
        };

        let input_mtl = match &data.metal_buffer {
            Some(mb) => mb.clone(),
            None => self.upload_f32(flat),
        };

        let result_mtl = self.alloc_f32(n * n);
        let n_buf = self.scalar_u32(n as u32);
        let dim_buf = self.scalar_u32(dim as u32);

        let cmd = self.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(pipeline);
        enc.set_buffer(0, Some(&input_mtl), 0);
        enc.set_buffer(1, Some(&result_mtl), 0);
        enc.set_buffer(2, Some(&n_buf), 0);
        enc.set_buffer(3, Some(&dim_buf), 0);

        let grid = MTLSize::new(n as u64, n as u64, 1);
        let tg = MTLSize::new(16.min(n as u64).max(1), 16.min(n as u64).max(1), 1);
        enc.dispatch_threads(grid, tg);
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();

        let result_data = Self::read_f32(&result_mtl, n * n);
        Ok(Buffer::from_metal(Some(result_data), result_mtl, n * n))
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

        let a_mtl = match &a.metal_buffer {
            Some(mb) => mb.clone(),
            None => self.upload_f32(a_data),
        };
        let b_mtl = match &b.metal_buffer {
            Some(mb) => mb.clone(),
            None => self.upload_f32(b_data),
        };

        let c_mtl = self.alloc_f32(m * n);
        let m_buf = self.scalar_u32(m as u32);
        let k_buf = self.scalar_u32(k as u32);
        let n_buf = self.scalar_u32(n as u32);

        let cmd = self.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&self.pipelines.matmul);
        enc.set_buffer(0, Some(&a_mtl), 0);
        enc.set_buffer(1, Some(&b_mtl), 0);
        enc.set_buffer(2, Some(&c_mtl), 0);
        enc.set_buffer(3, Some(&m_buf), 0);
        enc.set_buffer(4, Some(&k_buf), 0);
        enc.set_buffer(5, Some(&n_buf), 0);

        let grid = MTLSize::new(n as u64, m as u64, 1);
        let tg = MTLSize::new(16.min(n as u64).max(1), 16.min(m as u64).max(1), 1);
        enc.dispatch_threads(grid, tg);
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();

        let result_data = Self::read_f32(&c_mtl, m * n);
        Ok(Buffer::from_metal(Some(result_data), c_mtl, m * n))
    }

    fn batch_pairwise(
        &self,
        items: &Buffer,
        n: usize,
        item_len: usize,
        f: &dyn Fn(&[f64], &[f64]) -> f64,
    ) -> Result<Buffer> {
        // CPU fallback — closures can't be sent to GPU shaders.
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
        let mtl_buf = self.upload_f32(&result);
        Ok(Buffer::from_metal(Some(result), mtl_buf, n * n))
    }
}

#[cfg(all(test, feature = "metal"))]
mod tests {
    use super::*;

    fn backend() -> MetalBackend {
        MetalBackend::new().expect("Metal should be available on macOS")
    }

    #[test]
    fn device_info_kind() {
        let info = backend().device_info();
        assert_eq!(info.kind, BackendKind::Metal);
        assert!(info.max_parallelism >= 1);
        assert!(!info.name.is_empty());
    }

    #[test]
    fn buffer_round_trip() {
        let b = backend();
        let data = vec![1.0, 2.0, 3.0];
        let buf = b.buffer_from_slice(&data).unwrap();
        assert_eq!(b.buffer_len(&buf), 3);
        let read = b.read_buffer(&buf).unwrap();
        for (a, e) in read.iter().zip(&data) {
            assert!((a - e).abs() < 1e-6, "got {a}, expected {e}");
        }
    }

    #[test]
    fn buffer_zeros() {
        let b = backend();
        let buf = b.buffer_zeros(4).unwrap();
        let data = b.read_buffer(&buf).unwrap();
        for &v in &data {
            assert!(v.abs() < 1e-10);
        }
    }

    #[test]
    fn write_buffer_ok() {
        let b = backend();
        let mut buf = b.buffer_zeros(3).unwrap();
        b.write_buffer(&mut buf, &[7.0, 8.0, 9.0]).unwrap();
        let data = b.read_buffer(&buf).unwrap();
        assert!((data[0] - 7.0).abs() < 1e-6);
        assert!((data[1] - 8.0).abs() < 1e-6);
        assert!((data[2] - 9.0).abs() < 1e-6);
    }

    #[test]
    fn write_buffer_length_mismatch() {
        let b = backend();
        let mut buf = b.buffer_zeros(3).unwrap();
        assert!(b.write_buffer(&mut buf, &[1.0, 2.0]).is_err());
    }

    #[test]
    fn elementwise_map_double() {
        let b = backend();
        let input = b.buffer_from_slice(&[1.0, 2.0, 3.0]).unwrap();
        let mut output = b.buffer_zeros(3).unwrap();
        b.elementwise_map(&input, &mut output, &|x| x * 2.0)
            .unwrap();
        let data = b.read_buffer(&output).unwrap();
        assert!((data[0] - 2.0).abs() < 1e-6);
        assert!((data[1] - 4.0).abs() < 1e-6);
        assert!((data[2] - 6.0).abs() < 1e-6);
    }

    #[test]
    fn reduce_sum_ok() {
        let b = backend();
        let buf = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let sum = b.reduce_sum(&buf).unwrap();
        assert!((sum - 10.0).abs() < 1e-4, "got {sum}");
    }

    #[test]
    fn reduce_sum_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(b.reduce_sum(&buf).is_err());
    }

    #[test]
    fn reduce_min_max() {
        let b = backend();
        let buf = b.buffer_from_slice(&[3.0, 1.0, 4.0, 1.5]).unwrap();
        let min = b.reduce_min(&buf).unwrap();
        let max = b.reduce_max(&buf).unwrap();
        assert!((min - 1.0).abs() < 1e-6, "min: {min}");
        assert!((max - 4.0).abs() < 1e-6, "max: {max}");
    }

    #[test]
    fn reduce_min_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(b.reduce_min(&buf).is_err());
    }

    #[test]
    fn reduce_max_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(b.reduce_max(&buf).is_err());
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
        assert!((mat[0 * 3 + 1] - 3.0).abs() < 1e-4, "d(0,1) = {}", mat[1]);
        assert!((mat[0 * 3 + 2] - 4.0).abs() < 1e-4, "d(0,2) = {}", mat[2]);
        assert!((mat[1 * 3 + 2] - 5.0).abs() < 1e-4, "d(1,2) = {}", mat[5]);
        assert!(mat[0].abs() < 1e-6);
    }

    #[test]
    fn pairwise_manhattan() {
        let b = backend();
        let data = b.buffer_from_slice(&[0.0, 0.0, 3.0, 4.0]).unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 2, 2, DistanceMetricGpu::Manhattan)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[1] - 7.0).abs() < 1e-4, "got {}", mat[1]);
    }

    #[test]
    fn pairwise_cosine() {
        let b = backend();
        let data = b.buffer_from_slice(&[1.0, 0.0, 0.0, 1.0]).unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 2, 2, DistanceMetricGpu::Cosine)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[1] - 1.0).abs() < 1e-4, "got {}", mat[1]);
    }

    #[test]
    fn matrix_multiply_identity() {
        let b = backend();
        let eye = b.buffer_from_slice(&[1.0, 0.0, 0.0, 1.0]).unwrap();
        let mat = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let result = b.matrix_multiply(&eye, &mat, 2, 2, 2).unwrap();
        let data = b.read_buffer(&result).unwrap();
        assert!((data[0] - 1.0).abs() < 1e-4);
        assert!((data[1] - 2.0).abs() < 1e-4);
        assert!((data[2] - 3.0).abs() < 1e-4);
        assert!((data[3] - 4.0).abs() < 1e-4);
    }

    #[test]
    fn matrix_multiply_known() {
        let b = backend();
        let a = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let mb = b.buffer_from_slice(&[5.0, 6.0, 7.0, 8.0]).unwrap();
        let result = b.matrix_multiply(&a, &mb, 2, 2, 2).unwrap();
        let data = b.read_buffer(&result).unwrap();
        assert!((data[0] - 19.0).abs() < 1e-3, "C[0,0] = {}", data[0]);
        assert!((data[1] - 22.0).abs() < 1e-3, "C[0,1] = {}", data[1]);
        assert!((data[2] - 43.0).abs() < 1e-3, "C[1,0] = {}", data[2]);
        assert!((data[3] - 50.0).abs() < 1e-3, "C[1,1] = {}", data[3]);
    }

    #[test]
    fn batch_pairwise_abs_diff_sum() {
        let b = backend();
        let items = b
            .buffer_from_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
            .unwrap();
        let result = b
            .batch_pairwise(&items, 3, 2, &|a, bv| {
                a.iter().zip(bv).map(|(x, y)| (x - y).abs()).sum()
            })
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 3 + 1] - 4.0).abs() < 1e-6);
        assert!((mat[0 * 3 + 2] - 8.0).abs() < 1e-6);
        assert!((mat[1 * 3 + 0] - 4.0).abs() < 1e-6);
    }

    #[test]
    fn reduce_sum_large() {
        let b = backend();
        let data: Vec<f64> = (1..=1000).map(|i| i as f64).collect();
        let buf = b.buffer_from_slice(&data).unwrap();
        let sum = b.reduce_sum(&buf).unwrap();
        let expected = 500500.0;
        assert!(
            (sum - expected).abs() / expected < 1e-3,
            "got {sum}, expected {expected}"
        );
    }
}
