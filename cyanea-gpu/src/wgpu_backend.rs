//! WebGPU backend — GPU compute via wgpu for cross-platform and browser use.
//!
//! Uses WGSL compute shaders for reductions, pairwise distance matrices,
//! and matrix multiplication. All GPU computation uses `f32` (WebGPU has
//! no native f64 support); data is converted to/from `f64` at the
//! host-device boundary (same approach as the Metal backend).
//!
//! Feature-gated behind `wgpu`. Works on native platforms and in the browser
//! when compiled to `wasm32`.

use wgpu::util::DeviceExt;

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
use crate::buffer::Buffer;

/// WGSL shader sources.
const REDUCE_WGSL: &str = include_str!("shaders/reduce.wgsl");
const DISTANCE_WGSL: &str = include_str!("shaders/distance.wgsl");
const MATMUL_WGSL: &str = include_str!("shaders/matmul.wgsl");

/// WebGPU compute backend.
///
/// Compiles all WGSL shaders at construction time and caches pipeline
/// states. Uses `f32` on the GPU with f64↔f32 conversion at the boundary.
pub struct WgpuBackend {
    device: wgpu::Device,
    queue: wgpu::Queue,
    reduce_sum_pipeline: wgpu::ComputePipeline,
    reduce_min_pipeline: wgpu::ComputePipeline,
    reduce_max_pipeline: wgpu::ComputePipeline,
    pairwise_euclidean_pipeline: wgpu::ComputePipeline,
    pairwise_manhattan_pipeline: wgpu::ComputePipeline,
    pairwise_cosine_pipeline: wgpu::ComputePipeline,
    matmul_pipeline: wgpu::ComputePipeline,
}

impl WgpuBackend {
    /// Creates a WebGPU backend, requesting a device and compiling all shaders.
    ///
    /// # Errors
    ///
    /// Returns an error if no WebGPU adapter or device is available.
    pub fn new() -> Result<Self> {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        // Block on the async adapter/device request
        let (device, queue) = pollster::block_on(async {
            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    compatible_surface: None,
                    force_fallback_adapter: false,
                })
                .await
                .ok_or_else(|| CyaneaError::Other("no WebGPU adapter available".into()))?;

            adapter
                .request_device(
                    &wgpu::DeviceDescriptor {
                        label: Some("cyanea-gpu"),
                        required_features: wgpu::Features::empty(),
                        required_limits: wgpu::Limits::default(),
                        memory_hints: wgpu::MemoryHints::Performance,
                    },
                    None,
                )
                .await
                .map_err(|e| CyaneaError::Other(format!("WebGPU device request: {e}")))
        })?;

        // Compile shaders
        let reduce_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("reduce"),
            source: wgpu::ShaderSource::Wgsl(REDUCE_WGSL.into()),
        });
        let distance_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("distance"),
            source: wgpu::ShaderSource::Wgsl(DISTANCE_WGSL.into()),
        });
        let matmul_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("matmul"),
            source: wgpu::ShaderSource::Wgsl(MATMUL_WGSL.into()),
        });

        let make_pipeline =
            |module: &wgpu::ShaderModule, entry: &str| -> wgpu::ComputePipeline {
                device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some(entry),
                    layout: None,
                    module,
                    entry_point: Some(entry),
                    compilation_options: Default::default(),
                    cache: None,
                })
            };

        let reduce_sum_pipeline = make_pipeline(&reduce_module, "reduce_sum");
        let reduce_min_pipeline = make_pipeline(&reduce_module, "reduce_min");
        let reduce_max_pipeline = make_pipeline(&reduce_module, "reduce_max");
        let pairwise_euclidean_pipeline = make_pipeline(&distance_module, "pairwise_euclidean");
        let pairwise_manhattan_pipeline = make_pipeline(&distance_module, "pairwise_manhattan");
        let pairwise_cosine_pipeline = make_pipeline(&distance_module, "pairwise_cosine");
        let matmul_pipeline = make_pipeline(&matmul_module, "matmul");

        Ok(Self {
            device,
            queue,
            reduce_sum_pipeline,
            reduce_min_pipeline,
            reduce_max_pipeline,
            pairwise_euclidean_pipeline,
            pairwise_manhattan_pipeline,
            pairwise_cosine_pipeline,
            matmul_pipeline,
        })
    }

    /// Uploads f64 data as f32 to a GPU storage buffer.
    fn upload_f32(&self, data: &[f64]) -> wgpu::Buffer {
        let f32_data: Vec<f32> = data.iter().map(|&x| x as f32).collect();
        let bytes = bytemuck::cast_slice(&f32_data);
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: if bytes.is_empty() { &[0u8; 4] } else { bytes },
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_SRC
                    | wgpu::BufferUsages::COPY_DST,
            })
    }

    /// Creates a zero-initialized GPU storage buffer.
    fn alloc_f32(&self, count: usize) -> wgpu::Buffer {
        let size = ((count * 4) as u64).max(4);
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    /// Creates a uniform buffer from raw bytes.
    fn create_uniform(&self, data: &[u8]) -> wgpu::Buffer {
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: data,
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            })
    }

    /// Reads f32 data from a GPU buffer and converts to f64.
    fn read_f32_buffer(&self, buf: &wgpu::Buffer, count: usize) -> Vec<f64> {
        let size = (count * 4) as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: size.max(4),
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        encoder.copy_buffer_to_buffer(buf, 0, &staging, 0, size.max(4));
        self.queue.submit(Some(encoder.finish()));

        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = tx.send(result);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .unwrap_or(Err(wgpu::BufferAsyncError))
            .unwrap_or(());

        let data = slice.get_mapped_range();
        let f32_data: &[f32] = bytemuck::cast_slice(&data);
        f32_data[..count].iter().map(|&x| x as f64).collect()
    }

    /// Runs a two-phase reduction.
    fn run_reduce(
        &self,
        buf: &Buffer,
        pipeline: &wgpu::ComputePipeline,
        identity: f32,
        host_combine: fn(f32, f32) -> f32,
    ) -> Result<f64> {
        let host_data = buf
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".into()))?;
        if host_data.is_empty() {
            return Err(CyaneaError::InvalidInput("buffer is empty".into()));
        }

        let count = host_data.len();
        let num_groups = (count + 255) / 256;

        let input_buf = self.upload_f32(host_data);
        let output_buf = self.alloc_f32(num_groups);
        let params = (count as u32).to_le_bytes();
        let params_buf = self.create_uniform(&params);

        let bind_group_layout = pipeline.get_bind_group_layout(0);
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: input_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: output_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: params_buf.as_entire_binding(),
                },
            ],
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, Some(&bind_group), &[]);
            pass.dispatch_workgroups(num_groups as u32, 1, 1);
        }
        self.queue.submit(Some(encoder.finish()));

        let partials = self.read_f32_buffer(&output_buf, num_groups);
        let mut result = identity;
        for &p in &partials {
            result = host_combine(result, p as f32);
        }
        Ok(result as f64)
    }
}

// wgpu types are Send+Sync
unsafe impl Send for WgpuBackend {}
unsafe impl Sync for WgpuBackend {}

impl Backend for WgpuBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: "WebGPU".to_string(),
            kind: BackendKind::Wgpu,
            total_memory: 0,
            max_parallelism: 256,
        }
    }

    fn buffer_from_slice(&self, data: &[f64]) -> Result<Buffer> {
        let wgpu_buf = self.upload_f32(data);
        Ok(Buffer::from_wgpu(
            Some(data.to_vec()),
            wgpu_buf,
            data.len(),
        ))
    }

    fn buffer_zeros(&self, len: usize) -> Result<Buffer> {
        let wgpu_buf = self.alloc_f32(len);
        Ok(Buffer::from_wgpu(Some(vec![0.0; len]), wgpu_buf, len))
    }

    fn read_buffer(&self, buf: &Buffer) -> Result<Vec<f64>> {
        if let Some(ref data) = buf.host_data {
            return Ok(data.clone());
        }
        if let Some(ref wgpu_buf) = buf.wgpu_buffer {
            return Ok(self.read_f32_buffer(wgpu_buf, buf.len));
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
        buf.wgpu_buffer = Some(self.upload_f32(data));
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
        output.wgpu_buffer = Some(self.upload_f32(&mapped));
        output.host_data = Some(mapped);
        Ok(())
    }

    fn reduce_sum(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.reduce_sum_pipeline, 0.0, |a, b| a + b)
    }

    fn reduce_min(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(buf, &self.reduce_min_pipeline, f32::INFINITY, f32::min)
    }

    fn reduce_max(&self, buf: &Buffer) -> Result<f64> {
        self.run_reduce(
            buf,
            &self.reduce_max_pipeline,
            f32::NEG_INFINITY,
            f32::max,
        )
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
                "pairwise_distance_matrix: expected {} elements ({}×{}), got {}",
                n * dim,
                n,
                dim,
                flat.len()
            )));
        }

        let pipeline = match metric {
            DistanceMetricGpu::Euclidean => &self.pairwise_euclidean_pipeline,
            DistanceMetricGpu::Manhattan => &self.pairwise_manhattan_pipeline,
            DistanceMetricGpu::Cosine => &self.pairwise_cosine_pipeline,
        };

        let input_buf = self.upload_f32(flat);
        let result_buf = self.alloc_f32(n * n);

        // Params: n, dim as u32
        let mut params_bytes = Vec::with_capacity(8);
        params_bytes.extend_from_slice(&(n as u32).to_le_bytes());
        params_bytes.extend_from_slice(&(dim as u32).to_le_bytes());
        let params_buf = self.create_uniform(&params_bytes);

        let bind_group_layout = pipeline.get_bind_group_layout(0);
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: input_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: result_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: params_buf.as_entire_binding(),
                },
            ],
        });

        let wg_x = (n as u32 + 15) / 16;
        let wg_y = (n as u32 + 15) / 16;

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, Some(&bind_group), &[]);
            pass.dispatch_workgroups(wg_x, wg_y, 1);
        }
        self.queue.submit(Some(encoder.finish()));

        let result_data = self.read_f32_buffer(&result_buf, n * n);
        Ok(Buffer::from_wgpu(Some(result_data), result_buf, n * n))
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
                "matrix_multiply: a has {} elements, expected {}×{}={}",
                a_data.len(),
                m,
                k,
                m * k
            )));
        }
        if b_data.len() != k * n {
            return Err(CyaneaError::InvalidInput(format!(
                "matrix_multiply: b has {} elements, expected {}×{}={}",
                b_data.len(),
                k,
                n,
                k * n
            )));
        }

        let a_buf = self.upload_f32(a_data);
        let b_buf = self.upload_f32(b_data);
        let c_buf = self.alloc_f32(m * n);

        let mut params_bytes = Vec::with_capacity(12);
        params_bytes.extend_from_slice(&(m as u32).to_le_bytes());
        params_bytes.extend_from_slice(&(k as u32).to_le_bytes());
        params_bytes.extend_from_slice(&(n as u32).to_le_bytes());
        // Pad to 16-byte alignment for uniform buffer
        params_bytes.extend_from_slice(&0u32.to_le_bytes());
        let params_buf = self.create_uniform(&params_bytes);

        let bind_group_layout = self.matmul_pipeline.get_bind_group_layout(0);
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: a_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: b_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: c_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: params_buf.as_entire_binding(),
                },
            ],
        });

        let wg_x = (n as u32 + 15) / 16;
        let wg_y = (m as u32 + 15) / 16;

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.matmul_pipeline);
            pass.set_bind_group(0, Some(&bind_group), &[]);
            pass.dispatch_workgroups(wg_x, wg_y, 1);
        }
        self.queue.submit(Some(encoder.finish()));

        let result_data = self.read_f32_buffer(&c_buf, m * n);
        Ok(Buffer::from_wgpu(Some(result_data), c_buf, m * n))
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
                "batch_pairwise: expected {} elements ({}×{}), got {}",
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
        let wgpu_buf = self.upload_f32(&result);
        Ok(Buffer::from_wgpu(Some(result), wgpu_buf, n * n))
    }
}

#[cfg(all(test, feature = "wgpu"))]
mod tests {
    use super::*;

    fn try_backend() -> Option<WgpuBackend> {
        WgpuBackend::new().ok()
    }

    macro_rules! wgpu_test {
        ($name:ident, $body:expr) => {
            #[test]
            fn $name() {
                if let Some(b) = try_backend() {
                    $body(b);
                }
                // Skip if no WebGPU adapter available
            }
        };
    }

    wgpu_test!(device_info_kind, |b: WgpuBackend| {
        let info = b.device_info();
        assert_eq!(info.kind, BackendKind::Wgpu);
    });

    wgpu_test!(buffer_round_trip, |b: WgpuBackend| {
        let data = vec![1.0, 2.0, 3.0];
        let buf = b.buffer_from_slice(&data).unwrap();
        assert_eq!(b.buffer_len(&buf), 3);
        let read = b.read_buffer(&buf).unwrap();
        for (a, e) in read.iter().zip(&data) {
            assert!((a - e).abs() < 1e-5, "got {a}, expected {e}");
        }
    });

    wgpu_test!(buffer_zeros, |b: WgpuBackend| {
        let buf = b.buffer_zeros(4).unwrap();
        let data = b.read_buffer(&buf).unwrap();
        for &v in &data {
            assert!(v.abs() < 1e-10);
        }
    });

    wgpu_test!(write_buffer_ok, |b: WgpuBackend| {
        let mut buf = b.buffer_zeros(3).unwrap();
        b.write_buffer(&mut buf, &[7.0, 8.0, 9.0]).unwrap();
        let data = b.read_buffer(&buf).unwrap();
        assert!((data[0] - 7.0).abs() < 1e-5);
        assert!((data[1] - 8.0).abs() < 1e-5);
        assert!((data[2] - 9.0).abs() < 1e-5);
    });

    wgpu_test!(write_buffer_length_mismatch, |b: WgpuBackend| {
        let mut buf = b.buffer_zeros(3).unwrap();
        assert!(b.write_buffer(&mut buf, &[1.0, 2.0]).is_err());
    });

    wgpu_test!(reduce_sum_ok, |b: WgpuBackend| {
        let buf = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let sum = b.reduce_sum(&buf).unwrap();
        assert!((sum - 10.0).abs() < 1e-3, "got {sum}");
    });

    wgpu_test!(reduce_min_max, |b: WgpuBackend| {
        let buf = b.buffer_from_slice(&[3.0, 1.0, 4.0, 1.5]).unwrap();
        let min_val = b.reduce_min(&buf).unwrap();
        let max_val = b.reduce_max(&buf).unwrap();
        assert!((min_val - 1.0).abs() < 1e-5, "min: {min_val}");
        assert!((max_val - 4.0).abs() < 1e-5, "max: {max_val}");
    });

    wgpu_test!(pairwise_euclidean, |b: WgpuBackend| {
        let data = b
            .buffer_from_slice(&[0.0, 0.0, 3.0, 0.0, 0.0, 4.0])
            .unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 3, 2, DistanceMetricGpu::Euclidean)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 3 + 1] - 3.0).abs() < 1e-3);
        assert!((mat[0 * 3 + 2] - 4.0).abs() < 1e-3);
        assert!((mat[1 * 3 + 2] - 5.0).abs() < 1e-3);
    });

    wgpu_test!(matrix_multiply_known, |b: WgpuBackend| {
        let a = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let mb = b.buffer_from_slice(&[5.0, 6.0, 7.0, 8.0]).unwrap();
        let result = b.matrix_multiply(&a, &mb, 2, 2, 2).unwrap();
        let data = b.read_buffer(&result).unwrap();
        assert!((data[0] - 19.0).abs() < 1e-2, "C[0,0] = {}", data[0]);
        assert!((data[1] - 22.0).abs() < 1e-2, "C[0,1] = {}", data[1]);
        assert!((data[2] - 43.0).abs() < 1e-2, "C[1,0] = {}", data[2]);
        assert!((data[3] - 50.0).abs() < 1e-2, "C[1,1] = {}", data[3]);
    });

    wgpu_test!(batch_pairwise_abs_diff_sum, |b: WgpuBackend| {
        let items = b
            .buffer_from_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
            .unwrap();
        let result = b
            .batch_pairwise(&items, 3, 2, &|a, bv| {
                a.iter().zip(bv).map(|(x, y)| (x - y).abs()).sum()
            })
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 3 + 1] - 4.0).abs() < 1e-5);
    });
}
