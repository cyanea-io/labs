# cyanea-gpu

GPU compute abstraction layer with a trait-based backend system. Provides a uniform API for compute operations across CPU, CUDA, and Metal backends.

## Status: Complete

All three backends (CPU, CUDA, Metal) are fully implemented and tested. The CPU backend serves as the always-available fallback and reference implementation. CUDA and Metal backends are feature-gated and provide GPU-accelerated reductions, pairwise distance matrices, and matrix multiplication.

## Public API

### Backend abstraction (`backend.rs`)

| Type | Description |
|------|-------------|
| `BackendKind` | Enum: `Cpu`, `Cuda`, `Metal` |
| `DeviceInfo` | `name`, `kind`, `total_memory`, `max_parallelism` |
| `DistanceMetricGpu` | Enum: `Euclidean`, `Manhattan`, `Cosine` |
| `trait Backend: Send + Sync` | Abstract compute interface |

**Backend trait methods:**

| Method | Description |
|--------|-------------|
| `device_info() -> DeviceInfo` | Device information |
| `buffer_from_slice(data) -> Result<Buffer>` | Create buffer from data |
| `buffer_zeros(len) -> Result<Buffer>` | Create zero-initialized buffer |
| `read_buffer(buf) -> Result<Vec<f64>>` | Read buffer contents to host |
| `write_buffer(buf, data) -> Result<()>` | Overwrite buffer contents |
| `buffer_len(buf) -> usize` | Element count |
| `elementwise_map(input, output, f) -> Result<()>` | Apply function element-wise |
| `reduce_sum(buf) -> Result<f64>` | Sum reduction |
| `reduce_min(buf) -> Result<f64>` | Minimum |
| `reduce_max(buf) -> Result<f64>` | Maximum |
| `pairwise_distance_matrix(data, n, dim, metric) -> Result<Buffer>` | Pairwise distances |
| `matrix_multiply(a, b, m, k, n) -> Result<Buffer>` | Matrix multiplication |
| `batch_pairwise(items, n, item_len, f) -> Result<Buffer>` | Custom all-pairs |

### Buffers (`buffer.rs`)

| Type | Description |
|------|-------------|
| `Buffer` | Typed compute buffer — host-side `Vec<f64>` for CPU, Metal `Buffer` (f32) for Metal, `CudaSlice<f64>` for CUDA |

### Operations (`ops.rs`)

| Function | Description |
|----------|-------------|
| `reduce_sum(backend, buf) -> Result<f64>` | Sum reduction |
| `reduce_min(backend, buf) -> Result<f64>` | Minimum |
| `reduce_max(backend, buf) -> Result<f64>` | Maximum |
| `reduce_mean(backend, buf) -> Result<f64>` | Mean |
| `elementwise_map(backend, buf, fn) -> Result<Buffer>` | Apply function element-wise |
| `pairwise_distance_matrix(backend, data, n, dim, metric) -> Result<Buffer>` | Pairwise distances |
| `matrix_multiply(backend, a, b, m, k, n) -> Result<Buffer>` | Matrix multiplication |
| `batch_pairwise(backend, items, n, item_len, f) -> Result<Buffer>` | Custom all-pairs |
| `batch_z_score(backend, data, n_rows, n_cols) -> Result<Buffer>` | Column-wise z-score normalization |

### Backends

| Backend | Status | Description |
|---------|--------|-------------|
| `CpuBackend` | Complete | Full reference implementation, always available |
| `CudaBackend` | Complete | NVIDIA GPU via cudarc + NVRTC runtime compilation; uses f64 natively |
| `MetalBackend` | Complete | Apple Metal via metal-rs; MSL compute shaders, f32 on GPU with f64 conversion at boundary |

#### MetalBackend

- Compiles all MSL shaders at construction time and caches `ComputePipelineState` handles
- Uses `StorageModeShared` for unified memory on Apple Silicon
- Two-phase reduction: GPU threadgroup reduce + CPU final combine
- GPU kernels: `reduce_sum`, `reduce_min`, `reduce_max`, `pairwise_euclidean`, `pairwise_manhattan`, `pairwise_cosine`, `matmul`
- `elementwise_map` and `batch_pairwise` use CPU fallback (closures can't be sent to GPU shaders)
- Requires macOS with Apple Silicon or discrete AMD GPU

#### CudaBackend

- Compiles CUDA C kernels to PTX at construction via NVRTC
- Uses `cudarc` for safe CUDA driver API access (dynamically loads libcuda, no compile-time SDK dependency)
- All operations use `double` (f64) — NVIDIA GPUs have native f64 support
- Tiled 16x16 shared memory matrix multiplication
- 256-thread block tree reductions
- `elementwise_map` and `batch_pairwise` use CPU fallback (closures can't be sent to GPU)
- Exposes `context()` for reuse by other crates (e.g. cyanea-align GPU dispatch)
- Requires NVIDIA GPU with CUDA driver at runtime

### Auto-selection

| Function | Description |
|----------|-------------|
| `auto_backend() -> Box<dyn Backend>` | Returns best available backend: Metal > CUDA > CPU |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `cuda` | No | NVIDIA CUDA GPU backend (cudarc + NVRTC) |
| `metal` | No | Apple Metal GPU backend (metal-rs) |
| `serde` | No | Serialization support for `BackendKind`, `DeviceInfo`, `DistanceMetricGpu` |
| `parallel` | No | Rayon parallelism for CPU backend distance/matmul |

Note: No default features -- the crate is lightweight by default with only the CPU backend.

## Dependencies

- `cyanea-core` -- error types, `Summarizable` trait
- `metal-rs` (feature = "metal") -- Apple Metal bindings
- `cudarc` (feature = "cuda") -- Safe CUDA driver API
- `rayon` (feature = "parallel") -- CPU parallelism

## Tests

- 43 unit tests + 2 doc tests (default, CPU only)
- 61 total with `--features metal` (+18 Metal backend tests)
- CUDA tests (`--features cuda`) require an NVIDIA GPU at runtime

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 93 | Module declarations, `auto_backend()`, re-exports |
| `backend.rs` | 213 | `Backend` trait, `DeviceInfo`, `DistanceMetricGpu` |
| `buffer.rs` | 238 | `Buffer` type with Metal/CUDA device backing |
| `cpu.rs` | 498 | Full CPU backend implementation |
| `ops.rs` | 282 | High-level operation functions |
| `metal.rs` | 637 | Full Metal backend — MSL shaders, buffer management, all ops |
| `cuda/mod.rs` | 485 | Full CUDA backend — cudarc integration, NVRTC compilation |
| `cuda/kernels.rs` | 198 | CUDA C kernel source — reduce, pairwise distance, matmul |
| `shaders/reduce.metal` | — | MSL reduction kernels (sum, min, max) |
| `shaders/distance.metal` | — | MSL pairwise distance kernels (euclidean, manhattan, cosine) |
| `shaders/matmul.metal` | — | MSL matrix multiplication kernel |
| `shaders/mod.rs` | — | Shader source constants (`include_str!`) |
