# cyanea-gpu

GPU compute abstraction layer with a trait-based backend system. Provides a uniform API for compute operations across CPU, CUDA, and Metal backends.

## Status: Partial -- CPU Backend Only

The `Backend` trait abstraction and CPU reference implementation are complete and tested. CUDA and Metal backends are feature-gated stubs that panic on construction. The CPU backend serves as the always-available fallback and as the reference for future GPU implementations.

## Public API

### Backend abstraction (`backend.rs`)

| Type | Description |
|------|-------------|
| `BackendKind` | Enum: `Cpu`, `Cuda`, `Metal` |
| `DeviceInfo` | `name`, `kind`, `compute_capability`, `memory_mb` |
| `DistanceMetricGpu` | Enum: `Euclidean`, `Manhattan`, `Cosine` |
| `trait Backend: Send + Sync` | Abstract compute interface |

**Backend trait methods:**

| Method | Description |
|--------|-------------|
| `device_info() -> DeviceInfo` | Device information |
| `buffer_from_slice(data) -> Result<Buffer>` | Create buffer from data |
| Reduction, elementwise, matrix, distance ops | See Operations section |

### Buffers (`buffer.rs`)

| Type | Description |
|------|-------------|
| `Buffer` | Typed compute buffer (host-side for CPU, device-side for GPU) |

### Operations (`ops.rs`)

| Function | Description |
|----------|-------------|
| `reduce_sum(backend, buf) -> Result<f64>` | Sum reduction |
| `reduce_min(backend, buf) -> Result<f64>` | Minimum |
| `reduce_max(backend, buf) -> Result<f64>` | Maximum |
| `reduce_mean(backend, buf) -> Result<f64>` | Mean |
| `elementwise_map(backend, buf, fn) -> Result<Buffer>` | Apply function element-wise |
| `pairwise_distance_matrix(backend, data, metric) -> Result<Buffer>` | Pairwise distances |
| `matrix_multiply(backend, a, b) -> Result<Buffer>` | Matrix multiplication |
| `batch_pairwise(backend, queries, targets, metric) -> Result<Vec<Buffer>>` | Batch distances |
| `batch_z_score(backend, data) -> Result<Buffer>` | Batch z-score normalization |

### Backends

| Backend | Status | Description |
|---------|--------|-------------|
| `CpuBackend` | Complete | Full reference implementation, always available |
| `CudaBackend` | Stub | Panics on construction; gated behind `cuda` feature |
| `MetalBackend` | Stub | Panics on construction; gated behind `metal` feature |

### Auto-selection

| Function | Description |
|----------|-------------|
| `auto_backend() -> Box<dyn Backend>` | Returns best available backend (currently always CPU) |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `cuda` | No | CUDA GPU backend (future) |
| `metal` | No | Metal GPU backend (future) |
| `serde` | No | Serialization support |

Note: No default features -- the crate is lightweight by default with only the CPU backend.

## Dependencies

- `cyanea-core` -- error types

## Tests

43 unit tests + 2 doc tests covering the CPU backend, reduction operations, elementwise operations, distance matrices, and matrix multiplication.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 80 | Module declarations, `auto_backend()` |
| `backend.rs` | 212 | `Backend` trait, `DeviceInfo`, `DistanceMetricGpu` |
| `buffer.rs` | 146 | `Buffer` type |
| `cpu.rs` | 447 | Full CPU backend implementation |
| `ops.rs` | 281 | High-level operation functions |
| `cuda.rs` | 118 | CUDA stub |
| `metal.rs` | 117 | Metal stub |
