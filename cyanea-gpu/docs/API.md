# API Reference -- cyanea-gpu

GPU compute abstraction layer with a trait-based backend system. Provides a uniform API for compute operations across CPU, CUDA, Metal, and WebGPU backends, plus domain-specific GPU-accelerated bioinformatics operations.

## Public API

### Backend abstraction (`backend.rs`)

| Type | Description |
|------|-------------|
| `BackendKind` | Enum: `Cpu`, `Cuda`, `Metal`, `Wgpu` |
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
| `Buffer` | Typed compute buffer -- host-side `Vec<f64>` for CPU, Metal `Buffer` (f32) for Metal, `CudaSlice<f64>` for CUDA, `wgpu::Buffer` (f32) for WebGPU |

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
| `tiled_pairwise_distance(backend, data, n, dim, metric, tile_size) -> Result<Vec<f64>>` | Tiled pairwise distances for larger-than-memory matrices |

### GPU K-mer Counting (`kmer.rs`)

| Function/Type | Description |
|---------------|-------------|
| `KmerCountResult` | Struct with `counts: Vec<u32>`, `k: usize` |
| `gpu_kmer_count(sequences, k) -> Result<KmerCountResult>` | Auto-dispatch k-mer counting (Metal > CUDA > CPU) |
| `gpu_kmer_count_cpu(sequences, k) -> Result<KmerCountResult>` | CPU-only k-mer counting |

- 2-bit encoding (A=0, C=1, G=2, T=3), N positions skipped
- GPU limited to k <= 14 (4^14 = 268M entries)
- Metal shader: parallel per-position with atomic count increments
- CUDA kernel: same approach with `atomicAdd`

### GPU Smith-Waterman Protein (`smith_waterman.rs`)

| Function/Type | Description |
|---------------|-------------|
| `SwResult` | Struct with `score: i32`, `query_end: usize`, `target_end: usize` |
| `gpu_smith_waterman_batch(pairs, submat, gap_open, gap_extend) -> Result<Vec<SwResult>>` | Auto-dispatch batch SW (Metal > CUDA > CPU) |
| `gpu_sw_cpu(pairs, submat, gap_open, gap_extend) -> Result<Vec<SwResult>>` | CPU-only batch SW |
| `BLOSUM62_24` | Standard BLOSUM62 matrix extended to 24x24 |

- 24 amino acid alphabet (A..V + B, Z, X, *)
- Affine gap penalties (gap_open + gap_extend for opening)
- Substitution matrix passed as 24x24 flat array
- One GPU thread per sequence pair

### GPU MinHash Sketch (`minhash_gpu.rs`)

| Function/Type | Description |
|---------------|-------------|
| `MinHashSketch` | Struct with `hashes: Vec<u64>`, `k: usize` |
| `gpu_minhash(sequences, k, sketch_size) -> Result<Vec<MinHashSketch>>` | Auto-dispatch MinHash (Metal > CUDA > CPU) |
| `gpu_minhash_cpu(sequences, k, sketch_size) -> Result<Vec<MinHashSketch>>` | CPU-only MinHash |
| `gpu_minhash_jaccard(a, b) -> Result<f64>` | Jaccard similarity from two sketches |

- MurmurHash3 64-bit finalizer for k-mer hashing
- Bottom-sketch: retains `sketch_size` smallest hashes
- GPU parallelizes k-mer hashing, CPU selects bottom-k
- Jaccard estimation via merge-based intersection counting

### Backends

| Backend | Status | Description |
|---------|--------|-------------|
| `CpuBackend` | Complete | Full reference implementation, always available |
| `CudaBackend` | Complete | NVIDIA GPU via cudarc + NVRTC runtime compilation; uses f64 natively |
| `MetalBackend` | Complete | Apple Metal via metal-rs; MSL compute shaders, f32 on GPU with f64 conversion at boundary |
| `WgpuBackend` | Complete | WebGPU via wgpu; WGSL compute shaders, f32 on GPU with f64 conversion at boundary |

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
- All operations use `double` (f64) -- NVIDIA GPUs have native f64 support
- Tiled 16x16 shared memory matrix multiplication
- 256-thread block tree reductions
- `elementwise_map` and `batch_pairwise` use CPU fallback (closures can't be sent to GPU)
- Exposes `context()` for reuse by other crates (e.g. cyanea-align GPU dispatch)
- Requires NVIDIA GPU with CUDA driver at runtime

#### WgpuBackend

- Uses WGSL compute shaders equivalent to MSL and CUDA kernels
- Creates device/queue via pollster for synchronous initialization
- Uses `f32` on GPU (WebGPU has no native f64), with f64<->f32 conversion at host-device boundary
- Same two-phase reduction and tiled matmul patterns as Metal/CUDA
- Works on native platforms (Vulkan, Metal, DX12) and in browser (WebGPU)
- `elementwise_map` and `batch_pairwise` use CPU fallback
- Feature-gated behind `wgpu` with `wgpu`, `bytemuck`, and `pollster` dependencies

### Shaders (`shaders/`)

| File | Description |
|------|-------------|
| `shaders/mod.rs` | Shader source constants via `include_str!` for MSL and WGSL |
| `shaders/*.metal` | MSL compute shaders (reduce, distance, matmul, kmer, sw, minhash) |
| `shaders/*.wgsl` | WGSL compute shaders (reduce, distance, matmul) |

### Auto-selection

| Function | Description |
|----------|-------------|
| `auto_backend() -> Box<dyn Backend>` | Returns best available backend: Metal > CUDA > Wgpu > CPU |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `cuda` | No | NVIDIA CUDA GPU backend (cudarc + NVRTC) |
| `metal` | No | Apple Metal GPU backend (metal-rs) |
| `wgpu` | No | WebGPU GPU backend (wgpu + bytemuck + pollster) |
| `serde` | No | Serialization support for `BackendKind`, `DeviceInfo`, `DistanceMetricGpu` |
| `parallel` | No | Rayon parallelism for CPU backend distance/matmul |

Note: No default features -- the crate is lightweight by default with only the CPU backend.

## Dependencies

- `cyanea-core` -- error types, `Summarizable` trait
- `metal-rs` (feature = "metal") -- Apple Metal bindings
- `cudarc` (feature = "cuda") -- Safe CUDA driver API
- `wgpu` (feature = "wgpu") -- WebGPU compute
- `bytemuck` (feature = "wgpu") -- Safe transmutation for GPU buffer data
- `pollster` (feature = "wgpu") -- Synchronous async runtime for wgpu initialization
- `rayon` (feature = "parallel") -- CPU parallelism

## Tests

- 62 unit tests + 2 doc tests (default, CPU only)
- 80 total with `--features metal` (+18 Metal backend tests)
- +10 with `--features wgpu` (WebGPU backend tests, require GPU adapter)
- CUDA tests (`--features cuda`) require an NVIDIA GPU at runtime

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | ~100 | Module declarations, `auto_backend()`, re-exports |
| `backend.rs` | ~216 | `Backend` trait, `DeviceInfo`, `DistanceMetricGpu` |
| `buffer.rs` | ~260 | `Buffer` type with Metal/CUDA/wgpu device backing |
| `cpu.rs` | ~498 | Full CPU backend implementation |
| `ops.rs` | ~400 | High-level operation functions incl. tiled pairwise |
| `kmer.rs` | ~250 | GPU k-mer counting with Metal/CUDA/CPU dispatch |
| `smith_waterman.rs` | ~350 | GPU SW protein alignment with Metal/CUDA/CPU dispatch |
| `minhash_gpu.rs` | ~300 | GPU MinHash sketch computation with Metal/CUDA/CPU dispatch |
| `wgpu_backend.rs` | ~500 | Full WebGPU backend implementation |
| `metal.rs` | ~637 | Full Metal backend -- MSL shaders, buffer management, all ops |
| `cuda/mod.rs` | ~485 | Full CUDA backend -- cudarc integration, NVRTC compilation |
| `cuda/kernels.rs` | ~340 | CUDA C kernel source -- reduce, distance, matmul, kmer, sw, minhash |
| `shaders/mod.rs` | -- | Shader source constants (`include_str!`) for MSL and WGSL |
| `shaders/*.metal` | -- | MSL shaders (reduce, distance, matmul, kmer, sw, minhash) |
| `shaders/*.wgsl` | -- | WGSL shaders (reduce, distance, matmul) |
