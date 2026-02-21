# Architecture -- cyanea-gpu

Internal design documentation for the GPU compute abstraction crate.

## Backend Trait

The `Backend` trait provides a unified interface for compute operations across all hardware:

```
trait Backend: Send + Sync {
    fn device_info() -> DeviceInfo;
    fn buffer_from_slice(data: &[f64]) -> Result<Buffer>;
    fn buffer_zeros(len: usize) -> Result<Buffer>;
    fn read_buffer(buf: &Buffer) -> Result<Vec<f64>>;
    fn write_buffer(buf: &Buffer, data: &[f64]) -> Result<()>;
    fn buffer_len(buf: &Buffer) -> usize;
    fn reduce_sum/min/max(buf: &Buffer) -> Result<f64>;
    fn elementwise_map(input, output, f) -> Result<()>;
    fn pairwise_distance_matrix(data, n, dim, metric) -> Result<Buffer>;
    fn matrix_multiply(a, b, m, k, n) -> Result<Buffer>;
    fn batch_pairwise(items, n, item_len, f) -> Result<Buffer>;
}
```

Each backend implements the full trait. The `Send + Sync` bound allows backends to be shared across threads and used with `Box<dyn Backend>`.

## Buffer

`Buffer` is an enum with variants per backend:

- **Cpu**: wraps `Vec<f64>` directly in memory
- **Metal**: wraps `metal::Buffer` (f32 on device), with f64<->f32 conversion at the host-device boundary
- **Cuda**: wraps `CudaSlice<f64>` on the GPU (native f64)
- **Wgpu**: wraps `wgpu::Buffer` (f32 on device), same conversion pattern as Metal

This design avoids generics in the trait by using a single concrete Buffer type with runtime dispatch.

## CPU Backend (`cpu.rs`)

Full reference implementation used for testing and as an always-available fallback:

- Reductions: single-pass iteration over `Vec<f64>`
- Pairwise distances: nested loops over O(n^2) pairs, Euclidean/Manhattan/Cosine
- Matrix multiply: naive O(m*k*n) triple loop (with optional Rayon parallelism via `parallel` feature)
- `elementwise_map` and `batch_pairwise`: direct application of Rust closures
- When `parallel` feature is enabled, distance matrix and matmul use `rayon::par_iter` for multi-core speedup

## CUDA Backend (`cuda/`)

NVIDIA GPU compute via `cudarc` safe bindings:

**Initialization**: dynamically loads `libcuda.so`/`cuda.dll` at runtime via cudarc, no compile-time CUDA SDK required. Creates a context on device 0.

**Kernel compilation**: CUDA C source strings in `cuda/kernels.rs` are compiled to PTX at construction time using NVRTC (runtime compilation). Pipeline states are cached as `CudaFunction` handles.

**Memory**: `CudaSlice<f64>` for device-side buffers. Data transfer via `htod_copy` (host-to-device) and `dtoh_sync_copy` (device-to-host). All operations use `double` (f64) since NVIDIA GPUs have native f64 support.

**Kernels**:
- `reduce_sum/min/max`: 256-thread blocks with shared memory tree reduction, final combine on CPU
- `pairwise_euclidean/manhattan/cosine`: one thread per (i,j) pair in the output matrix
- `matmul`: 16x16 tiled multiplication with shared memory for coalesced access
- `kmer_count`: one thread per position, atomic increment into global count array
- `smith_waterman`: one thread per sequence pair, anti-diagonal wavefront DP
- `minhash`: one thread per k-mer position, MurmurHash3 finalization

**Context sharing**: the `context()` method exposes the `CudaDevice` for reuse by other crates (e.g., cyanea-align GPU dispatch).

## Metal Backend (`metal.rs`)

Apple GPU compute via `metal-rs`:

**Initialization**: creates a `Device` (system default), `CommandQueue`, compiles all MSL shader sources from `shaders/*.metal` into `ComputePipelineState` handles at construction time.

**Memory model**: uses `StorageModeShared` for unified memory on Apple Silicon. Both CPU and GPU can access the same physical memory without explicit transfers. Data marshaled as f32 (Metal has limited f64 support) with conversion at boundary.

**Dispatch pattern**:
1. Create command buffer from the queue
2. Create compute command encoder
3. Set pipeline state and bind buffers
4. Dispatch threadgroups (threads_per_group typically 256)
5. End encoding, commit command buffer, wait for completion
6. Read results from shared buffer

**Two-phase reduction**: GPU performs per-threadgroup reduction into a partial results buffer, CPU combines the partial results. This avoids multiple kernel launches for large inputs.

**Closures**: `elementwise_map` and `batch_pairwise` cannot run on GPU (closures are CPU-only), so they fall back to CPU execution.

## WebGPU Backend (`wgpu_backend.rs`)

Cross-platform GPU compute via the `wgpu` crate:

**Initialization**: uses `pollster::block_on` for synchronous adapter/device creation. Requests `wgpu::Features::empty()` for maximum compatibility. Works across Vulkan (Linux/Windows), Metal (macOS), DX12 (Windows), and WebGPU (browser).

**Shaders**: WGSL compute shaders in `shaders/*.wgsl`, functionally equivalent to the MSL and CUDA kernels. Compiled to `ComputePipeline` handles at construction.

**Memory**: device buffers created with `MAP_READ | COPY_DST | STORAGE` usage flags. Data marshaled as f32 (WebGPU standard lacks f64). Host-side reads use `buffer.slice(..).map_async()` with `pollster::block_on`.

**Dispatch**: bind groups, dispatch workgroups, submit command encoder. Same logical structure as Metal dispatch.

## Auto-Selection

`auto_backend()` returns the first available backend at runtime:

1. **Metal** (if `metal` feature enabled): `MetalBackend::new()` succeeds on macOS with compatible GPU
2. **CUDA** (if `cuda` feature enabled): `CudaBackend::new()` succeeds when NVIDIA driver is loaded
3. **WebGPU** (if `wgpu` feature enabled): `WgpuBackend::new()` succeeds when any GPU adapter is available
4. **CPU** (always available): `CpuBackend::new()` as guaranteed fallback

Selection is compile-time feature-gated first (only enabled backends are tried), then runtime-checked (constructor may fail if hardware is unavailable).

## Tiling for Large Matrices

The `tiled_pairwise_distance` function handles distance matrices that exceed available GPU memory:

1. Divide the n*n output matrix into tile_size * tile_size blocks
2. For each block (i_start..i_end, j_start..j_end):
   - Upload the relevant rows of input data
   - Compute the tile on GPU
   - Download the tile to host memory
3. Assemble the full distance matrix on the host

This allows processing datasets with millions of points where the n*n output would exceed GPU memory.

## Domain-Specific GPU Modules

These modules are standalone (not part of the Backend trait) and handle their own backend dispatch:

**K-mer counting** (`kmer.rs`): 2-bit encodes DNA sequences, dispatches to Metal/CUDA/CPU. GPU kernel processes one position per thread, using atomic increments into a count array of size 4^k.

**Smith-Waterman** (`smith_waterman.rs`): batch protein alignment with affine gap penalties. Substitution matrix (BLOSUM62) uploaded as a GPU buffer. One thread computes the full DP matrix for one (query, target) pair. Anti-diagonal parallelism within each pair on GPU.

**MinHash** (`minhash_gpu.rs`): parallel k-mer hashing using MurmurHash3 64-bit finalizer. GPU computes all k-mer hashes in parallel, then CPU selects the bottom-k smallest hashes per sequence. Jaccard similarity estimated from sketch intersection.

## Module Dependencies

```
backend  <--  cpu, metal, cuda, wgpu_backend, ops
buffer   <--  cpu, metal, cuda, wgpu_backend, ops, kmer, smith_waterman, minhash_gpu
cpu      <--  kmer, smith_waterman, minhash_gpu (as fallback)
shaders  <--  metal, wgpu_backend
```

All modules depend on `cyanea-core` for `CyaneaError` and `Result`.
