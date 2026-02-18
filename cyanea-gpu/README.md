# cyanea-gpu

> GPU compute abstraction layer with CPU, CUDA, Metal, and WebGPU backends.

## What's Inside

- **Backend trait** -- uniform API for compute operations across backends
- **CPU backend** -- always-available reference implementation
- **CUDA backend** -- NVIDIA GPU via cudarc + NVRTC runtime compilation, native f64
- **Metal backend** -- Apple GPU via metal-rs, MSL compute shaders, f32 with f64 boundary conversion
- **WebGPU backend** -- cross-platform via wgpu, WGSL shaders, works native and in-browser
- **Buffers** -- typed compute buffers backed by host memory or device memory
- **Operations** -- reduce (sum/min/max/mean), element-wise map, pairwise distance matrix, matrix multiply, batch z-score
- **Tiled pairwise distances** -- larger-than-memory matrices via tiling
- **GPU k-mer counting** -- 2-bit encoded parallel k-mer counting (k <= 14)
- **GPU Smith-Waterman** -- batch protein alignment with BLOSUM62, affine gaps
- **GPU MinHash** -- parallel bottom-k MinHash sketching with Jaccard estimation
- **Auto-selection** -- `auto_backend()` picks Metal > CUDA > WebGPU > CPU

## Quick Start

```toml
[dependencies]
cyanea-gpu = { version = "0.1", features = ["metal"] }  # or "cuda" or "wgpu"
```

```rust
use cyanea_gpu::{auto_backend, ops};

let backend = auto_backend();
println!("Using: {}", backend.device_info().name);

let buf = backend.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
let sum = ops::reduce_sum(&*backend, &buf).unwrap();
println!("Sum: {}", sum);
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `cuda` | No | NVIDIA CUDA backend (cudarc + NVRTC) |
| `metal` | No | Apple Metal backend (metal-rs) |
| `wgpu` | No | WebGPU backend (wgpu + bytemuck + pollster) |
| `serde` | No | Serialize/Deserialize derives |
| `parallel` | No | Rayon parallelism for CPU backend |

No default features -- the crate is lightweight with only the CPU backend.

## Modules

| Module | Description |
|--------|-------------|
| `backend` | `Backend` trait, `BackendKind`, `DeviceInfo`, `DistanceMetricGpu` |
| `buffer` | `Buffer` type (host-side or device-backed) |
| `cpu` | CPU reference backend |
| `metal` | Metal backend (feature-gated) |
| `cuda` | CUDA backend (feature-gated) |
| `wgpu_backend` | WebGPU backend (feature-gated) |
| `ops` | High-level operations (reduce, distance, matmul, z-score) |
| `kmer` | GPU k-mer counting |
| `smith_waterman` | GPU batch protein Smith-Waterman |
| `minhash_gpu` | GPU MinHash sketching |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
