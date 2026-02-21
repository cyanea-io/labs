# Building Cyanea Labs

## Prerequisites

- **Rust 1.93+** (edition 2021)
- Cargo workspace at the repository root

Optional system dependencies for specific features:

| Feature | System requirement |
|---------|-------------------|
| `h5ad` | HDF5 1.10.x (`brew install hdf5@1.10` / `apt install libhdf5-dev`) |
| `cuda` | CUDA toolkit with `nvcc` in PATH |
| `metal` | macOS with Xcode Command Line Tools |
| `wgpu` | Vulkan, Metal, or DX12 runtime |
| WASM | `wasm-pack` (`cargo install wasm-pack`) |
| Python | Python 3.13+, `maturin` (`pip install maturin`) |
| NIF | Erlang/OTP + Elixir + Rustler |

## Native Rust

### Check all crates

```bash
cargo check --workspace
```

### Run tests

```bash
# All crates except Python (needs special env var)
cargo test --workspace --exclude cyanea-py

# Single crate
cargo test -p cyanea-core
cargo test -p cyanea-seq
cargo test -p cyanea-align

# With specific features
cargo test -p cyanea-align --features simd,wfa
cargo test -p cyanea-io --features vcf,bed,gff,sam,bam
cargo test -p cyanea-omics --features zarr
```

### Run benchmarks

```bash
cargo bench -p cyanea-align
cargo bench -p cyanea-seq
cargo bench -p cyanea-ml
cargo bench -p cyanea-stats
cargo bench -p cyanea-chem
cargo bench -p cyanea-struct
cargo bench -p cyanea-gpu
```

### Build documentation

```bash
cargo doc --workspace --no-deps --open
```

### Clippy and formatting

```bash
cargo clippy --workspace -- -D warnings
cargo fmt --all -- --check
```

## WASM

The `cyanea-wasm` crate wraps 10 domain modules into a JSON-based WebAssembly API.

### Check

```bash
cargo check -p cyanea-wasm --features wasm
```

### Build with wasm-pack

```bash
wasm-pack build cyanea-wasm --features wasm
```

### TypeScript wrapper

```bash
cd cyanea-wasm
npm run build:ts
```

Published as `@cyanea/bio` on npm. TypeScript types are in `ts/types.ts` with a wrapper in `ts/index.ts`.

### Usage in JavaScript

```javascript
import { parseNewick, alignDna, describe } from '@cyanea/bio';

const tree = parseNewick('((A:0.1,B:0.2):0.3,C:0.4);');
const alignment = alignDna('ACGT', 'ACGTT', 'local');
const stats = describe([1.0, 2.0, 3.0, 4.0]);
```

## Python (PyO3 / maturin)

The `cyanea-py` crate provides Python bindings for 11 submodules via PyO3.

### Check

```bash
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py
```

### Build and install

```bash
cd cyanea-py
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

### With NumPy interop

```bash
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py --features numpy
cd cyanea-py
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release --features numpy
```

### Usage in Python

```python
from cyanea import seq, align, stats, ml

dna = seq.DnaSequence(b"ATCGATCG")
print(dna.gc_content())

result = align.align_dna("ACGT", "ACGTT", mode="local")
print(result.score, result.cigar_string)

summary = stats.describe([1.0, 2.0, 3.0, 4.0])
print(summary.mean, summary.median)
```

### Python ABI compatibility

PyO3 0.23 officially supports up to Python 3.13. For Python 3.14+, the `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1` environment variable is required.

## GPU -- CUDA

CUDA support is available in `cyanea-gpu` (compute primitives) and `cyanea-align` (batch alignment).

### Prerequisites

- NVIDIA GPU with CUDA driver
- CUDA toolkit (cudarc dynamically loads `libcuda`, no compile-time SDK dependency for `cyanea-gpu`)
- For `cyanea-align` GPU kernels: `nvcc` in PATH (NVRTC runtime compilation)

### Build and test

```bash
# GPU compute primitives
cargo test -p cyanea-gpu --features cuda

# GPU-accelerated alignment
cargo test -p cyanea-align --features cuda

# Benchmarks
cargo bench -p cyanea-gpu --features cuda
```

### Details

- Uses `cudarc` for safe CUDA driver API access
- Kernels compiled to PTX at runtime via NVRTC
- Native `f64` precision on NVIDIA GPUs
- Operations: pairwise distances, matrix multiply, reductions, k-mer counting, Smith-Waterman, MinHash

## GPU -- Metal

Metal support is available in `cyanea-gpu` and `cyanea-align`.

### Prerequisites

- macOS with Apple Silicon or discrete AMD GPU
- Xcode Command Line Tools

### Build and test

```bash
# GPU compute primitives
cargo test -p cyanea-gpu --features metal

# GPU-accelerated alignment
cargo test -p cyanea-align --features metal

# Benchmarks
cargo bench -p cyanea-gpu --features metal
```

### Details

- Uses `metal-rs` for Apple Metal bindings
- MSL compute shaders compiled at construction time
- `StorageModeShared` buffers for unified memory on Apple Silicon
- `f32` on GPU with `f64` conversion at host boundary

## GPU -- WebGPU

Cross-platform GPU via `wgpu`.

### Prerequisites

- A Vulkan, Metal, or DX12 capable GPU at runtime

### Build and test

```bash
cargo test -p cyanea-gpu --features wgpu
```

### Details

- WGSL compute shaders
- Same operations as Metal/CUDA backends
- Works on native platforms and in browsers (WebGPU API)
- `f32` on GPU with `f64` conversion at host boundary

## Elixir NIF

The NIF crate lives in `../cyanea/native/cyanea_native/` and depends on labs crates via path.

### Check

```bash
cargo check -p cyanea-native
```

### Build via Mix

```bash
cd ../cyanea
mix compile
```

### Prerequisites

- Erlang/OTP
- Elixir
- Rustler (added as a Mix dependency)

## Feature Combinations Reference

Common `--features` flags for specific use cases:

| Use case | Command |
|----------|---------|
| Full IO stack | `cargo test -p cyanea-io --features csv,vcf,bed,gff,gtf,sam,bam,blast,maf,genbank,bigwig` |
| IO + CRAM | `cargo test -p cyanea-io --features sam,bam,cram` |
| IO + Parquet | `cargo test -p cyanea-io --features parquet` |
| IO + variant calling | `cargo test -p cyanea-io --features variant-calling` |
| IO + Stockholm/Clustal/Phylip | `cargo test -p cyanea-io --features stockholm,clustal,phylip` |
| IO + indexed BAM | `cargo test -p cyanea-io --features indexed-bam` |
| IO + indexed VCF | `cargo test -p cyanea-io --features indexed-vcf` |
| Alignment + SIMD | `cargo test -p cyanea-align --features simd` |
| Alignment + WFA | `cargo test -p cyanea-align --features wfa` |
| Alignment + GPU (macOS) | `cargo test -p cyanea-align --features metal` |
| Alignment + GPU (NVIDIA) | `cargo test -p cyanea-align --features cuda` |
| Omics + h5ad | `HDF5_DIR="$(brew --prefix hdf5@1.10)" cargo test -p cyanea-omics --features h5ad` |
| Omics + Zarr | `cargo test -p cyanea-omics --features zarr` |
| Omics + single-cell | `cargo test -p cyanea-omics --features single-cell` |
| Phylo + tree construction | `cargo test -p cyanea-phylo --features ml` |
| Stats + BLAS PCA | `cargo test -p cyanea-stats --features blas` |
| ML + BLAS PCA | `cargo test -p cyanea-ml --features blas` |
| Parallel everything | `cargo test -p cyanea-align --features parallel` |
| WASM build | `cargo check -p cyanea-wasm --features wasm` |
| Python build | `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py` |
| Python + NumPy | `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py --features numpy` |
| All GPU backends | `cargo test -p cyanea-gpu --features cuda,metal,wgpu` |
| Seq + MinHash | `cargo test -p cyanea-seq --features minhash` |

## Troubleshooting

### HDF5 not found (`h5ad` feature)

```
error: could not find HDF5 library
```

Install HDF5 1.10.x and set the path:

```bash
# macOS
brew install hdf5@1.10
export HDF5_DIR="$(brew --prefix hdf5@1.10)"

# Linux
sudo apt install libhdf5-dev
```

Note: `hdf5-sys` 0.8 does not support HDF5 2.0.

### Python ABI mismatch

```
error: the configured Python interpreter version (3.14) is newer than the maximum supported
```

Set the compatibility flag:

```bash
export PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1
```

### CUDA not found

```
error: could not find CUDA driver
```

Ensure `libcuda.so` (Linux) or `libcuda.dylib` (macOS) is in the library path. `cudarc` uses dynamic loading -- no compile-time SDK is needed for `cyanea-gpu`, but NVRTC (`libnvrtc`) is needed for runtime kernel compilation.

### Metal shader compilation fails

Metal requires macOS. The `metal` feature will fail to compile on Linux/Windows.

### wgpu: no suitable adapter found

```
error: no suitable GPU adapter found
```

Ensure a Vulkan, Metal, or DX12 GPU driver is available. On headless servers, try setting `WGPU_BACKEND=vulkan`.

### CRAM reference path

CRAM files may require an external reference FASTA:

```rust
use cyanea_io::cram::{CramConfig, parse_cram};
let config = CramConfig { reference_path: Some("ref.fa".into()) };
let records = parse_cram("reads.cram", &config)?;
```

Use `parse_cram_default()` for CRAM files that don't need an external reference.

### Noodles version alignment

The CRAM feature uses `noodles-cram` 0.72, which pins specific versions of `noodles-sam` (0.67), `noodles-fasta` (0.46), and `noodles-core` (0.15). These are managed in the workspace `[workspace.dependencies]` section -- do not override them in individual crates.
