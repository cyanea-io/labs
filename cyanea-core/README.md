# cyanea-core

> Shared primitives, traits, and utilities for the Cyanea bioinformatics ecosystem.

## What's Inside

- **Error handling** -- `CyaneaError` enum used by every crate in the workspace
- **Core traits** -- `Sequence`, `ContentAddressable`, `Compressible`, `Scored`, `Annotated`, `Summarizable`
- **SHA-256 hashing** -- in-memory and streaming file hashing
- **Compression** -- zstd and gzip compress/decompress with auto-detection
- **Memory mapping** -- read-only memory-mapped file access via `memmap2`
- **Probability types** -- `LogProb` and `PhredProb` for numerically stable log-space computation
- **Rank/select bitvectors** -- O(1) rank, O(log n) select with superblock index
- **Wavelet matrix** -- O(log sigma) access/rank/select over integer alphabets
- **Fenwick tree** -- O(log n) prefix sums and point updates

## Quick Start

```toml
[dependencies]
cyanea-core = "0.1"
```

```rust
use cyanea_core::{sha256, LogProb};

let hash = sha256(b"ACGTACGT");

let p = LogProb::from_prob(0.001).unwrap();
let q = LogProb::from_prob(0.002).unwrap();
let combined = p.ln_add(q); // log-sum-exp
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Compression (zstd, flate2), memory mapping (memmap2) |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |

## Modules

| Module | Description |
|--------|-------------|
| `error` | `CyaneaError` enum and `Result<T>` alias |
| `traits` | Core trait definitions for the ecosystem |
| `hash` | SHA-256 hashing (in-memory and file) |
| `compress` | Zstd/gzip compression with auto-detection |
| `mmap` | Memory-mapped file access (`std` only) |
| `prob` | `LogProb` and `PhredProb` newtypes |
| `bitvec` | `RankSelectBitVec` and `WaveletMatrix` |
| `fenwick` | Generic `FenwickTree<T>` |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Internal Architecture](docs/ARCHITECTURE.md)
- [Workspace Architecture](../docs/ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
