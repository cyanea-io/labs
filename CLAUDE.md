# Cyanea Labs

Rust bioinformatics ecosystem — 13 crates, 659+ tests, targeting native, WASM, and Python.

## Workspace

Cargo workspace at `Cargo.toml`. Rust 1.93+, edition 2021, MIT/Apache-2.0.

```
cargo check --workspace                    # Check all crates
cargo test --workspace --exclude cyanea-py # Test all (py needs PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1)
cargo test -p cyanea-align                 # Test one crate
cargo check -p cyanea-wasm --features wasm # Check with wasm-bindgen
```

## Crate Map

All crates are complete. Each has `docs/STATUS.md` with full API docs.

| Crate | Purpose | Tests | Key deps |
|-------|---------|------:|----------|
| **cyanea-core** | Traits, errors, SHA-256, zstd, mmap | 14 | thiserror, sha2, zstd, flate2, memmap2 |
| **cyanea-seq** | DNA/RNA/protein, FASTA/FASTQ, k-mers, quality | 44 | cyanea-core, needletail |
| **cyanea-io** | CSV, VCF, BED, GFF3 (feature-gated) | 26 | cyanea-core, cyanea-omics, csv |
| **cyanea-align** | NW, SW, semi-global, MSA, banded, GPU dispatch | 77 | cyanea-core |
| **cyanea-omics** | Genomic coords, intervals, matrices, variants, AnnData | 89 | cyanea-core |
| **cyanea-stats** | Descriptive, correlation, t-tests, distributions, PCA | 77 | cyanea-core |
| **cyanea-ml** | Clustering, distances, embeddings, KNN, PCA, t-SNE | 109 | cyanea-core |
| **cyanea-chem** | SMILES/SDF, fingerprints, properties, substructure | 35 | cyanea-core, sha2 |
| **cyanea-struct** | PDB, geometry, DSSP, Kabsch, contact maps | 37 | cyanea-core, sha2 |
| **cyanea-phylo** | Newick/NEXUS, distances, UPGMA/NJ, Fitch/Sankoff | 73 | cyanea-core, cyanea-ml (optional) |
| **cyanea-gpu** | Backend trait, CPU impl, buffers, ops | 45 | cyanea-core |
| **cyanea-wasm** | JSON-based WASM bindings, wasm-bindgen | 56 | cyanea-core/seq/io/align/stats/ml, serde_json |
| **cyanea-py** | Python bindings via PyO3 | — | cyanea-core/seq/io/align/stats/ml, pyo3 |

## Conventions

### Error handling
- All crates use `CyaneaError` from `cyanea-core` (thiserror 2.x)
- Functions return `Result<T, CyaneaError>` (aliased as `Result<T>`)
- Binding crates convert: WASM → JSON `{"error": "msg"}`, Python → `PyErr`

### Feature flags
- `default = ["std"]` — every domain crate
- `std` gates: `memmap2`, `zstd`, `flate2`, file I/O
- `wasm` — WASM-specific code (wasm-bindgen annotations in cyanea-wasm)
- `serde` — serialization support
- `cuda` / `metal` — GPU backends in cyanea-gpu and cyanea-align
- `ml` — optional cyanea-ml dependency in cyanea-phylo

### Code patterns
- Each alignment module (NW, SW, semi-global) has its own private `push_cigar` helper
- WASM functions return `String` via `wasm_ok(val)` / `wasm_err(msg)` / `wasm_result(r)` JSON envelope
- Python bindings use `IntoPyResult` trait for `CyaneaError → PyErr` conversion
- Workspace shared deps in root `Cargo.toml` `[workspace.dependencies]`

### Testing
- Unit tests inline (`#[cfg(test)]` modules at bottom of each file)
- Doc tests for key public API functions
- Test data: inline strings/vecs, no external fixture files
- `tempfile` crate for tests needing filesystem

## Dependency Graph

```
cyanea-core (foundation — no internal deps)
├── cyanea-seq (+ needletail)
├── cyanea-align
├── cyanea-omics
├── cyanea-stats
├── cyanea-ml
├── cyanea-chem (+ sha2)
├── cyanea-struct (+ sha2)
├── cyanea-phylo (+ optional cyanea-ml)
├── cyanea-gpu
├── cyanea-io (+ cyanea-omics, csv)
├── cyanea-wasm (+ cyanea-seq/io/align/stats/ml, serde_json, wasm-bindgen)
└── cyanea-py (+ cyanea-seq/io/align/stats/ml, pyo3)
```

## Building Bindings

### WASM
```bash
cargo check -p cyanea-wasm --features wasm   # Check with wasm-bindgen
# Full build: wasm-pack build cyanea-wasm --features wasm
```

### Python
```bash
# Needs Python 3.13 or 3.14 with compat flag
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py
# Full build: cd cyanea-py && PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

### NIF (Elixir)
The NIF crate lives in `../cyanea/native/cyanea_native/` and depends on labs crates via path.
```bash
cargo check -p cyanea-native   # Check only (can't link without BEAM)
# Full build: cd ../cyanea && mix compile
```

## What's Not Done

- **GPU backends**: CUDA and Metal stubs exist but need hardware SDKs to implement
- **SIMD**: Banded alignment uses scalar code; true SIMD vectorization deferred pending profiling
- **Publishing**: Not yet on crates.io, PyPI, or npm
- **Advanced formats**: SAM/BAM/CRAM, Parquet, HDF5/Zarr not implemented in cyanea-io
- **Advanced stats**: ANOVA, chi-squared, Fisher's exact, UMAP not implemented
