# Cyanea Labs

Rust bioinformatics ecosystem — 13 crates, 1440+ tests, targeting native, WASM, Python, and Elixir NIFs.

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
| **cyanea-core** | Traits, errors, SHA-256, zstd, mmap, log-space probability types, rank/select bitvectors, wavelet matrix, Fenwick tree | 58 | thiserror, sha2, zstd, flate2, memmap2 |
| **cyanea-seq** | DNA/RNA/protein, FASTA/FASTQ, k-mers, 2-bit encoding, suffix array, FM-index, BWT, MinHash, pattern matching, PSSM/motif scanning, ORF finder, FASTA indexed reader (.fai), FMD-Index, quality trimming/filtering, codon tables (7 NCBI), codon usage/CAI, DUST/SEG/tandem repeat masking | 300 | cyanea-core, needletail |
| **cyanea-io** | CSV, VCF, BED, BEDPE, GFF3, SAM, BAM, CRAM, Parquet (feature-gated) | 81 | cyanea-core, cyanea-omics, csv, flate2, noodles, arrow/parquet |
| **cyanea-align** | NW, SW, semi-global, MSA, banded, seed-and-extend, minimizers, WFA, GPU dispatch, POA, LCSk++ sparse alignment, pair HMM, PAM40/120/200, BLOSUM30, CIGAR utilities, X-drop/Z-drop extension, spliced alignment | 290 | cyanea-core |
| **cyanea-omics** | Genomic coords, intervals, matrices, variants, AnnData, h5ad, zarr | 99 | cyanea-core, zarrs |
| **cyanea-stats** | Descriptive, correlation, hypothesis tests, distributions, PCA, effect sizes, Bayesian conjugate priors, combinatorics | 167 | cyanea-core |
| **cyanea-ml** | Clustering, distances, embeddings, KNN, PCA, t-SNE, UMAP, random forest, HMM | 161 | cyanea-core |
| **cyanea-chem** | SMILES/SDF V2000/V3000, fingerprints, MACCS keys, properties, substructure, stereochemistry, canonical SMILES | 79 | cyanea-core, sha2 |
| **cyanea-struct** | PDB, mmCIF, geometry, DSSP, Kabsch, contact maps, Ramachandran | 76 | cyanea-core, sha2 |
| **cyanea-phylo** | Newick/NEXUS, distances, UPGMA/NJ, Fitch/Sankoff, ML likelihood, bootstrap | 110 | cyanea-core, cyanea-ml (optional) |
| **cyanea-gpu** | Backend trait, CPU/CUDA/Metal backends, buffers, ops, benchmarks | 61 | cyanea-core, metal-rs, cudarc, criterion (bench) |
| **cyanea-wasm** | JSON-based WASM bindings (seq, io, align, stats, ml, chem, struct, phylo) | 96 | serde_json, wasm-bindgen |
| **cyanea-py** | Python bindings via PyO3 (seq, align, stats, ml, chem, struct, phylo, io) | — | pyo3 |

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
- `sam` — SAM format support in cyanea-io
- `bam` — BAM format support in cyanea-io (implies `sam`, depends on `flate2`)
- `wfa` — Wavefront alignment in cyanea-align
- `minhash` — MinHash/FracMinHash sketching in cyanea-seq
- `parallel` — Rayon parallelism in cyanea-align, cyanea-ml, cyanea-stats, cyanea-chem, cyanea-struct, cyanea-phylo, cyanea-io
- `blas` — Optional BLAS-backed PCA in cyanea-ml and cyanea-stats
- `simd` — SIMD-accelerated alignment in cyanea-align
- `h5ad` — HDF5-backed `.h5ad` file I/O in cyanea-omics (requires system HDF5 1.10.x; build with `HDF5_DIR`)
- `zarr` — Zarr v3 directory-based file I/O in cyanea-omics (pure Rust, zarrs 0.18)
- `cram` — CRAM alignment format in cyanea-io (implies `sam`, noodles-cram 0.72)
- `parquet` — Apache Parquet columnar format in cyanea-io (implies `vcf`, `bed`; arrow/parquet 54)

### Code patterns
- Each alignment module (NW, SW, semi-global, seed-extend) has its own private `push_cigar` helper
- WASM functions return `String` via `wasm_ok(val)` / `wasm_err(msg)` / `wasm_result(r)` JSON envelope
- Python bindings use `IntoPyResult` trait for `CyaneaError → PyErr` conversion
- Workspace shared deps in root `Cargo.toml` `[workspace.dependencies]`
- Seed-and-extend: `chain_seeds(seeds, k, max_gap)` scores each seed at `k` points with gap penalty
- Pattern matching algorithms in `cyanea-seq::pattern` work on `&[u8]` for generality; bitparallel methods (shift_and, bndm, myers) limited to patterns ≤ 64
- FMD-Index (`cyanea-seq::fmd_index`) indexes `text#revcomp(text)$` for bidirectional search and SMEM enumeration
- POA (`cyanea-align::poa`) uses Lee 2002 algorithm: DP alignment against DAG with topological sort
- CIGAR utilities (`cyanea-align::cigar`) — full SAM alphabet (9 ops: M/I/D/N/S/H/P/=/X), parsing, validation, coordinate queries, statistics, arithmetic (merge/reverse/split/collapse), alignment↔CIGAR conversion, MD tag generation

### Testing
- Unit tests inline (`#[cfg(test)]` modules at bottom of each file)
- Doc tests for key public API functions
- Test data: inline strings/vecs, no external fixture files
- `tempfile` crate for tests needing filesystem
- Property tests with `proptest` for round-trip and invariant verification
- Criterion benchmarks in cyanea-align, cyanea-seq, cyanea-ml, cyanea-stats, cyanea-chem, cyanea-struct, cyanea-gpu
- Fuzz targets with `cargo-fuzz` for all parsers (SMILES, SDF, PDB, Newick, NEXUS, VCF, BED, GFF3)

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
├── cyanea-io (+ cyanea-omics, csv, optional noodles)
├── cyanea-wasm (+ cyanea-seq/io/align/stats/ml/chem/struct/phylo, serde_json, wasm-bindgen)
└── cyanea-py (+ cyanea-seq/io/align/stats/ml/chem/struct/phylo, pyo3)
```

## Building Bindings

### WASM
```bash
cargo check -p cyanea-wasm --features wasm   # Check with wasm-bindgen
# Full build: wasm-pack build cyanea-wasm --features wasm
# TypeScript wrapper: npm run build:ts (from cyanea-wasm/)
```
Published as `@cyanea/bio` on npm (TypeScript types in `ts/types.ts`, wrapper in `ts/index.ts`).

### Python
```bash
# Needs Python 3.13 or 3.14 with compat flag
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py
# With NumPy interop:
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py --features numpy
# Full build: cd cyanea-py && PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```
Python bindings include: seq, align, stats, ml (pca/tsne/umap/kmeans/pairwise_distances + NumPy variants), chem, struct_bio, phylo, io (csv/vcf/bed/gff/sam/bam).

### NIF (Elixir)
The NIF crate lives in `../cyanea/native/cyanea_native/` and depends on labs crates via path.
```bash
cargo check -p cyanea-native   # Check only (can't link without BEAM)
# Full build: cd ../cyanea && mix compile
```

## CI & Benchmarking

- GitHub Actions CI: check, test, clippy, fmt, doc, WASM, Python (`.github/workflows/ci.yml`)
- Continuous benchmarking: Criterion on PRs with baseline comparison (`.github/workflows/bench.yml`)

## What's Not Done

- **Publishing**: Not yet on crates.io, PyPI, or npm (metadata ready, workflows pending)
- **Competitive benchmarks**: Head-to-head vs rust-bio not yet published
