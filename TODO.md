# Cyanea Labs — Open Work

> Best-in-class Rust bioinformatics ecosystem targeting CPU, GPU, WASM (browser), and BEAM (Elixir NIFs).

Last updated: 2026-02-09

---

## Status: All 13 crates complete (659+ tests)

Every crate has a `docs/STATUS.md` with full API documentation.

| Crate | Status | Tests |
|-------|--------|------:|
| cyanea-core | Complete | 14 |
| cyanea-seq | Complete | 44 |
| cyanea-io | Complete (CSV, VCF, BED, GFF3) | 26 |
| cyanea-align | Complete (NW, SW, semi-global, MSA, banded, GPU dispatch) | 77 |
| cyanea-omics | Complete (genomic coords, intervals, matrices, variants, AnnData) | 89 |
| cyanea-stats | Complete (descriptive, correlation, testing, distributions, PCA) | 77 |
| cyanea-ml | Complete (clustering, distances, embeddings, KNN, regression, PCA, t-SNE) | 109 |
| cyanea-chem | Complete (SMILES, SDF, fingerprints, properties, substructure) | 35 |
| cyanea-struct | Complete (PDB, geometry, DSSP, Kabsch, contact maps) | 37 |
| cyanea-phylo | Complete (Newick, NEXUS, distances, UPGMA/NJ, Fitch/Sankoff) | 73 |
| cyanea-gpu | CPU backend complete; CUDA/Metal stubs | 45 |
| cyanea-wasm | Complete (JSON API, wasm-bindgen annotations) | 56 |
| cyanea-py | Complete (seq, align, stats, core, ml submodules) | — |

---

## P0 — GPU backends (requires hardware SDKs)

### CUDA backend (`cyanea-gpu`)
- [ ] Metal compute shaders for matrix ops, pairwise distance, reductions
- [ ] CUDA kernels for same operation set
- [ ] GPU-accelerated Smith-Waterman batch alignment in `cyanea-align`
- [ ] GPU-accelerated pairwise distance matrices

### CPU parallelism (across crates)
- [ ] Add `rayon` as optional workspace dependency behind `parallel` feature
- [ ] Parallel batch alignment, distance matrices, fingerprint generation, etc.

---

## P1 — Algorithmic extensions

### Alignment (`cyanea-align`)
- [ ] SIMD-accelerated Smith-Waterman (striped, 8/16-way; Farrar 2007)
- [ ] Seed-and-extend heuristic alignment
- [ ] Minimizer-based seeding

### File formats (`cyanea-io`)
- [ ] SAM/BAM/CRAM parser
- [ ] Parquet reader/writer
- [ ] HDF5/Zarr reader

### Statistics (`cyanea-stats`)
- [ ] Chi-squared, F, Binomial distributions
- [ ] Fisher's exact test, chi-squared test of independence
- [ ] Effect sizes (Cohen's d, odds ratio)
- [ ] UMAP dimensionality reduction

### Phylogenetics (`cyanea-phylo`)
- [ ] Maximum likelihood tree inference (JC69 + Felsenstein pruning)
- [ ] Bootstrap support values

### Chemistry (`cyanea-chem`)
- [ ] Stereochemistry interpretation (R/S, E/Z)
- [ ] Canonical SMILES output
- [ ] SDF V3000 parser
- [ ] MACCS 166-key fingerprints

### Structural biology (`cyanea-struct`)
- [ ] mmCIF/PDBx parser
- [ ] Ramachandran validation and outlier detection
- [ ] B-factor analysis

### Sequences (`cyanea-seq`)
- [ ] Suffix arrays / FM-index
- [ ] 2-bit DNA encoding

### Linear algebra
- [ ] General dense matrix type with LU/eigen/SVD
- [ ] Decide: pure Rust vs optional BLAS/LAPACK backend

---

## P2 — Platform integration

### WASM (`cyanea-wasm`)
- [ ] Build and publish as npm package (@cyanea/*)
- [ ] TypeScript type definitions
- [ ] Web Worker-friendly API
- [ ] Add bindings for chem, struct, phylo modules
- [ ] Streaming FASTA/FASTQ parsing in browser (ReadableStream)

### Elixir NIFs (`cyanea/native/cyanea_native`)
- [ ] Add NIF bindings for chem, struct, phylo
- [ ] Dirty scheduler annotations for long-running computations
- [ ] Resource types for large objects (Molecule, Structure, PhyloTree)

### Python (`cyanea-py`)
- [ ] NumPy interop for matrices and distance results
- [ ] Publish to PyPI

---

## P3 — Quality and ecosystem

### Testing
- [ ] Integration tests across crate boundaries
- [ ] Property-based testing (proptest) for parsers
- [ ] Fuzzing targets for all parsers
- [ ] Benchmark suite (criterion)

### CI/CD
- [ ] `cargo check --target wasm32-unknown-unknown` for WASM
- [ ] Feature matrix: test with and without `std`, `serde`, `parallel`

### Documentation
- [ ] Crate-level doc examples for every public module
- [ ] Cross-crate usage guide
- [ ] Architecture diagram

### Publishing
- [ ] Publish to crates.io
- [ ] Publish to PyPI (maturin)
- [ ] Publish to npm (wasm-pack)
