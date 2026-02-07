# Cyanea Labs — Open Work

> Best-in-class Rust bioinformatics ecosystem targeting CPU, GPU, WASM (browser), and BEAM (Elixir NIFs).

Last updated: 2026-02-07

---

## P0 — Foundation gaps blocking real-world use

### ~~Sequence types (`cyanea-seq`)~~ ✓
- [x] `DnaSequence`, `RnaSequence`, `ProteinSequence` types implementing `Sequence` trait
- [x] Reverse complement, transcription (DNA→RNA), translation (RNA→Protein with codon table)
- [x] K-mer iterator (lazy, zero-alloc where possible)
- [x] Quality score representation (Phred33/64) and FASTQ record type
- [x] FASTQ parsing (via needletail, already a dep)
- [x] Sequence validation (IUPAC alphabets)

### File format I/O (`cyanea-io`)
- [ ] VCF parser (at minimum: header, INFO/FORMAT fields, genotype extraction)
- [ ] GFF3/GTF parser → `cyanea-omics::Gene`/`Transcript`/`Exon`
- [ ] BED parser (BED3/BED6/BED12) → `GenomicInterval`
- [ ] SAM parser (header + alignment records; BAM via optional `flate2`/`htslib`)
- [ ] Streaming/chunked parsing for large files (trait-based `RecordIterator`)
- [ ] Compressed file auto-detection (delegate to `cyanea-core::compress`)

### Linear algebra foundation
- [ ] General dense matrix type (row-major `Vec<f64>`, m×n) — could live in `cyanea-stats` or a new `cyanea-linalg`
- [ ] LU decomposition, eigenvalue decomposition, general SVD
- [ ] Decide: pure Rust vs optional BLAS/LAPACK backend behind feature flag
- [ ] This unblocks: PCA, t-SNE, UMAP, spectral clustering, GMM

---

## P1 — Algorithmic completions

### Dimensionality reduction (`cyanea-stats::reduction` + `cyanea-ml::reduction`)
- [ ] PCA (via eigendecomposition of covariance matrix)
- [ ] t-SNE (Barnes-Hut approximation for O(n log n))
- [ ] UMAP (nearest-neighbor graph + SGD layout)
- [ ] Consolidate into one crate — likely `cyanea-ml::reduction`

### Alignment (`cyanea-align`)
- [ ] Semi-global alignment (one sequence end-free)
- [ ] Banded Needleman-Wunsch (O(nd) for similar sequences)
- [ ] SIMD-accelerated Smith-Waterman (striped, 8/16-way; reference: Farrar 2007)
- [ ] Multiple sequence alignment — progressive (guide tree + profile alignment)

### Omics (`cyanea-omics`)
- [ ] Single-cell module: AnnData-like container (X matrix, obs, var, obsm layers)
- [ ] Interval tree for O(log n + k) overlap queries in `IntervalSet`
- [ ] Sparse matrix CSR/CSC formats (needed for scRNA-seq count matrices)

### Statistics (`cyanea-stats`)
- [ ] Chi-squared, F, Binomial distributions
- [ ] Fisher's exact test, chi-squared test of independence
- [ ] Effect sizes (Cohen's d, odds ratio)
- [ ] Survival analysis (Kaplan-Meier, log-rank test)

### Phylogenetics (`cyanea-phylo`)
- [ ] Maximum likelihood tree inference (at least JC69 + Felsenstein pruning)
- [ ] Bootstrap support values
- [ ] NEXUS format I/O
- [ ] Ancestral state reconstruction (marginal, via post-order traversal)

### Chemistry (`cyanea-chem`)
- [ ] Stereochemistry interpretation (R/S, E/Z from parsed @ and / tokens)
- [ ] Canonical SMILES output (for deduplication)
- [ ] SDF V3000 parser
- [ ] MACCS 166-key fingerprints

### Structural biology (`cyanea-struct`)
- [ ] mmCIF/PDBx parser (the modern replacement for PDB format)
- [ ] Ramachandran validation and outlier detection
- [ ] B-factor analysis and disorder prediction
- [ ] Mass-weighted superposition

---

## P2 — GPU and parallelism

### GPU backends (`cyanea-gpu`)
- [ ] Metal backend (macOS/iOS) — compute shaders for matrix ops, pairwise distance, reductions
- [ ] CUDA backend — same kernel set
- [ ] Unified dispatch: auto-select best available backend at runtime
- [ ] GPU-accelerated Smith-Waterman (batch short-read alignment)
- [ ] GPU-accelerated pairwise distance matrices (for large k-mer/fingerprint datasets)
- [ ] GPU-accelerated Morgan fingerprint generation (batch molecules)

### CPU parallelism (across crates)
- [ ] Add `rayon` as optional workspace dependency behind `parallel` feature
- [ ] `cyanea-align`: parallel batch alignment
- [ ] `cyanea-ml`: parallel k-means iterations, parallel pairwise distance matrix
- [ ] `cyanea-chem`: parallel SDF parsing, parallel fingerprint generation, parallel Tanimoto bulk
- [ ] `cyanea-stats`: parallel correlation matrix computation
- [ ] `cyanea-struct`: parallel contact map (all-atom), parallel RMSD over trajectory frames
- [ ] `cyanea-omics`: parallel expression matrix operations

---

## P3 — Platform integration

### WASM (`cyanea-wasm`)
- [ ] Add `wasm-bindgen` annotations to all existing functions
- [ ] Build and publish as npm package
- [ ] Streaming FASTA/FASTQ parsing in browser (ReadableStream → records)
- [ ] Web Worker-friendly API (transferable ArrayBuffers)
- [ ] Add WASM bindings for chem (SMILES→properties, fingerprints, Tanimoto)
- [ ] Add WASM bindings for struct (PDB parsing, geometry, contact maps)
- [ ] Add WASM bindings for phylo (Newick parse/write, tree construction)
- [ ] Benchmark and optimize: minimize allocations, use `wee_alloc`

### Elixir NIFs (`cyanea/native/cyanea_native`)
- [ ] Audit current NIF bridge — which labs crates are already exposed?
- [ ] Add NIF bindings for chem (SMILES parsing, properties, fingerprints, substructure search)
- [ ] Add NIF bindings for struct (PDB parsing, superposition, contact maps)
- [ ] Add NIF bindings for phylo (Newick I/O, tree construction, comparison)
- [ ] Add NIF bindings for seq types (DnaSequence, RnaSequence, ProteinSequence, FASTQ parsing)
- [ ] Dirty scheduler annotations for long-running computations (alignment, fingerprint bulk)
- [ ] Resource types for large objects (Molecule, Structure, PhyloTree) to avoid repeated serialization

### Python (`cyanea-py`) — low priority
- [ ] PyO3 + maturin setup
- [ ] Expose core functions: alignment, stats, k-mers, fingerprints
- [ ] NumPy interop for matrices and distance results

---

## P4 — Quality and ecosystem

### Testing
- [ ] Integration tests across crate boundaries (FASTA → sequence → alignment → stats)
- [ ] Property-based testing (proptest) for parsers (SMILES, PDB, Newick, VCF)
- [ ] Fuzzing targets for all parsers
- [ ] Benchmark suite (criterion) for alignment, fingerprints, distance matrices, PDB parsing

### CI/CD
- [ ] Workspace-wide `cargo check`, `cargo test`, `cargo clippy` in CI
- [ ] `cargo check --target wasm32-unknown-unknown` for WASM compatibility
- [ ] `cargo check -p cyanea-native` (NIF crate) in CI
- [ ] Feature matrix: test with and without `std`, `serde`, `parallel`

### Documentation
- [ ] Crate-level doc examples for every public module
- [ ] Cross-crate usage guide (how the pieces fit together)
- [ ] Architecture diagram (dependency graph of the 13 crates)

---

## Completed crates (for reference)

| Crate | Status | Tests |
|-------|--------|------:|
| cyanea-core | Complete | 14 |
| cyanea-align | Core complete (NW, SW, batch, scoring) | 42 |
| cyanea-omics | Core complete (intervals, expr, sparse, variant, annotation) | 75 |
| cyanea-stats | Core complete (descriptive, correlation, testing, correction, distributions) | 70 |
| cyanea-ml | Core complete (distance, k-mer, encoding, clustering, normalize, evaluate) | 77 |
| cyanea-chem | Complete | 34 |
| cyanea-struct | Complete | 36 |
| cyanea-phylo | Core complete (tree, newick, distance, compare, construct) | 53 |
| cyanea-gpu | CPU backend complete | 43 |
| cyanea-seq | Core complete (alphabet, sequence types, codon, k-mer, quality, FASTQ) | 44 |
| cyanea-wasm | Functions complete, needs wasm-bindgen | 35 |
