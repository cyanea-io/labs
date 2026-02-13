# Cyanea Labs — Open Work

> Best-in-class Rust bioinformatics ecosystem targeting CPU, GPU, WASM (browser), and BEAM (Elixir NIFs).

Last updated: 2026-02-13

---

## Status: All 13 crates complete (1440+ tests)

Every crate has a `docs/STATUS.md` with full API documentation.

| Crate | Status | Tests |
|-------|--------|------:|
| cyanea-core | Complete (traits, errors, hashing, compression, mmap, probability types, bitvectors, Fenwick tree) | 58 |
| cyanea-seq | Complete (FASTA, FASTQ, k-mers, 2-bit encoding, suffix arrays, FM-index, FMD-Index, BWT, MinHash, pattern matching, PSSM, ORF finder, FASTA indexed reader, codon tables, masking) | 300 |
| cyanea-io | Complete (CSV, VCF, BED, BEDPE, GFF3, SAM, BAM, CRAM, Parquet) | 81 |
| cyanea-align | Complete (NW, SW, semi-global, banded, MSA, SIMD SW, seed-and-extend, minimizers, WFA, GPU dispatch, POA, LCSk++, pair HMM, PAM40/120/200, BLOSUM30, CIGAR utilities, X-drop/Z-drop, spliced alignment) | 290 |
| cyanea-omics | Complete (genomic coords, intervals, matrices, variants, AnnData, h5ad, zarr) | 99 |
| cyanea-stats | Complete (descriptive, correlation, hypothesis tests, distributions, PCA, effect sizes, chi-squared, Fisher's exact, Bayesian conjugate priors, combinatorics) | 167 |
| cyanea-ml | Complete (clustering, distances, embeddings, KNN, PCA, t-SNE, UMAP, random forest, regression, HMM) | 161 |
| cyanea-chem | Complete (SMILES, SDF V2000/V3000, fingerprints, MACCS keys, properties, substructure, stereochemistry, canonical SMILES) | 79 |
| cyanea-struct | Complete (PDB, mmCIF, geometry, DSSP, Kabsch, contact maps, Ramachandran) | 76 |
| cyanea-phylo | Complete (Newick, NEXUS, distances, UPGMA/NJ, Fitch/Sankoff, ML likelihood, bootstrap) | 110 |
| cyanea-gpu | Complete (CPU + CUDA + Metal backends, buffers, ops) | 61 |
| cyanea-wasm | Complete (JSON API for seq, io, align, stats, ml, chem, struct, phylo) | 96 |
| cyanea-py | Complete (seq, align, stats, ml, chem, struct, phylo, io submodules) | — |

---

## P0 — Performance & SIMD

### AVX2 SIMD Smith-Waterman (`cyanea-align`)
- [x] NEON (aarch64) 8×i16 striped SW score-only
- [x] SSE4.1 (x86_64) 8×i16 striped SW score-only
- [x] AVX2 (x86_64) 16×i16 striped SW score-only
- [x] Runtime dispatch: AVX2 → SSE4.1 → scalar

### Wavefront Alignment (WFA)
- [x] Gap-affine WFA (O(ns) complexity) in `cyanea-align`
- [x] Feature-gated behind `wfa`
- [x] Asymmetric gap costs for correct score with match rewards
- [x] Proptest verification: WFA score equals NW score

### GPU backends (`cyanea-gpu`)
- [x] Metal compute shaders (matrix ops, pairwise distance, reductions)
- [x] CUDA kernels (cudarc/NVRTC, f64)
- [x] GPU-accelerated Smith-Waterman batch alignment dispatch in `cyanea-align`
- [x] Metal kernel integration testing on Apple Silicon (61 tests)
- [x] GPU alignment benchmarks (Criterion)

### Competitive benchmarks
- [x] Head-to-head vs rust-bio (SW, FM-index)
- [x] Benchmark matrix: scalar vs SSE4.1 vs AVX2 vs NEON
- [x] Results published in `benchmarks/README.md`

---

## P1 — Ecosystem completeness

### File formats (`cyanea-io`)
- [x] SAM text parser (feature-gated)
- [x] BAM parser with BGZF decompression (feature-gated behind `bam`)
- [x] CRAM parser (feature-gated behind `cram`, noodles-cram 0.72)
- [x] Parquet reader/writer (feature-gated behind `parquet`, arrow/parquet 54)

### File formats (`cyanea-omics`)
- [x] Zarr v3 directory-based I/O (feature-gated behind `zarr`, zarrs 0.18)
- [x] HDF5-backed AnnData `.h5ad` format (feature-gated behind `h5ad`)

### Sequences (`cyanea-seq`)
- [x] Suffix arrays
- [x] FM-index
- [x] 2-bit DNA encoding
- [x] MinHash / FracMinHash sketching (feature-gated behind `minhash`)

### Alignment (`cyanea-align`)
- [x] SIMD-accelerated Smith-Waterman (Farrar striped, NEON + SSE4.1 + AVX2)
- [x] Seed-and-extend heuristic alignment
- [x] Minimizer-based seeding
- [x] Wavefront Alignment (WFA) with gap-affine scoring

### Statistics (`cyanea-stats`)
- [x] Chi-squared, F, Binomial distributions
- [x] Fisher's exact test, chi-squared test of independence
- [x] Effect sizes (Cohen's d, odds ratio)

### ML (`cyanea-ml`)
- [x] UMAP dimensionality reduction
- [x] Random forests
- [x] Hidden Markov Models (forward, Viterbi, Baum-Welch)

### Phylogenetics (`cyanea-phylo`)
- [x] Maximum likelihood tree inference (JC69 + Felsenstein pruning)
- [x] Bootstrap support values

### Chemistry (`cyanea-chem`)
- [x] Stereochemistry interpretation (R/S, E/Z)
- [x] Canonical SMILES output
- [x] SDF V3000 parser (auto-detection in `parse_sdf`)
- [x] MACCS 166-key fingerprints

### Structural biology (`cyanea-struct`)
- [x] mmCIF/PDBx parser
- [x] Ramachandran validation and outlier detection

### Omics (`cyanea-omics`)
- [x] AnnData in-memory
- [x] HDF5-backed AnnData (`.h5ad` format, feature-gated behind `h5ad`)
- [x] Zarr v3 directory-based I/O (feature-gated behind `zarr`)

---

## P1.5 — Competitive Parity (rust-bio gaps)

### High-Value (implement now)

#### Pattern matching suite (`cyanea-seq`)
- [x] Myers bit-parallel approximate matching
- [x] Horspool (Boyer-Moore-Horspool)
- [x] KMP (Knuth-Morris-Pratt)
- [x] Shift-And bitparallel
- [x] BNDM bitparallel
- [x] BOM (Backward Oracle Matching)
- [x] Ukkonen cut-off approximate matching

#### Sequence analysis (`cyanea-seq`)
- [x] PSSM / Motif scanning (DNA + Protein)
- [x] ORF finder (configurable start/stop codons, min length, all 6 frames)
- [x] FASTA indexed reader (.fai) for random access
- [x] FMD-Index (bidirectional FM-index for strand-aware search)

#### Alignment algorithms (`cyanea-align`)
- [x] Partial Order Alignment (POA) for long-read consensus
- [x] Sparse alignment (LCSk++) for long sequences
- [x] Pair HMM (Match/Insert/Delete states, log-space)

### Medium-Value

- [x] Additional PAM matrices (PAM40, PAM120, PAM200) + BLOSUM30 → `cyanea-align`
- [x] BWT utilities (standalone construction/querying) → `cyanea-seq`
- [x] LogProb / PhredProb newtypes → `cyanea-core`
- [x] Bayesian statistics (conjugate priors: Beta, Gamma, NormalConjugate, Dirichlet) → `cyanea-stats`
- [x] Rank/select bitvectors, wavelet matrix → `cyanea-core`
- [x] Fenwick tree → `cyanea-core`
- [x] BEDPE format → `cyanea-io`
- [x] Combinatorics utilities (factorial, binomial, permutations, multinomial, combinations iterator) → `cyanea-stats`

---

## P2 — Platform integration & distribution

### WASM (`cyanea-wasm`)
- [x] JSON API with wasm-bindgen annotations for 8 modules
- [x] TypeScript type definitions (`ts/types.ts`, `ts/index.ts`)
- [x] npm package config (`@cyanea/bio` in `package.json`)
- [ ] npm package publishing (@cyanea/bio)
- [ ] Web Worker-friendly API

### Python (`cyanea-py`)
- [x] PyO3 bindings for 8 modules
- [x] NumPy interop for matrices and distance results (optional `numpy` feature)
- [x] SAM/BAM bindings (`parse_sam`, `sam_stats`, `parse_bam`, `bam_stats`)
- [x] PCA, t-SNE, k-means, pairwise_distances bindings
- [ ] Publish to PyPI

### Elixir NIFs (`cyanea/native/cyanea_native`)
- [x] Add NIF bindings for chem, struct, phylo (28 NIFs in native crate + Elixir domain modules)
- [ ] Dirty scheduler annotations for long-running computations

### Publishing
- [x] crates.io metadata: keywords + categories on all 13 crates
- [ ] Publish to crates.io (all 11 domain crates)
- [ ] Publish to PyPI (maturin)
- [ ] Publish to npm (wasm-pack)

---

## P3 — Quality and ecosystem

### Testing
- [x] Property-based testing (proptest) for parsers and alignment
- [x] Fuzz targets for all parsers (SMILES, SDF, PDB, Newick, NEXUS, VCF, BED, GFF3)
- [x] Criterion benchmark suite (align, seq, ml, stats, chem, struct, gpu)
- [ ] Integration tests across crate boundaries

### CI/CD
- [x] GitHub Actions CI: check, test, clippy, fmt, doc, WASM, Python
- [x] Continuous benchmarking with baseline comparison
- [ ] `cargo check --target wasm32-unknown-unknown` for WASM in CI
- [ ] Feature matrix: test with and without `std`, `serde`, `parallel`

### Documentation
- [x] Crate-level doc examples for public modules
- [ ] mdBook documentation site
- [ ] Cross-crate usage guide
- [ ] Architecture diagram

---

## New Features

See **[ROADMAP.md](ROADMAP.md)** for the full feature/capability roadmap (T1–T8 tiers).
