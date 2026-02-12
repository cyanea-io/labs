# Cyanea Labs — Open Work

> Best-in-class Rust bioinformatics ecosystem targeting CPU, GPU, WASM (browser), and BEAM (Elixir NIFs).

Last updated: 2026-02-12

---

## Status: All 13 crates complete (1055+ tests)

Every crate has a `docs/STATUS.md` with full API documentation.

| Crate | Status | Tests |
|-------|--------|------:|
| cyanea-core | Complete | 14 |
| cyanea-seq | Complete (FASTA, FASTQ, k-mers, 2-bit encoding, suffix arrays, FM-index, MinHash) | 111 |
| cyanea-io | Complete (CSV, VCF, BED, GFF3, SAM, BAM) | 56 |
| cyanea-align | Complete (NW, SW, semi-global, banded, MSA, SIMD SW, seed-and-extend, minimizers, WFA, GPU dispatch) | 138 |
| cyanea-omics | Complete (genomic coords, intervals, matrices, variants, AnnData, single-cell) | 87 |
| cyanea-stats | Complete (descriptive, correlation, hypothesis tests, distributions, PCA, effect sizes, chi-squared, Fisher's exact) | 127 |
| cyanea-ml | Complete (clustering, distances, embeddings, KNN, PCA, t-SNE, UMAP, random forest, regression, HMM) | 161 |
| cyanea-chem | Complete (SMILES, SDF V2000/V3000, fingerprints, MACCS keys, properties, substructure, stereochemistry, canonical SMILES) | 79 |
| cyanea-struct | Complete (PDB, mmCIF, geometry, DSSP, Kabsch, contact maps, Ramachandran) | 76 |
| cyanea-phylo | Complete (Newick, NEXUS, distances, UPGMA/NJ, Fitch/Sankoff, ML likelihood, bootstrap) | 101 |
| cyanea-gpu | CPU backend complete; CUDA/Metal kernel stubs | 43 |
| cyanea-wasm | Complete (JSON API for seq, io, align, stats, ml, chem, struct, phylo) | 83 |
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
- [x] Metal compute shader stubs (matrix ops, pairwise distance, reductions)
- [x] CUDA kernel stubs
- [x] GPU-accelerated Smith-Waterman batch alignment dispatch in `cyanea-align`
- [ ] Metal kernel integration testing on Apple Silicon
- [ ] GPU alignment benchmarks (Criterion)

### Competitive benchmarks
- [ ] Head-to-head vs rust-bio (SW, FM-index)
- [ ] Benchmark matrix: scalar vs SSE4.1 vs AVX2 vs NEON
- [ ] Results published in `benchmarks/README.md`

---

## P1 — Ecosystem completeness

### File formats (`cyanea-io`)
- [x] SAM text parser (feature-gated)
- [x] BAM parser with BGZF decompression (feature-gated behind `bam`)
- [ ] CRAM parser
- [ ] Parquet reader/writer
- [ ] HDF5/Zarr reader

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
- [ ] HDF5-backed AnnData (`.h5ad` format)

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
- [ ] Add NIF bindings for chem, struct, phylo
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
