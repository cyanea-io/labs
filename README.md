<div align="center">

<h1 align="center">
  Cyanea Labs
</h1>

<p align="center">
  <em>Modular bioinformatics in Rust.</em><br>
  <em>Sequences, alignments, structures, stats — native, WASM, Python.</em>
</p>

<p align="center">
  <a href="https://www.rust-lang.org/">
    <img alt="Rust" src="https://img.shields.io/badge/Rust-1.93+-000000?logo=rust&logoColor=white&style=for-the-badge">
  </a>
  <a href="Cargo.toml">
    <img alt="License" src="https://img.shields.io/badge/License-MIT%2FApache--2.0-c6a0f6?style=for-the-badge">
  </a>
</p>

<p align="center">
  <a href="docs/">Docs</a> &bull;
  <a href="ARCHITECTURE.md">Architecture</a> &bull;
  <a href="docs/BUILDING.md">Build Guide</a> &bull;
  <a href="https://github.com/cyanea-io/labs/issues">Issues</a> &bull;
  <a href="https://github.com/cyanea-io/labs/discussions">Discussions</a>
</p>

</div>

---

Cyanea Labs is a Cargo workspace of 13 crates covering the core primitives of computational biology — sequence analysis, alignment, genomic intervals, statistics, machine learning, cheminformatics, structural biology, and phylogenetics. Everything compiles to native, WebAssembly, and Python (via PyO3), with an Elixir NIF bridge for the Cyanea platform.

1570+ tests. Zero `unsafe`. No heavyweight C/C++ dependencies in the core path.

## Quick Start

```bash
# Check everything compiles
cargo check --workspace

# Run all tests (exclude cyanea-py unless Python env is configured)
cargo test --workspace --exclude cyanea-py

# Test a single crate
cargo test -p cyanea-seq

# Run benchmarks
cargo bench -p cyanea-align
```

### From Rust

```rust
use cyanea_seq::DnaSequence;
use cyanea_align::needleman_wunsch;
use cyanea_stats::describe;

// Sequence analysis
let seq = DnaSequence::new("ATCGATCG").unwrap();
assert_eq!(seq.gc_content(), 0.5);

// Global alignment
let result = needleman_wunsch(b"ACGT", b"ACGTT", 1, -1, -2).unwrap();
println!("{}", result.cigar); // 4=1I

// Descriptive statistics
let stats = describe(&[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
println!("mean={}, std={:.2}", stats.mean, stats.std_dev);
```

### From Python

```python
import cyanea

# Align two sequences
result = cyanea.smith_waterman("ACGTACGT", "CGTAC", 2, -1, -2)
print(result["score"], result["cigar"])

# Parse a VCF file
variants = cyanea.read_vcf("samples.vcf")

# PCA on a distance matrix
coords = cyanea.pca([[0,1,2],[1,0,3],[2,3,0]], n_components=2)
```

### From JavaScript / TypeScript

```typescript
import { align, seq, stats } from "@cyanea/bio";

const result = align.smithWaterman("ACGTACGT", "CGTAC", 2, -1, -2);
const gc = seq.gcContent("ATCGATCG");
const desc = stats.describe([1, 2, 3, 4, 5]);
```

## Crates

| Crate | What it does | Tests |
|-------|-------------|------:|
| **[cyanea-core](cyanea-core/)** | Error types, traits, SHA-256, zstd/gzip, mmap, log-space probability, rank/select bitvectors, wavelet matrix, Fenwick tree | 58 |
| **[cyanea-seq](cyanea-seq/)** | DNA/RNA/protein sequences, FASTA/FASTQ, k-mers, 2-bit encoding, suffix array, FM-index, BWT, FMD-index, MinHash, pattern matching (KMP, Boyer-Moore, Myers bit-parallel), PSSM/motif scanning, ORF finder, codon tables, sequence masking | 300 |
| **[cyanea-align](cyanea-align/)** | Needleman-Wunsch, Smith-Waterman, semi-global, banded, MSA, seed-and-extend, minimizers, WFA, POA, LCSk++, pair HMM, X-drop/Z-drop, spliced alignment, CIGAR utilities, substitution matrices | 290 |
| **[cyanea-omics](cyanea-omics/)** | Genomic coordinates, interval sets, genome arithmetic (intersect/subtract/complement/closest/Jaccard), expression matrices, sparse matrices, variants, gene annotations, coordinate liftover, AnnData/h5ad/zarr | 158 |
| **[cyanea-io](cyanea-io/)** | CSV, VCF, BED, BEDPE, GFF3, SAM, BAM, CRAM, Parquet | 81 |
| **[cyanea-stats](cyanea-stats/)** | Descriptive stats, correlation, hypothesis tests (t, chi-squared, Mann-Whitney, Fisher, KS), distributions, PCA, effect sizes, Bayesian conjugate priors, combinatorics | 167 |
| **[cyanea-ml](cyanea-ml/)** | K-means, DBSCAN, hierarchical clustering, pairwise distances, KNN, PCA, t-SNE, UMAP, random forest, HMM | 161 |
| **[cyanea-chem](cyanea-chem/)** | SMILES/SDF parsing, molecular fingerprints (Morgan, MACCS), substructure search, stereochemistry, canonical SMILES, physicochemical properties | 79 |
| **[cyanea-struct](cyanea-struct/)** | PDB/mmCIF parsing, 3D geometry, DSSP secondary structure, Kabsch superposition, contact maps, Ramachandran analysis | 76 |
| **[cyanea-phylo](cyanea-phylo/)** | Newick/NEXUS parsing, distance matrices, UPGMA/NJ tree building, Fitch/Sankoff parsimony, ML likelihood, bootstrap | 110 |
| **[cyanea-gpu](cyanea-gpu/)** | Backend trait with CPU, CUDA, and Metal implementations, GPU buffer management, benchmarks | 61 |
| **[cyanea-wasm](cyanea-wasm/)** | WebAssembly bindings via wasm-bindgen (seq, io, align, stats, ml, chem, struct, phylo) | 96 |
| **[cyanea-py](cyanea-py/)** | Python bindings via PyO3 (seq, align, stats, ml, chem, struct, phylo, io) with optional NumPy support | &mdash; |

## Architecture

> See [`ARCHITECTURE.md`](ARCHITECTURE.md) for the full architecture guide with ASCII diagrams, feature flag details, data flow pipelines, and platform support matrix.

```
cyanea-core (foundation)
├── cyanea-seq          Sequences, indexing, k-mers
├── cyanea-align        Pairwise & multiple alignment
├── cyanea-omics        Genomic intervals, matrices, variants
├── cyanea-io           File format I/O
├── cyanea-stats        Statistics & distributions
├── cyanea-ml           Machine learning & clustering
├── cyanea-chem         Chemical structure & fingerprints
├── cyanea-struct       Protein structure & geometry
├── cyanea-phylo        Phylogenetics (optional: cyanea-ml)
├── cyanea-gpu          GPU compute backends
├── cyanea-wasm         → WebAssembly (@cyanea/bio on npm)
└── cyanea-py           → Python (PyO3 + maturin)
```

Every domain crate depends only on `cyanea-core`. The binding crates (`wasm`, `py`) aggregate the domain crates to expose a unified API. The Elixir NIF bridge lives in a [separate repo](https://github.com/cyanea-io/cyanea) and depends on these crates via path.

## Feature Flags

All domain crates default to `std`. Opt into additional capabilities:

| Flag | Scope | What it enables |
|------|-------|-----------------|
| `parallel` | align, ml, stats, chem, struct, phylo, io | Rayon parallelism |
| `simd` | align | SIMD-accelerated alignment |
| `cuda` | gpu, align | CUDA GPU backend (requires CUDA toolkit) |
| `metal` | gpu, align | Metal GPU backend (macOS) |
| `blas` | ml, stats | BLAS-backed PCA |
| `wfa` | align | Wavefront alignment |
| `minhash` | seq | MinHash / FracMinHash sketching |
| `sam` | io | SAM format support |
| `bam` | io | BAM format (implies `sam`) |
| `cram` | io | CRAM format (implies `sam`) |
| `parquet` | io | Apache Parquet columnar format |
| `h5ad` | omics | HDF5-backed AnnData I/O |
| `zarr` | omics | Zarr v3 I/O |
| `serde` | all | Serialization support |
| `wasm` | wasm | wasm-bindgen annotations |
| `numpy` | py | NumPy array interop |

## Building Bindings

> See [`docs/BUILDING.md`](docs/BUILDING.md) for the complete build guide covering all targets, feature combinations, and troubleshooting.

### WebAssembly

```bash
cargo check -p cyanea-wasm --features wasm
wasm-pack build cyanea-wasm --features wasm
```

Published as `@cyanea/bio` on npm. TypeScript types and a Web Worker API are included in `cyanea-wasm/ts/`.

### Python

```bash
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo check -p cyanea-py
cd cyanea-py && PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

### Elixir NIF

```bash
cargo check -p cyanea-native   # Check only (linking requires BEAM)
cd ../cyanea && mix compile     # Full build via Rustler
```

## Contributing

```bash
# Full CI check
cargo fmt --all -- --check
cargo clippy --workspace --exclude cyanea-py
cargo test --workspace --exclude cyanea-py
```

**Using Claude Code?** See [`CLAUDE.md`](CLAUDE.md) for project context, architecture, and coding conventions.

## License

MIT OR Apache-2.0
