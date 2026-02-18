# Documentation

Cross-crate documentation index for the Cyanea Labs workspace.

## Getting Started

| Document | Description |
|----------|-------------|
| [README](../README.md) | Project overview, quick start, code examples (Rust / Python / TypeScript) |
| [Build Guide](BUILDING.md) | How to compile for every target: native, WASM, Python, CUDA, Metal, WebGPU, Elixir NIF |
| [Architecture](../ARCHITECTURE.md) | Dependency graph, feature flag architecture, data flow diagrams, platform support |

## Crate Documentation

Each crate has two documents: a **README** (user-facing overview, quick start, feature flags) and a **STATUS** (complete API reference with every public type and function).

### Foundation

| Crate | README | API Reference |
|-------|--------|---------------|
| **cyanea-core** | [README](../cyanea-core/README.md) | [STATUS](../cyanea-core/docs/STATUS.md) |

### Domain Crates

| Crate | README | API Reference | Area |
|-------|--------|---------------|------|
| **cyanea-seq** | [README](../cyanea-seq/README.md) | [STATUS](../cyanea-seq/docs/STATUS.md) | Sequences, FASTA/FASTQ, k-mers, indexing, trimming |
| **cyanea-align** | [README](../cyanea-align/README.md) | [STATUS](../cyanea-align/docs/STATUS.md) | NW/SW alignment, MSA, POA, HMMs, GPU dispatch |
| **cyanea-io** | [README](../cyanea-io/README.md) | [STATUS](../cyanea-io/docs/STATUS.md) | 15 file format parsers (VCF, BAM, Parquet, ...) |
| **cyanea-omics** | [README](../cyanea-omics/README.md) | [STATUS](../cyanea-omics/docs/STATUS.md) | Genomic intervals, variants, AnnData, h5ad/zarr |
| **cyanea-stats** | [README](../cyanea-stats/README.md) | [STATUS](../cyanea-stats/docs/STATUS.md) | Statistics, distributions, popgen, survival |
| **cyanea-ml** | [README](../cyanea-ml/README.md) | [STATUS](../cyanea-ml/docs/STATUS.md) | Clustering, trees, GBDT, HMM, PCA/t-SNE/UMAP |
| **cyanea-chem** | [README](../cyanea-chem/README.md) | [STATUS](../cyanea-chem/docs/STATUS.md) | SMILES/SDF, fingerprints, stereochemistry |
| **cyanea-struct** | [README](../cyanea-struct/README.md) | [STATUS](../cyanea-struct/docs/STATUS.md) | PDB/mmCIF, DSSP, Kabsch, Ramachandran |
| **cyanea-phylo** | [README](../cyanea-phylo/README.md) | [STATUS](../cyanea-phylo/docs/STATUS.md) | Newick/NEXUS, UPGMA/NJ, ML, bootstrap |
| **cyanea-gpu** | [README](../cyanea-gpu/README.md) | [STATUS](../cyanea-gpu/docs/STATUS.md) | CPU/CUDA/Metal/WebGPU backends |

### Binding Crates

| Crate | README | API Reference | Target |
|-------|--------|---------------|--------|
| **cyanea-wasm** | [README](../cyanea-wasm/README.md) | [STATUS](../cyanea-wasm/docs/STATUS.md) | WebAssembly (`@cyanea/bio` on npm) |
| **cyanea-py** | [README](../cyanea-py/README.md) | [STATUS](../cyanea-py/docs/STATUS.md) | Python (`pip install cyanea` via maturin) |

## Developer Reference

| Document | Description |
|----------|-------------|
| [CLAUDE.md](../CLAUDE.md) | Workspace conventions, crate map, dependency graph, coding patterns, testing strategy |
| [CI workflow](../.github/workflows/ci.yml) | GitHub Actions: check, test, clippy, fmt, doc, WASM, Python |
| [Bench workflow](../.github/workflows/bench.yml) | Criterion benchmarks on PRs with baseline comparison |

## Documentation Map

```
labs/
├── README.md              ← Project overview & quick start
├── ARCHITECTURE.md        ← Dependency graph, feature flags, data flow, platforms
├── CLAUDE.md              ← Conventions & patterns for contributors
├── docs/
│   ├── README.md          ← This file (documentation index)
│   └── BUILDING.md        ← Build guide for all targets
├── .github/workflows/
│   ├── ci.yml             ← CI pipeline
│   └── bench.yml          ← Benchmark pipeline
└── cyanea-*/
    ├── README.md          ← Per-crate overview, quick start, feature flags
    └── docs/
        └── STATUS.md      ← Per-crate full API reference
```

## Finding What You Need

| I want to... | Go to |
|--------------|-------|
| Get started quickly | [Root README](../README.md) |
| Understand the architecture | [ARCHITECTURE.md](../ARCHITECTURE.md) |
| Build for a specific target | [BUILDING.md](BUILDING.md) |
| Learn what a crate does | Crate README (linked above) |
| Look up a specific function | Crate STATUS.md (linked above) |
| Understand coding conventions | [CLAUDE.md](../CLAUDE.md) |
| See which features exist | [ARCHITECTURE.md -- Feature Flags](../ARCHITECTURE.md#feature-flag-architecture) |
| Check platform support | [ARCHITECTURE.md -- Platform Matrix](../ARCHITECTURE.md#platform-support-matrix) |
| Troubleshoot a build | [BUILDING.md -- Troubleshooting](BUILDING.md#troubleshooting) |
