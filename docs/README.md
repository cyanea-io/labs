# Documentation

## I want to...

| Goal | Document |
|------|----------|
| Get started quickly | [Root README](../README.md) |
| Understand the architecture | [Architecture](ARCHITECTURE.md) |
| Build for a specific target | [Build Guide](BUILDING.md) |
| Use Python/WASM/Elixir bindings | [Bindings Guide](BINDINGS.md) |
| Build a cross-crate pipeline | [Usage Guide](GUIDE.md) |
| See benchmark results | [Benchmarks](BENCHMARKS.md) |
| Learn a specific crate's API | Per-crate docs (below) |
| Understand coding conventions | [CLAUDE.md](../CLAUDE.md) |

## Per-Crate Documentation

Each crate has three documents in its `docs/` directory:
- **API.md** -- Full API reference (every public type and function)
- **GUIDE.md** -- Usage guide with code examples
- **ARCHITECTURE.md** -- Internal design and module map

### Foundation

| Crate | API | Guide | Architecture |
|-------|-----|-------|--------------|
| **cyanea-core** | [API.md](../cyanea-core/docs/API.md) | [GUIDE.md](../cyanea-core/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-core/docs/ARCHITECTURE.md) |

### Domain Crates

| Crate | Area | API | Guide | Architecture |
|-------|------|-----|-------|--------------|
| **cyanea-seq** | Sequences, FASTA/FASTQ, k-mers, FM-index, trimming, MinHash | [API.md](../cyanea-seq/docs/API.md) | [GUIDE.md](../cyanea-seq/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-seq/docs/ARCHITECTURE.md) |
| **cyanea-align** | NW/SW alignment, MSA, POA, HMMs, WFA, GPU dispatch, CIGAR | [API.md](../cyanea-align/docs/API.md) | [GUIDE.md](../cyanea-align/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-align/docs/ARCHITECTURE.md) |
| **cyanea-io** | 30+ file format parsers (VCF, BAM, Parquet, ...) | [API.md](../cyanea-io/docs/API.md) | [GUIDE.md](../cyanea-io/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-io/docs/ARCHITECTURE.md) |
| **cyanea-omics** | Genomic intervals, variants, single-cell, spatial, AnnData | [API.md](../cyanea-omics/docs/API.md) | [GUIDE.md](../cyanea-omics/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-omics/docs/ARCHITECTURE.md) |
| **cyanea-stats** | Statistics, distributions, popgen, survival, enrichment | [API.md](../cyanea-stats/docs/API.md) | [GUIDE.md](../cyanea-stats/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-stats/docs/ARCHITECTURE.md) |
| **cyanea-ml** | Clustering, trees, GBDT, HMM, PCA/t-SNE/UMAP | [API.md](../cyanea-ml/docs/API.md) | [GUIDE.md](../cyanea-ml/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-ml/docs/ARCHITECTURE.md) |
| **cyanea-chem** | SMILES/SDF, fingerprints, conformers, force fields, reactions | [API.md](../cyanea-chem/docs/API.md) | [GUIDE.md](../cyanea-chem/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-chem/docs/ARCHITECTURE.md) |
| **cyanea-struct** | PDB/mmCIF, DSSP, Kabsch, Ramachandran | [API.md](../cyanea-struct/docs/API.md) | [GUIDE.md](../cyanea-struct/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-struct/docs/ARCHITECTURE.md) |
| **cyanea-phylo** | Newick/NEXUS, UPGMA/NJ, ML, bootstrap, tree search | [API.md](../cyanea-phylo/docs/API.md) | [GUIDE.md](../cyanea-phylo/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-phylo/docs/ARCHITECTURE.md) |
| **cyanea-gpu** | CPU/CUDA/Metal/WebGPU backends | [API.md](../cyanea-gpu/docs/API.md) | [GUIDE.md](../cyanea-gpu/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-gpu/docs/ARCHITECTURE.md) |

### Binding Crates

| Crate | Target | API | Guide | Architecture |
|-------|--------|-----|-------|--------------|
| **cyanea-wasm** | WebAssembly (`@cyanea/bio` on npm, 10 modules) | [API.md](../cyanea-wasm/docs/API.md) | [GUIDE.md](../cyanea-wasm/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-wasm/docs/ARCHITECTURE.md) |
| **cyanea-py** | Python (`pip install cyanea`, 11 submodules) | [API.md](../cyanea-py/docs/API.md) | [GUIDE.md](../cyanea-py/docs/GUIDE.md) | [ARCHITECTURE.md](../cyanea-py/docs/ARCHITECTURE.md) |

## Developer Reference

| Document | Description |
|----------|-------------|
| [CLAUDE.md](../CLAUDE.md) | Workspace conventions, crate map, dependency graph, coding patterns, testing strategy |
| [CI workflow](../.github/workflows/ci.yml) | GitHub Actions: check, test, clippy, fmt, doc, WASM, Python |
| [Bench workflow](../.github/workflows/bench.yml) | Criterion benchmarks on PRs with baseline comparison |

## Documentation Map

```
labs/
+-- README.md                          <- Project overview & quick start
+-- CLAUDE.md                          <- Conventions & patterns for contributors
+-- docs/
|   +-- README.md                      <- This file (documentation index)
|   +-- ARCHITECTURE.md                <- Dependency graph, feature flags, data flow, platforms
|   +-- BUILDING.md                    <- Build guide for all targets
|   +-- BINDINGS.md                    <- Python, WASM, and Elixir binding guide
|   +-- GUIDE.md                       <- Cross-crate usage pipelines with examples
|   +-- BENCHMARKS.md                  <- Competitive benchmarks vs rust-bio
+-- .github/workflows/
|   +-- ci.yml                         <- CI pipeline
|   +-- bench.yml                      <- Benchmark pipeline
+-- cyanea-*/docs/
    +-- API.md                         <- Full API reference
    +-- GUIDE.md                       <- Usage guide with examples
    +-- ARCHITECTURE.md                <- Internal design and module map
```
