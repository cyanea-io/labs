# Architecture

Cyanea Labs is a Rust bioinformatics ecosystem comprising 13 crates in a Cargo workspace. It provides sequence analysis, alignment, statistics, machine learning, cheminformatics, structural biology, phylogenetics, and GPU acceleration -- all from a single dependency tree rooted in `cyanea-core`. The workspace compiles to native Rust, WebAssembly, Python (via PyO3), and Elixir NIFs.

## Dependency Graph

```
                          ┌─────────────┐
                          │ cyanea-core │
                          │             │
                          │  errors     │
                          │  traits     │
                          │  hashing    │
                          │  compress   │
                          │  probability│
                          │  bitvectors │
                          │  fenwick    │
                          └──────┬──────┘
                                 │
        ┌────────┬────────┬──────┼──────┬────────┬────────┬────────┐
        │        │        │      │      │        │        │        │
        ▼        ▼        ▼      ▼      ▼        ▼        ▼        ▼
     seq      align    omics   stats    ml     chem    struct   phylo
      │                  │                                │      │
      │                  │                                │      │
      │                  ├───── io ◄───────────(omics)────┘      │
      │                  │                                       │
      │                  │                          (optional)   │
      │                  │                        ┌──────────────┘
      │                  │                        ▼
      │                  │                    cyanea-ml
      │                  │
      │                  │
      ▼                  ▼
    ┌───────────────────────────────────┐
    │         Binding Crates           │
    ├───────────┬───────────┬──────────┤
    │  wasm     │    py     │  NIF*    │
    │           │           │          │
    │  8 domain │  9 domain │  path    │
    │  crates   │  crates   │  deps    │
    └───────────┴───────────┴──────────┘
                                * lives in ../cyanea/native/
```

**Key relationships:**

- Every domain crate depends only on `cyanea-core` (no cross-domain deps)
- `cyanea-io` depends on `cyanea-omics` for genomic types (`Variant`, `GenomicInterval`, `Gene`)
- `cyanea-phylo` optionally depends on `cyanea-ml` (behind `ml` feature) for distance matrices and tree construction
- `cyanea-io` optionally depends on `cyanea-stats` (behind `variant-calling` feature)
- Binding crates (`wasm`, `py`, NIF) aggregate multiple domain crates

## Feature Flag Architecture

Feature flags control optional dependencies and platform-specific code. Flags propagate through the dependency tree:

```
┌─────────────────────────────────────────────────────────────────────┐
│                       Workspace Feature Flags                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  std (default)                                                      │
│  ├── Enables: memmap2, zstd, flate2, file I/O                     │
│  └── Used by: core, seq, align, omics, stats, ml, chem, struct,   │
│               phylo                                                 │
│                                                                     │
│  serde                                                              │
│  ├── Enables: serde derives on public types                        │
│  └── Propagates: child crate serde → core/serde                   │
│                                                                     │
│  parallel                                                           │
│  ├── Enables: rayon parallelism                                    │
│  └── Used by: align, stats, ml, chem, struct, phylo, io           │
│                                                                     │
│  simd (default in align)                                            │
│  ├── Enables: NEON (aarch64), SSE4.1/AVX2 (x86_64)               │
│  └── Used by: align                                                │
│                                                                     │
│  GPU backends                                                       │
│  ├── cuda  → cudarc + NVRTC          (gpu, align)                  │
│  ├── metal → metal-rs                (gpu, align)                  │
│  └── wgpu  → wgpu + bytemuck        (gpu)                         │
│                                                                     │
│  Format flags (cyanea-io)                                           │
│  ├── csv (default) ─── csv, serde, serde_json                     │
│  ├── vcf ─── cyanea-omics                                         │
│  ├── bed ─── cyanea-omics                                         │
│  ├── gff ─── cyanea-omics                                         │
│  ├── gtf ─── cyanea-omics                                         │
│  ├── sam                                                           │
│  ├── bam ─── implies sam, adds flate2                              │
│  ├── cram ── implies sam, adds noodles                             │
│  ├── bcf ─── implies vcf, adds flate2                              │
│  ├── parquet ─ implies vcf + bed, adds arrow/parquet               │
│  ├── blast                                                         │
│  ├── maf                                                           │
│  ├── genbank                                                       │
│  └── bigwig ─ adds flate2                                          │
│                                                                     │
│  Storage backends (cyanea-omics)                                    │
│  ├── h5ad ── HDF5 .h5ad I/O (requires system HDF5)                │
│  └── zarr ── Zarr v3 directory I/O (pure Rust)                     │
│                                                                     │
│  minhash (cyanea-seq) ── MinHash/FracMinHash sketching             │
│  wfa (cyanea-align) ── Wavefront alignment                         │
│  ml (cyanea-phylo) ── cyanea-ml for tree construction              │
│  blas (stats, ml) ── ndarray-backed PCA                            │
│  numpy (cyanea-py) ── NumPy array interop                          │
│  wasm (cyanea-wasm) ── wasm-bindgen annotations                    │
│                                                                     │
│  Implication chains:                                                │
│    bam → sam                                                        │
│    cram → sam                                                       │
│    bcf → vcf                                                        │
│    parquet → vcf + bed                                              │
│    variant-calling → sam + vcf                                      │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

## Data Flow

### Raw reads to variant calling

```
  FASTQ files
      │
      ▼
  cyanea-seq          Parse, validate, quality scores
  ├── TrimPipeline    Adapter removal, quality trimming
  ├── PairedFastq     Mate-pair validation
  │
  ▼
  cyanea-align        Align reads to reference
  ├── banded_sw()     Banded local alignment
  ├── seed_extend()   Seed-and-extend for longer reads
  │
  ▼
  cyanea-io           Generate pileup, call variants
  ├── parse_sam()     Read alignments
  ├── pileup          Per-position coverage
  │
  ▼
  cyanea-omics        Genomic coordinates, intervals
  ├── Variant         VCF-style variant records
  └── IntervalSet     Region queries
```

### Sequence comparison pipeline

```
  Sequences (FASTA)
      │
      ▼
  cyanea-seq
  ├── kmers()            Extract k-mers
  ├── MinHash            Bottom-k sketching
  │
  ▼
  cyanea-ml
  ├── pairwise_distances()   Distance matrix
  ├── kmeans() / dbscan()    Clustering
  ├── umap() / tsne()        Visualization
  │
  ▼
  cyanea-phylo
  ├── neighbor_joining()     Tree construction
  └── bootstrap_support()    Statistical support
```

### Columnar analytics

```
  VCF / BED files
      │
      ▼
  cyanea-io
  ├── parse_vcf()                  Raw parsing
  ├── write_variants_parquet()     Columnar conversion
  │
  ▼
  cyanea-io (parquet)
  ├── read_variants_parquet_region()   Predicate pushdown
  │
  ▼
  cyanea-stats
  ├── popgen::fst()                Population genetics
  ├── popgen::tajimas_d()          Neutrality tests
  └── enrichment::gsea_preranked() Gene set enrichment
```

## Error Handling

All crates use `CyaneaError` from `cyanea-core`, built with `thiserror` 2.x:

```
CyaneaError
├── Io(#[from] std::io::Error)
├── Parse(String)
├── InvalidInput(String)
├── Compression(String)
├── Hash(String)
└── Other(String)
```

**Propagation at binding boundaries:**

```
                    Rust domain crates
                    Result<T, CyaneaError>
                           │
              ┌────────────┼────────────┐
              ▼            ▼            ▼
         cyanea-wasm   cyanea-py   cyanea-native
              │            │            │
              ▼            ▼            ▼
         JSON envelope  PyErr       ErlNifError
         {"ok": val}    IOError     {:error, msg}
         {"error": msg} ValueError
                        RuntimeError
```

- **WASM**: `wasm_ok(val)` / `wasm_err(msg)` / `wasm_result(r)` produce JSON strings
- **Python**: `IntoPyResult` trait converts `CyaneaError` variants to appropriate Python exceptions
- **NIF**: Errors become `{:error, reason}` tuples

## Testing Strategy

| Layer | Tool | What it tests |
|-------|------|---------------|
| Unit tests | `#[cfg(test)]` | Every module, inline at file bottom |
| Property tests | `proptest` | Round-trip encoding, invariant verification |
| Fuzz targets | `cargo-fuzz` | All parsers (SMILES, SDF, PDB, Newick, NEXUS, VCF, BED, GFF3) |
| Benchmarks | Criterion | Hot paths in align, seq, ml, stats, chem, struct, gpu |
| CI | GitHub Actions | check, test, clippy, fmt, doc, WASM, Python |

**Test counts by crate:**

| Crate | Tests |
|-------|------:|
| cyanea-core | 58 |
| cyanea-seq | 458 |
| cyanea-io | 158 |
| cyanea-align | 290 |
| cyanea-omics | 117 |
| cyanea-stats | 334 |
| cyanea-ml | 269 |
| cyanea-chem | 79 |
| cyanea-struct | 76 |
| cyanea-phylo | 134 |
| cyanea-gpu | 62+ |
| cyanea-wasm | 111 |
| **Total** | **2,100+** |

Test data is always inline (strings, vecs) -- no external fixture files. Tests needing the filesystem use the `tempfile` crate.

## Platform Support Matrix

| Crate | Native | WASM | Python | NIF | CUDA | Metal | WebGPU |
|-------|:------:|:----:|:------:|:---:|:----:|:-----:|:------:|
| cyanea-core | x | x | x | x | | | |
| cyanea-seq | x | x | x | x | | | |
| cyanea-io | x | x | x | x | | | |
| cyanea-align | x | x | x | x | x | x | |
| cyanea-omics | x | | x | x | | | |
| cyanea-stats | x | x | x | x | | | |
| cyanea-ml | x | x | x | x | | | |
| cyanea-chem | x | x | x | x | | | |
| cyanea-struct | x | x | x | x | | | |
| cyanea-phylo | x | x | x | x | | | |
| cyanea-gpu | x | | | | x | x | x |

## Crate Roles

| Crate | Role | Key Abstractions |
|-------|------|-----------------|
| **cyanea-core** | Foundation | `CyaneaError`, traits, SHA-256, compression, `LogProb`/`PhredProb`, bitvectors, Fenwick tree |
| **cyanea-seq** | Sequence analysis | `DnaSequence`/`RnaSequence`/`ProteinSequence`, FASTA/FASTQ, FM-index, MinHash, pattern matching, trimming |
| **cyanea-align** | Alignment | NW/SW/semi-global, banded, MSA, POA, pair/profile HMM, GPU dispatch, CIGAR |
| **cyanea-io** | File formats | 15 parsers (CSV through Parquet), feature-gated, streaming stats |
| **cyanea-omics** | Genomic data | Intervals, variants, AnnData, h5ad/zarr, OTU tables, networks, haplotypes |
| **cyanea-stats** | Statistics | Descriptive, hypothesis tests, distributions, PCA, Bayesian, popgen, survival |
| **cyanea-ml** | Machine learning | Clustering, KNN, trees/forests/GBDT, HMM, PCA/t-SNE/UMAP, metrics, CV |
| **cyanea-chem** | Cheminformatics | SMILES/SDF, Morgan/MACCS fingerprints, properties, substructure, stereochemistry |
| **cyanea-struct** | Structural biology | PDB/mmCIF, DSSP, Kabsch, contact maps, Ramachandran, B-factors |
| **cyanea-phylo** | Phylogenetics | Newick/NEXUS, UPGMA/NJ, parsimony, ML, bootstrap, consensus, dating, drawing |
| **cyanea-gpu** | GPU compute | Backend trait (CPU/CUDA/Metal/WebGPU), pairwise distances, k-mer counting, SW, MinHash |
| **cyanea-wasm** | WASM bindings | JSON API for 8 domain crates, `@cyanea/bio` npm package |
| **cyanea-py** | Python bindings | PyO3 classes for 9 modules, NumPy interop, `pip install cyanea` |
