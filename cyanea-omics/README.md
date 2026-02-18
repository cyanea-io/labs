# cyanea-omics

> Data structures for genomics, transcriptomics, and variant analysis.

## What's Inside

- **Genomic coordinates** -- `GenomicPosition`, `GenomicInterval`, `Strand` (0-based half-open)
- **Interval operations** -- `IntervalSet` with overlap queries, merging, coverage computation
- **Expression matrices** -- dense features x samples matrix with named rows/columns
- **Sparse matrices** -- COO-format with CSR conversion for high-dimensional data
- **Variant types** -- VCF-style `Variant` with `VariantType`, `Zygosity`, filters
- **Gene annotations** -- `Gene` / `Transcript` / `Exon` hierarchy with biotypes
- **AnnData container** -- single-cell data with obs/var metadata, obsm/varm embeddings, layers
- **HDF5 I/O** -- read/write `.h5ad` files (dense and CSR sparse, metadata, embeddings)
- **Zarr I/O** -- read/write Zarr v3 directories (pure Rust, same feature set as h5ad)
- **OTU/ASV tables** -- abundance tables with rarefaction, filtering, taxonomic collapse
- **Network biology** -- weighted graphs, degree/betweenness/closeness centrality, Louvain communities
- **Haplotype analysis** -- EM phasing, haplotype block detection, diversity statistics

## Quick Start

```toml
[dependencies]
cyanea-omics = { version = "0.1", features = ["zarr"] }
```

```rust
use cyanea_omics::{GenomicInterval, IntervalSet, Strand};

let intervals = vec![
    GenomicInterval::new("chr1", 100, 200, Strand::Plus),
    GenomicInterval::new("chr1", 150, 300, Strand::Plus),
];
let set = IntervalSet::from(intervals);
let merged = set.merge();
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `h5ad` | No | HDF5 `.h5ad` I/O (requires system HDF5 1.10.x) |
| `zarr` | No | Zarr v3 directory I/O (pure Rust) |

## Modules

| Module | Description |
|--------|-------------|
| `genomic` | `GenomicPosition`, `GenomicInterval`, `Strand` |
| `interval` | `IntervalSet` with overlap, merge, coverage |
| `expr` | Dense expression matrix |
| `sparse` | COO sparse matrix with CSR conversion |
| `variant` | `Variant`, `VariantType`, `Zygosity` |
| `annotation` | `Gene`, `Transcript`, `Exon` hierarchy |
| `single_cell` | `AnnData` container with typed metadata |
| `h5ad` | HDF5 `.h5ad` reader/writer (feature-gated) |
| `zarr` | Zarr v3 reader/writer (feature-gated) |
| `otu` | OTU/ASV abundance tables |
| `network` | Weighted graphs, centrality, Louvain |
| `haplotype` | EM phasing, haplotype blocks |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
