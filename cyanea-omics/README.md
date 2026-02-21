# cyanea-omics

> Data structures for genomics, transcriptomics, variant analysis, single-cell biology, and spatial transcriptomics.

434 unit tests + 2 doc tests.

## What's Inside

- **Genomic coordinates** -- `GenomicPosition`, `GenomicInterval`, `Strand` (0-based half-open)
- **Interval operations** -- `IntervalSet` with overlap queries, merging, coverage computation
- **Interval tree** -- augmented BST with O(log n + k) overlap, nearest/preceding/following queries
- **Coverage vectors** -- RLE-encoded genome-wide depth from intervals
- **Expression matrices** -- dense features x samples matrix with named rows/columns
- **Sparse matrices** -- COO-format with CSR conversion for high-dimensional data
- **Variant types** -- VCF-style `Variant` with `VariantType`, `Zygosity`, filters
- **Variant annotation** -- coding consequence prediction (missense/nonsense/frameshift/splice), HGVS notation, SIFT-style scoring
- **Gene annotations** -- `Gene` / `Transcript` / `Exon` hierarchy with biotypes
- **AnnData container** -- single-cell data with obs/var metadata, obsm/varm embeddings, layers, obsp, uns
- **HDF5 I/O** -- read/write `.h5ad` files (dense and CSR sparse, metadata, embeddings)
- **Zarr I/O** -- read/write Zarr v3 directories (pure Rust, same feature set as h5ad)
- **OTU/ASV tables** -- abundance tables with rarefaction, filtering, taxonomic collapse
- **Network biology** -- weighted graphs, degree/betweenness/closeness centrality, Louvain communities
- **Haplotype analysis** -- EM phasing, haplotype block detection, diversity statistics
- **Genome arithmetic** -- intersect, union, subtract, complement, closest, window, Jaccard
- **Liftover** -- UCSC chain file parsing, coordinate remapping between assemblies
- **CNV analysis** -- CBS segmentation, B-allele frequency, SV breakpoint detection, segment merging
- **Methylation** -- bisulfite conversion, CpG site calling, DMRs, CpG island detection
- **Spatial transcriptomics** -- Delaunay/kNN spatial graphs, Moran's I, Geary's C, co-occurrence, ligand-receptor
- **Single-cell pipeline** -- normalize, HVG, kNN graph, Leiden/Louvain, diffusion map, DPT, PAGA, RNA velocity, markers, Harmony/ComBat/MNN integration

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
| `single-cell` | No | Single-cell analysis pipeline (HVG, clustering, trajectory, markers, integration) |

## Modules

| Module | Description |
|--------|-------------|
| `genomic` | `GenomicPosition`, `GenomicInterval`, `Strand` |
| `interval` | `IntervalSet` with overlap, merge, coverage |
| `interval_tree` | Augmented BST with O(log n + k) overlap queries |
| `coverage` | RLE coverage vectors |
| `expr` | Dense expression matrix |
| `sparse` | COO sparse matrix with CSR conversion |
| `variant` | `Variant`, `VariantType`, `Zygosity` |
| `variant_annotation` | Consequence prediction, HGVS notation, splice scoring |
| `annotation` | `Gene`, `Transcript`, `Exon` hierarchy |
| `single_cell` | `AnnData` container with typed metadata |
| `h5ad` | HDF5 `.h5ad` reader/writer (feature-gated) |
| `zarr` | Zarr v3 reader/writer (feature-gated) |
| `otu` | OTU/ASV abundance tables |
| `network` | Weighted graphs, centrality, Louvain |
| `haplotype` | EM phasing, haplotype blocks |
| `genome_arithmetic` | Intersect, union, subtract, complement, closest, window, Jaccard |
| `liftover` | UCSC chain file parsing, coordinate liftover |
| `cnv` | CBS segmentation, BAF, SV breakpoints |
| `methylation` | CpG sites, DMRs, CpG islands, bisulfite conversion |
| `spatial` | Spatial neighbor graphs, Moran's I, Geary's C, co-occurrence, ligand-receptor |
| `sc_preprocess` | HVG, normalize, regress, doublet detection, gene scoring (feature-gated) |
| `sc_cluster` | kNN graph, Leiden, Louvain, NMI, ARI (feature-gated) |
| `sc_trajectory` | Diffusion map, DPT, PAGA, RNA velocity (feature-gated) |
| `sc_markers` | Marker gene detection: t-test, Wilcoxon, logistic regression (feature-gated) |
| `sc_integrate` | Harmony, ComBat, MNN, kBET/LISI metrics (feature-gated) |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Architecture](docs/ARCHITECTURE.md)
