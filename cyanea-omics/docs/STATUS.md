# cyanea-omics

Data structures for genomics, transcriptomics, and variant analysis. Provides the core types used across the Cyanea ecosystem for representing biological coordinates, expression data, variants, and gene annotations.

## Status: Complete

All omics data structures are implemented including genomic coordinates, interval operations, expression matrices, sparse matrices, variant types, gene annotations, and an AnnData-like single-cell container.

## Public API

### Genomic coordinates (`genomic.rs`)

| Type | Description |
|------|-------------|
| `Strand` | Enum: `Plus`, `Minus`, `Unstranded` |
| `GenomicPosition` | `chrom`, `position`, `strand` |
| `GenomicInterval` | `chrom`, `start`, `end`, `strand` (0-based half-open) |

### Interval operations (`interval.rs`)

| Type | Description |
|------|-------------|
| `IntervalSet` | Collection of `GenomicInterval` with overlap queries, merging, and coverage computation |

### Expression matrices (`expr.rs`)

| Type | Description |
|------|-------------|
| `ExpressionMatrix` | Dense features x samples matrix with named rows/columns |

**Methods:**

| Method | Description |
|--------|-------------|
| `new(data, feature_names, sample_names)` | Construct with validation |
| `shape() -> (usize, usize)` | (features, samples) |
| `get(feature_idx, sample_idx) -> Option<f64>` | Element access |
| Row/column selection, normalization | Subsetting and transforms |

### Sparse matrices (`sparse.rs`)

| Type | Description |
|------|-------------|
| `SparseMatrix` | COO-format sparse matrix for high-dimensional omics data |

### Variants (`variant.rs`)

| Type | Description |
|------|-------------|
| `VariantType` | Enum: `SNP`, `Indel`, `Structural`, `Other` |
| `Zygosity` | Enum: `Homozygous`, `Heterozygous`, `Missing` |
| `VariantFilter` | Filter enum for variant QC |
| `Variant` | VCF-style variant: `chrom`, `position`, `ref_allele`, `alt_alleles`, `quality`, `filter`, `info` |

### Gene annotations (`annotation.rs`)

| Type | Description |
|------|-------------|
| `GeneType` | Enum: `Protein`, `Lnc`, `Pseudogene`, `Other` |
| `Exon` | `start`, `end`, `phase` |
| `Transcript` | `id`, `exons`, `biotype` |
| `Gene` | `id`, `symbol`, `transcripts`, `gene_type`, `chrom`, `start`, `end`, `strand` |

### Single-cell container (`single_cell.rs`)

| Type | Description |
|------|-------------|
| `MatrixData` | Enum: `Dense(Vec<Vec<f64>>)`, `Sparse(SparseMatrix)` |
| `QcMetrics` | `total_counts`, `n_features` per observation |
| `AnnData` | AnnData-like container with obs, var, X, layers, obsm, varm |

**AnnData methods:**

| Method | Description |
|--------|-------------|
| `new(x, obs_names, var_names) -> Result<Self>` | Construct with dimension validation |
| `n_obs() -> usize` | Number of observations |
| `n_vars() -> usize` | Number of variables |
| `shape() -> (usize, usize)` | (observations, variables) |
| `x() -> &MatrixData` | Access expression matrix |
| `obs_names() / var_names()` | Access observation/variable names |
| `add_obs(key, values) / get_obs(key)` | Per-observation annotations |
| `add_var(key, values) / get_var(key)` | Per-variable annotations |
| `add_obsm(key, data) / get_obsm(key)` | Observation-level embeddings |
| `add_varm(key, data) / get_varm(key)` | Variable-level embeddings |
| `add_layer(key, layer) / get_layer(key)` | Alternative matrix layers |
| `subset_obs(indices) -> Result<AnnData>` | Subset by observation indices |
| `qc_metrics() -> QcMetrics` | Compute QC metrics (total counts, features per obs) |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |

## Dependencies

- `cyanea-core` -- error types, traits

## Tests

87 unit tests + 2 doc tests across 8 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 42 | Module declarations, re-exports |
| `genomic.rs` | 291 | Coordinates, strand, position, interval |
| `interval.rs` | 255 | IntervalSet with overlap/merge/coverage |
| `expr.rs` | 411 | Dense expression matrix |
| `sparse.rs` | 333 | COO sparse matrix |
| `variant.rs` | 310 | Variant types and filtering |
| `annotation.rs` | 297 | Gene/Transcript/Exon hierarchy |
| `single_cell.rs` | 513 | AnnData-like single-cell container |
