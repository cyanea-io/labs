# cyanea-omics

Data structures for genomics, transcriptomics, and variant analysis. Provides the core types used across the Cyanea ecosystem for representing biological coordinates, expression data, variants, and gene annotations.

## Status: Complete

All omics data structures are implemented including genomic coordinates, interval operations, expression matrices, sparse matrices, variant types, gene annotations, an AnnData-like single-cell container, and HDF5-backed `.h5ad` file I/O.

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

**Methods:**

| Method | Description |
|--------|-------------|
| `new(n_rows, n_cols)` | Empty matrix with dimensions |
| `from_triplets(rows, cols, values, n_rows, n_cols)` | Construct from COO triplets |
| `from_dense(data, threshold)` | Convert from dense, filtering by threshold |
| `from_csr(data, indices, indptr, n_rows, n_cols)` | Construct from CSR format |
| `insert(row, col, value)` | Insert a single entry |
| `get(row, col) -> f64` | Element access (O(nnz) scan) |
| `to_dense() -> Vec<Vec<f64>>` | Convert to dense 2D vector |
| `to_csr() -> (Vec<f64>, Vec<usize>, Vec<usize>)` | Convert to CSR `(data, indices, indptr)` |
| `nnz() / density() / shape()` | Matrix statistics |
| `row_nnz(row) / col_nnz(col)` | Per-row/column nonzero counts |
| `iter()` | Iterate over `(row, col, value)` triplets |

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

### Metadata columns (`single_cell.rs`)

| Type | Description |
|------|-------------|
| `ColumnData` | Enum: `Strings(Vec<String>)`, `Numeric(Vec<f64>)`, `Categorical { codes, categories }` |

**ColumnData methods:**

| Method | Description |
|--------|-------------|
| `len() / is_empty()` | Column length |
| `as_strings() -> Option<&Vec<String>>` | Try to access as string data |
| `as_numeric() -> Option<&Vec<f64>>` | Try to access as numeric data |

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
| `add_obs(key, Vec<String>)` | Add string observation annotation |
| `add_obs_numeric(key, Vec<f64>)` | Add numeric observation annotation |
| `add_obs_column(key, ColumnData)` | Add typed observation annotation |
| `get_obs(key) -> Option<&ColumnData>` | Get observation annotation (typed) |
| `get_obs_strings(key) -> Option<&Vec<String>>` | Get observation annotation (string convenience) |
| `obs_columns() -> &HashMap<String, ColumnData>` | All observation annotations |
| `add_var(key, Vec<String>)` | Add string variable annotation |
| `add_var_numeric(key, Vec<f64>)` | Add numeric variable annotation |
| `add_var_column(key, ColumnData)` | Add typed variable annotation |
| `get_var(key) -> Option<&ColumnData>` | Get variable annotation (typed) |
| `get_var_strings(key) -> Option<&Vec<String>>` | Get variable annotation (string convenience) |
| `var_columns() -> &HashMap<String, ColumnData>` | All variable annotations |
| `add_obsm(key, data) / get_obsm(key)` | Observation-level embeddings |
| `add_varm(key, data) / get_varm(key)` | Variable-level embeddings |
| `obsm_keys() / varm_keys()` | All embedding maps |
| `add_layer(key, layer) / get_layer(key)` | Alternative matrix layers |
| `layers_keys()` | All layer maps |
| `subset_obs(indices) -> Result<AnnData>` | Subset by observation indices |
| `qc_metrics() -> QcMetrics` | Compute QC metrics (total counts, features per obs) |

### HDF5 reader/writer (`h5ad.rs`, feature: `h5ad`)

| Function | Description |
|----------|-------------|
| `read_h5ad(path) -> Result<AnnData>` | Read `.h5ad` file into AnnData container |
| `write_h5ad(adata, path) -> Result<()>` | Write AnnData container to `.h5ad` file |

**Supported `.h5ad` components:**

- `X` — dense arrays and CSR sparse matrices
- `obs` / `var` — string, numeric, and categorical metadata columns
- `obs_names` / `var_names` — observation and variable index names
- `obsm` / `varm` — 2D embedding matrices (e.g. PCA, UMAP)
- `layers` — alternative dense or sparse matrices

**System requirement:** HDF5 1.10.x library (`brew install hdf5@1.10` on macOS, `apt install libhdf5-dev` on Linux). Build with `HDF5_DIR="$(brew --prefix hdf5@1.10)"`. Note: hdf5-sys 0.8 does not support HDF5 2.0.

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `h5ad` | No | HDF5-backed `.h5ad` file I/O (requires system HDF5 library) |

## Dependencies

- `cyanea-core` — error types, traits
- `hdf5` 0.8 (optional, `h5ad` feature) — HDF5 bindings
- `ndarray` 0.15 (optional, `h5ad` feature) — array types for HDF5 interop

## Tests

99 unit tests + 2 doc tests across 9 source files (with `h5ad` feature; 90 + 2 without).

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 47 | Module declarations, re-exports |
| `genomic.rs` | 291 | Coordinates, strand, position, interval |
| `interval.rs` | 255 | IntervalSet with overlap/merge/coverage |
| `expr.rs` | 411 | Dense expression matrix |
| `sparse.rs` | 481 | COO sparse matrix with CSR conversion |
| `variant.rs` | 310 | Variant types and filtering |
| `annotation.rs` | 297 | Gene/Transcript/Exon hierarchy |
| `single_cell.rs` | 631 | AnnData container with typed metadata |
| `h5ad.rs` | 735 | HDF5 `.h5ad` reader/writer |
