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
| `n_rows() / n_cols()` | Dimension accessors |
| `row_nnz(row) / col_nnz(col)` | Per-row/column nonzero counts |
| `column_sums() / column_means()` | Per-column aggregation |
| `row_sums()` | Per-row sums |
| `scale_rows(factors)` | Multiply rows by per-row factors |
| `map_values(f)` | Apply function to all stored values |
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
| `AnnData` | AnnData-like container with obs, var, X, layers, obsm, varm, obsp, uns |

**AnnData methods:**

| Method | Description |
|--------|-------------|
| `new(x, obs_names, var_names) -> Result<Self>` | Construct with dimension validation |
| `n_obs() -> usize` | Number of observations |
| `n_vars() -> usize` | Number of variables |
| `shape() -> (usize, usize)` | (observations, variables) |
| `x() -> &MatrixData` | Access expression matrix |
| `x_mut() -> &mut MatrixData` | Mutable access to expression matrix |
| `set_x(new_x) -> Result<()>` | Replace expression matrix (same shape) |
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
| `subset_vars(indices) -> Result<AnnData>` | Subset by variable indices |
| `add_obsp(key, SparseMatrix) / get_obsp(key)` | Pairwise observation annotations (e.g. kNN graphs) |
| `add_uns(key, String) / get_uns(key)` | Unstructured metadata |
| `get_layer_mut(key) -> Option<&mut MatrixData>` | Mutable layer access |
| `qc_metrics() -> QcMetrics` | Compute QC metrics (total counts, features per obs) |

### OTU/ASV tables (`otu.rs`)

OTU/ASV abundance table operations for metagenomics.

| Type | Description |
|------|-------------|
| `OtuTable` | Samples × OTUs abundance table |

**OtuTable methods:**

| Method | Description |
|--------|-------------|
| `new(counts, sample_names, otu_names) -> Result<Self>` | Create with dimension validation |
| `n_samples() / n_otus()` | Dimension accessors |
| `total_counts() -> Vec<usize>` | Per-sample totals |
| `relative_abundance() -> Vec<Vec<f64>>` | Normalized per sample |
| `filter_min_count(min_count) -> Self` | Keep OTUs with total ≥ threshold |
| `filter_min_prevalence(min_fraction) -> Self` | Keep OTUs present in ≥ fraction of samples |
| `rarefy(depth) -> Result<Self>` | Subsample to fixed depth |
| `collapse_taxonomy(level, taxonomy) -> Result<Self>` | Collapse by taxonomy prefix |
| `merge(other) -> Result<Self>` | Merge two tables with same OTUs |

### Network biology (`network.rs`)

Weighted graphs, centrality metrics, and community detection.

| Type | Description |
|------|-------------|
| `Graph` | Weighted directed or undirected graph |
| `CentralityScores` | Degree, betweenness, closeness centrality for all nodes |
| `Community` | Community assignments with modularity score |

**Graph methods:**

| Method | Description |
|--------|-------------|
| `new(n_nodes, directed) -> Self` | Create empty graph |
| `add_edge(from, to, weight) -> Result<()>` | Add a weighted edge |
| `from_correlation_matrix(matrix, threshold) -> Self` | Build from correlation matrix |
| `from_sparse_matrix(matrix) -> Self` | Build from symmetric sparse matrix |
| `n_nodes() / n_edges()` | Graph size |
| `neighbors(node) -> &[(usize, f64)]` | Neighbor list with weights |
| `degree_centrality() -> Vec<f64>` | Degree centrality |
| `betweenness_centrality() -> Vec<f64>` | Brandes' algorithm |
| `closeness_centrality() -> Vec<f64>` | BFS-based closeness |
| `centrality() -> CentralityScores` | All three centrality metrics |
| `louvain() -> Community` | Louvain community detection |
| `louvain_with_resolution(resolution) -> Community` | Louvain with resolution parameter (Phase 1+2) |
| `modularity(assignments) -> f64` | Compute modularity Q |
| `modularity_with_resolution(assignments, resolution) -> f64` | Modularity Q with resolution γ |

### Haplotype analysis (`haplotype.rs`)

EM phasing, haplotype block detection, and diversity statistics.

| Type | Description |
|------|-------------|
| `Haplotype` | Sequence of biallelic alleles (0/1) |
| `PhasedGenotypes` | EM phasing result with haplotype pairs and frequencies |
| `HaplotypeBlock` | Contiguous SNP block in high LD |

| Function | Description |
|----------|-------------|
| `phase_em(genotypes, max_iter) -> Result<PhasedGenotypes>` | EM-based haplotype phasing |
| `haplotype_blocks(genotypes, threshold) -> Result<Vec<HaplotypeBlock>>` | Detect blocks from LD structure |
| `haplotype_diversity(haplotypes) -> f64` | Haplotype diversity (expected heterozygosity analog) |

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

### Zarr reader/writer (`zarr.rs`, feature: `zarr`)

| Function | Description |
|----------|-------------|
| `read_zarr(path) -> Result<AnnData>` | Read Zarr directory into AnnData container |
| `write_zarr(adata, path) -> Result<()>` | Write AnnData container to Zarr directory |

**Supported Zarr components** (mirrors h5ad):

- `X` — dense arrays and CSR sparse matrices
- `obs` / `var` — string, numeric, and categorical metadata columns
- `obs_names` / `var_names` — observation and variable index names
- `obsm` / `varm` — 2D embedding matrices (e.g. PCA, UMAP)
- `layers` — alternative dense or sparse matrices

Pure Rust implementation using zarrs 0.18 (Zarr v3) — no system library requirements.

### Single-cell preprocessing (`sc_preprocess.rs`, feature: `single-cell`)

| Function | Description |
|----------|-------------|
| `highly_variable_genes(adata, config)` | HVG selection via Seurat v3 (VST) or CellRanger (mean/dispersion binning) |
| `normalize_total(adata, config)` | Per-cell library-size normalization + optional log1p |
| `regress_out(adata, keys)` | OLS regression per gene on confounders from obs |
| `scrublet_doublets(adata, config)` | Doublet detection: simulate doublets, PCA, kNN scoring |
| `score_genes(adata, gene_list, n_reference, name)` | Per-cell signature scoring (mean signature - mean reference) |

| Config type | Key fields |
|-------------|------------|
| `HvgConfig` | `n_top_genes`, `method: HvgMethod`, `min_mean`, `max_mean`, `min_disp`, `n_bins` |
| `NormalizeConfig` | `target_sum`, `log_transform`, `save_raw` |
| `ScrubletConfig` | `expected_doublet_rate`, `sim_doublet_ratio`, `n_pcs`, `k_neighbors`, `seed` |

### Single-cell clustering (`sc_cluster.rs`, feature: `single-cell`)

| Function | Description |
|----------|-------------|
| `neighbors(adata, config)` | Build kNN graph from PCA, UMAP-style fuzzy connectivities |
| `leiden(adata, config)` | Leiden algorithm (Traag 2019): local moves + refinement + aggregation |
| `louvain(adata, config)` | Louvain wrapper with resolution parameter |
| `nmi(a, b)` | Normalized Mutual Information between partitions |
| `adjusted_rand_index(a, b)` | Adjusted Rand Index between partitions |

| Config type | Key fields |
|-------------|------------|
| `NeighborsConfig` | `n_neighbors`, `n_pcs`, `metric: DistanceMetric`, `seed` |
| `ClusterConfig` | `resolution`, `n_iterations`, `seed`, `key_added` |

### Single-cell trajectory (`sc_trajectory.rs`, feature: `single-cell`)

| Function | Description |
|----------|-------------|
| `diffusion_map(adata, config) -> DiffusionResult` | Anisotropic diffusion map via power iteration |
| `dpt(adata, config)` | Diffusion pseudotime from root cell |
| `paga(adata, cluster_key) -> PagaResult` | Partition-based graph abstraction |
| `rna_velocity(adata, config)` | Steady-state RNA velocity (γ by regression) |

| Config/Result type | Key fields |
|---------------------|------------|
| `DiffusionConfig` | `n_components`, `alpha` |
| `DiffusionResult` | `eigenvalues`, `components` |
| `DptConfig` | `root_cell`, `n_branchings` |
| `PagaResult` | `connectivities`, `groups`, `cluster_sizes` |
| `VelocityConfig` | `mode: VelocityMode`, `n_neighbors`, `min_counts` |

### Single-cell markers (`sc_markers.rs`, feature: `single-cell`)

| Function | Description |
|----------|-------------|
| `rank_genes_groups(adata, config) -> MarkerResults` | One-vs-rest DE per cluster (t-test / Wilcoxon / logistic regression) |
| `filter_markers(results, config) -> MarkerResults` | Filter by log2FC, pct, padj thresholds |

| Config/Result type | Key fields |
|---------------------|------------|
| `MarkerConfig` | `method: MarkerMethod`, `cluster_key`, `log2fc_threshold`, `min_pct`, `padj_threshold`, `n_genes` |
| `MarkerGene` | `gene_name`, `gene_index`, `log2_fold_change`, `pct_in`, `pct_out`, `statistic`, `p_value`, `p_adjusted` |
| `MarkerResults` | `markers: HashMap<String, Vec<MarkerGene>>`, `method`, `n_clusters` |

### Single-cell integration (`sc_integrate.rs`, feature: `single-cell`)

| Function | Description |
|----------|-------------|
| `harmony(adata, config)` | Harmony iterative PCA correction (Korsunsky 2019) |
| `combat(adata, config)` | ComBat parametric empirical Bayes batch correction (Johnson 2007) |
| `mnn_correct(adata, config)` | Mutual Nearest Neighbors correction (Haghverdi 2018) |
| `integration_metrics(adata, config) -> IntegrationMetrics` | kBET + LISI integration quality metrics |

| Config/Result type | Key fields |
|---------------------|------------|
| `HarmonyConfig` | `batch_key`, `n_clusters`, `theta`, `sigma`, `max_iter` |
| `CombatConfig` | `batch_key`, `parametric` |
| `MnnConfig` | `batch_key`, `k`, `sigma`, `cos_norm` |
| `IntegrationMetrics` | `kbet_accept_rate`, `mean_ilisi`, `mean_clisi` |
| `MetricsConfig` | `batch_key`, `label_key`, `n_neighbors`, `alpha` |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `h5ad` | No | HDF5-backed `.h5ad` file I/O (requires system HDF5 library) |
| `zarr` | No | Zarr v3 directory-based file I/O (pure Rust) |
| `single-cell` | No | Single-cell analysis pipeline (HVG, clustering, trajectory, markers, integration) |

## Dependencies

- `cyanea-core` — error types, traits
- `cyanea-stats` (optional, `single-cell` feature) — statistical tests, BH correction
- `cyanea-ml` (optional, `single-cell` feature) — PCA, k-means, distance metrics
- `hdf5` 0.8 (optional, `h5ad` feature) — HDF5 bindings
- `ndarray` 0.15 (optional, `h5ad` feature) — array types for HDF5 interop
- `zarrs` 0.18 (optional, `zarr` feature) — Zarr v3 storage
- `serde_json` (optional, `zarr` feature) — metadata serialization

## Tests

245 unit tests + 2 doc tests (base); 366 + 2 with `single-cell` feature (151 new T10 tests).

## Source Files

| File | Purpose |
|------|---------|
| `lib.rs` | Module declarations, re-exports |
| `genomic.rs` | Coordinates, strand, position, interval |
| `interval.rs` | IntervalSet with overlap/merge/coverage |
| `expr.rs` | Dense expression matrix |
| `sparse.rs` | COO sparse matrix with CSR conversion, column/row ops |
| `variant.rs` | Variant types and filtering |
| `annotation.rs` | Gene/Transcript/Exon hierarchy |
| `single_cell.rs` | AnnData container with typed metadata, obsp, uns |
| `network.rs` | Graphs, centrality, Louvain with resolution |
| `h5ad.rs` | HDF5 `.h5ad` reader/writer |
| `zarr.rs` | Zarr v3 directory reader/writer |
| `sc_preprocess.rs` | HVG, normalize, regress, scrublet, score_genes |
| `sc_cluster.rs` | kNN graph, Leiden, Louvain, NMI, ARI |
| `sc_trajectory.rs` | Diffusion map, DPT, PAGA, RNA velocity |
| `sc_markers.rs` | Marker genes (t-test, Wilcoxon, logistic regression) |
| `sc_integrate.rs` | Harmony, ComBat, MNN, kBET/LISI metrics |
