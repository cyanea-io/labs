# cyanea-ml

Machine learning primitives for bioinformatics: clustering, distance metrics, sequence encoding, normalization, and evaluation.

## Status: Mostly Complete

Core ML primitives are implemented (clustering, distances, encoding, evaluation, normalization, k-mer counting). Neural embeddings, ONNX inference, and dimensionality reduction are stubbed.

## Public API

### Clustering (`cluster.rs`)

**K-Means:**

| Type | Description |
|------|-------------|
| `KMeansConfig` | `k`, `max_iterations`, `tolerance` |
| `KMeansResult` | `centroids`, `labels`, `iterations_run` |
| `kmeans(data, config) -> Result<KMeansResult>` | K-means clustering |

**DBSCAN:**

| Type | Description |
|------|-------------|
| `DbscanConfig` | `eps`, `min_samples` |
| `DbscanResult` | `labels`, `n_clusters`, `n_noise` |
| `dbscan(data, config) -> Result<DbscanResult>` | Density-based clustering |

**Hierarchical:**

| Type | Description |
|------|-------------|
| `Linkage` | Enum: `Single`, `Complete`, `Average`, `Ward` |
| `HierarchicalConfig` | `linkage`, `num_clusters` |
| `MergeStep` | Record of a hierarchical merge operation |
| `HierarchicalResult` | `dendrograms`, `labels`, `merge_steps` |
| `hierarchical(data, config) -> Result<HierarchicalResult>` | Agglomerative clustering |

### Distance metrics (`distance.rs`)

| Type/Function | Description |
|---------------|-------------|
| `DistanceMetric` | Enum: `Euclidean`, `Manhattan`, `Cosine`, `Hamming` |
| `euclidean(a, b) -> Result<f64>` | Euclidean distance |
| `manhattan(a, b) -> Result<f64>` | Manhattan distance |
| `cosine_similarity(a, b) -> Result<f64>` | Cosine similarity (0-1) |
| `cosine_distance(a, b) -> Result<f64>` | 1 - cosine similarity |
| `hamming(a, b) -> Result<usize>` | Byte-level Hamming distance |
| `compute_distance(a, b, metric) -> Result<f64>` | Dispatch by metric |
| `DistanceMatrix` | Precomputed pairwise distance matrix |
| `pairwise_distances(data, metric) -> Result<DistanceMatrix>` | Compute full distance matrix |

### Sequence encoding (`encoding.rs`)

| Function | Description |
|----------|-------------|
| `Alphabet` | Enum: `Dna`, `Rna`, `Protein` |
| `one_hot_encode(seq, alphabet) -> Vec<f64>` | One-hot encoding (n_bases x alphabet_size) |
| `label_encode(seq, alphabet) -> Vec<f64>` | Integer label encoding |

### Evaluation (`evaluate.rs`)

| Function | Description |
|----------|-------------|
| `silhouette_samples(data, labels) -> Result<Vec<f64>>` | Per-sample silhouette coefficients |
| `silhouette_score(data, labels) -> Result<f64>` | Mean silhouette score |

### K-mer counting (`kmer.rs`)

| Type | Description |
|------|-------------|
| `KmerCounter` | K-mer counter for DNA sequences |
| `KmerCounts` | Results with `total()`, `top_n(n)`, `to_frequency_vector(alphabet)` |

### Normalization (`normalize.rs`)

| Function | Description |
|----------|-------------|
| `min_max(data)` | In-place [0, 1] normalization |
| `z_score(data)` | In-place standardization |
| `l2_normalize(data)` | In-place unit norm |
| `min_max_columns(data, n_cols)` | Column-wise min-max |
| `z_score_columns(data, n_cols)` | Column-wise z-score |
| `l2_normalize_columns(data, n_cols)` | Column-wise L2 |

### Planned (stubbed)

| Module | Description |
|--------|-------------|
| `embedding` | Neural embeddings (ESM, ProtTrans, DNA/RNA language models) |
| `inference` | ONNX model inference |
| `reduction` | Dimensionality reduction (PCA, t-SNE, UMAP) |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |

## Dependencies

- `cyanea-core` -- error types

## Tests

77 tests across 6 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 44 | Module declarations, re-exports |
| `cluster.rs` | 812 | K-means, DBSCAN, hierarchical clustering |
| `distance.rs` | 309 | Distance metrics and pairwise matrices |
| `encoding.rs` | 141 | One-hot and label encoding |
| `evaluate.rs` | 186 | Silhouette score evaluation |
| `kmer.rs` | 228 | K-mer counting and frequency vectors |
| `normalize.rs` | 255 | Min-max, z-score, L2 normalization |
| `embedding.rs` | 8 | Stub |
| `inference.rs` | 8 | Stub |
| `reduction.rs` | 8 | Stub |
