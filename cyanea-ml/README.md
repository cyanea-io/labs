# cyanea-ml

> Machine learning primitives for bioinformatics: clustering, classification, dimensionality reduction, and sequence encoding.

## What's Inside

- **Clustering** -- K-means, DBSCAN, hierarchical (single/complete/average/Ward linkage)
- **Distance metrics** -- Euclidean, Manhattan, cosine, Hamming; pairwise distance matrices
- **KNN** -- k-nearest neighbors classification and regression
- **Linear regression** -- ordinary least squares via normal equation
- **Decision trees** -- CART classifier with Gini impurity
- **Random forests** -- bagged ensemble with bootstrap sampling and feature bagging
- **Gradient boosted trees** -- regression, binary/multiclass classification, early stopping, feature importance
- **Feature selection** -- variance threshold, mutual information, recursive feature elimination, Lasso L1
- **Hidden Markov models** -- forward, backward, Viterbi, Baum-Welch EM
- **PCA** -- principal component analysis (with optional BLAS acceleration)
- **t-SNE** -- t-distributed stochastic neighbor embedding
- **UMAP** -- uniform manifold approximation and projection
- **Sequence encoding** -- one-hot and label encoding for DNA/RNA/protein
- **K-mer counting** -- frequency vectors for ML features
- **Sequence embeddings** -- k-mer frequency and composition vector embeddings
- **Classification metrics** -- confusion matrix, precision/recall/F1, ROC/PR curves, AUC, MCC
- **Cross-validation** -- k-fold, stratified k-fold, leave-one-out
- **Normalization** -- min-max, z-score, L2 (sample-wise and column-wise)
- **Evaluation** -- silhouette score for cluster quality

## Quick Start

```toml
[dependencies]
cyanea-ml = "0.1"
```

```rust
use cyanea_ml::{
    cluster::{kmeans, KMeansConfig},
    distance::pairwise_distances,
    distance::DistanceMetric,
};

let data = vec![1.0, 2.0, 10.0, 11.0, 20.0, 21.0]; // 3 points, 2 features
let config = KMeansConfig { k: 2, max_iterations: 100, tolerance: 1e-4 };
let result = kmeans(&data, config).unwrap();
println!("Labels: {:?}", result.labels);
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `parallel` | No | Rayon parallelism |
| `blas` | No | BLAS-accelerated PCA via ndarray |

## Modules

| Module | Description |
|--------|-------------|
| `cluster` | K-means, DBSCAN, hierarchical clustering |
| `distance` | Distance metrics and pairwise distance matrices |
| `encoding` | One-hot and label encoding for sequences |
| `kmer` | K-mer counting and frequency vectors |
| `embedding` | K-mer frequency and composition vector embeddings |
| `normalize` | Min-max, z-score, L2 normalization |
| `inference` | KNN classification/regression, linear regression |
| `tree` | Decision tree classifier |
| `forest` | Random forest classifier |
| `gbdt` | Gradient boosted decision trees |
| `feature_selection` | Variance threshold, MI, RFE, Lasso |
| `hmm` | Hidden Markov Model (forward/backward/Viterbi/Baum-Welch) |
| `reduction` | PCA, t-SNE |
| `umap` | UMAP dimensionality reduction |
| `metrics` | Confusion matrix, ROC/PR curves, F1, MCC |
| `cross_validation` | K-fold, stratified k-fold, leave-one-out |
| `evaluate` | Silhouette score |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
