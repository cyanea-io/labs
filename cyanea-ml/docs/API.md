# cyanea-ml

Machine learning primitives for bioinformatics: clustering, distance metrics, sequence encoding, normalization, evaluation, classification metrics, cross-validation, decision trees, random forests, gradient boosted decision trees, feature selection, hidden Markov models, dimensionality reduction (PCA, t-SNE, UMAP).

## Status: Complete

All ML primitives are implemented across 17 modules covering clustering, distance metrics, sequence encoding, evaluation, normalization, k-mer counting, sequence embeddings, KNN/linear regression, dimensionality reduction (PCA, t-SNE, UMAP), decision tree and random forest classifiers, gradient boosted decision trees (regression, binary/multiclass classification, early stopping, feature importance), feature selection (variance threshold, mutual information, recursive feature elimination, Lasso L1), hidden Markov models (forward, backward, Viterbi, Baum-Welch), classification metrics (confusion matrix, ROC/PR curves, F1, MCC), and cross-validation (k-fold, stratified k-fold, leave-one-out).

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

### Classification metrics (`metrics.rs`)

**Confusion Matrix:**

| Type/Function | Description |
|---------------|-------------|
| `ConfusionMatrix` | Row-major confusion matrix for multi-class classification |
| `ConfusionMatrix::from_labels(actual, predicted, n_classes) -> Result<Self>` | Build from actual/predicted label vectors |
| `ConfusionMatrix::get(actual, predicted) -> usize` | Get count for (actual, predicted) pair |
| `ConfusionMatrix::total() -> usize` | Total number of samples |
| `ConfusionMatrix::true_positives(class) -> usize` | TP for a class |
| `ConfusionMatrix::false_positives(class) -> usize` | FP for a class |
| `ConfusionMatrix::true_negatives(class) -> usize` | TN for a class |
| `ConfusionMatrix::false_negatives(class) -> usize` | FN for a class |
| `ConfusionMatrix::accuracy() -> f64` | Overall accuracy |
| `ConfusionMatrix::precision(class) -> f64` | TP / (TP + FP) |
| `ConfusionMatrix::recall(class) -> f64` | TP / (TP + FN) |
| `ConfusionMatrix::f1(class) -> f64` | Harmonic mean of precision and recall |
| `ConfusionMatrix::specificity(class) -> f64` | TN / (TN + FP) |

**Standalone Scalar Metrics:**

| Function | Description |
|----------|-------------|
| `accuracy(actual, predicted) -> Result<f64>` | Overall accuracy |
| `f1_score(actual, predicted, class) -> Result<f64>` | F1 for a specific class |
| `f1_macro(actual, predicted) -> Result<f64>` | Macro-averaged F1 |
| `f1_weighted(actual, predicted) -> Result<f64>` | Weighted-averaged F1 |
| `matthews_corrcoef(actual, predicted) -> Result<f64>` | Matthews correlation coefficient |

**ROC Curve:**

| Type/Function | Description |
|---------------|-------------|
| `RocPoint` | `threshold`, `fpr`, `tpr` |
| `RocCurve` | `points: Vec<RocPoint>`, `auc: f64` |
| `roc_curve(scores, labels) -> Result<RocCurve>` | Compute ROC curve from scores and binary labels |
| `roc_auc(scores, labels) -> Result<f64>` | AUC shorthand |

**Precision-Recall Curve:**

| Type/Function | Description |
|---------------|-------------|
| `PrPoint` | `threshold`, `precision`, `recall` |
| `PrCurve` | `points: Vec<PrPoint>`, `auc: f64` |
| `pr_curve(scores, labels) -> Result<PrCurve>` | Compute PR curve from scores and binary labels |
| `pr_auc(scores, labels) -> Result<f64>` | AUC shorthand |

### Cross-validation (`cross_validation.rs`)

| Type/Function | Description |
|---------------|-------------|
| `FoldResult` | `fold`, `n_train`, `n_test`, `score` |
| `CvResult` | `folds: Vec<FoldResult>`, `mean_score`, `std_score` |
| `cross_validate_kfold(n_samples, k, seed, eval_fn) -> Result<CvResult>` | K-fold CV with shuffled splits |
| `cross_validate_stratified(labels, k, seed, eval_fn) -> Result<CvResult>` | Stratified k-fold preserving class proportions |
| `cross_validate_loo(n_samples, eval_fn) -> Result<CvResult>` | Leave-one-out CV |

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

### Sequence embeddings (`embedding.rs`)

| Type/Function | Description |
|---------------|-------------|
| `EmbeddingConfig` | `k`, `alphabet`, `normalize` |
| `SequenceEmbedding` | `vector`, `dim` |
| `kmer_embedding(seq, config) -> Result<SequenceEmbedding>` | K-mer frequency embedding |
| `composition_vector(seq, alphabet) -> Result<SequenceEmbedding>` | Nucleotide/amino acid composition |
| `batch_embed(sequences, config) -> Result<Vec<SequenceEmbedding>>` | Batch embedding |
| `pairwise_cosine_distances(embeddings) -> Result<Vec<f64>>` | Pairwise cosine distance matrix |

### Inference models (`inference.rs`)

**K-Nearest Neighbors:**

| Type/Function | Description |
|---------------|-------------|
| `KnnConfig` | `k`, `metric` |
| `KnnModel::fit(data, n_features, config) -> Result<Self>` | Fit KNN model |
| `KnnModel::neighbors(query) -> Result<Vec<(usize, f64)>>` | Find k nearest neighbors |
| `KnnModel::classify(query, labels) -> Result<i32>` | KNN classification (majority vote) |
| `KnnModel::regress(query, targets) -> Result<f64>` | KNN regression (weighted average) |

**Linear Regression:**

| Type/Function | Description |
|---------------|-------------|
| `LinearRegression::fit(data, n_features, targets) -> Result<Self>` | Fit via normal equation |
| `LinearRegression::predict(query) -> Result<f64>` | Predict single sample |
| `LinearRegression::predict_batch(queries) -> Result<Vec<f64>>` | Predict batch |
| Fields: `weights`, `bias`, `n_features`, `r_squared` | Model parameters |

### Decision tree (`tree.rs`)

| Type/Function | Description |
|---------------|-------------|
| `DecisionTree` | CART-style decision tree classifier using Gini impurity |
| `DecisionTree::fit(data, n_features, labels, max_depth) -> Result<Self>` | Fit a decision tree on flat row-major data |
| `DecisionTree::predict(sample) -> usize` | Predict class label for a single sample |
| `DecisionTree::predict_batch(data, n_features) -> Vec<usize>` | Predict class labels for multiple samples |

### Random forest (`forest.rs`)

| Type/Function | Description |
|---------------|-------------|
| `RandomForestConfig` | `n_trees`, `max_depth`, `max_features`, `seed` |
| `RandomForest` | Bagged ensemble of `DecisionTree` classifiers with bootstrap sampling and feature bagging |
| `RandomForest::fit(data, n_features, labels, config) -> Result<Self>` | Fit a random forest on flat row-major data |
| `RandomForest::predict(sample) -> usize` | Predict class label via majority vote |
| `RandomForest::predict_batch(data, n_features) -> Vec<usize>` | Predict class labels for multiple samples |
| `RandomForest::feature_importance(n_features) -> Vec<f64>` | Normalized split-frequency feature importance |
| `RandomForest::n_trees() -> usize` | Number of trees in the forest |
| `RandomForest::n_classes() -> usize` | Number of classes discovered during fitting |

### Gradient Boosted Decision Trees (`gbdt.rs`)

| Type/Function | Description |
|---------------|-------------|
| `GbdtConfig` | `n_estimators`, `learning_rate`, `max_depth`, `min_samples_leaf`, `subsample`, `max_features`, `seed`, `early_stopping_rounds`, `validation_fraction` |
| `GradientBoostedTrees` | Gradient boosted ensemble of internal regression trees |
| `GradientBoostedTrees::fit_regression(data, n_features, targets, config) -> Result<Self>` | Fit regression model (MSE loss) |
| `GradientBoostedTrees::fit_classification(data, n_features, labels, config) -> Result<Self>` | Fit classification model (binary log-loss or multiclass cross-entropy) |
| `GradientBoostedTrees::predict(sample) -> f64` | Raw prediction for a single sample |
| `GradientBoostedTrees::predict_batch(data, n_features) -> Vec<f64>` | Raw predictions for multiple samples |
| `GradientBoostedTrees::predict_class(sample) -> usize` | Predicted class label (classification only) |
| `GradientBoostedTrees::predict_class_batch(data, n_features) -> Vec<usize>` | Predicted class labels for multiple samples |
| `GradientBoostedTrees::predict_proba(sample) -> Vec<f64>` | Class probabilities for a single sample |
| `GradientBoostedTrees::predict_proba_batch(data, n_features) -> Vec<f64>` | Class probabilities for multiple samples (flat) |
| `GradientBoostedTrees::feature_importance() -> Vec<f64>` | Impurity-based feature importance (normalized, sums to 1.0) |
| `GradientBoostedTrees::permutation_importance_regression(data, n_features, targets, n_repeats, seed) -> Vec<f64>` | Permutation importance (MSE increase) |
| `GradientBoostedTrees::permutation_importance_classification(data, n_features, labels, n_repeats, seed) -> Vec<f64>` | Permutation importance (accuracy decrease) |
| `GradientBoostedTrees::n_estimators() -> usize` | Actual boosting rounds used |
| `GradientBoostedTrees::n_features() -> usize` | Number of features |
| `GradientBoostedTrees::is_regression() -> bool` | Whether model is regression |
| `GradientBoostedTrees::n_classes() -> usize` | Number of classes (0 for regression) |

### Feature Selection (`feature_selection.rs`)

**Variance Threshold:**

| Type/Function | Description |
|---------------|-------------|
| `FeatureSelection` | Result struct: `selected` (indices), `scores` (per-feature), `n_features` |
| `FeatureSelection::transform(data) -> Vec<f64>` | Filter data to keep only selected features |
| `FeatureSelection::n_selected() -> usize` | Number of selected features |
| `variance_threshold(data, n_features, threshold) -> Result<FeatureSelection>` | Select features with variance > threshold |

**Mutual Information:**

| Type/Function | Description |
|---------------|-------------|
| `mutual_information(data, n_features, labels, n_bins) -> Result<FeatureSelection>` | Select features with MI > 0 (discretized into equal-width bins) |
| `mutual_information_top_k(data, n_features, labels, n_bins, k) -> Result<FeatureSelection>` | Select top-k features by MI |

**Recursive Feature Elimination:**

| Type/Function | Description |
|---------------|-------------|
| `recursive_feature_elimination(data, n_features, n_select, importance_fn) -> Result<FeatureSelection>` | Iteratively remove least important features via user-supplied importance closure |

**Lasso (L1-Regularized) Selection:**

| Type/Function | Description |
|---------------|-------------|
| `LassoResult` | Result struct: `weights`, `bias`, `selected`, `n_iterations`, `n_features` |
| `LassoResult::transform(data) -> Vec<f64>` | Filter data to keep only non-zero-weight features |
| `LassoResult::predict(data) -> Vec<f64>` | Predict target values using fitted Lasso model |
| `LassoResult::n_selected() -> usize` | Number of non-zero-weight features |
| `lasso_selection(data, n_features, targets, alpha, max_iter, tol) -> Result<LassoResult>` | Coordinate descent Lasso for sparse feature selection |

### Hidden Markov Model (`hmm.rs`)

| Type/Function | Description |
|---------------|-------------|
| `HmmModel` | Discrete Hidden Markov Model with log-space arithmetic |
| `HmmModel::new(n_states, n_symbols, initial, transition, emission) -> Result<Self>` | Create and validate a new HMM |
| `HmmModel::forward(observations) -> Result<(Vec<Vec<f64>>, f64)>` | Forward algorithm (log-space alpha matrix + log-likelihood) |
| `HmmModel::backward(observations) -> Result<Vec<Vec<f64>>>` | Backward algorithm (log-space beta matrix) |
| `HmmModel::viterbi(observations) -> Result<(Vec<usize>, f64)>` | Viterbi decoding (most likely state path + log-probability) |
| `HmmModel::log_likelihood(observations) -> Result<f64>` | Log-likelihood of observation sequence |
| `HmmModel::baum_welch(observations, max_iter, tolerance) -> Result<f64>` | Baum-Welch (EM) parameter re-estimation |
| `HmmModel::n_states() -> usize` | Number of hidden states |
| `HmmModel::n_symbols() -> usize` | Number of observable symbols |

### Dimensionality reduction (`reduction.rs`)

**PCA:**

| Type/Function | Description |
|---------------|-------------|
| `PcaConfig` | `n_components`, `max_iter`, `tolerance` |
| `PcaResult` | `components`, `explained_variance`, `explained_variance_ratio`, `transformed`, `mean`, `n_features`, `n_components` |
| `pca(data, n_features, config) -> Result<PcaResult>` | Principal Component Analysis |

When the `blas` feature is enabled, PCA automatically dispatches to an ndarray-backed implementation (`blas_pca.rs`) that uses BLAS for covariance matrix computation and power iteration. The public API (`pca()`) is unchanged; the acceleration is transparent.

**t-SNE:**

| Type/Function | Description |
|---------------|-------------|
| `TsneConfig` | `n_components`, `perplexity`, `learning_rate`, `n_iter`, `seed` |
| `TsneResult` | `embedding`, `n_samples`, `n_components`, `kl_divergence` |
| `tsne(data, n_features, config) -> Result<TsneResult>` | t-distributed Stochastic Neighbor Embedding |

### UMAP (`umap.rs`)

| Type/Function | Description |
|---------------|-------------|
| `UmapInit` | Enum: `Random`, `Pca` |
| `UmapConfig` | `n_components`, `n_neighbors`, `min_dist`, `spread`, `learning_rate`, `n_epochs`, `negative_sample_rate`, `metric`, `init`, `seed` |
| `UmapResult` | `embedding`, `n_samples`, `n_components`, `n_epochs` |
| `umap(data, n_features, config) -> Result<UmapResult>` | Uniform Manifold Approximation and Projection |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `parallel` | No | Rayon parallelism |
| `blas` | No | BLAS-accelerated PCA via ndarray |

## Dependencies

- `cyanea-core` -- error types
- `ndarray` -- (optional, `blas` feature) matrix operations for BLAS-accelerated PCA
- `rayon` -- (optional, `parallel` feature) parallel iterators
- `serde` -- (optional, `serde` feature) serialization

## Tests

269 unit tests + 4 doc tests across 19 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 70 | Module declarations, re-exports |
| `cluster.rs` | 812 | K-means, DBSCAN, hierarchical clustering |
| `metrics.rs` | 798 | Confusion matrix, classification metrics, ROC/PR curves |
| `cross_validation.rs` | 472 | K-fold, stratified k-fold, leave-one-out CV |
| `distance.rs` | 309 | Distance metrics and pairwise matrices |
| `encoding.rs` | 141 | One-hot and label encoding |
| `evaluate.rs` | 186 | Silhouette score evaluation |
| `kmer.rs` | 228 | K-mer counting and frequency vectors |
| `normalize.rs` | 255 | Min-max, z-score, L2 normalization |
| `embedding.rs` | 262 | K-mer and composition vector embeddings |
| `inference.rs` | 546 | KNN and linear regression |
| `tree.rs` | 510 | Decision tree classifier (Gini impurity) |
| `forest.rs` | 447 | Random forest classifier (bagged ensemble) |
| `gbdt.rs` | 2165 | Gradient boosted decision trees (regression, classification, early stopping) |
| `feature_selection.rs` | 1226 | Feature selection (variance threshold, MI, RFE, Lasso) |
| `hmm.rs` | 770 | Hidden Markov Model (forward, backward, Viterbi, Baum-Welch) |
| `reduction.rs` | 773 | PCA and t-SNE dimensionality reduction |
| `blas_pca.rs` | 123 | BLAS-accelerated PCA via ndarray (feature `blas`) |
| `umap.rs` | 630 | UMAP dimensionality reduction |
