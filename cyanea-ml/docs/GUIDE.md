# cyanea-ml Usage Guide

Practical examples for machine learning in bioinformatics: clustering, dimensionality reduction, classification, HMMs, and model evaluation.

## K-Means Clustering

```rust
use cyanea_ml::cluster::{kmeans, KMeansConfig};

// 6 points in 2D (flat row-major)
let data = vec![
    1.0, 2.0,   // cluster A
    1.5, 1.8,
    1.2, 2.1,
    10.0, 10.0, // cluster B
    10.5, 9.8,
    9.8, 10.2,
];
let config = KMeansConfig { k: 2, max_iterations: 100, tolerance: 1e-4 };
let result = kmeans(&data, config).unwrap();

println!("Labels: {:?}", result.labels);
println!("Iterations: {}", result.iterations_run);
for (i, c) in result.centroids.chunks(2).enumerate() {
    println!("Centroid {}: ({:.1}, {:.1})", i, c[0], c[1]);
}
```

## DBSCAN

```rust
use cyanea_ml::cluster::{dbscan, DbscanConfig};

let data = vec![
    1.0, 1.0,
    1.1, 1.2,
    0.9, 0.8,
    10.0, 10.0,
    10.1, 10.2,
    50.0, 50.0, // noise point
];
let config = DbscanConfig { eps: 1.5, min_samples: 2 };
let result = dbscan(&data, config).unwrap();

println!("Clusters: {}, Noise: {}", result.n_clusters, result.n_noise);
println!("Labels: {:?}", result.labels); // -1 = noise
```

## Pairwise Distance Matrices

```rust
use cyanea_ml::distance::{pairwise_distances, DistanceMetric, DistanceMatrix};

let data = vec![
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
];
let dm = pairwise_distances(&data, DistanceMetric::Euclidean).unwrap();
println!("Distance 0-1: {:.3}", dm.get(0, 1));
println!("Distance 0-2: {:.3}", dm.get(0, 2));
```

## PCA Dimensionality Reduction

```rust
use cyanea_ml::reduction::{pca, PcaConfig};

// 10 samples, 5 features
let data: Vec<f64> = (0..50).map(|i| i as f64 * 0.1).collect();
let config = PcaConfig {
    n_components: 2,
    max_iter: 100,
    tolerance: 1e-6,
};
let result = pca(&data, 5, config).unwrap();

println!("Variance explained: {:?}", result.explained_variance_ratio);
// result.transformed contains the projected data (10 samples x 2 components)
```

## t-SNE Visualization

```rust
use cyanea_ml::reduction::{tsne, TsneConfig};

// High-dimensional data projected to 2D for visualization
let data: Vec<f64> = vec![/* 100 samples x 50 features */];
let config = TsneConfig {
    n_components: 2,
    perplexity: 30.0,
    learning_rate: 200.0,
    n_iter: 1000,
    seed: 42,
};
let result = tsne(&data, 50, config).unwrap();
println!("KL divergence: {:.4}", result.kl_divergence);
// result.embedding has shape 100 x 2
```

## UMAP Embedding

```rust
use cyanea_ml::umap::{umap, UmapConfig, UmapInit};
use cyanea_ml::distance::DistanceMetric;

let data: Vec<f64> = vec![/* samples x features */];
let config = UmapConfig {
    n_components: 2,
    n_neighbors: 15,
    min_dist: 0.1,
    spread: 1.0,
    learning_rate: 1.0,
    n_epochs: 200,
    negative_sample_rate: 5,
    metric: DistanceMetric::Euclidean,
    init: UmapInit::Pca,
    seed: 42,
};
let result = umap(&data, 50, config).unwrap();
```

## Random Forest Classification

```rust
use cyanea_ml::forest::{RandomForest, RandomForestConfig};

// Training data: 6 samples, 3 features
let train_data = vec![
    1.0, 0.0, 0.0,
    1.1, 0.1, 0.0,
    0.0, 1.0, 0.0,
    0.1, 1.1, 0.1,
    0.0, 0.0, 1.0,
    0.0, 0.1, 0.9,
];
let labels = vec![0, 0, 1, 1, 2, 2];

let config = RandomForestConfig {
    n_trees: 100,
    max_depth: Some(5),
    max_features: None, // sqrt(n_features) by default
    seed: 42,
};
let forest = RandomForest::fit(&train_data, 3, &labels, config).unwrap();

// Predict
let test_sample = vec![0.9, 0.2, 0.0];
let predicted = forest.predict(&test_sample);
println!("Predicted class: {}", predicted);

// Feature importance
let importance = forest.feature_importance(3);
println!("Feature importance: {:?}", importance);
```

## Gradient Boosted Trees

```rust
use cyanea_ml::gbdt::{GradientBoostedTrees, GbdtConfig};

// Regression
let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]; // 3 samples x 2 features
let targets = vec![1.5, 3.5, 5.5];

let config = GbdtConfig {
    n_estimators: 100,
    learning_rate: 0.1,
    max_depth: 3,
    early_stopping_rounds: Some(10),
    validation_fraction: Some(0.2),
    ..Default::default()
};
let model = GradientBoostedTrees::fit_regression(&data, 2, &targets, config).unwrap();
let prediction = model.predict(&[2.5, 3.5]);

// Classification
let labels = vec![0, 1, 0];
let clf = GradientBoostedTrees::fit_classification(&data, 2, &labels, Default::default()).unwrap();
let class = clf.predict_class(&[2.5, 3.5]);
let probas = clf.predict_proba(&[2.5, 3.5]);
println!("Class: {}, P(class=1): {:.3}", class, probas[1]);
```

## Hidden Markov Model

```rust
use cyanea_ml::hmm::HmmModel;

// 2 states (CpG island / non-CpG), 4 symbols (A, C, G, T)
let initial = vec![0.5, 0.5];
let transition = vec![
    0.7, 0.3,  // state 0 -> state 0, state 0 -> state 1
    0.4, 0.6,  // state 1 -> state 0, state 1 -> state 1
];
let emission = vec![
    0.15, 0.35, 0.35, 0.15, // state 0 emits: A, C, G, T (CpG-rich)
    0.30, 0.20, 0.20, 0.30, // state 1 emits: A, C, G, T (AT-rich)
];

let mut hmm = HmmModel::new(2, 4, initial, transition, emission).unwrap();

// Viterbi decoding: most likely state sequence
let observations = vec![1, 2, 1, 2, 0, 3, 0, 3]; // C, G, C, G, A, T, A, T
let (path, log_prob) = hmm.viterbi(&observations).unwrap();
println!("Path: {:?} (log P = {:.2})", path, log_prob);

// Forward algorithm: observation likelihood
let ll = hmm.log_likelihood(&observations).unwrap();
println!("Log-likelihood: {:.3}", ll);

// Baum-Welch: re-estimate parameters from data
let final_ll = hmm.baum_welch(&observations, 50, 1e-6).unwrap();
```

## Cross-Validation

```rust
use cyanea_ml::cross_validation::*;

let n_samples = 100;

// K-fold cross-validation
let result = cross_validate_kfold(n_samples, 5, 42, |train_idx, test_idx| {
    // Train model on train_idx, evaluate on test_idx
    // Return a score (e.g., accuracy)
    Ok(0.85) // placeholder
}).unwrap();

println!("Mean score: {:.3} +/- {:.3}", result.mean_score, result.std_score);

// Stratified k-fold (preserves class proportions)
let labels = vec![0; 50].into_iter().chain(vec![1; 50]).collect::<Vec<_>>();
let stratified = cross_validate_stratified(&labels, 5, 42, |train_idx, test_idx| {
    Ok(0.90)
}).unwrap();
```

## Feature Selection

```rust
use cyanea_ml::feature_selection::*;

// Data: 100 samples x 10 features
let data: Vec<f64> = vec![/* ... */];

// Variance threshold: remove low-variance features
let selected = variance_threshold(&data, 10, 0.01).unwrap();
println!("Selected {} of 10 features", selected.n_selected());

// Mutual information: top-k features by MI with labels
let labels = vec![0usize; 50].into_iter().chain(vec![1; 50]).collect::<Vec<_>>();
let mi = mutual_information_top_k(&data, 10, &labels, 10, 5).unwrap();
println!("Top 5 features: {:?}", mi.selected);

// Lasso L1 selection
let targets: Vec<f64> = (0..100).map(|i| i as f64).collect();
let lasso = lasso_selection(&data, 10, &targets, 0.1, 1000, 1e-4).unwrap();
println!("Non-zero features: {:?}", lasso.selected);
```

## Classification Metrics

```rust
use cyanea_ml::metrics::*;

let actual    = vec![0, 0, 1, 1, 2, 2, 0, 1];
let predicted = vec![0, 1, 1, 1, 2, 0, 0, 1];

// Confusion matrix
let cm = ConfusionMatrix::from_labels(&actual, &predicted, 3).unwrap();
println!("Accuracy: {:.3}", cm.accuracy());
println!("Class 1 precision: {:.3}", cm.precision(1));
println!("Class 1 recall: {:.3}", cm.recall(1));
println!("Class 1 F1: {:.3}", cm.f1(1));

// ROC curve (binary)
let scores = vec![0.9, 0.4, 0.8, 0.3, 0.7, 0.2];
let binary_labels = vec![1, 0, 1, 0, 1, 0];
let roc = roc_curve(&scores, &binary_labels).unwrap();
println!("AUC: {:.3}", roc.auc);

// Matthews Correlation Coefficient
let mcc = matthews_corrcoef(&actual, &predicted).unwrap();
println!("MCC: {:.3}", mcc);
```
