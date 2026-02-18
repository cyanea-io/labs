//! Gradient boosted decision trees (GBDT).
//!
//! Implements gradient boosting with internal regression trees for both
//! regression and classification (binary + multiclass) tasks. Features include
//! row/column subsampling, early stopping on a validation set, impurity-based
//! feature importance, and permutation importance.
//!
//! Data is flat row-major `&[f64]` with an `n_features` parameter, consistent
//! with the rest of the cyanea-ml crate.
//!
//! # Example
//!
//! ```
//! use cyanea_ml::gbdt::{GbdtConfig, GradientBoostedTrees};
//!
//! // Regression: fit y = 2*x
//! let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
//! let targets = vec![2.0, 4.0, 6.0, 8.0, 10.0];
//! let config = GbdtConfig { n_estimators: 50, ..Default::default() };
//! let model = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
//! let pred = model.predict(&[3.0]);
//! assert!((pred - 6.0).abs() < 1.0);
//! ```

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// PRNG (per-module pattern, same as forest.rs)
// ---------------------------------------------------------------------------

struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn next_bounded(&mut self, bound: u64) -> u64 {
        self.next_u64() % bound
    }

}

/// Shuffle a slice in-place using Fisher-Yates.
fn shuffle_indices(rng: &mut LcgRng, indices: &mut [usize]) {
    let n = indices.len();
    for i in (1..n).rev() {
        let j = rng.next_bounded((i + 1) as u64) as usize;
        indices.swap(i, j);
    }
}

/// Sample `count` indices from `0..n` without replacement.
fn subsample_indices(rng: &mut LcgRng, n: usize, count: usize) -> Vec<usize> {
    let count = count.min(n);
    let mut pool: Vec<usize> = (0..n).collect();
    for i in 0..count {
        let j = i + rng.next_bounded((n - i) as u64) as usize;
        pool.swap(i, j);
    }
    pool.truncate(count);
    pool
}

/// Select `count` distinct feature indices from `0..n_features`.
fn random_feature_subset(rng: &mut LcgRng, n_features: usize, count: usize) -> Vec<usize> {
    let count = count.min(n_features);
    if count == n_features {
        return (0..n_features).collect();
    }
    subsample_indices(rng, n_features, count)
}

// ---------------------------------------------------------------------------
// Internal regression tree
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
#[allow(dead_code)]
enum RegTreeNode {
    Split {
        feature_idx: usize,
        threshold: f64,
        left: usize,
        right: usize,
        impurity_decrease: f64,
        n_samples: usize,
    },
    Leaf {
        value: f64,
        n_samples: usize,
    },
}

#[derive(Debug, Clone)]
struct RegressionTree {
    nodes: Vec<RegTreeNode>,
}

impl RegressionTree {
    fn fit(
        data: &[f64],
        n_features: usize,
        targets: &[f64],
        indices: &[usize],
        candidate_features: &[usize],
        max_depth: usize,
        min_samples_leaf: usize,
    ) -> Self {
        let mut nodes = Vec::new();
        build_reg_tree(
            data,
            n_features,
            targets,
            indices,
            candidate_features,
            max_depth,
            min_samples_leaf,
            0,
            &mut nodes,
        );
        Self { nodes }
    }

    fn predict(&self, sample: &[f64]) -> f64 {
        let mut idx = 0;
        loop {
            match &self.nodes[idx] {
                RegTreeNode::Leaf { value, .. } => return *value,
                RegTreeNode::Split {
                    feature_idx,
                    threshold,
                    left,
                    right,
                    ..
                } => {
                    if sample[*feature_idx] <= *threshold {
                        idx = *left;
                    } else {
                        idx = *right;
                    }
                }
            }
        }
    }

    fn nodes(&self) -> &[RegTreeNode] {
        &self.nodes
    }
}

/// Compute MSE for a set of target values given by indices.
fn mse(targets: &[f64], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }
    let n = indices.len() as f64;
    let mean = indices.iter().map(|&i| targets[i]).sum::<f64>() / n;
    indices.iter().map(|&i| (targets[i] - mean).powi(2)).sum::<f64>() / n
}

/// Compute the mean of target values at the given indices.
fn mean_target(targets: &[f64], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }
    indices.iter().map(|&i| targets[i]).sum::<f64>() / indices.len() as f64
}

fn build_reg_tree(
    data: &[f64],
    n_features: usize,
    targets: &[f64],
    indices: &[usize],
    candidate_features: &[usize],
    max_depth: usize,
    min_samples_leaf: usize,
    depth: usize,
    nodes: &mut Vec<RegTreeNode>,
) -> usize {
    let leaf_value = mean_target(targets, indices);

    // Stop conditions
    if depth >= max_depth || indices.len() < 2 || indices.len() <= min_samples_leaf {
        let idx = nodes.len();
        nodes.push(RegTreeNode::Leaf {
            value: leaf_value,
            n_samples: indices.len(),
        });
        return idx;
    }

    // Check if all targets are the same
    let first_target = targets[indices[0]];
    if indices.iter().all(|&i| (targets[i] - first_target).abs() < 1e-15) {
        let idx = nodes.len();
        nodes.push(RegTreeNode::Leaf {
            value: leaf_value,
            n_samples: indices.len(),
        });
        return idx;
    }

    // Find best split
    if let Some((best_feature, best_threshold, best_decrease)) =
        find_best_mse_split(data, n_features, targets, indices, candidate_features, min_samples_leaf)
    {
        let (left_indices, right_indices) =
            partition(data, n_features, indices, best_feature, best_threshold);

        if left_indices.is_empty() || right_indices.is_empty() {
            let idx = nodes.len();
            nodes.push(RegTreeNode::Leaf {
                value: leaf_value,
                n_samples: indices.len(),
            });
            return idx;
        }

        let node_idx = nodes.len();
        nodes.push(RegTreeNode::Leaf {
            value: 0.0,
            n_samples: 0,
        }); // placeholder

        let left_child = build_reg_tree(
            data,
            n_features,
            targets,
            &left_indices,
            candidate_features,
            max_depth,
            min_samples_leaf,
            depth + 1,
            nodes,
        );
        let right_child = build_reg_tree(
            data,
            n_features,
            targets,
            &right_indices,
            candidate_features,
            max_depth,
            min_samples_leaf,
            depth + 1,
            nodes,
        );

        nodes[node_idx] = RegTreeNode::Split {
            feature_idx: best_feature,
            threshold: best_threshold,
            left: left_child,
            right: right_child,
            impurity_decrease: best_decrease,
            n_samples: indices.len(),
        };

        node_idx
    } else {
        let idx = nodes.len();
        nodes.push(RegTreeNode::Leaf {
            value: leaf_value,
            n_samples: indices.len(),
        });
        idx
    }
}

/// Find the best (feature, threshold) split by minimizing weighted MSE.
fn find_best_mse_split(
    data: &[f64],
    n_features: usize,
    targets: &[f64],
    indices: &[usize],
    candidate_features: &[usize],
    min_samples_leaf: usize,
) -> Option<(usize, f64, f64)> {
    let n = indices.len();
    let parent_mse = mse(targets, indices);
    let parent_total = parent_mse * n as f64;

    let mut best_decrease = 0.0;
    let mut best_feature = 0;
    let mut best_threshold = 0.0;
    let mut found = false;

    for &feat in candidate_features {
        // Collect (value, target) pairs and sort by value
        let mut pairs: Vec<(f64, f64)> = indices
            .iter()
            .map(|&i| (data[i * n_features + feat], targets[i]))
            .collect();
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Running sums for incremental MSE computation
        let total_sum: f64 = pairs.iter().map(|&(_, t)| t).sum();
        let total_sq_sum: f64 = pairs.iter().map(|&(_, t)| t * t).sum();

        let mut left_sum = 0.0;
        let mut left_sq_sum = 0.0;
        let mut left_count = 0usize;

        for i in 0..pairs.len() - 1 {
            left_sum += pairs[i].1;
            left_sq_sum += pairs[i].1 * pairs[i].1;
            left_count += 1;

            // Skip if same feature value as next (no valid threshold between them)
            if (pairs[i].0 - pairs[i + 1].0).abs() < 1e-15 {
                continue;
            }

            let right_count = n - left_count;

            // Check min_samples_leaf
            if left_count < min_samples_leaf || right_count < min_samples_leaf {
                continue;
            }

            let left_mean = left_sum / left_count as f64;
            let left_mse_val = left_sq_sum / left_count as f64 - left_mean * left_mean;

            let right_sum = total_sum - left_sum;
            let right_sq_sum = total_sq_sum - left_sq_sum;
            let right_mean = right_sum / right_count as f64;
            let right_mse_val = right_sq_sum / right_count as f64 - right_mean * right_mean;

            let weighted_mse =
                left_mse_val * left_count as f64 + right_mse_val * right_count as f64;
            let decrease = parent_total - weighted_mse;

            if decrease > best_decrease {
                best_decrease = decrease;
                best_feature = feat;
                best_threshold = (pairs[i].0 + pairs[i + 1].0) / 2.0;
                found = true;
            }
        }
    }

    if found {
        Some((best_feature, best_threshold, best_decrease))
    } else {
        None
    }
}

/// Partition sample indices by threshold on a given feature.
fn partition(
    data: &[f64],
    n_features: usize,
    indices: &[usize],
    feature_idx: usize,
    threshold: f64,
) -> (Vec<usize>, Vec<usize>) {
    let mut left = Vec::new();
    let mut right = Vec::new();
    for &i in indices {
        if data[i * n_features + feature_idx] <= threshold {
            left.push(i);
        } else {
            right.push(i);
        }
    }
    (left, right)
}

// ---------------------------------------------------------------------------
// Loss / activation helpers
// ---------------------------------------------------------------------------

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

fn softmax(raw_scores: &[f64]) -> Vec<f64> {
    // Log-sum-exp trick for numerical stability
    let max_val = raw_scores
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = raw_scores.iter().map(|&s| (s - max_val).exp()).collect();
    let sum: f64 = exps.iter().sum();
    exps.iter().map(|&e| e / sum).collect()
}

fn log_loss(y_true: f64, p: f64) -> f64 {
    let p = p.clamp(1e-15, 1.0 - 1e-15);
    -(y_true * p.ln() + (1.0 - y_true) * (1.0 - p).ln())
}

fn cross_entropy_loss(y_true_class: usize, probs: &[f64]) -> f64 {
    let p = probs[y_true_class].clamp(1e-15, 1.0);
    -p.ln()
}

// ---------------------------------------------------------------------------
// GbdtConfig
// ---------------------------------------------------------------------------

/// Configuration for gradient boosted decision trees.
#[derive(Debug, Clone)]
pub struct GbdtConfig {
    /// Number of boosting rounds (default: 100).
    pub n_estimators: usize,
    /// Learning rate / shrinkage factor (default: 0.1).
    pub learning_rate: f64,
    /// Maximum depth per regression tree (default: 5).
    pub max_depth: usize,
    /// Minimum number of samples in a leaf node (default: 1).
    pub min_samples_leaf: usize,
    /// Row subsampling ratio per round (default: 1.0 = no subsampling).
    pub subsample: f64,
    /// Number of features to consider at each split. `None` = all features.
    pub max_features: Option<usize>,
    /// Random seed for reproducibility (default: 42).
    pub seed: u64,
    /// If set, stop training after this many rounds with no improvement on
    /// a held-out validation set.
    pub early_stopping_rounds: Option<usize>,
    /// Fraction of training data held out for early stopping (default: 0.1).
    pub validation_fraction: f64,
}

impl Default for GbdtConfig {
    fn default() -> Self {
        Self {
            n_estimators: 100,
            learning_rate: 0.1,
            max_depth: 5,
            min_samples_leaf: 1,
            subsample: 1.0,
            max_features: None,
            seed: 42,
            early_stopping_rounds: None,
            validation_fraction: 0.1,
        }
    }
}

fn validate_config(config: &GbdtConfig) -> Result<()> {
    if config.n_estimators == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_estimators must be > 0".into(),
        ));
    }
    if config.learning_rate <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "learning_rate must be > 0".into(),
        ));
    }
    if config.subsample <= 0.0 || config.subsample > 1.0 {
        return Err(CyaneaError::InvalidInput(
            "subsample must be in (0, 1]".into(),
        ));
    }
    if config.early_stopping_rounds.is_some()
        && (config.validation_fraction <= 0.0 || config.validation_fraction >= 1.0)
    {
        return Err(CyaneaError::InvalidInput(
            "validation_fraction must be in (0, 1) when early stopping is enabled".into(),
        ));
    }
    Ok(())
}

fn validate_data(data: &[f64], n_features: usize) -> Result<usize> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    if n_features == 0 {
        return Err(CyaneaError::InvalidInput("n_features must be > 0".into()));
    }
    if data.len() % n_features != 0 {
        return Err(CyaneaError::InvalidInput(format!(
            "data length {} not divisible by n_features {}",
            data.len(),
            n_features
        )));
    }
    Ok(data.len() / n_features)
}

// ---------------------------------------------------------------------------
// Mode
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Regression,
    BinaryClassification,
    MulticlassClassification,
}

// ---------------------------------------------------------------------------
// GradientBoostedTrees
// ---------------------------------------------------------------------------

/// A gradient boosted decision tree model.
///
/// Supports regression, binary classification, and multiclass classification.
/// Each boosting round fits one or more internal regression trees to
/// pseudo-residuals of the current model.
#[derive(Debug, Clone)]
pub struct GradientBoostedTrees {
    /// Per-round trees. For multiclass, each round has `n_classes` trees.
    trees: Vec<RegressionTree>,
    mode: Mode,
    n_features: usize,
    n_classes: usize,
    /// Initial prediction (F_0). For multiclass, one per class.
    initial_prediction: Vec<f64>,
    learning_rate: f64,
    /// Pre-computed impurity-based feature importance (normalized).
    feature_importance_cache: Vec<f64>,
}

impl GradientBoostedTrees {
    // ── Fitting ─────────────────────────────────────────────────

    /// Fit a gradient boosted regression model.
    ///
    /// * `data` — flat row-major `n_samples x n_features`
    /// * `n_features` — number of features per sample
    /// * `targets` — continuous target value for each sample
    /// * `config` — GBDT hyper-parameters
    pub fn fit_regression(
        data: &[f64],
        n_features: usize,
        targets: &[f64],
        config: &GbdtConfig,
    ) -> Result<Self> {
        validate_config(config)?;
        let n_samples = validate_data(data, n_features)?;
        if targets.len() != n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "targets length {} != n_samples {}",
                targets.len(),
                n_samples
            )));
        }

        let mut rng = LcgRng::new(config.seed);
        let max_features = config.max_features.unwrap_or(n_features);

        // Split train/validation if early stopping
        let (train_idx, val_idx) = if config.early_stopping_rounds.is_some() {
            train_val_split(&mut rng, n_samples, config.validation_fraction)
        } else {
            ((0..n_samples).collect(), Vec::new())
        };

        // F_0 = mean(y_train)
        let f0 = mean_target(targets, &train_idx);
        let initial_prediction = vec![f0];

        // Current predictions for training samples
        let mut f_train: Vec<f64> = train_idx.iter().map(|_| f0).collect();
        let mut f_val: Vec<f64> = val_idx.iter().map(|_| f0).collect();

        let mut trees = Vec::new();
        let mut best_val_loss = f64::INFINITY;
        let mut rounds_without_improvement = 0usize;
        let patience = config.early_stopping_rounds.unwrap_or(0);

        for _ in 0..config.n_estimators {
            // Compute residuals on training set
            let residuals: Vec<f64> = train_idx
                .iter()
                .enumerate()
                .map(|(j, &i)| targets[i] - f_train[j])
                .collect();

            // Subsample rows
            let sub_count =
                ((train_idx.len() as f64 * config.subsample).round() as usize).max(1);
            let sub_indices = if config.subsample < 1.0 {
                subsample_indices(&mut rng, train_idx.len(), sub_count)
            } else {
                (0..train_idx.len()).collect()
            };

            // Feature subset
            let features = random_feature_subset(&mut rng, n_features, max_features);

            // Build tree on (train data, residuals) using sub_indices into train_idx
            // We need to create a temporary target array aligned to original indices
            // for the regression tree. Instead, build a local data/targets using
            // the training subset.
            let tree = fit_tree_on_subset(
                data,
                n_features,
                &residuals,
                &train_idx,
                &sub_indices,
                &features,
                config.max_depth,
                config.min_samples_leaf,
            );

            // Update predictions
            for (j, &ti) in train_idx.iter().enumerate() {
                let row = &data[ti * n_features..(ti + 1) * n_features];
                f_train[j] += config.learning_rate * tree.predict(row);
            }

            trees.push(tree);

            // Early stopping check
            if config.early_stopping_rounds.is_some() {
                // Update val predictions
                for (j, &vi) in val_idx.iter().enumerate() {
                    let row = &data[vi * n_features..(vi + 1) * n_features];
                    f_val[j] +=
                        config.learning_rate * trees.last().unwrap().predict(row);
                }

                // MSE on validation
                let val_loss: f64 = val_idx
                    .iter()
                    .enumerate()
                    .map(|(j, &i)| (targets[i] - f_val[j]).powi(2))
                    .sum::<f64>()
                    / val_idx.len() as f64;

                if val_loss < best_val_loss - 1e-10 {
                    best_val_loss = val_loss;
                    rounds_without_improvement = 0;
                } else {
                    rounds_without_improvement += 1;
                    if rounds_without_improvement >= patience {
                        // Trim trees to best round
                        let best_round = trees.len() - rounds_without_improvement;
                        trees.truncate(best_round.max(1));
                        break;
                    }
                }
            }
        }

        let fi = compute_feature_importance(&trees, n_features);

        Ok(Self {
            trees,
            mode: Mode::Regression,
            n_features,
            n_classes: 0,
            initial_prediction,
            learning_rate: config.learning_rate,
            feature_importance_cache: fi,
        })
    }

    /// Fit a gradient boosted classification model.
    ///
    /// * `data` — flat row-major `n_samples x n_features`
    /// * `n_features` — number of features per sample
    /// * `labels` — class label for each sample (`0..n_classes`)
    /// * `config` — GBDT hyper-parameters
    ///
    /// Automatically selects binary (log-loss) or multiclass (softmax
    /// cross-entropy) mode based on the number of distinct classes.
    pub fn fit_classification(
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        config: &GbdtConfig,
    ) -> Result<Self> {
        validate_config(config)?;
        let n_samples = validate_data(data, n_features)?;
        if labels.len() != n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "labels length {} != n_samples {}",
                labels.len(),
                n_samples
            )));
        }

        let n_classes = labels.iter().copied().max().map_or(0, |m| m + 1);
        if n_classes < 2 {
            // Degenerate: only one class
            return Self::fit_single_class(n_features, labels[0], config);
        }

        if n_classes == 2 {
            Self::fit_binary(data, n_features, labels, config)
        } else {
            Self::fit_multiclass(data, n_features, labels, n_classes, config)
        }
    }

    fn fit_single_class(
        n_features: usize,
        class: usize,
        config: &GbdtConfig,
    ) -> Result<Self> {
        // All samples are the same class — no trees needed
        let n_classes = class + 1;
        let mut initial = vec![f64::NEG_INFINITY; n_classes];
        initial[class] = 0.0; // will produce probability 1.0 for this class
        let fi = vec![0.0; n_features];
        Ok(Self {
            trees: Vec::new(),
            mode: if n_classes == 2 {
                Mode::BinaryClassification
            } else {
                Mode::MulticlassClassification
            },
            n_features,
            n_classes,
            initial_prediction: initial,
            learning_rate: config.learning_rate,
            feature_importance_cache: fi,
        })
    }

    fn fit_binary(
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        config: &GbdtConfig,
    ) -> Result<Self> {
        let n_samples = data.len() / n_features;
        let mut rng = LcgRng::new(config.seed);
        let max_features = config.max_features.unwrap_or(n_features);

        let (train_idx, val_idx) = if config.early_stopping_rounds.is_some() {
            train_val_split(&mut rng, n_samples, config.validation_fraction)
        } else {
            ((0..n_samples).collect(), Vec::new())
        };

        // F_0 = log(p / (1-p)) where p = mean(y_train)
        let p0 = train_idx.iter().map(|&i| labels[i] as f64).sum::<f64>()
            / train_idx.len() as f64;
        let p0 = p0.clamp(1e-8, 1.0 - 1e-8);
        let f0 = (p0 / (1.0 - p0)).ln();
        let initial_prediction = vec![f0];

        let mut f_train: Vec<f64> = train_idx.iter().map(|_| f0).collect();
        let mut f_val: Vec<f64> = val_idx.iter().map(|_| f0).collect();

        let mut trees = Vec::new();
        let mut best_val_loss = f64::INFINITY;
        let mut rounds_without_improvement = 0usize;
        let patience = config.early_stopping_rounds.unwrap_or(0);

        for _ in 0..config.n_estimators {
            // Pseudo-residuals: r_i = y_i - sigmoid(F(x_i))
            let residuals: Vec<f64> = train_idx
                .iter()
                .enumerate()
                .map(|(j, &i)| labels[i] as f64 - sigmoid(f_train[j]))
                .collect();

            let sub_count =
                ((train_idx.len() as f64 * config.subsample).round() as usize).max(1);
            let sub_indices = if config.subsample < 1.0 {
                subsample_indices(&mut rng, train_idx.len(), sub_count)
            } else {
                (0..train_idx.len()).collect()
            };

            let features = random_feature_subset(&mut rng, n_features, max_features);

            // Fit tree to pseudo-residuals
            let mut tree = fit_tree_on_subset(
                data,
                n_features,
                &residuals,
                &train_idx,
                &sub_indices,
                &features,
                config.max_depth,
                config.min_samples_leaf,
            );

            // Newton-Raphson leaf correction
            apply_binary_leaf_correction(
                &mut tree,
                data,
                n_features,
                &f_train,
                labels,
                &train_idx,
            );

            // Update predictions
            for (j, &ti) in train_idx.iter().enumerate() {
                let row = &data[ti * n_features..(ti + 1) * n_features];
                f_train[j] += config.learning_rate * tree.predict(row);
            }

            trees.push(tree);

            // Early stopping
            if config.early_stopping_rounds.is_some() {
                for (j, &vi) in val_idx.iter().enumerate() {
                    let row = &data[vi * n_features..(vi + 1) * n_features];
                    f_val[j] +=
                        config.learning_rate * trees.last().unwrap().predict(row);
                }

                let val_loss: f64 = val_idx
                    .iter()
                    .enumerate()
                    .map(|(j, &i)| log_loss(labels[i] as f64, sigmoid(f_val[j])))
                    .sum::<f64>()
                    / val_idx.len() as f64;

                if val_loss < best_val_loss - 1e-10 {
                    best_val_loss = val_loss;
                    rounds_without_improvement = 0;
                } else {
                    rounds_without_improvement += 1;
                    if rounds_without_improvement >= patience {
                        let best_round = trees.len() - rounds_without_improvement;
                        trees.truncate(best_round.max(1));
                        break;
                    }
                }
            }
        }

        let fi = compute_feature_importance(&trees, n_features);

        Ok(Self {
            trees,
            mode: Mode::BinaryClassification,
            n_features,
            n_classes: 2,
            initial_prediction,
            learning_rate: config.learning_rate,
            feature_importance_cache: fi,
        })
    }

    fn fit_multiclass(
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        n_classes: usize,
        config: &GbdtConfig,
    ) -> Result<Self> {
        let n_samples = data.len() / n_features;
        let mut rng = LcgRng::new(config.seed);
        let max_features = config.max_features.unwrap_or(n_features);

        let (train_idx, val_idx) = if config.early_stopping_rounds.is_some() {
            train_val_split(&mut rng, n_samples, config.validation_fraction)
        } else {
            ((0..n_samples).collect(), Vec::new())
        };

        // F_0^k = log(p_k) for each class k
        let mut class_counts = vec![0usize; n_classes];
        for &i in &train_idx {
            class_counts[labels[i]] += 1;
        }
        let n_train = train_idx.len() as f64;
        let initial_prediction: Vec<f64> = class_counts
            .iter()
            .map(|&c| {
                let p = (c as f64 / n_train).clamp(1e-8, 1.0 - 1e-8);
                p.ln()
            })
            .collect();

        // Current raw scores: f_train[j][k] for j-th training sample, class k
        let mut f_train: Vec<Vec<f64>> = train_idx
            .iter()
            .map(|_| initial_prediction.clone())
            .collect();
        let mut f_val: Vec<Vec<f64>> = val_idx
            .iter()
            .map(|_| initial_prediction.clone())
            .collect();

        let mut trees = Vec::new(); // K trees per round
        let mut best_val_loss = f64::INFINITY;
        let mut rounds_without_improvement = 0usize;
        let patience = config.early_stopping_rounds.unwrap_or(0);
        let mut actual_rounds = 0usize;

        for _ in 0..config.n_estimators {
            // For each class k, fit a tree
            for k in 0..n_classes {
                // Softmax probabilities for each training sample
                let probs: Vec<f64> = f_train.iter().map(|scores| softmax(scores)[k]).collect();

                // Pseudo-residuals: r_i^k = (y_i == k) - p_i^k
                let residuals: Vec<f64> = train_idx
                    .iter()
                    .enumerate()
                    .map(|(j, &i)| {
                        let indicator = if labels[i] == k { 1.0 } else { 0.0 };
                        indicator - probs[j]
                    })
                    .collect();

                let sub_count =
                    ((train_idx.len() as f64 * config.subsample).round() as usize).max(1);
                let sub_indices = if config.subsample < 1.0 {
                    subsample_indices(&mut rng, train_idx.len(), sub_count)
                } else {
                    (0..train_idx.len()).collect()
                };

                let features = random_feature_subset(&mut rng, n_features, max_features);

                let mut tree = fit_tree_on_subset(
                    data,
                    n_features,
                    &residuals,
                    &train_idx,
                    &sub_indices,
                    &features,
                    config.max_depth,
                    config.min_samples_leaf,
                );

                // Newton-Raphson multiclass leaf correction
                apply_multiclass_leaf_correction(
                    &mut tree,
                    data,
                    n_features,
                    &residuals,
                    &train_idx,
                    n_classes,
                );

                // Update training scores for class k
                for (j, &ti) in train_idx.iter().enumerate() {
                    let row = &data[ti * n_features..(ti + 1) * n_features];
                    f_train[j][k] += config.learning_rate * tree.predict(row);
                }

                trees.push(tree);
            }

            actual_rounds += 1;

            // Early stopping
            if config.early_stopping_rounds.is_some() {
                // Update val scores for all classes in this round
                let round_trees =
                    &trees[trees.len() - n_classes..];
                for (j, &vi) in val_idx.iter().enumerate() {
                    let row = &data[vi * n_features..(vi + 1) * n_features];
                    for (k, tree) in round_trees.iter().enumerate() {
                        f_val[j][k] += config.learning_rate * tree.predict(row);
                    }
                }

                let val_loss: f64 = val_idx
                    .iter()
                    .enumerate()
                    .map(|(j, &i)| cross_entropy_loss(labels[i], &softmax(&f_val[j])))
                    .sum::<f64>()
                    / val_idx.len() as f64;

                if val_loss < best_val_loss - 1e-10 {
                    best_val_loss = val_loss;
                    rounds_without_improvement = 0;
                } else {
                    rounds_without_improvement += 1;
                    if rounds_without_improvement >= patience {
                        let best_round = actual_rounds - rounds_without_improvement;
                        let keep = best_round.max(1) * n_classes;
                        trees.truncate(keep);
                        break;
                    }
                }
            }
        }

        let fi = compute_feature_importance(&trees, n_features);

        Ok(Self {
            trees,
            mode: Mode::MulticlassClassification,
            n_features,
            n_classes,
            initial_prediction,
            learning_rate: config.learning_rate,
            feature_importance_cache: fi,
        })
    }

    // ── Prediction ──────────────────────────────────────────────

    /// Predict the raw output for a single sample.
    ///
    /// For regression, this is the predicted value. For binary classification,
    /// this is the log-odds. For multiclass, this is the raw score for class 0
    /// (use [`predict_proba`](Self::predict_proba) for all class probabilities).
    pub fn predict(&self, sample: &[f64]) -> f64 {
        match self.mode {
            Mode::Regression => {
                let mut f = self.initial_prediction[0];
                for tree in &self.trees {
                    f += self.learning_rate * tree.predict(sample);
                }
                f
            }
            Mode::BinaryClassification => {
                let mut f = self.initial_prediction[0];
                for tree in &self.trees {
                    f += self.learning_rate * tree.predict(sample);
                }
                f
            }
            Mode::MulticlassClassification => {
                // Return raw score for class 0
                let scores = self.raw_scores(sample);
                scores[0]
            }
        }
    }

    /// Predict raw outputs for multiple samples.
    ///
    /// `data` is flat row-major with `n_features` columns.
    pub fn predict_batch(&self, data: &[f64], n_features: usize) -> Vec<f64> {
        let n_samples = data.len() / n_features;
        (0..n_samples)
            .map(|i| {
                let row = &data[i * n_features..(i + 1) * n_features];
                self.predict(row)
            })
            .collect()
    }

    /// Predict the class label for a single sample (classification only).
    ///
    /// For regression models, returns 0.
    pub fn predict_class(&self, sample: &[f64]) -> usize {
        match self.mode {
            Mode::Regression => 0,
            Mode::BinaryClassification => {
                let f = self.predict(sample);
                if sigmoid(f) >= 0.5 { 1 } else { 0 }
            }
            Mode::MulticlassClassification => {
                let scores = self.raw_scores(sample);
                scores
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|(i, _)| i)
                    .unwrap_or(0)
            }
        }
    }

    /// Predict class labels for multiple samples.
    pub fn predict_class_batch(&self, data: &[f64], n_features: usize) -> Vec<usize> {
        let n_samples = data.len() / n_features;
        (0..n_samples)
            .map(|i| {
                let row = &data[i * n_features..(i + 1) * n_features];
                self.predict_class(row)
            })
            .collect()
    }

    /// Predict class probabilities for a single sample.
    ///
    /// For binary classification, returns `[P(class=0), P(class=1)]`.
    /// For multiclass, returns `[P(class=0), ..., P(class=K-1)]`.
    /// For regression, returns an empty vector.
    pub fn predict_proba(&self, sample: &[f64]) -> Vec<f64> {
        match self.mode {
            Mode::Regression => Vec::new(),
            Mode::BinaryClassification => {
                let f = self.predict(sample);
                let p1 = sigmoid(f);
                vec![1.0 - p1, p1]
            }
            Mode::MulticlassClassification => {
                let scores = self.raw_scores(sample);
                softmax(&scores)
            }
        }
    }

    /// Predict class probabilities for multiple samples.
    ///
    /// Returns a flat vector of length `n_samples * n_classes`.
    pub fn predict_proba_batch(&self, data: &[f64], n_features: usize) -> Vec<f64> {
        let n_samples = data.len() / n_features;
        let mut result = Vec::with_capacity(n_samples * self.n_classes);
        for i in 0..n_samples {
            let row = &data[i * n_features..(i + 1) * n_features];
            result.extend(self.predict_proba(row));
        }
        result
    }

    /// Compute raw scores for all classes (multiclass).
    fn raw_scores(&self, sample: &[f64]) -> Vec<f64> {
        let mut scores = self.initial_prediction.clone();
        let k = self.n_classes;
        // Trees are stored K per round: [round0_class0, round0_class1, ..., round1_class0, ...]
        for tree_group in self.trees.chunks(k) {
            for (class_idx, tree) in tree_group.iter().enumerate() {
                scores[class_idx] += self.learning_rate * tree.predict(sample);
            }
        }
        scores
    }

    // ── Feature importance ──────────────────────────────────────

    /// Return impurity-based feature importance.
    ///
    /// The returned vector has `n_features` elements that sum to 1.0. The
    /// importance is the total weighted MSE decrease attributed to each
    /// feature across all split nodes in all trees.
    pub fn feature_importance(&self) -> Vec<f64> {
        self.feature_importance_cache.clone()
    }

    /// Compute permutation importance for regression.
    ///
    /// Measures the increase in MSE when each feature is randomly shuffled.
    /// Higher values indicate more important features.
    ///
    /// * `data` — flat row-major evaluation data
    /// * `n_features` — number of features per sample
    /// * `targets` — true target values
    /// * `n_repeats` — number of permutation repeats per feature
    /// * `seed` — PRNG seed
    pub fn permutation_importance_regression(
        &self,
        data: &[f64],
        n_features: usize,
        targets: &[f64],
        n_repeats: usize,
        seed: u64,
    ) -> Vec<f64> {
        let n_samples = data.len() / n_features;
        let mut rng = LcgRng::new(seed);

        // Baseline MSE
        let baseline_mse = {
            let preds = self.predict_batch(data, n_features);
            preds
                .iter()
                .zip(targets.iter())
                .map(|(&p, &t)| (p - t).powi(2))
                .sum::<f64>()
                / n_samples as f64
        };

        let mut importance = vec![0.0; n_features];
        let mut permuted_data = data.to_vec();

        for feat in 0..n_features {
            let mut total_increase = 0.0;

            for _ in 0..n_repeats {
                // Copy original data
                permuted_data.copy_from_slice(data);

                // Shuffle this feature column
                let mut col_indices: Vec<usize> = (0..n_samples).collect();
                shuffle_indices(&mut rng, &mut col_indices);
                for i in 0..n_samples {
                    permuted_data[i * n_features + feat] =
                        data[col_indices[i] * n_features + feat];
                }

                let preds = self.predict_batch(&permuted_data, n_features);
                let permuted_mse = preds
                    .iter()
                    .zip(targets.iter())
                    .map(|(&p, &t)| (p - t).powi(2))
                    .sum::<f64>()
                    / n_samples as f64;

                total_increase += permuted_mse - baseline_mse;
            }

            importance[feat] = total_increase / n_repeats as f64;
        }

        importance
    }

    /// Compute permutation importance for classification.
    ///
    /// Measures the decrease in accuracy when each feature is randomly shuffled.
    /// Higher values indicate more important features.
    ///
    /// * `data` — flat row-major evaluation data
    /// * `n_features` — number of features per sample
    /// * `labels` — true class labels
    /// * `n_repeats` — number of permutation repeats per feature
    /// * `seed` — PRNG seed
    pub fn permutation_importance_classification(
        &self,
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        n_repeats: usize,
        seed: u64,
    ) -> Vec<f64> {
        let n_samples = data.len() / n_features;
        let mut rng = LcgRng::new(seed);

        // Baseline accuracy
        let baseline_acc = {
            let preds = self.predict_class_batch(data, n_features);
            preds
                .iter()
                .zip(labels.iter())
                .filter(|(&p, &l)| p == l)
                .count() as f64
                / n_samples as f64
        };

        let mut importance = vec![0.0; n_features];
        let mut permuted_data = data.to_vec();

        for feat in 0..n_features {
            let mut total_decrease = 0.0;

            for _ in 0..n_repeats {
                permuted_data.copy_from_slice(data);

                let mut col_indices: Vec<usize> = (0..n_samples).collect();
                shuffle_indices(&mut rng, &mut col_indices);
                for i in 0..n_samples {
                    permuted_data[i * n_features + feat] =
                        data[col_indices[i] * n_features + feat];
                }

                let preds = self.predict_class_batch(&permuted_data, n_features);
                let permuted_acc = preds
                    .iter()
                    .zip(labels.iter())
                    .filter(|(&p, &l)| p == l)
                    .count() as f64
                    / n_samples as f64;

                total_decrease += baseline_acc - permuted_acc;
            }

            importance[feat] = total_decrease / n_repeats as f64;
        }

        importance
    }

    // ── Accessors ───────────────────────────────────────────────

    /// Number of boosting rounds actually used (may be less than configured
    /// if early stopping triggered).
    pub fn n_estimators(&self) -> usize {
        match self.mode {
            Mode::MulticlassClassification => {
                if self.n_classes == 0 {
                    0
                } else {
                    self.trees.len() / self.n_classes
                }
            }
            _ => self.trees.len(),
        }
    }

    /// Number of features expected per sample.
    pub fn n_features(&self) -> usize {
        self.n_features
    }

    /// Whether this model was fit for regression.
    pub fn is_regression(&self) -> bool {
        self.mode == Mode::Regression
    }

    /// Number of classes (0 for regression models).
    pub fn n_classes(&self) -> usize {
        self.n_classes
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a regression tree on a subset of training data.
///
/// `residuals` is indexed by position in `train_idx` (not original sample index).
/// `sub_indices` are indices into `train_idx` / `residuals`.
fn fit_tree_on_subset(
    data: &[f64],
    n_features: usize,
    residuals: &[f64],
    train_idx: &[usize],
    sub_indices: &[usize],
    candidate_features: &[usize],
    max_depth: usize,
    min_samples_leaf: usize,
) -> RegressionTree {
    // Build a view that the regression tree can work with:
    // We assemble flat data and targets aligned to sub_indices.
    let n_sub = sub_indices.len();
    let mut sub_data = Vec::with_capacity(n_sub * n_features);
    let mut sub_targets = Vec::with_capacity(n_sub);

    for &si in sub_indices {
        let orig_idx = train_idx[si];
        let row = &data[orig_idx * n_features..(orig_idx + 1) * n_features];
        sub_data.extend_from_slice(row);
        sub_targets.push(residuals[si]);
    }

    let indices: Vec<usize> = (0..n_sub).collect();
    RegressionTree::fit(
        &sub_data,
        n_features,
        &sub_targets,
        &indices,
        candidate_features,
        max_depth,
        min_samples_leaf,
    )
}

/// Apply Newton-Raphson leaf correction for binary classification.
///
/// For each leaf, replace the mean-residual value with
/// `sum(residuals) / sum(p * (1-p))` where p = sigmoid(F(x)) for the
/// training samples that land in that leaf.
fn apply_binary_leaf_correction(
    tree: &mut RegressionTree,
    data: &[f64],
    n_features: usize,
    f_train: &[f64],
    labels: &[usize],
    train_idx: &[usize],
) {
    // Route each training sample through the tree to find its leaf
    let mut leaf_residuals: std::collections::HashMap<usize, (f64, f64)> =
        std::collections::HashMap::new();

    for (j, &ti) in train_idx.iter().enumerate() {
        let row = &data[ti * n_features..(ti + 1) * n_features];
        let leaf_idx = find_leaf(tree, row);
        let p = sigmoid(f_train[j]);
        let r = labels[ti] as f64 - p;
        let w = (p * (1.0 - p)).max(1e-8);
        let entry = leaf_residuals.entry(leaf_idx).or_insert((0.0, 0.0));
        entry.0 += r;
        entry.1 += w;
    }

    for (leaf_idx, (sum_r, sum_w)) in &leaf_residuals {
        if let RegTreeNode::Leaf { value, .. } = &mut tree.nodes[*leaf_idx] {
            *value = sum_r / sum_w;
        }
    }
}

/// Apply Newton-Raphson leaf correction for multiclass classification.
///
/// For each leaf: `value = (K-1)/K * sum(r) / sum(|r| * (1 - |r|))`
fn apply_multiclass_leaf_correction(
    tree: &mut RegressionTree,
    data: &[f64],
    n_features: usize,
    residuals: &[f64],
    train_idx: &[usize],
    n_classes: usize,
) {
    let k = n_classes as f64;
    let factor = (k - 1.0) / k;

    let mut leaf_sums: std::collections::HashMap<usize, (f64, f64)> =
        std::collections::HashMap::new();

    for (j, &ti) in train_idx.iter().enumerate() {
        let row = &data[ti * n_features..(ti + 1) * n_features];
        let leaf_idx = find_leaf(tree, row);
        let r = residuals[j];
        let abs_r = r.abs();
        let w = (abs_r * (1.0 - abs_r)).max(1e-8);
        let entry = leaf_sums.entry(leaf_idx).or_insert((0.0, 0.0));
        entry.0 += r;
        entry.1 += w;
    }

    for (leaf_idx, (sum_r, sum_w)) in &leaf_sums {
        if let RegTreeNode::Leaf { value, .. } = &mut tree.nodes[*leaf_idx] {
            *value = factor * sum_r / sum_w;
        }
    }
}

/// Find the leaf node index that a sample routes to.
fn find_leaf(tree: &RegressionTree, sample: &[f64]) -> usize {
    let mut idx = 0;
    loop {
        match &tree.nodes[idx] {
            RegTreeNode::Leaf { .. } => return idx,
            RegTreeNode::Split {
                feature_idx,
                threshold,
                left,
                right,
                ..
            } => {
                if sample[*feature_idx] <= *threshold {
                    idx = *left;
                } else {
                    idx = *right;
                }
            }
        }
    }
}

/// Compute impurity-based feature importance across all trees.
fn compute_feature_importance(trees: &[RegressionTree], n_features: usize) -> Vec<f64> {
    let mut importance = vec![0.0; n_features];
    let mut total = 0.0;

    for tree in trees {
        for node in tree.nodes() {
            if let RegTreeNode::Split {
                feature_idx,
                impurity_decrease,
                n_samples,
                ..
            } = node
            {
                let weighted = *impurity_decrease * (*n_samples as f64);
                if *feature_idx < n_features {
                    importance[*feature_idx] += weighted;
                    total += weighted;
                }
            }
        }
    }

    if total > 0.0 {
        for v in &mut importance {
            *v /= total;
        }
    }

    importance
}

/// Split sample indices into training and validation sets.
fn train_val_split(
    rng: &mut LcgRng,
    n_samples: usize,
    val_fraction: f64,
) -> (Vec<usize>, Vec<usize>) {
    let mut indices: Vec<usize> = (0..n_samples).collect();
    shuffle_indices(rng, &mut indices);

    let val_size = ((n_samples as f64 * val_fraction).round() as usize).max(1);
    let val_size = val_size.min(n_samples - 1); // ensure at least 1 training sample
    let val_idx = indices[..val_size].to_vec();
    let train_idx = indices[val_size..].to_vec();
    (train_idx, val_idx)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Regression ──────────────────────────────────────────────

    #[test]
    fn regression_linear_function() {
        // y = 2*x + 1
        let mut data = Vec::new();
        let mut targets = Vec::new();
        for i in 0..50 {
            let x = i as f64 * 0.5;
            data.push(x);
            targets.push(2.0 * x + 1.0);
        }
        let config = GbdtConfig {
            n_estimators: 100,
            learning_rate: 0.1,
            max_depth: 4,
            seed: 42,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
        assert!(model.is_regression());
        assert_eq!(model.n_classes(), 0);

        // Test on a few points
        let pred = model.predict(&[5.0]);
        assert!(
            (pred - 11.0).abs() < 2.0,
            "expected ~11.0, got {}",
            pred
        );

        let pred = model.predict(&[10.0]);
        assert!(
            (pred - 21.0).abs() < 3.0,
            "expected ~21.0, got {}",
            pred
        );
    }

    #[test]
    fn regression_constant_target() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![7.0, 7.0, 7.0, 7.0, 7.0];
        let config = GbdtConfig {
            n_estimators: 10,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
        let pred = model.predict(&[99.0]);
        assert!(
            (pred - 7.0).abs() < 0.01,
            "expected ~7.0, got {}",
            pred
        );
    }

    #[test]
    fn regression_predict_batch_consistency() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let targets = vec![10.0, 20.0, 30.0];
        let config = GbdtConfig {
            n_estimators: 20,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 2, &targets, &config).unwrap();

        let batch = model.predict_batch(&data, 2);
        for i in 0..3 {
            let single = model.predict(&data[i * 2..(i + 1) * 2]);
            assert!(
                (batch[i] - single).abs() < 1e-12,
                "batch[{}]={} != single={}",
                i,
                batch[i],
                single
            );
        }
    }

    #[test]
    fn regression_feature_importance_informative() {
        // Feature 0 determines target, feature 1 is noise
        let mut data = Vec::new();
        let mut targets = Vec::new();
        for i in 0..40 {
            let x = i as f64;
            data.push(x);
            data.push(42.0); // constant noise feature
            targets.push(x * 3.0);
        }
        let config = GbdtConfig {
            n_estimators: 50,
            max_depth: 4,
            seed: 42,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 2, &targets, &config).unwrap();
        let fi = model.feature_importance();
        assert_eq!(fi.len(), 2);
        assert!(
            fi[0] > fi[1],
            "informative feature should have higher importance: {:?}",
            fi
        );
    }

    #[test]
    fn regression_subsample() {
        let data: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let targets: Vec<f64> = data.iter().map(|&x| x * 2.0).collect();
        let config = GbdtConfig {
            n_estimators: 50,
            subsample: 0.7,
            seed: 42,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
        let pred = model.predict(&[10.0]);
        assert!(
            (pred - 20.0).abs() < 5.0,
            "subsample pred {} too far from 20.0",
            pred
        );
    }

    #[test]
    fn regression_early_stopping() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let targets: Vec<f64> = data.iter().map(|&x| x * 2.0 + 1.0).collect();
        let config = GbdtConfig {
            n_estimators: 500,
            learning_rate: 0.3,
            max_depth: 3,
            early_stopping_rounds: Some(5),
            validation_fraction: 0.2,
            seed: 42,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
        // Should have stopped before 500 rounds
        assert!(
            model.n_estimators() < 500,
            "expected early stopping, got {} rounds",
            model.n_estimators()
        );
        assert!(model.n_estimators() > 0);
    }

    #[test]
    fn regression_learning_rate_effect() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let targets: Vec<f64> = data.iter().map(|&x| x * 2.0).collect();

        let config_slow = GbdtConfig {
            n_estimators: 10,
            learning_rate: 0.01,
            seed: 42,
            ..Default::default()
        };
        let config_fast = GbdtConfig {
            n_estimators: 10,
            learning_rate: 0.5,
            seed: 42,
            ..Default::default()
        };
        let model_slow =
            GradientBoostedTrees::fit_regression(&data, 1, &targets, &config_slow).unwrap();
        let model_fast =
            GradientBoostedTrees::fit_regression(&data, 1, &targets, &config_fast).unwrap();

        // Higher learning rate should converge faster (lower training error) in same rounds
        let err_slow: f64 = data
            .iter()
            .zip(targets.iter())
            .map(|(&x, &t)| (model_slow.predict(&[x]) - t).powi(2))
            .sum::<f64>();
        let err_fast: f64 = data
            .iter()
            .zip(targets.iter())
            .map(|(&x, &t)| (model_fast.predict(&[x]) - t).powi(2))
            .sum::<f64>();
        assert!(
            err_fast < err_slow,
            "fast lr error {} should be < slow lr error {}",
            err_fast,
            err_slow
        );
    }

    #[test]
    fn regression_single_sample() {
        let data = vec![5.0];
        let targets = vec![10.0];
        let config = GbdtConfig {
            n_estimators: 10,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
        let pred = model.predict(&[5.0]);
        assert!(
            (pred - 10.0).abs() < 0.5,
            "expected ~10.0, got {}",
            pred
        );
    }

    // ── Binary Classification ───────────────────────────────────

    #[test]
    fn binary_linearly_separable() {
        // Class 0: x < 5, Class 1: x >= 5
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..20 {
            let x = i as f64;
            data.push(x);
            data.push(0.0); // noise feature
            labels.push(if x < 10.0 { 0 } else { 1 });
        }
        let config = GbdtConfig {
            n_estimators: 50,
            max_depth: 3,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();
        assert!(!model.is_regression());
        assert_eq!(model.n_classes(), 2);

        let preds = model.predict_class_batch(&data, 2);
        let correct = preds.iter().zip(labels.iter()).filter(|(&p, &l)| p == l).count();
        let accuracy = correct as f64 / labels.len() as f64;
        assert!(
            accuracy > 0.9,
            "accuracy {:.2} too low",
            accuracy
        );
    }

    #[test]
    fn binary_predict_class_proba_consistency() {
        let data = vec![0.0, 1.0, 2.0, 10.0, 11.0, 12.0];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let config = GbdtConfig {
            n_estimators: 30,
            max_depth: 3,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 1, &labels, &config).unwrap();

        for i in 0..6 {
            let sample = &data[i..i + 1];
            let cls = model.predict_class(sample);
            let proba = model.predict_proba(sample);
            assert_eq!(proba.len(), 2);
            // Class should match argmax of proba
            let argmax = if proba[1] >= proba[0] { 1 } else { 0 };
            assert_eq!(cls, argmax, "class {} != argmax {} for sample {}", cls, argmax, i);
        }
    }

    #[test]
    fn binary_batch_consistency() {
        let data = vec![0.0, 1.0, 10.0, 11.0];
        let labels = vec![0, 0, 1, 1];
        let config = GbdtConfig {
            n_estimators: 30,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 1, &labels, &config).unwrap();

        let batch = model.predict_class_batch(&data, 1);
        for i in 0..4 {
            let single = model.predict_class(&data[i..i + 1]);
            assert_eq!(batch[i], single);
        }
    }

    #[test]
    fn binary_proba_sums_to_one() {
        let data = vec![0.0, 5.0, 10.0, 15.0];
        let labels = vec![0, 0, 1, 1];
        let config = GbdtConfig {
            n_estimators: 20,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 1, &labels, &config).unwrap();

        for i in 0..4 {
            let proba = model.predict_proba(&data[i..i + 1]);
            let sum: f64 = proba.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "proba sum {} != 1.0 for sample {}",
                sum,
                i
            );
        }
    }

    #[test]
    fn binary_feature_importance_sums_to_one() {
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..30 {
            data.push(i as f64);
            data.push(0.0);
            labels.push(if i < 15 { 0 } else { 1 });
        }
        let config = GbdtConfig {
            n_estimators: 30,
            max_depth: 3,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();
        let fi = model.feature_importance();
        assert_eq!(fi.len(), 2);
        let sum: f64 = fi.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "fi sum {} != 1.0",
            sum
        );
    }

    #[test]
    fn binary_early_stopping() {
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..50 {
            data.push(i as f64);
            labels.push(if i < 25 { 0 } else { 1 });
        }
        let config = GbdtConfig {
            n_estimators: 500,
            learning_rate: 0.3,
            max_depth: 3,
            early_stopping_rounds: Some(5),
            validation_fraction: 0.2,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 1, &labels, &config).unwrap();
        assert!(
            model.n_estimators() < 500,
            "expected early stopping, got {} rounds",
            model.n_estimators()
        );
    }

    #[test]
    fn binary_is_classification_checks() {
        let data = vec![0.0, 10.0];
        let labels = vec![0, 1];
        let config = GbdtConfig {
            n_estimators: 5,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 1, &labels, &config).unwrap();
        assert!(!model.is_regression());
        assert_eq!(model.n_classes(), 2);
        assert_eq!(model.n_features(), 1);
    }

    #[test]
    fn binary_single_class_edge() {
        // All samples belong to class 0
        let data = vec![1.0, 2.0, 3.0];
        let labels = vec![0, 0, 0];
        let config = GbdtConfig::default();
        let model =
            GradientBoostedTrees::fit_classification(&data, 1, &labels, &config).unwrap();
        assert_eq!(model.predict_class(&[1.0]), 0);
    }

    // ── Multiclass Classification ───────────────────────────────

    #[test]
    fn multiclass_three_clusters() {
        // Three well-separated clusters
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..15 {
            let offset = i as f64 * 0.1;
            data.push(0.0 + offset);
            data.push(0.0 + offset);
            labels.push(0);
        }
        for i in 0..15 {
            let offset = i as f64 * 0.1;
            data.push(10.0 + offset);
            data.push(10.0 + offset);
            labels.push(1);
        }
        for i in 0..15 {
            let offset = i as f64 * 0.1;
            data.push(20.0 + offset);
            data.push(20.0 + offset);
            labels.push(2);
        }

        let config = GbdtConfig {
            n_estimators: 50,
            max_depth: 4,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();
        assert_eq!(model.n_classes(), 3);

        let preds = model.predict_class_batch(&data, 2);
        let correct = preds.iter().zip(labels.iter()).filter(|(&p, &l)| p == l).count();
        let accuracy = correct as f64 / labels.len() as f64;
        assert!(
            accuracy > 0.9,
            "multiclass accuracy {:.2} too low",
            accuracy
        );
    }

    #[test]
    fn multiclass_proba_sums_to_one() {
        let data = vec![0.0, 0.0, 10.0, 10.0, 20.0, 20.0];
        let labels = vec![0, 1, 2];
        let config = GbdtConfig {
            n_estimators: 20,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();

        for i in 0..3 {
            let proba = model.predict_proba(&data[i * 2..(i + 1) * 2]);
            assert_eq!(proba.len(), 3);
            let sum: f64 = proba.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "proba sum {} != 1.0 for sample {}",
                sum,
                i
            );
        }
    }

    #[test]
    fn multiclass_class_matches_argmax_proba() {
        let data = vec![0.0, 0.0, 10.0, 10.0, 20.0, 20.0];
        let labels = vec![0, 1, 2];
        let config = GbdtConfig {
            n_estimators: 30,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();

        for i in 0..3 {
            let sample = &data[i * 2..(i + 1) * 2];
            let cls = model.predict_class(sample);
            let proba = model.predict_proba(sample);
            let argmax = proba
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(i, _)| i)
                .unwrap();
            assert_eq!(cls, argmax);
        }
    }

    #[test]
    fn multiclass_proba_batch_consistency() {
        let data = vec![0.0, 0.0, 10.0, 10.0, 20.0, 20.0];
        let labels = vec![0, 1, 2];
        let config = GbdtConfig {
            n_estimators: 20,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();

        let batch = model.predict_proba_batch(&data, 2);
        assert_eq!(batch.len(), 3 * 3); // 3 samples * 3 classes

        for i in 0..3 {
            let single = model.predict_proba(&data[i * 2..(i + 1) * 2]);
            for k in 0..3 {
                assert!(
                    (batch[i * 3 + k] - single[k]).abs() < 1e-12,
                    "mismatch at sample {} class {}",
                    i,
                    k
                );
            }
        }
    }

    #[test]
    fn multiclass_feature_importance_length() {
        let data = vec![0.0, 0.0, 10.0, 10.0, 20.0, 20.0];
        let labels = vec![0, 1, 2];
        let config = GbdtConfig {
            n_estimators: 20,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();
        let fi = model.feature_importance();
        assert_eq!(fi.len(), 2);
        let sum: f64 = fi.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10 || sum == 0.0,
            "fi sum {} unexpected",
            sum
        );
    }

    // ── Error Handling ──────────────────────────────────────────

    #[test]
    fn error_empty_data() {
        let config = GbdtConfig::default();
        assert!(GradientBoostedTrees::fit_regression(&[], 1, &[], &config).is_err());
        assert!(GradientBoostedTrees::fit_classification(&[], 1, &[], &config).is_err());
    }

    #[test]
    fn error_dimension_mismatch() {
        let data = vec![1.0, 2.0, 3.0]; // 3 elements, not divisible by 2
        let targets = vec![1.0];
        let config = GbdtConfig::default();
        assert!(GradientBoostedTrees::fit_regression(&data, 2, &targets, &config).is_err());
    }

    #[test]
    fn error_labels_length_mismatch() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let labels = vec![0]; // 2 samples but 1 label
        let config = GbdtConfig::default();
        assert!(GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).is_err());
    }

    #[test]
    fn error_invalid_config() {
        let data = vec![1.0];
        let targets = vec![1.0];

        // n_estimators = 0
        let config = GbdtConfig {
            n_estimators: 0,
            ..Default::default()
        };
        assert!(GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).is_err());

        // learning_rate <= 0
        let config = GbdtConfig {
            learning_rate: 0.0,
            ..Default::default()
        };
        assert!(GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).is_err());

        // subsample out of range
        let config = GbdtConfig {
            subsample: 0.0,
            ..Default::default()
        };
        assert!(GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).is_err());

        let config = GbdtConfig {
            subsample: 1.5,
            ..Default::default()
        };
        assert!(GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).is_err());
    }

    // ── Permutation Importance ──────────────────────────────────

    #[test]
    fn permutation_importance_regression_informative() {
        // Feature 0 determines target, feature 1 is constant noise
        let mut data = Vec::new();
        let mut targets = Vec::new();
        for i in 0..40 {
            let x = i as f64;
            data.push(x);
            data.push(42.0);
            targets.push(x * 3.0);
        }
        let config = GbdtConfig {
            n_estimators: 50,
            max_depth: 4,
            seed: 42,
            ..Default::default()
        };
        let model = GradientBoostedTrees::fit_regression(&data, 2, &targets, &config).unwrap();
        let pi = model.permutation_importance_regression(&data, 2, &targets, 5, 123);
        assert!(
            pi[0] > pi[1],
            "informative feature should have higher perm importance: {:?}",
            pi
        );
    }

    #[test]
    fn permutation_importance_classification_informative() {
        // Feature 0 determines class, feature 1 is noise
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..40 {
            let x = i as f64;
            data.push(x);
            data.push(42.0);
            labels.push(if x < 20.0 { 0 } else { 1 });
        }
        let config = GbdtConfig {
            n_estimators: 50,
            max_depth: 3,
            seed: 42,
            ..Default::default()
        };
        let model =
            GradientBoostedTrees::fit_classification(&data, 2, &labels, &config).unwrap();
        let pi = model.permutation_importance_classification(&data, 2, &labels, 5, 123);
        assert!(
            pi[0] > pi[1],
            "informative feature should have higher perm importance: {:?}",
            pi
        );
    }

    // ── Determinism ─────────────────────────────────────────────

    #[test]
    fn deterministic_with_same_seed() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let targets: Vec<f64> = data.iter().map(|&x| x * 2.0).collect();
        let config = GbdtConfig {
            n_estimators: 30,
            seed: 42,
            ..Default::default()
        };

        let model1 = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();
        let model2 = GradientBoostedTrees::fit_regression(&data, 1, &targets, &config).unwrap();

        let preds1 = model1.predict_batch(&data, 1);
        let preds2 = model2.predict_batch(&data, 1);
        assert_eq!(preds1.len(), preds2.len());
        for (a, b) in preds1.iter().zip(preds2.iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "predictions differ: {} vs {}",
                a,
                b
            );
        }
    }
}
