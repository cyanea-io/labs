//! Random forest classifier.
//!
//! Implements a bagged ensemble of [`DecisionTree`](crate::tree::DecisionTree)
//! classifiers with bootstrap sampling and feature bagging. Suitable for
//! small-to-medium biological datasets (expression matrices, k-mer profiles).
//!
//! Data is flat row-major `&[f64]` with an `n_features` parameter, consistent
//! with the rest of the cyanea-ml crate.

use cyanea_core::{CyaneaError, Result};

use crate::tree::{DecisionTree, TreeNode};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for random forest training.
#[derive(Debug, Clone)]
pub struct RandomForestConfig {
    /// Number of trees in the ensemble.
    pub n_trees: usize,
    /// Maximum depth per tree.
    pub max_depth: usize,
    /// Number of features to consider at each split. `None` defaults to
    /// `sqrt(n_features)`.
    pub max_features: Option<usize>,
    /// Random seed for reproducibility.
    pub seed: u64,
}

impl Default for RandomForestConfig {
    fn default() -> Self {
        Self {
            n_trees: 100,
            max_depth: 10,
            max_features: None,
            seed: 42,
        }
    }
}

// ---------------------------------------------------------------------------
// RandomForest
// ---------------------------------------------------------------------------

/// A random forest classifier (ensemble of decision trees).
#[derive(Debug, Clone)]
pub struct RandomForest {
    trees: Vec<DecisionTree>,
    n_classes: usize,
}

impl RandomForest {
    /// Fit a random forest on flat row-major data.
    ///
    /// * `data` — flat row-major `n_samples x n_features`
    /// * `n_features` — number of features per sample
    /// * `labels` — class label for each sample
    /// * `config` — forest hyper-parameters
    ///
    /// # Errors
    ///
    /// Returns an error if the data is empty, dimensions are inconsistent,
    /// or `n_trees` is 0.
    pub fn fit(
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        config: &RandomForestConfig,
    ) -> Result<Self> {
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
        let n_samples = data.len() / n_features;
        if labels.len() != n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "labels length {} != n_samples {}",
                labels.len(),
                n_samples
            )));
        }
        if config.n_trees == 0 {
            return Err(CyaneaError::InvalidInput("n_trees must be > 0".into()));
        }

        let n_classes = labels.iter().copied().max().map_or(0, |m| m + 1);
        let max_features = config
            .max_features
            .unwrap_or_else(|| isqrt(n_features).max(1));

        let mut rng = LcgRng::new(config.seed);
        let mut trees = Vec::with_capacity(config.n_trees);

        for _ in 0..config.n_trees {
            // Bootstrap sample: n_samples drawn with replacement
            let sample_indices: Vec<usize> = (0..n_samples)
                .map(|_| rng.next_bounded(n_samples as u64) as usize)
                .collect();

            // Feature bagging: random subset of features
            let candidate_features = random_feature_subset(&mut rng, n_features, max_features);

            let tree = DecisionTree::fit_with_features(
                data,
                n_features,
                labels,
                &sample_indices,
                &candidate_features,
                config.max_depth,
            )?;
            trees.push(tree);
        }

        Ok(Self { trees, n_classes })
    }

    /// Predict the class label for a single sample via majority vote.
    pub fn predict(&self, sample: &[f64]) -> usize {
        let mut votes = vec![0usize; self.n_classes.max(1)];
        for tree in &self.trees {
            let pred = tree.predict(sample);
            if pred < votes.len() {
                votes[pred] += 1;
            }
        }
        votes
            .iter()
            .enumerate()
            .max_by_key(|&(_, &c)| c)
            .map(|(cls, _)| cls)
            .unwrap_or(0)
    }

    /// Predict class labels for multiple samples.
    ///
    /// `data` is flat row-major with `n_features` columns.
    pub fn predict_batch(&self, data: &[f64], n_features: usize) -> Vec<usize> {
        let n_samples = data.len() / n_features;
        (0..n_samples)
            .map(|i| {
                let row = &data[i * n_features..(i + 1) * n_features];
                self.predict(row)
            })
            .collect()
    }

    /// Compute feature importance as the normalized frequency of each
    /// feature being used for splits across all trees.
    ///
    /// The returned vector has `n_features` elements that sum to 1.0.
    pub fn feature_importance(&self, n_features: usize) -> Vec<f64> {
        let mut counts = vec![0usize; n_features];
        let mut total = 0usize;

        for tree in &self.trees {
            for node in tree.nodes() {
                if let TreeNode::Split { feature_idx, .. } = node {
                    if *feature_idx < n_features {
                        counts[*feature_idx] += 1;
                        total += 1;
                    }
                }
            }
        }

        if total == 0 {
            return vec![0.0; n_features];
        }

        counts
            .iter()
            .map(|&c| c as f64 / total as f64)
            .collect()
    }

    /// Number of trees in the forest.
    pub fn n_trees(&self) -> usize {
        self.trees.len()
    }

    /// Number of classes discovered during fitting.
    pub fn n_classes(&self) -> usize {
        self.n_classes
    }
}

// ---------------------------------------------------------------------------
// PRNG and helpers
// ---------------------------------------------------------------------------

/// Simple LCG PRNG (linear congruential generator).
///
/// Uses the same constants as glibc: multiplier 6364136223846793005,
/// increment 1442695040888963407.
struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1), // avoid zero state
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

/// Integer square root (floor).
fn isqrt(n: usize) -> usize {
    (n as f64).sqrt() as usize
}

/// Select `count` distinct feature indices from `0..n_features` at random.
fn random_feature_subset(rng: &mut LcgRng, n_features: usize, count: usize) -> Vec<usize> {
    let count = count.min(n_features);
    if count == n_features {
        return (0..n_features).collect();
    }

    // Fisher-Yates partial shuffle
    let mut pool: Vec<usize> = (0..n_features).collect();
    for i in 0..count {
        let j = i + (rng.next_bounded((n_features - i) as u64) as usize);
        pool.swap(i, j);
    }
    pool.truncate(count);
    pool
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Create synthetic Iris-like data: 3 classes, 4 features, 10 samples each.
    fn iris_like_data() -> (Vec<f64>, Vec<usize>) {
        let mut data = Vec::new();
        let mut labels = Vec::new();

        // Class 0: low values
        for i in 0..10 {
            let offset = i as f64 * 0.1;
            data.extend_from_slice(&[1.0 + offset, 0.5 + offset, 0.2 + offset, 0.1 + offset]);
            labels.push(0);
        }
        // Class 1: medium values
        for i in 0..10 {
            let offset = i as f64 * 0.1;
            data.extend_from_slice(&[5.0 + offset, 3.0 + offset, 3.5 + offset, 1.0 + offset]);
            labels.push(1);
        }
        // Class 2: high values
        for i in 0..10 {
            let offset = i as f64 * 0.1;
            data.extend_from_slice(&[7.0 + offset, 4.0 + offset, 6.0 + offset, 2.0 + offset]);
            labels.push(2);
        }

        (data, labels)
    }

    #[test]
    fn fit_and_predict_iris_like() {
        let (data, labels) = iris_like_data();
        let config = RandomForestConfig {
            n_trees: 20,
            max_depth: 5,
            seed: 42,
            ..Default::default()
        };
        let forest = RandomForest::fit(&data, 4, &labels, &config).unwrap();

        // Should correctly classify training data (well-separated clusters)
        let preds = forest.predict_batch(&data, 4);
        let correct = preds
            .iter()
            .zip(labels.iter())
            .filter(|(&p, &l)| p == l)
            .count();
        let accuracy = correct as f64 / labels.len() as f64;
        assert!(
            accuracy > 0.9,
            "accuracy {:.2} too low on training data",
            accuracy
        );
    }

    #[test]
    fn single_sample_prediction() {
        let (data, labels) = iris_like_data();
        let config = RandomForestConfig {
            n_trees: 10,
            max_depth: 5,
            seed: 123,
            ..Default::default()
        };
        let forest = RandomForest::fit(&data, 4, &labels, &config).unwrap();

        // Predict a clearly class-0 sample
        let pred = forest.predict(&[1.0, 0.5, 0.2, 0.1]);
        assert_eq!(pred, 0);
    }

    #[test]
    fn feature_importance_sums_to_one() {
        let (data, labels) = iris_like_data();
        let config = RandomForestConfig {
            n_trees: 20,
            max_depth: 5,
            seed: 42,
            ..Default::default()
        };
        let forest = RandomForest::fit(&data, 4, &labels, &config).unwrap();

        let importance = forest.feature_importance(4);
        assert_eq!(importance.len(), 4);
        let total: f64 = importance.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-10,
            "importance sum {} != 1.0",
            total
        );
        // All values should be non-negative
        for &v in &importance {
            assert!(v >= 0.0);
        }
    }

    #[test]
    fn config_options_work() {
        let data = vec![0.0, 1.0, 10.0, 11.0, 0.1, 1.1, 10.1, 11.1];
        let labels = vec![0, 1, 0, 1];
        let config = RandomForestConfig {
            n_trees: 5,
            max_depth: 3,
            max_features: Some(1),
            seed: 99,
        };
        let forest = RandomForest::fit(&data, 2, &labels, &config).unwrap();
        assert_eq!(forest.n_trees(), 5);
        assert_eq!(forest.n_classes(), 2);
    }

    #[test]
    fn empty_data_error() {
        let config = RandomForestConfig::default();
        assert!(RandomForest::fit(&[], 2, &[], &config).is_err());
    }

    #[test]
    fn zero_trees_error() {
        let config = RandomForestConfig {
            n_trees: 0,
            ..Default::default()
        };
        assert!(RandomForest::fit(&[1.0, 2.0], 2, &[0], &config).is_err());
    }

    #[test]
    fn deterministic_with_seed() {
        let (data, labels) = iris_like_data();
        let config = RandomForestConfig {
            n_trees: 10,
            max_depth: 5,
            seed: 42,
            ..Default::default()
        };
        let forest1 = RandomForest::fit(&data, 4, &labels, &config).unwrap();
        let forest2 = RandomForest::fit(&data, 4, &labels, &config).unwrap();

        let preds1 = forest1.predict_batch(&data, 4);
        let preds2 = forest2.predict_batch(&data, 4);
        assert_eq!(preds1, preds2, "same seed should give same predictions");
    }

    #[test]
    fn binary_classification() {
        // Simple binary split on feature 0
        let data = vec![
            0.0, 0.0,
            1.0, 0.0,
            2.0, 0.0,
            10.0, 0.0,
            11.0, 0.0,
            12.0, 0.0,
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let config = RandomForestConfig {
            n_trees: 15,
            max_depth: 5,
            seed: 42,
            ..Default::default()
        };
        let forest = RandomForest::fit(&data, 2, &labels, &config).unwrap();
        assert_eq!(forest.predict(&[0.5, 0.0]), 0);
        assert_eq!(forest.predict(&[11.5, 0.0]), 1);
    }

    #[test]
    fn feature_importance_with_informative_feature() {
        // Feature 0 is informative, feature 1 is noise
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..20 {
            let x = i as f64;
            data.extend_from_slice(&[x, 0.0]); // feature 1 constant
            labels.push(if x < 10.0 { 0 } else { 1 });
        }

        let config = RandomForestConfig {
            n_trees: 30,
            max_depth: 5,
            max_features: Some(2), // consider all features
            seed: 42,
        };
        let forest = RandomForest::fit(&data, 2, &labels, &config).unwrap();
        let importance = forest.feature_importance(2);

        // Feature 0 should have higher importance since feature 1 is constant
        assert!(
            importance[0] > importance[1],
            "informative feature should have higher importance: {:?}",
            importance
        );
    }
}
