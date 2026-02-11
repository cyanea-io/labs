//! Decision tree classifier using Gini impurity.
//!
//! Provides a simple, dependency-free CART-style decision tree suitable for
//! classification tasks on biological data (e.g., expression profiles, k-mer
//! feature vectors).
//!
//! Data is flat row-major `&[f64]` with an `n_features` parameter, consistent
//! with the rest of the cyanea-ml crate.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Tree node representation (arena-allocated)
// ---------------------------------------------------------------------------

/// A single node in the decision tree.
#[derive(Debug, Clone)]
pub(crate) enum TreeNode {
    /// Internal split node.
    Split {
        feature_idx: usize,
        threshold: f64,
        left: usize,  // index into arena
        right: usize, // index into arena
    },
    /// Terminal leaf node.
    Leaf {
        class: usize,
    },
}

// ---------------------------------------------------------------------------
// DecisionTree
// ---------------------------------------------------------------------------

/// A decision tree classifier.
///
/// The tree is stored as a flat arena of [`TreeNode`] values, with index 0
/// as the root. This avoids recursive `Box` allocations and keeps the
/// structure cache-friendly.
#[derive(Debug, Clone)]
pub struct DecisionTree {
    nodes: Vec<TreeNode>,
}

impl DecisionTree {
    /// Fit a decision tree on flat row-major data.
    ///
    /// * `data` — flat row-major `n_samples x n_features`
    /// * `n_features` — number of features per sample
    /// * `labels` — class label for each sample (`0..n_classes`)
    /// * `max_depth` — maximum tree depth (0 = only root leaf)
    ///
    /// # Errors
    ///
    /// Returns an error if the data is empty, dimensions are inconsistent,
    /// or labels are empty.
    pub fn fit(
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        max_depth: usize,
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

        let indices: Vec<usize> = (0..n_samples).collect();
        let all_features: Vec<usize> = (0..n_features).collect();
        let mut nodes = Vec::new();

        build_tree(
            data,
            n_features,
            labels,
            &indices,
            &all_features,
            max_depth,
            0,
            &mut nodes,
        );

        Ok(Self { nodes })
    }

    /// Fit a decision tree considering only a subset of features at each
    /// split. Used internally by [`RandomForest`](super::forest::RandomForest).
    pub(crate) fn fit_with_features(
        data: &[f64],
        n_features: usize,
        labels: &[usize],
        sample_indices: &[usize],
        candidate_features: &[usize],
        max_depth: usize,
    ) -> Result<Self> {
        if sample_indices.is_empty() {
            return Err(CyaneaError::InvalidInput("empty sample set".into()));
        }
        let mut nodes = Vec::new();
        build_tree(
            data,
            n_features,
            labels,
            sample_indices,
            candidate_features,
            max_depth,
            0,
            &mut nodes,
        );
        Ok(Self { nodes })
    }

    /// Predict the class label for a single sample.
    ///
    /// `sample` must have exactly `n_features` elements (the same value used
    /// during fitting).
    pub fn predict(&self, sample: &[f64]) -> usize {
        let mut idx = 0;
        loop {
            match &self.nodes[idx] {
                TreeNode::Leaf { class } => return *class,
                TreeNode::Split {
                    feature_idx,
                    threshold,
                    left,
                    right,
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

    /// Return a reference to the internal node arena (used by forest for
    /// feature importance).
    pub(crate) fn nodes(&self) -> &[TreeNode] {
        &self.nodes
    }
}

// ---------------------------------------------------------------------------
// Tree building helpers
// ---------------------------------------------------------------------------

/// Recursively build the tree, returning the arena index of the created node.
fn build_tree(
    data: &[f64],
    n_features: usize,
    labels: &[usize],
    indices: &[usize],
    candidate_features: &[usize],
    max_depth: usize,
    depth: usize,
    nodes: &mut Vec<TreeNode>,
) -> usize {
    // Majority class for this subset
    let majority = majority_class(labels, indices);

    // Stop conditions: max depth, too few samples, or pure node
    if depth >= max_depth || indices.len() < 2 || is_pure(labels, indices) {
        let idx = nodes.len();
        nodes.push(TreeNode::Leaf { class: majority });
        return idx;
    }

    // Find best split across candidate features
    if let Some((best_feature, best_threshold)) =
        find_best_split(data, n_features, labels, indices, candidate_features)
    {
        // Partition indices
        let (left_indices, right_indices) =
            partition(data, n_features, indices, best_feature, best_threshold);

        // If split produces an empty child, make a leaf
        if left_indices.is_empty() || right_indices.is_empty() {
            let idx = nodes.len();
            nodes.push(TreeNode::Leaf { class: majority });
            return idx;
        }

        // Reserve a slot for this split node
        let node_idx = nodes.len();
        nodes.push(TreeNode::Leaf { class: 0 }); // placeholder

        let left_child = build_tree(
            data,
            n_features,
            labels,
            &left_indices,
            candidate_features,
            max_depth,
            depth + 1,
            nodes,
        );
        let right_child = build_tree(
            data,
            n_features,
            labels,
            &right_indices,
            candidate_features,
            max_depth,
            depth + 1,
            nodes,
        );

        nodes[node_idx] = TreeNode::Split {
            feature_idx: best_feature,
            threshold: best_threshold,
            left: left_child,
            right: right_child,
        };

        node_idx
    } else {
        // No valid split found
        let idx = nodes.len();
        nodes.push(TreeNode::Leaf { class: majority });
        idx
    }
}

/// Find the best (feature, threshold) split by minimizing weighted Gini
/// impurity.
fn find_best_split(
    data: &[f64],
    n_features: usize,
    labels: &[usize],
    indices: &[usize],
    candidate_features: &[usize],
) -> Option<(usize, f64)> {
    let n = indices.len();
    let parent_gini = gini_impurity(labels, indices);

    let mut best_gain = 0.0;
    let mut best_feature = 0;
    let mut best_threshold = 0.0;
    let mut found = false;

    for &feat in candidate_features {
        // Collect and sort unique values for this feature
        let mut values: Vec<f64> = indices
            .iter()
            .map(|&i| data[i * n_features + feat])
            .collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        values.dedup();

        if values.len() < 2 {
            continue;
        }

        // Try midpoints between consecutive unique values
        for w in values.windows(2) {
            let threshold = (w[0] + w[1]) / 2.0;

            // Count left/right and compute Gini
            let mut left_indices_buf = Vec::new();
            let mut right_indices_buf = Vec::new();
            for &i in indices {
                if data[i * n_features + feat] <= threshold {
                    left_indices_buf.push(i);
                } else {
                    right_indices_buf.push(i);
                }
            }

            if left_indices_buf.is_empty() || right_indices_buf.is_empty() {
                continue;
            }

            let left_gini = gini_impurity(labels, &left_indices_buf);
            let right_gini = gini_impurity(labels, &right_indices_buf);
            let weighted = (left_indices_buf.len() as f64 * left_gini
                + right_indices_buf.len() as f64 * right_gini)
                / n as f64;
            let gain = parent_gini - weighted;

            if gain > best_gain {
                best_gain = gain;
                best_feature = feat;
                best_threshold = threshold;
                found = true;
            }
        }
    }

    if found {
        Some((best_feature, best_threshold))
    } else {
        None
    }
}

/// Gini impurity for a subset of labels.
fn gini_impurity(labels: &[usize], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }
    // Count classes
    let max_label = indices.iter().map(|&i| labels[i]).max().unwrap_or(0);
    let mut counts = vec![0usize; max_label + 1];
    for &i in indices {
        counts[labels[i]] += 1;
    }
    let n = indices.len() as f64;
    let mut impurity = 1.0;
    for &c in &counts {
        let p = c as f64 / n;
        impurity -= p * p;
    }
    impurity
}

/// Check if all samples in the subset belong to the same class.
fn is_pure(labels: &[usize], indices: &[usize]) -> bool {
    if indices.is_empty() {
        return true;
    }
    let first = labels[indices[0]];
    indices.iter().all(|&i| labels[i] == first)
}

/// Return the most common class label among the given indices.
fn majority_class(labels: &[usize], indices: &[usize]) -> usize {
    if indices.is_empty() {
        return 0;
    }
    let max_label = indices.iter().map(|&i| labels[i]).max().unwrap_or(0);
    let mut counts = vec![0usize; max_label + 1];
    for &i in indices {
        counts[labels[i]] += 1;
    }
    counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, &c)| c)
        .map(|(cls, _)| cls)
        .unwrap_or(0)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_linearly_separable() {
        // 2D data: class 0 = low x, class 1 = high x
        let data = vec![
            0.0, 0.0,
            1.0, 0.0,
            2.0, 0.0,
            10.0, 0.0,
            11.0, 0.0,
            12.0, 0.0,
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let tree = DecisionTree::fit(&data, 2, &labels, 10).unwrap();

        // Should predict training data correctly
        for i in 0..6 {
            let row = &data[i * 2..(i + 1) * 2];
            assert_eq!(tree.predict(row), labels[i], "mismatch at sample {}", i);
        }
    }

    #[test]
    fn predict_batch_works() {
        let data = vec![
            0.0, 1.0,
            10.0, 1.0,
        ];
        let labels = vec![0, 1];
        let tree = DecisionTree::fit(&data, 2, &labels, 5).unwrap();

        let test_data = vec![0.5, 1.0, 9.5, 1.0];
        let preds = tree.predict_batch(&test_data, 2);
        assert_eq!(preds, vec![0, 1]);
    }

    #[test]
    fn pure_node_stops() {
        // All same class — tree should be a single leaf
        let data = vec![0.0, 1.0, 2.0, 3.0];
        let labels = vec![0, 0, 0, 0];
        let tree = DecisionTree::fit(&data, 1, &labels, 10).unwrap();
        assert_eq!(tree.nodes.len(), 1);
        assert_eq!(tree.predict(&[999.0]), 0);
    }

    #[test]
    fn single_sample() {
        let data = vec![5.0, 3.0];
        let labels = vec![2];
        let tree = DecisionTree::fit(&data, 2, &labels, 10).unwrap();
        assert_eq!(tree.predict(&[5.0, 3.0]), 2);
    }

    #[test]
    fn max_depth_zero() {
        // max_depth=0 should produce only a leaf (majority vote)
        let data = vec![0.0, 1.0, 10.0, 11.0];
        let labels = vec![0, 0, 1, 1];
        let tree = DecisionTree::fit(&data, 1, &labels, 0).unwrap();
        assert_eq!(tree.nodes.len(), 1);
        // With a tie, majority_class picks the first max
        let pred = tree.predict(&[5.0]);
        assert!(pred == 0 || pred == 1);
    }

    #[test]
    fn empty_data_error() {
        assert!(DecisionTree::fit(&[], 2, &[], 10).is_err());
    }

    #[test]
    fn dimension_mismatch_error() {
        let data = vec![1.0, 2.0, 3.0];
        let labels = vec![0];
        assert!(DecisionTree::fit(&data, 2, &labels, 10).is_err());
    }

    #[test]
    fn labels_length_mismatch() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let labels = vec![0]; // 2 samples but 1 label
        assert!(DecisionTree::fit(&data, 2, &labels, 10).is_err());
    }

    #[test]
    fn multiclass() {
        // Three classes, separated on a single feature
        let data = vec![0.0, 5.0, 10.0, 0.1, 5.1, 10.1];
        let labels = vec![0, 1, 2, 0, 1, 2];
        let tree = DecisionTree::fit(&data, 1, &labels, 10).unwrap();
        assert_eq!(tree.predict(&[0.05]), 0);
        assert_eq!(tree.predict(&[5.05]), 1);
        assert_eq!(tree.predict(&[10.05]), 2);
    }

    #[test]
    fn gini_pure() {
        let labels = vec![0, 0, 0];
        let indices = vec![0, 1, 2];
        assert!((gini_impurity(&labels, &indices) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn gini_half() {
        let labels = vec![0, 1, 0, 1];
        let indices = vec![0, 1, 2, 3];
        // Gini of 50/50 binary = 1 - 0.25 - 0.25 = 0.5
        assert!((gini_impurity(&labels, &indices) - 0.5).abs() < 1e-12);
    }
}
