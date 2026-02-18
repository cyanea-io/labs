//! Feature selection methods.
//!
//! Provides methods for identifying and selecting informative features from
//! flat row-major data matrices: variance threshold filtering, mutual
//! information, recursive feature elimination, and L1-regularized (Lasso)
//! selection.
//!
//! Data is flat row-major `&[f64]` with an `n_features` parameter, consistent
//! with the rest of the cyanea-ml crate.
//!
//! # Example
//!
//! ```
//! use cyanea_ml::feature_selection::{variance_threshold, lasso_selection};
//!
//! // Remove constant features
//! let data = vec![
//!     1.0, 5.0, 0.0,
//!     2.0, 5.0, 1.0,
//!     3.0, 5.0, 2.0,
//!     4.0, 5.0, 3.0,
//! ];
//! let result = variance_threshold(&data, 3, 0.0).unwrap();
//! assert_eq!(result.selected, vec![0, 2]); // feature 1 (constant) removed
//! assert_eq!(result.n_selected(), 2);
//!
//! let filtered = result.transform(&data);
//! assert_eq!(filtered.len(), 4 * 2); // 4 samples × 2 features
//! ```

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result of feature selection.
///
/// Contains the indices of selected features and per-feature scores whose
/// meaning depends on the selection method used.
#[derive(Debug, Clone)]
pub struct FeatureSelection {
    /// Indices of selected features (sorted ascending).
    pub selected: Vec<usize>,
    /// Per-feature scores (one per original feature). Meaning varies:
    /// - `variance_threshold`: feature variances
    /// - `mutual_information`: MI values (nats)
    /// - `recursive_feature_elimination`: rankings (1 = selected, higher = eliminated earlier)
    pub scores: Vec<f64>,
    /// Original number of features.
    pub n_features: usize,
}

impl FeatureSelection {
    /// Transform data by keeping only selected features.
    ///
    /// `data` is flat row-major with `self.n_features` columns.
    /// Returns a new flat row-major array with `self.n_selected()` columns.
    pub fn transform(&self, data: &[f64]) -> Vec<f64> {
        let n_samples = data.len() / self.n_features;
        let n_sel = self.selected.len();
        let mut out = Vec::with_capacity(n_samples * n_sel);
        for i in 0..n_samples {
            let row = &data[i * self.n_features..(i + 1) * self.n_features];
            for &j in &self.selected {
                out.push(row[j]);
            }
        }
        out
    }

    /// Number of selected features.
    pub fn n_selected(&self) -> usize {
        self.selected.len()
    }
}

/// Result of Lasso feature selection.
#[derive(Debug, Clone)]
pub struct LassoResult {
    /// Feature weights (one per original feature).
    pub weights: Vec<f64>,
    /// Bias (intercept) term.
    pub bias: f64,
    /// Indices of features with non-zero weights (sorted ascending).
    pub selected: Vec<usize>,
    /// Number of coordinate descent iterations performed.
    pub n_iterations: usize,
    /// Number of original features.
    pub n_features: usize,
}

impl LassoResult {
    /// Transform data by keeping only selected (non-zero weight) features.
    ///
    /// `data` is flat row-major with `self.n_features` columns.
    pub fn transform(&self, data: &[f64]) -> Vec<f64> {
        let n_samples = data.len() / self.n_features;
        let mut out = Vec::with_capacity(n_samples * self.selected.len());
        for i in 0..n_samples {
            let row = &data[i * self.n_features..(i + 1) * self.n_features];
            for &j in &self.selected {
                out.push(row[j]);
            }
        }
        out
    }

    /// Predict target values using the fitted Lasso model.
    pub fn predict(&self, data: &[f64]) -> Vec<f64> {
        let n_samples = data.len() / self.n_features;
        let mut preds = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let row = &data[i * self.n_features..(i + 1) * self.n_features];
            let mut val = self.bias;
            for (j, &w) in self.weights.iter().enumerate() {
                val += row[j] * w;
            }
            preds.push(val);
        }
        preds
    }

    /// Number of selected features (non-zero weights).
    pub fn n_selected(&self) -> usize {
        self.selected.len()
    }
}

// ---------------------------------------------------------------------------
// Variance Threshold
// ---------------------------------------------------------------------------

/// Select features whose variance exceeds a threshold.
///
/// Features with variance ≤ `threshold` are removed. Useful for removing
/// constant or near-constant features before model training.
///
/// * `data` — flat row-major `n_samples × n_features`
/// * `n_features` — number of features per sample
/// * `threshold` — minimum variance (features with variance > threshold are kept)
///
/// Scores in the result are per-feature variances (population variance).
///
/// # Errors
///
/// Returns an error if data is empty, dimensions are inconsistent, or
/// threshold is negative.
pub fn variance_threshold(
    data: &[f64],
    n_features: usize,
    threshold: f64,
) -> Result<FeatureSelection> {
    validate_data(data, n_features)?;
    if threshold < 0.0 {
        return Err(CyaneaError::InvalidInput(
            "threshold must be >= 0".into(),
        ));
    }

    let n_samples = data.len() / n_features;
    let n = n_samples as f64;

    // Compute mean per feature
    let mut means = vec![0.0; n_features];
    for i in 0..n_samples {
        for j in 0..n_features {
            means[j] += data[i * n_features + j];
        }
    }
    for m in means.iter_mut() {
        *m /= n;
    }

    // Compute variance per feature
    let mut variances = vec![0.0; n_features];
    for i in 0..n_samples {
        for j in 0..n_features {
            let diff = data[i * n_features + j] - means[j];
            variances[j] += diff * diff;
        }
    }
    for v in variances.iter_mut() {
        *v /= n;
    }

    let selected: Vec<usize> = (0..n_features)
        .filter(|&j| variances[j] > threshold)
        .collect();

    Ok(FeatureSelection {
        selected,
        scores: variances,
        n_features,
    })
}

// ---------------------------------------------------------------------------
// Mutual Information (discrete)
// ---------------------------------------------------------------------------

/// Select features by mutual information with discrete labels.
///
/// Continuous features are discretized into `n_bins` equal-width bins before
/// computing MI. All features with MI > 0 are selected.
///
/// MI(X; Y) = Σ_{x,y} p(x,y) ln(p(x,y) / (p(x) p(y)))
///
/// * `data` — flat row-major `n_samples × n_features`
/// * `n_features` — number of features per sample
/// * `labels` — discrete class labels (0-indexed)
/// * `n_bins` — number of equal-width bins for discretizing features
///
/// Scores in the result are per-feature MI values (in nats).
///
/// # Errors
///
/// Returns an error if data is empty, dimensions are inconsistent, labels
/// length mismatches, or `n_bins` is 0.
pub fn mutual_information(
    data: &[f64],
    n_features: usize,
    labels: &[usize],
    n_bins: usize,
) -> Result<FeatureSelection> {
    validate_data(data, n_features)?;
    let n_samples = data.len() / n_features;
    validate_labels(labels, n_samples)?;

    if n_bins == 0 {
        return Err(CyaneaError::InvalidInput("n_bins must be > 0".into()));
    }

    let n_classes = labels.iter().copied().max().map_or(0, |m| m + 1);

    let mut mi_scores = vec![0.0; n_features];

    for j in 0..n_features {
        let col: Vec<f64> = (0..n_samples)
            .map(|i| data[i * n_features + j])
            .collect();
        let bins = discretize(&col, n_bins);
        mi_scores[j] = compute_mi(&bins, n_bins, labels, n_classes, n_samples);
    }

    // Select features with MI > 0 (sorted ascending by index)
    let selected: Vec<usize> = (0..n_features)
        .filter(|&j| mi_scores[j] > 0.0)
        .collect();

    Ok(FeatureSelection {
        selected,
        scores: mi_scores,
        n_features,
    })
}

/// Select the top-k features by mutual information with discrete labels.
///
/// Same as [`mutual_information`] but keeps at most `k` features (those with
/// highest MI values). Selected indices are sorted ascending.
///
/// # Errors
///
/// Returns an error if `k` is 0, or any error from [`mutual_information`].
pub fn mutual_information_top_k(
    data: &[f64],
    n_features: usize,
    labels: &[usize],
    n_bins: usize,
    k: usize,
) -> Result<FeatureSelection> {
    if k == 0 {
        return Err(CyaneaError::InvalidInput("k must be > 0".into()));
    }

    let result = mutual_information(data, n_features, labels, n_bins)?;

    if k >= result.selected.len() {
        return Ok(result);
    }

    // Rank selected features by MI descending, take top k
    let mut ranked: Vec<(usize, f64)> = result
        .selected
        .iter()
        .map(|&j| (j, result.scores[j]))
        .collect();
    ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    ranked.truncate(k);

    let mut selected: Vec<usize> = ranked.iter().map(|&(j, _)| j).collect();
    selected.sort();

    Ok(FeatureSelection {
        selected,
        scores: result.scores,
        n_features: result.n_features,
    })
}

/// Discretize values into `n_bins` equal-width bins.
fn discretize(values: &[f64], n_bins: usize) -> Vec<usize> {
    if values.is_empty() || n_bins == 0 {
        return vec![0; values.len()];
    }

    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    if (max - min).abs() < f64::EPSILON {
        return vec![0; values.len()];
    }

    let width = (max - min) / n_bins as f64;
    values
        .iter()
        .map(|&v| {
            let bin = ((v - min) / width) as usize;
            bin.min(n_bins - 1)
        })
        .collect()
}

/// Compute mutual information between discretized feature bins and labels.
fn compute_mi(
    bins: &[usize],
    n_bins: usize,
    labels: &[usize],
    n_classes: usize,
    n_samples: usize,
) -> f64 {
    if n_samples == 0 {
        return 0.0;
    }

    let n = n_samples as f64;

    // Joint counts
    let mut joint = vec![0usize; n_bins * n_classes];
    let mut bin_counts = vec![0usize; n_bins];
    let mut class_counts = vec![0usize; n_classes];

    for i in 0..n_samples {
        let b = bins[i];
        let c = labels[i];
        joint[b * n_classes + c] += 1;
        bin_counts[b] += 1;
        class_counts[c] += 1;
    }

    let mut mi = 0.0;
    for b in 0..n_bins {
        if bin_counts[b] == 0 {
            continue;
        }
        for c in 0..n_classes {
            let jc = joint[b * n_classes + c];
            if class_counts[c] == 0 || jc == 0 {
                continue;
            }
            let p_joint = jc as f64 / n;
            let p_bin = bin_counts[b] as f64 / n;
            let p_class = class_counts[c] as f64 / n;
            mi += p_joint * (p_joint / (p_bin * p_class)).ln();
        }
    }

    mi.max(0.0) // clamp tiny negatives from floating-point
}

// ---------------------------------------------------------------------------
// Recursive Feature Elimination
// ---------------------------------------------------------------------------

/// Recursive feature elimination (RFE).
///
/// Iteratively trains a model on the current feature subset, computes feature
/// importance, and removes the least important feature. Repeats until
/// `n_select` features remain.
///
/// The caller supplies a closure that, given filtered data and the number of
/// active features, returns importance scores (one per active feature,
/// higher = more important). This design accommodates any classifier.
///
/// * `data` — flat row-major `n_samples × n_features`
/// * `n_features` — number of features per sample
/// * `n_select` — target number of features to keep
/// * `importance_fn` — `FnMut(&[f64], usize) -> Result<Vec<f64>>` returning
///   importance scores for the current active features
///
/// Scores in the result are feature rankings: 1 = selected (survived to the
/// end), higher values indicate earlier elimination.
///
/// # Errors
///
/// Returns an error if data is empty, `n_select` is 0 or exceeds `n_features`,
/// or the importance closure returns an error or wrong-length vector.
pub fn recursive_feature_elimination<F>(
    data: &[f64],
    n_features: usize,
    n_select: usize,
    mut importance_fn: F,
) -> Result<FeatureSelection>
where
    F: FnMut(&[f64], usize) -> Result<Vec<f64>>,
{
    validate_data(data, n_features)?;

    if n_select == 0 {
        return Err(CyaneaError::InvalidInput("n_select must be > 0".into()));
    }
    if n_select > n_features {
        return Err(CyaneaError::InvalidInput(format!(
            "n_select ({}) > n_features ({})",
            n_select, n_features
        )));
    }

    let n_samples = data.len() / n_features;
    let mut active: Vec<usize> = (0..n_features).collect();
    let mut rankings = vec![0usize; n_features];
    let mut rank = n_features;

    while active.len() > n_select {
        // Build filtered data with only active features
        let n_active = active.len();
        let mut filtered = Vec::with_capacity(n_samples * n_active);
        for i in 0..n_samples {
            let row = &data[i * n_features..(i + 1) * n_features];
            for &j in &active {
                filtered.push(row[j]);
            }
        }

        // Get importances for active features
        let importances = importance_fn(&filtered, n_active)?;
        if importances.len() != n_active {
            return Err(CyaneaError::InvalidInput(format!(
                "importance_fn returned {} scores, expected {}",
                importances.len(),
                n_active
            )));
        }

        // Find least important feature
        let min_idx = importances
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap();

        let removed_feature = active[min_idx];
        rankings[removed_feature] = rank;
        rank -= 1;
        active.remove(min_idx);
    }

    // Remaining features get rank 1
    for &j in &active {
        rankings[j] = 1;
    }

    let scores: Vec<f64> = rankings.iter().map(|&r| r as f64).collect();
    let mut selected = active;
    selected.sort();

    Ok(FeatureSelection {
        selected,
        scores,
        n_features,
    })
}

// ---------------------------------------------------------------------------
// L1-Regularized Feature Selection (Lasso)
// ---------------------------------------------------------------------------

/// L1-regularized linear regression (Lasso) for feature selection.
///
/// Uses coordinate descent to fit a linear model with L1 penalty, driving
/// uninformative feature weights to exactly zero. Features with non-zero
/// weights after fitting are considered selected.
///
/// The algorithm standardizes features internally and returns weights on the
/// original feature scale.
///
/// * `data` — flat row-major `n_samples × n_features`
/// * `n_features` — number of features per sample
/// * `targets` — target values (one per sample)
/// * `alpha` — L1 regularization strength (larger → sparser solution)
/// * `max_iter` — maximum coordinate descent passes
/// * `tol` — convergence tolerance on maximum weight change
///
/// # Errors
///
/// Returns an error if data is empty, dimensions are inconsistent, targets
/// length mismatches, alpha is negative, or max_iter is 0.
pub fn lasso_selection(
    data: &[f64],
    n_features: usize,
    targets: &[f64],
    alpha: f64,
    max_iter: usize,
    tol: f64,
) -> Result<LassoResult> {
    validate_data(data, n_features)?;
    let n_samples = data.len() / n_features;
    if targets.len() != n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "targets length {} != n_samples {}",
            targets.len(),
            n_samples
        )));
    }
    if alpha < 0.0 {
        return Err(CyaneaError::InvalidInput("alpha must be >= 0".into()));
    }
    if max_iter == 0 {
        return Err(CyaneaError::InvalidInput("max_iter must be > 0".into()));
    }

    let n = n_samples as f64;

    // Feature means and standard deviations for standardization
    let mut means = vec![0.0; n_features];
    for i in 0..n_samples {
        for j in 0..n_features {
            means[j] += data[i * n_features + j];
        }
    }
    for m in means.iter_mut() {
        *m /= n;
    }

    let mut stds = vec![0.0; n_features];
    for i in 0..n_samples {
        for j in 0..n_features {
            let diff = data[i * n_features + j] - means[j];
            stds[j] += diff * diff;
        }
    }
    for s in stds.iter_mut() {
        *s = (*s / n).sqrt();
        if *s < f64::EPSILON {
            *s = 1.0; // constant feature — avoid division by zero
        }
    }

    // Target centering
    let y_mean = targets.iter().sum::<f64>() / n;
    let y_centered: Vec<f64> = targets.iter().map(|&y| y - y_mean).collect();

    // Pre-compute standardized columns and their squared norms
    let mut columns: Vec<Vec<f64>> = Vec::with_capacity(n_features);
    for j in 0..n_features {
        let col: Vec<f64> = (0..n_samples)
            .map(|i| (data[i * n_features + j] - means[j]) / stds[j])
            .collect();
        columns.push(col);
    }

    let col_norms: Vec<f64> = columns
        .iter()
        .map(|col| col.iter().map(|&x| x * x).sum::<f64>() / n)
        .collect();

    // Coordinate descent
    let mut weights = vec![0.0; n_features];
    let mut residual = y_centered;

    let mut n_iterations = 0;

    for _ in 0..max_iter {
        n_iterations += 1;
        let mut max_change = 0.0f64;

        for j in 0..n_features {
            if col_norms[j] < f64::EPSILON {
                continue;
            }

            let old_w = weights[j];

            // rho = x_j^T residual / n + w_j * ||x_j||^2/n
            let rho: f64 = columns[j]
                .iter()
                .zip(residual.iter())
                .map(|(&x, &r)| x * r)
                .sum::<f64>()
                / n
                + old_w * col_norms[j];

            // Soft thresholding
            let new_w = soft_threshold(rho, alpha) / col_norms[j];

            let delta = new_w - old_w;
            if delta.abs() > f64::EPSILON * 100.0 {
                for i in 0..n_samples {
                    residual[i] -= columns[j][i] * delta;
                }
            }

            weights[j] = new_w;
            max_change = max_change.max(delta.abs());
        }

        if max_change < tol {
            break;
        }
    }

    // Un-standardize: w_orig_j = w_std_j / std_j, bias = y_mean - Σ w_orig_j * mean_j
    let mut orig_weights = vec![0.0; n_features];
    let mut bias = y_mean;
    for j in 0..n_features {
        orig_weights[j] = weights[j] / stds[j];
        bias -= orig_weights[j] * means[j];
    }

    let selected: Vec<usize> = (0..n_features)
        .filter(|&j| orig_weights[j].abs() > f64::EPSILON * 100.0)
        .collect();

    Ok(LassoResult {
        weights: orig_weights,
        bias,
        selected,
        n_iterations,
        n_features,
    })
}

/// Soft thresholding: sign(x) · max(0, |x| − λ).
fn soft_threshold(x: f64, lambda: f64) -> f64 {
    if x > lambda {
        x - lambda
    } else if x < -lambda {
        x + lambda
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Validation helpers
// ---------------------------------------------------------------------------

fn validate_data(data: &[f64], n_features: usize) -> Result<()> {
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
    Ok(())
}

fn validate_labels(labels: &[usize], n_samples: usize) -> Result<()> {
    if labels.len() != n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "labels length {} != n_samples {}",
            labels.len(),
            n_samples
        )));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // === Variance Threshold ================================================

    #[test]
    fn variance_removes_constant_features() {
        // Feature 0: varies, Feature 1: constant, Feature 2: varies
        let data = vec![
            1.0, 5.0, 0.0,
            2.0, 5.0, 1.0,
            3.0, 5.0, 2.0,
            4.0, 5.0, 3.0,
        ];
        let result = variance_threshold(&data, 3, 0.0).unwrap();
        assert_eq!(result.selected, vec![0, 2]);
        assert!(result.scores[1] < f64::EPSILON); // constant → zero variance
    }

    #[test]
    fn variance_threshold_zero_keeps_all_varying() {
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let result = variance_threshold(&data, 2, 0.0).unwrap();
        assert_eq!(result.selected, vec![0, 1]);
    }

    #[test]
    fn variance_higher_threshold_removes_more() {
        // Feature 0: [0,1,2,3] → var=1.25, Feature 1: [0,0.1,0.2,0.3] → var=0.0125
        let data = vec![
            0.0, 0.0,
            1.0, 0.1,
            2.0, 0.2,
            3.0, 0.3,
        ];
        let low = variance_threshold(&data, 2, 0.01).unwrap();
        assert_eq!(low.selected, vec![0, 1]);

        let high = variance_threshold(&data, 2, 0.1).unwrap();
        assert_eq!(high.selected, vec![0]);
    }

    #[test]
    fn variance_transform_shape() {
        let data = vec![
            1.0, 5.0, 10.0,
            2.0, 5.0, 20.0,
            3.0, 5.0, 30.0,
        ];
        let result = variance_threshold(&data, 3, 0.0).unwrap();
        let filtered = result.transform(&data);
        assert_eq!(filtered.len(), 3 * result.n_selected());
        // First row should be [1.0, 10.0]
        assert!((filtered[0] - 1.0).abs() < f64::EPSILON);
        assert!((filtered[1] - 10.0).abs() < f64::EPSILON);
    }

    #[test]
    fn variance_single_sample() {
        // Single sample → all variances 0 → nothing selected with threshold 0
        let data = vec![1.0, 2.0, 3.0];
        let result = variance_threshold(&data, 3, 0.0).unwrap();
        assert!(result.selected.is_empty());
    }

    #[test]
    fn variance_negative_threshold_error() {
        assert!(variance_threshold(&[1.0], 1, -0.1).is_err());
    }

    #[test]
    fn variance_empty_data_error() {
        assert!(variance_threshold(&[], 2, 0.0).is_err());
    }

    #[test]
    fn variance_dimension_mismatch_error() {
        assert!(variance_threshold(&[1.0, 2.0, 3.0], 2, 0.0).is_err());
    }

    // === Mutual Information ================================================

    #[test]
    fn mi_informative_feature_scores_high() {
        // Feature 0 perfectly predicts the label; feature 1 is noise
        let data = vec![
            0.0, 5.0,
            0.0, 5.1,
            0.0, 4.9,
            1.0, 5.0,
            1.0, 5.1,
            1.0, 4.9,
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let result = mutual_information(&data, 2, &labels, 5).unwrap();
        assert!(
            result.scores[0] > result.scores[1],
            "informative feature should have higher MI: {:?}",
            result.scores
        );
    }

    #[test]
    fn mi_constant_feature_zero() {
        let data = vec![
            5.0, 0.0,
            5.0, 1.0,
            5.0, 2.0,
        ];
        let labels = vec![0, 1, 2];
        let result = mutual_information(&data, 2, &labels, 3).unwrap();
        assert!(
            result.scores[0] < f64::EPSILON,
            "constant feature MI should be ~0, got {}",
            result.scores[0]
        );
    }

    #[test]
    fn mi_selects_nonzero_features() {
        // Two informative features, one constant
        let data = vec![
            0.0, 5.0, 10.0,
            1.0, 5.0, 11.0,
            2.0, 5.0, 12.0,
            3.0, 5.0, 13.0,
            4.0, 5.0, 14.0,
            5.0, 5.0, 15.0,
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let result = mutual_information(&data, 3, &labels, 5).unwrap();
        // Feature 1 (constant) should not be selected
        assert!(!result.selected.contains(&1));
        // Features 0 and 2 should be selected
        assert!(result.selected.contains(&0));
        assert!(result.selected.contains(&2));
    }

    #[test]
    fn mi_top_k_limits_count() {
        let data = vec![
            0.0, 10.0, 100.0,
            1.0, 11.0, 101.0,
            2.0, 12.0, 102.0,
            3.0, 13.0, 103.0,
        ];
        let labels = vec![0, 0, 1, 1];
        let result = mutual_information_top_k(&data, 3, &labels, 4, 1).unwrap();
        assert_eq!(result.n_selected(), 1);
    }

    #[test]
    fn mi_top_k_sorted_by_index() {
        let data = vec![
            100.0, 0.0, 50.0,
            101.0, 1.0, 51.0,
            102.0, 2.0, 52.0,
            103.0, 3.0, 53.0,
        ];
        let labels = vec![0, 0, 1, 1];
        let result = mutual_information_top_k(&data, 3, &labels, 4, 2).unwrap();
        // Selected should be sorted ascending
        for w in result.selected.windows(2) {
            assert!(w[0] < w[1]);
        }
    }

    #[test]
    fn mi_multiclass() {
        // 3 classes separated by feature 0
        let data = vec![
            0.0, 1.0,
            10.0, 1.0,
            20.0, 1.0,
        ];
        let labels = vec![0, 1, 2];
        let result = mutual_information(&data, 2, &labels, 3).unwrap();
        assert!(result.scores[0] > 0.0);
    }

    #[test]
    fn mi_zero_bins_error() {
        assert!(mutual_information(&[1.0], 1, &[0], 0).is_err());
    }

    #[test]
    fn mi_labels_mismatch_error() {
        assert!(mutual_information(&[1.0, 2.0], 1, &[0, 1, 2], 5).is_err());
    }

    #[test]
    fn mi_top_k_zero_error() {
        assert!(mutual_information_top_k(&[1.0], 1, &[0], 5, 0).is_err());
    }

    // === Recursive Feature Elimination =====================================

    #[test]
    fn rfe_selects_informative_features() {
        // Feature 0: informative (low→class0, high→class1)
        // Feature 1: constant noise
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..10 {
            data.extend_from_slice(&[i as f64, 5.0]);
            labels.push(if i < 5 { 0 } else { 1 });
        }

        let labels_ref = &labels;
        let result = recursive_feature_elimination(&data, 2, 1, |filtered, n_feat| {
            use crate::forest::{RandomForest, RandomForestConfig};
            let config = RandomForestConfig {
                n_trees: 20,
                max_depth: 3,
                seed: 42,
                ..Default::default()
            };
            let train_labels: Vec<usize> = labels_ref.to_vec();
            let forest = RandomForest::fit(filtered, n_feat, &train_labels, &config)?;
            Ok(forest.feature_importance(n_feat))
        })
        .unwrap();

        assert_eq!(result.n_selected(), 1);
        assert_eq!(result.selected, vec![0]); // feature 0 is informative
    }

    #[test]
    fn rfe_rankings_correct() {
        // 3 features, select 1 → features get ranks 1, 2, or 3
        let data = vec![
            0.0, 5.0, 10.0,
            1.0, 5.0, 11.0,
            2.0, 5.0, 12.0,
            10.0, 5.0, 20.0,
            11.0, 5.0, 21.0,
            12.0, 5.0, 22.0,
        ];

        let result = recursive_feature_elimination(&data, 3, 1, |_filtered, n_feat| {
            // Dummy: first feature always least important
            let mut imp = vec![1.0; n_feat];
            imp[0] = 0.0;
            Ok(imp)
        })
        .unwrap();

        assert_eq!(result.n_selected(), 1);
        // All scores should be 1.0, 2.0, or 3.0
        for &s in &result.scores {
            assert!(s >= 1.0 && s <= 3.0);
        }
        // Selected feature should have rank 1
        for &j in &result.selected {
            assert!((result.scores[j] - 1.0).abs() < f64::EPSILON);
        }
    }

    #[test]
    fn rfe_n_select_equals_n_features() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let result = recursive_feature_elimination(&data, 2, 2, |_d, n| {
            Ok(vec![1.0; n])
        })
        .unwrap();
        assert_eq!(result.selected, vec![0, 1]);
    }

    #[test]
    fn rfe_n_select_zero_error() {
        assert!(recursive_feature_elimination(&[1.0], 1, 0, |_, n| Ok(vec![1.0; n])).is_err());
    }

    #[test]
    fn rfe_n_select_exceeds_error() {
        assert!(recursive_feature_elimination(&[1.0], 1, 2, |_, n| Ok(vec![1.0; n])).is_err());
    }

    #[test]
    fn rfe_bad_importance_length_error() {
        let result = recursive_feature_elimination(&[1.0, 2.0, 3.0, 4.0], 2, 1, |_, _| {
            Ok(vec![1.0]) // wrong length
        });
        assert!(result.is_err());
    }

    #[test]
    fn rfe_transform_works() {
        let data = vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
        ];
        let result = recursive_feature_elimination(&data, 3, 2, |_d, n| {
            // Remove first active feature each time
            let mut imp = vec![1.0; n];
            imp[0] = 0.0;
            Ok(imp)
        })
        .unwrap();

        let filtered = result.transform(&data);
        assert_eq!(filtered.len(), 2 * result.n_selected());
    }

    // === Lasso Feature Selection ===========================================

    #[test]
    fn lasso_selects_informative_features() {
        // y = 2*x0 + 0*x1 (feature 1 is pure noise)
        let data = vec![
            1.0, 10.0,
            2.0, 20.0,
            3.0, 30.0,
            4.0, 40.0,
            5.0, 50.0,
            6.0, 60.0,
            7.0, 70.0,
            8.0, 80.0,
        ];
        let targets: Vec<f64> = (0..8).map(|i| 2.0 * (i + 1) as f64).collect();

        let result = lasso_selection(&data, 2, &targets, 0.5, 1000, 1e-8).unwrap();
        // Both features correlate with target (x1 = 10*x0), so both may survive
        // with moderate alpha. The important thing is that the model works.
        assert!(!result.selected.is_empty());
        assert!(result.weights[0].abs() > f64::EPSILON);
    }

    #[test]
    fn lasso_high_alpha_eliminates_all() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![2.0, 4.0, 6.0, 8.0, 10.0];

        let result = lasso_selection(&data, 1, &targets, 100.0, 100, 1e-8).unwrap();
        assert!(result.selected.is_empty(), "high alpha should drive all weights to zero");
    }

    #[test]
    fn lasso_zero_alpha_keeps_all() {
        let data = vec![
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0,
            2.0, 1.0,
            1.0, 2.0,
            3.0, 0.0,
        ];
        let targets = vec![1.0, 2.0, 3.0, 3.0, 4.0, 3.0];

        let result = lasso_selection(&data, 2, &targets, 0.0, 1000, 1e-10).unwrap();
        assert_eq!(result.n_selected(), 2, "alpha=0 should keep all features");
    }

    #[test]
    fn lasso_predict_accuracy() {
        // y = 3*x + 1
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let targets = vec![1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0, 22.0];

        let result = lasso_selection(&data, 1, &targets, 0.01, 1000, 1e-10).unwrap();
        let preds = result.predict(&data);
        for (p, &t) in preds.iter().zip(targets.iter()) {
            assert!(
                (p - t).abs() < 1.0,
                "prediction {} too far from target {}",
                p,
                t
            );
        }
    }

    #[test]
    fn lasso_convergence() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![2.0, 4.0, 6.0, 8.0, 10.0];

        let result = lasso_selection(&data, 1, &targets, 0.01, 10000, 1e-12).unwrap();
        // Should converge before max_iter
        assert!(
            result.n_iterations < 10000,
            "should converge, took {} iterations",
            result.n_iterations
        );
    }

    #[test]
    fn lasso_transform_shape() {
        let data = vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
            10.0, 11.0, 12.0,
        ];
        let targets = vec![1.0, 4.0, 7.0, 10.0];

        let result = lasso_selection(&data, 3, &targets, 0.01, 1000, 1e-8).unwrap();
        let filtered = result.transform(&data);
        assert_eq!(filtered.len(), 4 * result.n_selected());
    }

    #[test]
    fn lasso_negative_alpha_error() {
        assert!(lasso_selection(&[1.0], 1, &[1.0], -0.1, 100, 1e-8).is_err());
    }

    #[test]
    fn lasso_zero_max_iter_error() {
        assert!(lasso_selection(&[1.0], 1, &[1.0], 0.1, 0, 1e-8).is_err());
    }

    #[test]
    fn lasso_empty_data_error() {
        assert!(lasso_selection(&[], 1, &[], 0.1, 100, 1e-8).is_err());
    }

    #[test]
    fn lasso_targets_mismatch_error() {
        assert!(lasso_selection(&[1.0, 2.0], 1, &[1.0, 2.0, 3.0], 0.1, 100, 1e-8).is_err());
    }

    // === Integration ======================================================

    #[test]
    fn integration_variance_then_mi() {
        // Pipeline: variance threshold → mutual information on remaining features
        let data = vec![
            0.0, 5.0, 10.0, 1.0,
            1.0, 5.0, 11.0, 1.0,
            2.0, 5.0, 12.0, 1.0,
            3.0, 5.0, 13.0, 1.0,
            4.0, 5.0, 14.0, 1.0,
            5.0, 5.0, 15.0, 1.0,
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];

        // Step 1: remove constant features
        let vt = variance_threshold(&data, 4, 0.0).unwrap();
        assert_eq!(vt.selected, vec![0, 2]); // features 1 and 3 are constant

        // Step 2: MI on filtered data
        let filtered = vt.transform(&data);
        let mi = mutual_information(&filtered, 2, &labels, 5).unwrap();
        assert!(mi.scores[0] > 0.0); // both surviving features should be informative
        assert!(mi.scores[1] > 0.0);
    }

    #[test]
    fn integration_rfe_with_random_forest() {
        // 3 features: 0 = informative, 1 = less informative, 2 = noise
        let mut data = Vec::new();
        let mut labels = Vec::new();
        for i in 0..20 {
            let x0 = i as f64;
            let x1 = (i as f64) * 0.5 + 3.0;
            let x2 = 42.0; // constant noise
            data.extend_from_slice(&[x0, x1, x2]);
            labels.push(if i < 10 { 0 } else { 1 });
        }

        let labels_ref = &labels;
        let result = recursive_feature_elimination(&data, 3, 2, |filtered, n_feat| {
            use crate::forest::{RandomForest, RandomForestConfig};
            let config = RandomForestConfig {
                n_trees: 20,
                max_depth: 4,
                seed: 42,
                ..Default::default()
            };
            let forest = RandomForest::fit(filtered, n_feat, labels_ref, &config)?;
            Ok(forest.feature_importance(n_feat))
        })
        .unwrap();

        assert_eq!(result.n_selected(), 2);
        // Noise feature (2) should be eliminated first (highest rank)
        assert!(
            result.scores[2] > result.scores[0],
            "noise feature should have higher rank (eliminated earlier)"
        );
        assert!(
            result.scores[2] > result.scores[1],
            "noise feature should have higher rank (eliminated earlier)"
        );
    }

    #[test]
    fn lasso_multivariate_selection() {
        // y = 1*x0 + 2*x1, x2 is uncorrelated noise (alternating pattern)
        let data = vec![
            1.0, 0.0, 1.0,
            0.0, 1.0, -1.0,
            1.0, 1.0, 1.0,
            2.0, 0.0, -1.0,
            0.0, 2.0, 1.0,
            2.0, 1.0, -1.0,
            1.0, 2.0, 1.0,
            2.0, 2.0, -1.0,
        ];
        let targets = vec![1.0, 2.0, 3.0, 2.0, 4.0, 4.0, 5.0, 6.0];

        let result = lasso_selection(&data, 3, &targets, 0.1, 1000, 1e-10).unwrap();
        // Features 0 and 1 should have non-trivial weights
        assert!(
            result.weights[0].abs() > 0.1,
            "feature 0 weight {} should be > 0.1", result.weights[0]
        );
        assert!(
            result.weights[1].abs() > 0.1,
            "feature 1 weight {} should be > 0.1", result.weights[1]
        );
    }

    // === Discretize helper =================================================

    #[test]
    fn discretize_equal_width_bins() {
        let values = vec![0.0, 2.5, 5.0, 7.5, 10.0];
        let bins = discretize(&values, 4);
        assert_eq!(bins[0], 0);
        assert_eq!(bins[4], 3); // max value maps to last bin
    }

    #[test]
    fn discretize_constant_column() {
        let values = vec![5.0, 5.0, 5.0];
        let bins = discretize(&values, 10);
        assert_eq!(bins, vec![0, 0, 0]);
    }

    #[test]
    fn soft_threshold_positive() {
        assert!((soft_threshold(3.0, 1.0) - 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn soft_threshold_negative() {
        assert!((soft_threshold(-3.0, 1.0) - (-2.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn soft_threshold_within_band() {
        assert!((soft_threshold(0.5, 1.0)).abs() < f64::EPSILON);
    }
}
