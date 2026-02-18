//! Classification evaluation metrics.
//!
//! Provides confusion matrix computation, per-class precision / recall / F1,
//! macro/weighted averaging, Matthews correlation coefficient, ROC curves,
//! precision-recall curves, and AUC (trapezoidal).

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Confusion Matrix
// ---------------------------------------------------------------------------

/// Row-major confusion matrix for multi-class classification.
///
/// Entry `(i, j)` counts samples whose **actual** class is `i` and
/// **predicted** class is `j`.
#[derive(Debug, Clone)]
pub struct ConfusionMatrix {
    /// Row-major storage: `matrix[actual * n_classes + predicted]`.
    pub matrix: Vec<usize>,
    /// Number of classes.
    pub n_classes: usize,
}

impl ConfusionMatrix {
    /// Build a confusion matrix from actual and predicted label vectors.
    ///
    /// `n_classes` is inferred from the maximum label + 1 when `None`.
    ///
    /// # Errors
    ///
    /// Returns an error if the slices are empty or have different lengths.
    pub fn from_labels(
        actual: &[usize],
        predicted: &[usize],
        n_classes: Option<usize>,
    ) -> Result<Self> {
        if actual.is_empty() {
            return Err(CyaneaError::InvalidInput("empty label vectors".into()));
        }
        if actual.len() != predicted.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "actual length {} != predicted length {}",
                actual.len(),
                predicted.len()
            )));
        }

        let nc = n_classes.unwrap_or_else(|| {
            let max_a = actual.iter().copied().max().unwrap_or(0);
            let max_p = predicted.iter().copied().max().unwrap_or(0);
            max_a.max(max_p) + 1
        });

        let mut matrix = vec![0usize; nc * nc];
        for (&a, &p) in actual.iter().zip(predicted.iter()) {
            if a < nc && p < nc {
                matrix[a * nc + p] += 1;
            }
        }

        Ok(Self {
            matrix,
            n_classes: nc,
        })
    }

    /// Get the count for a specific (actual, predicted) pair.
    #[inline]
    pub fn get(&self, actual: usize, predicted: usize) -> usize {
        self.matrix[actual * self.n_classes + predicted]
    }

    /// Total number of samples.
    pub fn total(&self) -> usize {
        self.matrix.iter().sum()
    }

    /// True positives for a given class.
    pub fn true_positives(&self, class: usize) -> usize {
        self.get(class, class)
    }

    /// False positives for a given class (predicted as `class` but actually
    /// something else).
    pub fn false_positives(&self, class: usize) -> usize {
        let mut fp = 0;
        for i in 0..self.n_classes {
            if i != class {
                fp += self.get(i, class);
            }
        }
        fp
    }

    /// True negatives for a given class.
    pub fn true_negatives(&self, class: usize) -> usize {
        let mut tn = 0;
        for i in 0..self.n_classes {
            for j in 0..self.n_classes {
                if i != class && j != class {
                    tn += self.get(i, j);
                }
            }
        }
        tn
    }

    /// False negatives for a given class (actually `class` but predicted as
    /// something else).
    pub fn false_negatives(&self, class: usize) -> usize {
        let mut fne = 0;
        for j in 0..self.n_classes {
            if j != class {
                fne += self.get(class, j);
            }
        }
        fne
    }

    /// Overall accuracy (correct predictions / total).
    pub fn accuracy(&self) -> f64 {
        let total = self.total();
        if total == 0 {
            return 0.0;
        }
        let correct: usize = (0..self.n_classes).map(|c| self.get(c, c)).sum();
        correct as f64 / total as f64
    }

    /// Precision for a given class: `TP / (TP + FP)`.
    ///
    /// Returns 0.0 if `TP + FP == 0`.
    pub fn precision(&self, class: usize) -> f64 {
        let tp = self.true_positives(class);
        let fp = self.false_positives(class);
        let denom = tp + fp;
        if denom == 0 {
            0.0
        } else {
            tp as f64 / denom as f64
        }
    }

    /// Recall (sensitivity) for a given class: `TP / (TP + FN)`.
    ///
    /// Returns 0.0 if `TP + FN == 0`.
    pub fn recall(&self, class: usize) -> f64 {
        let tp = self.true_positives(class);
        let fne = self.false_negatives(class);
        let denom = tp + fne;
        if denom == 0 {
            0.0
        } else {
            tp as f64 / denom as f64
        }
    }

    /// F1 score for a given class (harmonic mean of precision and recall).
    ///
    /// Returns 0.0 if both precision and recall are 0.
    pub fn f1(&self, class: usize) -> f64 {
        let p = self.precision(class);
        let r = self.recall(class);
        if p + r == 0.0 {
            0.0
        } else {
            2.0 * p * r / (p + r)
        }
    }

    /// Specificity for a given class: `TN / (TN + FP)`.
    ///
    /// Returns 0.0 if `TN + FP == 0`.
    pub fn specificity(&self, class: usize) -> f64 {
        let tn = self.true_negatives(class);
        let fp = self.false_positives(class);
        let denom = tn + fp;
        if denom == 0 {
            0.0
        } else {
            tn as f64 / denom as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Standalone scalar metrics
// ---------------------------------------------------------------------------

/// Overall accuracy: fraction of correct predictions.
///
/// # Errors
///
/// Returns an error if the slices are empty or have different lengths.
pub fn accuracy(actual: &[usize], predicted: &[usize]) -> Result<f64> {
    let cm = ConfusionMatrix::from_labels(actual, predicted, None)?;
    Ok(cm.accuracy())
}

/// F1 score for a specific class.
///
/// # Errors
///
/// Returns an error if the slices are empty or have different lengths.
pub fn f1_score(actual: &[usize], predicted: &[usize], class: usize) -> Result<f64> {
    let cm = ConfusionMatrix::from_labels(actual, predicted, None)?;
    Ok(cm.f1(class))
}

/// Macro-averaged F1 score (unweighted mean of per-class F1).
///
/// # Errors
///
/// Returns an error if the slices are empty or have different lengths.
pub fn f1_macro(actual: &[usize], predicted: &[usize]) -> Result<f64> {
    let cm = ConfusionMatrix::from_labels(actual, predicted, None)?;
    let sum: f64 = (0..cm.n_classes).map(|c| cm.f1(c)).sum();
    Ok(sum / cm.n_classes as f64)
}

/// Weighted-averaged F1 score (weighted by class support).
///
/// # Errors
///
/// Returns an error if the slices are empty or have different lengths.
pub fn f1_weighted(actual: &[usize], predicted: &[usize]) -> Result<f64> {
    let cm = ConfusionMatrix::from_labels(actual, predicted, None)?;
    let total = cm.total() as f64;
    if total == 0.0 {
        return Ok(0.0);
    }
    let mut weighted_sum = 0.0;
    for c in 0..cm.n_classes {
        let support = (cm.true_positives(c) + cm.false_negatives(c)) as f64;
        weighted_sum += cm.f1(c) * support;
    }
    Ok(weighted_sum / total)
}

/// Matthews correlation coefficient (MCC) for binary or multiclass.
///
/// For binary (2-class) case:
/// `MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))`
///
/// For multiclass, uses the generalized formula based on the confusion matrix:
/// `MCC = (c*s - sum(p_k * t_k)) / sqrt((s^2 - sum(p_k^2)) * (s^2 - sum(t_k^2)))`
/// where c = correct predictions, s = total samples, p_k = column sum for k,
/// t_k = row sum for k.
///
/// Returns 0.0 if the denominator is zero.
///
/// # Errors
///
/// Returns an error if the slices are empty or have different lengths.
pub fn matthews_corrcoef(actual: &[usize], predicted: &[usize]) -> Result<f64> {
    let cm = ConfusionMatrix::from_labels(actual, predicted, None)?;
    let nc = cm.n_classes;

    // Generalized MCC via the multiclass formula (reduces to binary MCC for nc=2)
    let s = cm.total() as f64;
    let c: f64 = (0..nc).map(|k| cm.get(k, k) as f64).sum();

    // p_k = sum of column k (total predicted as k)
    // t_k = sum of row k (total actual k)
    let mut sum_pk_sq = 0.0;
    let mut sum_tk_sq = 0.0;
    let mut sum_pk_tk = 0.0;
    for k in 0..nc {
        let pk: f64 = (0..nc).map(|i| cm.get(i, k) as f64).sum();
        let tk: f64 = (0..nc).map(|j| cm.get(k, j) as f64).sum();
        sum_pk_sq += pk * pk;
        sum_tk_sq += tk * tk;
        sum_pk_tk += pk * tk;
    }

    let numer = c * s - sum_pk_tk;
    let denom_sq = (s * s - sum_pk_sq) * (s * s - sum_tk_sq);
    if denom_sq <= 0.0 {
        return Ok(0.0);
    }
    Ok(numer / denom_sq.sqrt())
}

// ---------------------------------------------------------------------------
// ROC Curve
// ---------------------------------------------------------------------------

/// A single point on the ROC curve.
#[derive(Debug, Clone)]
pub struct RocPoint {
    /// Score threshold at which this point is computed.
    pub threshold: f64,
    /// False positive rate: FP / (FP + TN).
    pub fpr: f64,
    /// True positive rate (recall): TP / (TP + FN).
    pub tpr: f64,
}

/// ROC curve with AUC.
#[derive(Debug, Clone)]
pub struct RocCurve {
    /// Points on the curve, from (0, 0) to (1, 1).
    pub points: Vec<RocPoint>,
    /// Area under the ROC curve (trapezoidal rule).
    pub auc: f64,
}

/// Compute the ROC curve from predicted scores and binary labels.
///
/// Sorts by descending score and walks thresholds to compute (FPR, TPR)
/// at each distinct score. Includes endpoints (0,0) and (1,1).
///
/// # Errors
///
/// Returns an error if the slices are empty, have different lengths, or
/// contain no positive / no negative samples.
pub fn roc_curve(scores: &[f64], labels: &[bool]) -> Result<RocCurve> {
    if scores.is_empty() {
        return Err(CyaneaError::InvalidInput("empty input".into()));
    }
    if scores.len() != labels.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "scores length {} != labels length {}",
            scores.len(),
            labels.len()
        )));
    }

    let total_pos = labels.iter().filter(|&&l| l).count();
    let total_neg = labels.len() - total_pos;
    if total_pos == 0 {
        return Err(CyaneaError::InvalidInput("no positive samples".into()));
    }
    if total_neg == 0 {
        return Err(CyaneaError::InvalidInput("no negative samples".into()));
    }

    // Sort by descending score; ties: negatives first (pessimistic)
    let mut indices: Vec<usize> = (0..scores.len()).collect();
    indices.sort_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| labels[a].cmp(&labels[b])) // false < true, so negatives first
    });

    let p = total_pos as f64;
    let n = total_neg as f64;

    let mut points = Vec::new();
    // Start at origin (threshold = +inf, nothing predicted positive)
    points.push(RocPoint {
        threshold: f64::INFINITY,
        fpr: 0.0,
        tpr: 0.0,
    });

    let mut tp = 0usize;
    let mut fp = 0usize;

    let mut i = 0;
    while i < indices.len() {
        // Accumulate all samples with the same score
        let current_score = scores[indices[i]];
        while i < indices.len() && scores[indices[i]] == current_score {
            if labels[indices[i]] {
                tp += 1;
            } else {
                fp += 1;
            }
            i += 1;
        }

        points.push(RocPoint {
            threshold: current_score,
            fpr: fp as f64 / n,
            tpr: tp as f64 / p,
        });
    }

    // Compute AUC via trapezoidal rule
    let auc = trapezoidal_auc(
        &points.iter().map(|p| p.fpr).collect::<Vec<_>>(),
        &points.iter().map(|p| p.tpr).collect::<Vec<_>>(),
    );

    Ok(RocCurve { points, auc })
}

/// Compute only the AUC of the ROC curve.
///
/// Shorthand for `roc_curve(scores, labels)?.auc`.
pub fn roc_auc(scores: &[f64], labels: &[bool]) -> Result<f64> {
    Ok(roc_curve(scores, labels)?.auc)
}

// ---------------------------------------------------------------------------
// Precision-Recall Curve
// ---------------------------------------------------------------------------

/// A single point on the precision-recall curve.
#[derive(Debug, Clone)]
pub struct PrPoint {
    /// Score threshold at which this point is computed.
    pub threshold: f64,
    /// Precision: TP / (TP + FP).
    pub precision: f64,
    /// Recall: TP / (TP + FN).
    pub recall: f64,
}

/// Precision-recall curve with AUC.
#[derive(Debug, Clone)]
pub struct PrCurve {
    /// Points on the curve.
    pub points: Vec<PrPoint>,
    /// Area under the PR curve (trapezoidal rule).
    pub auc: f64,
}

/// Compute the precision-recall curve from predicted scores and binary labels.
///
/// # Errors
///
/// Returns an error if the slices are empty, have different lengths, or
/// contain no positive samples.
pub fn pr_curve(scores: &[f64], labels: &[bool]) -> Result<PrCurve> {
    if scores.is_empty() {
        return Err(CyaneaError::InvalidInput("empty input".into()));
    }
    if scores.len() != labels.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "scores length {} != labels length {}",
            scores.len(),
            labels.len()
        )));
    }

    let total_pos = labels.iter().filter(|&&l| l).count();
    if total_pos == 0 {
        return Err(CyaneaError::InvalidInput("no positive samples".into()));
    }

    // Sort by descending score; ties: negatives first
    let mut indices: Vec<usize> = (0..scores.len()).collect();
    indices.sort_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| labels[a].cmp(&labels[b]))
    });

    let p = total_pos as f64;

    let mut points = Vec::new();
    // Start at (recall=0, precision=1)
    points.push(PrPoint {
        threshold: f64::INFINITY,
        precision: 1.0,
        recall: 0.0,
    });

    let mut tp = 0usize;
    let mut fp = 0usize;

    let mut i = 0;
    while i < indices.len() {
        let current_score = scores[indices[i]];
        while i < indices.len() && scores[indices[i]] == current_score {
            if labels[indices[i]] {
                tp += 1;
            } else {
                fp += 1;
            }
            i += 1;
        }

        let precision = tp as f64 / (tp + fp) as f64;
        let recall = tp as f64 / p;
        points.push(PrPoint {
            threshold: current_score,
            precision,
            recall,
        });
    }

    // AUC via trapezoidal rule over (recall, precision)
    let auc = trapezoidal_auc(
        &points.iter().map(|p| p.recall).collect::<Vec<_>>(),
        &points.iter().map(|p| p.precision).collect::<Vec<_>>(),
    );

    Ok(PrCurve { points, auc })
}

/// Compute only the AUC of the precision-recall curve.
///
/// Shorthand for `pr_curve(scores, labels)?.auc`.
pub fn pr_auc(scores: &[f64], labels: &[bool]) -> Result<f64> {
    Ok(pr_curve(scores, labels)?.auc)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Trapezoidal AUC: sum of trapezoids between consecutive (x, y) points.
fn trapezoidal_auc(x: &[f64], y: &[f64]) -> f64 {
    let mut auc = 0.0;
    for i in 1..x.len() {
        auc += (x[i] - x[i - 1]).abs() * (y[i] + y[i - 1]) / 2.0;
    }
    auc
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Confusion Matrix ────────────────────────────────────────

    #[test]
    fn cm_binary_basic() {
        // actual:    [1, 1, 0, 0, 1, 0]
        // predicted: [1, 0, 0, 1, 1, 0]
        let actual = [1, 1, 0, 0, 1, 0];
        let predicted = [1, 0, 0, 1, 1, 0];
        let cm = ConfusionMatrix::from_labels(&actual, &predicted, None).unwrap();
        assert_eq!(cm.n_classes, 2);
        assert_eq!(cm.total(), 6);

        // Class 1: TP=2, FP=1, FN=1, TN=2
        assert_eq!(cm.true_positives(1), 2);
        assert_eq!(cm.false_positives(1), 1);
        assert_eq!(cm.false_negatives(1), 1);
        assert_eq!(cm.true_negatives(1), 2);
    }

    #[test]
    fn cm_multiclass_3class() {
        let actual = [0, 0, 1, 1, 2, 2];
        let predicted = [0, 1, 1, 2, 2, 0];
        let cm = ConfusionMatrix::from_labels(&actual, &predicted, None).unwrap();
        assert_eq!(cm.n_classes, 3);
        assert_eq!(cm.total(), 6);
        // Class 0: predicted as 0 correctly once
        assert_eq!(cm.true_positives(0), 1);
        // Class 1: predicted as 1 correctly once
        assert_eq!(cm.true_positives(1), 1);
        // Class 2: predicted as 2 correctly once
        assert_eq!(cm.true_positives(2), 1);
    }

    #[test]
    fn cm_perfect_accuracy() {
        let labels = [0, 1, 2, 0, 1, 2];
        let cm = ConfusionMatrix::from_labels(&labels, &labels, None).unwrap();
        assert!((cm.accuracy() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cm_half_accuracy() {
        let actual = [0, 0, 1, 1];
        let predicted = [0, 1, 0, 1];
        let cm = ConfusionMatrix::from_labels(&actual, &predicted, None).unwrap();
        assert!((cm.accuracy() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn cm_precision_recall_f1() {
        // Binary: actual [1,1,1,0,0], predicted [1,1,0,1,0]
        // Class 1: TP=2, FP=1, FN=1
        let actual = [1, 1, 1, 0, 0];
        let predicted = [1, 1, 0, 1, 0];
        let cm = ConfusionMatrix::from_labels(&actual, &predicted, None).unwrap();
        let prec = cm.precision(1);
        let rec = cm.recall(1);
        let f1 = cm.f1(1);
        assert!((prec - 2.0 / 3.0).abs() < 1e-12);
        assert!((rec - 2.0 / 3.0).abs() < 1e-12);
        assert!((f1 - 2.0 / 3.0).abs() < 1e-12); // precision == recall → F1 == same
    }

    #[test]
    fn cm_specificity() {
        // actual [1,1,0,0], predicted [1,0,0,0]
        // Class 1: TN=2, FP=0 → specificity = 1.0
        let actual = [1, 1, 0, 0];
        let predicted = [1, 0, 0, 0];
        let cm = ConfusionMatrix::from_labels(&actual, &predicted, None).unwrap();
        assert!((cm.specificity(1) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cm_empty_error() {
        let empty: &[usize] = &[];
        assert!(ConfusionMatrix::from_labels(empty, empty, None).is_err());
    }

    #[test]
    fn cm_length_mismatch_error() {
        assert!(ConfusionMatrix::from_labels(&[0, 1], &[0], None).is_err());
    }

    // ── Scalar metrics ──────────────────────────────────────────

    #[test]
    fn accuracy_known() {
        let actual = [0, 0, 1, 1, 2, 2];
        let predicted = [0, 0, 1, 1, 2, 2];
        assert!((accuracy(&actual, &predicted).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn f1_binary_known() {
        // Class 1: TP=2, FP=1, FN=1 → P=2/3, R=2/3, F1=2/3
        let actual = [1, 1, 1, 0, 0];
        let predicted = [1, 1, 0, 1, 0];
        let f1 = f1_score(&actual, &predicted, 1).unwrap();
        assert!((f1 - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn f1_macro_multiclass() {
        let actual = [0, 0, 1, 1, 2, 2];
        let predicted = [0, 0, 1, 1, 2, 2];
        let f1 = f1_macro(&actual, &predicted).unwrap();
        assert!((f1 - 1.0).abs() < 1e-12); // perfect
    }

    #[test]
    fn f1_weighted_multiclass() {
        let actual = [0, 0, 1, 1, 2, 2];
        let predicted = [0, 0, 1, 1, 2, 2];
        let f1 = f1_weighted(&actual, &predicted).unwrap();
        assert!((f1 - 1.0).abs() < 1e-12); // perfect
    }

    #[test]
    fn mcc_known_values() {
        // Binary: actual [1,1,0,0], predicted [1,0,0,1]
        // TP=1, TN=1, FP=1, FN=1 → MCC = 0
        let actual = [1, 1, 0, 0];
        let predicted = [1, 0, 0, 1];
        let mcc = matthews_corrcoef(&actual, &predicted).unwrap();
        assert!(mcc.abs() < 1e-12);
    }

    // ── MCC edge cases ──────────────────────────────────────────

    #[test]
    fn mcc_perfect() {
        let actual = [0, 0, 1, 1];
        let predicted = [0, 0, 1, 1];
        let mcc = matthews_corrcoef(&actual, &predicted).unwrap();
        assert!((mcc - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mcc_inverse() {
        let actual = [0, 0, 1, 1];
        let predicted = [1, 1, 0, 0];
        let mcc = matthews_corrcoef(&actual, &predicted).unwrap();
        assert!((mcc - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn mcc_zero_denominator() {
        // All predicted as same class → column for class 1 is all zeros
        let actual = [0, 0, 1, 1];
        let predicted = [0, 0, 0, 0];
        let mcc = matthews_corrcoef(&actual, &predicted).unwrap();
        assert!((mcc - 0.0).abs() < 1e-12);
    }

    // ── ROC Curve ───────────────────────────────────────────────

    #[test]
    fn roc_perfect_auc() {
        // Perfect separation: all positives have higher scores than negatives
        let scores = vec![0.9, 0.8, 0.3, 0.1];
        let labels = vec![true, true, false, false];
        let roc = roc_curve(&scores, &labels).unwrap();
        assert!((roc.auc - 1.0).abs() < 1e-12);
    }

    #[test]
    fn roc_random_auc() {
        // Anti-perfect: all positives scored lower → AUC = 0
        let scores = vec![0.1, 0.2, 0.8, 0.9];
        let labels = vec![true, true, false, false];
        let roc = roc_curve(&scores, &labels).unwrap();
        assert!((roc.auc - 0.0).abs() < 1e-12);
    }

    #[test]
    fn roc_known_curve_points() {
        // Scores: [0.9, 0.7, 0.5, 0.3]  Labels: [T, F, T, F]
        // After sorting: (0.9,T), (0.7,F), (0.5,T), (0.3,F)
        // Thresholds at 0.9: TP=1, FP=0 → (0.0, 0.5)
        // Thresholds at 0.7: TP=1, FP=1 → (0.5, 0.5)
        // Thresholds at 0.5: TP=2, FP=1 → (0.5, 1.0)
        // Thresholds at 0.3: TP=2, FP=2 → (1.0, 1.0)
        let scores = vec![0.9, 0.7, 0.5, 0.3];
        let labels = vec![true, false, true, false];
        let roc = roc_curve(&scores, &labels).unwrap();
        assert!(roc.points.len() >= 3);
        // AUC should be 0.75
        assert!((roc.auc - 0.75).abs() < 1e-12);
    }

    #[test]
    fn roc_endpoints() {
        let scores = vec![0.9, 0.1];
        let labels = vec![true, false];
        let roc = roc_curve(&scores, &labels).unwrap();
        // First point should be (0, 0)
        assert!((roc.points[0].fpr - 0.0).abs() < 1e-12);
        assert!((roc.points[0].tpr - 0.0).abs() < 1e-12);
        // Last point should be (1, 1)
        let last = roc.points.last().unwrap();
        assert!((last.fpr - 1.0).abs() < 1e-12);
        assert!((last.tpr - 1.0).abs() < 1e-12);
    }

    #[test]
    fn roc_all_same_label_error() {
        let scores = vec![0.9, 0.8, 0.7];
        let labels = vec![true, true, true];
        assert!(roc_curve(&scores, &labels).is_err());
    }

    #[test]
    fn roc_empty_error() {
        let scores: Vec<f64> = vec![];
        let labels: Vec<bool> = vec![];
        assert!(roc_curve(&scores, &labels).is_err());
    }

    // ── PR Curve ────────────────────────────────────────────────

    #[test]
    fn pr_perfect_auc() {
        let scores = vec![0.9, 0.8, 0.3, 0.1];
        let labels = vec![true, true, false, false];
        let pr = pr_curve(&scores, &labels).unwrap();
        assert!((pr.auc - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pr_known_curve_points() {
        // Scores: [0.9, 0.7, 0.5, 0.3]  Labels: [T, F, T, F]
        // After sorting: (0.9,T), (0.7,F), (0.5,T), (0.3,F)
        // At 0.9: TP=1, FP=0 → P=1.0, R=0.5
        // At 0.7: TP=1, FP=1 → P=0.5, R=0.5
        // At 0.5: TP=2, FP=1 → P=2/3, R=1.0
        // At 0.3: TP=2, FP=2 → P=0.5, R=1.0
        let scores = vec![0.9, 0.7, 0.5, 0.3];
        let labels = vec![true, false, true, false];
        let pr = pr_curve(&scores, &labels).unwrap();
        assert!(pr.points.len() >= 3);
        // Check a few points (skipping the initial anchor at recall=0, precision=1)
        let p1 = &pr.points[1]; // threshold=0.9
        assert!((p1.precision - 1.0).abs() < 1e-12);
        assert!((p1.recall - 0.5).abs() < 1e-12);
    }

    #[test]
    fn pr_all_positive_edge() {
        // All labels positive — no negative class
        let scores = vec![0.9, 0.8, 0.7];
        let labels = vec![true, true, true];
        let pr = pr_curve(&scores, &labels).unwrap();
        // Precision is always 1.0 when all are positive
        assert!((pr.auc - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pr_no_positives_error() {
        let scores = vec![0.5, 0.3];
        let labels = vec![false, false];
        assert!(pr_curve(&scores, &labels).is_err());
    }

    #[test]
    fn pr_auc_matches_curve() {
        let scores = vec![0.9, 0.7, 0.5, 0.3];
        let labels = vec![true, false, true, false];
        let from_curve = pr_curve(&scores, &labels).unwrap().auc;
        let shortcut = pr_auc(&scores, &labels).unwrap();
        assert!((from_curve - shortcut).abs() < 1e-12);
    }
}
