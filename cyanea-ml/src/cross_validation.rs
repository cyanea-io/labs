//! Cross-validation utilities.
//!
//! Provides k-fold, stratified k-fold, and leave-one-out cross-validation.
//! All functions are closure-based: the caller supplies a function
//! `Fn(&[usize], &[usize]) -> Result<f64>` that receives train/test index
//! splits and returns a score (e.g., accuracy). This design accommodates
//! classifiers with different `fit()` signatures.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result for a single fold.
#[derive(Debug, Clone)]
pub struct FoldResult {
    /// Fold number (0-indexed).
    pub fold: usize,
    /// Number of training samples.
    pub n_train: usize,
    /// Number of test samples.
    pub n_test: usize,
    /// Evaluation score returned by the user's closure.
    pub score: f64,
}

/// Aggregated cross-validation result.
#[derive(Debug, Clone)]
pub struct CvResult {
    /// Per-fold results.
    pub folds: Vec<FoldResult>,
    /// Mean score across folds.
    pub mean_score: f64,
    /// Standard deviation of scores across folds.
    pub std_score: f64,
}

// ---------------------------------------------------------------------------
// PRNG (copied from forest.rs — private per-module, matches push_cigar pattern)
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

/// Fisher-Yates shuffle using LCG.
fn shuffle(rng: &mut LcgRng, data: &mut [usize]) {
    let n = data.len();
    for i in (1..n).rev() {
        let j = rng.next_bounded((i + 1) as u64) as usize;
        data.swap(i, j);
    }
}

// ---------------------------------------------------------------------------
// K-Fold Cross-Validation
// ---------------------------------------------------------------------------

/// K-fold cross-validation.
///
/// Splits `n_samples` indices into `k` folds (after shuffling with `seed`),
/// then calls `eval_fn(train_indices, test_indices)` for each fold.
///
/// # Errors
///
/// Returns an error if `k < 2`, `k > n_samples`, or the evaluation closure
/// returns an error.
pub fn cross_validate_kfold<F>(
    n_samples: usize,
    k: usize,
    seed: u64,
    mut eval_fn: F,
) -> Result<CvResult>
where
    F: FnMut(&[usize], &[usize]) -> Result<f64>,
{
    if k < 2 {
        return Err(CyaneaError::InvalidInput(
            "k must be at least 2".into(),
        ));
    }
    if k > n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "k ({}) > n_samples ({})",
            k, n_samples
        )));
    }

    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = LcgRng::new(seed);
    shuffle(&mut rng, &mut indices);

    let folds = build_folds(&indices, k);
    run_folds(&folds, k, &mut eval_fn)
}

// ---------------------------------------------------------------------------
// Stratified K-Fold Cross-Validation
// ---------------------------------------------------------------------------

/// Stratified k-fold cross-validation.
///
/// Preserves the class distribution in each fold. Groups indices by label,
/// shuffles each group, and distributes round-robin across `k` folds.
///
/// # Errors
///
/// Returns an error if `k < 2`, `k > n_samples`, or the evaluation closure
/// returns an error.
pub fn cross_validate_stratified<F>(
    labels: &[usize],
    k: usize,
    seed: u64,
    mut eval_fn: F,
) -> Result<CvResult>
where
    F: FnMut(&[usize], &[usize]) -> Result<f64>,
{
    let n_samples = labels.len();
    if k < 2 {
        return Err(CyaneaError::InvalidInput(
            "k must be at least 2".into(),
        ));
    }
    if k > n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "k ({}) > n_samples ({})",
            k, n_samples
        )));
    }

    let n_classes = labels.iter().copied().max().map_or(0, |m| m + 1);

    // Group indices by class
    let mut class_indices: Vec<Vec<usize>> = vec![Vec::new(); n_classes];
    for (i, &label) in labels.iter().enumerate() {
        class_indices[label].push(i);
    }

    // Shuffle each class's indices
    let mut rng = LcgRng::new(seed);
    for group in &mut class_indices {
        shuffle(&mut rng, group);
    }

    // Distribute round-robin across folds
    let mut fold_indices: Vec<Vec<usize>> = vec![Vec::new(); k];
    for group in &class_indices {
        for (i, &idx) in group.iter().enumerate() {
            fold_indices[i % k].push(idx);
        }
    }

    run_folds(&fold_indices, k, &mut eval_fn)
}

// ---------------------------------------------------------------------------
// Leave-One-Out Cross-Validation
// ---------------------------------------------------------------------------

/// Leave-one-out cross-validation (LOO-CV).
///
/// Each sample is used once as the test set, with all remaining samples as
/// training. Equivalent to k-fold with `k = n_samples`.
///
/// # Errors
///
/// Returns an error if `n_samples < 2` or the evaluation closure returns
/// an error.
pub fn cross_validate_loo<F>(n_samples: usize, mut eval_fn: F) -> Result<CvResult>
where
    F: FnMut(&[usize], &[usize]) -> Result<f64>,
{
    if n_samples < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 samples for LOO".into(),
        ));
    }

    let mut folds = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let train: Vec<usize> = (0..n_samples).filter(|&j| j != i).collect();
        let test = vec![i];
        let score = eval_fn(&train, &test)?;
        folds.push(FoldResult {
            fold: i,
            n_train: n_samples - 1,
            n_test: 1,
            score,
        });
    }

    let (mean, std) = score_stats(&folds);
    Ok(CvResult {
        folds,
        mean_score: mean,
        std_score: std,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Split a flat list of indices into `k` fold groups.
fn build_folds(indices: &[usize], k: usize) -> Vec<Vec<usize>> {
    let n = indices.len();
    let base_size = n / k;
    let remainder = n % k;

    let mut folds = Vec::with_capacity(k);
    let mut start = 0;
    for i in 0..k {
        let size = base_size + if i < remainder { 1 } else { 0 };
        folds.push(indices[start..start + size].to_vec());
        start += size;
    }
    folds
}

/// Run the eval function for each fold (test = fold i, train = rest).
fn run_folds<F>(fold_indices: &[Vec<usize>], k: usize, eval_fn: &mut F) -> Result<CvResult>
where
    F: FnMut(&[usize], &[usize]) -> Result<f64>,
{
    let mut results = Vec::with_capacity(k);
    for i in 0..k {
        let test = &fold_indices[i];
        let train: Vec<usize> = fold_indices
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .flat_map(|(_, v)| v.iter().copied())
            .collect();

        let score = eval_fn(&train, test)?;
        results.push(FoldResult {
            fold: i,
            n_train: train.len(),
            n_test: test.len(),
            score,
        });
    }

    let (mean, std) = score_stats(&results);
    Ok(CvResult {
        folds: results,
        mean_score: mean,
        std_score: std,
    })
}

/// Mean and standard deviation of fold scores.
fn score_stats(folds: &[FoldResult]) -> (f64, f64) {
    let n = folds.len() as f64;
    let mean = folds.iter().map(|f| f.score).sum::<f64>() / n;
    let var = folds.iter().map(|f| (f.score - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── K-Fold ──────────────────────────────────────────────────

    #[test]
    fn kfold_sizes_correct() {
        let result = cross_validate_kfold(10, 5, 42, |train, test| {
            assert_eq!(train.len(), 8);
            assert_eq!(test.len(), 2);
            Ok(1.0)
        })
        .unwrap();
        assert_eq!(result.folds.len(), 5);
    }

    #[test]
    fn kfold_uneven_split() {
        // 7 samples, 3 folds → sizes 3, 2, 2
        let result = cross_validate_kfold(7, 3, 42, |train, test| {
            assert_eq!(train.len() + test.len(), 7);
            Ok(1.0)
        })
        .unwrap();
        assert_eq!(result.folds.len(), 3);
        let test_sizes: Vec<usize> = result.folds.iter().map(|f| f.n_test).collect();
        assert!(test_sizes.contains(&3) || test_sizes.contains(&2));
    }

    #[test]
    fn kfold_no_overlap() {
        let mut all_test: Vec<usize> = Vec::new();
        cross_validate_kfold(10, 5, 42, |_train, test| {
            all_test.extend_from_slice(test);
            Ok(1.0)
        })
        .unwrap();
        all_test.sort();
        assert_eq!(all_test, (0..10).collect::<Vec<_>>());
    }

    #[test]
    fn kfold_k_gt_n_error() {
        let result = cross_validate_kfold(3, 5, 42, |_, _| Ok(1.0));
        assert!(result.is_err());
    }

    // ── Stratified K-Fold ───────────────────────────────────────

    #[test]
    fn stratified_class_proportions() {
        // 6 class-0, 4 class-1, k=2 → each fold should have ~3 class-0, ~2 class-1
        let labels = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1];
        cross_validate_stratified(&labels, 2, 42, |_train, test| {
            let n0 = test.iter().filter(|&&i| labels[i] == 0).count();
            let n1 = test.iter().filter(|&&i| labels[i] == 1).count();
            assert_eq!(n0, 3);
            assert_eq!(n1, 2);
            Ok(1.0)
        })
        .unwrap();
    }

    #[test]
    fn stratified_all_classes_in_each_fold() {
        let labels = vec![0, 0, 1, 1, 2, 2, 0, 1, 2];
        cross_validate_stratified(&labels, 3, 42, |_train, test| {
            let classes: std::collections::HashSet<usize> =
                test.iter().map(|&i| labels[i]).collect();
            assert_eq!(classes.len(), 3, "not all classes in fold");
            Ok(1.0)
        })
        .unwrap();
    }

    #[test]
    fn stratified_deterministic_with_seed() {
        let labels = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
        let dummy = |_train: &[usize], test: &[usize]| Ok(test.len() as f64);
        let r1 = cross_validate_stratified(&labels, 3, 42, dummy).unwrap();
        let r2 = cross_validate_stratified(&labels, 3, 42, dummy).unwrap();
        for (f1, f2) in r1.folds.iter().zip(r2.folds.iter()) {
            assert!((f1.score - f2.score).abs() < 1e-12);
        }
    }

    #[test]
    fn stratified_covers_all_samples() {
        let labels = vec![0, 0, 1, 1, 2, 2];
        let mut all_test: Vec<usize> = Vec::new();
        cross_validate_stratified(&labels, 3, 42, |_train, test| {
            all_test.extend_from_slice(test);
            Ok(1.0)
        })
        .unwrap();
        all_test.sort();
        assert_eq!(all_test, (0..6).collect::<Vec<_>>());
    }

    // ── Leave-One-Out ───────────────────────────────────────────

    #[test]
    fn loo_fold_count_and_sizes() {
        let result = cross_validate_loo(5, |train, test| {
            assert_eq!(train.len(), 4);
            assert_eq!(test.len(), 1);
            Ok(1.0)
        })
        .unwrap();
        assert_eq!(result.folds.len(), 5);
    }

    #[test]
    fn loo_trivially_separable() {
        // Two distinct classes, well separated — predict by nearest neighbor logic
        // Data: class 0 at x=0, class 1 at x=100
        let data = vec![0.0, 0.0, 100.0, 100.0];
        let labels = vec![0, 0, 1, 1];
        let n_features = 1;

        let result = cross_validate_loo(4, |train, test| {
            use crate::tree::DecisionTree;
            let train_data: Vec<f64> = train
                .iter()
                .map(|&i| data[i * n_features..(i + 1) * n_features].to_vec())
                .flatten()
                .collect();
            let train_labels: Vec<usize> = train.iter().map(|&i| labels[i]).collect();
            let tree =
                DecisionTree::fit(&train_data, n_features, &train_labels, 3)?;

            let mut correct = 0;
            for &i in test {
                let sample = &data[i * n_features..(i + 1) * n_features];
                if tree.predict(sample) == labels[i] {
                    correct += 1;
                }
            }
            Ok(correct as f64 / test.len() as f64)
        })
        .unwrap();

        assert!((result.mean_score - 1.0).abs() < 1e-12);
    }

    // ── Integration ─────────────────────────────────────────────

    #[test]
    fn integration_decision_tree_kfold() {
        // Simple linearly separable 2-class dataset
        // Class 0: features around [0, 0], Class 1: features around [10, 10]
        let data: Vec<f64> = vec![
            0.0, 0.1, 0.1, 0.0, 0.2, 0.1, 0.0, 0.2, 0.1, 0.1,
            10.0, 10.1, 10.1, 10.0, 10.2, 10.1, 10.0, 10.2, 10.1, 10.1,
        ];
        let labels: Vec<usize> = vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
        let n_features = 2;

        let result = cross_validate_kfold(10, 5, 42, |train, test| {
            use crate::tree::DecisionTree;
            let train_data: Vec<f64> = train
                .iter()
                .flat_map(|&i| data[i * n_features..(i + 1) * n_features].iter().copied())
                .collect();
            let train_labels: Vec<usize> = train.iter().map(|&i| labels[i]).collect();
            let tree =
                DecisionTree::fit(&train_data, n_features, &train_labels, 5)?;

            let correct: usize = test
                .iter()
                .filter(|&&i| {
                    let sample = &data[i * n_features..(i + 1) * n_features];
                    tree.predict(sample) == labels[i]
                })
                .count();
            Ok(correct as f64 / test.len() as f64)
        })
        .unwrap();

        assert!(
            result.mean_score > 0.8,
            "expected high accuracy on separable data, got {}",
            result.mean_score
        );
    }
}
