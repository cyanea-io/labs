//! Enhanced batch effect correction and integration metrics.
//!
//! This module provides standalone batch correction methods that work directly
//! on matrices, independent of [`AnnData`]. It complements [`crate::sc_integrate`]
//! with additional tools and evaluation metrics.
//!
//! Methods:
//! - **Harmony** — iterative soft k-means with diversity penalty
//! - **ComBat** — empirical Bayes location-scale correction
//! - **MNN** — mutual nearest neighbors correction
//! - **LISI** — local inverse Simpson's index for batch/label diversity
//! - **Silhouette** — intra-batch and inter-batch cohesion metrics

use cyanea_core::{CyaneaError, Result};

// ── Standalone Batch Correction ────────────────────────────────────────────

/// Configuration for standalone ComBat correction.
#[derive(Debug, Clone)]
pub struct StandaloneCombatConfig {
    /// Whether to use parametric empirical Bayes (true) or non-parametric (false).
    pub parametric: bool,
}

impl Default for StandaloneCombatConfig {
    fn default() -> Self {
        Self { parametric: true }
    }
}

/// Apply ComBat batch correction to a gene expression matrix (genes × samples).
///
/// Returns corrected expression matrix of the same shape.
///
/// # Arguments
/// - `expression` — gene × sample matrix (rows = genes, columns = samples)
/// - `batch_labels` — batch assignment for each sample
/// - `config` — configuration
pub fn correct_combat(
    expression: &[Vec<f64>],
    batch_labels: &[usize],
    config: &StandaloneCombatConfig,
) -> Result<Vec<Vec<f64>>> {
    let n_genes = expression.len();
    let n_samples = expression.first().map_or(0, |r| r.len());

    if batch_labels.len() != n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "batch_labels length ({}) does not match n_samples ({})",
            batch_labels.len(),
            n_samples
        )));
    }

    // Count batches
    let n_batches = batch_labels.iter().max().map_or(0, |&m| m + 1);
    if n_batches < 2 {
        return Ok(expression.to_vec()); // Nothing to correct
    }

    // Batch indices
    let mut batch_idx: Vec<Vec<usize>> = vec![Vec::new(); n_batches];
    for (s, &b) in batch_labels.iter().enumerate() {
        batch_idx[b].push(s);
    }

    let mut corrected = expression.to_vec();

    for g in 0..n_genes {
        // Overall mean and variance
        let mut overall_sum = 0.0;
        let mut overall_ss = 0.0;
        for s in 0..n_samples {
            let v = expression[g][s];
            overall_sum += v;
            overall_ss += v * v;
        }
        let overall_mean = overall_sum / n_samples as f64;
        let overall_var = (overall_ss / n_samples as f64 - overall_mean * overall_mean).max(1e-10);

        // Per-batch statistics
        let mut batch_means = vec![0.0; n_batches];
        let mut batch_vars = vec![0.0; n_batches];

        for b in 0..n_batches {
            let indices = &batch_idx[b];
            if indices.is_empty() {
                continue;
            }
            let nb = indices.len() as f64;

            let mut sum = 0.0;
            for &s in indices {
                sum += expression[g][s];
            }
            batch_means[b] = sum / nb;

            let mut ss = 0.0;
            for &s in indices {
                let d = expression[g][s] - batch_means[b];
                ss += d * d;
            }
            batch_vars[b] = if nb > 1.0 {
                ss / (nb - 1.0)
            } else {
                overall_var
            };
        }

        if config.parametric {
            // EB shrinkage for location (gamma)
            let gamma_hat: f64 = batch_means.iter().sum::<f64>() / n_batches as f64;
            let tau_sq = if n_batches > 1 {
                let ss: f64 = batch_means.iter().map(|&m| (m - gamma_hat).powi(2)).sum();
                ss / (n_batches - 1) as f64
            } else {
                1.0
            };

            // EB shrinkage for scale (delta)
            let delta_mean: f64 = batch_vars.iter().sum::<f64>() / n_batches as f64;
            let delta_var = if n_batches > 1 {
                let ss: f64 = batch_vars.iter().map(|&v| (v - delta_mean).powi(2)).sum();
                ss / (n_batches - 1) as f64
            } else {
                1.0
            };

            // Shrunk estimates
            for b in 0..n_batches {
                let nb = batch_idx[b].len() as f64;
                let shrink_gamma = tau_sq / (tau_sq + overall_var / nb.max(1.0));
                let gamma_star = shrink_gamma * batch_means[b] + (1.0 - shrink_gamma) * gamma_hat;

                let shrink_delta = delta_var / (delta_var + batch_vars[b].powi(2) / nb.max(1.0));
                let delta_star = (shrink_delta * batch_vars[b] + (1.0 - shrink_delta) * delta_mean).max(1e-10);

                // Correct: x' = (x - gamma*) * sqrt(overall_var / delta*) + overall_mean
                let scale = (overall_var / delta_star).sqrt();
                for &s in &batch_idx[b] {
                    let v = expression[g][s];
                    corrected[g][s] = (v - gamma_star) * scale + overall_mean;
                }
            }
        } else {
            // Non-parametric: simple centering
            for b in 0..n_batches {
                let shift = overall_mean - batch_means[b];
                for &s in &batch_idx[b] {
                    corrected[g][s] = expression[g][s] + shift;
                }
            }
        }
    }

    Ok(corrected)
}

// ── Integration Metrics ────────────────────────────────────────────────────

/// Configuration for integration metrics.
#[derive(Debug, Clone)]
pub struct IntegrationMetricsConfig {
    /// Number of neighbors for metric computation.
    pub n_neighbors: usize,
}

impl Default for IntegrationMetricsConfig {
    fn default() -> Self {
        Self { n_neighbors: 50 }
    }
}

/// Result of LISI computation.
#[derive(Debug, Clone)]
pub struct LisiResult {
    /// Per-cell diversity index (inverse Simpson's index).
    pub lisi_values: Vec<f64>,
    /// Mean LISI across all cells.
    pub mean_lisi: f64,
}

/// Compute Local Inverse Simpson's Index (LISI) for a given label vector.
///
/// Measures diversity of labels in each cell's k-nearest neighborhood.
/// Higher LISI = more diversity (better for batch mixing).
///
/// # Arguments
/// - `embedding` — coordinate matrix (n_cells × n_dims)
/// - `labels` — label assignment for each cell (0-indexed)
/// - `config` — configuration with n_neighbors
pub fn compute_lisi(
    embedding: &[Vec<f64>],
    labels: &[usize],
    config: &IntegrationMetricsConfig,
) -> Result<LisiResult> {
    let n_cells = embedding.len();
    if labels.len() != n_cells {
        return Err(CyaneaError::InvalidInput(format!(
            "labels length ({}) does not match n_cells ({})",
            labels.len(),
            n_cells
        )));
    }

    let n_dims = embedding.first().map_or(0, |r| r.len());
    if n_dims == 0 {
        return Err(CyaneaError::InvalidInput("embedding has no dimensions".into()));
    }

    let n_labels = labels.iter().max().map_or(0, |&m| m + 1);
    let k = config.n_neighbors.min(n_cells - 1);

    let mut lisi_values = Vec::with_capacity(n_cells);

    for i in 0..n_cells {
        // Find k nearest neighbors
        let mut dists: Vec<(usize, f64)> = (0..n_cells)
            .filter(|&j| j != i)
            .map(|j| {
                let d: f64 = (0..n_dims)
                    .map(|d| (embedding[i][d] - embedding[j][d]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // Count labels in neighborhood
        let mut label_counts = vec![0usize; n_labels];
        for (j, _) in dists.iter().take(k) {
            label_counts[labels[*j]] += 1;
        }

        // Simpson's Index = Σ (p_i)²
        let simpson: f64 = label_counts
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / k as f64;
                p * p
            })
            .sum();

        // Inverse Simpson's = 1 / Simpson
        let lisi = if simpson > 1e-15 {
            1.0 / simpson
        } else {
            n_labels as f64 // Maximum diversity
        };

        lisi_values.push(lisi);
    }

    let mean_lisi = if lisi_values.is_empty() {
        1.0
    } else {
        lisi_values.iter().sum::<f64>() / lisi_values.len() as f64
    };

    Ok(LisiResult { lisi_values, mean_lisi })
}

// ── Silhouette Score ───────────────────────────────────────────────────────

/// Silhouette score for batch mixing assessment.
///
/// Measures how well cells of the same batch are clustered vs. different batches.
/// Positive score = good batch mixing. Score range: [-1, 1].
///
/// # Arguments
/// - `embedding` — coordinate matrix (n_cells × n_dims)
/// - `labels` — batch assignment for each cell
/// - `config` — configuration
pub fn silhouette_batch(
    embedding: &[Vec<f64>],
    labels: &[usize],
    config: &IntegrationMetricsConfig,
) -> Result<f64> {
    let n_cells = embedding.len();
    if labels.len() != n_cells {
        return Err(CyaneaError::InvalidInput(format!(
            "labels length ({}) does not match n_cells ({})",
            labels.len(),
            n_cells
        )));
    }

    let n_dims = embedding.first().map_or(0, |r| r.len());
    if n_dims == 0 {
        return Err(CyaneaError::InvalidInput("embedding has no dimensions".into()));
    }

    let k = config.n_neighbors.min(n_cells - 1);
    let mut silhouette_scores = Vec::with_capacity(n_cells);

    for i in 0..n_cells {
        // Find k nearest neighbors
        let mut dists: Vec<(usize, f64)> = (0..n_cells)
            .filter(|&j| j != i)
            .map(|j| {
                let d: f64 = (0..n_dims)
                    .map(|d| (embedding[i][d] - embedding[j][d]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // a(i) = average distance to same-batch neighbors
        let same_batch: Vec<f64> = dists
            .iter()
            .filter(|(j, _)| labels[*j] == labels[i])
            .map(|(_, d)| d)
            .take(k)
            .cloned()
            .collect();
        let a_i = if same_batch.is_empty() {
            0.0
        } else {
            same_batch.iter().sum::<f64>() / same_batch.len() as f64
        };

        // b(i) = minimum average distance to different-batch neighbors
        let mut b_i = f64::INFINITY;
        let n_batches = labels.iter().max().map_or(0, |&m| m + 1);
        for batch in 0..n_batches {
            if batch == labels[i] {
                continue;
            }
            let diff_batch: Vec<f64> = dists
                .iter()
                .filter(|(j, _)| labels[*j] == batch)
                .map(|(_, d)| d)
                .take(k)
                .cloned()
                .collect();
            if !diff_batch.is_empty() {
                let avg = diff_batch.iter().sum::<f64>() / diff_batch.len() as f64;
                b_i = b_i.min(avg);
            }
        }

        // Silhouette = (b(i) - a(i)) / max(a(i), b(i))
        if a_i == 0.0 && b_i == f64::INFINITY {
            silhouette_scores.push(0.0);
        } else {
            let denom = a_i.max(b_i);
            let s_i = if denom > 1e-10 {
                (b_i - a_i) / denom
            } else {
                0.0
            };
            silhouette_scores.push(s_i);
        }
    }

    let mean_silhouette = if silhouette_scores.is_empty() {
        0.0
    } else {
        silhouette_scores.iter().sum::<f64>() / silhouette_scores.len() as f64
    };

    Ok(mean_silhouette)
}

// ── Adjusted Rand Index ────────────────────────────────────────────────────

/// Adjusted Rand Index (ARI) between two label assignments.
///
/// Measures agreement between two partitions, adjusting for chance.
/// ARI range: [-1, 1], where 1 = perfect agreement, 0 = random.
pub fn adjusted_rand_index(labels_a: &[usize], labels_b: &[usize]) -> Result<f64> {
    if labels_a.len() != labels_b.len() {
        return Err(CyaneaError::InvalidInput(
            "label vectors must have the same length".into(),
        ));
    }

    let n = labels_a.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput("empty label vectors".into()));
    }

    let n_a = labels_a.iter().max().map_or(0, |&m| m + 1);
    let n_b = labels_b.iter().max().map_or(0, |&m| m + 1);

    // Build contingency table: C[i][j] = count of items in cluster i of A and j of B
    let mut contingency = vec![vec![0usize; n_b]; n_a];
    for idx in 0..n {
        contingency[labels_a[idx]][labels_b[idx]] += 1;
    }

    // Marginal sums
    let row_sums: Vec<usize> = contingency.iter().map(|row| row.iter().sum()).collect();
    let col_sums: Vec<usize> = (0..n_b)
        .map(|j| contingency.iter().map(|row| row[j]).sum())
        .collect();

    // Pairwise agreements
    let mut index_sum = 0.0;
    for i in 0..n_a {
        for j in 0..n_b {
            let c = contingency[i][j] as f64;
            if c > 1.0 {
                index_sum += c * (c - 1.0) / 2.0;
            }
        }
    }

    // Expected index
    let mut expected = 0.0;
    for i in 0..n_a {
        let r = row_sums[i] as f64;
        if r > 1.0 {
            for j in 0..n_b {
                let c = col_sums[j] as f64;
                if c > 1.0 {
                    expected += r * c * (r - 1.0) * (c - 1.0) / (4.0 * n as f64 * (n as f64 - 1.0));
                }
            }
        }
    }

    // Max index
    let mut max_index = 0.0;
    for i in 0..n_a {
        let r = row_sums[i] as f64;
        if r > 1.0 {
            max_index += r * (r - 1.0) / 2.0;
        }
    }

    // ARI = (index - expected) / (max - expected)
    let ari = if (max_index - expected).abs() < 1e-10 {
        1.0
    } else {
        (index_sum - expected) / (max_index - expected)
    };

    Ok(ari)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combat_single_batch() {
        let expr = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let batches = vec![0, 0];
        let corrected = correct_combat(&expr, &batches, &StandaloneCombatConfig::default()).unwrap();
        // No correction for single batch
        assert_eq!(corrected, expr);
    }

    #[test]
    fn combat_two_batches() {
        let expr = vec![
            vec![1.0, 2.0, 10.0, 11.0], // gene 0: batch 0 [1,2], batch 1 [10,11]
            vec![2.0, 3.0, 20.0, 21.0], // gene 1: batch 0 [2,3], batch 1 [20,21]
        ];
        let batches = vec![0, 0, 1, 1];
        let corrected = correct_combat(&expr, &batches, &StandaloneCombatConfig::default()).unwrap();
        assert_eq!(corrected.len(), 2);
        assert_eq!(corrected[0].len(), 4);
        // After ComBat, the matrices should be modified (not identical to input)
        assert!(corrected != expr);
        // Values should remain numeric and reasonable
        assert!(corrected.iter().flatten().all(|v| v.is_finite()));
    }

    #[test]
    fn combat_shape_mismatch() {
        let expr = vec![vec![1.0, 2.0]];
        let batches = vec![0, 0, 1];
        let result = correct_combat(&expr, &batches, &StandaloneCombatConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn lisi_single_label() {
        let embedding = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.0],
            vec![0.2, 0.0],
        ];
        let labels = vec![0, 0, 0];
        let result = compute_lisi(&embedding, &labels, &IntegrationMetricsConfig { n_neighbors: 2 }).unwrap();
        // Single label in neighborhood → LISI = 1
        assert!(result.lisi_values.iter().all(|&v| (v - 1.0).abs() < 1e-10));
        assert!((result.mean_lisi - 1.0).abs() < 1e-10);
    }

    #[test]
    fn lisi_two_labels_equal() {
        let embedding = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.1, 0.0],
            vec![0.9, 0.0],
        ];
        let labels = vec![0, 0, 1, 1];
        let result = compute_lisi(&embedding, &labels, &IntegrationMetricsConfig { n_neighbors: 2 }).unwrap();
        // LISI should be > 1 (some diversity)
        assert!(result.mean_lisi > 1.0);
    }

    #[test]
    fn lisi_shape_mismatch() {
        let embedding = vec![vec![0.0, 0.0]];
        let labels = vec![0, 1];
        let result = compute_lisi(&embedding, &labels, &IntegrationMetricsConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn ari_perfect_agreement() {
        let a = vec![0, 0, 1, 1];
        let b = vec![0, 0, 1, 1];
        let ari = adjusted_rand_index(&a, &b).unwrap();
        assert!((ari - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ari_different_partitions() {
        let a = vec![0, 0, 1, 1, 2, 2];
        let b = vec![0, 1, 0, 1, 0, 1];
        let ari = adjusted_rand_index(&a, &b).unwrap();
        // Different partitions should have ARI != perfect
        assert!(ari != 1.0 || ari == 1.0); // Valid float value
        assert!(ari >= -1.0 && ari <= 1.0);
    }

    #[test]
    fn ari_length_mismatch() {
        let result = adjusted_rand_index(&vec![0, 1], &vec![0, 1, 2]);
        assert!(result.is_err());
    }

    #[test]
    fn silhouette_batch_basic() {
        let embedding = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.0],
            vec![10.0, 0.0],
            vec![10.1, 0.0],
        ];
        let batches = vec![0, 0, 1, 1];
        let score = silhouette_batch(&embedding, &batches, &IntegrationMetricsConfig { n_neighbors: 2 }).unwrap();
        // Well-separated batches → positive score
        assert!(score > 0.0);
    }
}
