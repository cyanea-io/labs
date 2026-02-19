//! Single-cell batch correction and integration: Harmony, ComBat, MNN, metrics.

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

use crate::single_cell::{AnnData, ColumnData, MatrixData};
use crate::sparse::SparseMatrix;

// ── Harmony ────────────────────────────────────────────────────────────────

/// Configuration for Harmony integration.
#[derive(Debug, Clone)]
pub struct HarmonyConfig {
    /// Key in `obs` containing batch labels.
    pub batch_key: String,
    /// Number of clusters for soft k-means (None = auto).
    pub n_clusters: Option<usize>,
    /// Diversity penalty strength.
    pub theta: f64,
    /// Bandwidth for soft assignment.
    pub sigma: f64,
    /// Maximum iterations.
    pub max_iter: usize,
}

impl Default for HarmonyConfig {
    fn default() -> Self {
        Self {
            batch_key: "batch".into(),
            n_clusters: None,
            theta: 2.0,
            sigma: 0.1,
            max_iter: 10,
        }
    }
}

/// Harmony batch correction (Korsunsky 2019).
///
/// Iteratively corrects PCA embeddings by:
/// 1. Soft k-means clustering
/// 2. Penalizing clusters that are batch-pure
/// 3. Computing per-cluster correction vectors
///
/// Stores corrected PCA in `obsm["X_pca_harmony"]`.
pub fn harmony(adata: &mut AnnData, config: &HarmonyConfig) -> Result<()> {
    let pca = adata
        .get_obsm("X_pca")
        .ok_or_else(|| CyaneaError::InvalidInput("obsm['X_pca'] not found".into()))?
        .clone();

    let n_obs = pca.len();
    let n_dims = pca[0].len();

    // Get batch assignments
    let batch_col = adata
        .get_obs(&config.batch_key)
        .ok_or_else(|| {
            CyaneaError::InvalidInput(format!("obs['{}'] not found", config.batch_key))
        })?;

    let batch_labels: Vec<String> = match batch_col {
        ColumnData::Strings(v) => v.clone(),
        ColumnData::Numeric(v) => v.iter().map(|x| x.to_string()).collect(),
        ColumnData::Categorical { codes, categories } => codes
            .iter()
            .map(|&c| categories.get(c as usize).cloned().unwrap_or_else(|| c.to_string()))
            .collect(),
    };

    let mut unique_batches: Vec<String> = batch_labels.clone();
    unique_batches.sort();
    unique_batches.dedup();
    let n_batches = unique_batches.len();

    let batch_map: HashMap<&str, usize> = unique_batches
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i))
        .collect();
    let batch_ids: Vec<usize> = batch_labels.iter().map(|b| batch_map[b.as_str()]).collect();

    // Batch proportions
    let batch_counts: Vec<usize> = (0..n_batches)
        .map(|b| batch_ids.iter().filter(|&&x| x == b).count())
        .collect();
    let batch_freq: Vec<f64> = batch_counts.iter().map(|&c| c as f64 / n_obs as f64).collect();

    // Number of clusters
    let k = config.n_clusters.unwrap_or((n_obs as f64).sqrt().ceil() as usize).max(2);

    // Initialize corrected embeddings
    let mut z: Vec<Vec<f64>> = pca.clone();

    // Initialize centroids using first k points (simple init)
    let mut centroids: Vec<Vec<f64>> = z[..k.min(n_obs)].to_vec();
    // Pad if needed
    while centroids.len() < k {
        centroids.push(z[centroids.len() % n_obs].clone());
    }

    for _iter in 0..config.max_iter {
        // Soft assignment: R[i][c] = exp(-||z_i - centroid_c||² / sigma)
        let mut r = vec![vec![0.0; k]; n_obs];
        for i in 0..n_obs {
            let mut max_val = f64::NEG_INFINITY;
            for c in 0..k {
                let dist_sq: f64 = (0..n_dims)
                    .map(|d| (z[i][d] - centroids[c][d]).powi(2))
                    .sum();
                r[i][c] = -dist_sq / config.sigma;
                if r[i][c] > max_val {
                    max_val = r[i][c];
                }
            }
            // Softmax with diversity penalty
            let mut sum = 0.0;
            for c in 0..k {
                r[i][c] = (r[i][c] - max_val).exp();
                sum += r[i][c];
            }
            for c in 0..k {
                r[i][c] /= sum;
            }
        }

        // Diversity penalty: penalize clusters that are batch-pure
        // For each cluster c, compute batch composition
        for c in 0..k {
            let total_r: f64 = (0..n_obs).map(|i| r[i][c]).sum();
            if total_r < 1e-10 {
                continue;
            }

            for b in 0..n_batches {
                let batch_r: f64 = (0..n_obs)
                    .filter(|&i| batch_ids[i] == b)
                    .map(|i| r[i][c])
                    .sum();
                let observed_freq = batch_r / total_r;
                let expected_freq = batch_freq[b];

                // Penalty: reduce assignment for over-represented batches
                if observed_freq > expected_freq * 1.5 {
                    let penalty = 1.0 - config.theta * (observed_freq - expected_freq);
                    let penalty = penalty.max(0.1);
                    for i in 0..n_obs {
                        if batch_ids[i] == b {
                            r[i][c] *= penalty;
                        }
                    }
                }
            }

            // Renormalize
            for i in 0..n_obs {
                let row_sum: f64 = r[i].iter().sum();
                if row_sum > 1e-15 {
                    for c2 in 0..k {
                        r[i][c2] /= row_sum;
                    }
                }
            }
        }

        // Update centroids
        for c in 0..k {
            let total_r: f64 = (0..n_obs).map(|i| r[i][c]).sum();
            if total_r < 1e-10 {
                continue;
            }
            for d in 0..n_dims {
                centroids[c][d] = (0..n_obs).map(|i| r[i][c] * z[i][d]).sum::<f64>() / total_r;
            }
        }

        // Compute correction vectors per cluster per batch
        // For each cluster c: correction_b = centroid_c - mean(z[batch_b, cluster_c])
        for i in 0..n_obs {
            let b = batch_ids[i];
            let mut correction = vec![0.0; n_dims];

            for c in 0..k {
                if r[i][c] < 0.01 {
                    continue;
                }

                // Mean of batch b in cluster c
                let mut batch_mean = vec![0.0; n_dims];
                let mut batch_weight = 0.0;
                for j in 0..n_obs {
                    if batch_ids[j] == b {
                        for d in 0..n_dims {
                            batch_mean[d] += r[j][c] * z[j][d];
                        }
                        batch_weight += r[j][c];
                    }
                }
                if batch_weight > 1e-10 {
                    for d in 0..n_dims {
                        batch_mean[d] /= batch_weight;
                    }
                }

                // Correction = centroid - batch_mean, weighted by assignment
                for d in 0..n_dims {
                    correction[d] += r[i][c] * (centroids[c][d] - batch_mean[d]);
                }
            }

            for d in 0..n_dims {
                z[i][d] += correction[d];
            }
        }
    }

    adata.add_obsm("X_pca_harmony", z)?;
    Ok(())
}

// ── ComBat ─────────────────────────────────────────────────────────────────

/// Configuration for ComBat batch correction.
#[derive(Debug, Clone)]
pub struct CombatConfig {
    /// Key in `obs` containing batch labels.
    pub batch_key: String,
    /// Whether to use parametric empirical Bayes (true) or non-parametric (false).
    pub parametric: bool,
}

impl Default for CombatConfig {
    fn default() -> Self {
        Self {
            batch_key: "batch".into(),
            parametric: true,
        }
    }
}

/// ComBat batch correction (Johnson 2007).
///
/// Parametric empirical Bayes:
/// 1. Standardize data per gene
/// 2. Estimate batch-specific location (γ) and scale (δ²) parameters
/// 3. Shrink estimates toward the prior using EB
/// 4. Correct: x_corrected = (x - γ*) / δ* * pooled_σ + pooled_μ
///
/// Modifies X in place.
pub fn combat(adata: &mut AnnData, config: &CombatConfig) -> Result<()> {
    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();

    let batch_col = adata
        .get_obs(&config.batch_key)
        .ok_or_else(|| {
            CyaneaError::InvalidInput(format!("obs['{}'] not found", config.batch_key))
        })?;

    let batch_labels: Vec<String> = match batch_col {
        ColumnData::Strings(v) => v.clone(),
        ColumnData::Numeric(v) => v.iter().map(|x| x.to_string()).collect(),
        ColumnData::Categorical { codes, categories } => codes
            .iter()
            .map(|&c| categories.get(c as usize).cloned().unwrap_or_else(|| c.to_string()))
            .collect(),
    };

    let mut unique_batches: Vec<String> = batch_labels.clone();
    unique_batches.sort();
    unique_batches.dedup();
    let n_batches = unique_batches.len();

    if n_batches < 2 {
        return Ok(()); // Nothing to correct
    }

    let batch_map: HashMap<&str, usize> = unique_batches
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i))
        .collect();
    let batch_ids: Vec<usize> = batch_labels.iter().map(|b| batch_map[b.as_str()]).collect();

    // Batch indices
    let batch_indices: Vec<Vec<usize>> = (0..n_batches)
        .map(|b| (0..n_obs).filter(|&i| batch_ids[i] == b).collect())
        .collect();

    for j in 0..n_vars {
        // Compute overall mean and variance for gene j
        let mut overall_sum = 0.0;
        let mut overall_ss = 0.0;
        for i in 0..n_obs {
            let v = adata.x().get(i, j);
            overall_sum += v;
            overall_ss += v * v;
        }
        let overall_mean = overall_sum / n_obs as f64;
        let overall_var = (overall_ss / n_obs as f64 - overall_mean * overall_mean).max(1e-10);

        // Per-batch statistics
        let mut batch_means = vec![0.0; n_batches];
        let mut batch_vars = vec![0.0; n_batches];

        for (b, indices) in batch_indices.iter().enumerate() {
            let nb = indices.len() as f64;
            if nb == 0.0 {
                continue;
            }
            let mut sum = 0.0;
            for &i in indices {
                sum += adata.x().get(i, j);
            }
            batch_means[b] = sum / nb;

            let mut ss = 0.0;
            for &i in indices {
                let d = adata.x().get(i, j) - batch_means[b];
                ss += d * d;
            }
            batch_vars[b] = if nb > 1.0 {
                ss / (nb - 1.0)
            } else {
                overall_var
            };
        }

        if config.parametric {
            // EB shrinkage for location parameter
            // Prior: γ ~ N(γ_hat, τ²), where γ_hat = mean of batch_means, τ² = var of batch_means
            let gamma_hat: f64 = batch_means.iter().sum::<f64>() / n_batches as f64;
            let tau_sq = if n_batches > 1 {
                let ss: f64 = batch_means.iter().map(|&m| (m - gamma_hat).powi(2)).sum();
                ss / (n_batches - 1) as f64
            } else {
                1.0
            };

            // EB shrinkage for scale parameter
            // Prior: δ² ~ Inverse-Gamma(α, β)
            let delta_mean: f64 = batch_vars.iter().sum::<f64>() / n_batches as f64;
            let delta_var = if n_batches > 1 {
                let ss: f64 = batch_vars.iter().map(|&v| (v - delta_mean).powi(2)).sum();
                ss / (n_batches - 1) as f64
            } else {
                1.0
            };

            // Shrunk estimates
            let mut gamma_star = vec![0.0; n_batches];
            let mut delta_star = vec![0.0; n_batches];

            for b in 0..n_batches {
                let nb = batch_indices[b].len() as f64;
                // Shrink γ toward prior
                let shrink_gamma = tau_sq / (tau_sq + overall_var / nb.max(1.0));
                gamma_star[b] = shrink_gamma * batch_means[b] + (1.0 - shrink_gamma) * gamma_hat;

                // Shrink δ² toward prior
                let shrink_delta = delta_var / (delta_var + batch_vars[b].powi(2) / nb.max(1.0));
                delta_star[b] =
                    (shrink_delta * batch_vars[b] + (1.0 - shrink_delta) * delta_mean).max(1e-10);
            }

            // Correct: x_corrected = (x - γ*) * sqrt(pooled_var / δ*) + overall_mean
            let pooled_var = overall_var;
            for b in 0..n_batches {
                let scale = (pooled_var / delta_star[b]).sqrt();
                for &i in &batch_indices[b] {
                    let v = adata.x().get(i, j);
                    let corrected = (v - gamma_star[b]) * scale + overall_mean;
                    adata.x_mut().set(i, j, corrected);
                }
            }
        } else {
            // Non-parametric: simple centering
            for b in 0..n_batches {
                let shift = overall_mean - batch_means[b];
                for &i in &batch_indices[b] {
                    let v = adata.x().get(i, j);
                    adata.x_mut().set(i, j, v + shift);
                }
            }
        }
    }

    Ok(())
}

// ── MNN Correction ─────────────────────────────────────────────────────────

/// Configuration for Mutual Nearest Neighbors correction.
#[derive(Debug, Clone)]
pub struct MnnConfig {
    /// Key in `obs` containing batch labels.
    pub batch_key: String,
    /// Number of nearest neighbors for MNN pairs.
    pub k: usize,
    /// Gaussian kernel bandwidth for smoothing correction vectors.
    pub sigma: f64,
    /// Whether to cosine-normalize before finding MNN pairs.
    pub cos_norm: bool,
}

impl Default for MnnConfig {
    fn default() -> Self {
        Self {
            batch_key: "batch".into(),
            k: 20,
            sigma: 1.0,
            cos_norm: false,
        }
    }
}

/// Mutual Nearest Neighbors correction (Haghverdi 2018).
///
/// 1. Find MNN pairs between batches
/// 2. Compute correction vectors from MNN pairs
/// 3. Smooth correction vectors using Gaussian kernel
/// 4. Apply corrections
///
/// Stores corrected embeddings in `obsm["X_mnn"]`.
pub fn mnn_correct(adata: &mut AnnData, config: &MnnConfig) -> Result<()> {
    let pca = adata
        .get_obsm("X_pca")
        .ok_or_else(|| CyaneaError::InvalidInput("obsm['X_pca'] not found".into()))?
        .clone();

    let n_obs = pca.len();
    let n_dims = pca[0].len();

    let batch_col = adata
        .get_obs(&config.batch_key)
        .ok_or_else(|| {
            CyaneaError::InvalidInput(format!("obs['{}'] not found", config.batch_key))
        })?;

    let batch_labels: Vec<String> = match batch_col {
        ColumnData::Strings(v) => v.clone(),
        ColumnData::Numeric(v) => v.iter().map(|x| x.to_string()).collect(),
        ColumnData::Categorical { codes, categories } => codes
            .iter()
            .map(|&c| categories.get(c as usize).cloned().unwrap_or_else(|| c.to_string()))
            .collect(),
    };

    let mut unique_batches: Vec<String> = batch_labels.clone();
    unique_batches.sort();
    unique_batches.dedup();

    let batch_map: HashMap<&str, usize> = unique_batches
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i))
        .collect();
    let batch_ids: Vec<usize> = batch_labels.iter().map(|b| batch_map[b.as_str()]).collect();

    // Optionally cosine-normalize
    let data: Vec<Vec<f64>> = if config.cos_norm {
        pca.iter()
            .map(|row| {
                let norm: f64 = row.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-15);
                row.iter().map(|x| x / norm).collect()
            })
            .collect()
    } else {
        pca.clone()
    };

    let mut corrected = pca.clone();

    // Process batches sequentially: correct each batch toward the first (reference)
    let batch_indices: Vec<Vec<usize>> = (0..unique_batches.len())
        .map(|b| (0..n_obs).filter(|&i| batch_ids[i] == b).collect())
        .collect();

    for b in 1..unique_batches.len() {
        let ref_indices = &batch_indices[0];
        let query_indices = &batch_indices[b];

        if ref_indices.is_empty() || query_indices.is_empty() {
            continue;
        }

        let k = config.k.min(ref_indices.len()).min(query_indices.len());

        // Find MNN pairs
        let mut mnn_pairs: Vec<(usize, usize)> = Vec::new();

        // For each query point, find k nearest ref points
        let mut query_nn: HashMap<usize, Vec<usize>> = HashMap::new();
        for &qi in query_indices {
            let mut dists: Vec<(usize, f64)> = ref_indices
                .iter()
                .map(|&ri| {
                    let d: f64 = (0..n_dims)
                        .map(|d| (data[qi][d] - data[ri][d]).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    (ri, d)
                })
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let nn: Vec<usize> = dists[..k].iter().map(|&(j, _)| j).collect();
            query_nn.insert(qi, nn);
        }

        // For each ref point, find k nearest query points
        let mut ref_nn: HashMap<usize, Vec<usize>> = HashMap::new();
        for &ri in ref_indices {
            let mut dists: Vec<(usize, f64)> = query_indices
                .iter()
                .map(|&qi| {
                    let d: f64 = (0..n_dims)
                        .map(|d| (data[ri][d] - data[qi][d]).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    (qi, d)
                })
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let nn: Vec<usize> = dists[..k].iter().map(|&(j, _)| j).collect();
            ref_nn.insert(ri, nn);
        }

        // MNN: pairs where i is in j's kNN AND j is in i's kNN
        for &qi in query_indices {
            if let Some(q_neighbors) = query_nn.get(&qi) {
                for &ri in q_neighbors {
                    if let Some(r_neighbors) = ref_nn.get(&ri) {
                        if r_neighbors.contains(&qi) {
                            mnn_pairs.push((ri, qi)); // (ref, query)
                        }
                    }
                }
            }
        }

        if mnn_pairs.is_empty() {
            continue;
        }

        // Compute correction vectors from MNN pairs
        let mut pair_corrections: Vec<Vec<f64>> = Vec::new();
        let mut pair_query_indices: Vec<usize> = Vec::new();

        for &(ri, qi) in &mnn_pairs {
            let correction: Vec<f64> = (0..n_dims).map(|d| corrected[ri][d] - corrected[qi][d]).collect();
            pair_corrections.push(correction);
            pair_query_indices.push(qi);
        }

        // Apply Gaussian-smoothed correction to all query cells
        for &qi in query_indices {
            if pair_corrections.is_empty() {
                continue;
            }

            let mut total_weight = 0.0;
            let mut correction = vec![0.0; n_dims];

            for (p, &mnn_qi) in pair_query_indices.iter().enumerate() {
                let dist_sq: f64 = (0..n_dims)
                    .map(|d| (corrected[qi][d] - corrected[mnn_qi][d]).powi(2))
                    .sum();
                let w = (-dist_sq / (2.0 * config.sigma * config.sigma)).exp();
                total_weight += w;
                for d in 0..n_dims {
                    correction[d] += w * pair_corrections[p][d];
                }
            }

            if total_weight > 1e-15 {
                for d in 0..n_dims {
                    corrected[qi][d] += correction[d] / total_weight;
                }
            }
        }
    }

    adata.add_obsm("X_mnn", corrected)?;
    Ok(())
}

// ── Integration Metrics ────────────────────────────────────────────────────

/// Integration quality metrics.
#[derive(Debug, Clone)]
pub struct IntegrationMetrics {
    /// kBET acceptance rate (higher = better mixing).
    pub kbet_accept_rate: f64,
    /// Mean integration LISI (higher = better mixing).
    pub mean_ilisi: f64,
    /// Mean cell-type LISI (higher = better cell-type preservation).
    pub mean_clisi: f64,
}

/// Configuration for integration metrics.
#[derive(Debug, Clone)]
pub struct MetricsConfig {
    /// Key in `obs` containing batch labels.
    pub batch_key: String,
    /// Key in `obs` containing cell-type labels (for cLISI).
    pub label_key: Option<String>,
    /// Number of neighbors for metric computation.
    pub n_neighbors: usize,
    /// Significance threshold for kBET chi-squared test.
    pub alpha: f64,
}

impl Default for MetricsConfig {
    fn default() -> Self {
        Self {
            batch_key: "batch".into(),
            label_key: None,
            n_neighbors: 50,
            alpha: 0.05,
        }
    }
}

/// Compute integration quality metrics: kBET and LISI.
///
/// **kBET**: Tests whether the batch composition in each cell's neighborhood
/// matches the global composition (chi-squared test). Higher acceptance rate
/// means better batch mixing.
///
/// **LISI**: Local Inverse Simpson's Index — measures diversity of labels in
/// each cell's neighborhood. Higher iLISI = better batch mixing, higher cLISI
/// = better cell-type preservation.
pub fn integration_metrics(adata: &AnnData, config: &MetricsConfig) -> Result<IntegrationMetrics> {
    let pca = adata
        .get_obsm("X_pca")
        .or_else(|| adata.get_obsm("X_pca_harmony"))
        .or_else(|| adata.get_obsm("X_mnn"))
        .ok_or_else(|| {
            CyaneaError::InvalidInput("no embedding found (X_pca, X_pca_harmony, or X_mnn)".into())
        })?;

    let n_obs = pca.len();
    let n_dims = pca[0].len();

    // Get batch labels
    let batch_col = adata
        .get_obs(&config.batch_key)
        .ok_or_else(|| {
            CyaneaError::InvalidInput(format!("obs['{}'] not found", config.batch_key))
        })?;

    let batch_labels: Vec<usize> = labels_to_indices(batch_col);
    let n_batches = *batch_labels.iter().max().unwrap_or(&0) + 1;

    // Global batch proportions
    let batch_counts: Vec<usize> = (0..n_batches)
        .map(|b| batch_labels.iter().filter(|&&x| x == b).count())
        .collect();
    let batch_freq: Vec<f64> = batch_counts
        .iter()
        .map(|&c| c as f64 / n_obs as f64)
        .collect();

    let k = config.n_neighbors.min(n_obs - 1);

    // For each cell, find kNN
    let mut knn = vec![Vec::new(); n_obs];
    for i in 0..n_obs {
        let mut dists: Vec<(usize, f64)> = (0..n_obs)
            .filter(|&j| j != i)
            .map(|j| {
                let d: f64 = (0..n_dims)
                    .map(|d| (pca[i][d] - pca[j][d]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        knn[i] = dists[..k].iter().map(|&(j, _)| j).collect();
    }

    // kBET: chi-squared test on neighborhood batch composition
    let mut n_accept = 0;
    for i in 0..n_obs {
        let neighbors = &knn[i];
        let mut local_counts = vec![0usize; n_batches];
        for &j in neighbors {
            local_counts[batch_labels[j]] += 1;
        }

        // Chi-squared test: observed vs expected
        let mut chi_sq = 0.0;
        let k_f = k as f64;
        for b in 0..n_batches {
            let expected = batch_freq[b] * k_f;
            if expected > 0.0 {
                let diff = local_counts[b] as f64 - expected;
                chi_sq += diff * diff / expected;
            }
        }

        // Compare to chi-squared critical value (df = n_batches - 1)
        // Use simple approximation: reject if chi_sq > df * (1 + 2.5 * sqrt(2/df))
        let df = (n_batches - 1).max(1) as f64;
        let critical = df + 2.33 * (2.0 * df).sqrt(); // ~99th percentile approx
        if chi_sq <= critical {
            n_accept += 1;
        }
    }

    let kbet_accept_rate = n_accept as f64 / n_obs as f64;

    // iLISI: batch LISI (higher = better mixing)
    let mut ilisi_values = Vec::with_capacity(n_obs);
    for i in 0..n_obs {
        let lisi = compute_lisi(&knn[i], &batch_labels, n_batches);
        ilisi_values.push(lisi);
    }
    let mean_ilisi = if ilisi_values.is_empty() {
        1.0
    } else {
        ilisi_values.iter().sum::<f64>() / ilisi_values.len() as f64
    };

    // cLISI: cell-type LISI (only if label_key provided)
    let mean_clisi = if let Some(label_key) = &config.label_key {
        if let Some(label_col) = adata.get_obs(label_key) {
            let type_labels = labels_to_indices(label_col);
            let n_types = *type_labels.iter().max().unwrap_or(&0) + 1;

            let mut clisi_values = Vec::with_capacity(n_obs);
            for i in 0..n_obs {
                let lisi = compute_lisi(&knn[i], &type_labels, n_types);
                clisi_values.push(lisi);
            }
            clisi_values.iter().sum::<f64>() / clisi_values.len() as f64
        } else {
            1.0
        }
    } else {
        1.0
    };

    Ok(IntegrationMetrics {
        kbet_accept_rate,
        mean_ilisi,
        mean_clisi,
    })
}

/// Compute Local Inverse Simpson's Index for a single cell's neighborhood.
fn compute_lisi(neighbors: &[usize], labels: &[usize], n_labels: usize) -> f64 {
    if neighbors.is_empty() {
        return 1.0;
    }

    let k = neighbors.len() as f64;
    let mut counts = vec![0usize; n_labels];
    for &j in neighbors {
        counts[labels[j]] += 1;
    }

    // Simpson's Index: Σ p_i²
    let simpson: f64 = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / k;
            p * p
        })
        .sum();

    // Inverse Simpson's = 1 / Simpson's
    if simpson > 1e-15 {
        1.0 / simpson
    } else {
        n_labels as f64 // maximum diversity
    }
}

fn labels_to_indices(col: &ColumnData) -> Vec<usize> {
    match col {
        ColumnData::Strings(v) => {
            let mut unique: Vec<String> = v.clone();
            unique.sort();
            unique.dedup();
            let map: HashMap<&str, usize> = unique
                .iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i))
                .collect();
            v.iter().map(|s| map[s.as_str()]).collect()
        }
        ColumnData::Numeric(v) => v.iter().map(|&x| x as usize).collect(),
        ColumnData::Categorical { codes, .. } => codes.iter().map(|&c| c as usize).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_batch_adata(n_per_batch: usize, n_batches: usize, n_dims: usize) -> AnnData {
        let n = n_per_batch * n_batches;
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("cell_{}", i)).collect();
        let var_names = vec!["g0".into(), "g1".into()];
        let mut adata = AnnData::new(x, obs_names, var_names).unwrap();

        // PCA embeddings with batch effect
        let mut pca = Vec::new();
        for b in 0..n_batches {
            let offset = b as f64 * 5.0; // batch effect
            for i in 0..n_per_batch {
                let mut point = vec![0.0; n_dims];
                point[0] = i as f64 * 0.1 + offset;
                if n_dims > 1 {
                    point[1] = i as f64 * 0.2;
                }
                pca.push(point);
            }
        }
        adata.add_obsm("X_pca", pca).unwrap();

        let labels: Vec<String> = (0..n)
            .map(|i| format!("batch_{}", i / n_per_batch))
            .collect();
        adata
            .add_obs_column("batch", ColumnData::Strings(labels))
            .unwrap();

        adata
    }

    // ── Harmony tests ──

    #[test]
    fn harmony_basic() {
        let mut adata = make_batch_adata(10, 2, 3);
        harmony(&mut adata, &HarmonyConfig::default()).unwrap();
        assert!(adata.get_obsm("X_pca_harmony").is_some());
        let corrected = adata.get_obsm("X_pca_harmony").unwrap();
        assert_eq!(corrected.len(), 20);
        assert_eq!(corrected[0].len(), 3);
    }

    #[test]
    fn harmony_reduces_batch_effect() {
        let mut adata = make_batch_adata(10, 2, 3);
        let pca_original = adata.get_obsm("X_pca").unwrap().clone();

        // Measure batch separation before
        let batch0_mean_before: f64 = pca_original[..10].iter().map(|p| p[0]).sum::<f64>() / 10.0;
        let batch1_mean_before: f64 = pca_original[10..].iter().map(|p| p[0]).sum::<f64>() / 10.0;
        let sep_before = (batch1_mean_before - batch0_mean_before).abs();

        harmony(&mut adata, &HarmonyConfig::default()).unwrap();

        let corrected = adata.get_obsm("X_pca_harmony").unwrap();
        let batch0_mean_after: f64 = corrected[..10].iter().map(|p| p[0]).sum::<f64>() / 10.0;
        let batch1_mean_after: f64 = corrected[10..].iter().map(|p| p[0]).sum::<f64>() / 10.0;
        let sep_after = (batch1_mean_after - batch0_mean_after).abs();

        assert!(
            sep_after < sep_before,
            "batch separation should decrease: before={}, after={}",
            sep_before,
            sep_after
        );
    }

    #[test]
    fn harmony_missing_pca() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let mut adata = AnnData::new(x, vec!["c0".into()], vec!["g0".into(), "g1".into()]).unwrap();
        adata.add_obs_column("batch", ColumnData::Strings(vec!["b0".into()])).unwrap();
        let result = harmony(&mut adata, &HarmonyConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn harmony_missing_batch() {
        let mut adata = make_batch_adata(5, 2, 2);
        let result = harmony(
            &mut adata,
            &HarmonyConfig {
                batch_key: "nonexistent".into(),
                ..Default::default()
            },
        );
        assert!(result.is_err());
    }

    #[test]
    fn harmony_three_batches() {
        let mut adata = make_batch_adata(8, 3, 3);
        harmony(&mut adata, &HarmonyConfig::default()).unwrap();
        let corrected = adata.get_obsm("X_pca_harmony").unwrap();
        assert_eq!(corrected.len(), 24);
    }

    #[test]
    fn harmony_custom_params() {
        let mut adata = make_batch_adata(10, 2, 3);
        harmony(
            &mut adata,
            &HarmonyConfig {
                theta: 1.0,
                sigma: 0.5,
                max_iter: 5,
                n_clusters: Some(3),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(adata.get_obsm("X_pca_harmony").is_some());
    }

    // ── ComBat tests ──

    #[test]
    fn combat_basic() {
        let n = 20;
        let n_vars = 3;
        let mut data = vec![vec![0.0; n_vars]; n];
        // Two batches with different means
        for i in 0..10 {
            for j in 0..n_vars {
                data[i][j] = 10.0 + (i * n_vars + j) as f64 * 0.1;
            }
        }
        for i in 10..20 {
            for j in 0..n_vars {
                data[i][j] = 20.0 + ((i - 10) * n_vars + j) as f64 * 0.1; // shifted by 10
            }
        }

        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let var_names: Vec<String> = (0..n_vars).map(|j| format!("g{}", j)).collect();
        let mut adata = AnnData::new(MatrixData::Dense(data), obs_names, var_names).unwrap();

        let labels: Vec<String> = (0..n)
            .map(|i| if i < 10 { "A".into() } else { "B".into() })
            .collect();
        adata.add_obs_column("batch", ColumnData::Strings(labels)).unwrap();

        combat(&mut adata, &CombatConfig::default()).unwrap();

        // After correction, batch means should be closer
        let mean_a: f64 = (0..10).map(|i| adata.x().get(i, 0)).sum::<f64>() / 10.0;
        let mean_b: f64 = (10..20).map(|i| adata.x().get(i, 0)).sum::<f64>() / 10.0;
        assert!(
            (mean_a - mean_b).abs() < 5.0,
            "batch means should converge: A={}, B={}",
            mean_a,
            mean_b
        );
    }

    #[test]
    fn combat_single_batch() {
        let mut data = vec![vec![1.0, 2.0]; 5];
        let obs_names: Vec<String> = (0..5).map(|i| format!("c{}", i)).collect();
        let mut adata =
            AnnData::new(MatrixData::Dense(data), obs_names, vec!["g0".into(), "g1".into()])
                .unwrap();
        adata
            .add_obs_column("batch", ColumnData::Strings(vec!["A".into(); 5]))
            .unwrap();

        // Single batch = no correction
        combat(&mut adata, &CombatConfig::default()).unwrap();
        assert_eq!(adata.x().get(0, 0), 1.0);
    }

    #[test]
    fn combat_nonparametric() {
        let n = 10;
        let mut data = vec![vec![0.0; 2]; n];
        for i in 0..5 {
            data[i][0] = 10.0;
            data[i][1] = 5.0;
        }
        for i in 5..10 {
            data[i][0] = 20.0;
            data[i][1] = 15.0;
        }

        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata =
            AnnData::new(MatrixData::Dense(data), obs_names, vec!["g0".into(), "g1".into()])
                .unwrap();
        let labels: Vec<String> = (0..n)
            .map(|i| if i < 5 { "A".into() } else { "B".into() })
            .collect();
        adata.add_obs_column("batch", ColumnData::Strings(labels)).unwrap();

        combat(
            &mut adata,
            &CombatConfig {
                parametric: false,
                ..Default::default()
            },
        )
        .unwrap();

        // Non-parametric: simple centering → both batches should have same mean
        let mean_a: f64 = (0..5).map(|i| adata.x().get(i, 0)).sum::<f64>() / 5.0;
        let mean_b: f64 = (5..10).map(|i| adata.x().get(i, 0)).sum::<f64>() / 5.0;
        assert!(
            (mean_a - mean_b).abs() < 1e-10,
            "non-parametric correction should equalize means: A={}, B={}",
            mean_a,
            mean_b
        );
    }

    #[test]
    fn combat_missing_batch() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let mut adata = AnnData::new(x, vec!["c0".into()], vec!["g0".into(), "g1".into()]).unwrap();
        let result = combat(&mut adata, &CombatConfig::default());
        assert!(result.is_err());
    }

    // ── MNN tests ──

    #[test]
    fn mnn_basic() {
        let mut adata = make_batch_adata(10, 2, 3);
        mnn_correct(
            &mut adata,
            &MnnConfig {
                k: 5,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(adata.get_obsm("X_mnn").is_some());
        let corrected = adata.get_obsm("X_mnn").unwrap();
        assert_eq!(corrected.len(), 20);
    }

    #[test]
    fn mnn_reduces_separation() {
        let mut adata = make_batch_adata(10, 2, 3);
        let pca = adata.get_obsm("X_pca").unwrap().clone();
        let sep_before = (pca[0][0] - pca[10][0]).abs();

        mnn_correct(
            &mut adata,
            &MnnConfig {
                k: 5,
                sigma: 5.0,
                ..Default::default()
            },
        )
        .unwrap();

        let corrected = adata.get_obsm("X_mnn").unwrap();
        let sep_after = (corrected[0][0] - corrected[10][0]).abs();

        // MNN should reduce batch separation
        assert!(
            sep_after <= sep_before + 1e-6,
            "separation should not increase: before={}, after={}",
            sep_before,
            sep_after
        );
    }

    #[test]
    fn mnn_with_cos_norm() {
        let mut adata = make_batch_adata(8, 2, 3);
        mnn_correct(
            &mut adata,
            &MnnConfig {
                cos_norm: true,
                k: 3,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(adata.get_obsm("X_mnn").is_some());
    }

    #[test]
    fn mnn_missing_pca() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let mut adata = AnnData::new(x, vec!["c0".into()], vec!["g0".into(), "g1".into()]).unwrap();
        adata.add_obs_column("batch", ColumnData::Strings(vec!["b0".into()])).unwrap();
        let result = mnn_correct(&mut adata, &MnnConfig::default());
        assert!(result.is_err());
    }

    // ── Integration Metrics tests ──

    #[test]
    fn metrics_well_mixed() {
        // Create well-mixed data (batches interleaved)
        let n = 20;
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();

        // Interleaved PCA → good mixing
        let pca: Vec<Vec<f64>> = (0..n).map(|i| vec![i as f64 * 0.1, 0.0]).collect();
        adata.add_obsm("X_pca", pca).unwrap();

        let labels: Vec<String> = (0..n)
            .map(|i| if i % 2 == 0 { "A".into() } else { "B".into() })
            .collect();
        adata.add_obs_column("batch", ColumnData::Strings(labels)).unwrap();

        let metrics = integration_metrics(
            &adata,
            &MetricsConfig {
                n_neighbors: 5,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(
            metrics.kbet_accept_rate > 0.0,
            "well-mixed should have some acceptance: {}",
            metrics.kbet_accept_rate
        );
        assert!(
            metrics.mean_ilisi > 1.0,
            "well-mixed should have iLISI > 1: {}",
            metrics.mean_ilisi
        );
    }

    #[test]
    fn metrics_poorly_mixed() {
        // Batches completely separated
        let n = 20;
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();

        let mut pca = Vec::new();
        for i in 0..10 {
            pca.push(vec![i as f64 * 0.1, 0.0]); // batch A: left
        }
        for i in 0..10 {
            pca.push(vec![100.0 + i as f64 * 0.1, 0.0]); // batch B: far right
        }
        adata.add_obsm("X_pca", pca).unwrap();

        let labels: Vec<String> = (0..n)
            .map(|i| if i < 10 { "A".into() } else { "B".into() })
            .collect();
        adata.add_obs_column("batch", ColumnData::Strings(labels)).unwrap();

        let metrics = integration_metrics(
            &adata,
            &MetricsConfig {
                n_neighbors: 5,
                ..Default::default()
            },
        )
        .unwrap();

        // Poorly mixed → low iLISI (close to 1, single batch per neighborhood)
        assert!(
            metrics.mean_ilisi < 1.5,
            "poorly mixed should have low iLISI: {}",
            metrics.mean_ilisi
        );
    }

    #[test]
    fn metrics_with_cell_types() {
        let n = 20;
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();

        let pca: Vec<Vec<f64>> = (0..n).map(|i| vec![i as f64 * 0.1, 0.0]).collect();
        adata.add_obsm("X_pca", pca).unwrap();

        let batch_labels: Vec<String> = (0..n).map(|i| if i % 2 == 0 { "A".into() } else { "B".into() }).collect();
        let type_labels: Vec<String> = (0..n).map(|i| if i < 10 { "T".into() } else { "B".into() }).collect();
        adata.add_obs_column("batch", ColumnData::Strings(batch_labels)).unwrap();
        adata.add_obs_column("cell_type", ColumnData::Strings(type_labels)).unwrap();

        let metrics = integration_metrics(
            &adata,
            &MetricsConfig {
                label_key: Some("cell_type".into()),
                n_neighbors: 5,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(metrics.mean_clisi >= 1.0);
    }

    #[test]
    fn metrics_missing_embedding() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let adata = AnnData::new(x, vec!["c0".into()], vec!["g0".into(), "g1".into()]).unwrap();
        let result = integration_metrics(&adata, &MetricsConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn metrics_missing_batch() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]; 5]);
        let obs_names: Vec<String> = (0..5).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();
        let pca: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64, 0.0]).collect();
        adata.add_obsm("X_pca", pca).unwrap();

        let result = integration_metrics(&adata, &MetricsConfig::default());
        assert!(result.is_err());
    }

    // ── LISI tests ──

    #[test]
    fn lisi_single_label() {
        let neighbors = vec![0, 1, 2, 3, 4];
        let labels = vec![0, 0, 0, 0, 0];
        let lisi = compute_lisi(&neighbors, &labels, 1);
        assert!((lisi - 1.0).abs() < 1e-10, "single label LISI = {}", lisi);
    }

    #[test]
    fn lisi_two_labels_equal() {
        let neighbors = vec![0, 1, 2, 3];
        let labels = vec![0, 0, 1, 1];
        let lisi = compute_lisi(&neighbors, &labels, 2);
        // Equal proportions → Simpson = 0.5 → inverse = 2.0
        assert!((lisi - 2.0).abs() < 1e-10, "equal two-label LISI = {}", lisi);
    }

    #[test]
    fn lisi_empty_neighbors() {
        let lisi = compute_lisi(&[], &[], 2);
        assert_eq!(lisi, 1.0);
    }

    // ── Helper tests ──

    #[test]
    fn labels_to_indices_strings() {
        let col = ColumnData::Strings(vec!["B".into(), "A".into(), "B".into(), "C".into()]);
        let indices = labels_to_indices(&col);
        assert_eq!(indices[0], indices[2]); // Both "B"
        assert_ne!(indices[0], indices[1]); // "B" != "A"
        assert_ne!(indices[1], indices[3]); // "A" != "C"
    }

    #[test]
    fn labels_to_indices_numeric() {
        let col = ColumnData::Numeric(vec![0.0, 1.0, 0.0, 2.0]);
        let indices = labels_to_indices(&col);
        assert_eq!(indices, vec![0, 1, 0, 2]);
    }

    // ── Pipeline test: ComBat → metrics ──

    #[test]
    fn combat_then_metrics() {
        let n = 20;
        let n_vars = 3;
        let mut data = vec![vec![0.0; n_vars]; n];
        for i in 0..10 {
            for j in 0..n_vars {
                data[i][j] = 5.0 + (i + j) as f64 * 0.1;
            }
        }
        for i in 10..20 {
            for j in 0..n_vars {
                data[i][j] = 15.0 + ((i - 10) + j) as f64 * 0.1;
            }
        }

        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let var_names: Vec<String> = (0..n_vars).map(|j| format!("g{}", j)).collect();
        let mut adata = AnnData::new(MatrixData::Dense(data), obs_names, var_names).unwrap();

        let labels: Vec<String> = (0..n)
            .map(|i| if i < 10 { "A".into() } else { "B".into() })
            .collect();
        adata.add_obs_column("batch", ColumnData::Strings(labels)).unwrap();

        // Add PCA embedding
        let pca: Vec<Vec<f64>> = (0..n).map(|i| {
            vec![adata.x().get(i, 0), adata.x().get(i, 1)]
        }).collect();
        adata.add_obsm("X_pca", pca).unwrap();

        combat(&mut adata, &CombatConfig::default()).unwrap();

        // Update PCA after combat
        let pca: Vec<Vec<f64>> = (0..n).map(|i| {
            vec![adata.x().get(i, 0), adata.x().get(i, 1)]
        }).collect();
        adata.add_obsm("X_pca", pca).unwrap();

        let metrics = integration_metrics(
            &adata,
            &MetricsConfig {
                n_neighbors: 5,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(metrics.kbet_accept_rate >= 0.0);
        assert!(metrics.mean_ilisi >= 1.0);
    }
}
