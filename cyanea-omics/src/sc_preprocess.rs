//! Single-cell preprocessing: normalization, HVG selection, doublet detection, gene scoring.
//!
//! All functions operate on [`AnnData`], storing results in obs/var/layers
//! following scanpy conventions.

use cyanea_core::{CyaneaError, Result};

use crate::single_cell::{AnnData, ColumnData, MatrixData};
use crate::sparse::SparseMatrix;

// ── HVG Selection ──────────────────────────────────────────────────────────

/// Method for highly variable gene selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HvgMethod {
    /// Seurat v3 variance-stabilizing transformation.
    SeuratV3,
    /// Cell Ranger mean/dispersion binning method.
    CellRanger,
}

/// Configuration for highly variable gene selection.
#[derive(Debug, Clone)]
pub struct HvgConfig {
    /// Number of top genes to select.
    pub n_top_genes: usize,
    /// Selection method.
    pub method: HvgMethod,
    /// Minimum mean expression for CellRanger method.
    pub min_mean: f64,
    /// Maximum mean expression for CellRanger method.
    pub max_mean: f64,
    /// Minimum normalized dispersion for CellRanger method.
    pub min_disp: f64,
    /// Number of bins for CellRanger dispersion normalization.
    pub n_bins: usize,
}

impl Default for HvgConfig {
    fn default() -> Self {
        Self {
            n_top_genes: 2000,
            method: HvgMethod::SeuratV3,
            min_mean: 0.0125,
            max_mean: 3.0,
            min_disp: 0.5,
            n_bins: 20,
        }
    }
}

/// Select highly variable genes and annotate `var["highly_variable"]` and
/// `var["dispersions_norm"]`.
///
/// **Seurat v3 (VST):** For each gene, fits a loess-like curve of variance vs mean,
/// then selects genes with highest standardized variance.
///
/// **CellRanger:** Bins genes by mean expression, normalizes dispersion within bins,
/// then filters by mean/dispersion thresholds.
pub fn highly_variable_genes(adata: &mut AnnData, config: &HvgConfig) -> Result<()> {
    let n_vars = adata.n_vars();
    let n_obs = adata.n_obs();
    if n_obs < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 observations for HVG selection".into(),
        ));
    }

    // Compute per-gene mean and variance
    let mut means = vec![0.0; n_vars];
    let mut variances = vec![0.0; n_vars];

    for j in 0..n_vars {
        let mut sum = 0.0;
        for i in 0..n_obs {
            sum += adata.x().get(i, j);
        }
        means[j] = sum / n_obs as f64;
    }

    for j in 0..n_vars {
        let mut ss = 0.0;
        for i in 0..n_obs {
            let d = adata.x().get(i, j) - means[j];
            ss += d * d;
        }
        variances[j] = ss / (n_obs - 1) as f64;
    }

    let dispersions_norm = match config.method {
        HvgMethod::SeuratV3 => hvg_seurat_v3(&means, &variances, n_vars),
        HvgMethod::CellRanger => {
            hvg_cell_ranger(&means, &variances, n_vars, config.n_bins)
        }
    };

    // Select top genes
    let n_select = config.n_top_genes.min(n_vars);
    let mut ranked: Vec<(usize, f64)> = dispersions_norm.iter().copied().enumerate().collect();
    ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut highly_variable = vec![0.0; n_vars];
    for &(idx, _) in ranked.iter().take(n_select) {
        highly_variable[idx] = 1.0;
    }

    // For CellRanger, also apply mean/dispersion filters
    if config.method == HvgMethod::CellRanger {
        for j in 0..n_vars {
            if means[j] < config.min_mean
                || means[j] > config.max_mean
                || dispersions_norm[j] < config.min_disp
            {
                highly_variable[j] = 0.0;
            }
        }
    }

    adata.add_var_column("highly_variable", ColumnData::Numeric(highly_variable))?;
    adata.add_var_column("dispersions_norm", ColumnData::Numeric(dispersions_norm))?;
    adata.add_var_column("means", ColumnData::Numeric(means))?;
    adata.add_var_column("variances", ColumnData::Numeric(variances))?;

    Ok(())
}

/// Seurat v3 VST: standardized variance after fitting mean-variance relationship.
fn hvg_seurat_v3(means: &[f64], variances: &[f64], n_vars: usize) -> Vec<f64> {
    // Clip means to avoid log(0)
    let clipped_means: Vec<f64> = means.iter().map(|&m| m.max(1e-12)).collect();
    let log_means: Vec<f64> = clipped_means.iter().map(|m| m.ln()).collect();

    // Sort genes by mean for local fitting
    let mut order: Vec<usize> = (0..n_vars).collect();
    order.sort_by(|&a, &b| log_means[a].partial_cmp(&log_means[b]).unwrap_or(std::cmp::Ordering::Equal));

    // Expected variance: use local window to estimate mean-variance trend
    let window = (n_vars / 10).max(10).min(n_vars);
    let mut expected_var = vec![0.0; n_vars];

    for (rank, &idx) in order.iter().enumerate() {
        let start = rank.saturating_sub(window / 2);
        let end = (rank + window / 2 + 1).min(n_vars);
        let mut sum = 0.0;
        let mut count = 0;
        for &j in &order[start..end] {
            sum += variances[j];
            count += 1;
        }
        expected_var[idx] = if count > 0 { sum / count as f64 } else { variances[idx] };
    }

    // Standardized variance: variance / expected_variance
    (0..n_vars)
        .map(|j| {
            let exp = expected_var[j].max(1e-12);
            variances[j] / exp
        })
        .collect()
}

/// CellRanger: bin genes by mean, normalize dispersion within bins.
fn hvg_cell_ranger(
    means: &[f64],
    variances: &[f64],
    n_vars: usize,
    n_bins: usize,
) -> Vec<f64> {
    // Dispersion = variance / mean (coefficient of variation squared)
    let dispersions: Vec<f64> = (0..n_vars)
        .map(|j| {
            let m = means[j].max(1e-12);
            variances[j] / m
        })
        .collect();

    // Bin genes by mean expression
    let mean_max = means.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean_min = means.iter().cloned().fold(f64::INFINITY, f64::min);
    let range = (mean_max - mean_min).max(1e-12);
    let bin_width = range / n_bins as f64;

    // Assign each gene to a bin
    let bins: Vec<usize> = means
        .iter()
        .map(|&m| {
            let b = ((m - mean_min) / bin_width) as usize;
            b.min(n_bins - 1)
        })
        .collect();

    // Compute median and MAD of log-dispersions within each bin
    let log_disps: Vec<f64> = dispersions.iter().map(|&d| d.max(1e-12).ln()).collect();

    let mut bin_medians = vec![0.0; n_bins];
    let mut bin_mads = vec![1.0; n_bins];

    for b in 0..n_bins {
        let vals: Vec<f64> = (0..n_vars)
            .filter(|&j| bins[j] == b)
            .map(|j| log_disps[j])
            .collect();
        if vals.is_empty() {
            continue;
        }
        let med = simple_median(&vals);
        bin_medians[b] = med;
        let deviations: Vec<f64> = vals.iter().map(|&v| (v - med).abs()).collect();
        bin_mads[b] = simple_median(&deviations).max(1e-12);
    }

    // Normalized dispersion = (log_disp - bin_median) / bin_mad
    (0..n_vars)
        .map(|j| {
            let b = bins[j];
            (log_disps[j] - bin_medians[b]) / bin_mads[b]
        })
        .collect()
}

fn simple_median(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

// ── Normalization ──────────────────────────────────────────────────────────

/// Configuration for total-count normalization.
#[derive(Debug, Clone)]
pub struct NormalizeConfig {
    /// Target sum per cell after normalization.
    pub target_sum: f64,
    /// Whether to apply log1p transformation after scaling.
    pub log_transform: bool,
    /// Whether to save raw counts in `layers["counts"]`.
    pub save_raw: bool,
}

impl Default for NormalizeConfig {
    fn default() -> Self {
        Self {
            target_sum: 1e4,
            log_transform: true,
            save_raw: true,
        }
    }
}

/// Normalize per-cell total counts to `target_sum`, optionally log-transforming.
///
/// Follows the scanpy `normalize_total` + `log1p` pattern. If `save_raw` is true,
/// original counts are saved to `layers["counts"]`.
pub fn normalize_total(adata: &mut AnnData, config: &NormalizeConfig) -> Result<()> {
    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();

    // Save raw if requested
    if config.save_raw {
        let raw = adata.x().clone();
        adata.add_layer("counts", raw)?;
    }

    // Compute per-cell size factors
    let row_sums = adata.x().row_sums();
    let factors: Vec<f64> = row_sums
        .iter()
        .map(|&s| {
            if s > 0.0 {
                config.target_sum / s
            } else {
                0.0
            }
        })
        .collect();

    // Apply normalization
    match adata.x_mut() {
        MatrixData::Dense(rows) => {
            for i in 0..n_obs {
                for j in 0..n_vars {
                    rows[i][j] *= factors[i];
                    if config.log_transform {
                        rows[i][j] = (rows[i][j] + 1.0).ln();
                    }
                }
            }
        }
        MatrixData::Sparse(s) => {
            s.scale_rows(&factors);
            if config.log_transform {
                s.map_values(|v| (v + 1.0).ln());
            }
        }
    }

    Ok(())
}

// ── Regress Out ────────────────────────────────────────────────────────────

/// Regress out confounding variables from the expression matrix.
///
/// For each gene, fits OLS on the specified `obs` keys and subtracts the
/// predicted values. Expects numeric columns in `obs`.
pub fn regress_out(adata: &mut AnnData, keys: &[&str]) -> Result<()> {
    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();

    if keys.is_empty() {
        return Ok(());
    }

    // Build design matrix: intercept + covariates
    let n_covariates = keys.len();
    let n_cols = n_covariates + 1; // +1 for intercept
    let mut design = vec![vec![0.0; n_cols]; n_obs];
    for i in 0..n_obs {
        design[i][0] = 1.0; // intercept
    }

    for (k, &key) in keys.iter().enumerate() {
        let col = adata
            .get_obs(key)
            .ok_or_else(|| CyaneaError::InvalidInput(format!("obs key '{}' not found", key)))?;
        let values = col
            .as_numeric()
            .ok_or_else(|| {
                CyaneaError::InvalidInput(format!("obs '{}' must be numeric for regression", key))
            })?;
        for i in 0..n_obs {
            design[i][k + 1] = values[i];
        }
    }

    // Solve OLS for each gene: beta = (X^T X)^{-1} X^T y
    // X^T X is small (n_cols × n_cols), so direct solve is fine
    let xtx = mat_mul_ata(&design, n_obs, n_cols);
    let xtx_inv = match invert_small_matrix(&xtx, n_cols) {
        Some(inv) => inv,
        None => return Err(CyaneaError::InvalidInput("singular design matrix".into())),
    };

    for j in 0..n_vars {
        // y = gene expression column
        let y: Vec<f64> = (0..n_obs).map(|i| adata.x().get(i, j)).collect();

        // X^T y
        let xty: Vec<f64> = (0..n_cols)
            .map(|k| (0..n_obs).map(|i| design[i][k] * y[i]).sum())
            .collect();

        // beta = (X^T X)^{-1} X^T y
        let beta: Vec<f64> = (0..n_cols)
            .map(|k| (0..n_cols).map(|l| xtx_inv[k * n_cols + l] * xty[l]).sum())
            .collect();

        // residual = y - X * beta
        for i in 0..n_obs {
            let predicted: f64 = (0..n_cols).map(|k| design[i][k] * beta[k]).sum();
            let residual = y[i] - predicted;
            adata.x_mut().set(i, j, residual);
        }
    }

    Ok(())
}

/// Compute X^T X for a design matrix stored as Vec<Vec<f64>>.
fn mat_mul_ata(x: &[Vec<f64>], n_rows: usize, n_cols: usize) -> Vec<f64> {
    let mut result = vec![0.0; n_cols * n_cols];
    for i in 0..n_rows {
        for j in 0..n_cols {
            for k in j..n_cols {
                let v = x[i][j] * x[i][k];
                result[j * n_cols + k] += v;
                if j != k {
                    result[k * n_cols + j] += v;
                }
            }
        }
    }
    result
}

/// Invert a small symmetric positive-definite matrix via Gauss-Jordan.
fn invert_small_matrix(m: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            aug[i * 2 * n + j] = m[i * n + j];
        }
        aug[i * 2 * n + n + i] = 1.0;
    }

    for col in 0..n {
        // Find pivot
        let mut max_val = aug[col * 2 * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let val = aug[row * 2 * n + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_val < 1e-15 {
            return None;
        }

        // Swap rows
        if max_row != col {
            for j in 0..2 * n {
                let tmp = aug[col * 2 * n + j];
                aug[col * 2 * n + j] = aug[max_row * 2 * n + j];
                aug[max_row * 2 * n + j] = tmp;
            }
        }

        // Scale pivot row
        let pivot = aug[col * 2 * n + col];
        for j in 0..2 * n {
            aug[col * 2 * n + j] /= pivot;
        }

        // Eliminate
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row * 2 * n + col];
            for j in 0..2 * n {
                aug[row * 2 * n + j] -= factor * aug[col * 2 * n + j];
            }
        }
    }

    let mut result = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            result[i * n + j] = aug[i * 2 * n + n + j];
        }
    }
    Some(result)
}

// ── Scrublet Doublet Detection ─────────────────────────────────────────────

/// Configuration for Scrublet doublet detection.
#[derive(Debug, Clone)]
pub struct ScrubletConfig {
    /// Expected doublet rate (default 0.06).
    pub expected_doublet_rate: f64,
    /// Ratio of simulated doublets to observed cells.
    pub sim_doublet_ratio: f64,
    /// Number of PCs to use.
    pub n_pcs: usize,
    /// Number of neighbors for kNN scoring.
    pub k_neighbors: usize,
    /// Random seed.
    pub seed: u64,
}

impl Default for ScrubletConfig {
    fn default() -> Self {
        Self {
            expected_doublet_rate: 0.06,
            sim_doublet_ratio: 2.0,
            n_pcs: 30,
            k_neighbors: 20,
            seed: 42,
        }
    }
}

/// Detect doublets using the Scrublet algorithm (Wolock 2019).
///
/// 1. Simulate doublets by averaging random pairs of cells
/// 2. Run PCA on combined (observed + simulated)
/// 3. Score each cell by the fraction of kNN that are simulated doublets
/// 4. Store scores in `obs["doublet_score"]` and calls in `obs["predicted_doublet"]`
pub fn scrublet_doublets(adata: &mut AnnData, config: &ScrubletConfig) -> Result<()> {
    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();

    if n_obs < 3 {
        return Err(CyaneaError::InvalidInput(
            "need at least 3 cells for doublet detection".into(),
        ));
    }

    let n_sim = (n_obs as f64 * config.sim_doublet_ratio).ceil() as usize;

    // Build observed matrix as flat row-major
    let obs_flat = adata.x().to_flat_row_major();

    // Simulate doublets by averaging random pairs
    let mut sim_flat = vec![0.0; n_sim * n_vars];
    let mut rng_state = config.seed;
    for s in 0..n_sim {
        let i = lcg_next(&mut rng_state) as usize % n_obs;
        let j = lcg_next(&mut rng_state) as usize % n_obs;
        for v in 0..n_vars {
            sim_flat[s * n_vars + v] =
                (obs_flat[i * n_vars + v] + obs_flat[j * n_vars + v]) / 2.0;
        }
    }

    // Combine observed + simulated
    let n_total = n_obs + n_sim;
    let mut combined = Vec::with_capacity(n_total * n_vars);
    combined.extend_from_slice(&obs_flat);
    combined.extend_from_slice(&sim_flat);

    // PCA
    let n_pcs = config.n_pcs.min(n_vars).min(n_total);
    let pca_config = cyanea_ml::reduction::PcaConfig {
        n_components: n_pcs,
        max_iter: 1000,
        tolerance: 1e-10,
    };
    let pca_result = cyanea_ml::reduction::pca(&combined, n_vars, &pca_config)?;

    // Score each observed cell: fraction of kNN that are simulated
    let mut scores = vec![0.0; n_obs];
    for i in 0..n_obs {
        let query = &pca_result.transformed[i * n_pcs..(i + 1) * n_pcs];
        let mut dists: Vec<(usize, f64)> = (0..n_total)
            .filter(|&j| j != i)
            .map(|j| {
                let other = &pca_result.transformed[j * n_pcs..(j + 1) * n_pcs];
                let d: f64 = query
                    .iter()
                    .zip(other.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let k = config.k_neighbors.min(dists.len());
        let sim_count = dists[..k]
            .iter()
            .filter(|&&(idx, _)| idx >= n_obs)
            .count();
        scores[i] = sim_count as f64 / k as f64;
    }

    // Bimodal threshold detection: use simple Otsu-like method
    let threshold = otsu_threshold(&scores);

    let predicted: Vec<f64> = scores.iter().map(|&s| if s > threshold { 1.0 } else { 0.0 }).collect();

    adata.add_obs_column("doublet_score", ColumnData::Numeric(scores))?;
    adata.add_obs_column("predicted_doublet", ColumnData::Numeric(predicted))?;

    Ok(())
}

/// Simple Otsu threshold for bimodal score distribution.
fn otsu_threshold(scores: &[f64]) -> f64 {
    if scores.is_empty() {
        return 0.5;
    }
    let n = scores.len() as f64;
    let mut sorted = scores.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let total_sum: f64 = sorted.iter().sum();
    let mut best_threshold = sorted[sorted.len() / 2];
    let mut best_variance = f64::NEG_INFINITY;

    // Test each unique value as threshold
    let mut sum_bg = 0.0;
    let mut count_bg = 0.0;

    for i in 0..sorted.len() - 1 {
        sum_bg += sorted[i];
        count_bg += 1.0;
        let count_fg = n - count_bg;
        if count_fg == 0.0 {
            break;
        }
        let mean_bg = sum_bg / count_bg;
        let mean_fg = (total_sum - sum_bg) / count_fg;
        let between_var = count_bg * count_fg * (mean_bg - mean_fg).powi(2);
        if between_var > best_variance {
            best_variance = between_var;
            best_threshold = (sorted[i] + sorted[i + 1]) / 2.0;
        }
    }

    best_threshold
}

/// Simple LCG random number generator.
fn lcg_next(state: &mut u64) -> u64 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *state >> 33
}

// ── Gene Scoring ───────────────────────────────────────────────────────────

/// Score cells by expression of a gene set (like scanpy's `score_genes`).
///
/// For each cell, the score = mean(signature genes) - mean(reference genes).
/// Reference genes are expression-matched from the remaining genes.
///
/// Stores result in `obs[score_name]`.
pub fn score_genes(
    adata: &mut AnnData,
    gene_list: &[usize],
    n_reference: usize,
    score_name: &str,
) -> Result<()> {
    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();

    for &g in gene_list {
        if g >= n_vars {
            return Err(CyaneaError::InvalidInput(format!(
                "gene index {} out of bounds (n_vars={})",
                g, n_vars
            )));
        }
    }

    // Compute mean expression of each gene across cells (for matching)
    let gene_means = adata.x().column_means();

    // Sort all genes by mean expression
    let mut all_genes: Vec<usize> = (0..n_vars).collect();
    all_genes.sort_by(|&a, &b| {
        gene_means[a]
            .partial_cmp(&gene_means[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // For each signature gene, pick `n_reference` expression-matched reference genes
    let gene_set: std::collections::HashSet<usize> = gene_list.iter().copied().collect();
    let mut reference_genes = std::collections::HashSet::new();

    for &g in gene_list {
        let rank = all_genes.iter().position(|&x| x == g).unwrap_or(0);
        let window = n_reference;
        let start = rank.saturating_sub(window / 2);
        let end = (rank + window / 2 + 1).min(all_genes.len());
        for &candidate in &all_genes[start..end] {
            if !gene_set.contains(&candidate) {
                reference_genes.insert(candidate);
            }
        }
    }

    let ref_genes: Vec<usize> = reference_genes.into_iter().collect();

    // Score each cell
    let mut scores = vec![0.0; n_obs];
    for i in 0..n_obs {
        let sig_mean = if gene_list.is_empty() {
            0.0
        } else {
            gene_list.iter().map(|&g| adata.x().get(i, g)).sum::<f64>() / gene_list.len() as f64
        };
        let ref_mean = if ref_genes.is_empty() {
            0.0
        } else {
            ref_genes.iter().map(|&g| adata.x().get(i, g)).sum::<f64>() / ref_genes.len() as f64
        };
        scores[i] = sig_mean - ref_mean;
    }

    adata.add_obs_column(score_name, ColumnData::Numeric(scores))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_adata(data: Vec<Vec<f64>>) -> AnnData {
        let n_obs = data.len();
        let n_vars = data[0].len();
        let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{}", i)).collect();
        let var_names: Vec<String> = (0..n_vars).map(|j| format!("gene_{}", j)).collect();
        AnnData::new(MatrixData::Dense(data), obs_names, var_names).unwrap()
    }

    // ── HVG tests ──

    #[test]
    fn hvg_seurat_basic() {
        // 10 cells, 5 genes. First 2 genes high variance, rest low.
        let mut data = vec![vec![0.0; 5]; 10];
        for i in 0..10 {
            data[i][0] = (i as f64) * 10.0; // high variance gene
            data[i][1] = (i as f64) * 8.0;  // high variance gene
            data[i][2] = 1.0;               // constant
            data[i][3] = 1.1;               // near-constant
            data[i][4] = 1.2;               // near-constant
        }
        let mut adata = make_adata(data);

        let config = HvgConfig {
            n_top_genes: 2,
            method: HvgMethod::SeuratV3,
            ..Default::default()
        };
        highly_variable_genes(&mut adata, &config).unwrap();

        let hv = adata.get_var("highly_variable").unwrap().as_numeric().unwrap();
        assert_eq!(hv[0], 1.0); // high variance
        assert_eq!(hv[1], 1.0); // high variance
        assert_eq!(hv[2], 0.0); // constant
    }

    #[test]
    fn hvg_cellranger_basic() {
        let mut data = vec![vec![0.0; 5]; 10];
        for i in 0..10 {
            data[i][0] = (i as f64) * 5.0;
            data[i][1] = (i as f64) * 4.0;
            data[i][2] = 1.0;
            data[i][3] = 1.0;
            data[i][4] = 1.0;
        }
        let mut adata = make_adata(data);
        let config = HvgConfig {
            n_top_genes: 2,
            method: HvgMethod::CellRanger,
            min_mean: 0.0,
            max_mean: 100.0,
            min_disp: -10.0,
            n_bins: 5,
        };
        highly_variable_genes(&mut adata, &config).unwrap();
        let disp = adata.get_var("dispersions_norm").unwrap().as_numeric().unwrap();
        assert_eq!(disp.len(), 5);
    }

    #[test]
    fn hvg_annotates_means_variances() {
        let mut data = vec![vec![0.0; 3]; 5];
        for i in 0..5 {
            data[i][0] = i as f64;
            data[i][1] = 10.0;
            data[i][2] = 0.0;
        }
        let mut adata = make_adata(data);
        highly_variable_genes(&mut adata, &HvgConfig::default()).unwrap();
        assert!(adata.get_var("means").is_some());
        assert!(adata.get_var("variances").is_some());
    }

    #[test]
    fn hvg_too_few_cells() {
        let mut adata = make_adata(vec![vec![1.0, 2.0]]);
        let result = highly_variable_genes(&mut adata, &HvgConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn hvg_n_top_exceeds_n_vars() {
        let mut data = vec![vec![0.0; 3]; 5];
        for i in 0..5 {
            data[i][0] = i as f64;
            data[i][1] = (i * 2) as f64;
            data[i][2] = (i * 3) as f64;
        }
        let mut adata = make_adata(data);
        let config = HvgConfig {
            n_top_genes: 100, // more than n_vars
            ..Default::default()
        };
        highly_variable_genes(&mut adata, &config).unwrap();
        let hv = adata.get_var("highly_variable").unwrap().as_numeric().unwrap();
        // All genes should be selected
        assert!(hv.iter().all(|&v| v == 1.0));
    }

    // ── Normalize tests ──

    #[test]
    fn normalize_total_basic() {
        let mut adata = make_adata(vec![
            vec![1.0, 2.0, 3.0], // sum = 6
            vec![4.0, 5.0, 6.0], // sum = 15
        ]);
        let config = NormalizeConfig {
            target_sum: 100.0,
            log_transform: false,
            save_raw: true,
        };
        normalize_total(&mut adata, &config).unwrap();

        // Check scaling
        let sum0: f64 = (0..3).map(|j| adata.x().get(0, j)).sum();
        assert!((sum0 - 100.0).abs() < 1e-6);
        let sum1: f64 = (0..3).map(|j| adata.x().get(1, j)).sum();
        assert!((sum1 - 100.0).abs() < 1e-6);

        // Raw counts saved
        let raw = adata.get_layer("counts").unwrap();
        assert_eq!(raw.get(0, 0), 1.0);
    }

    #[test]
    fn normalize_total_with_log() {
        let mut adata = make_adata(vec![vec![10.0, 0.0]]);
        let config = NormalizeConfig {
            target_sum: 10.0,
            log_transform: true,
            save_raw: false,
        };
        normalize_total(&mut adata, &config).unwrap();
        // 10 * (10/10) = 10, log1p(10) = ln(11)
        assert!((adata.x().get(0, 0) - (11.0_f64).ln()).abs() < 1e-10);
        // 0 * (10/10) = 0, log1p(0) = ln(1) = 0
        assert!((adata.x().get(0, 1) - 0.0).abs() < 1e-10);
        assert!(adata.get_layer("counts").is_none());
    }

    #[test]
    fn normalize_total_zero_cell() {
        let mut adata = make_adata(vec![vec![0.0, 0.0], vec![1.0, 1.0]]);
        let config = NormalizeConfig {
            target_sum: 100.0,
            log_transform: false,
            save_raw: false,
        };
        normalize_total(&mut adata, &config).unwrap();
        // Zero cell stays zero
        assert_eq!(adata.x().get(0, 0), 0.0);
    }

    #[test]
    fn normalize_total_sparse() {
        let s = SparseMatrix::from_triplets(
            vec![0, 0, 1],
            vec![0, 1, 0],
            vec![3.0, 7.0, 5.0],
            2,
            2,
        )
        .unwrap();
        let mut adata = AnnData::new(
            MatrixData::Sparse(s),
            vec!["c0".into(), "c1".into()],
            vec!["g0".into(), "g1".into()],
        )
        .unwrap();
        let config = NormalizeConfig {
            target_sum: 10.0,
            log_transform: false,
            save_raw: false,
        };
        normalize_total(&mut adata, &config).unwrap();
        // cell 0: sum=10, factor=1.0; cell 1: sum=5, factor=2.0
        assert!((adata.x().get(0, 0) - 3.0).abs() < 1e-10);
        assert!((adata.x().get(1, 0) - 10.0).abs() < 1e-10);
    }

    // ── Regress out tests ──

    #[test]
    fn regress_out_basic() {
        // Gene expression correlated with a confounder
        let n = 20;
        let mut data = vec![vec![0.0; 2]; n];
        let mut depth = vec![0.0; n];
        for i in 0..n {
            let d = (i + 1) as f64;
            depth[i] = d;
            data[i][0] = d * 2.0 + 1.0; // perfectly correlated with depth
            data[i][1] = 5.0; // constant
        }
        let mut adata = make_adata(data);
        adata
            .add_obs_column("depth", ColumnData::Numeric(depth))
            .unwrap();

        regress_out(&mut adata, &["depth"]).unwrap();

        // After regressing out depth, gene 0 residuals should be ~0
        let residuals: Vec<f64> = (0..n).map(|i| adata.x().get(i, 0)).collect();
        let max_resid = residuals.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        assert!(max_resid < 1e-10, "max residual = {}", max_resid);

        // Gene 1 was constant, so residuals should be ~0
        let residuals1: Vec<f64> = (0..n).map(|i| adata.x().get(i, 1)).collect();
        let max_resid1 = residuals1.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        assert!(max_resid1 < 1e-10, "max residual = {}", max_resid1);
    }

    #[test]
    fn regress_out_no_keys() {
        let mut adata = make_adata(vec![vec![1.0, 2.0]]);
        regress_out(&mut adata, &[]).unwrap();
        assert_eq!(adata.x().get(0, 0), 1.0);
    }

    #[test]
    fn regress_out_missing_key() {
        let mut adata = make_adata(vec![vec![1.0, 2.0]]);
        let result = regress_out(&mut adata, &["missing"]);
        assert!(result.is_err());
    }

    #[test]
    fn regress_out_non_numeric() {
        let mut adata = make_adata(vec![vec![1.0]; 3]);
        adata
            .add_obs("batch", vec!["a".into(), "b".into(), "c".into()])
            .unwrap();
        let result = regress_out(&mut adata, &["batch"]);
        assert!(result.is_err());
    }

    // ── Scrublet tests ──

    #[test]
    fn scrublet_basic() {
        // Create a dataset where some cells are obvious doublets
        let n = 30;
        let n_vars = 10;
        let mut data = vec![vec![0.0; n_vars]; n];
        for i in 0..n {
            for j in 0..n_vars {
                // Two distinct clusters
                if i < n / 2 {
                    data[i][j] = if j < n_vars / 2 { 10.0 } else { 0.0 };
                } else {
                    data[i][j] = if j < n_vars / 2 { 0.0 } else { 10.0 };
                }
            }
        }
        let mut adata = make_adata(data);
        scrublet_doublets(&mut adata, &ScrubletConfig::default()).unwrap();

        let scores = adata
            .get_obs("doublet_score")
            .unwrap()
            .as_numeric()
            .unwrap();
        assert_eq!(scores.len(), n);
        // All scores should be between 0 and 1
        assert!(scores.iter().all(|&s| (0.0..=1.0).contains(&s)));
    }

    #[test]
    fn scrublet_too_few_cells() {
        let mut adata = make_adata(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        let result = scrublet_doublets(&mut adata, &ScrubletConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn scrublet_stores_predictions() {
        let n = 20;
        let data: Vec<Vec<f64>> = (0..n)
            .map(|i| vec![i as f64, (i * 2) as f64, (i * 3) as f64])
            .collect();
        let mut adata = make_adata(data);
        scrublet_doublets(&mut adata, &ScrubletConfig { n_pcs: 2, ..Default::default() }).unwrap();
        assert!(adata.get_obs("predicted_doublet").is_some());
        let preds = adata.get_obs("predicted_doublet").unwrap().as_numeric().unwrap();
        assert!(preds.iter().all(|&p| p == 0.0 || p == 1.0));
    }

    // ── Score genes tests ──

    #[test]
    fn score_genes_basic() {
        let data = vec![
            vec![10.0, 0.0, 0.0, 0.0],
            vec![0.0, 10.0, 0.0, 0.0],
            vec![0.0, 0.0, 10.0, 0.0],
        ];
        let mut adata = make_adata(data);
        score_genes(&mut adata, &[0], 2, "sig_score").unwrap();
        let scores = adata.get_obs("sig_score").unwrap().as_numeric().unwrap();
        assert_eq!(scores.len(), 3);
        // Cell 0 highly expresses gene 0 → high score
        assert!(scores[0] > scores[1]);
        assert!(scores[0] > scores[2]);
    }

    #[test]
    fn score_genes_empty_list() {
        let mut adata = make_adata(vec![vec![1.0, 2.0]]);
        score_genes(&mut adata, &[], 5, "empty").unwrap();
        let scores = adata.get_obs("empty").unwrap().as_numeric().unwrap();
        assert_eq!(scores[0], 0.0);
    }

    #[test]
    fn score_genes_invalid_index() {
        let mut adata = make_adata(vec![vec![1.0, 2.0]]);
        let result = score_genes(&mut adata, &[10], 5, "bad");
        assert!(result.is_err());
    }

    #[test]
    fn score_genes_all_genes_as_signature() {
        let data = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        let mut adata = make_adata(data);
        // All genes in signature → no reference genes, reference mean = 0
        score_genes(&mut adata, &[0, 1], 5, "all_sig").unwrap();
        let scores = adata.get_obs("all_sig").unwrap().as_numeric().unwrap();
        // Score = mean(all genes) - 0
        assert!((scores[0] - 1.5).abs() < 1e-10);
        assert!((scores[1] - 3.5).abs() < 1e-10);
    }

    // ── Helper tests ──

    #[test]
    fn lcg_produces_different_values() {
        let mut state = 12345u64;
        let a = lcg_next(&mut state);
        let b = lcg_next(&mut state);
        let c = lcg_next(&mut state);
        assert_ne!(a, b);
        assert_ne!(b, c);
    }

    #[test]
    fn otsu_threshold_bimodal() {
        let mut scores = vec![0.0; 100];
        for s in scores.iter_mut().take(80) {
            *s = 0.1;
        }
        for s in scores.iter_mut().skip(80) {
            *s = 0.9;
        }
        let t = otsu_threshold(&scores);
        assert!(t > 0.1 && t < 0.9, "threshold = {}", t);
    }

    #[test]
    fn simple_median_odd() {
        assert_eq!(simple_median(&[3.0, 1.0, 2.0]), 2.0);
    }

    #[test]
    fn simple_median_even() {
        assert_eq!(simple_median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
    }

    #[test]
    fn simple_median_empty() {
        assert_eq!(simple_median(&[]), 0.0);
    }

    #[test]
    fn invert_small_matrix_2x2() {
        // [[2, 1], [1, 3]] → inverse [[3/5, -1/5], [-1/5, 2/5]]
        let m = vec![2.0, 1.0, 1.0, 3.0];
        let inv = invert_small_matrix(&m, 2).unwrap();
        assert!((inv[0] - 0.6).abs() < 1e-10);
        assert!((inv[1] - (-0.2)).abs() < 1e-10);
        assert!((inv[2] - (-0.2)).abs() < 1e-10);
        assert!((inv[3] - 0.4).abs() < 1e-10);
    }

    #[test]
    fn invert_singular_matrix() {
        let m = vec![1.0, 2.0, 2.0, 4.0]; // singular
        assert!(invert_small_matrix(&m, 2).is_none());
    }

    // ── Normalize + HVG pipeline ──

    #[test]
    fn normalize_then_hvg() {
        let mut data = vec![vec![0.0; 5]; 10];
        for i in 0..10 {
            data[i][0] = (i as f64 + 1.0) * 100.0; // high expression, high variance
            data[i][1] = (i as f64 + 1.0) * 50.0;
            data[i][2] = 1.0;
            data[i][3] = 1.0;
            data[i][4] = 1.0;
        }
        let mut adata = make_adata(data);
        normalize_total(&mut adata, &NormalizeConfig::default()).unwrap();
        highly_variable_genes(
            &mut adata,
            &HvgConfig {
                n_top_genes: 2,
                ..Default::default()
            },
        )
        .unwrap();
        let hv = adata.get_var("highly_variable").unwrap().as_numeric().unwrap();
        // Should have exactly 2 HVGs
        let n_hvg: usize = hv.iter().filter(|&&v| v == 1.0).count();
        assert_eq!(n_hvg, 2);
    }

    #[test]
    fn regress_out_multiple_covariates() {
        let n = 30;
        let mut data = vec![vec![0.0; 2]; n];
        let mut cov1 = vec![0.0; n];
        let mut cov2 = vec![0.0; n];
        for i in 0..n {
            cov1[i] = i as f64;
            cov2[i] = (i as f64).powi(2) / 100.0;
            data[i][0] = 3.0 * cov1[i] + 2.0 * cov2[i] + 5.0;
            data[i][1] = 10.0;
        }
        let mut adata = make_adata(data);
        adata.add_obs_column("cov1", ColumnData::Numeric(cov1)).unwrap();
        adata.add_obs_column("cov2", ColumnData::Numeric(cov2)).unwrap();

        regress_out(&mut adata, &["cov1", "cov2"]).unwrap();

        let max_resid = (0..n).map(|i| adata.x().get(i, 0).abs()).fold(0.0f64, f64::max);
        assert!(max_resid < 1e-8, "max residual = {}", max_resid);
    }
}
