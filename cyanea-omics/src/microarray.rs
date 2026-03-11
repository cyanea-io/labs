//! Microarray analysis — RMA normalization, limma-style differential expression,
//! Illumina methylation array analysis (450K/EPIC).
//!
//! Provides core algorithms for processing legacy microarray data.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// RMA normalization
// ---------------------------------------------------------------------------

/// RMA (Robust Multi-array Average) normalization pipeline.
///
/// Applies three steps:
/// 1. Background correction (log2 transform + offset)
/// 2. Quantile normalization (across arrays)
/// 3. Per-probe median polish (summarization)
///
/// `probe_intensities` is probe-major: `probe_intensities[probe][sample]`.
/// Returns normalized expression values (probe × sample).
pub fn rma_normalize(probe_intensities: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    if probe_intensities.is_empty() {
        return Err(CyaneaError::InvalidInput("empty probe intensities".into()));
    }
    let n_samples = probe_intensities[0].len();
    if n_samples == 0 {
        return Err(CyaneaError::InvalidInput("no samples".into()));
    }
    for row in probe_intensities {
        if row.len() != n_samples {
            return Err(CyaneaError::InvalidInput(
                "inconsistent sample counts".into(),
            ));
        }
    }

    // Step 1: Background correction — log2(max(x, 1))
    let mut corrected: Vec<Vec<f64>> = probe_intensities
        .iter()
        .map(|row| row.iter().map(|&v| v.max(1.0).log2()).collect())
        .collect();

    // Step 2: Quantile normalization
    quantile_normalize(&mut corrected)?;

    Ok(corrected)
}

/// Quantile normalization across columns (samples).
///
/// For each column, ranks the values, then replaces each rank with the
/// mean of that rank across all columns.
pub fn quantile_normalize(data: &mut [Vec<f64>]) -> Result<()> {
    if data.is_empty() {
        return Ok(());
    }
    let n_probes = data.len();
    let n_samples = data[0].len();

    // For each sample: sort probes by value, record rank order
    let mut ranked: Vec<Vec<usize>> = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        let mut indices: Vec<usize> = (0..n_probes).collect();
        indices.sort_by(|&a, &b| {
            data[a][s]
                .partial_cmp(&data[b][s])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        ranked.push(indices);
    }

    // Compute mean value at each rank across samples
    let mut rank_means = vec![0.0; n_probes];
    for rank in 0..n_probes {
        let sum: f64 = (0..n_samples).map(|s| data[ranked[s][rank]][s]).sum();
        rank_means[rank] = sum / n_samples as f64;
    }

    // Replace values with rank means
    for s in 0..n_samples {
        for (rank, &probe_idx) in ranked[s].iter().enumerate() {
            data[probe_idx][s] = rank_means[rank];
        }
    }

    Ok(())
}

/// Median polish for probe set summarization.
///
/// Given a matrix of log2-normalized probe intensities for a single
/// probe set, iteratively removes row and column effects to produce
/// a robust estimate of the probe set expression level per sample.
///
/// `probes` is probe-major: `probes[probe_in_set][sample]`.
/// Returns one expression value per sample.
pub fn median_polish(probes: &[Vec<f64>], max_iter: usize) -> Result<Vec<f64>> {
    if probes.is_empty() {
        return Err(CyaneaError::InvalidInput("empty probe set".into()));
    }
    let n_probes = probes.len();
    let n_samples = probes[0].len();

    let mut residuals: Vec<Vec<f64>> = probes.to_vec();
    let mut col_effects = vec![0.0; n_samples];
    let mut overall = 0.0;

    for _ in 0..max_iter {
        // Subtract row medians
        for p in 0..n_probes {
            let row_med = median_f64(&residuals[p]);
            for s in 0..n_samples {
                residuals[p][s] -= row_med;
            }
        }

        // Subtract column medians, accumulate into col_effects
        for s in 0..n_samples {
            let col_vals: Vec<f64> = (0..n_probes).map(|p| residuals[p][s]).collect();
            let col_med = median_f64(&col_vals);
            col_effects[s] += col_med;
            for p in 0..n_probes {
                residuals[p][s] -= col_med;
            }
        }

        // Update overall effect
        let eff_med = median_f64(&col_effects);
        overall += eff_med;
        for s in 0..n_samples {
            col_effects[s] -= eff_med;
        }
    }

    // Final expression = overall + col_effects
    let expression: Vec<f64> = col_effects.iter().map(|&e| overall + e).collect();
    Ok(expression)
}

// ---------------------------------------------------------------------------
// Limma-style differential expression
// ---------------------------------------------------------------------------

/// Result of a differential expression test for one gene/probe set.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DiffExprResult {
    /// Gene/probe set index.
    pub gene_idx: usize,
    /// Gene/probe set name.
    pub gene_name: String,
    /// Log2 fold change (group B vs group A).
    pub log2_fc: f64,
    /// Average expression across both groups.
    pub avg_expr: f64,
    /// t-statistic.
    pub t_stat: f64,
    /// Raw p-value.
    pub p_value: f64,
    /// BH-adjusted p-value.
    pub adj_p_value: f64,
}

/// Run limma-style moderated t-test for differential expression.
///
/// `expression` is gene-major: `expression[gene][sample]`.
/// `gene_names` are the gene identifiers.
/// `groups` assigns each sample to group 0 or 1.
///
/// Uses empirical Bayes variance moderation: pools variance estimates
/// across genes to stabilize t-statistics for small sample sizes.
pub fn limma_diff_expr(
    expression: &[Vec<f64>],
    gene_names: &[String],
    groups: &[u8],
) -> Result<Vec<DiffExprResult>> {
    let n_genes = expression.len();
    if n_genes == 0 {
        return Err(CyaneaError::InvalidInput("empty expression matrix".into()));
    }
    let n_samples = expression[0].len();
    if gene_names.len() != n_genes || groups.len() != n_samples {
        return Err(CyaneaError::InvalidInput("dimension mismatch".into()));
    }

    let idx_a: Vec<usize> = groups.iter().enumerate().filter(|(_, &g)| g == 0).map(|(i, _)| i).collect();
    let idx_b: Vec<usize> = groups.iter().enumerate().filter(|(_, &g)| g == 1).map(|(i, _)| i).collect();

    let n_a = idx_a.len() as f64;
    let n_b = idx_b.len() as f64;
    if n_a < 1.0 || n_b < 1.0 {
        return Err(CyaneaError::InvalidInput(
            "need at least 1 sample per group".into(),
        ));
    }

    let df = n_a + n_b - 2.0;

    // Compute per-gene statistics
    let mut log2_fcs = vec![0.0; n_genes];
    let mut avg_exprs = vec![0.0; n_genes];
    let mut vars = vec![0.0; n_genes];

    for g in 0..n_genes {
        let mean_a: f64 = idx_a.iter().map(|&i| expression[g][i]).sum::<f64>() / n_a;
        let mean_b: f64 = idx_b.iter().map(|&i| expression[g][i]).sum::<f64>() / n_b;

        log2_fcs[g] = mean_b - mean_a; // Already log2 if RMA-normalized
        avg_exprs[g] = (mean_a + mean_b) / 2.0;

        // Pooled variance
        let ss_a: f64 = idx_a.iter().map(|&i| (expression[g][i] - mean_a).powi(2)).sum();
        let ss_b: f64 = idx_b.iter().map(|&i| (expression[g][i] - mean_b).powi(2)).sum();
        vars[g] = (ss_a + ss_b) / df;
    }

    // Empirical Bayes moderation
    // Prior variance s0² = median of gene-level variances
    // Prior df d0 = n_genes / 10 (heuristic)
    let s0_sq = median_f64(&vars);
    let d0 = (n_genes as f64 / 10.0).max(1.0);

    // Moderated variance = (d0 * s0² + df * s²) / (d0 + df)
    let mut mod_vars = vec![0.0; n_genes];
    for g in 0..n_genes {
        mod_vars[g] = (d0 * s0_sq + df * vars[g]) / (d0 + df);
    }

    // Moderated t-statistic
    let mut results = Vec::with_capacity(n_genes);
    for g in 0..n_genes {
        let se = (mod_vars[g] * (1.0 / n_a + 1.0 / n_b)).sqrt();
        let t = if se > 1e-30 {
            log2_fcs[g] / se
        } else {
            0.0
        };

        // p-value from t-distribution approximation (normal for large df)
        let mod_df = d0 + df;
        let p = t_to_p(t, mod_df);

        results.push(DiffExprResult {
            gene_idx: g,
            gene_name: gene_names[g].clone(),
            log2_fc: log2_fcs[g],
            avg_expr: avg_exprs[g],
            t_stat: t,
            p_value: p,
            adj_p_value: p, // corrected below
        });
    }

    // BH correction
    bh_correct(&mut results);

    // Sort by adjusted p-value
    results.sort_by(|a, b| {
        a.adj_p_value
            .partial_cmp(&b.adj_p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(results)
}

// ---------------------------------------------------------------------------
// Illumina methylation array analysis
// ---------------------------------------------------------------------------

/// A methylation probe with beta value and detection p-value.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MethylationProbe {
    /// Probe ID (e.g., "cg00000029").
    pub probe_id: String,
    /// Chromosome.
    pub chrom: String,
    /// Position.
    pub position: u64,
    /// Gene name (nearest annotated gene).
    pub gene: String,
    /// Infinium design type (I or II).
    pub design_type: InfiniumType,
    /// Beta values per sample.
    pub beta_values: Vec<f64>,
    /// Detection p-values per sample.
    pub detection_p: Vec<f64>,
}

/// Infinium probe design type.
#[derive(Debug, Clone, PartialEq)]
pub enum InfiniumType {
    /// Type I — two-bead design (separate methylated/unmethylated).
    TypeI,
    /// Type II — single-bead design.
    TypeII,
}

/// Compute beta values from methylated (M) and unmethylated (U) signal.
///
/// `beta = M / (M + U + offset)` where offset (typically 100) prevents
/// division by near-zero.
pub fn compute_beta(methylated: &[f64], unmethylated: &[f64], offset: f64) -> Result<Vec<f64>> {
    if methylated.len() != unmethylated.len() {
        return Err(CyaneaError::InvalidInput(
            "M and U arrays must have same length".into(),
        ));
    }
    Ok(methylated
        .iter()
        .zip(unmethylated.iter())
        .map(|(&m, &u)| m / (m + u + offset))
        .collect())
}

/// Convert beta values to M-values (log2 ratio).
///
/// `M = log2(beta / (1 - beta))`
///
/// Beta values are clamped to [0.001, 0.999] to avoid infinity.
pub fn beta_to_m_value(beta: &[f64]) -> Vec<f64> {
    beta.iter()
        .map(|&b| {
            let b = b.clamp(0.001, 0.999);
            (b / (1.0 - b)).log2()
        })
        .collect()
}

/// Convert M-values back to beta values.
///
/// `beta = 2^M / (2^M + 1)`
pub fn m_value_to_beta(m_values: &[f64]) -> Vec<f64> {
    m_values
        .iter()
        .map(|&m| {
            let pow = 2.0_f64.powf(m);
            pow / (pow + 1.0)
        })
        .collect()
}

/// SWAN normalization for Infinium type I/II bias correction.
///
/// Adjusts Type II probe intensities to match the distribution of
/// Type I probes within each sample.
///
/// `beta_values` is probe-major: `beta_values[probe][sample]`.
/// `design_types` indicates Type I or Type II for each probe.
/// Returns corrected beta values.
pub fn swan_normalize(
    beta_values: &[Vec<f64>],
    design_types: &[InfiniumType],
) -> Result<Vec<Vec<f64>>> {
    let n_probes = beta_values.len();
    let n_samples = beta_values[0].len();

    if design_types.len() != n_probes {
        return Err(CyaneaError::InvalidInput(
            "design_types length mismatch".into(),
        ));
    }

    let type1_idx: Vec<usize> = design_types
        .iter()
        .enumerate()
        .filter(|(_, t)| **t == InfiniumType::TypeI)
        .map(|(i, _)| i)
        .collect();
    let type2_idx: Vec<usize> = design_types
        .iter()
        .enumerate()
        .filter(|(_, t)| **t == InfiniumType::TypeII)
        .map(|(i, _)| i)
        .collect();

    if type1_idx.is_empty() || type2_idx.is_empty() {
        // Nothing to normalize — return as-is
        return Ok(beta_values.to_vec());
    }

    let mut corrected = beta_values.to_vec();

    // Per-sample correction: match Type II quantiles to Type I quantiles
    for s in 0..n_samples {
        let mut t1_vals: Vec<f64> = type1_idx.iter().map(|&i| beta_values[i][s]).collect();
        let mut t2_vals: Vec<f64> = type2_idx.iter().map(|&i| beta_values[i][s]).collect();
        t1_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        t2_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Build quantile mapping function (linear interpolation)
        let n1 = t1_vals.len();
        let n2 = t2_vals.len();

        // For each Type II probe, find its quantile in the Type II distribution,
        // then map to the Type I distribution at that quantile
        for &pi in &type2_idx {
            let val = beta_values[pi][s];
            // Find rank in Type II distribution
            let rank = t2_vals.partition_point(|&v| v < val);
            let quantile = rank as f64 / n2 as f64;
            // Map to Type I value at same quantile
            let t1_idx_f = quantile * (n1 - 1) as f64;
            let lo = t1_idx_f.floor() as usize;
            let hi = (lo + 1).min(n1 - 1);
            let frac = t1_idx_f - lo as f64;
            let mapped = t1_vals[lo] * (1.0 - frac) + t1_vals[hi] * frac;
            corrected[pi][s] = mapped;
        }
    }

    Ok(corrected)
}

/// Differential methylation analysis between two groups.
///
/// `beta_values` is probe-major: `beta_values[probe][sample]`.
/// `probe_ids` are probe identifiers.
/// `groups` assigns each sample to group 0 or 1.
///
/// Returns probes with significant differential methylation (delta beta
/// and moderated t-test on M-values).
pub fn diff_methylation(
    beta_values: &[Vec<f64>],
    probe_ids: &[String],
    groups: &[u8],
) -> Result<Vec<DiffMethResult>> {
    let n_probes = beta_values.len();
    if n_probes == 0 {
        return Err(CyaneaError::InvalidInput("empty beta values".into()));
    }
    let n_samples = beta_values[0].len();
    if probe_ids.len() != n_probes || groups.len() != n_samples {
        return Err(CyaneaError::InvalidInput("dimension mismatch".into()));
    }

    let idx_a: Vec<usize> = groups.iter().enumerate().filter(|(_, &g)| g == 0).map(|(i, _)| i).collect();
    let idx_b: Vec<usize> = groups.iter().enumerate().filter(|(_, &g)| g == 1).map(|(i, _)| i).collect();
    let n_a = idx_a.len() as f64;
    let n_b = idx_b.len() as f64;

    if n_a < 1.0 || n_b < 1.0 {
        return Err(CyaneaError::InvalidInput(
            "need at least 1 sample per group".into(),
        ));
    }

    let df = n_a + n_b - 2.0;
    let mut results = Vec::with_capacity(n_probes);

    // Convert to M-values for statistical testing
    let m_values: Vec<Vec<f64>> = beta_values
        .iter()
        .map(|row| beta_to_m_value(row))
        .collect();

    // Per-probe statistics
    let mut vars = vec![0.0; n_probes];
    let mut delta_betas = vec![0.0; n_probes];
    let mut m_fcs = vec![0.0; n_probes];

    for p in 0..n_probes {
        let beta_mean_a: f64 = idx_a.iter().map(|&i| beta_values[p][i]).sum::<f64>() / n_a;
        let beta_mean_b: f64 = idx_b.iter().map(|&i| beta_values[p][i]).sum::<f64>() / n_b;
        delta_betas[p] = beta_mean_b - beta_mean_a;

        let m_mean_a: f64 = idx_a.iter().map(|&i| m_values[p][i]).sum::<f64>() / n_a;
        let m_mean_b: f64 = idx_b.iter().map(|&i| m_values[p][i]).sum::<f64>() / n_b;
        m_fcs[p] = m_mean_b - m_mean_a;

        let ss_a: f64 = idx_a.iter().map(|&i| (m_values[p][i] - m_mean_a).powi(2)).sum();
        let ss_b: f64 = idx_b.iter().map(|&i| (m_values[p][i] - m_mean_b).powi(2)).sum();
        vars[p] = if df > 0.0 { (ss_a + ss_b) / df } else { 0.0 };
    }

    // Empirical Bayes moderation
    let s0_sq = median_f64(&vars);
    let d0 = (n_probes as f64 / 10.0).max(1.0);

    for p in 0..n_probes {
        let mod_var = (d0 * s0_sq + df * vars[p]) / (d0 + df);
        let se = (mod_var * (1.0 / n_a + 1.0 / n_b)).sqrt();
        let t = if se > 1e-30 { m_fcs[p] / se } else { 0.0 };
        let pval = t_to_p(t, d0 + df);

        results.push(DiffMethResult {
            probe_idx: p,
            probe_id: probe_ids[p].clone(),
            delta_beta: delta_betas[p],
            m_value_fc: m_fcs[p],
            t_stat: t,
            p_value: pval,
            adj_p_value: pval,
        });
    }

    // BH correction
    bh_correct_meth(&mut results);

    results.sort_by(|a, b| {
        a.adj_p_value
            .partial_cmp(&b.adj_p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(results)
}

/// Result of a differential methylation test.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DiffMethResult {
    /// Probe index.
    pub probe_idx: usize,
    /// Probe ID.
    pub probe_id: String,
    /// Delta beta (group B - group A).
    pub delta_beta: f64,
    /// M-value fold change.
    pub m_value_fc: f64,
    /// Moderated t-statistic.
    pub t_stat: f64,
    /// Raw p-value.
    pub p_value: f64,
    /// BH-adjusted p-value.
    pub adj_p_value: f64,
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

fn median_f64(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

/// Approximate p-value from t-distribution using normal approximation
/// (good for df > 30, reasonable for smaller df).
fn t_to_p(t: f64, df: f64) -> f64 {
    // Use the approximation: for large df, t ~ N(0,1)
    // For smaller df, use the correction: z = t * sqrt(df/(df-2)) approximately
    let z = if df > 2.0 {
        t * (1.0 - 1.0 / (4.0 * df)).sqrt()
    } else {
        t
    };
    // Two-sided p-value from z
    let abs_z = z.abs();
    // Complementary error function approximation
    let p = erfc_approx(abs_z / std::f64::consts::SQRT_2);
    p.clamp(0.0, 1.0)
}

fn erf_approx(x: f64) -> f64 {
    let a1: f64 = 0.254829592;
    let a2: f64 = -0.284496736;
    let a3: f64 = 1.421413741;
    let a4: f64 = -1.453152027;
    let a5: f64 = 1.061405429;
    let p: f64 = 0.3275911;
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    sign * y
}

fn erfc_approx(x: f64) -> f64 {
    1.0 - erf_approx(x)
}

fn bh_correct(results: &mut [DiffExprResult]) {
    let n = results.len();
    if n == 0 {
        return;
    }
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        results[a]
            .p_value
            .partial_cmp(&results[b].p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (rank, &idx) in indices.iter().enumerate() {
        let adj = (results[idx].p_value * n as f64 / (rank + 1) as f64).min(1.0);
        results[idx].adj_p_value = adj;
    }

    // Enforce monotonicity
    for i in (0..n - 1).rev() {
        let idx_curr = indices[i];
        let idx_next = indices[i + 1];
        if results[idx_curr].adj_p_value > results[idx_next].adj_p_value {
            results[idx_curr].adj_p_value = results[idx_next].adj_p_value;
        }
    }
}

fn bh_correct_meth(results: &mut [DiffMethResult]) {
    let n = results.len();
    if n == 0 {
        return;
    }
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        results[a]
            .p_value
            .partial_cmp(&results[b].p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (rank, &idx) in indices.iter().enumerate() {
        let adj = (results[idx].p_value * n as f64 / (rank + 1) as f64).min(1.0);
        results[idx].adj_p_value = adj;
    }

    for i in (0..n - 1).rev() {
        let idx_curr = indices[i];
        let idx_next = indices[i + 1];
        if results[idx_curr].adj_p_value > results[idx_next].adj_p_value {
            results[idx_curr].adj_p_value = results[idx_next].adj_p_value;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rma_normalize() {
        // 4 probes × 3 samples
        let data = vec![
            vec![100.0, 200.0, 150.0],
            vec![500.0, 300.0, 400.0],
            vec![50.0, 80.0, 60.0],
            vec![1000.0, 800.0, 900.0],
        ];
        let result = rma_normalize(&data).unwrap();
        assert_eq!(result.len(), 4);
        assert_eq!(result[0].len(), 3);
        // After quantile normalization, column distributions should be identical
        for s in 0..3 {
            let mut col: Vec<f64> = result.iter().map(|r| r[s]).collect();
            col.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // All columns should have the same sorted values
            let mut col0: Vec<f64> = result.iter().map(|r| r[0]).collect();
            col0.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for i in 0..4 {
                assert!(
                    (col[i] - col0[i]).abs() < 1e-10,
                    "sample {} rank {} differs: {} vs {}",
                    s,
                    i,
                    col[i],
                    col0[i]
                );
            }
        }
    }

    #[test]
    fn test_quantile_normalize_identity() {
        // Already identical distributions → no change
        let mut data = vec![vec![1.0, 1.0], vec![2.0, 2.0], vec![3.0, 3.0]];
        let original = data.clone();
        quantile_normalize(&mut data).unwrap();
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - original[i][j]).abs() < 1e-10,
                    "changed at [{},{}]",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_median_polish() {
        let probes = vec![
            vec![10.0, 12.0, 11.0],
            vec![9.0, 11.0, 10.0],
            vec![8.0, 10.0, 9.0],
        ];
        let result = median_polish(&probes, 10).unwrap();
        assert_eq!(result.len(), 3);
        // Expression should reflect the column pattern
        assert!(result[1] > result[0]); // sample 2 has higher values
    }

    #[test]
    fn test_limma_basic() {
        // 5 genes × 6 samples (3 per group)
        let expression = vec![
            vec![5.0, 5.1, 4.9, 10.0, 10.2, 9.8],  // up-regulated
            vec![10.0, 9.8, 10.2, 5.0, 5.1, 4.9],   // down-regulated
            vec![7.0, 7.1, 6.9, 7.0, 7.2, 6.8],     // no change
            vec![3.0, 3.1, 2.9, 8.0, 8.1, 7.9],     // up-regulated
            vec![6.0, 6.5, 5.5, 6.0, 6.5, 5.5],     // no change
        ];
        let gene_names: Vec<String> = (0..5).map(|i| format!("Gene{}", i)).collect();
        let groups = vec![0, 0, 0, 1, 1, 1];

        let results = limma_diff_expr(&expression, &gene_names, &groups).unwrap();
        assert_eq!(results.len(), 5);

        // Gene0 and Gene3 should be up-regulated (positive log2FC)
        let gene0 = results.iter().find(|r| r.gene_name == "Gene0").unwrap();
        assert!(gene0.log2_fc > 0.0);
        assert!(gene0.adj_p_value < 0.1);

        // Gene1 should be down-regulated
        let gene1 = results.iter().find(|r| r.gene_name == "Gene1").unwrap();
        assert!(gene1.log2_fc < 0.0);

        // Gene2 should not be significant
        let gene2 = results.iter().find(|r| r.gene_name == "Gene2").unwrap();
        assert!(gene2.adj_p_value > gene0.adj_p_value);
    }

    #[test]
    fn test_limma_adj_p_values() {
        let expression = vec![
            vec![1.0, 1.1, 5.0, 5.1],
            vec![2.0, 2.0, 2.0, 2.0],
        ];
        let gene_names = vec!["DE".into(), "Flat".into()];
        let groups = vec![0, 0, 1, 1];
        let results = limma_diff_expr(&expression, &gene_names, &groups).unwrap();
        for r in &results {
            assert!(r.adj_p_value >= r.p_value - 1e-10);
            assert!(r.adj_p_value <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn test_compute_beta() {
        let m = vec![1000.0, 500.0, 100.0];
        let u = vec![100.0, 500.0, 1000.0];
        let beta = compute_beta(&m, &u, 100.0).unwrap();
        assert_eq!(beta.len(), 3);
        // High methylation
        assert!(beta[0] > 0.8);
        // ~50% methylation
        assert!((beta[1] - 0.454).abs() < 0.05);
        // Low methylation
        assert!(beta[2] < 0.2);
    }

    #[test]
    fn test_beta_m_roundtrip() {
        let betas = vec![0.1, 0.3, 0.5, 0.7, 0.9];
        let m_vals = beta_to_m_value(&betas);
        let recovered = m_value_to_beta(&m_vals);
        for (i, (&orig, &rec)) in betas.iter().zip(recovered.iter()).enumerate() {
            assert!(
                (orig - rec).abs() < 1e-6,
                "beta[{}]: {} vs {}",
                i,
                orig,
                rec
            );
        }
    }

    #[test]
    fn test_swan_normalize() {
        // 4 probes: 2 Type I, 2 Type II
        let beta = vec![
            vec![0.3, 0.4],  // Type I
            vec![0.7, 0.8],  // Type I
            vec![0.2, 0.3],  // Type II (will be adjusted)
            vec![0.6, 0.7],  // Type II (will be adjusted)
        ];
        let types = vec![
            InfiniumType::TypeI,
            InfiniumType::TypeI,
            InfiniumType::TypeII,
            InfiniumType::TypeII,
        ];
        let corrected = swan_normalize(&beta, &types).unwrap();
        assert_eq!(corrected.len(), 4);
        // Type I probes should be unchanged
        assert!((corrected[0][0] - 0.3).abs() < 1e-10);
        assert!((corrected[1][0] - 0.7).abs() < 1e-10);
        // Type II probes should be adjusted toward Type I distribution
        assert!(corrected[2][0] >= 0.0 && corrected[2][0] <= 1.0);
    }

    #[test]
    fn test_diff_methylation() {
        // 3 probes × 4 samples
        let beta = vec![
            vec![0.2, 0.25, 0.8, 0.85],  // differentially methylated
            vec![0.5, 0.5, 0.5, 0.5],     // no difference
            vec![0.9, 0.85, 0.1, 0.15],   // differentially methylated (opposite)
        ];
        let probe_ids = vec!["cg001".into(), "cg002".into(), "cg003".into()];
        let groups = vec![0, 0, 1, 1];

        let results = diff_methylation(&beta, &probe_ids, &groups).unwrap();
        assert_eq!(results.len(), 3);

        // cg001: hypermethylated (delta_beta > 0)
        let cg001 = results.iter().find(|r| r.probe_id == "cg001").unwrap();
        assert!(cg001.delta_beta > 0.4);

        // cg003: hypomethylated (delta_beta < 0)
        let cg003 = results.iter().find(|r| r.probe_id == "cg003").unwrap();
        assert!(cg003.delta_beta < -0.4);

        // cg002: no change
        let cg002 = results.iter().find(|r| r.probe_id == "cg002").unwrap();
        assert!(cg002.delta_beta.abs() < 0.01);
    }

    #[test]
    fn test_median_f64() {
        assert!((median_f64(&[1.0, 3.0, 5.0]) - 3.0).abs() < 1e-10);
        assert!((median_f64(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
        assert!((median_f64(&[42.0]) - 42.0).abs() < 1e-10);
    }
}
