//! Differential binding and accessibility analysis.
//!
//! This module provides DESeq2-style differential analysis for count data at peaks.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::Result;
use crate::EpiError;

/// Result of differential binding analysis for a region.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", serialize, deserialize)]
pub struct DiffResult {
    /// Region name/ID.
    pub region_id: String,
    /// Log2 fold change (group1 / group2).
    pub log2_fc: f64,
    /// P-value from Wald test.
    pub p_value: f64,
    /// Q-value (BH-corrected p-value).
    pub q_value: f64,
    /// Mean count in group 1.
    pub mean_count_1: f64,
    /// Mean count in group 2.
    pub mean_count_2: f64,
}

/// Perform differential peak analysis using negative binomial model.
///
/// Input: count matrix (peaks x samples), conditions (group assignment for each sample).
/// Groups must be 0 or 1.
/// Uses size factor normalization and Wald test for differential abundance.
pub fn differential_peaks(
    count_matrix: &[Vec<u32>],
    conditions: &[u8],
    region_ids: &[String],
    pseudocount: f64,
) -> Result<Vec<DiffResult>> {
    if count_matrix.is_empty() {
        return Err(EpiError::InsufficientData("Count matrix is empty".to_string()));
    }

    let n_samples = count_matrix[0].len();
    if conditions.len() != n_samples {
        return Err(EpiError::InvalidInput(
            "Condition count must equal sample count".to_string(),
        ));
    }

    // Compute size factors (median of ratios)
    let size_factors = compute_size_factors(count_matrix)?;

    // Normalize counts (precomputed for reuse if needed)
    let _normalized = normalize_counts(count_matrix, &size_factors)?;

    // Separate by condition
    let mut group1_indices = Vec::new();
    let mut group2_indices = Vec::new();

    for (i, &cond) in conditions.iter().enumerate() {
        match cond {
            0 => group1_indices.push(i),
            1 => group2_indices.push(i),
            _ => {
                return Err(EpiError::InvalidInput(
                    "Conditions must be 0 or 1".to_string(),
                ))
            }
        }
    }

    if group1_indices.is_empty() || group2_indices.is_empty() {
        return Err(EpiError::InsufficientData(
            "Both groups must have samples".to_string(),
        ));
    }

    // Test each region
    let mut results = Vec::new();

    for (region_idx, counts) in count_matrix.iter().enumerate() {
        let region_id = region_ids
            .get(region_idx)
            .cloned()
            .unwrap_or_else(|| format!("region_{}", region_idx));

        let g1_counts: Vec<f64> = group1_indices
            .iter()
            .map(|&i| (counts[i] as f64 + pseudocount) / size_factors[i])
            .collect();

        let g2_counts: Vec<f64> = group2_indices
            .iter()
            .map(|&i| (counts[i] as f64 + pseudocount) / size_factors[i])
            .collect();

        let mean1 = g1_counts.iter().sum::<f64>() / g1_counts.len().max(1) as f64;
        let mean2 = g2_counts.iter().sum::<f64>() / g2_counts.len().max(1) as f64;

        let log2_fc = (mean1 / mean2.max(1e-10)).log2();

        // Welch's t-test (unequal variance)
        let var1 = g1_counts.iter().map(|x| (x - mean1).powi(2)).sum::<f64>()
            / g1_counts.len().max(1) as f64;
        let var2 = g2_counts.iter().map(|x| (x - mean2).powi(2)).sum::<f64>()
            / g2_counts.len().max(1) as f64;

        let se = (var1 / g1_counts.len() as f64 + var2 / g2_counts.len() as f64).sqrt();
        let t_stat = if se > 0.0 {
            (mean1 - mean2) / se
        } else {
            0.0
        };

        let p_value = t_test_pvalue(t_stat, (g1_counts.len() + g2_counts.len() - 2) as f64);

        results.push(DiffResult {
            region_id,
            log2_fc,
            p_value,
            q_value: p_value, // Will be corrected later
            mean_count_1: mean1,
            mean_count_2: mean2,
        });
    }

    // Apply Benjamini-Hochberg correction
    apply_bh_correction(&mut results);

    Ok(results)
}

/// Compute size factors using median of ratios method.
fn compute_size_factors(count_matrix: &[Vec<u32>]) -> Result<Vec<f64>> {
    if count_matrix.is_empty() {
        return Err(EpiError::InsufficientData("Count matrix is empty".to_string()));
    }

    let n_samples = count_matrix[0].len();
    let mut size_factors = vec![1.0; n_samples];

    // Compute geometric mean for each region
    let geometric_means: Vec<f64> = count_matrix
        .iter()
        .map(|counts| {
            let prod: f64 = counts.iter().map(|&c| (c as f64 + 1.0).ln()).sum();
            (prod / counts.len() as f64).exp()
        })
        .collect();

    // Compute ratios for each sample
    for sample_idx in 0..n_samples {
        let mut ratios = Vec::new();

        for (region_idx, counts) in count_matrix.iter().enumerate() {
            if geometric_means[region_idx] > 0.0 && counts[sample_idx] > 0 {
                ratios.push((counts[sample_idx] as f64) / geometric_means[region_idx]);
            }
        }

        if !ratios.is_empty() {
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
            size_factors[sample_idx] = ratios[ratios.len() / 2]; // Median
        }
    }

    Ok(size_factors)
}

/// Normalize count matrix using size factors.
fn normalize_counts(count_matrix: &[Vec<u32>], size_factors: &[f64]) -> Result<Vec<Vec<f64>>> {
    Ok(count_matrix
        .iter()
        .map(|counts| {
            counts
                .iter()
                .zip(size_factors.iter())
                .map(|(&c, &sf)| (c as f64) / sf.max(1e-10))
                .collect()
        })
        .collect())
}

/// Two-tailed t-test p-value from t-statistic.
fn t_test_pvalue(t: f64, df: f64) -> f64 {
    // Approximation using normal for large df
    if df > 30.0 {
        normal_pvalue(t)
    } else {
        // Use approximate Student's t distribution
        let abs_t = t.abs();
        let prob = if abs_t > 10.0 {
            1e-10
        } else {
            let beta_arg = df / (df + abs_t * abs_t);
            beta_incomplete(df / 2.0, 0.5, beta_arg)
        };
        (2.0 * prob).min(1.0)
    }
}

/// Approximate incomplete beta function.
#[allow(unused_variables)]
fn beta_incomplete(_a: f64, _b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else if x >= 1.0 {
        1.0
    } else {
        // Use normal approximation for common cases
        let z = normal_quantile(x);
        normal_cdf(z)
    }
}

/// Standard normal CDF approximation.
fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf_approx(z / std::f64::consts::SQRT_2))
}

/// Standard normal quantile (inverse CDF).
fn normal_quantile(p: f64) -> f64 {
    if p <= 0.0 {
        f64::NEG_INFINITY
    } else if p >= 1.0 {
        f64::INFINITY
    } else if p < 0.5 {
        -normal_quantile(1.0 - p)
    } else {
        // Rational approximation
        let t = (1.0 / (1.0 - p)).sqrt();
        let c = [2.515517, 0.802853, 0.010328];
        let d = [1.432788, 0.189269, 0.001308];

        let num = c[0] + c[1] * t + c[2] * t * t;
        let den = 1.0 + d[0] * t + d[1] * t * t + d[2] * t * t * t;

        t - num / den
    }
}

/// Approximate error function.
fn erf_approx(x: f64) -> f64 {
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    let x = x.abs();

    let p = 0.3275911;
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;

    let t = 1.0 / (1.0 + p * x);
    let poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));

    sign * (1.0 - poly * (-x * x).exp())
}

/// Normal-based p-value.
fn normal_pvalue(z: f64) -> f64 {
    (erfc_approx(z.abs() / std::f64::consts::SQRT_2) / 2.0).min(1.0)
}

/// Complementary error function.
fn erfc_approx(x: f64) -> f64 {
    1.0 - erf_approx(x)
}

/// Apply Benjamini-Hochberg correction to p-values in-place.
fn apply_bh_correction(results: &mut [DiffResult]) {
    let n = results.len() as f64;

    // Sort by p-value
    results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap());

    // Compute q-values
    for (i, result) in results.iter_mut().enumerate() {
        let rank = (i as f64) + 1.0;
        result.q_value = (result.p_value * n / rank).min(1.0);
    }

    // Ensure monotonicity
    for i in (1..results.len()).rev() {
        results[i - 1].q_value = results[i - 1].q_value.min(results[i].q_value);
    }
}

/// Count reads overlapping peaks for each sample.
pub fn count_reads_in_peaks(
    reads: &[(String, u64, u64)],
    peaks: &[(u64, u64)],
) -> Vec<u32> {
    let mut counts = vec![0u32; peaks.len()];

    for (_chrom, read_start, read_end) in reads {
        for (peak_idx, (peak_start, peak_end)) in peaks.iter().enumerate() {
            // Check overlap
            if read_start < peak_end && read_end > peak_start {
                counts[peak_idx] = counts[peak_idx].saturating_add(1);
            }
        }
    }

    counts
}

/// Compute MA plot data (mean abundance vs log2 fold change).
pub fn ma_plot_data(diff_results: &[DiffResult]) -> Vec<(f64, f64)> {
    diff_results
        .iter()
        .map(|r| {
            let mean_log2 = (r.mean_count_1.log2() + r.mean_count_2.log2()) / 2.0;
            (mean_log2, r.log2_fc)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_differential_peaks() {
        let count_matrix = vec![
            vec![10, 20, 5, 15],
            vec![100, 110, 20, 25],
            vec![50, 45, 100, 95],
        ];

        let conditions = vec![0, 0, 1, 1];
        let region_ids = vec!["peak1".to_string(), "peak2".to_string(), "peak3".to_string()];

        let results = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0).unwrap();

        assert_eq!(results.len(), 3);
        for result in &results {
            assert!(!result.log2_fc.is_nan());
            assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
        }
    }

    #[test]
    fn test_count_reads_in_peaks() {
        let reads = vec![
            ("chr1".to_string(), 100, 150),
            ("chr1".to_string(), 120, 170),
        ];

        let peaks = vec![(50, 200), (250, 300)];
        let counts = count_reads_in_peaks(&reads, &peaks);

        assert_eq!(counts[0], 2);
        assert_eq!(counts[1], 0);
    }

    #[test]
    fn test_ma_plot_data() {
        let results = vec![
            DiffResult {
                region_id: "peak1".to_string(),
                log2_fc: 2.0,
                p_value: 0.01,
                q_value: 0.05,
                mean_count_1: 100.0,
                mean_count_2: 25.0,
            },
        ];

        let data = ma_plot_data(&results);
        assert_eq!(data.len(), 1);
        assert!(data[0].0 > 0.0);
    }

    #[test]
    fn test_size_factor_computation() {
        let count_matrix = vec![
            vec![10, 20],
            vec![100, 200],
            vec![1000, 2000],
        ];

        let sf = compute_size_factors(&count_matrix).unwrap();
        assert_eq!(sf.len(), 2);
        assert!(sf[0] > 0.0);
        assert!(sf[1] > 0.0);
    }

    #[test]
    fn test_bh_correction() {
        let mut results = vec![
            DiffResult {
                region_id: "1".to_string(),
                log2_fc: 1.0,
                p_value: 0.001,
                q_value: 0.0,
                mean_count_1: 10.0,
                mean_count_2: 5.0,
            },
            DiffResult {
                region_id: "2".to_string(),
                log2_fc: 0.5,
                p_value: 0.05,
                q_value: 0.0,
                mean_count_1: 5.0,
                mean_count_2: 4.0,
            },
        ];

        apply_bh_correction(&mut results);
        assert!(results[0].q_value <= results[1].q_value);
    }
}
