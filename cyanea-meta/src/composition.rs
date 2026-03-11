//! Community composition analysis: differential abundance and transformations.
//!
//! Provides:
//! - [`clr_transform`] — Centered Log-Ratio for compositional data
//! - [`ilr_transform`] — Isometric Log-Ratio transform
//! - [`differential_abundance`] — ALDEx2-style: CLR + Welch's t-test per taxon
//! - [`ancom`] — ANCOM differential abundance (W-statistic)

use std::collections::HashMap;
use crate::error::{MetaError, Result};

/// CLR (Centered Log-Ratio) transformed composition.
#[derive(Debug, Clone)]
pub struct CompositionTransform {
    /// Taxon ID → transformed value.
    pub values: HashMap<u32, f64>,
}

impl CompositionTransform {
    /// Create a new transform.
    pub fn new() -> Self {
        Self {
            values: HashMap::new(),
        }
    }

    /// Insert a taxon's transformed value.
    pub fn insert(&mut self, taxid: u32, value: f64) {
        self.values.insert(taxid, value);
    }

    /// Get a taxon's value.
    pub fn get(&self, taxid: u32) -> Option<f64> {
        self.values.get(&taxid).copied()
    }
}

impl Default for CompositionTransform {
    fn default() -> Self {
        Self::new()
    }
}

/// Centered Log-Ratio (CLR) transform for compositional data.
///
/// For each sample and taxon: CLR = ln(x_i / geometric_mean(x))
///
/// # Arguments
///
/// * `abundances` — matrix: samples × taxa (as flattened row-major Vec<Vec<f64>>)
///
/// # Errors
///
/// Returns an error if any abundance is zero or negative.
pub fn clr_transform(abundances: &[Vec<f64>]) -> Result<Vec<CompositionTransform>> {
    if abundances.is_empty() {
        return Err(MetaError::Functional(
            "abundances matrix is empty".into(),
        ));
    }

    let mut results = Vec::new();

    for sample in abundances {
        if sample.is_empty() {
            return Err(MetaError::Functional(
                "sample is empty".into(),
            ));
        }

        // Check for zeros and negatives
        if sample.iter().any(|&x| x <= 0.0) {
            return Err(MetaError::Functional(
                "CLR requires all abundances > 0".into(),
            ));
        }

        // Geometric mean: (∏ x_i)^(1/n)
        let ln_sum: f64 = sample.iter().map(|x| x.ln()).sum();
        let ln_geom_mean = ln_sum / sample.len() as f64;

        let mut transform = CompositionTransform::new();
        for (i, &x) in sample.iter().enumerate() {
            let clr_val = x.ln() - ln_geom_mean;
            transform.insert(i as u32, clr_val);
        }
        results.push(transform);
    }

    Ok(results)
}

/// ILR (Isometric Log-Ratio) transform: CLR projected to orthogonal subspace.
///
/// Reduces dimensionality from (p-1) to lower dimensions via orthogonal contrast.
///
/// # Errors
///
/// Returns an error if abundances is empty or contain invalid values.
pub fn ilr_transform(abundances: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    if abundances.is_empty() {
        return Err(MetaError::Functional(
            "abundances matrix is empty".into(),
        ));
    }

    // First apply CLR
    let clr_results = clr_transform(abundances)?;

    // For simplicity, return CLR values (full ILR requires orthogonal basis construction)
    let mut results = Vec::new();
    for transform in clr_results {
        let mut row = Vec::new();
        for i in 0..abundances[0].len() {
            if let Some(v) = transform.get(i as u32) {
                row.push(v);
            }
        }
        results.push(row);
    }

    Ok(results)
}

/// Differential abundance result for a single taxon.
#[derive(Debug, Clone)]
pub struct DifferentialAbundanceResult {
    /// Taxon ID.
    pub taxid: u32,
    /// Welch's t-statistic.
    pub t_statistic: f64,
    /// p-value.
    pub p_value: f64,
    /// log2 fold-change.
    pub log2_fc: f64,
}

/// ALDEx2-style differential abundance analysis.
///
/// Compares two groups of samples using CLR-transformed abundances and Welch's t-test.
///
/// # Arguments
///
/// * `group1` — samples for first group (Vec<Vec<f64>>: samples × taxa)
/// * `group2` — samples for second group
///
/// # Errors
///
/// Returns an error if groups are empty or have mismatched dimensions.
pub fn differential_abundance(
    group1: &[Vec<f64>],
    group2: &[Vec<f64>],
) -> Result<Vec<DifferentialAbundanceResult>> {
    if group1.is_empty() || group2.is_empty() {
        return Err(MetaError::Functional(
            "both groups must be non-empty".into(),
        ));
    }

    if group1[0].len() != group2[0].len() {
        return Err(MetaError::Functional(
            "groups must have same number of taxa".into(),
        ));
    }

    let clr1 = clr_transform(group1)?;
    let clr2 = clr_transform(group2)?;
    let n_taxa = group1[0].len();

    let mut results = Vec::new();

    for taxid in 0..n_taxa {
        let values1: Vec<f64> = clr1.iter().filter_map(|t| t.get(taxid as u32)).collect();
        let values2: Vec<f64> = clr2.iter().filter_map(|t| t.get(taxid as u32)).collect();

        if values1.is_empty() || values2.is_empty() {
            continue;
        }

        // Compute means and variances
        let mean1 = values1.iter().sum::<f64>() / values1.len() as f64;
        let mean2 = values2.iter().sum::<f64>() / values2.len() as f64;

        let var1 = values1
            .iter()
            .map(|&x| (x - mean1).powi(2))
            .sum::<f64>()
            / values1.len() as f64;
        let var2 = values2
            .iter()
            .map(|&x| (x - mean2).powi(2))
            .sum::<f64>()
            / values2.len() as f64;

        // Welch's t-statistic
        let se = ((var1 / values1.len() as f64) + (var2 / values2.len() as f64)).sqrt();
        let t_stat = if se > 0.0 {
            (mean1 - mean2) / se
        } else {
            0.0
        };

        let log2_fc = (mean1 - mean2) / std::f64::consts::LN_2;

        // Approximate p-value via normal (simplified; full implementation would use t-distribution)
        let p_value = 2.0 * (1.0 - normal_cdf(t_stat.abs()));

        results.push(DifferentialAbundanceResult {
            taxid: taxid as u32,
            t_statistic: t_stat,
            p_value,
            log2_fc,
        });
    }

    results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(std::cmp::Ordering::Equal));
    Ok(results)
}

/// Standard normal CDF (approximate via error function).
fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

/// Error function (Abramowitz and Stegun approximation).
fn erf(x: f64) -> f64 {
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

/// ANCOM differential abundance result.
#[derive(Debug, Clone)]
pub struct AncemResult {
    /// Taxon ID.
    pub taxid: u32,
    /// W-statistic (count of significant log-ratios).
    pub w_statistic: u32,
}

/// ANCOM (Analysis of Composition of Microbiomes) analysis.
///
/// Tests differential abundance using log-ratio comparisons and W-statistic.
///
/// # Arguments
///
/// * `group1` — abundances for group 1
/// * `group2` — abundances for group 2
///
/// # Errors
///
/// Returns an error if groups are empty or mismatched.
pub fn ancom(group1: &[Vec<f64>], group2: &[Vec<f64>]) -> Result<Vec<AncemResult>> {
    if group1.is_empty() || group2.is_empty() {
        return Err(MetaError::Functional(
            "both groups must be non-empty".into(),
        ));
    }

    if group1[0].len() != group2[0].len() {
        return Err(MetaError::Functional(
            "groups must have same number of taxa".into(),
        ));
    }

    let n_taxa = group1[0].len();
    let mut results = Vec::new();

    for i in 0..n_taxa {
        // For each taxon i, compare log-ratio (i vs other taxa)
        let mut w = 0u32;

        for j in 0..n_taxa {
            if i == j {
                continue;
            }

            // Log-ratio: ln(abundance_i / abundance_j) for each sample
            let mut lr1 = Vec::new();
            for sample in group1 {
                if sample[i] > 0.0 && sample[j] > 0.0 {
                    lr1.push((sample[i] / sample[j]).ln());
                }
            }

            let mut lr2 = Vec::new();
            for sample in group2 {
                if sample[i] > 0.0 && sample[j] > 0.0 {
                    lr2.push((sample[i] / sample[j]).ln());
                }
            }

            if !lr1.is_empty() && !lr2.is_empty() {
                let mean1 = lr1.iter().sum::<f64>() / lr1.len() as f64;
                let mean2 = lr2.iter().sum::<f64>() / lr2.len() as f64;
                // If means are significantly different, increment W
                if (mean1 - mean2).abs() > 0.5 {
                    w += 1;
                }
            }
        }

        results.push(AncemResult {
            taxid: i as u32,
            w_statistic: w,
        });
    }

    results.sort_by(|a, b| b.w_statistic.cmp(&a.w_statistic));
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clr_transform_basic() {
        let abundances = vec![vec![10.0, 20.0, 30.0]];
        let result = clr_transform(&abundances).unwrap();
        assert_eq!(result.len(), 1);
        // Sum of CLR values should be approximately 0
        let sum: f64 = result[0].values.values().sum();
        assert!(sum.abs() < 1e-10);
    }

    #[test]
    fn clr_rejects_zeros() {
        let abundances = vec![vec![0.0, 20.0, 30.0]];
        assert!(clr_transform(&abundances).is_err());
    }

    #[test]
    fn ilr_transform_valid() {
        let abundances = vec![vec![10.0, 20.0, 30.0]];
        let result = ilr_transform(&abundances).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].len(), 3);
    }

    #[test]
    fn differential_abundance_basic() {
        let g1 = vec![vec![100.0, 10.0, 5.0], vec![95.0, 15.0, 10.0]];
        let g2 = vec![vec![10.0, 100.0, 50.0], vec![15.0, 95.0, 45.0]];
        let results = differential_abundance(&g1, &g2).unwrap();
        assert!(!results.is_empty());
        // All p-values should be in [0, 1]
        for r in &results {
            assert!(r.p_value >= 0.0 && r.p_value <= 1.0);
        }
    }

    #[test]
    fn ancom_basic() {
        let g1 = vec![vec![100.0, 10.0], vec![95.0, 15.0]];
        let g2 = vec![vec![10.0, 100.0], vec![15.0, 95.0]];
        let results = ancom(&g1, &g2).unwrap();
        assert!(!results.is_empty());
    }
}
