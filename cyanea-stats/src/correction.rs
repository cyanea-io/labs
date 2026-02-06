//! Multiple testing correction.
//!
//! When running many hypothesis tests simultaneously, p-values must be
//! adjusted to control the family-wise error rate or false discovery rate.

use cyanea_core::{CyaneaError, Result};

/// Multiple testing correction method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrectionMethod {
    /// Bonferroni correction — controls family-wise error rate (FWER).
    Bonferroni,
    /// Benjamini-Hochberg procedure — controls false discovery rate (FDR).
    BenjaminiHochberg,
}

/// Apply a multiple testing correction to `p_values`.
///
/// Returns a new `Vec<f64>` of adjusted p-values in the same order as the
/// input.
pub fn correct(p_values: &[f64], method: CorrectionMethod) -> Result<Vec<f64>> {
    match method {
        CorrectionMethod::Bonferroni => bonferroni(p_values),
        CorrectionMethod::BenjaminiHochberg => benjamini_hochberg(p_values),
    }
}

/// Bonferroni correction: `p_adj = min(p * n, 1.0)`.
pub fn bonferroni(p_values: &[f64]) -> Result<Vec<f64>> {
    validate_p_values(p_values)?;
    let n = p_values.len() as f64;
    Ok(p_values.iter().map(|&p| (p * n).min(1.0)).collect())
}

/// Benjamini-Hochberg procedure for controlling the false discovery rate.
///
/// Sorts p-values, adjusts as `p * n / rank`, enforces monotonicity
/// from right to left, and clamps to [0, 1].
pub fn benjamini_hochberg(p_values: &[f64]) -> Result<Vec<f64>> {
    validate_p_values(p_values)?;
    let n = p_values.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Sort indices by p-value.
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| p_values[a].total_cmp(&p_values[b]));

    let n_f = n as f64;
    let mut adjusted = vec![0.0; n];

    // Compute adjusted p-values and enforce monotonicity (right to left).
    let mut prev = f64::INFINITY;
    for i in (0..n).rev() {
        let rank = (i + 1) as f64;
        let adj = (p_values[indices[i]] * n_f / rank).min(1.0);
        let adj = adj.min(prev);
        adjusted[indices[i]] = adj;
        prev = adj;
    }

    Ok(adjusted)
}

fn validate_p_values(p_values: &[f64]) -> Result<()> {
    for (i, &p) in p_values.iter().enumerate() {
        if !(0.0..=1.0).contains(&p) {
            return Err(CyaneaError::InvalidInput(format!(
                "p-value at index {} is out of range [0, 1]: {}",
                i, p,
            )));
        }
    }
    Ok(())
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn bonferroni_basic() {
        let p = [0.01, 0.04, 0.03, 0.005];
        let adj = bonferroni(&p).unwrap();
        assert!((adj[0] - 0.04).abs() < TOL);
        assert!((adj[1] - 0.16).abs() < TOL);
        assert!((adj[2] - 0.12).abs() < TOL);
        assert!((adj[3] - 0.02).abs() < TOL);
    }

    #[test]
    fn bonferroni_clamp() {
        let p = [0.5, 0.8];
        let adj = bonferroni(&p).unwrap();
        assert!((adj[0] - 1.0).abs() < TOL);
        assert!((adj[1] - 1.0).abs() < TOL);
    }

    #[test]
    fn bh_known() {
        // Classic BH example
        let p = [0.01, 0.04, 0.03, 0.005];
        let adj = benjamini_hochberg(&p).unwrap();
        // Sorted: 0.005(idx3), 0.01(idx0), 0.03(idx2), 0.04(idx1)
        // Ranks:    1            2            3            4
        // Raw adj: 0.005*4/1=0.02, 0.01*4/2=0.02, 0.03*4/3=0.04, 0.04*4/4=0.04
        // Monotonicity (R-to-L): 0.04, 0.04, 0.02, 0.02
        assert!((adj[3] - 0.02).abs() < TOL);
        assert!((adj[0] - 0.02).abs() < TOL);
        assert!((adj[2] - 0.04).abs() < TOL);
        assert!((adj[1] - 0.04).abs() < TOL);
    }

    #[test]
    fn bh_monotonicity() {
        // Ensure sorted adjusted p-values are non-decreasing.
        let p = [0.1, 0.001, 0.05, 0.01, 0.5];
        let adj = benjamini_hochberg(&p).unwrap();
        let mut sorted_adj: Vec<(f64, f64)> = p.iter().copied().zip(adj.iter().copied()).collect();
        sorted_adj.sort_by(|a, b| a.0.total_cmp(&b.0));
        for w in sorted_adj.windows(2) {
            assert!(
                w[1].1 >= w[0].1 - TOL,
                "monotonicity violated: {} > {}",
                w[0].1,
                w[1].1
            );
        }
    }

    #[test]
    fn bh_clamp() {
        let p = [0.9, 0.95];
        let adj = benjamini_hochberg(&p).unwrap();
        assert!(adj[0] <= 1.0);
        assert!(adj[1] <= 1.0);
    }

    #[test]
    fn correction_empty() {
        assert_eq!(bonferroni(&[]).unwrap(), Vec::<f64>::new());
        assert_eq!(benjamini_hochberg(&[]).unwrap(), Vec::<f64>::new());
    }

    #[test]
    fn correction_single() {
        assert!((bonferroni(&[0.05]).unwrap()[0] - 0.05).abs() < TOL);
        assert!((benjamini_hochberg(&[0.05]).unwrap()[0] - 0.05).abs() < TOL);
    }

    #[test]
    fn correction_invalid_p() {
        assert!(bonferroni(&[0.5, 1.5]).is_err());
        assert!(benjamini_hochberg(&[-0.1, 0.5]).is_err());
    }
}
