//! Correlation analysis.
//!
//! Provides Pearson and Spearman correlation coefficients, and a
//! [`CorrelationMatrix`] for pairwise analysis of multiple variables.

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::rank::{rank, RankMethod};

/// Pearson product-moment correlation coefficient between `x` and `y`.
///
/// Returns 0.0 if either series is constant (zero variance).
pub fn pearson(x: &[f64], y: &[f64]) -> Result<f64> {
    validate_paired(x, y)?;

    let n = x.len() as f64;
    let mean_x: f64 = x.iter().sum::<f64>() / n;
    let mean_y: f64 = y.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;
    for (xi, yi) in x.iter().zip(y.iter()) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 {
        return Ok(0.0);
    }
    Ok(cov / denom)
}

/// Spearman rank correlation coefficient between `x` and `y`.
///
/// Ranks both series with [`RankMethod::Average`], then computes Pearson
/// correlation on the ranks.
pub fn spearman(x: &[f64], y: &[f64]) -> Result<f64> {
    validate_paired(x, y)?;
    let rx = rank(x, RankMethod::Average);
    let ry = rank(y, RankMethod::Average);
    pearson(&rx, &ry)
}

fn validate_paired(x: &[f64], y: &[f64]) -> Result<()> {
    if x.len() != y.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "correlation: x and y must have the same length ({} vs {})",
            x.len(),
            y.len(),
        )));
    }
    if x.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "correlation: need at least 2 observations".into(),
        ));
    }
    Ok(())
}

// ── Correlation matrix ─────────────────────────────────────────────────────

/// Pairwise Pearson correlation matrix for a set of variables.
#[derive(Debug, Clone)]
pub struct CorrelationMatrix {
    /// Flat storage (row-major, n×n).
    data: Vec<f64>,
    /// Number of variables.
    size: usize,
    /// Optional variable labels.
    labels: Option<Vec<String>>,
}

impl CorrelationMatrix {
    /// Build a correlation matrix from rows of observations.
    ///
    /// Each inner slice is one variable's observations (all must have the same
    /// length and at least 2 elements).
    pub fn from_rows(rows: &[&[f64]]) -> Result<Self> {
        Self::build(rows, None)
    }

    /// Build a labeled correlation matrix.
    pub fn from_rows_labeled(rows: &[&[f64]], labels: &[&str]) -> Result<Self> {
        if labels.len() != rows.len() {
            return Err(CyaneaError::InvalidInput(
                "CorrelationMatrix: labels length must match rows length".into(),
            ));
        }
        let labels: Vec<String> = labels.iter().map(|s| s.to_string()).collect();
        Self::build(rows, Some(labels))
    }

    fn build(rows: &[&[f64]], labels: Option<Vec<String>>) -> Result<Self> {
        if rows.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "CorrelationMatrix: need at least one variable".into(),
            ));
        }
        let obs_len = rows[0].len();
        for (i, row) in rows.iter().enumerate() {
            if row.len() != obs_len {
                return Err(CyaneaError::InvalidInput(format!(
                    "CorrelationMatrix: row {} has {} observations, expected {}",
                    i,
                    row.len(),
                    obs_len,
                )));
            }
        }

        let n = rows.len();
        #[cfg(feature = "parallel")]
        let data = {
            use rayon::prelude::*;
            let upper: Vec<Vec<(usize, f64)>> = (0..n)
                .into_par_iter()
                .map(|i| {
                    ((i + 1)..n)
                        .map(|j| {
                            let r = pearson(rows[i], rows[j]).unwrap();
                            (j, r)
                        })
                        .collect()
                })
                .collect();
            let mut data = vec![0.0; n * n];
            for i in 0..n {
                data[i * n + i] = 1.0;
                for &(j, r) in &upper[i] {
                    data[i * n + j] = r;
                    data[j * n + i] = r;
                }
            }
            data
        };
        #[cfg(not(feature = "parallel"))]
        let data = {
            let mut data = vec![0.0; n * n];
            for i in 0..n {
                data[i * n + i] = 1.0;
                for j in (i + 1)..n {
                    let r = pearson(rows[i], rows[j])?;
                    data[i * n + j] = r;
                    data[j * n + i] = r;
                }
            }
            data
        };

        Ok(Self {
            data,
            size: n,
            labels,
        })
    }

    /// Get the correlation between variable `i` and variable `j`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.size + j]
    }

    /// Number of variables.
    pub fn n(&self) -> usize {
        self.size
    }

    /// Variable labels, if provided.
    pub fn labels(&self) -> Option<&[String]> {
        self.labels.as_deref()
    }
}

impl Summarizable for CorrelationMatrix {
    fn summary(&self) -> String {
        format!("CorrelationMatrix: {}x{}", self.size, self.size)
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn pearson_perfect_positive() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((pearson(&x, &y).unwrap() - 1.0).abs() < TOL);
    }

    #[test]
    fn pearson_perfect_negative() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [10.0, 8.0, 6.0, 4.0, 2.0];
        assert!((pearson(&x, &y).unwrap() - (-1.0)).abs() < TOL);
    }

    #[test]
    fn pearson_zero_correlation() {
        // Orthogonal pattern
        let x = [1.0, 0.0, -1.0, 0.0];
        let y = [0.0, 1.0, 0.0, -1.0];
        assert!((pearson(&x, &y).unwrap()).abs() < TOL);
    }

    #[test]
    fn pearson_constant_series() {
        let x = [3.0, 3.0, 3.0];
        let y = [1.0, 2.0, 3.0];
        assert!((pearson(&x, &y).unwrap()).abs() < TOL);
    }

    #[test]
    fn pearson_length_mismatch() {
        assert!(pearson(&[1.0, 2.0], &[1.0]).is_err());
    }

    #[test]
    fn pearson_too_short() {
        assert!(pearson(&[1.0], &[2.0]).is_err());
    }

    #[test]
    fn spearman_monotonic() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 8.0, 27.0, 64.0, 125.0]; // x^3 — monotonically increasing
        assert!((spearman(&x, &y).unwrap() - 1.0).abs() < TOL);
    }

    #[test]
    fn spearman_reverse() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [5.0, 4.0, 3.0, 2.0, 1.0];
        assert!((spearman(&x, &y).unwrap() - (-1.0)).abs() < TOL);
    }

    #[test]
    fn correlation_matrix_diagonal() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = [7.0, 8.0, 9.0];
        let cm = CorrelationMatrix::from_rows(&[&a[..], &b[..], &c[..]]).unwrap();
        assert_eq!(cm.n(), 3);
        assert!((cm.get(0, 0) - 1.0).abs() < TOL);
        assert!((cm.get(1, 1) - 1.0).abs() < TOL);
        assert!((cm.get(2, 2) - 1.0).abs() < TOL);
    }

    #[test]
    fn correlation_matrix_symmetric() {
        let a = [1.0, 2.0, 3.0, 4.0];
        let b = [4.0, 3.0, 2.0, 1.0];
        let cm = CorrelationMatrix::from_rows(&[&a[..], &b[..]]).unwrap();
        assert!((cm.get(0, 1) - cm.get(1, 0)).abs() < TOL);
    }

    #[test]
    fn correlation_matrix_summary() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let cm = CorrelationMatrix::from_rows(&[&a[..], &b[..]]).unwrap();
        assert_eq!(cm.summary(), "CorrelationMatrix: 2x2");
    }
}
