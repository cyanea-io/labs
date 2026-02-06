//! Feature normalization for ML pipelines.

use cyanea_core::{CyaneaError, Result};

/// Scale values in-place to `[0, 1]` using min-max normalization.
///
/// Constant data (max == min) is set to 0.0.
pub fn min_max(data: &mut [f64]) -> Result<()> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    let min = data.iter().copied().fold(f64::INFINITY, f64::min);
    let max = data.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    if range == 0.0 {
        data.iter_mut().for_each(|v| *v = 0.0);
    } else {
        data.iter_mut().for_each(|v| *v = (*v - min) / range);
    }
    Ok(())
}

/// Standardize values in-place to zero mean and unit variance (z-score).
///
/// Constant data (std == 0) is set to 0.0.
pub fn z_score(data: &mut [f64]) -> Result<()> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    let n = data.len() as f64;
    let mean = data.iter().sum::<f64>() / n;
    let var = data.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    let std = var.sqrt();
    if std == 0.0 {
        data.iter_mut().for_each(|v| *v = 0.0);
    } else {
        data.iter_mut().for_each(|v| *v = (*v - mean) / std);
    }
    Ok(())
}

/// Normalize values in-place to unit L2 norm.
///
/// Zero vectors are left unchanged.
pub fn l2_normalize(data: &mut [f64]) -> Result<()> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    let norm: f64 = data.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm > 0.0 {
        data.iter_mut().for_each(|v| *v /= norm);
    }
    Ok(())
}

// --- Column-wise variants for flat row-major matrices ---

fn validate_matrix(data: &[f64], n_cols: usize) -> Result<()> {
    if n_cols == 0 {
        return Err(CyaneaError::InvalidInput("n_cols must be > 0".into()));
    }
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    if data.len() % n_cols != 0 {
        return Err(CyaneaError::InvalidInput(format!(
            "data length {} not divisible by n_cols {}",
            data.len(),
            n_cols
        )));
    }
    Ok(())
}

/// Min-max normalize each column of a flat row-major matrix in-place.
pub fn min_max_columns(data: &mut [f64], n_cols: usize) -> Result<()> {
    validate_matrix(data, n_cols)?;
    let n_rows = data.len() / n_cols;
    for col in 0..n_cols {
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        for row in 0..n_rows {
            let v = data[row * n_cols + col];
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
        }
        let range = max - min;
        for row in 0..n_rows {
            let idx = row * n_cols + col;
            if range == 0.0 {
                data[idx] = 0.0;
            } else {
                data[idx] = (data[idx] - min) / range;
            }
        }
    }
    Ok(())
}

/// Z-score normalize each column of a flat row-major matrix in-place.
pub fn z_score_columns(data: &mut [f64], n_cols: usize) -> Result<()> {
    validate_matrix(data, n_cols)?;
    let n_rows = data.len() / n_cols;
    let nf = n_rows as f64;
    for col in 0..n_cols {
        let mean: f64 = (0..n_rows)
            .map(|r| data[r * n_cols + col])
            .sum::<f64>()
            / nf;
        let var: f64 = (0..n_rows)
            .map(|r| (data[r * n_cols + col] - mean).powi(2))
            .sum::<f64>()
            / nf;
        let std = var.sqrt();
        for row in 0..n_rows {
            let idx = row * n_cols + col;
            if std == 0.0 {
                data[idx] = 0.0;
            } else {
                data[idx] = (data[idx] - mean) / std;
            }
        }
    }
    Ok(())
}

/// L2-normalize each column of a flat row-major matrix in-place.
pub fn l2_normalize_columns(data: &mut [f64], n_cols: usize) -> Result<()> {
    validate_matrix(data, n_cols)?;
    let n_rows = data.len() / n_cols;
    for col in 0..n_cols {
        let norm: f64 = (0..n_rows)
            .map(|r| data[r * n_cols + col].powi(2))
            .sum::<f64>()
            .sqrt();
        if norm > 0.0 {
            for row in 0..n_rows {
                data[row * n_cols + col] /= norm;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_max_basic() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        min_max(&mut data).unwrap();
        assert!((data[0] - 0.0).abs() < 1e-12);
        assert!((data[4] - 1.0).abs() < 1e-12);
        assert!((data[2] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn min_max_constant() {
        let mut data = vec![3.0, 3.0, 3.0];
        min_max(&mut data).unwrap();
        assert!(data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn min_max_empty_error() {
        assert!(min_max(&mut []).is_err());
    }

    #[test]
    fn z_score_basic() {
        let mut data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        z_score(&mut data).unwrap();
        let mean: f64 = data.iter().sum::<f64>() / data.len() as f64;
        assert!(mean.abs() < 1e-10);
    }

    #[test]
    fn z_score_constant() {
        let mut data = vec![5.0, 5.0, 5.0];
        z_score(&mut data).unwrap();
        assert!(data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn z_score_empty_error() {
        assert!(z_score(&mut []).is_err());
    }

    #[test]
    fn l2_normalize_basic() {
        let mut data = vec![3.0, 4.0];
        l2_normalize(&mut data).unwrap();
        let norm: f64 = data.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-12);
    }

    #[test]
    fn l2_normalize_zero_vector() {
        let mut data = vec![0.0, 0.0, 0.0];
        l2_normalize(&mut data).unwrap();
        assert!(data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn l2_normalize_empty_error() {
        assert!(l2_normalize(&mut []).is_err());
    }

    #[test]
    fn min_max_columns_basic() {
        // 3x2 matrix: [[1,10], [2,20], [3,30]]
        let mut data = vec![1.0, 10.0, 2.0, 20.0, 3.0, 30.0];
        min_max_columns(&mut data, 2).unwrap();
        // col 0: [0.0, 0.5, 1.0], col 1: [0.0, 0.5, 1.0]
        assert!((data[0] - 0.0).abs() < 1e-12);
        assert!((data[2] - 0.5).abs() < 1e-12);
        assert!((data[4] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn z_score_columns_basic() {
        let mut data = vec![1.0, 100.0, 3.0, 300.0, 5.0, 500.0];
        z_score_columns(&mut data, 2).unwrap();
        // each column should have mean ~0
        let col0_mean = (data[0] + data[2] + data[4]) / 3.0;
        let col1_mean = (data[1] + data[3] + data[5]) / 3.0;
        assert!(col0_mean.abs() < 1e-10);
        assert!(col1_mean.abs() < 1e-10);
    }

    #[test]
    fn l2_normalize_columns_basic() {
        let mut data = vec![3.0, 0.0, 4.0, 0.0];
        l2_normalize_columns(&mut data, 2).unwrap();
        // col 0: [3/5, 4/5], col 1: [0, 0] (zero col unchanged)
        assert!((data[0] - 0.6).abs() < 1e-12);
        assert!((data[2] - 0.8).abs() < 1e-12);
        assert!((data[1] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn column_validate_bad_ncols() {
        assert!(min_max_columns(&mut [1.0, 2.0], 0).is_err());
    }

    #[test]
    fn column_validate_not_divisible() {
        assert!(min_max_columns(&mut [1.0, 2.0, 3.0], 2).is_err());
    }
}
