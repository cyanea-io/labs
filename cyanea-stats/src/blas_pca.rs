//! BLAS-accelerated PCA using ndarray for matrix operations.
//!
//! When the `blas` feature is enabled, this module provides a PCA implementation
//! that uses ndarray's matrix multiply (which can dispatch to BLAS) for the
//! covariance matrix computation and power iteration.

use cyanea_core::Result;
use ndarray::{Array1, Array2};

use super::PcaResult;

/// PCA implementation using ndarray for matrix operations.
pub(super) fn pca_ndarray(
    data: &[&[f64]],
    n_features: usize,
    n_samples: usize,
    n_components: usize,
) -> Result<PcaResult> {
    // Build the data matrix
    let mut flat = Vec::with_capacity(n_samples * n_features);
    for row in data {
        flat.extend_from_slice(row);
    }
    let x = Array2::from_shape_vec((n_samples, n_features), flat)
        .expect("shape already validated");

    // Compute and subtract mean
    let mean_arr = x.mean_axis(ndarray::Axis(0)).unwrap();
    let centered = &x - &mean_arr;

    // Covariance matrix: C = X^T X / (n-1)
    let scale = if n_samples > 1 {
        (n_samples - 1) as f64
    } else {
        1.0
    };
    let cov = centered.t().dot(&centered) / scale;

    let total_variance: f64 = cov.diag().sum();

    // Power iteration with deflation
    let mut cov_work = cov;
    let mut components = Vec::with_capacity(n_components);
    let mut eigenvalues = Vec::with_capacity(n_components);

    for _ in 0..n_components {
        let (eigenvalue, eigenvector) = ndarray_power_iteration(&cov_work, 300);
        eigenvalues.push(eigenvalue);
        components.push(eigenvector.to_vec());

        // Deflate
        let v = eigenvector.view().insert_axis(ndarray::Axis(1));
        cov_work = cov_work - eigenvalue * v.dot(&v.t());
    }

    let explained_variance_ratio: Vec<f64> = eigenvalues
        .iter()
        .map(|&e| if total_variance > 0.0 { e / total_variance } else { 0.0 })
        .collect();

    // Project data onto components
    let transformed: Vec<Vec<f64>> = (0..n_samples)
        .map(|i| {
            let row = centered.row(i);
            components
                .iter()
                .map(|comp| {
                    row.iter()
                        .zip(comp.iter())
                        .map(|(a, b)| a * b)
                        .sum()
                })
                .collect()
        })
        .collect();

    let mean = mean_arr.to_vec();

    Ok(PcaResult {
        components,
        explained_variance: eigenvalues,
        explained_variance_ratio,
        transformed,
        mean,
    })
}

/// Power iteration using ndarray matrix-vector multiply.
fn ndarray_power_iteration(matrix: &Array2<f64>, max_iter: usize) -> (f64, Array1<f64>) {
    let n = matrix.nrows();

    let mut v = Array1::from_vec((0..n).map(|i| 1.0 / ((i + 1) as f64)).collect());
    let norm = v.dot(&v).sqrt();
    if norm > 0.0 {
        v /= norm;
    }

    let mut eigenvalue = 0.0;

    for _ in 0..max_iter {
        let w = matrix.dot(&v);
        let new_eigenvalue = v.dot(&w);
        let wnorm = w.dot(&w).sqrt();
        if wnorm < 1e-15 {
            break;
        }
        let w_normalized = &w / wnorm;
        let diff: f64 = (&v - &w_normalized).mapv(|x| x * x).sum();
        v = w_normalized;
        eigenvalue = new_eigenvalue;
        if diff < 1e-12 {
            break;
        }
    }

    (eigenvalue.max(0.0), v)
}
