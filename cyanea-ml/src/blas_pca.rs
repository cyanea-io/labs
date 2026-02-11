//! BLAS-accelerated PCA using ndarray for matrix operations.
//!
//! When the `blas` feature is enabled, this module provides a PCA implementation
//! that uses ndarray's matrix multiply (which can dispatch to BLAS) for the
//! covariance matrix computation and power iteration, giving significant speedups
//! for large feature counts.

use cyanea_core::Result;
use ndarray::{Array1, Array2};

use super::{PcaConfig, PcaResult};

/// PCA implementation using ndarray for matrix operations.
pub(super) fn pca_ndarray(
    data: &[f64],
    n_features: usize,
    n_samples: usize,
    n_components: usize,
    config: &PcaConfig,
) -> Result<PcaResult> {
    // Build the data matrix (n_samples × n_features)
    let x = Array2::from_shape_vec((n_samples, n_features), data.to_vec())
        .expect("shape already validated");

    // Compute and subtract mean
    let mean_arr = x.mean_axis(ndarray::Axis(0)).unwrap();
    let centered = &x - &mean_arr;

    // Covariance matrix: C = X^T X / (n-1), using ndarray gemm (BLAS when available)
    let scale = (n_samples - 1) as f64;
    let cov = centered.t().dot(&centered) / scale;

    // Power iteration with deflation
    let mut cov_work = cov.clone();
    let mut components = Vec::with_capacity(n_components * n_features);
    let mut eigenvalues = Vec::with_capacity(n_components);

    for _ in 0..n_components {
        let (eigenvalue, eigenvector) =
            ndarray_power_iteration(&cov_work, config.max_iter, config.tolerance);
        components.extend(eigenvector.iter());
        eigenvalues.push(eigenvalue);

        // Deflate: C = C - lambda * v * v^T
        let v = eigenvector.view().insert_axis(ndarray::Axis(1)); // column vector
        cov_work = cov_work - eigenvalue * v.dot(&v.t());
    }

    // Explained variance ratio
    let orig_trace: f64 = cov.diag().sum();
    let explained_variance_ratio = if orig_trace > 0.0 {
        eigenvalues.iter().map(|&ev| ev / orig_trace).collect()
    } else {
        vec![0.0; n_components]
    };

    // Project: transformed = centered * components^T
    let comp_matrix =
        Array2::from_shape_vec((n_components, n_features), components.clone()).unwrap();
    let transformed = centered.dot(&comp_matrix.t());

    let mean = mean_arr.to_vec();

    Ok(PcaResult {
        components,
        explained_variance: eigenvalues,
        explained_variance_ratio,
        transformed: transformed
            .into_raw_vec_and_offset()
            .0,
        mean,
        n_features,
        n_components,
    })
}

/// Power iteration using ndarray matrix-vector multiply.
fn ndarray_power_iteration(
    matrix: &Array2<f64>,
    max_iter: usize,
    tol: f64,
) -> (f64, Array1<f64>) {
    let n = matrix.nrows();

    // Initialize with deterministic non-zero vector
    let mut v = Array1::from_vec((0..n).map(|i| 1.0 / ((i + 1) as f64)).collect());
    let norm = v.dot(&v).sqrt();
    if norm > 0.0 {
        v /= norm;
    }

    let mut eigenvalue = 0.0;

    for _ in 0..max_iter {
        // w = M * v (BLAS gemv when available)
        let w = matrix.dot(&v);

        // Eigenvalue estimate
        let new_eigenvalue = v.dot(&w);

        // Normalize
        let wnorm = w.dot(&w).sqrt();
        if wnorm == 0.0 {
            break;
        }
        let w_normalized = &w / wnorm;

        // Check convergence
        let diff = (&v - &w_normalized)
            .mapv(|x| x * x)
            .sum()
            .sqrt();

        v = w_normalized;
        eigenvalue = new_eigenvalue;

        if diff < tol {
            break;
        }
    }

    (eigenvalue.abs(), v)
}
