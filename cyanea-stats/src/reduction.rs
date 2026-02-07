//! Dimensionality reduction — Principal Component Analysis (PCA).
//!
//! Provides a pure-Rust PCA implementation using eigendecomposition of the
//! covariance matrix via power iteration. No external linear algebra
//! dependencies are required.

use cyanea_core::{CyaneaError, Result};

/// Result of a PCA computation.
#[derive(Debug, Clone)]
pub struct PcaResult {
    /// Principal component vectors (directions of maximum variance).
    /// Shape: `n_components × n_features`.
    pub components: Vec<Vec<f64>>,
    /// Variance explained by each component.
    pub explained_variance: Vec<f64>,
    /// Fraction of total variance explained by each component.
    pub explained_variance_ratio: Vec<f64>,
    /// Projected data in PC space. Shape: `n_samples × n_components`.
    pub transformed: Vec<Vec<f64>>,
    /// Mean of each feature (used for centering).
    pub mean: Vec<f64>,
}

/// Perform PCA on a data matrix.
///
/// `data` is a slice of sample vectors (each sample is a row of feature values).
/// `n_components` is the number of principal components to extract.
///
/// # Errors
///
/// Returns an error if data is empty, samples have inconsistent lengths,
/// or `n_components` is zero or exceeds the number of features.
pub fn pca(data: &[&[f64]], n_components: usize) -> Result<PcaResult> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    let n_samples = data.len();
    let n_features = data[0].len();
    if n_features == 0 {
        return Err(CyaneaError::InvalidInput("zero features".into()));
    }
    for (i, row) in data.iter().enumerate() {
        if row.len() != n_features {
            return Err(CyaneaError::InvalidInput(format!(
                "sample {} has {} features, expected {}",
                i,
                row.len(),
                n_features
            )));
        }
    }
    if n_components == 0 || n_components > n_features {
        return Err(CyaneaError::InvalidInput(format!(
            "n_components ({}) must be in [1, {}]",
            n_components, n_features
        )));
    }

    // Center data
    let mean = compute_mean(data, n_features);
    let centered: Vec<Vec<f64>> = data
        .iter()
        .map(|row| row.iter().zip(&mean).map(|(x, m)| x - m).collect())
        .collect();

    // Covariance matrix (n_features × n_features)
    let cov = covariance_matrix(&centered, n_features, n_samples);

    let total_variance: f64 = (0..n_features).map(|i| cov[i * n_features + i]).sum();

    // Extract eigenvectors via power iteration with deflation
    let mut components = Vec::with_capacity(n_components);
    let mut eigenvalues = Vec::with_capacity(n_components);
    let mut deflated = cov;

    for _ in 0..n_components {
        let (eigenvalue, eigenvector) = power_iteration(&deflated, n_features, 300);
        eigenvalues.push(eigenvalue);
        components.push(eigenvector.clone());
        deflate(&mut deflated, eigenvalue, &eigenvector, n_features);
    }

    let explained_variance_ratio: Vec<f64> = eigenvalues
        .iter()
        .map(|&e| {
            if total_variance > 0.0 {
                e / total_variance
            } else {
                0.0
            }
        })
        .collect();

    // Project data onto components
    let transformed: Vec<Vec<f64>> = centered
        .iter()
        .map(|sample| components.iter().map(|comp| dot(sample, comp)).collect())
        .collect();

    Ok(PcaResult {
        components,
        explained_variance: eigenvalues,
        explained_variance_ratio,
        transformed,
        mean,
    })
}

fn compute_mean(data: &[&[f64]], n_features: usize) -> Vec<f64> {
    let n = data.len() as f64;
    let mut mean = vec![0.0; n_features];
    for row in data {
        for (i, &val) in row.iter().enumerate() {
            mean[i] += val;
        }
    }
    for m in mean.iter_mut() {
        *m /= n;
    }
    mean
}

fn covariance_matrix(centered: &[Vec<f64>], n_features: usize, n_samples: usize) -> Vec<f64> {
    let denom = if n_samples > 1 {
        (n_samples - 1) as f64
    } else {
        1.0
    };
    let mut cov = vec![0.0; n_features * n_features];
    for sample in centered {
        for i in 0..n_features {
            for j in i..n_features {
                let v = sample[i] * sample[j];
                cov[i * n_features + j] += v;
                if i != j {
                    cov[j * n_features + i] += v;
                }
            }
        }
    }
    for c in cov.iter_mut() {
        *c /= denom;
    }
    cov
}

fn power_iteration(matrix: &[f64], n: usize, max_iter: usize) -> (f64, Vec<f64>) {
    // Initialize with a non-degenerate vector
    let mut v: Vec<f64> = (0..n).map(|i| 1.0 / ((i + 1) as f64)).collect();
    normalize(&mut v);

    let mut eigenvalue = 0.0;

    for _ in 0..max_iter {
        let mut new_v = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                new_v[i] += matrix[i * n + j] * v[j];
            }
        }
        eigenvalue = dot(&new_v, &v);
        let norm = l2_norm(&new_v);
        if norm < 1e-15 {
            break;
        }
        for val in new_v.iter_mut() {
            *val /= norm;
        }
        let diff: f64 = new_v.iter().zip(&v).map(|(a, b)| (a - b).powi(2)).sum();
        v = new_v;
        if diff < 1e-12 {
            break;
        }
    }

    (eigenvalue.max(0.0), v)
}

fn deflate(cov: &mut [f64], eigenvalue: f64, eigenvector: &[f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            cov[i * n + j] -= eigenvalue * eigenvector[i] * eigenvector[j];
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn normalize(v: &mut [f64]) {
    let n = l2_norm(v);
    if n > 0.0 {
        for val in v.iter_mut() {
            *val /= n;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pca_basic_2d() {
        let data: Vec<Vec<f64>> = vec![
            vec![1.0, 0.1],
            vec![2.0, 0.2],
            vec![3.0, 0.3],
            vec![4.0, 0.4],
            vec![5.0, 0.5],
        ];
        let refs: Vec<&[f64]> = data.iter().map(|v| v.as_slice()).collect();
        let result = pca(&refs, 2).unwrap();
        assert_eq!(result.components.len(), 2);
        assert_eq!(result.transformed.len(), 5);
        assert_eq!(result.transformed[0].len(), 2);
        assert!(result.explained_variance_ratio[0] > 0.9);
    }

    #[test]
    fn pca_single_component() {
        let data: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let refs: Vec<&[f64]> = data.iter().map(|v| v.as_slice()).collect();
        let result = pca(&refs, 1).unwrap();
        assert_eq!(result.components.len(), 1);
        assert_eq!(result.components[0].len(), 3);
    }

    #[test]
    fn pca_mean_centering() {
        let data: Vec<Vec<f64>> = vec![vec![10.0, 20.0], vec![12.0, 22.0], vec![14.0, 24.0]];
        let refs: Vec<&[f64]> = data.iter().map(|v| v.as_slice()).collect();
        let result = pca(&refs, 1).unwrap();
        assert!((result.mean[0] - 12.0).abs() < 1e-10);
        assert!((result.mean[1] - 22.0).abs() < 1e-10);
    }

    #[test]
    fn pca_empty_data_error() {
        let data: Vec<&[f64]> = vec![];
        assert!(pca(&data, 1).is_err());
    }

    #[test]
    fn pca_too_many_components_error() {
        let data: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let refs: Vec<&[f64]> = data.iter().map(|v| v.as_slice()).collect();
        assert!(pca(&refs, 3).is_err());
    }

    #[test]
    fn pca_variance_ratio_sums_near_one() {
        let data: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.5],
            vec![2.0, 1.0, 0.3],
            vec![3.0, 2.0, 0.1],
            vec![4.0, 3.0, 0.7],
        ];
        let refs: Vec<&[f64]> = data.iter().map(|v| v.as_slice()).collect();
        let result = pca(&refs, 3).unwrap();
        let total: f64 = result.explained_variance_ratio.iter().sum();
        assert!((total - 1.0).abs() < 0.15);
    }
}
