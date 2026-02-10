//! Dimensionality reduction: PCA and t-SNE.
//!
//! Provides unsupervised methods for projecting high-dimensional data into
//! lower-dimensional spaces for visualization and feature extraction.
//!
//! - **PCA** — Principal Component Analysis via power iteration on the
//!   covariance matrix (no LAPACK dependency).
//! - **t-SNE** — t-distributed Stochastic Neighbor Embedding for
//!   visualization of high-dimensional data in 2D/3D.

use cyanea_core::{CyaneaError, Result, Summarizable};

// ---------------------------------------------------------------------------
// PCA
// ---------------------------------------------------------------------------

/// Configuration for PCA.
#[derive(Debug, Clone)]
pub struct PcaConfig {
    /// Number of principal components to compute.
    pub n_components: usize,
    /// Maximum iterations for power iteration per component.
    pub max_iter: usize,
    /// Convergence tolerance for eigenvector.
    pub tolerance: f64,
}

impl Default for PcaConfig {
    fn default() -> Self {
        Self {
            n_components: 2,
            max_iter: 1000,
            tolerance: 1e-10,
        }
    }
}

/// Result of PCA computation.
#[derive(Debug, Clone)]
pub struct PcaResult {
    /// Principal components (eigenvectors), row-major: `n_components × n_features`.
    pub components: Vec<f64>,
    /// Eigenvalues (variance explained) for each component.
    pub explained_variance: Vec<f64>,
    /// Fraction of total variance explained by each component.
    pub explained_variance_ratio: Vec<f64>,
    /// Projected data, row-major: `n_samples × n_components`.
    pub transformed: Vec<f64>,
    /// Per-feature mean of the input data.
    pub mean: Vec<f64>,
    /// Number of features in the original data.
    pub n_features: usize,
    /// Number of components computed.
    pub n_components: usize,
}

impl Summarizable for PcaResult {
    fn summary(&self) -> String {
        let total: f64 = self.explained_variance_ratio.iter().sum();
        format!(
            "PCA: {} components, {:.1}% variance explained",
            self.n_components,
            total * 100.0,
        )
    }
}

/// Run PCA on a flat row-major data matrix.
///
/// `data` has shape `n_samples × n_features`, stored as a flat vector.
///
/// # Errors
///
/// Returns an error if the data is empty, dimensions are inconsistent,
/// or `n_components` exceeds `n_features`.
pub fn pca(data: &[f64], n_features: usize, config: &PcaConfig) -> Result<PcaResult> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    if n_features == 0 {
        return Err(CyaneaError::InvalidInput("n_features must be > 0".into()));
    }
    if data.len() % n_features != 0 {
        return Err(CyaneaError::InvalidInput(format!(
            "data length {} not divisible by n_features {}",
            data.len(),
            n_features
        )));
    }
    let n_samples = data.len() / n_features;
    if n_samples < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 samples".into(),
        ));
    }
    let n_components = config.n_components.min(n_features).min(n_samples);
    if n_components == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_components must be > 0".into(),
        ));
    }

    // Compute mean per feature
    let mut mean = vec![0.0; n_features];
    for row in 0..n_samples {
        for col in 0..n_features {
            mean[col] += data[row * n_features + col];
        }
    }
    let nf = n_samples as f64;
    for m in mean.iter_mut() {
        *m /= nf;
    }

    // Center the data
    let mut centered = data.to_vec();
    for row in 0..n_samples {
        for col in 0..n_features {
            centered[row * n_features + col] -= mean[col];
        }
    }

    // Compute covariance matrix: C = X^T X / (n - 1)
    let mut cov = vec![0.0; n_features * n_features];
    for row in 0..n_samples {
        let r = &centered[row * n_features..(row + 1) * n_features];
        for i in 0..n_features {
            for j in i..n_features {
                let val = r[i] * r[j];
                cov[i * n_features + j] += val;
                if i != j {
                    cov[j * n_features + i] += val;
                }
            }
        }
    }
    let scale = (n_samples - 1) as f64;
    for v in cov.iter_mut() {
        *v /= scale;
    }

    // Power iteration with deflation for each component
    let mut components = Vec::with_capacity(n_components * n_features);
    let mut eigenvalues = Vec::with_capacity(n_components);

    for _ in 0..n_components {
        let (eigenvalue, eigenvector) =
            power_iteration(&cov, n_features, config.max_iter, config.tolerance);

        components.extend_from_slice(&eigenvector);
        eigenvalues.push(eigenvalue);

        // Deflate: C = C - lambda * v * v^T
        for i in 0..n_features {
            for j in 0..n_features {
                cov[i * n_features + j] -= eigenvalue * eigenvector[i] * eigenvector[j];
            }
        }
    }

    // Compute explained variance ratio
    let total_variance: f64 = eigenvalues.iter().sum::<f64>();
    let explained_variance_ratio = if total_variance > 0.0 {
        // Total variance = trace of original covariance
        let mut orig_trace = 0.0;
        for row in 0..n_samples {
            let r = &centered[row * n_features..(row + 1) * n_features];
            for i in 0..n_features {
                orig_trace += r[i] * r[i];
            }
        }
        orig_trace /= scale;
        eigenvalues.iter().map(|&ev| ev / orig_trace).collect()
    } else {
        vec![0.0; n_components]
    };

    // Project data: transformed = centered * components^T
    let mut transformed = vec![0.0; n_samples * n_components];
    for row in 0..n_samples {
        let r = &centered[row * n_features..(row + 1) * n_features];
        for comp in 0..n_components {
            let c = &components[comp * n_features..(comp + 1) * n_features];
            let mut dot = 0.0;
            for k in 0..n_features {
                dot += r[k] * c[k];
            }
            transformed[row * n_components + comp] = dot;
        }
    }

    Ok(PcaResult {
        components,
        explained_variance: eigenvalues,
        explained_variance_ratio,
        transformed,
        mean,
        n_features,
        n_components,
    })
}

/// Power iteration: find the dominant eigenvector of a symmetric matrix.
fn power_iteration(
    matrix: &[f64],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> (f64, Vec<f64>) {
    // Initialize with [1, 0, 0, ...] then normalize
    let mut v = vec![0.0; n];
    // Use a deterministic non-zero init
    for (i, val) in v.iter_mut().enumerate() {
        *val = 1.0 / ((i + 1) as f64);
    }
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 0.0 {
        for val in v.iter_mut() {
            *val /= norm;
        }
    }

    let mut eigenvalue = 0.0;

    for _ in 0..max_iter {
        // w = M * v
        let mut w = vec![0.0; n];
        for i in 0..n {
            let mut sum = 0.0;
            for j in 0..n {
                sum += matrix[i * n + j] * v[j];
            }
            w[i] = sum;
        }

        // New eigenvalue estimate = v^T * w
        let new_eigenvalue: f64 = v.iter().zip(w.iter()).map(|(a, b)| a * b).sum();

        // Normalize w
        let wnorm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if wnorm == 0.0 {
            break;
        }
        for val in w.iter_mut() {
            *val /= wnorm;
        }

        // Check convergence
        let diff: f64 = v
            .iter()
            .zip(w.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();

        v = w;
        eigenvalue = new_eigenvalue;

        if diff < tol {
            break;
        }
    }

    (eigenvalue.abs(), v)
}

// ---------------------------------------------------------------------------
// t-SNE
// ---------------------------------------------------------------------------

/// Configuration for t-SNE.
#[derive(Debug, Clone)]
pub struct TsneConfig {
    /// Output dimensionality (typically 2 or 3).
    pub n_components: usize,
    /// Perplexity: effective number of neighbors (5-50 typical).
    pub perplexity: f64,
    /// Learning rate.
    pub learning_rate: f64,
    /// Number of gradient descent iterations.
    pub n_iter: usize,
    /// Random seed.
    pub seed: u64,
}

impl Default for TsneConfig {
    fn default() -> Self {
        Self {
            n_components: 2,
            perplexity: 30.0,
            learning_rate: 200.0,
            n_iter: 1000,
            seed: 42,
        }
    }
}

/// Result of t-SNE computation.
#[derive(Debug, Clone)]
pub struct TsneResult {
    /// Embedded coordinates, row-major: `n_samples × n_components`.
    pub embedding: Vec<f64>,
    /// Number of samples.
    pub n_samples: usize,
    /// Number of output components.
    pub n_components: usize,
    /// Final KL divergence.
    pub kl_divergence: f64,
}

impl Summarizable for TsneResult {
    fn summary(&self) -> String {
        format!(
            "t-SNE: {} samples in {}D, KL={:.4}",
            self.n_samples, self.n_components, self.kl_divergence,
        )
    }
}

/// Run t-SNE on a flat row-major data matrix.
///
/// `data` has shape `n_samples × n_features`.
///
/// # Errors
///
/// Returns an error if the data is empty or dimensions are inconsistent.
pub fn tsne(data: &[f64], n_features: usize, config: &TsneConfig) -> Result<TsneResult> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    if n_features == 0 {
        return Err(CyaneaError::InvalidInput("n_features must be > 0".into()));
    }
    if data.len() % n_features != 0 {
        return Err(CyaneaError::InvalidInput(format!(
            "data length {} not divisible by n_features {}",
            data.len(),
            n_features
        )));
    }
    let n = data.len() / n_features;
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 samples".into(),
        ));
    }
    if config.perplexity >= n as f64 {
        return Err(CyaneaError::InvalidInput(format!(
            "perplexity ({}) must be less than n_samples ({})",
            config.perplexity, n
        )));
    }
    let out_dim = config.n_components;

    // Compute pairwise squared distances
    #[cfg(feature = "parallel")]
    let sq_dists = {
        use rayon::prelude::*;
        let rows: Vec<Vec<(usize, f64)>> = (0..n)
            .into_par_iter()
            .map(|i| {
                let ri = &data[i * n_features..(i + 1) * n_features];
                ((i + 1)..n)
                    .map(|j| {
                        let rj = &data[j * n_features..(j + 1) * n_features];
                        let d: f64 = ri.iter().zip(rj).map(|(a, b)| (a - b).powi(2)).sum();
                        (j, d)
                    })
                    .collect()
            })
            .collect();
        let mut sq_dists = vec![0.0; n * n];
        for (i, row) in rows.into_iter().enumerate() {
            for (j, d) in row {
                sq_dists[i * n + j] = d;
                sq_dists[j * n + i] = d;
            }
        }
        sq_dists
    };
    #[cfg(not(feature = "parallel"))]
    let sq_dists = {
        let mut sq_dists = vec![0.0; n * n];
        for i in 0..n {
            let ri = &data[i * n_features..(i + 1) * n_features];
            for j in (i + 1)..n {
                let rj = &data[j * n_features..(j + 1) * n_features];
                let d: f64 = ri.iter().zip(rj).map(|(a, b)| (a - b).powi(2)).sum();
                sq_dists[i * n + j] = d;
                sq_dists[j * n + i] = d;
            }
        }
        sq_dists
    };

    // Compute pairwise affinities P (symmetrized)
    let p = compute_joint_probabilities(&sq_dists, n, config.perplexity);

    // Initialize embedding with small random values
    let mut rng = Xorshift64(config.seed.max(1));
    let mut y = vec![0.0; n * out_dim];
    for val in y.iter_mut() {
        *val = rng.next_f64() * 0.01 - 0.005;
    }
    let mut gains = vec![1.0; n * out_dim];
    let mut velocities = vec![0.0; n * out_dim];

    let momentum_switch = 250;
    let initial_momentum = 0.5;
    let final_momentum = 0.8;

    let mut kl_div = 0.0;

    for iter in 0..config.n_iter {
        let momentum = if iter < momentum_switch {
            initial_momentum
        } else {
            final_momentum
        };

        // Compute Q distribution (Student-t with 1 DOF)
        let mut q_num = vec![0.0; n * n];
        let mut q_sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let mut d2 = 0.0;
                for d in 0..out_dim {
                    let diff = y[i * out_dim + d] - y[j * out_dim + d];
                    d2 += diff * diff;
                }
                let val = 1.0 / (1.0 + d2);
                q_num[i * n + j] = val;
                q_num[j * n + i] = val;
                q_sum += 2.0 * val;
            }
        }
        if q_sum == 0.0 {
            q_sum = 1.0;
        }

        // Compute gradients
        let mut grad = vec![0.0; n * out_dim];
        kl_div = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let q_ij = (q_num[i * n + j] / q_sum).max(1e-12);
                let p_ij = p[i * n + j];
                let mult = 4.0 * (p_ij - q_ij) * q_num[i * n + j];
                for d in 0..out_dim {
                    grad[i * out_dim + d] += mult * (y[i * out_dim + d] - y[j * out_dim + d]);
                }
                if p_ij > 1e-12 {
                    kl_div += p_ij * (p_ij / q_ij).ln();
                }
            }
        }
        // KL is double-counted for i<j pairs
        kl_div /= 2.0;

        // Update with momentum and adaptive gains
        for idx in 0..y.len() {
            let sign_match = (grad[idx] > 0.0) == (velocities[idx] > 0.0);
            gains[idx] = if sign_match {
                (gains[idx] * 0.8_f64).max(0.01)
            } else {
                gains[idx] + 0.2
            };
            velocities[idx] =
                momentum * velocities[idx] - config.learning_rate * gains[idx] * grad[idx];
            y[idx] += velocities[idx];
        }

        // Re-center
        for d in 0..out_dim {
            let mean: f64 = (0..n).map(|i| y[i * out_dim + d]).sum::<f64>() / n as f64;
            for i in 0..n {
                y[i * out_dim + d] -= mean;
            }
        }
    }

    Ok(TsneResult {
        embedding: y,
        n_samples: n,
        n_components: out_dim,
        kl_divergence: kl_div,
    })
}

/// Compute symmetrized joint probabilities P from squared distances.
fn compute_joint_probabilities(sq_dists: &[f64], n: usize, perplexity: f64) -> Vec<f64> {
    let target_entropy = perplexity.ln();
    let mut p = vec![0.0; n * n];

    // For each point, find sigma via binary search on perplexity
    for i in 0..n {
        let mut lo = 1e-10_f64;
        let mut hi = 1e4_f64;
        let mut sigma = 1.0;

        for _ in 0..50 {
            sigma = (lo + hi) / 2.0;
            let beta = 1.0 / (2.0 * sigma * sigma);

            // Compute conditional probabilities p(j|i)
            let mut sum_exp = 0.0;
            for j in 0..n {
                if j == i {
                    continue;
                }
                sum_exp += (-beta * sq_dists[i * n + j]).exp();
            }

            if sum_exp == 0.0 {
                lo = sigma;
                continue;
            }

            // Compute entropy
            let mut entropy = 0.0;
            for j in 0..n {
                if j == i {
                    continue;
                }
                let pj = (-beta * sq_dists[i * n + j]).exp() / sum_exp;
                if pj > 1e-12 {
                    entropy -= pj * pj.ln();
                }
            }

            if entropy > target_entropy {
                hi = sigma;
            } else {
                lo = sigma;
            }
        }

        // Store conditional probabilities
        let beta = 1.0 / (2.0 * sigma * sigma);
        let mut sum_exp = 0.0;
        for j in 0..n {
            if j == i {
                continue;
            }
            sum_exp += (-beta * sq_dists[i * n + j]).exp();
        }
        if sum_exp > 0.0 {
            for j in 0..n {
                if j == i {
                    continue;
                }
                p[i * n + j] = (-beta * sq_dists[i * n + j]).exp() / sum_exp;
            }
        }
    }

    // Symmetrize: P_ij = (p(j|i) + p(i|j)) / (2n)
    let scale = 1.0 / (2.0 * n as f64);
    let mut sym = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let val = (p[i * n + j] + p[j * n + i]) * scale;
            let val = val.max(1e-12);
            sym[i * n + j] = val;
            sym[j * n + i] = val;
        }
    }
    sym
}

/// Minimal xorshift64 PRNG.
struct Xorshift64(u64);

impl Xorshift64 {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- PCA tests ---

    #[test]
    fn pca_basic_2d() {
        // 4 points in 2D, main variance along x
        let data = vec![
            0.0, 0.0,
            1.0, 0.1,
            2.0, 0.2,
            3.0, 0.3,
        ];
        let config = PcaConfig {
            n_components: 2,
            ..Default::default()
        };
        let result = pca(&data, 2, &config).unwrap();
        assert_eq!(result.n_components, 2);
        assert_eq!(result.transformed.len(), 4 * 2);
        // First component should capture most variance
        assert!(result.explained_variance[0] > result.explained_variance[1]);
    }

    #[test]
    fn pca_one_component() {
        let data = vec![
            0.0, 0.0,
            1.0, 1.0,
            2.0, 2.0,
            3.0, 3.0,
        ];
        let config = PcaConfig {
            n_components: 1,
            ..Default::default()
        };
        let result = pca(&data, 2, &config).unwrap();
        assert_eq!(result.n_components, 1);
        assert_eq!(result.transformed.len(), 4);
    }

    #[test]
    fn pca_mean_centering() {
        let data = vec![
            100.0, 200.0,
            101.0, 201.0,
            102.0, 202.0,
            103.0, 203.0,
        ];
        let config = PcaConfig::default();
        let result = pca(&data, 2, &config).unwrap();
        assert!((result.mean[0] - 101.5).abs() < 1e-10);
        assert!((result.mean[1] - 201.5).abs() < 1e-10);
    }

    #[test]
    fn pca_empty_error() {
        let config = PcaConfig::default();
        assert!(pca(&[], 2, &config).is_err());
    }

    #[test]
    fn pca_dimension_mismatch() {
        let data = vec![1.0, 2.0, 3.0];
        let config = PcaConfig::default();
        assert!(pca(&data, 2, &config).is_err()); // 3 not divisible by 2
    }

    #[test]
    fn pca_summary() {
        let data = vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
        let config = PcaConfig::default();
        let result = pca(&data, 2, &config).unwrap();
        let s = result.summary();
        assert!(s.contains("PCA"));
        assert!(s.contains("2 components"));
    }

    #[test]
    fn pca_variance_ratio_sums_leq_one() {
        let data = vec![
            0.0, 0.0, 0.0,
            1.0, 0.1, 0.0,
            2.0, 0.2, 0.1,
            3.0, 0.3, 0.1,
        ];
        let config = PcaConfig {
            n_components: 3,
            ..Default::default()
        };
        let result = pca(&data, 3, &config).unwrap();
        let ratio_sum: f64 = result.explained_variance_ratio.iter().sum();
        assert!(ratio_sum <= 1.0 + 1e-6);
        assert!(ratio_sum > 0.0);
    }

    // --- t-SNE tests ---

    #[test]
    fn tsne_basic_clusters() {
        // Two well-separated clusters in 3D
        let mut data = Vec::new();
        for _ in 0..5 {
            data.extend_from_slice(&[0.0, 0.0, 0.0]);
        }
        for _ in 0..5 {
            data.extend_from_slice(&[10.0, 10.0, 10.0]);
        }
        let config = TsneConfig {
            n_components: 2,
            perplexity: 3.0,
            n_iter: 200,
            ..Default::default()
        };
        let result = tsne(&data, 3, &config).unwrap();
        assert_eq!(result.n_samples, 10);
        assert_eq!(result.n_components, 2);
        assert_eq!(result.embedding.len(), 20);
    }

    #[test]
    fn tsne_3d_output() {
        let data = vec![
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0,
            5.0, 5.0,
            6.0, 5.0,
        ];
        let config = TsneConfig {
            n_components: 3,
            perplexity: 2.0,
            n_iter: 100,
            ..Default::default()
        };
        let result = tsne(&data, 2, &config).unwrap();
        assert_eq!(result.n_components, 3);
        assert_eq!(result.embedding.len(), 6 * 3);
    }

    #[test]
    fn tsne_empty_error() {
        let config = TsneConfig::default();
        assert!(tsne(&[], 2, &config).is_err());
    }

    #[test]
    fn tsne_perplexity_too_high() {
        let data = vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0];
        let config = TsneConfig {
            perplexity: 10.0, // >= n_samples=3
            ..Default::default()
        };
        assert!(tsne(&data, 2, &config).is_err());
    }

    #[test]
    fn tsne_summary() {
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let config = TsneConfig {
            perplexity: 1.5,
            n_iter: 50,
            ..Default::default()
        };
        let result = tsne(&data, 1, &config).unwrap();
        let s = result.summary();
        assert!(s.contains("t-SNE"));
        assert!(s.contains("2D"));
    }
}
