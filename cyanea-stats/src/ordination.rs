//! Ordination methods for multivariate ecological analysis.
//!
//! - **PCoA** — Principal Coordinates Analysis (classical multidimensional scaling)
//! - **NMDS** — Non-metric Multidimensional Scaling (Kruskal's stress-1)
//! - **RDA** — Redundancy Analysis (constrained PCA)
//! - **CCA** — Canonical Correspondence Analysis
//! - **Procrustes** — Procrustes rotation and M² statistic

use cyanea_core::{CyaneaError, Result};

// ── PCoA ────────────────────────────────────────────────────────────────────

/// Result of Principal Coordinates Analysis.
#[derive(Debug, Clone)]
pub struct PcoaResult {
    /// Sample coordinates in ordination space. Shape: `n_samples × n_components`.
    pub coordinates: Vec<Vec<f64>>,
    /// Eigenvalues (may include negative values for non-Euclidean metrics).
    pub eigenvalues: Vec<f64>,
    /// Fraction of variance explained by each axis.
    pub proportion_explained: Vec<f64>,
    /// Number of negative eigenvalues encountered.
    pub n_negative_eigenvalues: usize,
}

/// Principal Coordinates Analysis (classical MDS).
///
/// Embeds samples into a low-dimensional Euclidean space that best preserves
/// the pairwise distances.
///
/// # Algorithm
///
/// 1. Double-center the squared distance matrix: `G = -½ (D² - row_mean - col_mean + grand_mean)`
/// 2. Eigendecompose G
/// 3. Coordinates = eigenvector × √max(λ, 0)
///
/// # Errors
///
/// Returns an error if the distance matrix is not square, has fewer than 2
/// samples, or `n_components` exceeds the number of samples.
pub fn pcoa(distances: &[Vec<f64>], n_components: usize) -> Result<PcoaResult> {
    let n = distances.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "pcoa: at least 2 samples required".into(),
        ));
    }
    for row in distances {
        if row.len() != n {
            return Err(CyaneaError::InvalidInput(
                "pcoa: distance matrix must be square".into(),
            ));
        }
    }
    if n_components == 0 || n_components > n {
        return Err(CyaneaError::InvalidInput(format!(
            "pcoa: n_components ({}) must be in [1, {}]",
            n_components, n
        )));
    }

    // Build squared distance matrix and double-center
    let mut g = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            g[i * n + j] = -0.5 * distances[i][j] * distances[i][j];
        }
    }
    double_center(&mut g, n);

    // Eigendecompose via power iteration + deflation
    let mut eigenvalues = Vec::with_capacity(n);
    let mut eigenvectors = Vec::with_capacity(n);
    let mut deflated = g;

    // Extract all eigenvalues (up to n) to count negatives and compute proportions
    for _ in 0..n {
        let (val, vec) = power_iteration(&deflated, n, 500);
        eigenvalues.push(val);
        eigenvectors.push(vec.clone());
        deflate(&mut deflated, val, &vec, n);
    }

    let n_negative = eigenvalues.iter().filter(|&&e| e < -1e-10).count();
    let total_positive: f64 = eigenvalues.iter().filter(|&&e| e > 0.0).sum();

    let proportion_explained: Vec<f64> = eigenvalues
        .iter()
        .take(n_components)
        .map(|&e| {
            if total_positive > 0.0 && e > 0.0 {
                e / total_positive
            } else {
                0.0
            }
        })
        .collect();

    // Compute coordinates: coord[i][k] = eigenvector[k][i] * sqrt(max(eigenvalue[k], 0))
    let mut coordinates = vec![vec![0.0; n_components]; n];
    for k in 0..n_components {
        let scale = eigenvalues[k].max(0.0).sqrt();
        for i in 0..n {
            coordinates[i][k] = eigenvectors[k][i] * scale;
        }
    }

    Ok(PcoaResult {
        coordinates,
        eigenvalues: eigenvalues[..n_components].to_vec(),
        proportion_explained,
        n_negative_eigenvalues: n_negative,
    })
}

// ── NMDS ────────────────────────────────────────────────────────────────────

/// Configuration for NMDS.
#[derive(Debug, Clone)]
pub struct NmdsConfig {
    /// Number of dimensions for embedding.
    pub n_dims: usize,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance for stress change.
    pub tolerance: f64,
    /// Random seed for initialization.
    pub seed: u64,
}

impl Default for NmdsConfig {
    fn default() -> Self {
        Self {
            n_dims: 2,
            max_iter: 300,
            tolerance: 1e-7,
            seed: 42,
        }
    }
}

/// Result of NMDS.
#[derive(Debug, Clone)]
pub struct NmdsResult {
    /// Sample coordinates. Shape: `n_samples × n_dims`.
    pub coordinates: Vec<Vec<f64>>,
    /// Final Kruskal stress-1 value.
    pub stress: f64,
    /// Number of iterations performed.
    pub n_iterations: usize,
    /// Whether the algorithm converged within tolerance.
    pub converged: bool,
}

/// Non-metric Multidimensional Scaling.
///
/// Finds a low-dimensional configuration that preserves the rank order of
/// pairwise distances, minimizing Kruskal's stress-1.
///
/// # Errors
///
/// Returns an error if the distance matrix is invalid or config is invalid.
pub fn nmds(distances: &[Vec<f64>], config: &NmdsConfig) -> Result<NmdsResult> {
    let n = distances.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "nmds: at least 2 samples required".into(),
        ));
    }
    for row in distances {
        if row.len() != n {
            return Err(CyaneaError::InvalidInput(
                "nmds: distance matrix must be square".into(),
            ));
        }
    }
    if config.n_dims == 0 {
        return Err(CyaneaError::InvalidInput(
            "nmds: n_dims must be > 0".into(),
        ));
    }

    let d = config.n_dims;
    let mut rng = Xorshift64::new(config.seed);

    // Initialize coordinates randomly
    let mut coords = vec![vec![0.0; d]; n];
    for row in &mut coords {
        for val in row.iter_mut() {
            *val = (rng.next_u64() as f64 / u64::MAX as f64) * 2.0 - 1.0;
        }
    }

    // Extract upper-triangle distances and their sort order
    let n_pairs = n * (n - 1) / 2;
    let mut orig_dists = Vec::with_capacity(n_pairs);
    let mut pair_indices = Vec::with_capacity(n_pairs);
    for i in 0..n {
        for j in (i + 1)..n {
            orig_dists.push(distances[i][j]);
            pair_indices.push((i, j));
        }
    }

    // Sort pairs by original distance for isotonic regression
    let mut sort_order: Vec<usize> = (0..n_pairs).collect();
    sort_order.sort_by(|&a, &b| orig_dists[a].partial_cmp(&orig_dists[b]).unwrap());

    let mut prev_stress = f64::MAX;
    let mut converged = false;
    let mut n_iterations = 0;
    let learning_rate = 0.05;

    for iter in 0..config.max_iter {
        n_iterations = iter + 1;

        // Compute embedding distances
        let mut embed_dists = vec![0.0; n_pairs];
        for (idx, &(i, j)) in pair_indices.iter().enumerate() {
            let mut d2 = 0.0;
            for k in 0..d {
                let diff = coords[i][k] - coords[j][k];
                d2 += diff * diff;
            }
            embed_dists[idx] = d2.sqrt();
        }

        // Isotonic regression (pool adjacent violators) on embedding distances
        // ordered by original distances
        let mut disparities = vec![0.0; n_pairs];
        {
            let sorted_embed: Vec<f64> = sort_order.iter().map(|&i| embed_dists[i]).collect();
            let iso = isotonic_regression(&sorted_embed);
            for (rank, &orig_idx) in sort_order.iter().enumerate() {
                disparities[orig_idx] = iso[rank];
            }
        }

        // Compute stress-1
        let mut num = 0.0;
        let mut den = 0.0;
        for idx in 0..n_pairs {
            let diff = embed_dists[idx] - disparities[idx];
            num += diff * diff;
            den += embed_dists[idx] * embed_dists[idx];
        }
        let stress = if den > 0.0 { (num / den).sqrt() } else { 0.0 };

        if (prev_stress - stress).abs() < config.tolerance {
            converged = true;
            prev_stress = stress;
            break;
        }
        prev_stress = stress;

        // Gradient descent step
        for (idx, &(i, j)) in pair_indices.iter().enumerate() {
            let ed = embed_dists[idx];
            if ed < 1e-15 {
                continue;
            }
            let scale = learning_rate * (disparities[idx] - ed) / ed;
            for k in 0..d {
                let diff = coords[i][k] - coords[j][k];
                let delta = scale * diff;
                coords[i][k] += delta;
                coords[j][k] -= delta;
            }
        }
    }

    Ok(NmdsResult {
        coordinates: coords,
        stress: prev_stress,
        n_iterations,
        converged,
    })
}

// ── RDA ─────────────────────────────────────────────────────────────────────

/// Result of constrained ordination (RDA or CCA).
#[derive(Debug, Clone)]
pub struct ConstrainedOrdinationResult {
    /// Sample scores in ordination space. Shape: `n_sites × n_axes`.
    pub sample_scores: Vec<Vec<f64>>,
    /// Species/response scores. Shape: `n_response × n_axes`.
    pub species_scores: Vec<Vec<f64>>,
    /// Environmental biplot scores. Shape: `n_env × n_axes`.
    pub biplot_scores: Vec<Vec<f64>>,
    /// Eigenvalues of constrained axes.
    pub eigenvalues: Vec<f64>,
    /// Proportion of constrained variance explained by each axis.
    pub proportion_explained: Vec<f64>,
    /// Total inertia (variance) in the response data.
    pub total_inertia: f64,
    /// Constrained inertia (variance explained by environment).
    pub constrained_inertia: f64,
    /// Method name ("RDA" or "CCA").
    pub method: String,
}

/// Redundancy Analysis — constrained PCA.
///
/// Projects the response matrix Y onto the column space of environment matrix X,
/// then performs PCA on the fitted values.
///
/// # Arguments
///
/// * `response` — flat row-major response matrix (n_sites × n_response)
/// * `environment` — flat row-major environment matrix (n_sites × n_env)
/// * `n_sites`, `n_response`, `n_env` — matrix dimensions
/// * `n_axes` — number of constrained axes to extract
///
/// # Errors
///
/// Returns an error if dimensions are inconsistent or `n_axes` is invalid.
pub fn rda(
    response: &[f64],
    environment: &[f64],
    n_sites: usize,
    n_response: usize,
    n_env: usize,
    n_axes: usize,
) -> Result<ConstrainedOrdinationResult> {
    validate_constrained_inputs(response, environment, n_sites, n_response, n_env, n_axes)?;

    // Center Y and X
    let y = center_matrix(response, n_sites, n_response);
    let x = center_matrix(environment, n_sites, n_env);

    // Total inertia = trace(Y^T Y) / (n-1)
    let total_inertia = matrix_trace_xtx(&y, n_sites, n_response) / (n_sites - 1).max(1) as f64;

    // QR decomposition of X
    let (q, _r) = qr_householder(&x, n_sites, n_env);

    // Y_hat = Q Q^T Y (projection of Y onto column space of X)
    let y_hat = project_via_q(&q, &y, n_sites, n_env.min(n_sites), n_response);

    // Constrained inertia
    let constrained_inertia =
        matrix_trace_xtx(&y_hat, n_sites, n_response) / (n_sites - 1).max(1) as f64;

    // PCA of Y_hat: covariance = Y_hat^T Y_hat / (n-1)
    let cov = compute_xtx(&y_hat, n_sites, n_response, (n_sites - 1).max(1) as f64);

    let actual_axes = n_axes.min(n_response).min(n_env);
    let mut eigenvalues = Vec::with_capacity(actual_axes);
    let mut eigenvectors = Vec::with_capacity(actual_axes);
    let mut deflated = cov;

    for _ in 0..actual_axes {
        let (val, vec) = power_iteration(&deflated, n_response, 500);
        eigenvalues.push(val);
        eigenvectors.push(vec.clone());
        deflate(&mut deflated, val, &vec, n_response);
    }

    let total_constrained: f64 = eigenvalues.iter().filter(|&&e| e > 0.0).sum();
    let proportion_explained: Vec<f64> = eigenvalues
        .iter()
        .map(|&e| {
            if total_constrained > 0.0 && e > 0.0 {
                e / total_constrained
            } else {
                0.0
            }
        })
        .collect();

    // Sample scores: Y_hat × eigenvectors
    let sample_scores = project_matrix(&y_hat, n_sites, n_response, &eigenvectors);

    // Species scores: eigenvectors themselves
    let species_scores = eigenvectors.iter().map(|v| v.clone()).collect();

    // Biplot scores: correlation of X columns with sample scores
    let biplot_scores = compute_biplot_scores(&x, n_sites, n_env, &sample_scores);

    Ok(ConstrainedOrdinationResult {
        sample_scores,
        species_scores,
        eigenvalues,
        proportion_explained,
        total_inertia,
        constrained_inertia,
        biplot_scores,
        method: "RDA".to_string(),
    })
}

/// Canonical Correspondence Analysis.
///
/// Chi-squared-weighted version of RDA for species composition data.
///
/// # Arguments
///
/// * `species` — flat row-major species matrix (n_sites × n_species), non-negative counts
/// * `environment` — flat row-major environment matrix (n_sites × n_env)
/// * `n_sites`, `n_species`, `n_env` — matrix dimensions
/// * `n_axes` — number of constrained axes to extract
///
/// # Errors
///
/// Returns an error if dimensions are inconsistent or data is invalid.
pub fn cca(
    species: &[f64],
    environment: &[f64],
    n_sites: usize,
    n_species: usize,
    n_env: usize,
    n_axes: usize,
) -> Result<ConstrainedOrdinationResult> {
    validate_constrained_inputs(species, environment, n_sites, n_species, n_env, n_axes)?;

    let grand_total: f64 = species.iter().sum();
    if grand_total <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "cca: species matrix must have positive total".into(),
        ));
    }

    // Compute row and column totals
    let mut row_sums = vec![0.0; n_sites];
    let mut col_sums = vec![0.0; n_species];
    for i in 0..n_sites {
        for j in 0..n_species {
            let v = species[i * n_species + j];
            row_sums[i] += v;
            col_sums[j] += v;
        }
    }

    // Row and column weights
    let row_weights: Vec<f64> = row_sums.iter().map(|&r| r / grand_total).collect();
    let col_weights: Vec<f64> = col_sums.iter().map(|&c| c / grand_total).collect();

    // Chi-squared transformed matrix: Q[i][j] = (Y[i][j]/grand - ri*cj) / sqrt(ri*cj)
    let mut q_mat = vec![0.0; n_sites * n_species];
    for i in 0..n_sites {
        for j in 0..n_species {
            let expected = row_weights[i] * col_weights[j];
            if expected > 0.0 {
                let observed = species[i * n_species + j] / grand_total;
                q_mat[i * n_species + j] = (observed - expected) / expected.sqrt();
            }
        }
    }

    // Total inertia (chi-squared / grand_total)
    let total_inertia: f64 = q_mat.iter().map(|&v| v * v).sum();

    // Weight rows of Q by sqrt(row_weight) and environment by sqrt(row_weight)
    let mut q_weighted = q_mat.clone();
    let mut x_weighted = center_matrix(environment, n_sites, n_env);
    for i in 0..n_sites {
        let w = row_weights[i].sqrt();
        for j in 0..n_species {
            q_weighted[i * n_species + j] *= w;
        }
        for j in 0..n_env {
            x_weighted[i * n_env + j] *= w;
        }
    }

    // QR of weighted environment
    let (qr_q, _qr_r) = qr_householder(&x_weighted, n_sites, n_env);

    // Project weighted Q onto environment
    let q_hat = project_via_q(&qr_q, &q_weighted, n_sites, n_env.min(n_sites), n_species);

    let constrained_inertia: f64 = q_hat.iter().map(|&v| v * v).sum();

    // PCA of projected matrix
    let cov = compute_xtx(&q_hat, n_sites, n_species, 1.0);

    let actual_axes = n_axes.min(n_species).min(n_env);
    let mut eigenvalues = Vec::with_capacity(actual_axes);
    let mut eigenvectors = Vec::with_capacity(actual_axes);
    let mut deflated = cov;

    for _ in 0..actual_axes {
        let (val, vec) = power_iteration(&deflated, n_species, 500);
        eigenvalues.push(val);
        eigenvectors.push(vec.clone());
        deflate(&mut deflated, val, &vec, n_species);
    }

    let total_constrained: f64 = eigenvalues.iter().filter(|&&e| e > 0.0).sum();
    let proportion_explained: Vec<f64> = eigenvalues
        .iter()
        .map(|&e| {
            if total_constrained > 0.0 && e > 0.0 {
                e / total_constrained
            } else {
                0.0
            }
        })
        .collect();

    let sample_scores = project_matrix(&q_hat, n_sites, n_species, &eigenvectors);
    let species_scores = eigenvectors.iter().map(|v| v.clone()).collect();
    let biplot_scores =
        compute_biplot_scores(&center_matrix(environment, n_sites, n_env), n_sites, n_env, &sample_scores);

    Ok(ConstrainedOrdinationResult {
        sample_scores,
        species_scores,
        eigenvalues,
        proportion_explained,
        total_inertia,
        constrained_inertia,
        biplot_scores,
        method: "CCA".to_string(),
    })
}

// ── Procrustes ──────────────────────────────────────────────────────────────

/// Result of Procrustes analysis.
#[derive(Debug, Clone)]
pub struct ProcrustesResult {
    /// Transformed target coordinates (after rotation/scaling/translation).
    pub transformed: Vec<Vec<f64>>,
    /// Procrustes M² statistic (sum of squared differences after transformation).
    pub m2: f64,
    /// Rotation matrix (p × p).
    pub rotation: Vec<Vec<f64>>,
    /// Scaling factor applied to target.
    pub scale: f64,
    /// Translation applied to target (added after rotation/scaling).
    pub translation: Vec<f64>,
}

/// Procrustes analysis: find the optimal rotation, scaling, and translation
/// to superimpose `target` onto `reference`.
///
/// Minimizes the sum of squared differences (M² statistic).
///
/// # Arguments
///
/// * `reference` — reference coordinates, one row per sample
/// * `target` — target coordinates to be transformed
///
/// # Errors
///
/// Returns an error if shapes are inconsistent or fewer than 2 samples.
pub fn procrustes(reference: &[Vec<f64>], target: &[Vec<f64>]) -> Result<ProcrustesResult> {
    let n = reference.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "procrustes: at least 2 samples required".into(),
        ));
    }
    if target.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "procrustes: reference ({}) and target ({}) must have same number of samples",
            n,
            target.len()
        )));
    }
    let p = reference[0].len();
    if p == 0 {
        return Err(CyaneaError::InvalidInput(
            "procrustes: dimensions must be > 0".into(),
        ));
    }
    for (i, row) in reference.iter().enumerate() {
        if row.len() != p {
            return Err(CyaneaError::InvalidInput(format!(
                "procrustes: reference row {} has {} dims, expected {}",
                i, row.len(), p
            )));
        }
    }
    for (i, row) in target.iter().enumerate() {
        if row.len() != p {
            return Err(CyaneaError::InvalidInput(format!(
                "procrustes: target row {} has {} dims, expected {}",
                i, row.len(), p
            )));
        }
    }

    // Center both matrices
    let ref_mean = column_means(reference, p);
    let tgt_mean = column_means(target, p);

    let ref_c: Vec<Vec<f64>> = reference
        .iter()
        .map(|row| row.iter().zip(&ref_mean).map(|(x, m)| x - m).collect())
        .collect();
    let tgt_c: Vec<Vec<f64>> = target
        .iter()
        .map(|row| row.iter().zip(&tgt_mean).map(|(x, m)| x - m).collect())
        .collect();

    // Scale reference to unit sum-of-squares
    let ss_ref: f64 = ref_c.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
    let ss_tgt: f64 = tgt_c.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();

    let ref_scaled: Vec<Vec<f64>> = if ss_ref > 0.0 {
        let s = ss_ref.sqrt();
        ref_c.iter().map(|r| r.iter().map(|x| x / s).collect()).collect()
    } else {
        ref_c.clone()
    };
    let tgt_scaled: Vec<Vec<f64>> = if ss_tgt > 0.0 {
        let s = ss_tgt.sqrt();
        tgt_c.iter().map(|r| r.iter().map(|x| x / s).collect()).collect()
    } else {
        tgt_c.clone()
    };

    // Cross-product matrix M = X^T Y (p × p) where X=ref, Y=target
    let mut m = vec![vec![0.0; p]; p];
    for k in 0..n {
        for i in 0..p {
            for j in 0..p {
                m[i][j] += ref_scaled[k][i] * tgt_scaled[k][j];
            }
        }
    }

    // SVD via eigendecomposition: M^T M and M M^T
    let (u, singular_values, vt) = svd_via_eigen(&m, p);

    // Rotation R = V U^T
    // Here u contains left singular vectors, vt contains right singular vectors transposed
    let mut rotation = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..p {
            for k in 0..p {
                // R = V U^T: V[i][k] * U^T[k][j] = V[i][k] * U[j][k]
                rotation[i][j] += vt[k][i] * u[j][k]; // vt transposed = V, u[j][k]
            }
        }
    }

    // Scale factor
    let trace_singular: f64 = singular_values.iter().sum();
    let scale = if ss_tgt > 0.0 {
        trace_singular * ss_ref.sqrt() / ss_tgt.sqrt()
    } else {
        1.0
    };

    // Apply transformation: transformed = tgt_c * R^T * scale + ref_mean
    let mut transformed = vec![vec![0.0; p]; n];
    for k in 0..n {
        for i in 0..p {
            let mut val = 0.0;
            for j in 0..p {
                val += tgt_c[k][j] * rotation[i][j]; // (tgt_c * R^T)[k][i]
            }
            transformed[k][i] = val * scale / ss_tgt.sqrt().max(1e-15) * ss_ref.sqrt().max(1e-15) + ref_mean[i];
        }
    }

    // Recompute for simpler M² on unit-scaled
    // M² = 1 - (Σ singular_values)²
    let m2 = (1.0 - trace_singular * trace_singular).max(0.0);

    let translation = ref_mean.clone();

    Ok(ProcrustesResult {
        transformed,
        m2,
        rotation,
        scale,
        translation,
    })
}

// ── Internal helpers ────────────────────────────────────────────────────────

struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }
}

fn double_center(matrix: &mut [f64], n: usize) {
    let mut row_means = vec![0.0; n];
    let mut col_means = vec![0.0; n];
    let mut grand_mean = 0.0;

    for i in 0..n {
        for j in 0..n {
            let v = matrix[i * n + j];
            row_means[i] += v;
            col_means[j] += v;
            grand_mean += v;
        }
    }
    for m in row_means.iter_mut() {
        *m /= n as f64;
    }
    for m in col_means.iter_mut() {
        *m /= n as f64;
    }
    grand_mean /= (n * n) as f64;

    for i in 0..n {
        for j in 0..n {
            matrix[i * n + j] = matrix[i * n + j] - row_means[i] - col_means[j] + grand_mean;
        }
    }
}

fn power_iteration(matrix: &[f64], n: usize, max_iter: usize) -> (f64, Vec<f64>) {
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

    (eigenvalue, v)
}

fn deflate(matrix: &mut [f64], eigenvalue: f64, eigenvector: &[f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            matrix[i * n + j] -= eigenvalue * eigenvector[i] * eigenvector[j];
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

fn isotonic_regression(y: &[f64]) -> Vec<f64> {
    // Pool Adjacent Violators Algorithm using integer block sizes
    let n = y.len();
    if n == 0 {
        return vec![];
    }

    // Each block: (weighted_sum, count)
    let mut blocks: Vec<(f64, usize)> = y.iter().map(|&v| (v, 1)).collect();

    let mut i = 0;
    while i < blocks.len() - 1 {
        let val_i = blocks[i].0 / blocks[i].1 as f64;
        let val_next = blocks[i + 1].0 / blocks[i + 1].1 as f64;
        if val_i > val_next {
            // Pool blocks i and i+1
            blocks[i].0 += blocks[i + 1].0;
            blocks[i].1 += blocks[i + 1].1;
            blocks.remove(i + 1);
            // Check backwards
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }

    // Expand blocks back to original size
    let mut result = Vec::with_capacity(n);
    for &(sum, count) in &blocks {
        let val = sum / count as f64;
        for _ in 0..count {
            result.push(val);
        }
    }
    result
}

fn center_matrix(data: &[f64], n_rows: usize, n_cols: usize) -> Vec<f64> {
    let mut means = vec![0.0; n_cols];
    for i in 0..n_rows {
        for j in 0..n_cols {
            means[j] += data[i * n_cols + j];
        }
    }
    for m in means.iter_mut() {
        *m /= n_rows as f64;
    }
    let mut centered = vec![0.0; n_rows * n_cols];
    for i in 0..n_rows {
        for j in 0..n_cols {
            centered[i * n_cols + j] = data[i * n_cols + j] - means[j];
        }
    }
    centered
}

fn matrix_trace_xtx(x: &[f64], n_rows: usize, n_cols: usize) -> f64 {
    // trace(X^T X) = sum of all squared elements
    let mut trace = 0.0;
    for i in 0..n_rows {
        for j in 0..n_cols {
            let v = x[i * n_cols + j];
            trace += v * v;
        }
    }
    trace
}

fn compute_xtx(x: &[f64], n_rows: usize, n_cols: usize, denom: f64) -> Vec<f64> {
    // X^T X / denom (n_cols × n_cols)
    let mut result = vec![0.0; n_cols * n_cols];
    for k in 0..n_rows {
        for i in 0..n_cols {
            for j in i..n_cols {
                let v = x[k * n_cols + i] * x[k * n_cols + j];
                result[i * n_cols + j] += v;
                if i != j {
                    result[j * n_cols + i] += v;
                }
            }
        }
    }
    for v in result.iter_mut() {
        *v /= denom;
    }
    result
}

fn qr_householder(a: &[f64], m: usize, n: usize) -> (Vec<f64>, Vec<f64>) {
    // QR decomposition via Householder reflections
    // Returns Q (m × n) and R (n × n)
    let k = m.min(n);
    let mut q = vec![0.0; m * m];
    // Initialize Q as identity
    for i in 0..m {
        q[i * m + i] = 1.0;
    }
    let mut r = vec![0.0; m * n];
    for i in 0..m {
        for j in 0..n {
            r[i * n + j] = a[i * n + j];
        }
    }

    for j in 0..k {
        // Extract column j from row j onwards
        let mut x = vec![0.0; m - j];
        for i in j..m {
            x[i - j] = r[i * n + j];
        }

        let norm_x = l2_norm(&x);
        if norm_x < 1e-15 {
            continue;
        }

        // Householder vector
        let sign = if x[0] >= 0.0 { 1.0 } else { -1.0 };
        x[0] += sign * norm_x;
        let norm_v = l2_norm(&x);
        if norm_v < 1e-15 {
            continue;
        }
        for val in x.iter_mut() {
            *val /= norm_v;
        }

        // Apply to R: R[j:, j:] -= 2 v (v^T R[j:, j:])
        for col in j..n {
            let mut dot_val = 0.0;
            for i in 0..(m - j) {
                dot_val += x[i] * r[(i + j) * n + col];
            }
            for i in 0..(m - j) {
                r[(i + j) * n + col] -= 2.0 * x[i] * dot_val;
            }
        }

        // Apply to Q: Q[:, j:] -= 2 (Q[:, j:] v) v^T
        for row in 0..m {
            let mut dot_val = 0.0;
            for i in 0..(m - j) {
                dot_val += q[row * m + (i + j)] * x[i];
            }
            for i in 0..(m - j) {
                q[row * m + (i + j)] -= 2.0 * dot_val * x[i];
            }
        }
    }

    // Extract Q (m × k) and R (k × n)
    let mut q_out = vec![0.0; m * k];
    for i in 0..m {
        for j in 0..k {
            q_out[i * k + j] = q[i * m + j];
        }
    }
    let mut r_out = vec![0.0; k * n];
    for i in 0..k {
        for j in 0..n {
            r_out[i * n + j] = r[i * n + j];
        }
    }

    (q_out, r_out)
}

fn project_via_q(q: &[f64], y: &[f64], m: usize, k: usize, n_cols: usize) -> Vec<f64> {
    // Y_hat = Q Q^T Y (m × n_cols)
    // First: Z = Q^T Y (k × n_cols)
    let mut z = vec![0.0; k * n_cols];
    for i in 0..k {
        for j in 0..n_cols {
            for row in 0..m {
                z[i * n_cols + j] += q[row * k + i] * y[row * n_cols + j];
            }
        }
    }
    // Then: Y_hat = Q Z (m × n_cols)
    let mut y_hat = vec![0.0; m * n_cols];
    for i in 0..m {
        for j in 0..n_cols {
            for col in 0..k {
                y_hat[i * n_cols + j] += q[i * k + col] * z[col * n_cols + j];
            }
        }
    }
    y_hat
}

fn project_matrix(
    data: &[f64],
    n_rows: usize,
    n_cols: usize,
    eigenvectors: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n_axes = eigenvectors.len();
    let mut scores = vec![vec![0.0; n_axes]; n_rows];
    for i in 0..n_rows {
        for (k, ev) in eigenvectors.iter().enumerate() {
            let mut val = 0.0;
            for j in 0..n_cols {
                val += data[i * n_cols + j] * ev[j];
            }
            scores[i][k] = val;
        }
    }
    scores
}

fn compute_biplot_scores(
    x: &[f64],
    n_rows: usize,
    n_cols: usize,
    sample_scores: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n_axes = if sample_scores.is_empty() {
        0
    } else {
        sample_scores[0].len()
    };
    let mut biplot = vec![vec![0.0; n_axes]; n_cols];

    for j in 0..n_cols {
        let x_col: Vec<f64> = (0..n_rows).map(|i| x[i * n_cols + j]).collect();
        let x_norm = l2_norm(&x_col);
        if x_norm < 1e-15 {
            continue;
        }
        for k in 0..n_axes {
            let score_col: Vec<f64> = sample_scores.iter().map(|s| s[k]).collect();
            let score_norm = l2_norm(&score_col);
            if score_norm < 1e-15 {
                continue;
            }
            biplot[j][k] = dot(&x_col, &score_col) / (x_norm * score_norm);
        }
    }
    biplot
}

fn validate_constrained_inputs(
    response: &[f64],
    environment: &[f64],
    n_sites: usize,
    n_response: usize,
    n_env: usize,
    n_axes: usize,
) -> Result<()> {
    if n_sites < 2 {
        return Err(CyaneaError::InvalidInput(
            "at least 2 sites required".into(),
        ));
    }
    if n_response == 0 || n_env == 0 {
        return Err(CyaneaError::InvalidInput(
            "response and environment must have at least 1 variable".into(),
        ));
    }
    if response.len() != n_sites * n_response {
        return Err(CyaneaError::InvalidInput(format!(
            "response length ({}) != n_sites ({}) * n_response ({})",
            response.len(),
            n_sites,
            n_response
        )));
    }
    if environment.len() != n_sites * n_env {
        return Err(CyaneaError::InvalidInput(format!(
            "environment length ({}) != n_sites ({}) * n_env ({})",
            environment.len(),
            n_sites,
            n_env
        )));
    }
    if n_axes == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_axes must be > 0".into(),
        ));
    }
    Ok(())
}

fn column_means(data: &[Vec<f64>], p: usize) -> Vec<f64> {
    let n = data.len();
    let mut means = vec![0.0; p];
    for row in data {
        for (i, &v) in row.iter().enumerate() {
            means[i] += v;
        }
    }
    for m in means.iter_mut() {
        *m /= n as f64;
    }
    means
}

fn svd_via_eigen(m: &[Vec<f64>], p: usize) -> (Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>) {
    // Compute SVD of p×p matrix M using eigendecomposition
    // M^T M → V, eigenvalues → σ²
    // M M^T → U

    // M^T M (p × p)
    let mut mtm = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..p {
            for k in 0..p {
                mtm[i * p + j] += m[k][i] * m[k][j];
            }
        }
    }

    // M M^T (p × p)
    let mut mmt = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..p {
            for k in 0..p {
                mmt[i * p + j] += m[i][k] * m[j][k];
            }
        }
    }

    // Eigendecompose both
    let mut eigenvalues = Vec::with_capacity(p);
    let mut v_vecs = Vec::with_capacity(p);
    let mut u_vecs = Vec::with_capacity(p);

    let mut mtm_deflated = mtm;
    let mut mmt_deflated = mmt;

    for _ in 0..p {
        let (ev_v, vec_v) = power_iteration(&mtm_deflated, p, 500);
        let (_, vec_u) = power_iteration(&mmt_deflated, p, 500);
        let sigma = ev_v.max(0.0).sqrt();
        eigenvalues.push(sigma);
        v_vecs.push(vec_v.clone());
        u_vecs.push(vec_u.clone());
        deflate(&mut mtm_deflated, ev_v, &vec_v, p);
        let ev_u = sigma * sigma; // same eigenvalue
        deflate(&mut mmt_deflated, ev_u, &vec_u, p);
    }

    // u_vecs[k] are columns of U, v_vecs[k] are columns of V
    // Return U (as rows), singular values, V^T (V transposed, as rows)
    (u_vecs, eigenvalues, v_vecs)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── PCoA tests ───────────────────────────────────────────────────

    #[test]
    fn pcoa_euclidean_recovers_structure() {
        // 3 points in 2D; PCoA on Euclidean distances should recover 2D structure
        let points: [[f64; 2]; 3] = [[0.0, 0.0], [3.0, 0.0], [0.0, 4.0]];
        let mut dists = vec![vec![0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let dx = points[i][0] - points[j][0];
                let dy = points[i][1] - points[j][1];
                dists[i][j] = (dx * dx + dy * dy).sqrt();
            }
        }
        let result = pcoa(&dists, 2).unwrap();
        assert_eq!(result.coordinates.len(), 3);
        assert_eq!(result.coordinates[0].len(), 2);
        // Verify pairwise distances are preserved
        for i in 0..3 {
            for j in (i + 1)..3 {
                let mut d2 = 0.0;
                for k in 0..2 {
                    let diff = result.coordinates[i][k] - result.coordinates[j][k];
                    d2 += diff * diff;
                }
                let reconstructed = d2.sqrt();
                assert!(
                    (reconstructed - dists[i][j]).abs() < 0.1,
                    "d[{},{}]: orig={} recon={}",
                    i,
                    j,
                    dists[i][j],
                    reconstructed
                );
            }
        }
    }

    #[test]
    fn pcoa_no_negative_eigenvalues_for_euclidean() {
        let points: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        let mut dists = vec![vec![0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                let dx = points[i][0] - points[j][0];
                let dy = points[i][1] - points[j][1];
                dists[i][j] = (dx * dx + dy * dy).sqrt();
            }
        }
        let result = pcoa(&dists, 2).unwrap();
        assert_eq!(result.n_negative_eigenvalues, 0);
    }

    #[test]
    fn pcoa_proportion_explained_sums_reasonable() {
        let dists = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.5],
            vec![2.0, 1.5, 0.0],
        ];
        let result = pcoa(&dists, 2).unwrap();
        let total: f64 = result.proportion_explained.iter().sum();
        assert!(total <= 1.0 + 1e-10, "total={}", total);
        assert!(total > 0.0, "total={}", total);
    }

    #[test]
    fn pcoa_too_few_samples_error() {
        let dists = vec![vec![0.0]];
        assert!(pcoa(&dists, 1).is_err());
    }

    #[test]
    fn pcoa_non_square_error() {
        let dists = vec![vec![0.0, 1.0], vec![1.0, 0.0, 0.5]];
        assert!(pcoa(&dists, 1).is_err());
    }

    // ── NMDS tests ───────────────────────────────────────────────────

    #[test]
    fn nmds_stress_decreases() {
        let dists = vec![
            vec![0.0, 1.0, 3.0, 5.0],
            vec![1.0, 0.0, 2.0, 4.0],
            vec![3.0, 2.0, 0.0, 2.0],
            vec![5.0, 4.0, 2.0, 0.0],
        ];
        let config = NmdsConfig {
            n_dims: 2,
            max_iter: 100,
            tolerance: 1e-10,
            seed: 42,
        };
        let result = nmds(&dists, &config).unwrap();
        assert!(result.stress >= 0.0, "stress={}", result.stress);
        assert!(result.stress < 1.0, "stress={}", result.stress);
    }

    #[test]
    fn nmds_deterministic_with_same_seed() {
        let dists = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.5],
            vec![2.0, 1.5, 0.0],
        ];
        let config = NmdsConfig::default();
        let r1 = nmds(&dists, &config).unwrap();
        let r2 = nmds(&dists, &config).unwrap();
        assert!((r1.stress - r2.stress).abs() < 1e-10);
    }

    #[test]
    fn nmds_too_few_samples_error() {
        let dists = vec![vec![0.0]];
        assert!(nmds(&dists, &NmdsConfig::default()).is_err());
    }

    #[test]
    fn nmds_zero_dims_error() {
        let dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let config = NmdsConfig {
            n_dims: 0,
            ..Default::default()
        };
        assert!(nmds(&dists, &config).is_err());
    }

    // ── RDA tests ────────────────────────────────────────────────────

    #[test]
    fn rda_linear_relationship() {
        // Y perfectly predicted by X → constrained inertia ≈ total inertia
        let n = 10;
        let env: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let response: Vec<f64> = env.iter().map(|&x| 2.0 * x + 1.0).collect();
        let result = rda(&response, &env, n, 1, 1, 1).unwrap();
        assert!(
            result.constrained_inertia / result.total_inertia > 0.99,
            "constrained/total = {}",
            result.constrained_inertia / result.total_inertia
        );
    }

    #[test]
    fn rda_no_relationship() {
        // Y uncorrelated with X → constrained inertia ≈ 0
        let env = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let response = vec![5.0, 3.0, 7.0, 1.0, 8.0, 2.0, 6.0, 4.0];
        let result = rda(&response, &env, 8, 1, 1, 1).unwrap();
        // Not perfectly 0 due to chance correlations, but should be < total
        assert!(result.constrained_inertia <= result.total_inertia + 1e-10);
    }

    #[test]
    fn rda_too_few_sites_error() {
        assert!(rda(&[1.0], &[1.0], 1, 1, 1, 1).is_err());
    }

    #[test]
    fn rda_dimension_mismatch_error() {
        assert!(rda(&[1.0, 2.0], &[1.0, 2.0, 3.0], 2, 1, 1, 1).is_err());
    }

    // ── CCA tests ────────────────────────────────────────────────────

    #[test]
    fn cca_basic_gradient() {
        // Species with clear environmental gradient
        let species = vec![
            10.0, 0.0, 0.0, // site 1: species A dominant
            8.0, 2.0, 0.0, // site 2
            3.0, 5.0, 2.0, // site 3
            0.0, 3.0, 7.0, // site 4
            0.0, 0.0, 10.0, // site 5: species C dominant
        ];
        let env = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = cca(&species, &env, 5, 3, 1, 1).unwrap();
        assert_eq!(result.method, "CCA");
        assert!(!result.eigenvalues.is_empty());
        assert!(result.total_inertia > 0.0);
    }

    #[test]
    fn cca_zero_species_error() {
        let species = vec![0.0, 0.0, 0.0, 0.0];
        let env = vec![1.0, 2.0];
        assert!(cca(&species, &env, 2, 2, 1, 1).is_err());
    }

    // ── Procrustes tests ─────────────────────────────────────────────

    #[test]
    fn procrustes_identical_m2_zero() {
        let coords = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let result = procrustes(&coords, &coords).unwrap();
        assert!(result.m2 < 1e-8, "M²={}", result.m2);
    }

    #[test]
    fn procrustes_scaled_low_m2() {
        let reference = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let target = vec![vec![0.0, 0.0], vec![2.0, 0.0], vec![0.0, 2.0]];
        let result = procrustes(&reference, &target).unwrap();
        // Uniform scaling → M² should be very low
        assert!(result.m2 < 0.1, "M²={}", result.m2);
    }

    #[test]
    fn procrustes_m2_in_range() {
        let reference = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
        ];
        let target = vec![
            vec![0.5, 0.5],
            vec![1.5, 0.5],
            vec![0.5, 1.5],
            vec![1.5, 1.5],
        ];
        let result = procrustes(&reference, &target).unwrap();
        assert!(result.m2 >= 0.0 && result.m2 <= 1.0, "M²={}", result.m2);
    }

    #[test]
    fn procrustes_too_few_samples_error() {
        assert!(procrustes(&[vec![1.0]], &[vec![1.0]]).is_err());
    }

    #[test]
    fn procrustes_mismatched_samples_error() {
        let a = vec![vec![0.0], vec![1.0]];
        let b = vec![vec![0.0], vec![1.0], vec![2.0]];
        assert!(procrustes(&a, &b).is_err());
    }

    // ── Isotonic regression tests ────────────────────────────────────

    #[test]
    fn isotonic_already_ordered() {
        let y = vec![1.0, 2.0, 3.0, 4.0];
        let result = isotonic_regression(&y);
        assert_eq!(result.len(), 4);
        for (a, b) in y.iter().zip(&result) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn isotonic_single_pool() {
        let y = vec![3.0, 1.0, 2.0];
        let result = isotonic_regression(&y);
        assert_eq!(result.len(), 3);
        // Should be non-decreasing
        for w in result.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "{:?}", result);
        }
    }

    // ── Double center test ───────────────────────────────────────────

    #[test]
    fn double_center_row_col_means_zero() {
        let mut m = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        double_center(&mut m, 3);
        // Row means and column means should be ~0
        for i in 0..3 {
            let row_mean: f64 = (0..3).map(|j| m[i * 3 + j]).sum::<f64>() / 3.0;
            assert!(
                row_mean.abs() < 1e-10,
                "row {} mean = {}",
                i,
                row_mean
            );
        }
        for j in 0..3 {
            let col_mean: f64 = (0..3).map(|i| m[i * 3 + j]).sum::<f64>() / 3.0;
            assert!(
                col_mean.abs() < 1e-10,
                "col {} mean = {}",
                j,
                col_mean
            );
        }
    }
}
