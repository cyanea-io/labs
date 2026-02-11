//! UMAP (Uniform Manifold Approximation and Projection).
//!
//! Nonlinear dimensionality reduction that preserves both local and global
//! structure better than t-SNE, with faster convergence. Particularly suited
//! for biological data (single-cell, k-mer embeddings, expression matrices).

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::distance::{compute_distance, DistanceMetric};

// ---------------------------------------------------------------------------
// Config & Result types
// ---------------------------------------------------------------------------

/// Initialization strategy for UMAP embedding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmapInit {
    /// Random initialization using Xorshift64.
    Random,
    /// PCA-based initialization (reuses `crate::reduction::pca`).
    Pca,
}

/// Configuration for UMAP.
#[derive(Debug, Clone)]
pub struct UmapConfig {
    /// Number of output dimensions.
    pub n_components: usize,
    /// Number of nearest neighbors for graph construction.
    pub n_neighbors: usize,
    /// Minimum distance between points in the embedding.
    pub min_dist: f64,
    /// Effective scale of embedded points.
    pub spread: f64,
    /// Learning rate for SGD.
    pub learning_rate: f64,
    /// Number of optimization epochs.
    pub n_epochs: usize,
    /// Number of negative samples per positive edge.
    pub negative_sample_rate: usize,
    /// Distance metric for KNN computation.
    pub metric: DistanceMetric,
    /// Embedding initialization method.
    pub init: UmapInit,
    /// Random seed for reproducibility.
    pub seed: u64,
}

impl Default for UmapConfig {
    fn default() -> Self {
        Self {
            n_components: 2,
            n_neighbors: 15,
            min_dist: 0.1,
            spread: 1.0,
            learning_rate: 1.0,
            n_epochs: 200,
            negative_sample_rate: 5,
            metric: DistanceMetric::Euclidean,
            init: UmapInit::Pca,
            seed: 42,
        }
    }
}

/// Result of UMAP computation.
#[derive(Debug, Clone)]
pub struct UmapResult {
    /// Embedded coordinates, row-major: `n_samples * n_components`.
    pub embedding: Vec<f64>,
    /// Number of samples.
    pub n_samples: usize,
    /// Number of output components.
    pub n_components: usize,
    /// Number of optimization epochs run.
    pub n_epochs: usize,
}

impl Summarizable for UmapResult {
    fn summary(&self) -> String {
        format!(
            "UMAP: {} samples in {}D, {} epochs",
            self.n_samples, self.n_components, self.n_epochs,
        )
    }
}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// Sparse edge with symmetrized weight.
#[derive(Debug, Clone)]
struct Edge {
    i: usize,
    j: usize,
    weight: f64,
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

    fn next_bounded(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Run UMAP on a flat row-major data matrix.
///
/// `data` has shape `n_samples * n_features`, stored as a flat vector.
///
/// # Errors
///
/// Returns an error if the data is empty, dimensions are inconsistent,
/// or `n_neighbors` is out of range.
pub fn umap(data: &[f64], n_features: usize, config: &UmapConfig) -> Result<UmapResult> {
    // --- Phase 1: Input validation ---
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
    if config.n_neighbors < 2 {
        return Err(CyaneaError::InvalidInput(
            "n_neighbors must be >= 2".into(),
        ));
    }
    if config.n_neighbors >= n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "n_neighbors ({}) must be < n_samples ({})",
            config.n_neighbors, n_samples
        )));
    }
    if config.n_components == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_components must be > 0".into(),
        ));
    }

    let k = config.n_neighbors;
    let out_dim = config.n_components;

    // --- Phase 2: KNN graph ---
    let (knn_indices, knn_dists) =
        compute_knn_graph(data, n_features, n_samples, k, config.metric)?;

    // --- Phase 3: Smooth KNN distances ---
    let sigmas = smooth_knn_distances(&knn_dists, n_samples, k);

    // --- Phase 4: Fuzzy graph + symmetrize ---
    let edges = compute_fuzzy_graph(&knn_indices, &knn_dists, &sigmas, n_samples, k);

    // --- Phase 5: Fit a/b params ---
    let (a, b) = fit_ab_params(config.min_dist, config.spread);

    // --- Phase 6: Initialize embedding ---
    let mut embedding = initialize_embedding(
        data,
        n_features,
        n_samples,
        out_dim,
        config.init,
        config.seed,
    )?;

    // --- Phase 7: SGD optimization ---
    optimize_embedding(
        &mut embedding,
        &edges,
        n_samples,
        out_dim,
        a,
        b,
        config.learning_rate,
        config.n_epochs,
        config.negative_sample_rate,
        config.seed,
    );

    Ok(UmapResult {
        embedding,
        n_samples,
        n_components: out_dim,
        n_epochs: config.n_epochs,
    })
}

// ---------------------------------------------------------------------------
// Phase 2: KNN graph
// ---------------------------------------------------------------------------

/// Brute-force KNN graph. Returns (indices, distances), each n_samples * k.
fn compute_knn_graph(
    data: &[f64],
    n_features: usize,
    n_samples: usize,
    k: usize,
    metric: DistanceMetric,
) -> Result<(Vec<Vec<usize>>, Vec<Vec<f64>>)> {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let results: Vec<Result<(Vec<usize>, Vec<f64>)>> = (0..n_samples)
            .into_par_iter()
            .map(|i| {
                let ri = &data[i * n_features..(i + 1) * n_features];
                let mut dists: Vec<(usize, f64)> = Vec::with_capacity(n_samples - 1);
                for j in 0..n_samples {
                    if j == i {
                        continue;
                    }
                    let rj = &data[j * n_features..(j + 1) * n_features];
                    let d = compute_distance(ri, rj, metric)?;
                    dists.push((j, d));
                }
                dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
                dists.truncate(k);
                let (idx, dst): (Vec<_>, Vec<_>) = dists.into_iter().unzip();
                Ok((idx, dst))
            })
            .collect();
        let mut all_indices = Vec::with_capacity(n_samples);
        let mut all_dists = Vec::with_capacity(n_samples);
        for r in results {
            let (idx, dst) = r?;
            all_indices.push(idx);
            all_dists.push(dst);
        }
        return Ok((all_indices, all_dists));
    }

    #[cfg(not(feature = "parallel"))]
    {
        let mut all_indices = Vec::with_capacity(n_samples);
        let mut all_dists = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let ri = &data[i * n_features..(i + 1) * n_features];
            let mut dists: Vec<(usize, f64)> = Vec::with_capacity(n_samples - 1);
            for j in 0..n_samples {
                if j == i {
                    continue;
                }
                let rj = &data[j * n_features..(j + 1) * n_features];
                let d = compute_distance(ri, rj, metric)?;
                dists.push((j, d));
            }
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            dists.truncate(k);
            let (idx, dst): (Vec<_>, Vec<_>) = dists.into_iter().unzip();
            all_indices.push(idx);
            all_dists.push(dst);
        }
        Ok((all_indices, all_dists))
    }
}

// ---------------------------------------------------------------------------
// Phase 3: Smooth KNN distances (binary search for sigma)
// ---------------------------------------------------------------------------

/// For each point, find sigma_i such that sum of exp(-d / sigma) over k
/// neighbors yields the target log2(n_neighbors).
fn smooth_knn_distances(
    knn_dists: &[Vec<f64>],
    n_samples: usize,
    k: usize,
) -> Vec<f64> {
    let target = (k as f64).ln();
    let mut sigmas = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let dists = &knn_dists[i];
        let rho = dists[0].max(0.0); // distance to nearest neighbor

        let mut lo = 1e-10_f64;
        let mut hi = 1e4_f64;
        let mut sigma = 1.0;

        for _ in 0..64 {
            sigma = (lo + hi) / 2.0;
            let mut sum = 0.0;
            for &d in dists.iter() {
                let shifted = (d - rho).max(0.0);
                sum += (-shifted / sigma).exp();
            }
            // We subtract the self-contribution (d=0 gives exp(0)=1 after shift)
            // but rho is the nearest-neighbor dist so the first neighbor's shifted
            // dist is 0, giving exp(0)=1. The entropy target doesn't include self.
            let entropy = sum.ln();

            if (entropy - target).abs() < 1e-5 {
                break;
            }
            if entropy > target {
                hi = sigma;
            } else {
                lo = sigma;
            }
        }

        sigmas.push(sigma);
    }
    sigmas
}

// ---------------------------------------------------------------------------
// Phase 4: Fuzzy simplicial set (directed → symmetrized undirected graph)
// ---------------------------------------------------------------------------

/// Build directed graph from KNN + sigmas, then symmetrize.
fn compute_fuzzy_graph(
    knn_indices: &[Vec<usize>],
    knn_dists: &[Vec<f64>],
    sigmas: &[f64],
    n_samples: usize,
    k: usize,
) -> Vec<Edge> {
    // Build directed weights: w(i→j) = exp(-(d(i,j) - rho_i) / sigma_i)
    // Store in a flat map: (i, j) → w_directed
    let mut directed: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_samples];

    for i in 0..n_samples {
        let rho = knn_dists[i][0].max(0.0);
        let sigma = sigmas[i];
        for idx in 0..k {
            let j = knn_indices[i][idx];
            let d = knn_dists[i][idx];
            let shifted = (d - rho).max(0.0);
            let w = if sigma > 0.0 {
                (-shifted / sigma).exp()
            } else {
                1.0
            };
            directed[i].push((j, w));
        }
    }

    // Symmetrize: w_sym(i,j) = w(i→j) + w(j→i) - w(i→j) * w(j→i)
    // Collect all unique edges
    use std::collections::BTreeMap;
    let mut pair_weights: BTreeMap<(usize, usize), (f64, f64)> = BTreeMap::new();

    for i in 0..n_samples {
        for &(j, w) in &directed[i] {
            let key = if i < j { (i, j) } else { (j, i) };
            let entry = pair_weights.entry(key).or_insert((0.0, 0.0));
            if i < j {
                entry.0 = w; // w(i→j)
            } else {
                entry.1 = w; // w(j→i) stored as second
            }
        }
    }

    // For pairs where only one direction exists, the other is 0
    let mut edges = Vec::with_capacity(pair_weights.len());
    for (&(i, j), &(w_ij, w_ji)) in &pair_weights {
        let sym = w_ij + w_ji - w_ij * w_ji;
        if sym > 0.0 {
            edges.push(Edge {
                i,
                j,
                weight: sym,
            });
        }
    }

    edges
}

// ---------------------------------------------------------------------------
// Phase 5: Fit a/b parameters
// ---------------------------------------------------------------------------

/// Fit the a, b parameters for the low-dimensional distance curve.
///
/// For default min_dist=0.1, spread=1.0, returns precomputed values.
/// Otherwise runs a simple Gauss-Newton fit.
fn fit_ab_params(min_dist: f64, spread: f64) -> (f64, f64) {
    // Precomputed defaults for the common case
    if (min_dist - 0.1).abs() < 1e-10 && (spread - 1.0).abs() < 1e-10 {
        return (1.929, 0.7915);
    }

    // Gauss-Newton least squares fit of:
    //   f(d) = 1 / (1 + a * d^(2b))
    // to the target piecewise function:
    //   target(d) = 1.0 if d <= min_dist
    //   target(d) = exp(-(d - min_dist) / spread) if d > min_dist

    let n_points = 300;
    let max_d = 3.0 * spread + min_dist;

    let mut a = 1.0;
    let mut b = 1.0;

    for _ in 0..100 {
        let mut jt_j_00 = 0.0;
        let mut jt_j_01 = 0.0;
        let mut jt_j_11 = 0.0;
        let mut jt_r_0 = 0.0;
        let mut jt_r_1 = 0.0;

        for k in 0..n_points {
            let d = (k as f64 + 0.5) * max_d / n_points as f64;
            let target = if d <= min_dist {
                1.0
            } else {
                (-(d - min_dist) / spread).exp()
            };

            let d2b = d.powf(2.0 * b);
            let denom = 1.0 + a * d2b;
            let f = 1.0 / denom;
            let residual = f - target;

            // Partial derivatives
            let df_da = -d2b / (denom * denom);
            let df_db = if d > 1e-12 {
                -a * d2b * 2.0 * d.ln() / (denom * denom)
            } else {
                0.0
            };

            // Accumulate J^T J and J^T r
            jt_j_00 += df_da * df_da;
            jt_j_01 += df_da * df_db;
            jt_j_11 += df_db * df_db;
            jt_r_0 += df_da * residual;
            jt_r_1 += df_db * residual;
        }

        // Solve 2x2 system (J^T J) * delta = -J^T r
        let det = jt_j_00 * jt_j_11 - jt_j_01 * jt_j_01;
        if det.abs() < 1e-20 {
            break;
        }
        let da = -(jt_j_11 * jt_r_0 - jt_j_01 * jt_r_1) / det;
        let db = -(jt_j_00 * jt_r_1 - jt_j_01 * jt_r_0) / det;

        a += da;
        b += db;

        // Clamp to reasonable range
        a = a.max(0.001);
        b = b.max(0.001);

        if da.abs() < 1e-8 && db.abs() < 1e-8 {
            break;
        }
    }

    (a, b)
}

// ---------------------------------------------------------------------------
// Phase 6: Initialize embedding
// ---------------------------------------------------------------------------

fn initialize_embedding(
    data: &[f64],
    n_features: usize,
    n_samples: usize,
    n_components: usize,
    init: UmapInit,
    seed: u64,
) -> Result<Vec<f64>> {
    match init {
        UmapInit::Pca => {
            let pca_config = crate::reduction::PcaConfig {
                n_components,
                max_iter: 100,
                tolerance: 1e-6,
            };
            let pca_result = crate::reduction::pca(data, n_features, &pca_config)?;
            let pca_dims = pca_result.n_components;
            // If PCA produced fewer components (n_features < n_components),
            // pad remaining dimensions with small random values
            let mut embedding = vec![0.0; n_samples * n_components];
            let mut rng = Xorshift64(seed.max(1));
            for i in 0..n_samples {
                for d in 0..n_components {
                    if d < pca_dims {
                        embedding[i * n_components + d] =
                            pca_result.transformed[i * pca_dims + d];
                    } else {
                        embedding[i * n_components + d] = rng.next_f64() * 0.01 - 0.005;
                    }
                }
            }
            // Scale PCA embedding to small range for SGD stability
            let max_abs = embedding
                .iter()
                .map(|x| x.abs())
                .fold(0.0_f64, f64::max);
            if max_abs > 0.0 {
                let scale = 10.0 / max_abs;
                for v in embedding.iter_mut() {
                    *v *= scale;
                }
            }
            Ok(embedding)
        }
        UmapInit::Random => {
            let mut rng = Xorshift64(seed.max(1));
            let mut embedding = vec![0.0; n_samples * n_components];
            for val in embedding.iter_mut() {
                *val = rng.next_f64() * 20.0 - 10.0;
            }
            Ok(embedding)
        }
    }
}

// ---------------------------------------------------------------------------
// Phase 7: SGD optimization
// ---------------------------------------------------------------------------

fn optimize_embedding(
    embedding: &mut [f64],
    edges: &[Edge],
    n_samples: usize,
    n_components: usize,
    a: f64,
    b: f64,
    initial_lr: f64,
    n_epochs: usize,
    negative_sample_rate: usize,
    seed: u64,
) {
    if edges.is_empty() || n_epochs == 0 {
        return;
    }

    // Compute epochs_per_sample: how often each edge is sampled
    let max_weight = edges
        .iter()
        .map(|e| e.weight)
        .fold(0.0_f64, f64::max);
    if max_weight <= 0.0 {
        return;
    }

    let epochs_per_sample: Vec<f64> = edges
        .iter()
        .map(|e| {
            let w = e.weight / max_weight;
            if w > 0.0 {
                1.0 / w
            } else {
                n_epochs as f64 + 1.0
            }
        })
        .collect();

    let mut epochs_of_next_sample: Vec<f64> = epochs_per_sample.clone();
    let mut rng = Xorshift64(seed.max(1));

    let clip = 4.0;

    for epoch in 0..n_epochs {
        let lr = initial_lr * (1.0 - epoch as f64 / n_epochs as f64);

        for (edge_idx, edge) in edges.iter().enumerate() {
            if epochs_of_next_sample[edge_idx] > epoch as f64 {
                continue;
            }

            let i = edge.i;
            let j = edge.j;

            // Compute squared distance in embedding space
            let mut dist_sq = 0.0;
            for d in 0..n_components {
                let diff = embedding[i * n_components + d] - embedding[j * n_components + d];
                dist_sq += diff * diff;
            }
            let dist_sq = dist_sq.max(1e-10);

            // Attractive force
            let grad_coeff = -2.0 * a * b * dist_sq.powf(b - 1.0)
                / (1.0 + a * dist_sq.powf(b));

            for d in 0..n_components {
                let diff = embedding[i * n_components + d] - embedding[j * n_components + d];
                let grad = grad_coeff * diff;
                let grad = grad.clamp(-clip, clip);
                embedding[i * n_components + d] += lr * grad;
                embedding[j * n_components + d] -= lr * grad;
            }

            // Repulsive forces (negative sampling)
            for _ in 0..negative_sample_rate {
                let neg = rng.next_bounded(n_samples);
                if neg == i {
                    continue;
                }

                let mut neg_dist_sq = 0.0;
                for d in 0..n_components {
                    let diff =
                        embedding[i * n_components + d] - embedding[neg * n_components + d];
                    neg_dist_sq += diff * diff;
                }
                let neg_dist_sq = neg_dist_sq.max(1e-10);

                let grad_coeff = 2.0 * b
                    / ((0.001 + neg_dist_sq) * (1.0 + a * neg_dist_sq.powf(b)));

                for d in 0..n_components {
                    let diff =
                        embedding[i * n_components + d] - embedding[neg * n_components + d];
                    let grad = grad_coeff * diff;
                    let grad = grad.clamp(-clip, clip);
                    embedding[i * n_components + d] += lr * grad;
                }
            }

            epochs_of_next_sample[edge_idx] += epochs_per_sample[edge_idx];
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn umap_basic_clusters() {
        // Two well-separated clusters in 3D
        let mut data = Vec::new();
        // Cluster A: near origin
        for i in 0..8 {
            data.push(i as f64 * 0.1);
            data.push(i as f64 * 0.05);
            data.push(0.0);
        }
        // Cluster B: far away
        for i in 0..8 {
            data.push(50.0 + i as f64 * 0.1);
            data.push(50.0 + i as f64 * 0.05);
            data.push(50.0);
        }
        let config = UmapConfig {
            n_neighbors: 5,
            n_epochs: 100,
            ..Default::default()
        };
        let result = umap(&data, 3, &config).unwrap();
        assert_eq!(result.n_samples, 16);
        assert_eq!(result.n_components, 2);
        assert_eq!(result.embedding.len(), 32);

        // Clusters should remain separated: centroid distance should be nonzero
        let mut centroid_a = [0.0; 2];
        let mut centroid_b = [0.0; 2];
        for i in 0..8 {
            centroid_a[0] += result.embedding[i * 2];
            centroid_a[1] += result.embedding[i * 2 + 1];
        }
        for i in 8..16 {
            centroid_b[0] += result.embedding[i * 2];
            centroid_b[1] += result.embedding[i * 2 + 1];
        }
        centroid_a[0] /= 8.0;
        centroid_a[1] /= 8.0;
        centroid_b[0] /= 8.0;
        centroid_b[1] /= 8.0;
        let dist = ((centroid_a[0] - centroid_b[0]).powi(2)
            + (centroid_a[1] - centroid_b[1]).powi(2))
        .sqrt();
        assert!(dist > 0.1, "clusters should be separated, got dist={dist}");
    }

    #[test]
    fn umap_3d_output() {
        let data = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0, 5.0, 6.0, 5.0, 5.0, 6.0, 6.0, 6.0,
        ];
        let config = UmapConfig {
            n_components: 3,
            n_neighbors: 3,
            n_epochs: 50,
            ..Default::default()
        };
        let result = umap(&data, 2, &config).unwrap();
        assert_eq!(result.n_components, 3);
        assert_eq!(result.embedding.len(), 8 * 3);
    }

    #[test]
    fn umap_empty_error() {
        let config = UmapConfig::default();
        assert!(umap(&[], 2, &config).is_err());
    }

    #[test]
    fn umap_dimension_mismatch() {
        let data = vec![1.0, 2.0, 3.0];
        let config = UmapConfig {
            n_neighbors: 2,
            ..Default::default()
        };
        // 3 not divisible by 2 → error
        assert!(umap(&data, 2, &config).is_err());
    }

    #[test]
    fn umap_n_neighbors_too_small() {
        let data = vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0];
        let config = UmapConfig {
            n_neighbors: 1,
            ..Default::default()
        };
        assert!(umap(&data, 2, &config).is_err());
    }

    #[test]
    fn umap_n_neighbors_too_large() {
        let data = vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0];
        let config = UmapConfig {
            n_neighbors: 3, // == n_samples
            ..Default::default()
        };
        assert!(umap(&data, 2, &config).is_err());
    }

    #[test]
    fn umap_summary() {
        let data = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0, 5.0, 6.0, 5.0,
        ];
        let config = UmapConfig {
            n_neighbors: 3,
            n_epochs: 50,
            ..Default::default()
        };
        let result = umap(&data, 2, &config).unwrap();
        let s = result.summary();
        assert!(s.contains("UMAP"));
        assert!(s.contains("2D"));
        assert!(s.contains("50 epochs"));
    }

    #[test]
    fn umap_random_init() {
        let data = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0, 5.0, 6.0, 5.0,
        ];
        let config = UmapConfig {
            n_neighbors: 3,
            n_epochs: 50,
            init: UmapInit::Random,
            ..Default::default()
        };
        let result = umap(&data, 2, &config).unwrap();
        assert_eq!(result.embedding.len(), 6 * 2);
        // Check all values are finite
        assert!(result.embedding.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn umap_pca_init() {
        let data = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0, 5.0, 6.0, 5.0,
        ];
        let config = UmapConfig {
            n_neighbors: 3,
            n_epochs: 50,
            init: UmapInit::Pca,
            ..Default::default()
        };
        let result = umap(&data, 2, &config).unwrap();
        assert_eq!(result.embedding.len(), 6 * 2);
        assert!(result.embedding.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn umap_deterministic() {
        let data = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0, 5.0, 6.0, 5.0,
        ];
        let config = UmapConfig {
            n_neighbors: 3,
            n_epochs: 50,
            seed: 123,
            ..Default::default()
        };
        let r1 = umap(&data, 2, &config).unwrap();
        let r2 = umap(&data, 2, &config).unwrap();
        assert_eq!(r1.embedding, r2.embedding);
    }

    #[test]
    fn umap_custom_metric() {
        let data = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0, 5.0, 6.0, 5.0,
        ];
        let config_manhattan = UmapConfig {
            n_neighbors: 3,
            n_epochs: 50,
            metric: DistanceMetric::Manhattan,
            ..Default::default()
        };
        let r = umap(&data, 2, &config_manhattan).unwrap();
        assert_eq!(r.embedding.len(), 6 * 2);
        assert!(r.embedding.iter().all(|v| v.is_finite()));

        let config_cosine = UmapConfig {
            n_neighbors: 3,
            n_epochs: 50,
            metric: DistanceMetric::Cosine,
            init: UmapInit::Random,
            ..Default::default()
        };
        let r = umap(&data, 2, &config_cosine).unwrap();
        assert_eq!(r.embedding.len(), 6 * 2);
        assert!(r.embedding.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn smooth_knn_convergence() {
        // Verify that smooth_knn_distances produces positive sigma values
        let knn_dists = vec![
            vec![0.5, 1.0, 1.5, 2.0, 3.0],
            vec![0.3, 0.8, 1.2, 1.8, 2.5],
            vec![0.1, 0.4, 0.9, 1.3, 2.0],
        ];
        let sigmas = smooth_knn_distances(&knn_dists, 3, 5);
        assert_eq!(sigmas.len(), 3);
        for &s in &sigmas {
            assert!(s > 0.0, "sigma should be positive, got {s}");
            assert!(s.is_finite(), "sigma should be finite");
        }
    }

    #[test]
    fn fit_ab_default() {
        let (a, b) = fit_ab_params(0.1, 1.0);
        assert!((a - 1.929).abs() < 0.001, "a={a}");
        assert!((b - 0.7915).abs() < 0.001, "b={b}");
    }

    #[test]
    fn fit_ab_custom() {
        let (a, b) = fit_ab_params(0.5, 1.0);
        assert!(a > 0.0, "a should be positive, got {a}");
        assert!(b > 0.0, "b should be positive, got {b}");
        assert!(a.is_finite());
        assert!(b.is_finite());
        // Different min_dist should give different a/b
        assert!((a - 1.929).abs() > 0.01 || (b - 0.7915).abs() > 0.01);
    }

    #[test]
    fn symmetrize_weights() {
        // Build a small KNN graph and verify symmetrization
        let knn_indices = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let knn_dists = vec![vec![1.0, 2.0], vec![1.0, 1.5], vec![2.0, 1.5]];
        let sigmas = smooth_knn_distances(&knn_dists, 3, 2);
        let edges = compute_fuzzy_graph(&knn_indices, &knn_dists, &sigmas, 3, 2);

        // All edges should have positive weight
        for edge in &edges {
            assert!(edge.weight > 0.0);
            assert!(edge.weight <= 1.0, "symmetrized weight should be <= 1.0");
        }

        // Should have undirected edges (i < j convention)
        for edge in &edges {
            assert!(edge.i < edge.j, "edges should be stored with i < j");
        }
    }
}
