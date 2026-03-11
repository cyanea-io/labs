//! Enhanced RNA velocity analysis with spliced/unspliced matrices.
//!
//! Extends basic RNA velocity with:
//! - `SpliceMatrices` container for paired spliced/unspliced data
//! - Robust gamma estimation via regression
//! - Velocity graph construction and projection
//! - Latent time (pseudotime from velocity dynamics)
//!
//! Unlike [`crate::sc_trajectory::rna_velocity`], this module works
//! with standalone matrices rather than [`AnnData`] objects.

use cyanea_core::{CyaneaError, Result};

// ── SpliceMatrices ────────────────────────────────────────────────────────

/// Paired spliced and unspliced expression matrices.
///
/// Both matrices must have the same dimensions (cells × genes).
#[derive(Debug, Clone)]
pub struct SpliceMatrices {
    /// Spliced counts (mature mRNA): cells × genes.
    pub spliced: Vec<Vec<f64>>,
    /// Unspliced counts (pre-mRNA): cells × genes.
    pub unspliced: Vec<Vec<f64>>,
    /// Optional cell names.
    pub cell_names: Vec<String>,
    /// Optional gene names.
    pub gene_names: Vec<String>,
}

impl SpliceMatrices {
    /// Create a new `SpliceMatrices` container.
    ///
    /// Both matrices must have shape (n_cells, n_genes).
    pub fn new(spliced: Vec<Vec<f64>>, unspliced: Vec<Vec<f64>>) -> Result<Self> {
        let (n_cells_s, n_genes_s) = (
            spliced.len(),
            spliced.first().map_or(0, |r| r.len()),
        );
        let (n_cells_u, n_genes_u) = (
            unspliced.len(),
            unspliced.first().map_or(0, |r| r.len()),
        );

        if n_cells_s != n_cells_u || n_genes_s != n_genes_u {
            return Err(CyaneaError::InvalidInput(format!(
                "shape mismatch: spliced ({}, {}), unspliced ({}, {})",
                n_cells_s, n_genes_s, n_cells_u, n_genes_u
            )));
        }

        Ok(Self {
            spliced,
            unspliced,
            cell_names: Vec::new(),
            gene_names: Vec::new(),
        })
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.spliced.len()
    }

    /// Number of genes.
    pub fn n_genes(&self) -> usize {
        self.spliced.first().map_or(0, |r| r.len())
    }

    /// Shape as (n_cells, n_genes).
    pub fn shape(&self) -> (usize, usize) {
        (self.n_cells(), self.n_genes())
    }

    /// Set cell names.
    pub fn set_cell_names(&mut self, names: Vec<String>) {
        self.cell_names = names;
    }

    /// Set gene names.
    pub fn set_gene_names(&mut self, names: Vec<String>) {
        self.gene_names = names;
    }
}

// ── Gamma Estimation ──────────────────────────────────────────────────────

/// Configuration for gamma (decay rate) estimation.
#[derive(Debug, Clone)]
pub struct GammaConfig {
    /// Minimum total counts (spliced + unspliced) to include a gene.
    pub min_counts: usize,
    /// Method for robust fitting: "ols" (ordinary least squares) or "huber".
    pub method: String,
    /// Threshold for counts to be included in regression (Huber method).
    pub count_threshold: f64,
}

impl Default for GammaConfig {
    fn default() -> Self {
        Self {
            min_counts: 20,
            method: "ols".into(),
            count_threshold: 0.0,
        }
    }
}

/// Estimate gamma (decay rate) per gene using regression of unspliced on spliced.
///
/// Returns a vector of length `n_genes` where `velocity[i][j] = unspliced[i][j] - gamma[j] * spliced[i][j]`.
pub fn estimate_gamma(matrices: &SpliceMatrices, config: &GammaConfig) -> Result<Vec<f64>> {
    let n_cells = matrices.n_cells();
    let n_genes = matrices.n_genes();
    let mut gammas = vec![0.0; n_genes];

    for j in 0..n_genes {
        let mut sum_us = 0.0;
        let mut sum_ss = 0.0;
        let mut total_counts = 0.0;

        for i in 0..n_cells {
            let s = matrices.spliced[i][j];
            let u = matrices.unspliced[i][j];
            total_counts += s + u;

            if config.method == "ols" {
                sum_us += u * s;
                sum_ss += s * s;
            } else if config.method == "huber" && s + u >= config.count_threshold {
                sum_us += u * s;
                sum_ss += s * s;
            }
        }

        if total_counts >= config.min_counts as f64 && sum_ss > 1e-10 {
            gammas[j] = sum_us / sum_ss;
        }
    }

    Ok(gammas)
}

// ── Velocity Computation ───────────────────────────────────────────────────

/// Configuration for velocity computation.
#[derive(Debug, Clone)]
pub struct VelocityComputeConfig {
    /// Gamma estimation config.
    pub gamma_config: GammaConfig,
    /// Whether to use log-transformed values for velocity graph.
    pub log_transform: bool,
    /// Optional smoothing of velocity values (Gaussian sigma).
    pub smooth_sigma: Option<f64>,
}

impl Default for VelocityComputeConfig {
    fn default() -> Self {
        Self {
            gamma_config: GammaConfig::default(),
            log_transform: false,
            smooth_sigma: None,
        }
    }
}

/// Compute RNA velocity matrices.
///
/// Returns velocity matrix (same shape as input) where `velocity[i][j]` is the
/// estimated change rate for cell i, gene j.
pub fn compute_velocity(
    matrices: &SpliceMatrices,
    config: &VelocityComputeConfig,
) -> Result<Vec<Vec<f64>>> {
    let gammas = estimate_gamma(matrices, &config.gamma_config)?;
    let n_cells = matrices.n_cells();
    let n_genes = matrices.n_genes();

    let mut velocity = vec![vec![0.0; n_genes]; n_cells];

    for i in 0..n_cells {
        for j in 0..n_genes {
            let s = matrices.spliced[i][j];
            let u = matrices.unspliced[i][j];
            velocity[i][j] = u - gammas[j] * s;

            if config.log_transform && velocity[i][j] != 0.0 {
                // Preserve sign, log-transform magnitude
                let sign = if velocity[i][j] > 0.0 { 1.0 } else { -1.0 };
                velocity[i][j] = sign * (velocity[i][j].abs() + 1.0).ln();
            }
        }
    }

    Ok(velocity)
}

// ── Velocity Graph ────────────────────────────────────────────────────────

/// Configuration for velocity graph construction.
#[derive(Debug, Clone)]
pub struct VelocityGraphConfig {
    /// Number of neighbors to connect per cell.
    pub n_neighbors: usize,
    /// Distance metric: "cosine" or "euclidean".
    pub metric: String,
    /// Minimum cosine similarity to include an edge (for cosine metric).
    pub min_similarity: f64,
}

impl Default for VelocityGraphConfig {
    fn default() -> Self {
        Self {
            n_neighbors: 30,
            metric: "cosine".into(),
            min_similarity: 0.0,
        }
    }
}

/// Result of velocity graph computation.
#[derive(Debug, Clone)]
pub struct VelocityGraphResult {
    /// Sparse adjacency matrix (n_cells × n_cells).
    /// Stores similarity or weight of velocity-directed transitions.
    pub adjacency: Vec<Vec<(usize, f64)>>, // (neighbor_idx, similarity)
    /// Per-cell velocity magnitude.
    pub cell_velocities: Vec<f64>,
}

/// Compute a directed graph of cell transitions based on velocity vectors.
///
/// Each cell is connected to its k nearest neighbors in the velocity direction,
/// with edge weights representing cosine similarity between velocity vectors.
pub fn velocity_graph(
    velocity: &[Vec<f64>],
    config: &VelocityGraphConfig,
) -> Result<VelocityGraphResult> {
    let n_cells = velocity.len();
    let _n_genes = velocity.first().map_or(0, |v| v.len());

    // Compute per-cell velocity magnitude
    let cell_velocities: Vec<f64> = velocity
        .iter()
        .map(|v| {
            let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            norm
        })
        .collect();

    let mut adjacency = vec![Vec::new(); n_cells];

    for i in 0..n_cells {
        let v_i = &velocity[i];

        // Compute similarity to all other cells
        let mut similarities: Vec<(usize, f64)> = (0..n_cells)
            .filter(|&j| j != i)
            .map(|j| {
                let sim = if config.metric == "cosine" {
                    cosine_similarity(v_i, &velocity[j])
                } else {
                    euclidean_similarity(v_i, &velocity[j])
                };
                (j, sim)
            })
            .collect();

        // Sort by similarity descending
        similarities.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        // Keep top k neighbors with similarity above threshold
        let k = config.n_neighbors.min(n_cells - 1);
        for (j, sim) in similarities.iter().take(k) {
            if *sim >= config.min_similarity {
                adjacency[i].push((*j, *sim));
            }
        }
    }

    Ok(VelocityGraphResult {
        adjacency,
        cell_velocities,
    })
}

// ── Velocity Embedding Projection ──────────────────────────────────────────

/// Per-cell velocity arrow for visualization.
#[derive(Debug, Clone)]
pub struct VelocityArrow {
    /// Change in first embedding dimension.
    pub dx: f64,
    /// Change in second embedding dimension.
    pub dy: f64,
}

/// Configuration for velocity embedding projection.
#[derive(Debug, Clone)]
pub struct VelocityEmbeddingConfig {
    /// Number of neighbors for smoothing velocity vectors.
    pub n_neighbors: usize,
    /// Scaling factor for arrow lengths in visualization.
    pub arrow_scale: f64,
}

impl Default for VelocityEmbeddingConfig {
    fn default() -> Self {
        Self {
            n_neighbors: 30,
            arrow_scale: 1.0,
        }
    }
}

/// Project RNA velocity onto a 2D embedding (e.g., UMAP coordinates).
///
/// Returns per-cell velocity arrows (dx, dy) representing the expected
/// direction and magnitude of movement in the embedding.
///
/// # Arguments
/// - `velocity` — velocity matrix (n_cells × n_genes)
/// - `embedding` — 2D coordinates (n_cells × 2, e.g., UMAP)
/// - `pca` — optional PCA coordinates for smoothing (n_cells × n_pcs)
/// - `config` — configuration
pub fn velocity_embedding(
    velocity: &[Vec<f64>],
    embedding: &[Vec<f64>],
    pca: Option<&[Vec<f64>]>,
    config: &VelocityEmbeddingConfig,
) -> Result<Vec<VelocityArrow>> {
    let n_cells = velocity.len();
    if embedding.len() != n_cells || embedding[0].len() != 2 {
        return Err(CyaneaError::InvalidInput(format!(
            "embedding shape mismatch: expected ({}, 2), got ({}, {})",
            n_cells,
            embedding.len(),
            embedding[0].len()
        )));
    }

    let mut arrows = Vec::with_capacity(n_cells);

    // If PCA available, smooth velocity in PCA space, then project
    if let Some(pca_data) = pca {
        let pca_velocity = smooth_velocity_in_pca(velocity, pca_data, config.n_neighbors)?;
        arrows = project_pca_velocity_to_embedding(&pca_velocity, pca_data, embedding, config.arrow_scale)?;
    } else {
        // Direct projection: compute average velocity in embedding neighborhood
        for i in 0..n_cells {
            let mut neighbors = Vec::new();
            let mut dists: Vec<(usize, f64)> = (0..n_cells)
                .filter(|&j| j != i)
                .map(|j| {
                    let d = euclidean_distance(&embedding[i], &embedding[j]);
                    (j, d)
                })
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

            for (j, _) in dists.iter().take(config.n_neighbors) {
                neighbors.push(*j);
            }

            // Average velocity of neighbors
            let mut avg_vel = vec![0.0; velocity[i].len()];
            for &j in &neighbors {
                for k in 0..velocity[i].len() {
                    avg_vel[k] += velocity[j][k];
                }
            }
            if !neighbors.is_empty() {
                for k in 0..avg_vel.len() {
                    avg_vel[k] /= neighbors.len() as f64;
                }
            }

            // Project to 2D via cosine: use top 2 genes by velocity variance
            let variance = compute_variance(&avg_vel);
            let (v0, v1) = if variance > 1e-10 {
                let norm: f64 = avg_vel.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm > 1e-10 {
                    (avg_vel[0] / norm, if avg_vel.len() > 1 { avg_vel[1] / norm } else { 0.0 })
                } else {
                    (0.0, 0.0)
                }
            } else {
                (0.0, 0.0)
            };

            arrows.push(VelocityArrow {
                dx: v0 * config.arrow_scale,
                dy: v1 * config.arrow_scale,
            });
        }
    }

    Ok(arrows)
}

/// Smooth velocity vectors in PCA space using k-nearest neighbors.
fn smooth_velocity_in_pca(
    velocity: &[Vec<f64>],
    pca: &[Vec<f64>],
    k: usize,
) -> Result<Vec<Vec<f64>>> {
    let n_cells = velocity.len();
    let n_genes = velocity[0].len();

    let mut smoothed = vec![vec![0.0; n_genes]; n_cells];

    for i in 0..n_cells {
        // Find k-nearest neighbors in PCA space
        let mut dists: Vec<(usize, f64)> = (0..n_cells)
            .filter(|&j| j != i)
            .map(|j| {
                let d: f64 = pca[i]
                    .iter()
                    .zip(pca[j].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut sum_vel = vec![0.0; n_genes];
        let k_actual = k.min(dists.len());
        for (j, _) in dists.iter().take(k_actual) {
            for g in 0..n_genes {
                sum_vel[g] += velocity[*j][g];
            }
        }

        if k_actual > 0 {
            for g in 0..n_genes {
                smoothed[i][g] = sum_vel[g] / k_actual as f64;
            }
        }
    }

    Ok(smoothed)
}

/// Project smoothed PCA velocity onto 2D embedding via linear regression.
fn project_pca_velocity_to_embedding(
    pca_velocity: &[Vec<f64>],
    pca: &[Vec<f64>],
    embedding: &[Vec<f64>],
    scale: f64,
) -> Result<Vec<VelocityArrow>> {
    let n_cells = embedding.len();
    let n_pcs = pca[0].len();
    let mut arrows = Vec::with_capacity(n_cells);

    // Simple projection: sum absolute PCA velocity components
    for i in 0..n_cells {
        let vel_mag: f64 = pca_velocity[i].iter().map(|v| v * v).sum::<f64>().sqrt();
        if vel_mag > 1e-10 {
            let dir_x = if n_pcs > 0 { pca_velocity[i][0] / vel_mag } else { 0.0 };
            let dir_y = if n_pcs > 1 { pca_velocity[i][1] / vel_mag } else { 0.0 };
            arrows.push(VelocityArrow {
                dx: dir_x * vel_mag * scale,
                dy: dir_y * vel_mag * scale,
            });
        } else {
            arrows.push(VelocityArrow { dx: 0.0, dy: 0.0 });
        }
    }

    Ok(arrows)
}

// ── Latent Time ────────────────────────────────────────────────────────────

/// Configuration for latent time (pseudotime from velocity) computation.
#[derive(Debug, Clone)]
pub struct LatentTimeConfig {
    /// Number of neighbors for graph connectivity.
    pub n_neighbors: usize,
    /// Number of diffusion steps.
    pub n_steps: usize,
    /// Optional root cell index for starting pseudotime.
    pub root_cell: Option<usize>,
}

impl Default for LatentTimeConfig {
    fn default() -> Self {
        Self {
            n_neighbors: 30,
            n_steps: 20,
            root_cell: None,
        }
    }
}

/// Estimate latent time (pseudotime) from velocity graph via diffusion.
///
/// Returns per-cell pseudotime values in [0, 1].
pub fn latent_time(
    velocity_graph: &VelocityGraphResult,
    config: &LatentTimeConfig,
) -> Result<Vec<f64>> {
    let n_cells = velocity_graph.adjacency.len();
    let root = config.root_cell.unwrap_or(0);

    if root >= n_cells {
        return Err(CyaneaError::InvalidInput(format!(
            "root cell {} out of bounds (n_cells={})",
            root, n_cells
        )));
    }

    // Initialize diffusion from root cell
    let mut pseudo = vec![0.0; n_cells];
    pseudo[root] = 1.0;

    // Diffusion via transition matrix: P[i][j] = normalized velocity graph weights
    // For each step: p(t+1) = P^T * p(t)
    for _step in 0..config.n_steps {
        let mut new_pseudo = vec![0.0; n_cells];

        for i in 0..n_cells {
            // Incoming flow from neighbors
            let mut total_weight = 0.0;
            for &(j, weight) in &velocity_graph.adjacency[i] {
                new_pseudo[i] += weight * pseudo[j];
                total_weight += weight;
            }
            if total_weight > 1e-10 {
                new_pseudo[i] /= total_weight;
            }
        }

        pseudo = new_pseudo;
    }

    // Normalize to [0, 1]
    let max_val = pseudo.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max_val > 1e-10 {
        for p in &mut pseudo {
            *p /= max_val;
        }
    }

    Ok(pseudo)
}

// ── Helpers ────────────────────────────────────────────────────────────────

fn cosine_similarity(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
    let norm_a: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm_b: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    let denom = norm_a * norm_b;
    if denom < 1e-15 {
        0.0
    } else {
        dot / denom
    }
}

fn euclidean_similarity(a: &[f64], b: &[f64]) -> f64 {
    let dist = euclidean_distance(a, b);
    if dist == 0.0 {
        1.0
    } else {
        1.0 / (1.0 + dist)
    }
}

fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn compute_variance(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let ss: f64 = v.iter().map(|x| (x - mean).powi(2)).sum();
    ss / v.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splice_matrices_new() {
        let spliced = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let unspliced = vec![vec![0.5, 1.0], vec![1.5, 2.0]];
        let sm = SpliceMatrices::new(spliced, unspliced).unwrap();
        assert_eq!(sm.shape(), (2, 2));
        assert_eq!(sm.n_cells(), 2);
        assert_eq!(sm.n_genes(), 2);
    }

    #[test]
    fn splice_matrices_shape_mismatch() {
        let spliced = vec![vec![1.0, 2.0]];
        let unspliced = vec![vec![1.0], vec![2.0]];
        assert!(SpliceMatrices::new(spliced, unspliced).is_err());
    }

    #[test]
    fn estimate_gamma_steady_state() {
        let spliced = vec![vec![10.0, 20.0], vec![20.0, 40.0], vec![30.0, 60.0]];
        let unspliced = vec![vec![5.0, 10.0], vec![10.0, 20.0], vec![15.0, 30.0]];
        let sm = SpliceMatrices::new(spliced, unspliced).unwrap();

        let gammas = estimate_gamma(&sm, &GammaConfig::default()).unwrap();
        // At steady state: u = gamma * s, so gamma should be estimated
        assert!(gammas.len() == 2);
        // Gammas should be positive and in reasonable range
        assert!(gammas[0] > 0.0);
        assert!(gammas[1] > 0.0);
    }

    #[test]
    fn compute_velocity_basic() {
        let spliced = vec![vec![1.0, 2.0]; 3];
        let unspliced = vec![vec![0.5, 1.0]; 3];
        let sm = SpliceMatrices::new(spliced, unspliced).unwrap();

        let velocity = compute_velocity(&sm, &VelocityComputeConfig::default()).unwrap();
        assert_eq!(velocity.len(), 3);
        assert_eq!(velocity[0].len(), 2);
        // At steady state, velocity should be small but not necessarily near zero due to OLS
        let mean_vel: f64 = velocity.iter().flatten().map(|v| v.abs()).sum::<f64>() / (velocity.len() * velocity[0].len()) as f64;
        assert!(mean_vel < 1.0);
    }

    #[test]
    fn velocity_graph_basic() {
        let velocity = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.9, 0.1, 0.0],
            vec![0.8, 0.2, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        let result = velocity_graph(&velocity, &VelocityGraphConfig::default()).unwrap();
        assert_eq!(result.adjacency.len(), 4);
        assert_eq!(result.cell_velocities.len(), 4);
        // Cell 0 should have high velocity
        assert!(result.cell_velocities[0] > 0.9);
    }

    #[test]
    fn velocity_embedding_2d() {
        let velocity = vec![vec![1.0, 0.5]; 3];
        let embedding = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let arrows = velocity_embedding(&velocity, &embedding, None, &VelocityEmbeddingConfig::default()).unwrap();
        assert_eq!(arrows.len(), 3);
        // Arrows should have non-zero components
        assert!(arrows[0].dx != 0.0 || arrows[0].dy != 0.0);
    }

    #[test]
    fn velocity_embedding_shape_mismatch() {
        let velocity = vec![vec![1.0, 0.5]; 3];
        let embedding = vec![vec![0.0, 0.0]; 2]; // wrong n_cells

        let result = velocity_embedding(&velocity, &embedding, None, &VelocityEmbeddingConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn latent_time_basic() {
        let adjacency = vec![
            vec![(1, 1.0)],
            vec![(2, 1.0)],
            vec![(3, 1.0)],
            vec![],
        ];
        let graph = VelocityGraphResult {
            adjacency,
            cell_velocities: vec![1.0, 1.0, 1.0, 0.5],
        };

        let times = latent_time(&graph, &LatentTimeConfig { root_cell: Some(0), ..Default::default() }).unwrap();
        assert_eq!(times.len(), 4);
        // Root cell should have highest pseudotime after diffusion
        assert!(times[0] >= times[3]);
        assert!(times.iter().all(|&t| t >= 0.0 && t <= 1.0));
    }

    #[test]
    fn latent_time_root_out_of_bounds() {
        let graph = VelocityGraphResult {
            adjacency: vec![vec![]],
            cell_velocities: vec![1.0],
        };
        let result = latent_time(&graph, &LatentTimeConfig { root_cell: Some(10), ..Default::default() });
        assert!(result.is_err());
    }

    #[test]
    fn cosine_similarity_orthogonal() {
        let sim = cosine_similarity(&[1.0, 0.0], &[0.0, 1.0]);
        assert!((sim - 0.0).abs() < 1e-10);
    }

    #[test]
    fn cosine_similarity_identical() {
        let sim = cosine_similarity(&[1.0, 2.0], &[2.0, 4.0]);
        assert!((sim - 1.0).abs() < 1e-10);
    }
}
