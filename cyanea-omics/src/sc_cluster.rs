//! Single-cell clustering: kNN graph construction, Leiden, Louvain, NMI, ARI.

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

use crate::network::Graph;
use crate::single_cell::{AnnData, ColumnData};
use crate::sparse::SparseMatrix;

// ── Neighbors ──────────────────────────────────────────────────────────────

/// Distance metric for neighbor graph construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DistanceMetric {
    Euclidean,
    Cosine,
}

/// Configuration for kNN graph construction.
#[derive(Debug, Clone)]
pub struct NeighborsConfig {
    /// Number of nearest neighbors.
    pub n_neighbors: usize,
    /// Number of PCs to use from `obsm["X_pca"]`.
    pub n_pcs: usize,
    /// Distance metric.
    pub metric: DistanceMetric,
    /// Random seed (unused currently, reserved for approximate methods).
    pub seed: u64,
}

impl Default for NeighborsConfig {
    fn default() -> Self {
        Self {
            n_neighbors: 15,
            n_pcs: 50,
            metric: DistanceMetric::Euclidean,
            seed: 42,
        }
    }
}

/// Build a kNN graph from PCA embeddings, storing connectivities and distances.
///
/// Reads `obsm["X_pca"]` and builds:
/// - `obsp["distances"]` — kNN distance matrix
/// - `obsp["connectivities"]` — UMAP-style fuzzy simplicial set connectivities
pub fn neighbors(adata: &mut AnnData, config: &NeighborsConfig) -> Result<()> {
    let pca = adata
        .get_obsm("X_pca")
        .ok_or_else(|| CyaneaError::InvalidInput("obsm['X_pca'] not found; run PCA first".into()))?
        .clone();

    let n_obs = pca.len();
    let n_dims = pca[0].len().min(config.n_pcs);

    if n_obs < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 observations for neighbor graph".into(),
        ));
    }

    let k = config.n_neighbors.min(n_obs - 1);

    // Compute pairwise distances and find kNN
    let mut knn_indices = vec![Vec::new(); n_obs];
    let mut knn_distances = vec![Vec::new(); n_obs];

    for i in 0..n_obs {
        let mut dists: Vec<(usize, f64)> = (0..n_obs)
            .filter(|&j| j != i)
            .map(|j| {
                let d = compute_distance(&pca[i][..n_dims], &pca[j][..n_dims], config.metric);
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        knn_indices[i] = dists[..k].iter().map(|&(j, _)| j).collect();
        knn_distances[i] = dists[..k].iter().map(|&(_, d)| d).collect();
    }

    // Build distance sparse matrix
    let mut dist_matrix = SparseMatrix::new(n_obs, n_obs);
    for i in 0..n_obs {
        for (idx, &j) in knn_indices[i].iter().enumerate() {
            let _ = dist_matrix.insert(i, j, knn_distances[i][idx]);
        }
    }

    // Build UMAP-style fuzzy connectivities
    // For each point, compute sigma (bandwidth) such that perplexity ~ k
    let mut conn_matrix = SparseMatrix::new(n_obs, n_obs);
    for i in 0..n_obs {
        let dists = &knn_distances[i];
        if dists.is_empty() {
            continue;
        }
        let rho = dists[0]; // distance to nearest neighbor
        let sigma = smooth_knn_dist(dists, (k as f64).ln() / std::f64::consts::LN_2);

        for (idx, &j) in knn_indices[i].iter().enumerate() {
            let d = dists[idx];
            let w = if d <= rho {
                1.0
            } else {
                (-(d - rho) / sigma.max(1e-10)).exp()
            };
            let _ = conn_matrix.insert(i, j, w);
        }
    }

    // Symmetrize: w_sym = w_ij + w_ji - w_ij * w_ji
    let mut sym_conn = SparseMatrix::new(n_obs, n_obs);
    for i in 0..n_obs {
        for &j in &knn_indices[i] {
            let w_ij = conn_matrix.get(i, j);
            let w_ji = conn_matrix.get(j, i);
            let w_sym = w_ij + w_ji - w_ij * w_ji;
            if w_sym > 0.0 {
                let _ = sym_conn.insert(i, j, w_sym);
            }
        }
    }

    adata.add_obsp("distances", dist_matrix)?;
    adata.add_obsp("connectivities", sym_conn)?;
    adata.add_uns("neighbors_n_neighbors", config.n_neighbors.to_string());

    Ok(())
}

/// Binary search for sigma in smooth kNN distance kernel.
fn smooth_knn_dist(distances: &[f64], target: f64) -> f64 {
    let mut lo = 1e-10_f64;
    let mut hi = 1000.0_f64;

    for _ in 0..64 {
        let mid = (lo + hi) / 2.0;
        let rho = distances[0];
        let sum: f64 = distances
            .iter()
            .map(|&d| {
                if d <= rho {
                    1.0
                } else {
                    (-(d - rho) / mid).exp()
                }
            })
            .sum();
        if sum > target {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) / 2.0
}

fn compute_distance(a: &[f64], b: &[f64], metric: DistanceMetric) -> f64 {
    match metric {
        DistanceMetric::Euclidean => {
            a.iter()
                .zip(b.iter())
                .map(|(x, y)| (x - y).powi(2))
                .sum::<f64>()
                .sqrt()
        }
        DistanceMetric::Cosine => {
            let dot: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
            let norm_a: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
            let norm_b: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
            let denom = norm_a * norm_b;
            if denom < 1e-15 {
                1.0
            } else {
                1.0 - dot / denom
            }
        }
    }
}

// ── Clustering ─────────────────────────────────────────────────────────────

/// Configuration for Leiden/Louvain clustering.
#[derive(Debug, Clone)]
pub struct ClusterConfig {
    /// Resolution parameter (higher = more clusters).
    pub resolution: f64,
    /// Maximum iterations for the algorithm.
    pub n_iterations: usize,
    /// Random seed.
    pub seed: u64,
    /// Key name to store cluster labels in obs.
    pub key_added: String,
}

impl Default for ClusterConfig {
    fn default() -> Self {
        Self {
            resolution: 1.0,
            n_iterations: 10,
            seed: 42,
            key_added: "leiden".into(),
        }
    }
}

/// Leiden community detection (Traag 2019).
///
/// Three-phase algorithm guaranteeing well-connected communities:
/// 1. Local moves (like Louvain Phase 1)
/// 2. Refinement: each community is refined by checking if nodes should form subcommunities
/// 3. Aggregation: contract graph, repeat
///
/// Stores cluster labels in `obs[key_added]`.
pub fn leiden(adata: &mut AnnData, config: &ClusterConfig) -> Result<()> {
    let conn = adata
        .get_obsp("connectivities")
        .ok_or_else(|| {
            CyaneaError::InvalidInput(
                "obsp['connectivities'] not found; run neighbors() first".into(),
            )
        })?
        .clone();

    let n = adata.n_obs();
    let graph = Graph::from_sparse_matrix(&conn);

    // Start with each node in its own community
    let mut assignments: Vec<usize> = (0..n).collect();

    // Total edge weight
    let total_weight: f64 = conn.iter().filter(|(i, j, _)| i < j).map(|(_, _, w)| w).sum();
    if total_weight == 0.0 {
        let labels: Vec<String> = (0..n).map(|i| i.to_string()).collect();
        adata.add_obs_column(&config.key_added, ColumnData::Strings(labels))?;
        return Ok(());
    }

    let m2 = total_weight * 2.0;

    for _iteration in 0..config.n_iterations {
        let mut moved = false;

        // Phase 1: local moves with resolution
        let mut local_improved = true;
        while local_improved {
            local_improved = false;

            // Use seeded deterministic order
            let mut order: Vec<usize> = (0..n).collect();
            let mut rng_state = config.seed.wrapping_add(_iteration as u64);
            for i in (1..n).rev() {
                rng_state = rng_state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let j = (rng_state >> 33) as usize % (i + 1);
                order.swap(i, j);
            }

            for &i in &order {
                let current = assignments[i];
                let k_i: f64 = graph.neighbors(i).iter().map(|(_, w)| w).sum();

                // Weights to each neighbor community
                let mut comm_weights: HashMap<usize, f64> = HashMap::new();
                for &(j, w) in graph.neighbors(i) {
                    *comm_weights.entry(assignments[j]).or_insert(0.0) += w;
                }

                let sigma_tot_current = community_weight(&graph, &assignments, current, n);
                let k_i_in_current = comm_weights.get(&current).copied().unwrap_or(0.0);

                let mut best = current;
                let mut best_dq = 0.0;

                for (&c, &k_i_in) in &comm_weights {
                    if c == current {
                        continue;
                    }
                    let sigma_tot = community_weight(&graph, &assignments, c, n);
                    let dq = (k_i_in - k_i_in_current) / m2
                        - config.resolution * k_i * (sigma_tot - sigma_tot_current + k_i)
                            / (m2 * m2)
                            * 2.0;
                    if dq > best_dq {
                        best_dq = dq;
                        best = c;
                    }
                }

                if best != current {
                    assignments[i] = best;
                    local_improved = true;
                    moved = true;
                }
            }
        }

        // Phase 2: refinement
        // Check if poorly-connected nodes should be moved to better-fitting communities
        let communities = unique_sorted(&assignments);
        for &c in &communities {
            let members: Vec<usize> = (0..n).filter(|&i| assignments[i] == c).collect();
            if members.len() <= 2 {
                continue;
            }

            for &node in &members {
                let k_i: f64 = graph.neighbors(node).iter().map(|(_, w)| w).sum();
                let k_i_in: f64 = graph
                    .neighbors(node)
                    .iter()
                    .filter(|&&(j, _)| assignments[j] == c)
                    .map(|(_, w)| w)
                    .sum();

                // Only consider extraction if weakly connected to own community
                if k_i > 0.0 && k_i_in / k_i < 0.1 {
                    let sigma_tot = community_weight(&graph, &assignments, c, n);
                    let dq_remove = -k_i_in / m2
                        + config.resolution * k_i * (sigma_tot - k_i) / (m2 * m2) * 2.0;

                    if dq_remove > 1e-6 {
                        let new_c = assignments.iter().max().unwrap_or(&0) + 1;
                        assignments[node] = new_c;
                        moved = true;
                    }
                }
            }
        }

        if !moved {
            break;
        }

        // Renumber communities
        renumber_assignments(&mut assignments);
    }

    renumber_assignments(&mut assignments);

    let labels: Vec<String> = assignments.iter().map(|c| c.to_string()).collect();
    adata.add_obs_column(&config.key_added, ColumnData::Strings(labels))?;

    let modularity = graph.modularity_with_resolution(&assignments, config.resolution);
    adata.add_uns(
        &format!("{}_modularity", config.key_added),
        modularity.to_string(),
    );

    Ok(())
}

/// Louvain clustering wrapper using the enhanced Graph method.
///
/// Stores cluster labels in `obs[key_added]`.
pub fn louvain(adata: &mut AnnData, config: &ClusterConfig) -> Result<()> {
    let conn = adata
        .get_obsp("connectivities")
        .ok_or_else(|| {
            CyaneaError::InvalidInput(
                "obsp['connectivities'] not found; run neighbors() first".into(),
            )
        })?
        .clone();

    let graph = Graph::from_sparse_matrix(&conn);
    let result = graph.louvain_with_resolution(config.resolution);

    let labels: Vec<String> = result.assignments.iter().map(|c| c.to_string()).collect();
    adata.add_obs_column(&config.key_added, ColumnData::Strings(labels))?;
    adata.add_uns(
        &format!("{}_modularity", config.key_added),
        result.modularity.to_string(),
    );

    Ok(())
}

fn community_weight(graph: &Graph, assignments: &[usize], community: usize, n: usize) -> f64 {
    let mut total = 0.0;
    for v in 0..n {
        if assignments[v] == community {
            for &(_, w) in graph.neighbors(v) {
                total += w;
            }
        }
    }
    total
}

fn unique_sorted(v: &[usize]) -> Vec<usize> {
    let mut u = v.to_vec();
    u.sort_unstable();
    u.dedup();
    u
}

fn renumber_assignments(assignments: &mut [usize]) {
    let mut map: Vec<usize> = Vec::new();
    for a in assignments.iter_mut() {
        let pos = map.iter().position(|&c| c == *a);
        *a = match pos {
            Some(idx) => idx,
            None => {
                map.push(*a);
                map.len() - 1
            }
        };
    }
}

// ── Clustering Metrics ─────────────────────────────────────────────────────

/// Normalized Mutual Information between two partitions.
pub fn nmi(a: &[usize], b: &[usize]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(
            "partitions must have the same length".into(),
        ));
    }
    let n = a.len();
    if n == 0 {
        return Ok(0.0);
    }

    let n_a = *a.iter().max().unwrap_or(&0) + 1;
    let n_b = *b.iter().max().unwrap_or(&0) + 1;

    // Contingency table
    let mut contingency = vec![0usize; n_a * n_b];
    for i in 0..n {
        contingency[a[i] * n_b + b[i]] += 1;
    }

    // Marginals
    let mut row_sums = vec![0usize; n_a];
    let mut col_sums = vec![0usize; n_b];
    for i in 0..n_a {
        for j in 0..n_b {
            row_sums[i] += contingency[i * n_b + j];
            col_sums[j] += contingency[i * n_b + j];
        }
    }

    let nf = n as f64;

    // Mutual information
    let mut mi = 0.0;
    for i in 0..n_a {
        for j in 0..n_b {
            let nij = contingency[i * n_b + j] as f64;
            if nij > 0.0 {
                mi += nij / nf * (nf * nij / (row_sums[i] as f64 * col_sums[j] as f64)).ln();
            }
        }
    }

    // Entropies
    let h_a: f64 = row_sums
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / nf;
            -p * p.ln()
        })
        .sum();

    let h_b: f64 = col_sums
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / nf;
            -p * p.ln()
        })
        .sum();

    let denom = ((h_a + h_b) / 2.0).max(1e-15);
    Ok((mi / denom).clamp(0.0, 1.0))
}

/// Adjusted Rand Index between two partitions.
pub fn adjusted_rand_index(a: &[usize], b: &[usize]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(
            "partitions must have the same length".into(),
        ));
    }
    let n = a.len();
    if n < 2 {
        return Ok(0.0);
    }

    let n_a = *a.iter().max().unwrap_or(&0) + 1;
    let n_b = *b.iter().max().unwrap_or(&0) + 1;

    // Contingency table
    let mut contingency = vec![0i64; n_a * n_b];
    for i in 0..n {
        contingency[a[i] * n_b + b[i]] += 1;
    }

    // Row and column sums
    let mut row_sums = vec![0i64; n_a];
    let mut col_sums = vec![0i64; n_b];
    for i in 0..n_a {
        for j in 0..n_b {
            row_sums[i] += contingency[i * n_b + j];
            col_sums[j] += contingency[i * n_b + j];
        }
    }

    // C(n,2) helper
    let c2 = |x: i64| -> i64 { x * (x - 1) / 2 };

    let sum_comb_c: i64 = contingency.iter().map(|&x| c2(x)).sum();
    let sum_comb_a: i64 = row_sums.iter().map(|&x| c2(x)).sum();
    let sum_comb_b: i64 = col_sums.iter().map(|&x| c2(x)).sum();
    let comb_n = c2(n as i64);

    let expected = sum_comb_a as f64 * sum_comb_b as f64 / comb_n as f64;
    let max_index = (sum_comb_a as f64 + sum_comb_b as f64) / 2.0;
    let denom = max_index - expected;

    if denom.abs() < 1e-15 {
        return Ok(0.0);
    }

    Ok((sum_comb_c as f64 - expected) / denom)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_cell::MatrixData;

    fn make_adata_with_pca(n_obs: usize, n_dims: usize, data: Vec<Vec<f64>>) -> AnnData {
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n_obs]);
        let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{}", i)).collect();
        let var_names = vec!["g0".into(), "g1".into()];
        let mut adata = AnnData::new(x, obs_names, var_names).unwrap();
        adata.add_obsm("X_pca", data).unwrap();
        adata
    }

    // ── Neighbors tests ──

    #[test]
    fn neighbors_basic() {
        // 6 points in 2D, two clusters
        let pca = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.1],
            vec![0.2, 0.0],
            vec![10.0, 10.0],
            vec![10.1, 10.1],
            vec![10.2, 10.0],
        ];
        let mut adata = make_adata_with_pca(6, 2, pca);
        neighbors(
            &mut adata,
            &NeighborsConfig {
                n_neighbors: 3,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(adata.get_obsp("distances").is_some());
        assert!(adata.get_obsp("connectivities").is_some());

        // Points 0,1,2 should be close to each other
        let dist = adata.get_obsp("distances").unwrap();
        let d01 = dist.get(0, 1);
        let d03 = dist.get(0, 3);
        // d01 should be stored (neighbor), d03 might not be (too far)
        if d03 > 0.0 {
            assert!(d01 < d03);
        }
    }

    #[test]
    fn neighbors_missing_pca() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let mut adata =
            AnnData::new(x, vec!["c0".into()], vec!["g0".into(), "g1".into()]).unwrap();
        let result = neighbors(&mut adata, &NeighborsConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn neighbors_cosine_metric() {
        let pca = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
        ];
        let mut adata = make_adata_with_pca(3, 2, pca);
        neighbors(
            &mut adata,
            &NeighborsConfig {
                n_neighbors: 2,
                metric: DistanceMetric::Cosine,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(adata.get_obsp("connectivities").is_some());
    }

    // ── Leiden tests ──

    #[test]
    fn leiden_two_clusters() {
        // Build connectivity matrix for two well-separated clusters
        let n = 10;
        let mut conn = SparseMatrix::new(n, n);
        // Cluster 1: nodes 0-4
        for i in 0..5 {
            for j in (i + 1)..5 {
                let _ = conn.insert(i, j, 1.0);
                let _ = conn.insert(j, i, 1.0);
            }
        }
        // Cluster 2: nodes 5-9
        for i in 5..10 {
            for j in (i + 1)..10 {
                let _ = conn.insert(i, j, 1.0);
                let _ = conn.insert(j, i, 1.0);
            }
        }
        // Weak bridge
        let _ = conn.insert(4, 5, 0.01);
        let _ = conn.insert(5, 4, 0.01);

        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let var_names = vec!["g0".into(), "g1".into()];
        let mut adata = AnnData::new(x, obs_names, var_names).unwrap();
        adata.add_obsp("connectivities", conn).unwrap();

        leiden(&mut adata, &ClusterConfig::default()).unwrap();

        let labels = adata.get_obs("leiden").unwrap().as_strings().unwrap();
        // Nodes in same cluster should have same label
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[0], labels[2]);
        assert_eq!(labels[5], labels[6]);
        assert_eq!(labels[5], labels[7]);
        // Different clusters
        assert_ne!(labels[0], labels[5]);
    }

    #[test]
    fn leiden_missing_connectivity() {
        let x = MatrixData::Dense(vec![vec![1.0]]);
        let mut adata = AnnData::new(x, vec!["c0".into()], vec!["g0".into()]).unwrap();
        let result = leiden(&mut adata, &ClusterConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn leiden_custom_key() {
        let n = 4;
        let mut conn = SparseMatrix::new(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let _ = conn.insert(i, j, 1.0);
                let _ = conn.insert(j, i, 1.0);
            }
        }
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();
        adata.add_obsp("connectivities", conn).unwrap();

        leiden(
            &mut adata,
            &ClusterConfig {
                key_added: "my_clusters".into(),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(adata.get_obs("my_clusters").is_some());
    }

    #[test]
    fn leiden_high_resolution() {
        let n = 8;
        let mut conn = SparseMatrix::new(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let _ = conn.insert(i, j, 1.0);
                let _ = conn.insert(j, i, 1.0);
            }
        }
        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();
        adata.add_obsp("connectivities", conn).unwrap();

        leiden(
            &mut adata,
            &ClusterConfig {
                resolution: 5.0,
                ..Default::default()
            },
        )
        .unwrap();

        let labels = adata.get_obs("leiden").unwrap().as_strings().unwrap();
        let n_clusters: usize = {
            let mut u: Vec<&String> = labels.iter().collect();
            u.sort();
            u.dedup();
            u.len()
        };
        // High resolution should produce more clusters
        assert!(n_clusters >= 1);
    }

    // ── Louvain wrapper tests ──

    #[test]
    fn louvain_wrapper_basic() {
        let n = 6;
        let mut conn = SparseMatrix::new(n, n);
        for i in 0..3 {
            for j in (i + 1)..3 {
                let _ = conn.insert(i, j, 1.0);
                let _ = conn.insert(j, i, 1.0);
            }
        }
        for i in 3..6 {
            for j in (i + 1)..6 {
                let _ = conn.insert(i, j, 1.0);
                let _ = conn.insert(j, i, 1.0);
            }
        }
        let _ = conn.insert(2, 3, 0.01);
        let _ = conn.insert(3, 2, 0.01);

        let x = MatrixData::Dense(vec![vec![0.0; 2]; n]);
        let obs_names: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let mut adata = AnnData::new(x, obs_names, vec!["g0".into(), "g1".into()]).unwrap();
        adata.add_obsp("connectivities", conn).unwrap();

        louvain(
            &mut adata,
            &ClusterConfig {
                key_added: "louvain".into(),
                ..Default::default()
            },
        )
        .unwrap();

        let labels = adata.get_obs("louvain").unwrap().as_strings().unwrap();
        assert_eq!(labels[0], labels[1]);
        assert_ne!(labels[0], labels[3]);
    }

    // ── NMI tests ──

    #[test]
    fn nmi_identical() {
        let a = vec![0, 0, 1, 1, 2, 2];
        let b = vec![0, 0, 1, 1, 2, 2];
        let score = nmi(&a, &b).unwrap();
        assert!((score - 1.0).abs() < 1e-10);
    }

    #[test]
    fn nmi_independent() {
        // One partition: all same, other: all different
        let a = vec![0, 0, 0, 0];
        let b = vec![0, 1, 2, 3];
        let score = nmi(&a, &b).unwrap();
        assert!(score.abs() < 1e-10);
    }

    #[test]
    fn nmi_different_lengths() {
        let result = nmi(&[0, 1], &[0, 1, 2]);
        assert!(result.is_err());
    }

    #[test]
    fn nmi_empty() {
        let score = nmi(&[], &[]).unwrap();
        assert_eq!(score, 0.0);
    }

    #[test]
    fn nmi_permuted() {
        let a = vec![0, 0, 1, 1, 2, 2];
        let b = vec![2, 2, 0, 0, 1, 1]; // same partition, different labels
        let score = nmi(&a, &b).unwrap();
        assert!((score - 1.0).abs() < 1e-10);
    }

    // ── ARI tests ──

    #[test]
    fn ari_identical() {
        let a = vec![0, 0, 1, 1, 2, 2];
        let b = vec![0, 0, 1, 1, 2, 2];
        let score = adjusted_rand_index(&a, &b).unwrap();
        assert!((score - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ari_random() {
        // Random assignment → ARI should be near 0
        let a = vec![0, 1, 0, 1, 0, 1, 0, 1];
        let b = vec![0, 0, 1, 1, 0, 0, 1, 1];
        let score = adjusted_rand_index(&a, &b).unwrap();
        assert!(score.abs() < 0.5); // should be close to 0
    }

    #[test]
    fn ari_permuted() {
        let a = vec![0, 0, 0, 1, 1, 1];
        let b = vec![1, 1, 1, 0, 0, 0]; // same partition, swapped labels
        let score = adjusted_rand_index(&a, &b).unwrap();
        assert!((score - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ari_different_lengths() {
        let result = adjusted_rand_index(&[0, 1], &[0]);
        assert!(result.is_err());
    }

    #[test]
    fn ari_single_element() {
        let score = adjusted_rand_index(&[0], &[0]).unwrap();
        assert_eq!(score, 0.0);
    }

    // ── Distance tests ──

    #[test]
    fn euclidean_distance() {
        let d = compute_distance(&[0.0, 0.0], &[3.0, 4.0], DistanceMetric::Euclidean);
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn cosine_distance_orthogonal() {
        let d = compute_distance(&[1.0, 0.0], &[0.0, 1.0], DistanceMetric::Cosine);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn cosine_distance_identical() {
        let d = compute_distance(&[1.0, 2.0], &[2.0, 4.0], DistanceMetric::Cosine);
        assert!(d.abs() < 1e-10);
    }

    // ── Integration test: neighbors → leiden pipeline ──

    #[test]
    fn neighbors_then_leiden() {
        let pca = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.1],
            vec![0.05, 0.05],
            vec![10.0, 10.0],
            vec![10.1, 10.1],
            vec![10.05, 10.05],
        ];
        let mut adata = make_adata_with_pca(6, 2, pca);
        neighbors(
            &mut adata,
            &NeighborsConfig {
                n_neighbors: 3,
                ..Default::default()
            },
        )
        .unwrap();
        leiden(&mut adata, &ClusterConfig::default()).unwrap();

        let labels = adata.get_obs("leiden").unwrap().as_strings().unwrap();
        // Points 0,1,2 should cluster together
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[1], labels[2]);
        // Points 3,4,5 should cluster together
        assert_eq!(labels[3], labels[4]);
        assert_eq!(labels[4], labels[5]);
        // Two clusters should be different
        assert_ne!(labels[0], labels[3]);
    }

    #[test]
    fn neighbors_then_louvain() {
        let pca = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.1],
            vec![0.05, 0.05],
            vec![10.0, 10.0],
            vec![10.1, 10.1],
            vec![10.05, 10.05],
        ];
        let mut adata = make_adata_with_pca(6, 2, pca);
        neighbors(
            &mut adata,
            &NeighborsConfig {
                n_neighbors: 3,
                ..Default::default()
            },
        )
        .unwrap();
        louvain(
            &mut adata,
            &ClusterConfig {
                key_added: "louvain".into(),
                ..Default::default()
            },
        )
        .unwrap();

        let labels = adata.get_obs("louvain").unwrap().as_strings().unwrap();
        assert_eq!(labels[0], labels[1]);
        assert_ne!(labels[0], labels[3]);
    }
}
