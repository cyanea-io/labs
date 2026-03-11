//! Spatial domain detection — identifying tissue regions with coherent
//! transcriptomic profiles from spatial transcriptomics data.
//!
//! Implements Leiden-like clustering on a spatial+expression combined graph,
//! Hidden Markov Random Field (HMRF) smoothing, and SVG (spatially variable
//! gene) detection.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A spatial domain assignment.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpatialDomain {
    /// Domain ID (0-based).
    pub domain_id: usize,
    /// Indices of spots/cells in this domain.
    pub members: Vec<usize>,
    /// Mean expression profile (gene-major) for this domain.
    pub mean_expression: Vec<f64>,
}

/// Result of spatial domain detection.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DomainResult {
    /// Per-spot domain assignment.
    pub labels: Vec<usize>,
    /// Spatial domains with member lists and profiles.
    pub domains: Vec<SpatialDomain>,
    /// Number of domains detected.
    pub n_domains: usize,
}

/// A spatially variable gene with significance.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpatiallyVariableGene {
    /// Gene index.
    pub gene_idx: usize,
    /// Gene name.
    pub gene_name: String,
    /// Moran's I statistic.
    pub morans_i: f64,
    /// p-value (permutation or analytical).
    pub p_value: f64,
    /// Adjusted p-value (BH-corrected).
    pub adjusted_p_value: f64,
}

/// Parameters for spatial domain detection.
#[derive(Debug, Clone)]
pub struct DomainParams {
    /// Weight for spatial distance (0 = expression only, 1 = space only).
    pub spatial_weight: f64,
    /// Number of domains (if known). If `None`, determined automatically.
    pub n_domains: Option<usize>,
    /// Number of k-nearest neighbors for the spatial graph.
    pub n_neighbors: usize,
    /// Maximum iterations for HMRF optimization.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

impl Default for DomainParams {
    fn default() -> Self {
        Self {
            spatial_weight: 0.5,
            n_domains: None,
            n_neighbors: 6,
            max_iter: 50,
            tolerance: 1e-4,
        }
    }
}

// ---------------------------------------------------------------------------
// Spatial domain detection
// ---------------------------------------------------------------------------

/// Detect spatial domains using spatially-aware k-means.
///
/// Combines expression similarity and spatial proximity to cluster spots/cells
/// into coherent tissue domains.
///
/// `expression` is spot-major: `expression[spot][gene]`.
/// `coords` are `(x, y)` per spot.
///
/// The algorithm:
/// 1. Normalize expression and spatial coordinates to [0, 1].
/// 2. Build combined distance: `d = (1 - α) * expr_dist + α * spatial_dist`
///    where `α = spatial_weight`.
/// 3. Run k-means on the combined distances.
pub fn detect_domains(
    expression: &[Vec<f64>],
    coords: &[(f64, f64)],
    params: &DomainParams,
) -> Result<DomainResult> {
    let n_spots = expression.len();
    if n_spots == 0 {
        return Err(CyaneaError::InvalidInput("empty expression matrix".into()));
    }
    if coords.len() != n_spots {
        return Err(CyaneaError::InvalidInput(
            "coords length must match expression rows".into(),
        ));
    }
    let n_genes = expression[0].len();
    if n_genes == 0 {
        return Err(CyaneaError::InvalidInput("no genes".into()));
    }
    for row in expression {
        if row.len() != n_genes {
            return Err(CyaneaError::InvalidInput(
                "inconsistent gene counts".into(),
            ));
        }
    }

    let k = params.n_domains.unwrap_or_else(|| {
        // Heuristic: sqrt(n_spots / 2), clamped to [2, 20]
        let est = ((n_spots as f64) / 2.0).sqrt().round() as usize;
        est.max(2).min(20)
    });

    // Normalize expression: z-score per gene
    let mut expr_norm: Vec<Vec<f64>> = expression.to_vec();
    for g in 0..n_genes {
        let mean = expression.iter().map(|r| r[g]).sum::<f64>() / n_spots as f64;
        let var = expression.iter().map(|r| (r[g] - mean).powi(2)).sum::<f64>() / n_spots as f64;
        let std = var.sqrt().max(1e-10);
        for s in 0..n_spots {
            expr_norm[s][g] = (expression[s][g] - mean) / std;
        }
    }

    // Normalize coordinates to [0, 1]
    let x_min = coords.iter().map(|c| c.0).fold(f64::MAX, f64::min);
    let x_max = coords.iter().map(|c| c.0).fold(f64::MIN, f64::max);
    let y_min = coords.iter().map(|c| c.1).fold(f64::MAX, f64::min);
    let y_max = coords.iter().map(|c| c.1).fold(f64::MIN, f64::max);
    let x_range = (x_max - x_min).max(1e-10);
    let y_range = (y_max - y_min).max(1e-10);
    let coords_norm: Vec<(f64, f64)> = coords
        .iter()
        .map(|&(x, y)| ((x - x_min) / x_range, (y - y_min) / y_range))
        .collect();

    let alpha = params.spatial_weight;

    // Combined distance function
    let combined_dist = |a: usize, b: usize| -> f64 {
        // Expression distance (Euclidean on z-scored)
        let expr_d: f64 = expr_norm[a]
            .iter()
            .zip(expr_norm[b].iter())
            .map(|(va, vb)| (va - vb).powi(2))
            .sum::<f64>()
            .sqrt()
            / (n_genes as f64).sqrt();

        // Spatial distance (Euclidean on normalized coords)
        let sp_d = ((coords_norm[a].0 - coords_norm[b].0).powi(2)
            + (coords_norm[a].1 - coords_norm[b].1).powi(2))
        .sqrt();

        (1.0 - alpha) * expr_d + alpha * sp_d
    };

    // K-means initialization: spread seeds evenly
    let mut centroids_expr: Vec<Vec<f64>> = Vec::with_capacity(k);
    let mut centroids_coord: Vec<(f64, f64)> = Vec::with_capacity(k);
    let step = n_spots.max(1) / k.max(1);
    for ki in 0..k {
        let idx = (ki * step).min(n_spots - 1);
        centroids_expr.push(expr_norm[idx].clone());
        centroids_coord.push(coords_norm[idx]);
    }

    let mut labels = vec![0usize; n_spots];

    for _iter in 0..params.max_iter {
        // Assignment step
        let mut changed = false;
        for s in 0..n_spots {
            let mut best_k = 0;
            let mut best_d = f64::MAX;
            for ki in 0..k {
                let expr_d: f64 = expr_norm[s]
                    .iter()
                    .zip(centroids_expr[ki].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt()
                    / (n_genes as f64).sqrt();
                let sp_d = ((coords_norm[s].0 - centroids_coord[ki].0).powi(2)
                    + (coords_norm[s].1 - centroids_coord[ki].1).powi(2))
                .sqrt();
                let d = (1.0 - alpha) * expr_d + alpha * sp_d;
                if d < best_d {
                    best_d = d;
                    best_k = ki;
                }
            }
            if labels[s] != best_k {
                labels[s] = best_k;
                changed = true;
            }
        }

        if !changed {
            break;
        }

        // Update step
        for ki in 0..k {
            let members: Vec<usize> = (0..n_spots).filter(|&s| labels[s] == ki).collect();
            if members.is_empty() {
                continue;
            }
            let count = members.len() as f64;
            for g in 0..n_genes {
                centroids_expr[ki][g] = members.iter().map(|&s| expr_norm[s][g]).sum::<f64>() / count;
            }
            centroids_coord[ki] = (
                members.iter().map(|&s| coords_norm[s].0).sum::<f64>() / count,
                members.iter().map(|&s| coords_norm[s].1).sum::<f64>() / count,
            );
        }
    }

    // Build domain result
    let mut domains = Vec::new();
    for ki in 0..k {
        let members: Vec<usize> = (0..n_spots).filter(|&s| labels[s] == ki).collect();
        if members.is_empty() {
            continue;
        }
        let count = members.len() as f64;
        let mean_expr: Vec<f64> = (0..n_genes)
            .map(|g| members.iter().map(|&s| expression[s][g]).sum::<f64>() / count)
            .collect();
        domains.push(SpatialDomain {
            domain_id: ki,
            members,
            mean_expression: mean_expr,
        });
    }

    // Suppress unused combined_dist warning
    let _ = combined_dist;

    let n_domains = domains.len();
    Ok(DomainResult {
        labels,
        domains,
        n_domains,
    })
}

// ---------------------------------------------------------------------------
// HMRF smoothing
// ---------------------------------------------------------------------------

/// Smooth domain labels using a Hidden Markov Random Field (HMRF).
///
/// Iteratively refines labels by considering both expression fit and
/// spatial neighborhood consistency (Ising-like prior).
///
/// `labels` is the initial per-spot assignment (e.g. from k-means).
/// `expression` is spot-major `[spot][gene]`.
/// `neighbors` maps each spot to its neighbor indices.
/// `beta` controls spatial smoothness (higher = smoother).
/// `max_iter` is the maximum number of ICM iterations.
pub fn hmrf_smooth(
    labels: &[usize],
    expression: &[Vec<f64>],
    neighbors: &[Vec<usize>],
    beta: f64,
    max_iter: usize,
) -> Result<Vec<usize>> {
    let n = labels.len();
    if expression.len() != n || neighbors.len() != n {
        return Err(CyaneaError::InvalidInput(
            "labels, expression, and neighbors must have same length".into(),
        ));
    }
    if n == 0 {
        return Ok(Vec::new());
    }

    let n_genes = expression[0].len();
    let k = *labels.iter().max().unwrap_or(&0) + 1;

    let mut current_labels = labels.to_vec();

    for _iter in 0..max_iter {
        // Compute cluster means
        let mut means: Vec<Vec<f64>> = vec![vec![0.0; n_genes]; k];
        let mut counts = vec![0usize; k];
        for i in 0..n {
            let l = current_labels[i];
            counts[l] += 1;
            for g in 0..n_genes {
                means[l][g] += expression[i][g];
            }
        }
        for ki in 0..k {
            if counts[ki] > 0 {
                for g in 0..n_genes {
                    means[ki][g] /= counts[ki] as f64;
                }
            }
        }

        // ICM update
        let mut changed = false;
        for i in 0..n {
            let mut best_label = current_labels[i];
            let mut best_energy = f64::MAX;

            for ki in 0..k {
                // Expression cost (squared Euclidean distance to cluster mean)
                let expr_cost: f64 = expression[i]
                    .iter()
                    .zip(means[ki].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum();

                // Spatial cost (number of neighbors with different label)
                let spatial_cost: f64 = neighbors[i]
                    .iter()
                    .filter(|&&j| current_labels[j] != ki)
                    .count() as f64;

                let energy = expr_cost + beta * spatial_cost;
                if energy < best_energy {
                    best_energy = energy;
                    best_label = ki;
                }
            }

            if best_label != current_labels[i] {
                current_labels[i] = best_label;
                changed = true;
            }
        }

        if !changed {
            break;
        }
    }

    Ok(current_labels)
}

// ---------------------------------------------------------------------------
// Spatially variable gene detection
// ---------------------------------------------------------------------------

/// Detect spatially variable genes using Moran's I.
///
/// For each gene, computes Moran's I on the spatial neighbor graph and
/// assesses significance via permutation testing.
///
/// Returns genes sorted by adjusted p-value (ascending).
pub fn find_spatially_variable_genes(
    expression: &[Vec<f64>],
    gene_names: &[String],
    neighbors: &[Vec<usize>],
    n_permutations: usize,
    seed: u64,
) -> Result<Vec<SpatiallyVariableGene>> {
    let n_spots = expression.len();
    if n_spots == 0 {
        return Err(CyaneaError::InvalidInput("empty expression".into()));
    }
    let n_genes = expression[0].len();
    if gene_names.len() != n_genes {
        return Err(CyaneaError::InvalidInput(
            "gene_names length must match gene count".into(),
        ));
    }
    if neighbors.len() != n_spots {
        return Err(CyaneaError::InvalidInput(
            "neighbors length must match spot count".into(),
        ));
    }

    // Build edges for weight computation
    let mut w_total = 0.0_f64;
    for nb in neighbors {
        w_total += nb.len() as f64;
    }

    let mut rng_state = if seed == 0 { 0x5EED_DEAD_BEEF_CAFE } else { seed };
    let xorshift = |state: &mut u64| -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    };

    let mut results = Vec::with_capacity(n_genes);

    for g in 0..n_genes {
        let values: Vec<f64> = expression.iter().map(|r| r[g]).collect();
        let mean = values.iter().sum::<f64>() / n_spots as f64;
        let devs: Vec<f64> = values.iter().map(|&v| v - mean).collect();
        let sum_sq: f64 = devs.iter().map(|d| d * d).sum();

        if sum_sq < 1e-30 {
            // Zero-variance gene
            results.push(SpatiallyVariableGene {
                gene_idx: g,
                gene_name: gene_names[g].clone(),
                morans_i: 0.0,
                p_value: 1.0,
                adjusted_p_value: 1.0,
            });
            continue;
        }

        // Compute Moran's I
        let mut numerator = 0.0_f64;
        for i in 0..n_spots {
            for &j in &neighbors[i] {
                numerator += devs[i] * devs[j];
            }
        }
        let morans_i = (n_spots as f64 / w_total) * numerator / sum_sq;

        // Permutation test
        let mut perm_ge = 0usize;
        let mut perm_indices: Vec<usize> = (0..n_spots).collect();
        for _ in 0..n_permutations {
            // Fisher-Yates shuffle
            for i in (1..n_spots).rev() {
                let j = (xorshift(&mut rng_state) as usize) % (i + 1);
                perm_indices.swap(i, j);
            }
            let mut perm_num = 0.0_f64;
            for i in 0..n_spots {
                for &j in &neighbors[i] {
                    perm_num += devs[perm_indices[i]] * devs[perm_indices[j]];
                }
            }
            let perm_i = (n_spots as f64 / w_total) * perm_num / sum_sq;
            if perm_i >= morans_i {
                perm_ge += 1;
            }
        }

        let p_value = if n_permutations > 0 {
            (perm_ge as f64 + 1.0) / (n_permutations as f64 + 1.0)
        } else {
            1.0
        };

        results.push(SpatiallyVariableGene {
            gene_idx: g,
            gene_name: gene_names[g].clone(),
            morans_i,
            p_value,
            adjusted_p_value: p_value, // corrected below
        });
    }

    // BH correction
    let n_tests = results.len();
    if n_tests > 0 {
        // Sort by p-value
        let mut indices: Vec<usize> = (0..n_tests).collect();
        indices.sort_by(|&a, &b| {
            results[a]
                .p_value
                .partial_cmp(&results[b].p_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut adj_p = vec![1.0_f64; n_tests];
        for (rank, &idx) in indices.iter().enumerate() {
            let bh = results[idx].p_value * n_tests as f64 / (rank + 1) as f64;
            adj_p[idx] = bh.min(1.0);
        }

        // Enforce monotonicity (from highest rank down)
        for i in (0..n_tests - 1).rev() {
            let idx_curr = indices[i];
            let idx_next = indices[i + 1];
            if adj_p[idx_curr] > adj_p[idx_next] {
                adj_p[idx_curr] = adj_p[idx_next];
            }
        }

        for (i, adj) in adj_p.into_iter().enumerate() {
            results[i].adjusted_p_value = adj;
        }
    }

    // Sort by adjusted p-value
    results.sort_by(|a, b| {
        a.adjusted_p_value
            .partial_cmp(&b.adjusted_p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(results)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_domains_basic() {
        // Two clusters in expression and space
        let mut expression = Vec::new();
        let mut coords = Vec::new();
        // Cluster 1: high expr near origin
        for i in 0..10 {
            expression.push(vec![10.0, 0.0]);
            coords.push((i as f64, 0.0));
        }
        // Cluster 2: low expr far away
        for i in 0..10 {
            expression.push(vec![0.0, 10.0]);
            coords.push((100.0 + i as f64, 0.0));
        }

        let params = DomainParams {
            n_domains: Some(2),
            ..Default::default()
        };
        let result = detect_domains(&expression, &coords, &params).unwrap();
        assert_eq!(result.n_domains, 2);
        assert_eq!(result.labels.len(), 20);

        // Check that the two clusters are in different domains
        let label0 = result.labels[0];
        let label10 = result.labels[10];
        assert_ne!(label0, label10);
    }

    #[test]
    fn test_detect_domains_auto_k() {
        let expression = vec![vec![1.0, 2.0]; 50];
        let coords: Vec<(f64, f64)> = (0..50).map(|i| (i as f64, 0.0)).collect();
        let params = DomainParams::default(); // n_domains = None
        let result = detect_domains(&expression, &coords, &params).unwrap();
        assert!(result.n_domains >= 1);
        assert_eq!(result.labels.len(), 50);
    }

    #[test]
    fn test_hmrf_smooth() {
        // Two spatial clusters, initial labels mostly correct with some noise
        let expression = vec![
            vec![10.0, 0.0], // cluster 0
            vec![10.0, 0.0], // cluster 0
            vec![10.0, 0.0], // cluster 0 (mislabeled as 1)
            vec![0.0, 10.0], // cluster 1
            vec![0.0, 10.0], // cluster 1
        ];
        let labels = vec![0, 0, 1, 1, 1]; // spot 2 mislabeled
        let neighbors = vec![
            vec![1],       // 0 -> 1
            vec![0, 2],    // 1 -> 0, 2
            vec![1, 3],    // 2 -> 1, 3
            vec![2, 4],    // 3 -> 2, 4
            vec![3],       // 4 -> 3
        ];
        let smoothed = hmrf_smooth(&labels, &expression, &neighbors, 1.0, 20).unwrap();
        assert_eq!(smoothed.len(), 5);
        // With beta=1 and expression cost, spot 2 should flip to 0
        assert_eq!(smoothed[2], 0);
    }

    #[test]
    fn test_hmrf_no_change() {
        // Already optimal labels
        let expression = vec![vec![10.0], vec![10.0], vec![0.0], vec![0.0]];
        let labels = vec![0, 0, 1, 1];
        let neighbors = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let smoothed = hmrf_smooth(&labels, &expression, &neighbors, 1.0, 10).unwrap();
        assert_eq!(smoothed, labels);
    }

    #[test]
    fn test_svg_detection() {
        // Gene 0: spatially clustered (left high, right low)
        // Gene 1: random (no spatial pattern)
        let n = 20;
        let mut expression = Vec::new();
        let mut coords = Vec::new();
        for i in 0..n {
            let g0 = if i < 10 { 10.0 } else { 0.0 }; // spatial pattern
            let g1 = if i % 2 == 0 { 5.0 } else { 5.1 }; // no pattern
            expression.push(vec![g0, g1]);
            coords.push((i as f64, 0.0));
        }

        // Build neighbors (linear chain)
        let neighbors: Vec<Vec<usize>> = (0..n)
            .map(|i| {
                let mut nb = Vec::new();
                if i > 0 { nb.push(i - 1); }
                if i < n - 1 { nb.push(i + 1); }
                nb
            })
            .collect();

        let gene_names = vec!["SpatialGene".into(), "RandomGene".into()];
        let results = find_spatially_variable_genes(
            &expression,
            &gene_names,
            &neighbors,
            199,
            42,
        )
        .unwrap();

        assert_eq!(results.len(), 2);
        // SpatialGene should have higher Moran's I
        let spatial = results.iter().find(|r| r.gene_name == "SpatialGene").unwrap();
        let random = results.iter().find(|r| r.gene_name == "RandomGene").unwrap();
        assert!(spatial.morans_i > random.morans_i);
    }

    #[test]
    fn test_svg_zero_variance() {
        let expression = vec![vec![5.0, 0.0]; 10];
        let neighbors: Vec<Vec<usize>> = (0..10)
            .map(|i| {
                let mut nb = Vec::new();
                if i > 0 { nb.push(i - 1); }
                if i < 9 { nb.push(i + 1); }
                nb
            })
            .collect();
        let gene_names = vec!["Constant".into(), "Zero".into()];
        let results = find_spatially_variable_genes(
            &expression,
            &gene_names,
            &neighbors,
            50,
            1,
        )
        .unwrap();
        // Constant gene should have p_value = 1.0
        let constant = results.iter().find(|r| r.gene_name == "Constant").unwrap();
        assert!((constant.p_value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bh_correction() {
        // With small sample, check that adjusted p-values are >= raw p-values
        let expression = vec![
            vec![10.0, 0.0, 5.0],
            vec![10.0, 0.0, 3.0],
            vec![10.0, 0.0, 7.0],
            vec![0.0, 10.0, 4.0],
            vec![0.0, 10.0, 6.0],
            vec![0.0, 10.0, 5.0],
        ];
        let neighbors = vec![
            vec![1, 2],
            vec![0, 2],
            vec![0, 1, 3],
            vec![2, 4, 5],
            vec![3, 5],
            vec![3, 4],
        ];
        let gene_names = vec!["G1".into(), "G2".into(), "G3".into()];
        let results = find_spatially_variable_genes(
            &expression,
            &gene_names,
            &neighbors,
            99,
            42,
        )
        .unwrap();
        for r in &results {
            assert!(r.adjusted_p_value >= r.p_value - 1e-10);
            assert!(r.adjusted_p_value <= 1.0 + 1e-10);
        }
    }
}
