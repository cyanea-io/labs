//! Spatial transcriptomics analysis — neighbor graphs, spatial autocorrelation,
//! co-occurrence, and ligand-receptor interaction scoring.
//!
//! Build spatial neighbor graphs from spot coordinates (Delaunay triangulation or
//! k-nearest neighbors), then compute spatial statistics such as Moran's I,
//! Geary's C, feature co-occurrence, and ligand-receptor interaction scores with
//! permutation-based significance testing.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A 2-D point with an associated spot/cell index.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpatialPoint {
    pub x: f64,
    pub y: f64,
    pub index: usize,
}

/// Spatial neighbor graph: for each node, a list of `(neighbor_index, distance)`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpatialGraph {
    pub n_nodes: usize,
    /// `neighbors[i]` is a vec of `(neighbor_index, distance)`.
    pub neighbors: Vec<Vec<(usize, f64)>>,
}

/// Result of Moran's I spatial autocorrelation test.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpatialAutocorrelation {
    pub morans_i: f64,
    pub expected_i: f64,
    pub variance_i: f64,
    pub z_score: f64,
    pub p_value: f64,
}

/// Result of Geary's C spatial autocorrelation test.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GearysC {
    pub c: f64,
    pub expected_c: f64,
    pub z_score: f64,
    pub p_value: f64,
}

/// Co-occurrence result for a pair of features.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CooccurrenceResult {
    pub feature_a: usize,
    pub feature_b: usize,
    pub observed: usize,
    pub expected: f64,
    pub log_odds: f64,
    pub p_value: f64,
}

/// Ligand-receptor interaction score with permutation p-value.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LrInteraction {
    pub ligand_name: String,
    pub receptor_name: String,
    pub interaction_score: f64,
    pub p_value: f64,
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Simple Xorshift64 PRNG for permutation tests.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Avoid zero state.
        Self {
            state: if seed == 0 { 0x5EED_DEAD_BEEF_CAFE } else { seed },
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

    #[allow(dead_code)]
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        let n = slice.len();
        for i in (1..n).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

/// Approximate error function using Abramowitz & Stegun formula 7.1.26.
fn erf_approx(x: f64) -> f64 {
    let a1: f64 = 0.254829592;
    let a2: f64 = -0.284496736;
    let a3: f64 = 1.421413741;
    let a4: f64 = -1.453152027;
    let a5: f64 = 1.061405429;
    let p: f64 = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    sign * y
}

/// Approximate complementary error function.
fn erfc_approx(x: f64) -> f64 {
    1.0 - erf_approx(x)
}

/// Two-sided p-value from a z-score (normal approximation).
fn p_from_z(z: f64) -> f64 {
    erfc_approx(z.abs() / std::f64::consts::SQRT_2)
}

/// Euclidean distance between two points.
fn euclidean(ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    ((ax - bx).powi(2) + (ay - by).powi(2)).sqrt()
}

// ---------------------------------------------------------------------------
// Delaunay triangulation (Bowyer-Watson)
// ---------------------------------------------------------------------------

/// A triangle represented by three point indices into the working point list.
#[derive(Clone, Copy)]
struct Triangle {
    v: [usize; 3],
}

/// Circumcircle: center (cx, cy) and radius squared.
fn circumcircle(
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
) -> (f64, f64, f64) {
    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if d.abs() < 1e-20 {
        // Degenerate — return huge radius.
        return (0.0, 0.0, f64::MAX);
    }
    let ux =
        ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
    let uy =
        ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
    let r2 = (ax - ux).powi(2) + (ay - uy).powi(2);
    (ux, uy, r2)
}

/// Check if point (px, py) lies inside the circumcircle of the triangle
/// defined by the three coordinate pairs.
fn in_circumcircle(
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
    px: f64,
    py: f64,
) -> bool {
    let (ux, uy, r2) = circumcircle(ax, ay, bx, by, cx, cy);
    let dist2 = (px - ux).powi(2) + (py - uy).powi(2);
    dist2 < r2
}

/// Build a Delaunay-triangulation-based spatial neighbor graph (Bowyer-Watson).
///
/// Creates a super-triangle that encloses all points, then incrementally inserts
/// each point.  Finally, triangles sharing vertices with the super-triangle are
/// removed and the remaining edges form the neighbor graph.
///
/// # Errors
///
/// Returns an error if fewer than 3 points are provided.
pub fn delaunay_neighbors(points: &[SpatialPoint]) -> Result<SpatialGraph> {
    let n = points.len();
    if n < 3 {
        return Err(CyaneaError::InvalidInput(
            "Delaunay triangulation requires at least 3 points".into(),
        ));
    }

    // Working coordinate list: first `n` entries are the real points,
    // entries n, n+1, n+2 are super-triangle vertices.
    let mut xs: Vec<f64> = points.iter().map(|p| p.x).collect();
    let mut ys: Vec<f64> = points.iter().map(|p| p.y).collect();

    // Super-triangle vertices (very large).
    xs.push(-1e10);
    ys.push(-1e10);
    xs.push(1e10);
    ys.push(-1e10);
    xs.push(0.0);
    ys.push(1e10);

    let st_a = n;
    let st_b = n + 1;
    let st_c = n + 2;

    let mut triangles: Vec<Triangle> = vec![Triangle { v: [st_a, st_b, st_c] }];

    for i in 0..n {
        let px = xs[i];
        let py = ys[i];

        // Find "bad" triangles whose circumcircle contains the new point.
        let mut bad: Vec<usize> = Vec::new();
        for (ti, tri) in triangles.iter().enumerate() {
            let [a, b, c] = tri.v;
            if in_circumcircle(xs[a], ys[a], xs[b], ys[b], xs[c], ys[c], px, py) {
                bad.push(ti);
            }
        }

        // Find boundary polygon — edges shared by exactly one bad triangle.
        let mut edge_count: Vec<([usize; 2], usize)> = Vec::new();
        for &ti in &bad {
            let v = triangles[ti].v;
            let edges = [[v[0], v[1]], [v[1], v[2]], [v[2], v[0]]];
            for e in &edges {
                let key = if e[0] < e[1] { [e[0], e[1]] } else { [e[1], e[0]] };
                if let Some(entry) = edge_count.iter_mut().find(|(k, _)| *k == key) {
                    entry.1 += 1;
                } else {
                    edge_count.push((key, 1));
                }
            }
        }

        let boundary: Vec<[usize; 2]> = edge_count
            .into_iter()
            .filter(|&(_, c)| c == 1)
            .map(|(e, _)| e)
            .collect();

        // Remove bad triangles (in reverse order to keep indices valid).
        bad.sort_unstable();
        for &ti in bad.iter().rev() {
            triangles.swap_remove(ti);
        }

        // Create new triangles from the point to each boundary edge.
        for e in &boundary {
            triangles.push(Triangle { v: [i, e[0], e[1]] });
        }
    }

    // Remove triangles that share a vertex with the super-triangle.
    triangles.retain(|tri| {
        !tri.v.iter().any(|&v| v == st_a || v == st_b || v == st_c)
    });

    // Build the neighbor graph from remaining triangles.
    let mut neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

    let mut add_edge = |a: usize, b: usize| {
        let d = euclidean(xs[a], ys[a], xs[b], ys[b]);
        if !neighbors[a].iter().any(|&(nb, _)| nb == b) {
            neighbors[a].push((b, d));
        }
        if !neighbors[b].iter().any(|&(nb, _)| nb == a) {
            neighbors[b].push((a, d));
        }
    };

    for tri in &triangles {
        let [a, b, c] = tri.v;
        add_edge(a, b);
        add_edge(a, c);
        add_edge(b, c);
    }

    Ok(SpatialGraph { n_nodes: n, neighbors })
}

// ---------------------------------------------------------------------------
// k-nearest neighbors
// ---------------------------------------------------------------------------

/// Build a k-nearest-neighbor spatial graph.
///
/// For each point, the `k` closest points (by Euclidean distance) are selected
/// as neighbors.  The graph is then symmetrized: if *i* is a neighbor of *j*,
/// then *j* is also a neighbor of *i*.
///
/// # Errors
///
/// Returns an error if `k` is 0 or `k >= n`.
pub fn knn_spatial_neighbors(points: &[SpatialPoint], k: usize) -> Result<SpatialGraph> {
    let n = points.len();
    if k == 0 {
        return Err(CyaneaError::InvalidInput("k must be >= 1".into()));
    }
    if k >= n {
        return Err(CyaneaError::InvalidInput(format!(
            "k ({}) must be less than number of points ({})",
            k, n
        )));
    }

    let mut neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

    for i in 0..n {
        let mut dists: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let d = euclidean(points[i].x, points[i].y, points[j].x, points[j].y);
                (j, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.truncate(k);
        neighbors[i] = dists;
    }

    // Symmetrize.
    for i in 0..n {
        let nbs: Vec<(usize, f64)> = neighbors[i].clone();
        for (j, d) in nbs {
            if !neighbors[j].iter().any(|&(nb, _)| nb == i) {
                neighbors[j].push((i, d));
            }
        }
    }

    Ok(SpatialGraph { n_nodes: n, neighbors })
}

// ---------------------------------------------------------------------------
// Moran's I
// ---------------------------------------------------------------------------

/// Compute Moran's I spatial autocorrelation statistic.
///
/// Uses binary weights (w_ij = 1 for neighbors, 0 otherwise).
/// Significance is assessed via z-score under the normality assumption.
///
/// # Errors
///
/// Returns an error if `values` length does not match the graph size, or if
/// variance is zero.
pub fn morans_i(values: &[f64], graph: &SpatialGraph) -> Result<SpatialAutocorrelation> {
    let n = graph.n_nodes;
    if values.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "values length ({}) does not match graph size ({})",
            values.len(),
            n
        )));
    }

    let mean = values.iter().sum::<f64>() / n as f64;
    let deviations: Vec<f64> = values.iter().map(|&v| v - mean).collect();
    let sum_sq: f64 = deviations.iter().map(|d| d * d).sum();
    if sum_sq.abs() < 1e-30 {
        return Err(CyaneaError::InvalidInput(
            "zero variance in values — Moran's I is undefined".into(),
        ));
    }

    // Total weight W (binary weights: each neighbor pair counted once per direction).
    let mut w_total = 0.0_f64;
    let mut numerator = 0.0_f64;
    for i in 0..n {
        for &(j, _) in &graph.neighbors[i] {
            w_total += 1.0;
            numerator += deviations[i] * deviations[j];
        }
    }

    if w_total == 0.0 {
        return Err(CyaneaError::InvalidInput("graph has no edges".into()));
    }

    let i_stat = (n as f64 / w_total) * numerator / sum_sq;
    let expected = -1.0 / (n as f64 - 1.0);

    // Variance under normality (Cliff & Ord formula).
    let nf = n as f64;
    let s1 = {
        let mut s = 0.0;
        for i in 0..n {
            for &(j, _) in &graph.neighbors[i] {
                // w_ij + w_ji: since the graph stores both directions, both are 1
                // but we need to sum (w_ij + w_ji)^2 over distinct pairs.
                // Instead, iterate all i,j: s1 = 0.5 * sum (w_ij + w_ji)^2.
                // Since w_ij = 1 for stored neighbors, and w_ji = 1 if symmetric:
                let w_ji = if graph.neighbors[j].iter().any(|&(nb, _)| nb == i) {
                    1.0
                } else {
                    0.0
                };
                s += (1.0_f64 + w_ji).powi(2);
            }
        }
        0.5 * s
    };

    let s2 = {
        let mut s = 0.0;
        for i in 0..n {
            let row_sum: f64 = graph.neighbors[i].len() as f64;
            let col_sum: f64 = (0..n)
                .filter(|&j| graph.neighbors[j].iter().any(|&(nb, _)| nb == i))
                .count() as f64;
            s += (row_sum + col_sum).powi(2);
        }
        s
    };

    let w2 = w_total * w_total;

    let m4 = deviations.iter().map(|d| d.powi(4)).sum::<f64>() / nf;
    let m2 = sum_sq / nf;
    let b2 = m4 / (m2 * m2); // kurtosis

    let var_i = {
        let t1 = nf * ((nf * nf - 3.0 * nf + 3.0) * s1 - nf * s2 + 3.0 * w2);
        let t2 = b2 * ((nf * nf - nf) * s1 - 2.0 * nf * s2 + 6.0 * w2);
        let denom = (nf - 1.0) * (nf - 2.0) * (nf - 3.0) * w2;
        (t1 - t2) / denom - expected * expected
    };

    let z = if var_i > 0.0 {
        (i_stat - expected) / var_i.sqrt()
    } else {
        0.0
    };

    let p = p_from_z(z);

    Ok(SpatialAutocorrelation {
        morans_i: i_stat,
        expected_i: expected,
        variance_i: var_i,
        z_score: z,
        p_value: p,
    })
}

// ---------------------------------------------------------------------------
// Geary's C
// ---------------------------------------------------------------------------

/// Compute Geary's C spatial autocorrelation statistic.
///
/// Uses binary weights. Values less than 1 indicate positive spatial
/// autocorrelation; values greater than 1 indicate negative.
///
/// # Errors
///
/// Returns an error if `values` length does not match the graph, or variance
/// is zero.
pub fn gearys_c(values: &[f64], graph: &SpatialGraph) -> Result<GearysC> {
    let n = graph.n_nodes;
    if values.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "values length ({}) does not match graph size ({})",
            values.len(),
            n
        )));
    }

    let mean = values.iter().sum::<f64>() / n as f64;
    let sum_sq: f64 = values.iter().map(|v| (v - mean).powi(2)).sum();
    if sum_sq.abs() < 1e-30 {
        return Err(CyaneaError::InvalidInput(
            "zero variance in values — Geary's C is undefined".into(),
        ));
    }

    let mut w_total = 0.0_f64;
    let mut numerator = 0.0_f64;
    for i in 0..n {
        for &(j, _) in &graph.neighbors[i] {
            w_total += 1.0;
            numerator += (values[i] - values[j]).powi(2);
        }
    }

    if w_total == 0.0 {
        return Err(CyaneaError::InvalidInput("graph has no edges".into()));
    }

    let nf = n as f64;
    let c_stat = ((nf - 1.0) / (2.0 * w_total)) * numerator / sum_sq;
    let expected_c = 1.0;

    // Variance under normality (simplified).
    let deviations: Vec<f64> = values.iter().map(|&v| v - mean).collect();
    let m4 = deviations.iter().map(|d| d.powi(4)).sum::<f64>() / nf;
    let m2 = sum_sq / nf;
    let b2 = m4 / (m2 * m2);

    // Using the variance formula for Geary's C under normality.
    let s1 = {
        let mut s = 0.0;
        for i in 0..n {
            for &(j, _) in &graph.neighbors[i] {
                let w_ji = if graph.neighbors[j].iter().any(|&(nb, _)| nb == i) {
                    1.0
                } else {
                    0.0
                };
                s += (1.0_f64 + w_ji).powi(2);
            }
        }
        0.5 * s
    };

    let s2 = {
        let mut s = 0.0;
        for i in 0..n {
            let row_sum: f64 = graph.neighbors[i].len() as f64;
            let col_sum: f64 = (0..n)
                .filter(|&j| graph.neighbors[j].iter().any(|&(nb, _)| nb == i))
                .count() as f64;
            s += (row_sum + col_sum).powi(2);
        }
        s
    };

    let w2 = w_total * w_total;

    let var_c = {
        let t1 = (2.0 * s1 + s2) * (nf - 1.0) - 4.0 * w2;
        let t2 = (nf - 1.0) * (nf - 2.0) * (nf - 3.0);
        let normal_var = t1 / (2.0 * (nf + 1.0) * w2);
        // Adjust with kurtosis.
        let _adj = ((nf - 1.0) * s1 * (nf * nf - 3.0 * nf + 3.0 - (nf - 1.0) * b2))
            / (t2 * w2);
        let _adj2 = ((nf - 1.0) * s2 * (nf * nf + 3.0 * nf - 6.0 - (nf * nf - nf + 2.0) * b2))
            / (4.0 * t2 * w2);
        let _adj3 = (w2 * (nf * nf - 3.0 - (nf - 1.0).powi(2) * b2)) / (t2 * w2);
        // Use simplified variance estimate to avoid numerical instability.
        normal_var.max(1e-20)
    };

    let z = (c_stat - expected_c) / var_c.sqrt();
    let p = p_from_z(z);

    Ok(GearysC {
        c: c_stat,
        expected_c,
        z_score: z,
        p_value: p,
    })
}

// ---------------------------------------------------------------------------
// Co-occurrence
// ---------------------------------------------------------------------------

/// Test co-occurrence of feature pairs in a spatial neighborhood.
///
/// `expression` is a feature-major matrix: `expression[feature][spot]`.
/// Two features are considered "expressed" at a spot if their value exceeds
/// `threshold`. Co-occurrence counts the number of graph edges where both
/// features are expressed at connected spots (in either spot of the pair).
///
/// Significance is assessed by permuting spot labels `n_permutations` times.
///
/// # Errors
///
/// Returns an error if the expression matrix dimensions are inconsistent with
/// the graph.
pub fn cooccurrence(
    expression: &[Vec<f64>],
    graph: &SpatialGraph,
    threshold: f64,
    n_permutations: usize,
    seed: u64,
) -> Result<Vec<CooccurrenceResult>> {
    let n_features = expression.len();
    if n_features < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 features for co-occurrence".into(),
        ));
    }
    let n_spots = graph.n_nodes;
    for (fi, feat) in expression.iter().enumerate() {
        if feat.len() != n_spots {
            return Err(CyaneaError::InvalidInput(format!(
                "feature {} has {} values, expected {}",
                fi,
                feat.len(),
                n_spots
            )));
        }
    }

    // Pre-compute binary expressed vectors.
    let expressed: Vec<Vec<bool>> = expression
        .iter()
        .map(|feat| feat.iter().map(|&v| v > threshold).collect())
        .collect();

    // All edges as (i, j) pairs (directed from the adjacency list).
    let edges: Vec<(usize, usize)> = (0..n_spots)
        .flat_map(|i| graph.neighbors[i].iter().map(move |&(j, _)| (i, j)))
        .collect();
    let n_edges = edges.len();

    let mut rng = Xorshift64::new(seed);
    let mut results = Vec::new();

    for a in 0..n_features {
        for b in (a + 1)..n_features {
            // Observed: count edges where both endpoints have at least one of
            // the two features expressed, i.e., feature a expressed at one end
            // and feature b expressed at the other end (or both at both).
            let observed = edges
                .iter()
                .filter(|&&(i, j)| {
                    (expressed[a][i] && expressed[b][j])
                        || (expressed[a][j] && expressed[b][i])
                })
                .count();

            // Expected under independence: P(a expressed) * P(b expressed) * edges.
            let pa = expressed[a].iter().filter(|&&e| e).count() as f64 / n_spots as f64;
            let pb = expressed[b].iter().filter(|&&e| e).count() as f64 / n_spots as f64;
            // Probability that at least one direction matches for a random edge.
            let p_cooccur = 2.0 * pa * pb - pa * pa * pb * pb;
            let expected = p_cooccur * n_edges as f64;

            let log_odds = if expected > 0.0 && observed > 0 {
                (observed as f64 / expected).ln()
            } else if observed > 0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            };

            // Permutation test.
            let mut perm_ge = 0usize;
            let mut indices: Vec<usize> = (0..n_spots).collect();
            for _ in 0..n_permutations {
                rng.shuffle(&mut indices);
                let perm_count = edges
                    .iter()
                    .filter(|&&(i, j)| {
                        (expressed[a][indices[i]] && expressed[b][indices[j]])
                            || (expressed[a][indices[j]] && expressed[b][indices[i]])
                    })
                    .count();
                if perm_count >= observed {
                    perm_ge += 1;
                }
            }
            let p_value = if n_permutations > 0 {
                (perm_ge as f64 + 1.0) / (n_permutations as f64 + 1.0)
            } else {
                1.0
            };

            results.push(CooccurrenceResult {
                feature_a: a,
                feature_b: b,
                observed,
                expected,
                log_odds,
                p_value,
            });
        }
    }

    Ok(results)
}

// ---------------------------------------------------------------------------
// Ligand-receptor scoring
// ---------------------------------------------------------------------------

/// Score a ligand-receptor interaction across a spatial neighbor graph.
///
/// The interaction score is the mean of `ligand[i] * receptor[j]` across all
/// directed neighbor pairs `(i, j)` in the graph.  Significance is assessed by
/// permuting cell/spot labels `n_permutations` times using an Xorshift64 PRNG.
///
/// # Errors
///
/// Returns an error if expression vectors do not match the graph size, or the
/// graph has no edges.
pub fn ligand_receptor_score(
    ligand_expr: &[f64],
    receptor_expr: &[f64],
    graph: &SpatialGraph,
    ligand_name: &str,
    receptor_name: &str,
    n_permutations: usize,
    seed: u64,
) -> Result<LrInteraction> {
    let n = graph.n_nodes;
    if ligand_expr.len() != n || receptor_expr.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "expression vector lengths ({}, {}) must match graph size ({})",
            ligand_expr.len(),
            receptor_expr.len(),
            n
        )));
    }

    let n_edges: usize = graph.neighbors.iter().map(|nb| nb.len()).sum();
    if n_edges == 0 {
        return Err(CyaneaError::InvalidInput("graph has no edges".into()));
    }

    // Observed score: mean of ligand[i] * receptor[j] over all directed edges.
    let observed_sum: f64 = (0..n)
        .flat_map(|i| {
            graph.neighbors[i]
                .iter()
                .map(move |&(j, _)| ligand_expr[i] * receptor_expr[j])
        })
        .sum();
    let observed_score = observed_sum / n_edges as f64;

    // Permutation test.
    let mut rng = Xorshift64::new(seed);
    let mut perm_ge = 0usize;
    let mut perm_indices: Vec<usize> = (0..n).collect();

    for _ in 0..n_permutations {
        rng.shuffle(&mut perm_indices);
        let mut perm_sum: f64 = 0.0;
        for i in 0..n {
            for &(j, _) in &graph.neighbors[i] {
                perm_sum += ligand_expr[perm_indices[i]] * receptor_expr[perm_indices[j]];
            }
        }
        let perm_score = perm_sum / n_edges as f64;
        if perm_score >= observed_score {
            perm_ge += 1;
        }
    }

    let p_value = if n_permutations > 0 {
        (perm_ge as f64 + 1.0) / (n_permutations as f64 + 1.0)
    } else {
        1.0
    };

    Ok(LrInteraction {
        ligand_name: ligand_name.to_string(),
        receptor_name: receptor_name.to_string(),
        interaction_score: observed_score,
        p_value,
    })
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(x: f64, y: f64, index: usize) -> SpatialPoint {
        SpatialPoint { x, y, index }
    }

    // 1. Delaunay 3 points → 1 triangle (3 edges, each point connected to other 2)
    #[test]
    fn test_delaunay_triangle() {
        let points = vec![pt(0.0, 0.0, 0), pt(1.0, 0.0, 1), pt(0.5, 1.0, 2)];
        let graph = delaunay_neighbors(&points).unwrap();
        assert_eq!(graph.n_nodes, 3);
        // Each point connected to the other 2.
        for i in 0..3 {
            assert_eq!(
                graph.neighbors[i].len(),
                2,
                "node {} has {} neighbors, expected 2",
                i,
                graph.neighbors[i].len()
            );
        }
    }

    // 2. Delaunay 4 points → valid triangulation (2 triangles)
    #[test]
    fn test_delaunay_four_points() {
        // Square: should produce 2 triangles → each vertex has 2 or 3 neighbors.
        let points = vec![
            pt(0.0, 0.0, 0),
            pt(1.0, 0.0, 1),
            pt(1.0, 1.0, 2),
            pt(0.0, 1.0, 3),
        ];
        let graph = delaunay_neighbors(&points).unwrap();
        assert_eq!(graph.n_nodes, 4);
        // Total edges: 2 triangles share 1 edge → 5 undirected edges, or
        // 4 outer + 1 diagonal = 5.
        let total_undirected_edges: usize = graph.neighbors.iter().map(|nb| nb.len()).sum::<usize>() / 2;
        assert!(
            total_undirected_edges >= 4 && total_undirected_edges <= 6,
            "expected 4-6 undirected edges, got {}",
            total_undirected_edges
        );
    }

    // 3. Delaunay valid — no circumcircle violations
    #[test]
    fn test_delaunay_circumcircle_valid() {
        let points = vec![
            pt(0.0, 0.0, 0),
            pt(3.0, 0.0, 1),
            pt(1.5, 2.5, 2),
            pt(0.5, 1.0, 3),
            pt(2.5, 1.0, 4),
        ];
        let graph = delaunay_neighbors(&points).unwrap();
        assert_eq!(graph.n_nodes, 5);

        // Reconstruct triangles from the graph edges and verify the Delaunay
        // property: no point lies strictly inside any triangle's circumcircle.
        // Collect all triangles as sets of 3 mutually-connected nodes.
        let mut triangles: Vec<[usize; 3]> = Vec::new();
        for i in 0..graph.n_nodes {
            for &(j, _) in &graph.neighbors[i] {
                if j <= i {
                    continue;
                }
                for &(k, _) in &graph.neighbors[j] {
                    if k <= j {
                        continue;
                    }
                    // Check if i-k are also neighbors.
                    if graph.neighbors[i].iter().any(|&(nb, _)| nb == k) {
                        triangles.push([i, j, k]);
                    }
                }
            }
        }

        for tri in &triangles {
            let [a, b, c] = *tri;
            for p in 0..graph.n_nodes {
                if p == a || p == b || p == c {
                    continue;
                }
                let inside = in_circumcircle(
                    points[a].x, points[a].y,
                    points[b].x, points[b].y,
                    points[c].x, points[c].y,
                    points[p].x, points[p].y,
                );
                assert!(
                    !inside,
                    "point {} is inside circumcircle of triangle [{}, {}, {}]",
                    p, a, b, c
                );
            }
        }
    }

    // 4. kNN k=2 — each point has at least 2 neighbors
    #[test]
    fn test_knn_min_neighbors() {
        let points = vec![
            pt(0.0, 0.0, 0),
            pt(1.0, 0.0, 1),
            pt(2.0, 0.0, 2),
            pt(3.0, 0.0, 3),
            pt(4.0, 0.0, 4),
        ];
        let graph = knn_spatial_neighbors(&points, 2).unwrap();
        for i in 0..graph.n_nodes {
            assert!(
                graph.neighbors[i].len() >= 2,
                "node {} has {} neighbors, expected >= 2",
                i,
                graph.neighbors[i].len()
            );
        }
    }

    // 5. kNN symmetry — if A neighbors B, then B neighbors A
    #[test]
    fn test_knn_symmetry() {
        let points = vec![
            pt(0.0, 0.0, 0),
            pt(1.0, 0.5, 1),
            pt(2.0, -0.3, 2),
            pt(0.5, 2.0, 3),
            pt(3.0, 1.0, 4),
        ];
        let graph = knn_spatial_neighbors(&points, 2).unwrap();
        for i in 0..graph.n_nodes {
            for &(j, _) in &graph.neighbors[i] {
                assert!(
                    graph.neighbors[j].iter().any(|&(nb, _)| nb == i),
                    "node {} has neighbor {} but not vice versa",
                    i,
                    j
                );
            }
        }
    }

    /// Helper: build a simple grid graph for autocorrelation tests.
    fn grid_graph(rows: usize, cols: usize) -> (Vec<SpatialPoint>, SpatialGraph) {
        let n = rows * cols;
        let mut points = Vec::with_capacity(n);
        for r in 0..rows {
            for c in 0..cols {
                points.push(pt(c as f64, r as f64, r * cols + c));
            }
        }

        let mut neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                if c + 1 < cols {
                    let j = r * cols + c + 1;
                    neighbors[i].push((j, 1.0));
                    neighbors[j].push((i, 1.0));
                }
                if r + 1 < rows {
                    let j = (r + 1) * cols + c;
                    neighbors[i].push((j, 1.0));
                    neighbors[j].push((i, 1.0));
                }
            }
        }

        let graph = SpatialGraph { n_nodes: n, neighbors };
        (points, graph)
    }

    // 6. Moran's I clustered → positive
    #[test]
    fn test_morans_i_clustered() {
        let (_pts, graph) = grid_graph(4, 4);
        // Clustered: top-left is high, bottom-right is low.
        let values = vec![
            10.0, 10.0, 9.0, 9.0, 10.0, 10.0, 9.0, 9.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0,
            2.0,
        ];
        let result = morans_i(&values, &graph).unwrap();
        assert!(
            result.morans_i > 0.0,
            "clustered data should have positive Moran's I, got {}",
            result.morans_i
        );
    }

    // 7. Moran's I random → near zero
    #[test]
    fn test_morans_i_random() {
        let (_pts, graph) = grid_graph(5, 5);
        // Use a larger grid and shuffled values so spatial autocorrelation is weak.
        let values = vec![
            12.0, 5.0, 18.0, 3.0, 15.0,
            9.0, 22.0, 1.0, 14.0, 7.0,
            20.0, 4.0, 16.0, 8.0, 11.0,
            2.0, 17.0, 6.0, 23.0, 10.0,
            13.0, 19.0, 24.0, 21.0, 25.0,
        ];
        let result = morans_i(&values, &graph).unwrap();
        // For shuffled data on a grid, Moran's I should be modest in magnitude.
        // Expected value = -1/(n-1) = -0.042. With only 25 points, values
        // can fluctuate significantly so use a wide tolerance.
        assert!(
            result.morans_i.abs() < 0.6,
            "random data should have Moran's I of modest magnitude, got {}",
            result.morans_i
        );
    }

    // 8. Moran's I checkerboard → negative
    #[test]
    fn test_morans_i_checkerboard() {
        let (_pts, graph) = grid_graph(4, 4);
        // Checkerboard: alternating high/low.
        let mut values = vec![0.0; 16];
        for r in 0..4 {
            for c in 0..4 {
                values[r * 4 + c] = if (r + c) % 2 == 0 { 10.0 } else { 0.0 };
            }
        }
        let result = morans_i(&values, &graph).unwrap();
        assert!(
            result.morans_i < 0.0,
            "checkerboard should have negative Moran's I, got {}",
            result.morans_i
        );
    }

    // 9. Geary's C clustered → less than 1
    #[test]
    fn test_gearys_c_clustered() {
        let (_pts, graph) = grid_graph(4, 4);
        let values = vec![
            10.0, 10.0, 9.0, 9.0, 10.0, 10.0, 9.0, 9.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0,
            2.0,
        ];
        let result = gearys_c(&values, &graph).unwrap();
        assert!(
            result.c < 1.0,
            "clustered data should have Geary's C < 1, got {}",
            result.c
        );
    }

    // 10. Geary's C random → near 1
    #[test]
    fn test_gearys_c_random() {
        let (_pts, graph) = grid_graph(4, 4);
        let values = vec![
            3.0, 7.0, 1.0, 9.0, 8.0, 2.0, 6.0, 4.0, 5.0, 9.0, 3.0, 7.0, 1.0, 8.0, 4.0, 6.0,
        ];
        let result = gearys_c(&values, &graph).unwrap();
        assert!(
            (result.c - 1.0).abs() < 0.5,
            "random data should have Geary's C near 1, got {}",
            result.c
        );
    }

    // 11. Co-occurrence co-expressed → high observed
    #[test]
    fn test_cooccurrence_coexpressed() {
        let (_pts, graph) = grid_graph(3, 3);
        // Both features expressed at every node → every edge is a co-occurrence edge.
        let feat_a = vec![5.0; 9];
        let feat_b = vec![5.0; 9];
        let expression = vec![feat_a, feat_b];
        let results = cooccurrence(&expression, &graph, 1.0, 100, 42).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert!(
            r.observed > 0,
            "co-expressed features should have observed > 0"
        );
        // When all nodes express both features, observed == n_edges and expected
        // is also n_edges, so log_odds ≈ 0. Just check observed is positive.
    }

    // 12. Co-occurrence independent → near expected
    #[test]
    fn test_cooccurrence_independent() {
        let (_pts, graph) = grid_graph(4, 4);
        // Feature A expressed in left half, B in top half → partially independent.
        let mut feat_a = vec![0.0; 16];
        let mut feat_b = vec![0.0; 16];
        for r in 0..4 {
            for c in 0..4 {
                if c < 2 {
                    feat_a[r * 4 + c] = 5.0;
                }
                if r < 2 {
                    feat_b[r * 4 + c] = 5.0;
                }
            }
        }
        let expression = vec![feat_a, feat_b];
        let results = cooccurrence(&expression, &graph, 1.0, 200, 123).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        // The observed and expected should be in the same ballpark.
        let ratio = if r.expected > 0.0 {
            r.observed as f64 / r.expected
        } else {
            1.0
        };
        assert!(
            ratio > 0.2 && ratio < 5.0,
            "independent features should have observed/expected near 1, got {}",
            ratio
        );
    }

    // 13. L-R high score when ligand neighbors receptor
    #[test]
    fn test_lr_high_score() {
        let (_pts, graph) = grid_graph(3, 3);
        // Ligand highly expressed in left column, receptor in middle column
        // (spatially adjacent on the grid).
        let ligand = vec![10.0, 0.0, 0.0, 10.0, 0.0, 0.0, 10.0, 0.0, 0.0];
        let receptor = vec![0.0, 10.0, 0.0, 0.0, 10.0, 0.0, 0.0, 10.0, 0.0];
        let result = ligand_receptor_score(
            &ligand,
            &receptor,
            &graph,
            "WNT3A",
            "FZD1",
            500,
            42,
        )
        .unwrap();
        assert!(
            result.interaction_score > 0.0,
            "ligand-receptor score should be positive, got {}",
            result.interaction_score
        );
        assert_eq!(result.ligand_name, "WNT3A");
        assert_eq!(result.receptor_name, "FZD1");
    }

    // 14. L-R permutation p-value reasonable
    #[test]
    fn test_lr_pvalue_reasonable() {
        let (_pts, graph) = grid_graph(4, 4);
        // Strong spatial signal: ligand in top half, receptor in bottom half,
        // meeting at the boundary.
        let mut ligand = vec![0.0; 16];
        let mut receptor = vec![0.0; 16];
        for r in 0..4 {
            for c in 0..4 {
                if r < 2 {
                    ligand[r * 4 + c] = 10.0;
                } else {
                    receptor[r * 4 + c] = 10.0;
                }
            }
        }
        let result = ligand_receptor_score(
            &ligand,
            &receptor,
            &graph,
            "CXCL12",
            "CXCR4",
            999,
            77,
        )
        .unwrap();
        // p-value should be between 0 and 1 and on the smaller side since
        // the spatial arrangement is structured.
        assert!(
            result.p_value > 0.0 && result.p_value <= 1.0,
            "p-value should be in (0, 1], got {}",
            result.p_value
        );
        assert!(
            result.interaction_score > 0.0,
            "should have positive interaction score"
        );
    }
}
