//! Clustering algorithms: k-means, DBSCAN, and hierarchical.

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::distance::{compute_distance, DistanceMatrix, DistanceMetric};

// ---------------------------------------------------------------------------
// K-Means
// ---------------------------------------------------------------------------

/// Configuration for k-means clustering.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KMeansConfig {
    pub n_clusters: usize,
    pub max_iter: usize,
    pub tolerance: f64,
    pub seed: u64,
}

impl Default for KMeansConfig {
    fn default() -> Self {
        Self {
            n_clusters: 2,
            max_iter: 300,
            tolerance: 1e-4,
            seed: 42,
        }
    }
}

/// Result of k-means clustering.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KMeansResult {
    /// Flat centroid data: `n_clusters * n_features` values.
    pub centroids: Vec<f64>,
    /// Cluster label for each data point.
    pub labels: Vec<usize>,
    /// Sum of squared distances to nearest centroid.
    pub inertia: f64,
    /// Number of iterations run.
    pub n_iter: usize,
    /// Dimensionality of the data.
    pub n_features: usize,
}

impl Summarizable for KMeansResult {
    fn summary(&self) -> String {
        let k = if self.n_features > 0 {
            self.centroids.len() / self.n_features
        } else {
            0
        };
        format!(
            "KMeans: k={}, inertia={:.4}, iterations={}",
            k, self.inertia, self.n_iter,
        )
    }
}

/// Run k-means clustering on the given data points.
///
/// Uses k-means++ initialization and Lloyd's algorithm.
pub fn kmeans(data: &[&[f64]], config: &KMeansConfig) -> Result<KMeansResult> {
    let n = data.len();
    let k = config.n_clusters;

    if n == 0 {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    if k == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_clusters must be > 0".into(),
        ));
    }
    if k > n {
        return Err(CyaneaError::InvalidInput(format!(
            "n_clusters ({}) > n_samples ({})",
            k, n
        )));
    }

    let dim = data[0].len();
    if dim == 0 {
        return Err(CyaneaError::InvalidInput("zero-dimensional data".into()));
    }
    for (i, row) in data.iter().enumerate() {
        if row.len() != dim {
            return Err(CyaneaError::InvalidInput(format!(
                "point {} has dimension {}, expected {}",
                i,
                row.len(),
                dim
            )));
        }
    }

    // k-means++ init
    let mut rng = Xorshift64(config.seed.max(1));
    let mut centroids = vec![0.0; k * dim];
    // First centroid: random point
    let first = rng.next_bounded(n as u64) as usize;
    centroids[..dim].copy_from_slice(data[first]);

    for c in 1..k {
        // Compute min distance from each point to existing centroids
        let mut dists = vec![f64::INFINITY; n];
        for i in 0..n {
            for prev in 0..c {
                let cent = &centroids[prev * dim..(prev + 1) * dim];
                let d = sq_euclidean(data[i], cent);
                if d < dists[i] {
                    dists[i] = d;
                }
            }
        }
        // Weighted random selection proportional to dist^2
        let total: f64 = dists.iter().sum();
        if total == 0.0 {
            // All points identical; just pick next point
            centroids[c * dim..(c + 1) * dim].copy_from_slice(data[c % n]);
        } else {
            let threshold = rng.next_f64() * total;
            let mut cumulative = 0.0;
            let mut chosen = n - 1;
            for (i, &d) in dists.iter().enumerate() {
                cumulative += d;
                if cumulative >= threshold {
                    chosen = i;
                    break;
                }
            }
            centroids[c * dim..(c + 1) * dim].copy_from_slice(data[chosen]);
        }
    }

    // Lloyd's iterations
    let mut labels = vec![0usize; n];
    let mut n_iter = 0;

    for _iter in 0..config.max_iter {
        n_iter += 1;

        // Assign each point to nearest centroid
        for i in 0..n {
            let mut best_dist = f64::INFINITY;
            let mut best_c = 0;
            for c in 0..k {
                let cent = &centroids[c * dim..(c + 1) * dim];
                let d = sq_euclidean(data[i], cent);
                if d < best_dist {
                    best_dist = d;
                    best_c = c;
                }
            }
            labels[i] = best_c;
        }

        // Update centroids
        let mut new_centroids = vec![0.0; k * dim];
        let mut counts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i];
            counts[c] += 1;
            for d in 0..dim {
                new_centroids[c * dim + d] += data[i][d];
            }
        }
        for c in 0..k {
            if counts[c] > 0 {
                let cnt = counts[c] as f64;
                for d in 0..dim {
                    new_centroids[c * dim + d] /= cnt;
                }
            } else {
                // Empty cluster: keep old centroid
                new_centroids[c * dim..(c + 1) * dim]
                    .copy_from_slice(&centroids[c * dim..(c + 1) * dim]);
            }
        }

        // Check convergence
        let mut max_shift = 0.0_f64;
        for c in 0..k {
            let shift = sq_euclidean(
                &centroids[c * dim..(c + 1) * dim],
                &new_centroids[c * dim..(c + 1) * dim],
            )
            .sqrt();
            if shift > max_shift {
                max_shift = shift;
            }
        }

        centroids = new_centroids;

        if max_shift < config.tolerance {
            break;
        }
    }

    // Compute inertia
    let mut inertia = 0.0;
    for i in 0..n {
        let c = labels[i];
        inertia += sq_euclidean(data[i], &centroids[c * dim..(c + 1) * dim]);
    }

    Ok(KMeansResult {
        centroids,
        labels,
        inertia,
        n_iter,
        n_features: dim,
    })
}

// ---------------------------------------------------------------------------
// DBSCAN
// ---------------------------------------------------------------------------

/// Configuration for DBSCAN clustering.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DbscanConfig {
    pub eps: f64,
    pub min_samples: usize,
    pub metric: DistanceMetric,
}

/// Result of DBSCAN clustering.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DbscanResult {
    /// Cluster label for each point. -1 indicates noise.
    pub labels: Vec<i32>,
    /// Number of clusters found (not counting noise).
    pub n_clusters: usize,
}

impl Summarizable for DbscanResult {
    fn summary(&self) -> String {
        let noise = self.labels.iter().filter(|&&l| l < 0).count();
        format!(
            "DBSCAN: {} clusters, {} noise points",
            self.n_clusters, noise,
        )
    }
}

/// Run DBSCAN clustering on the given data points.
pub fn dbscan(data: &[&[f64]], config: &DbscanConfig) -> Result<DbscanResult> {
    let n = data.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput("empty data".into()));
    }
    if config.eps <= 0.0 {
        return Err(CyaneaError::InvalidInput("eps must be > 0".into()));
    }
    if config.min_samples == 0 {
        return Err(CyaneaError::InvalidInput(
            "min_samples must be > 0".into(),
        ));
    }

    let dim = data[0].len();
    for (i, row) in data.iter().enumerate() {
        if row.len() != dim {
            return Err(CyaneaError::InvalidInput(format!(
                "point {} has dimension {}, expected {}",
                i,
                row.len(),
                dim
            )));
        }
    }

    // Precompute neighborhoods
    let mut neighborhoods: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut neighbors = Vec::new();
        for j in 0..n {
            if i == j {
                neighbors.push(j);
                continue;
            }
            let d = compute_distance(data[i], data[j], config.metric)?;
            if d <= config.eps {
                neighbors.push(j);
            }
        }
        neighborhoods.push(neighbors);
    }

    let mut labels = vec![-1i32; n];
    let mut cluster_id = 0i32;

    for i in 0..n {
        if labels[i] != -1 {
            continue; // already assigned
        }
        if neighborhoods[i].len() < config.min_samples {
            continue; // not a core point
        }

        // Expand cluster
        labels[i] = cluster_id;
        let mut queue: Vec<usize> = neighborhoods[i]
            .iter()
            .copied()
            .filter(|&j| j != i)
            .collect();
        let mut head = 0;

        while head < queue.len() {
            let j = queue[head];
            head += 1;

            if labels[j] == -1 {
                labels[j] = cluster_id;
            } else {
                continue; // already in a cluster
            }

            // If j is also a core point, add its neighbors
            if neighborhoods[j].len() >= config.min_samples {
                for &nb in &neighborhoods[j] {
                    if labels[nb] == -1 {
                        queue.push(nb);
                    }
                }
            }
        }

        cluster_id += 1;
    }

    Ok(DbscanResult {
        labels,
        n_clusters: cluster_id as usize,
    })
}

// ---------------------------------------------------------------------------
// Hierarchical (agglomerative)
// ---------------------------------------------------------------------------

/// Linkage criterion for hierarchical clustering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Linkage {
    Single,
    Complete,
    Average,
}

/// Configuration for hierarchical clustering.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HierarchicalConfig {
    pub n_clusters: usize,
    pub linkage: Linkage,
}

/// A single merge step in the dendrogram.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MergeStep {
    pub cluster_a: usize,
    pub cluster_b: usize,
    pub distance: f64,
    pub size: usize,
}

/// Result of hierarchical clustering.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HierarchicalResult {
    /// Cluster label for each point.
    pub labels: Vec<usize>,
    /// Ordered merge history.
    pub merge_history: Vec<MergeStep>,
}

impl Summarizable for HierarchicalResult {
    fn summary(&self) -> String {
        let k = if self.labels.is_empty() {
            0
        } else {
            *self.labels.iter().max().unwrap() + 1
        };
        format!(
            "Hierarchical: {} clusters, {} merges",
            k,
            self.merge_history.len(),
        )
    }
}

/// Run agglomerative hierarchical clustering on a precomputed distance matrix.
pub fn hierarchical(
    distances: &DistanceMatrix,
    config: &HierarchicalConfig,
) -> Result<HierarchicalResult> {
    let n = distances.n();
    if config.n_clusters == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_clusters must be > 0".into(),
        ));
    }
    if config.n_clusters > n {
        return Err(CyaneaError::InvalidInput(format!(
            "n_clusters ({}) > n_points ({})",
            config.n_clusters, n
        )));
    }

    // Each point starts in its own cluster.
    // cluster_members[i] = list of original point indices in cluster i.
    let mut cluster_members: Vec<Option<Vec<usize>>> =
        (0..n).map(|i| Some(vec![i])).collect();
    let mut active: Vec<usize> = (0..n).collect();
    let mut merge_history = Vec::with_capacity(n - config.n_clusters);

    // Working distance matrix (full n x n, but we only use active pairs)
    let mut dist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = distances.get(i, j);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    let target_clusters = config.n_clusters;

    while active.len() > target_clusters {
        // Find the two closest active clusters
        let mut best_dist = f64::INFINITY;
        let mut best_a = 0;
        let mut best_b = 0;
        for (ai, &a) in active.iter().enumerate() {
            for &b in &active[ai + 1..] {
                if dist[a][b] < best_dist {
                    best_dist = dist[a][b];
                    best_a = a;
                    best_b = b;
                }
            }
        }

        // Merge b into a
        let members_b = cluster_members[best_b].take().unwrap();
        let members_a = cluster_members[best_a].as_ref().unwrap();
        let new_size = members_a.len() + members_b.len();

        merge_history.push(MergeStep {
            cluster_a: best_a,
            cluster_b: best_b,
            distance: best_dist,
            size: new_size,
        });

        // Update distances from merged cluster to all other active clusters
        for &c in &active {
            if c == best_a || c == best_b {
                continue;
            }
            let d_ac = dist[best_a][c];
            let d_bc = dist[best_b][c];
            let new_d = match config.linkage {
                Linkage::Single => d_ac.min(d_bc),
                Linkage::Complete => d_ac.max(d_bc),
                Linkage::Average => {
                    let size_a = members_a.len() as f64;
                    let size_b = members_b.len() as f64;
                    (d_ac * size_a + d_bc * size_b) / (size_a + size_b)
                }
            };
            dist[best_a][c] = new_d;
            dist[c][best_a] = new_d;
        }

        // Merge member lists
        let mut merged = cluster_members[best_a].take().unwrap();
        merged.extend(members_b);
        cluster_members[best_a] = Some(merged);

        // Remove best_b from active
        active.retain(|&c| c != best_b);
    }

    // Assign labels
    let mut labels = vec![0usize; n];
    for (label, &cluster_idx) in active.iter().enumerate() {
        if let Some(members) = &cluster_members[cluster_idx] {
            for &pt in members {
                labels[pt] = label;
            }
        }
    }

    Ok(HierarchicalResult {
        labels,
        merge_history,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Squared Euclidean distance (no sqrt for speed).
fn sq_euclidean(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum()
}

/// Minimal xorshift64 PRNG — no external dependency needed.
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

    fn next_bounded(&mut self, bound: u64) -> u64 {
        self.next_u64() % bound
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_refs(data: &[Vec<f64>]) -> Vec<&[f64]> {
        data.iter().map(|v| v.as_slice()).collect()
    }

    // --- K-Means ---

    #[test]
    fn kmeans_two_clusters() {
        let data = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.1],
            vec![0.2, 0.0],
            vec![10.0, 10.0],
            vec![10.1, 10.1],
            vec![10.2, 10.0],
        ];
        let refs = make_refs(&data);
        let config = KMeansConfig {
            n_clusters: 2,
            ..Default::default()
        };
        let result = kmeans(&refs, &config).unwrap();
        assert_eq!(result.labels.len(), 6);
        // First 3 should be in the same cluster, last 3 in another
        assert_eq!(result.labels[0], result.labels[1]);
        assert_eq!(result.labels[0], result.labels[2]);
        assert_eq!(result.labels[3], result.labels[4]);
        assert_eq!(result.labels[3], result.labels[5]);
        assert_ne!(result.labels[0], result.labels[3]);
    }

    #[test]
    fn kmeans_single_cluster() {
        let data = vec![vec![1.0, 2.0], vec![1.1, 2.1], vec![0.9, 1.9]];
        let refs = make_refs(&data);
        let config = KMeansConfig {
            n_clusters: 1,
            ..Default::default()
        };
        let result = kmeans(&refs, &config).unwrap();
        assert!(result.labels.iter().all(|&l| l == 0));
    }

    #[test]
    fn kmeans_too_many_clusters() {
        let data = vec![vec![1.0], vec![2.0]];
        let refs = make_refs(&data);
        let config = KMeansConfig {
            n_clusters: 3,
            ..Default::default()
        };
        assert!(kmeans(&refs, &config).is_err());
    }

    #[test]
    fn kmeans_empty_data() {
        let refs: Vec<&[f64]> = vec![];
        let config = KMeansConfig::default();
        assert!(kmeans(&refs, &config).is_err());
    }

    #[test]
    fn kmeans_summary() {
        let data = vec![vec![0.0], vec![1.0], vec![10.0], vec![11.0]];
        let refs = make_refs(&data);
        let config = KMeansConfig {
            n_clusters: 2,
            ..Default::default()
        };
        let result = kmeans(&refs, &config).unwrap();
        let s = result.summary();
        assert!(s.contains("k=2"));
    }

    #[test]
    fn kmeans_convergence() {
        let data = vec![vec![0.0], vec![100.0]];
        let refs = make_refs(&data);
        let config = KMeansConfig {
            n_clusters: 2,
            ..Default::default()
        };
        let result = kmeans(&refs, &config).unwrap();
        // Should converge in very few iterations
        assert!(result.n_iter <= 5);
    }

    // --- DBSCAN ---

    #[test]
    fn dbscan_two_clusters() {
        let data = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.0],
            vec![0.0, 0.1],
            vec![10.0, 10.0],
            vec![10.1, 10.0],
            vec![10.0, 10.1],
        ];
        let refs = make_refs(&data);
        let config = DbscanConfig {
            eps: 1.0,
            min_samples: 2,
            metric: DistanceMetric::Euclidean,
        };
        let result = dbscan(&refs, &config).unwrap();
        assert_eq!(result.n_clusters, 2);
        assert_eq!(result.labels[0], result.labels[1]);
        assert_eq!(result.labels[3], result.labels[4]);
        assert_ne!(result.labels[0], result.labels[3]);
    }

    #[test]
    fn dbscan_all_noise() {
        let data = vec![vec![0.0], vec![10.0], vec![20.0]];
        let refs = make_refs(&data);
        let config = DbscanConfig {
            eps: 0.1,
            min_samples: 2,
            metric: DistanceMetric::Euclidean,
        };
        let result = dbscan(&refs, &config).unwrap();
        assert_eq!(result.n_clusters, 0);
        assert!(result.labels.iter().all(|&l| l == -1));
    }

    #[test]
    fn dbscan_one_cluster() {
        let data = vec![vec![0.0], vec![0.1], vec![0.2]];
        let refs = make_refs(&data);
        let config = DbscanConfig {
            eps: 1.0,
            min_samples: 2,
            metric: DistanceMetric::Euclidean,
        };
        let result = dbscan(&refs, &config).unwrap();
        assert_eq!(result.n_clusters, 1);
        assert!(result.labels.iter().all(|&l| l == 0));
    }

    #[test]
    fn dbscan_empty_data() {
        let refs: Vec<&[f64]> = vec![];
        let config = DbscanConfig {
            eps: 1.0,
            min_samples: 2,
            metric: DistanceMetric::Euclidean,
        };
        assert!(dbscan(&refs, &config).is_err());
    }

    #[test]
    fn dbscan_summary() {
        let data = vec![vec![0.0], vec![0.1], vec![10.0]];
        let refs = make_refs(&data);
        let config = DbscanConfig {
            eps: 0.5,
            min_samples: 2,
            metric: DistanceMetric::Euclidean,
        };
        let result = dbscan(&refs, &config).unwrap();
        let s = result.summary();
        assert!(s.contains("DBSCAN"));
    }

    // --- Hierarchical ---

    #[test]
    fn hierarchical_single_linkage() {
        let pts = vec![vec![0.0], vec![1.0], vec![2.0], vec![10.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 2,
            linkage: Linkage::Single,
        };
        let result = hierarchical(&dm, &config).unwrap();
        // 0, 1, 2 should cluster together; 10 is separate
        assert_eq!(result.labels[0], result.labels[1]);
        assert_eq!(result.labels[1], result.labels[2]);
        assert_ne!(result.labels[0], result.labels[3]);
    }

    #[test]
    fn hierarchical_complete_linkage() {
        let pts = vec![vec![0.0], vec![1.0], vec![10.0], vec![11.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 2,
            linkage: Linkage::Complete,
        };
        let result = hierarchical(&dm, &config).unwrap();
        assert_eq!(result.labels[0], result.labels[1]);
        assert_eq!(result.labels[2], result.labels[3]);
        assert_ne!(result.labels[0], result.labels[2]);
    }

    #[test]
    fn hierarchical_average_linkage() {
        let pts = vec![vec![0.0], vec![1.0], vec![10.0], vec![11.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 2,
            linkage: Linkage::Average,
        };
        let result = hierarchical(&dm, &config).unwrap();
        assert_eq!(result.labels[0], result.labels[1]);
        assert_eq!(result.labels[2], result.labels[3]);
        assert_ne!(result.labels[0], result.labels[2]);
    }

    #[test]
    fn hierarchical_single_cluster() {
        let pts = vec![vec![0.0], vec![1.0], vec![2.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 1,
            linkage: Linkage::Single,
        };
        let result = hierarchical(&dm, &config).unwrap();
        assert!(result.labels.iter().all(|&l| l == 0));
    }

    #[test]
    fn hierarchical_too_many_clusters() {
        let pts = vec![vec![0.0], vec![1.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 5,
            linkage: Linkage::Single,
        };
        assert!(hierarchical(&dm, &config).is_err());
    }

    #[test]
    fn hierarchical_summary() {
        let pts = vec![vec![0.0], vec![1.0], vec![10.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 2,
            linkage: Linkage::Single,
        };
        let result = hierarchical(&dm, &config).unwrap();
        let s = result.summary();
        assert!(s.contains("Hierarchical"));
        assert!(s.contains("2 clusters"));
    }

    #[test]
    fn hierarchical_merge_history() {
        let pts = vec![vec![0.0], vec![1.0], vec![10.0]];
        let refs = make_refs(&pts);
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        let config = HierarchicalConfig {
            n_clusters: 1,
            linkage: Linkage::Single,
        };
        let result = hierarchical(&dm, &config).unwrap();
        assert_eq!(result.merge_history.len(), 2); // n-1 merges for n=3 -> 1 cluster
        // First merge should be distance 1.0 (between 0.0 and 1.0)
        assert!((result.merge_history[0].distance - 1.0).abs() < 1e-12);
    }
}
