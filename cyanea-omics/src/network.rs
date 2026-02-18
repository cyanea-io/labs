//! Network biology — weighted graphs, centrality metrics, and community detection.
//!
//! Build graphs from correlation matrices, compute degree/betweenness/closeness
//! centrality, and detect communities with the Louvain algorithm.

use std::collections::VecDeque;

use cyanea_core::{CyaneaError, Result};

/// A weighted graph (directed or undirected).
#[derive(Debug, Clone)]
pub struct Graph {
    n_nodes: usize,
    edges: Vec<(usize, usize, f64)>,
    adjacency: Vec<Vec<(usize, f64)>>,
    directed: bool,
}

/// Centrality scores for all nodes in a graph.
#[derive(Debug, Clone)]
pub struct CentralityScores {
    /// Degree centrality for each node.
    pub degree: Vec<f64>,
    /// Betweenness centrality for each node.
    pub betweenness: Vec<f64>,
    /// Closeness centrality for each node.
    pub closeness: Vec<f64>,
}

/// Community detection result.
#[derive(Debug, Clone)]
pub struct Community {
    /// Community assignment for each node.
    pub assignments: Vec<usize>,
    /// Modularity Q of the partition.
    pub modularity: f64,
}

impl Graph {
    /// Create an empty graph with `n_nodes` nodes.
    pub fn new(n_nodes: usize, directed: bool) -> Self {
        Self {
            n_nodes,
            edges: Vec::new(),
            adjacency: vec![Vec::new(); n_nodes],
            directed,
        }
    }

    /// Add a weighted edge.
    ///
    /// # Errors
    ///
    /// Returns an error if either node index is out of range.
    pub fn add_edge(&mut self, from: usize, to: usize, weight: f64) -> Result<()> {
        if from >= self.n_nodes || to >= self.n_nodes {
            return Err(CyaneaError::InvalidInput(format!(
                "node index out of range: from={}, to={}, n_nodes={}",
                from, to, self.n_nodes
            )));
        }
        self.edges.push((from, to, weight));
        self.adjacency[from].push((to, weight));
        if !self.directed {
            self.adjacency[to].push((from, weight));
        }
        Ok(())
    }

    /// Build an undirected graph from a correlation matrix with a threshold.
    ///
    /// Adds an edge between nodes i and j if |correlation[i][j]| ≥ threshold.
    /// Edge weight is the absolute correlation.
    pub fn from_correlation_matrix(matrix: &[Vec<f64>], threshold: f64) -> Self {
        let n = matrix.len();
        let mut graph = Self::new(n, false);
        for i in 0..n {
            for j in (i + 1)..n {
                let corr = matrix[i][j].abs();
                if corr >= threshold {
                    let _ = graph.add_edge(i, j, corr);
                }
            }
        }
        graph
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.n_nodes
    }

    /// Number of edges.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Neighbors of a node with edge weights.
    pub fn neighbors(&self, node: usize) -> &[(usize, f64)] {
        if node < self.n_nodes {
            &self.adjacency[node]
        } else {
            &[]
        }
    }

    /// Degree centrality: degree(v) / (n - 1).
    pub fn degree_centrality(&self) -> Vec<f64> {
        let n = self.n_nodes;
        if n <= 1 {
            return vec![0.0; n];
        }
        let denom = (n - 1) as f64;
        (0..n)
            .map(|v| self.adjacency[v].len() as f64 / denom)
            .collect()
    }

    /// Betweenness centrality using Brandes' algorithm.
    pub fn betweenness_centrality(&self) -> Vec<f64> {
        let n = self.n_nodes;
        let mut cb = vec![0.0f64; n];

        for s in 0..n {
            // BFS from s (unweighted shortest paths).
            let mut stack = Vec::new();
            let mut pred: Vec<Vec<usize>> = vec![Vec::new(); n];
            let mut sigma = vec![0.0f64; n];
            sigma[s] = 1.0;
            let mut dist = vec![-1i64; n];
            dist[s] = 0;

            let mut queue = VecDeque::new();
            queue.push_back(s);

            while let Some(v) = queue.pop_front() {
                stack.push(v);
                for &(w, _) in &self.adjacency[v] {
                    // First visit?
                    if dist[w] < 0 {
                        dist[w] = dist[v] + 1;
                        queue.push_back(w);
                    }
                    // Shortest path through v?
                    if dist[w] == dist[v] + 1 {
                        sigma[w] += sigma[v];
                        pred[w].push(v);
                    }
                }
            }

            // Accumulate dependencies.
            let mut delta = vec![0.0f64; n];
            while let Some(w) = stack.pop() {
                for &v in &pred[w] {
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
                if w != s {
                    cb[w] += delta[w];
                }
            }
        }

        // For undirected graphs, each pair is counted twice.
        if !self.directed {
            for v in &mut cb {
                *v /= 2.0;
            }
        }

        // Normalize by (n-1)(n-2) for comparability.
        let norm = if n > 2 {
            ((n - 1) * (n - 2)) as f64
        } else {
            1.0
        };
        for v in &mut cb {
            *v /= norm;
        }

        cb
    }

    /// Closeness centrality: (n-1) / Σ d(v, u).
    pub fn closeness_centrality(&self) -> Vec<f64> {
        let n = self.n_nodes;
        if n <= 1 {
            return vec![0.0; n];
        }

        (0..n)
            .map(|v| {
                let distances = self.bfs_distances(v);
                let sum_dist: usize = distances.iter().filter(|&&d| d > 0).sum();
                let reachable = distances.iter().filter(|&&d| d > 0).count();
                if sum_dist > 0 && reachable > 0 {
                    reachable as f64 / sum_dist as f64
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Compute all three centrality metrics.
    pub fn centrality(&self) -> CentralityScores {
        CentralityScores {
            degree: self.degree_centrality(),
            betweenness: self.betweenness_centrality(),
            closeness: self.closeness_centrality(),
        }
    }

    /// Louvain community detection.
    ///
    /// Iteratively moves nodes between communities to maximize modularity.
    pub fn louvain(&self) -> Community {
        if self.n_nodes == 0 {
            return Community {
                assignments: Vec::new(),
                modularity: 0.0,
            };
        }

        let n = self.n_nodes;
        let mut assignments: Vec<usize> = (0..n).collect();

        // Total edge weight (2m for undirected).
        let m2: f64 = if self.directed {
            self.edges.iter().map(|(_, _, w)| w).sum()
        } else {
            self.edges.iter().map(|(_, _, w)| 2.0 * w).sum()
        };

        if m2 == 0.0 {
            return Community {
                assignments,
                modularity: 0.0,
            };
        }

        // Phase 1: local moves.
        let mut improved = true;
        while improved {
            improved = false;

            for i in 0..n {
                let current_community = assignments[i];

                // Compute k_i (weighted degree of node i).
                let k_i: f64 = self.adjacency[i].iter().map(|(_, w)| w).sum();

                // Sum of weights to each neighbor community.
                let mut community_weights: Vec<(usize, f64)> = Vec::new();
                for &(j, w) in &self.adjacency[i] {
                    let cj = assignments[j];
                    if let Some(entry) = community_weights.iter_mut().find(|(c, _)| *c == cj) {
                        entry.1 += w;
                    } else {
                        community_weights.push((cj, w));
                    }
                }

                // Σ_tot for current community.
                let sigma_tot_current = self.community_total_weight(&assignments, current_community);
                let k_i_in_current = community_weights
                    .iter()
                    .find(|(c, _)| *c == current_community)
                    .map_or(0.0, |(_, w)| *w);

                let mut best_community = current_community;
                let mut best_delta_q = 0.0;

                for &(cj, k_i_in) in &community_weights {
                    if cj == current_community {
                        continue;
                    }
                    let sigma_tot = self.community_total_weight(&assignments, cj);

                    // ΔQ for moving i from current to cj.
                    let delta_q = (k_i_in - k_i_in_current) / m2
                        - k_i * (sigma_tot - sigma_tot_current + k_i) / (m2 * m2) * 2.0;

                    if delta_q > best_delta_q {
                        best_delta_q = delta_q;
                        best_community = cj;
                    }
                }

                if best_community != current_community {
                    assignments[i] = best_community;
                    improved = true;
                }
            }
        }

        // Renumber communities to be contiguous.
        let mut community_map: Vec<usize> = Vec::new();
        for a in &mut assignments {
            let pos = community_map.iter().position(|&c| c == *a);
            *a = match pos {
                Some(idx) => idx,
                None => {
                    community_map.push(*a);
                    community_map.len() - 1
                }
            };
        }

        let modularity = self.modularity(&assignments);
        Community {
            assignments,
            modularity,
        }
    }

    /// Compute modularity Q for a given community assignment.
    ///
    /// Uses the community-level formula:
    /// Q = Σ_c [L_c / m - (d_c / 2m)²]
    /// where L_c = sum of edge weights within community c, m = total edge weight,
    /// and d_c = sum of node degrees in community c.
    pub fn modularity(&self, assignments: &[usize]) -> f64 {
        let m: f64 = self.edges.iter().map(|(_, _, w)| w).sum();
        if m == 0.0 {
            return 0.0;
        }

        // Find unique communities.
        let mut communities: Vec<usize> = assignments.to_vec();
        communities.sort_unstable();
        communities.dedup();

        let mut q = 0.0;
        for &c in &communities {
            // L_c: sum of edge weights where both endpoints are in community c.
            let l_c: f64 = self
                .edges
                .iter()
                .filter(|&&(i, j, _)| {
                    assignments.get(i) == Some(&c) && assignments.get(j) == Some(&c)
                })
                .map(|(_, _, w)| w)
                .sum();

            // d_c: sum of weighted degrees of nodes in community c.
            let d_c: f64 = (0..self.n_nodes)
                .filter(|&v| assignments.get(v) == Some(&c))
                .map(|v| -> f64 { self.adjacency[v].iter().map(|(_, w)| w).sum() })
                .sum();

            q += l_c / m - (d_c / (2.0 * m)).powi(2);
        }
        q
    }

    fn bfs_distances(&self, start: usize) -> Vec<usize> {
        let n = self.n_nodes;
        let mut dist = vec![usize::MAX; n];
        dist[start] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(v) = queue.pop_front() {
            for &(w, _) in &self.adjacency[v] {
                if dist[w] == usize::MAX {
                    dist[w] = dist[v] + 1;
                    queue.push_back(w);
                }
            }
        }
        dist
    }

    fn community_total_weight(&self, assignments: &[usize], community: usize) -> f64 {
        let mut total = 0.0;
        for (v, neighbors) in self.adjacency.iter().enumerate() {
            if assignments[v] == community {
                for &(_, w) in neighbors {
                    total += w;
                }
            }
        }
        total
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn graph_creation() {
        let mut g = Graph::new(5, false);
        g.add_edge(0, 1, 1.0).unwrap();
        g.add_edge(1, 2, 1.0).unwrap();
        assert_eq!(g.n_nodes(), 5);
        assert_eq!(g.n_edges(), 2);
    }

    #[test]
    fn degree_centrality_star() {
        // Star graph: center (0) connected to 1,2,3,4.
        let mut g = Graph::new(5, false);
        for i in 1..5 {
            g.add_edge(0, i, 1.0).unwrap();
        }
        let dc = g.degree_centrality();
        // Center has degree 4, centrality = 4/4 = 1.0.
        assert!((dc[0] - 1.0).abs() < 1e-10);
        // Leaves have degree 1, centrality = 1/4 = 0.25.
        for i in 1..5 {
            assert!((dc[i] - 0.25).abs() < 1e-10);
        }
    }

    #[test]
    fn betweenness_centrality_line() {
        // Line: 0-1-2-3-4. Node 2 should have highest betweenness.
        let mut g = Graph::new(5, false);
        for i in 0..4 {
            g.add_edge(i, i + 1, 1.0).unwrap();
        }
        let bc = g.betweenness_centrality();
        // Node 2 (center of line) has highest betweenness.
        let max_node = bc
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(max_node, 2);
    }

    #[test]
    fn closeness_centrality_complete() {
        // Complete graph: all nodes equidistant → equal closeness.
        let mut g = Graph::new(4, false);
        for i in 0..4 {
            for j in (i + 1)..4 {
                g.add_edge(i, j, 1.0).unwrap();
            }
        }
        let cc = g.closeness_centrality();
        for i in 0..4 {
            assert!((cc[i] - cc[0]).abs() < 1e-10);
        }
        // Closeness should be 1.0: (n-1) / sum(distances) = 3/3 = 1.0.
        assert!((cc[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn louvain_two_cliques() {
        // Two cliques of 3 nodes each, connected by a single weak edge.
        let mut g = Graph::new(6, false);
        // Clique 1: 0-1-2.
        g.add_edge(0, 1, 1.0).unwrap();
        g.add_edge(1, 2, 1.0).unwrap();
        g.add_edge(0, 2, 1.0).unwrap();
        // Clique 2: 3-4-5.
        g.add_edge(3, 4, 1.0).unwrap();
        g.add_edge(4, 5, 1.0).unwrap();
        g.add_edge(3, 5, 1.0).unwrap();
        // Weak bridge.
        g.add_edge(2, 3, 0.1).unwrap();

        let community = g.louvain();
        // Nodes in the same clique should be in the same community.
        assert_eq!(community.assignments[0], community.assignments[1]);
        assert_eq!(community.assignments[1], community.assignments[2]);
        assert_eq!(community.assignments[3], community.assignments[4]);
        assert_eq!(community.assignments[4], community.assignments[5]);
        // The two cliques should be in different communities.
        assert_ne!(community.assignments[0], community.assignments[3]);
        assert!(community.modularity > 0.0);
    }

    #[test]
    fn modularity_single_community() {
        // All nodes in one community → Q should be 0.
        let mut g = Graph::new(3, false);
        g.add_edge(0, 1, 1.0).unwrap();
        g.add_edge(1, 2, 1.0).unwrap();
        let assignments = vec![0, 0, 0];
        let q = g.modularity(&assignments);
        assert!(q.abs() < 1e-10, "Q should be ~0 for single community, got {}", q);
    }

    #[test]
    fn from_correlation_matrix_threshold() {
        let matrix = vec![
            vec![1.0, 0.9, 0.2],
            vec![0.9, 1.0, 0.3],
            vec![0.2, 0.3, 1.0],
        ];
        let g = Graph::from_correlation_matrix(&matrix, 0.5);
        // Only the 0.9 edge should pass the 0.5 threshold.
        assert_eq!(g.n_edges(), 1);
        assert_eq!(g.n_nodes(), 3);
    }
}
