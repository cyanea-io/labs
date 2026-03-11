//! Community detection algorithms for biological networks.
//!
//! Implements Louvain modularity optimization and label propagation for
//! identifying functional modules in PPI networks and gene co-expression networks.

use crate::graph::{Graph, GraphType};
use cyanea_core::Result;
use std::collections::HashMap;

/// Result of community detection.
#[derive(Debug, Clone)]
pub struct CommunityResult {
    /// Node ID → community label.
    pub assignments: HashMap<String, usize>,
    /// Modularity score (Q) of the partition.
    pub modularity: f64,
    /// Number of communities found.
    pub num_communities: usize,
}

impl CommunityResult {
    /// Get all nodes in a given community.
    pub fn members(&self, community: usize) -> Vec<&str> {
        self.assignments
            .iter()
            .filter(|(_, &c)| c == community)
            .map(|(id, _)| id.as_str())
            .collect()
    }

    /// Get community sizes as (community_id, size) pairs, sorted by size descending.
    pub fn community_sizes(&self) -> Vec<(usize, usize)> {
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for &c in self.assignments.values() {
            *counts.entry(c).or_insert(0) += 1;
        }
        let mut sizes: Vec<(usize, usize)> = counts.into_iter().collect();
        sizes.sort_by(|a, b| b.1.cmp(&a.1));
        sizes
    }
}

/// Compute Newman-Girvan modularity for a given partition.
///
/// Q = (1/2m) * Σ_ij [A_ij - k_i*k_j/(2m)] * δ(c_i, c_j)
///
/// where m = total edge weight, k_i = weighted degree, c_i = community of node i.
pub fn modularity(graph: &Graph, assignments: &HashMap<String, usize>) -> f64 {
    let (ids, matrix) = graph.adjacency_matrix();
    let n = ids.len();
    if n == 0 {
        return 0.0;
    }

    let _id_to_idx: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

    // Total edge weight
    let m: f64 = match graph.graph_type {
        GraphType::Undirected => graph.edges.iter().map(|e| e.weight).sum::<f64>(),
        GraphType::Directed => graph.edges.iter().map(|e| e.weight).sum(),
    };

    if m == 0.0 {
        return 0.0;
    }

    let two_m = 2.0 * m;

    // Weighted degree of each node
    let mut k = vec![0.0; n];
    for i in 0..n {
        k[i] = matrix[i].iter().sum();
    }

    let mut q = 0.0;
    for i in 0..n {
        let ci = assignments.get(ids[i].as_str()).copied().unwrap_or(0);
        for j in 0..n {
            let cj = assignments.get(ids[j].as_str()).copied().unwrap_or(0);
            if ci == cj {
                q += matrix[i][j] - (k[i] * k[j]) / two_m;
            }
        }
    }

    q / two_m
}

/// Louvain community detection algorithm.
///
/// Greedy modularity optimization: iteratively move nodes to the community
/// that gives the largest modularity gain, then aggregate communities.
///
/// # Arguments
///
/// * `graph` - The network to partition
/// * `resolution` - Resolution parameter (1.0 = standard modularity, >1.0 = smaller communities)
pub fn louvain(graph: &Graph, resolution: f64) -> Result<CommunityResult> {
    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let n = ids.len();

    if n == 0 {
        return Ok(CommunityResult {
            assignments: HashMap::new(),
            modularity: 0.0,
            num_communities: 0,
        });
    }

    // Initialize: each node in its own community
    let mut community: HashMap<String, usize> = ids
        .iter()
        .enumerate()
        .map(|(i, id)| (id.clone(), i))
        .collect();

    let (_, matrix) = graph.adjacency_matrix();

    // Total edge weight
    let m: f64 = match graph.graph_type {
        GraphType::Undirected => graph.edges.iter().map(|e| e.weight).sum::<f64>(),
        GraphType::Directed => graph.edges.iter().map(|e| e.weight).sum(),
    };

    if m == 0.0 {
        let num_c = community.values().copied().max().map(|m| m + 1).unwrap_or(0);
        return Ok(CommunityResult {
            modularity: 0.0,
            num_communities: num_c,
            assignments: community,
        });
    }

    let two_m = 2.0 * m;

    // Sorted IDs for consistent indexing
    let mut sorted_ids: Vec<String> = ids.clone();
    sorted_ids.sort();
    let id_to_idx: HashMap<&str, usize> = sorted_ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

    // Weighted degree
    let k: Vec<f64> = (0..n).map(|i| matrix[i].iter().sum()).collect();

    // Community total weight
    let mut sigma_tot: HashMap<usize, f64> = HashMap::new();
    for (id, &c) in &community {
        let idx = id_to_idx[id.as_str()];
        *sigma_tot.entry(c).or_insert(0.0) += k[idx];
    }

    let max_passes = 20;
    for _ in 0..max_passes {
        let mut improved = false;

        for id in &sorted_ids {
            let i = id_to_idx[id.as_str()];
            let current_c = community[id];

            // Compute weights to each neighboring community
            let mut comm_weights: HashMap<usize, f64> = HashMap::new();
            for j in 0..n {
                if matrix[i][j] > 0.0 {
                    let cj = community[&sorted_ids[j]];
                    *comm_weights.entry(cj).or_insert(0.0) += matrix[i][j];
                }
            }

            // Remove node from current community
            *sigma_tot.get_mut(&current_c).unwrap() -= k[i];

            // Find best community
            let mut best_c = current_c;
            let mut best_gain = 0.0;

            for (&c, &w_ic) in &comm_weights {
                let st = sigma_tot.get(&c).copied().unwrap_or(0.0);
                let gain = w_ic - resolution * k[i] * st / two_m;
                if gain > best_gain {
                    best_gain = gain;
                    best_c = c;
                }
            }

            // Also check staying in current community
            let w_current = comm_weights.get(&current_c).copied().unwrap_or(0.0);
            let st_current = sigma_tot.get(&current_c).copied().unwrap_or(0.0);
            let stay_gain = w_current - resolution * k[i] * st_current / two_m;
            if stay_gain >= best_gain {
                best_c = current_c;
            }

            // Move node to best community
            *sigma_tot.entry(best_c).or_insert(0.0) += k[i];
            if best_c != current_c {
                community.insert(id.clone(), best_c);
                improved = true;
            }
        }

        if !improved {
            break;
        }
    }

    // Renumber communities to be contiguous
    let mut label_map: HashMap<usize, usize> = HashMap::new();
    let mut next_label = 0;
    for c in community.values() {
        if !label_map.contains_key(c) {
            label_map.insert(*c, next_label);
            next_label += 1;
        }
    }

    let assignments: HashMap<String, usize> = community
        .into_iter()
        .map(|(id, c)| (id, label_map[&c]))
        .collect();

    let mod_score = modularity(graph, &assignments);

    Ok(CommunityResult {
        num_communities: next_label,
        modularity: mod_score,
        assignments,
    })
}

/// Label propagation community detection.
///
/// Each node adopts the most frequent label among its neighbors.
/// Simple, fast, and doesn't require a resolution parameter.
/// Non-deterministic due to iteration order.
pub fn label_propagation(graph: &Graph, max_iter: usize) -> Result<CommunityResult> {
    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let n = ids.len();

    if n == 0 {
        return Ok(CommunityResult {
            assignments: HashMap::new(),
            modularity: 0.0,
            num_communities: 0,
        });
    }

    // Initialize each node with unique label
    let mut labels: HashMap<String, usize> = ids
        .iter()
        .enumerate()
        .map(|(i, id)| (id.clone(), i))
        .collect();

    for _ in 0..max_iter {
        let mut changed = false;

        for id in &ids {
            let neighbors = graph.neighbors(id).unwrap_or_default();
            if neighbors.is_empty() {
                continue;
            }

            // Count weighted labels among neighbors
            let mut label_weights: HashMap<usize, f64> = HashMap::new();
            for (neighbor, weight) in &neighbors {
                let label = labels[*neighbor];
                *label_weights.entry(label).or_insert(0.0) += weight;
            }

            // Find most common label
            let best_label = label_weights
                .into_iter()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(l, _)| l)
                .unwrap();

            if best_label != labels[id] {
                labels.insert(id.clone(), best_label);
                changed = true;
            }
        }

        if !changed {
            break;
        }
    }

    // Renumber labels
    let mut label_map: HashMap<usize, usize> = HashMap::new();
    let mut next = 0;
    for l in labels.values() {
        if !label_map.contains_key(l) {
            label_map.insert(*l, next);
            next += 1;
        }
    }

    let assignments: HashMap<String, usize> = labels
        .into_iter()
        .map(|(id, l)| (id, label_map[&l]))
        .collect();

    let mod_score = modularity(graph, &assignments);

    Ok(CommunityResult {
        num_communities: next,
        modularity: mod_score,
        assignments,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_cliques() -> Graph {
        // Two triangles connected by a bridge
        let mut g = Graph::new(GraphType::Undirected);
        for id in &["A", "B", "C", "D", "E", "F"] {
            g.add_node(id, id).unwrap();
        }
        // Clique 1: A-B-C
        g.add_edge("A", "B", 1.0).unwrap();
        g.add_edge("B", "C", 1.0).unwrap();
        g.add_edge("A", "C", 1.0).unwrap();
        // Clique 2: D-E-F
        g.add_edge("D", "E", 1.0).unwrap();
        g.add_edge("E", "F", 1.0).unwrap();
        g.add_edge("D", "F", 1.0).unwrap();
        // Bridge
        g.add_edge("C", "D", 1.0).unwrap();
        g
    }

    #[test]
    fn test_modularity_trivial() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        g.add_edge("A", "B", 1.0).unwrap();
        let assign: HashMap<String, usize> =
            vec![("A".into(), 0), ("B".into(), 0)].into_iter().collect();
        let q = modularity(&g, &assign);
        // Single community: Q = 0 for complete graph
        assert!(q.abs() < 1e-10);
    }

    #[test]
    fn test_modularity_two_cliques() {
        let g = two_cliques();
        // Optimal partition: {A,B,C} and {D,E,F}
        let assign: HashMap<String, usize> = vec![
            ("A".into(), 0),
            ("B".into(), 0),
            ("C".into(), 0),
            ("D".into(), 1),
            ("E".into(), 1),
            ("F".into(), 1),
        ]
        .into_iter()
        .collect();
        let q = modularity(&g, &assign);
        assert!(q > 0.3); // Should be significantly positive
    }

    #[test]
    fn test_louvain_two_cliques() {
        let g = two_cliques();
        let result = louvain(&g, 1.0).unwrap();
        assert_eq!(result.num_communities, 2);
        assert!(result.modularity > 0.3);
        // A, B, C should be in one community; D, E, F in another
        assert_eq!(result.assignments["A"], result.assignments["B"]);
        assert_eq!(result.assignments["B"], result.assignments["C"]);
        assert_eq!(result.assignments["D"], result.assignments["E"]);
        assert_eq!(result.assignments["E"], result.assignments["F"]);
        assert_ne!(result.assignments["A"], result.assignments["D"]);
    }

    #[test]
    fn test_louvain_single_node() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        let result = louvain(&g, 1.0).unwrap();
        assert_eq!(result.num_communities, 1);
    }

    #[test]
    fn test_louvain_empty() {
        let g = Graph::new(GraphType::Undirected);
        let result = louvain(&g, 1.0).unwrap();
        assert_eq!(result.num_communities, 0);
    }

    #[test]
    fn test_label_propagation_two_cliques() {
        let g = two_cliques();
        let result = label_propagation(&g, 100).unwrap();
        assert!(result.num_communities >= 1);
        // Within each clique, nodes should share a label
        assert_eq!(result.assignments["A"], result.assignments["B"]);
        assert_eq!(result.assignments["D"], result.assignments["E"]);
    }

    #[test]
    fn test_community_sizes() {
        let g = two_cliques();
        let result = louvain(&g, 1.0).unwrap();
        let sizes = result.community_sizes();
        assert_eq!(sizes.len(), 2);
        assert_eq!(sizes[0].1, 3);
        assert_eq!(sizes[1].1, 3);
    }

    #[test]
    fn test_community_members() {
        let g = two_cliques();
        let result = louvain(&g, 1.0).unwrap();
        let c0 = result.members(result.assignments["A"]);
        assert_eq!(c0.len(), 3);
    }
}
