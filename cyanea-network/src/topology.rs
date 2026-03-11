//! Network topology metrics: centrality, shortest paths, clustering.
//!
//! Computes standard graph-theoretic measures useful for analyzing biological
//! networks including protein-protein interaction networks and gene regulatory
//! networks.

use crate::graph::Graph;
use cyanea_core::Result;
use std::collections::{HashMap, HashSet, VecDeque};

/// Centrality scores for all nodes in a graph.
#[derive(Debug, Clone)]
pub struct CentralityResult {
    /// Node ID → centrality score.
    pub scores: HashMap<String, f64>,
}

impl CentralityResult {
    /// Get the top-k nodes by centrality score.
    pub fn top_k(&self, k: usize) -> Vec<(&str, f64)> {
        let mut sorted: Vec<(&str, f64)> = self
            .scores
            .iter()
            .map(|(id, &score)| (id.as_str(), score))
            .collect();
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        sorted.truncate(k);
        sorted
    }
}

/// Degree centrality: C_d(v) = deg(v) / (n - 1).
pub fn degree_centrality(graph: &Graph) -> Result<CentralityResult> {
    let n = graph.node_count();
    if n <= 1 {
        let scores = graph
            .node_ids()
            .into_iter()
            .map(|id| (id.to_string(), 0.0))
            .collect();
        return Ok(CentralityResult { scores });
    }

    let denom = (n - 1) as f64;
    let scores = graph
        .node_ids()
        .into_iter()
        .map(|id| {
            let deg = graph.degree(id).unwrap_or(0) as f64;
            (id.to_string(), deg / denom)
        })
        .collect();

    Ok(CentralityResult { scores })
}

/// Closeness centrality: C_c(v) = (n-1) / sum(d(v,u) for all u).
///
/// Uses BFS shortest paths (unweighted). Nodes in disconnected components
/// get closeness 0.0.
pub fn closeness_centrality(graph: &Graph) -> Result<CentralityResult> {
    let n = graph.node_count();
    if n <= 1 {
        let scores = graph
            .node_ids()
            .into_iter()
            .map(|id| (id.to_string(), 0.0))
            .collect();
        return Ok(CentralityResult { scores });
    }

    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let mut scores = HashMap::new();

    for id in &ids {
        let dists = bfs_distances(graph, id);
        let reachable = dists.values().filter(|&&d| d > 0).count();
        if reachable == 0 {
            scores.insert(id.clone(), 0.0);
        } else {
            let sum_dist: usize = dists.values().sum();
            // Wasserman-Faust normalization for disconnected graphs
            let closeness = if sum_dist > 0 {
                (reachable as f64) / (sum_dist as f64)
                    * (reachable as f64 / (n - 1) as f64)
            } else {
                0.0
            };
            scores.insert(id.clone(), closeness);
        }
    }

    Ok(CentralityResult { scores })
}

/// Betweenness centrality: C_b(v) = sum(σ_st(v)/σ_st for all s≠v≠t).
///
/// Uses Brandes' algorithm (BFS-based, unweighted). Normalized by
/// 1/((n-1)(n-2)) for directed graphs and 2/((n-1)(n-2)) for undirected.
pub fn betweenness_centrality(graph: &Graph) -> Result<CentralityResult> {
    let n = graph.node_count();
    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let mut cb: HashMap<String, f64> = ids.iter().map(|id| (id.clone(), 0.0)).collect();

    if n <= 2 {
        return Ok(CentralityResult { scores: cb });
    }

    // Brandes' algorithm
    for s in &ids {
        let mut stack: Vec<String> = Vec::new();
        let mut predecessors: HashMap<String, Vec<String>> = HashMap::new();
        let mut sigma: HashMap<String, f64> = ids.iter().map(|id| (id.clone(), 0.0)).collect();
        *sigma.get_mut(s).unwrap() = 1.0;
        let mut dist: HashMap<String, i64> = ids.iter().map(|id| (id.clone(), -1)).collect();
        *dist.get_mut(s).unwrap() = 0;

        let mut queue = VecDeque::new();
        queue.push_back(s.clone());

        while let Some(v) = queue.pop_front() {
            stack.push(v.clone());
            let neighbors = graph.neighbors(&v).unwrap_or_default();
            for (w, _) in neighbors {
                let w = w.to_string();
                if dist[&w] < 0 {
                    queue.push_back(w.clone());
                    *dist.get_mut(&w).unwrap() = dist[&v] + 1;
                }
                if dist[&w] == dist[&v] + 1 {
                    *sigma.get_mut(&w).unwrap() += sigma[&v];
                    predecessors.entry(w).or_default().push(v.clone());
                }
            }
        }

        let mut delta: HashMap<String, f64> = ids.iter().map(|id| (id.clone(), 0.0)).collect();
        while let Some(w) = stack.pop() {
            if let Some(preds) = predecessors.get(&w) {
                for v in preds {
                    let d = (sigma[v] / sigma[&w]) * (1.0 + delta[&w]);
                    *delta.get_mut(v).unwrap() += d;
                }
            }
            if w != *s {
                *cb.get_mut(&w).unwrap() += delta[&w];
            }
        }
    }

    // Normalize
    let norm = match graph.graph_type {
        crate::graph::GraphType::Directed => {
            ((n - 1) * (n - 2)) as f64
        }
        crate::graph::GraphType::Undirected => {
            ((n - 1) * (n - 2)) as f64 / 2.0
        }
    };

    if norm > 0.0 {
        for score in cb.values_mut() {
            *score /= norm;
        }
    }

    Ok(CentralityResult { scores: cb })
}

/// PageRank centrality.
///
/// Iterative power method with damping factor (default 0.85).
/// Converges when max change < tolerance (default 1e-6).
pub fn pagerank(
    graph: &Graph,
    damping: f64,
    max_iter: usize,
    tolerance: f64,
) -> Result<CentralityResult> {
    let n = graph.node_count();
    if n == 0 {
        return Ok(CentralityResult {
            scores: HashMap::new(),
        });
    }

    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let mut rank: HashMap<String, f64> = ids.iter().map(|id| (id.clone(), 1.0 / n as f64)).collect();

    for _ in 0..max_iter {
        let mut new_rank: HashMap<String, f64> =
            ids.iter().map(|id| (id.clone(), (1.0 - damping) / n as f64)).collect();

        // Accumulate dangling node mass
        let mut dangling_sum = 0.0;
        for id in &ids {
            let deg = graph.degree(id).unwrap_or(0);
            if deg == 0 {
                dangling_sum += rank[id];
            }
        }

        // Distribute dangling mass evenly
        let dangling_contrib = damping * dangling_sum / n as f64;
        for val in new_rank.values_mut() {
            *val += dangling_contrib;
        }

        // Add edge contributions
        for id in &ids {
            let deg = graph.degree(id).unwrap_or(0);
            if deg > 0 {
                let neighbors = graph.neighbors(id).unwrap_or_default();
                let contrib = damping * rank[id] / deg as f64;
                for (neighbor, _) in &neighbors {
                    if let Some(r) = new_rank.get_mut(*neighbor) {
                        *r += contrib;
                    }
                }
            }
        }

        // Check convergence
        let max_diff = ids
            .iter()
            .map(|id| (new_rank[id] - rank[id]).abs())
            .fold(0.0f64, f64::max);

        rank = new_rank;

        if max_diff < tolerance {
            break;
        }
    }

    Ok(CentralityResult { scores: rank })
}

/// Local clustering coefficient for a node.
///
/// The fraction of pairs of a node's neighbors that are connected.
/// Returns 0.0 for nodes with degree < 2.
pub fn clustering_coefficient(graph: &Graph, node_id: &str) -> Result<f64> {
    let neighbors = graph.neighbors(node_id)?;
    let k = neighbors.len();
    if k < 2 {
        return Ok(0.0);
    }

    let neighbor_set: HashSet<&str> = neighbors.iter().map(|(id, _)| *id).collect();
    let mut triangles = 0;

    for (n1, _) in &neighbors {
        let n1_neighbors = graph.neighbors(n1).unwrap_or_default();
        for (n2, _) in &n1_neighbors {
            if neighbor_set.contains(n2) {
                triangles += 1;
            }
        }
    }

    // Each triangle is counted twice in undirected graphs
    let pairs = k * (k - 1);
    Ok(triangles as f64 / pairs as f64)
}

/// Average clustering coefficient across all nodes.
pub fn average_clustering(graph: &Graph) -> Result<f64> {
    let ids = graph.node_ids();
    if ids.is_empty() {
        return Ok(0.0);
    }

    let sum: f64 = ids
        .iter()
        .map(|id| clustering_coefficient(graph, id).unwrap_or(0.0))
        .sum();

    Ok(sum / ids.len() as f64)
}

/// BFS shortest-path distances from a source node.
pub fn bfs_distances(graph: &Graph, source: &str) -> HashMap<String, usize> {
    let mut dist = HashMap::new();
    dist.insert(source.to_string(), 0usize);

    let mut queue = VecDeque::new();
    queue.push_back(source.to_string());

    while let Some(v) = queue.pop_front() {
        let d = dist[&v];
        if let Ok(neighbors) = graph.neighbors(&v) {
            for (w, _) in neighbors {
                if !dist.contains_key(w) {
                    dist.insert(w.to_string(), d + 1);
                    queue.push_back(w.to_string());
                }
            }
        }
    }

    dist
}

/// Diameter of the graph (longest shortest path).
///
/// Returns None if the graph is disconnected.
pub fn diameter(graph: &Graph) -> Result<Option<usize>> {
    let ids = graph.node_ids();
    let n = ids.len();
    if n <= 1 {
        return Ok(Some(0));
    }

    let mut max_dist = 0;
    for id in &ids {
        let dists = bfs_distances(graph, id);
        if dists.len() < n {
            return Ok(None); // disconnected
        }
        let local_max = dists.values().copied().max().unwrap_or(0);
        max_dist = max_dist.max(local_max);
    }

    Ok(Some(max_dist))
}

/// Connected components (for undirected graphs).
///
/// Returns a list of components, each being a set of node IDs.
pub fn connected_components(graph: &Graph) -> Vec<HashSet<String>> {
    let mut visited = HashSet::new();
    let mut components = Vec::new();

    for id in graph.node_ids() {
        if visited.contains(id) {
            continue;
        }
        let mut component = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(id.to_string());

        while let Some(v) = queue.pop_front() {
            if !visited.insert(v.clone()) {
                continue;
            }
            component.insert(v.clone());
            if let Ok(neighbors) = graph.neighbors(&v) {
                for (w, _) in neighbors {
                    if !visited.contains(w) {
                        queue.push_back(w.to_string());
                    }
                }
            }
        }

        components.push(component);
    }

    components
}

/// Shortest path between two nodes using BFS (unweighted).
///
/// Returns None if no path exists.
pub fn shortest_path(graph: &Graph, source: &str, target: &str) -> Result<Option<Vec<String>>> {
    if !graph.nodes.contains_key(source) {
        return Err(cyanea_core::CyaneaError::InvalidInput(
            format!("source '{}' not found", source),
        ));
    }
    if !graph.nodes.contains_key(target) {
        return Err(cyanea_core::CyaneaError::InvalidInput(
            format!("target '{}' not found", target),
        ));
    }

    if source == target {
        return Ok(Some(vec![source.to_string()]));
    }

    let mut visited = HashSet::new();
    let mut parent: HashMap<String, String> = HashMap::new();
    let mut queue = VecDeque::new();
    queue.push_back(source.to_string());
    visited.insert(source.to_string());

    while let Some(v) = queue.pop_front() {
        if let Ok(neighbors) = graph.neighbors(&v) {
            for (w, _) in neighbors {
                if !visited.contains(w) {
                    visited.insert(w.to_string());
                    parent.insert(w.to_string(), v.clone());
                    if w == target {
                        // Reconstruct path
                        let mut path = vec![target.to_string()];
                        let mut current = target.to_string();
                        while let Some(p) = parent.get(&current) {
                            path.push(p.clone());
                            current = p.clone();
                        }
                        path.reverse();
                        return Ok(Some(path));
                    }
                    queue.push_back(w.to_string());
                }
            }
        }
    }

    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::GraphType;

    fn triangle_graph() -> Graph {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        g.add_node("C", "C").unwrap();
        g.add_edge("A", "B", 1.0).unwrap();
        g.add_edge("B", "C", 1.0).unwrap();
        g.add_edge("A", "C", 1.0).unwrap();
        g
    }

    fn line_graph() -> Graph {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        g.add_node("C", "C").unwrap();
        g.add_node("D", "D").unwrap();
        g.add_edge("A", "B", 1.0).unwrap();
        g.add_edge("B", "C", 1.0).unwrap();
        g.add_edge("C", "D", 1.0).unwrap();
        g
    }

    #[test]
    fn test_degree_centrality() {
        let g = triangle_graph();
        let dc = degree_centrality(&g).unwrap();
        // All nodes have degree 2 / (3-1) = 1.0
        for (_, score) in &dc.scores {
            assert!((score - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_degree_centrality_line() {
        let g = line_graph();
        let dc = degree_centrality(&g).unwrap();
        // A and D have degree 1/(4-1) = 0.333
        assert!((dc.scores["A"] - 1.0 / 3.0).abs() < 0.01);
        // B and C have degree 2/(4-1) = 0.667
        assert!((dc.scores["B"] - 2.0 / 3.0).abs() < 0.01);
    }

    #[test]
    fn test_closeness_centrality() {
        let g = line_graph();
        let cc = closeness_centrality(&g).unwrap();
        // B: distances [0,1,1,2] → sum=4, reachable=3 → (3/4)*(3/3) = 0.75
        // Wait, it's to B: A=1, C=1, D=2 → sum=4
        assert!(cc.scores["B"] > cc.scores["A"]);
    }

    #[test]
    fn test_betweenness_centrality() {
        let g = line_graph();
        let bc = betweenness_centrality(&g).unwrap();
        // B and C should have highest betweenness
        assert!(bc.scores["B"] > bc.scores["A"]);
        assert!(bc.scores["C"] > bc.scores["D"]);
    }

    #[test]
    fn test_pagerank() {
        let g = triangle_graph();
        let pr = pagerank(&g, 0.85, 100, 1e-8).unwrap();
        // Symmetric graph: all ranks should be equal
        let vals: Vec<f64> = pr.scores.values().copied().collect();
        for v in &vals {
            assert!((v - 1.0 / 3.0).abs() < 0.01);
        }
    }

    #[test]
    fn test_pagerank_star() {
        // Star graph: center has highest PageRank
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("center", "center").unwrap();
        for i in 0..5 {
            let id = format!("leaf{}", i);
            g.add_node(&id, &id).unwrap();
            g.add_edge("center", &id, 1.0).unwrap();
        }
        let pr = pagerank(&g, 0.85, 100, 1e-8).unwrap();
        assert!(pr.scores["center"] > pr.scores["leaf0"]);
    }

    #[test]
    fn test_clustering_coefficient_triangle() {
        let g = triangle_graph();
        let cc = clustering_coefficient(&g, "A").unwrap();
        assert!((cc - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clustering_coefficient_line() {
        let g = line_graph();
        let cc = clustering_coefficient(&g, "B").unwrap();
        assert!((cc - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_average_clustering() {
        let g = triangle_graph();
        let avg = average_clustering(&g).unwrap();
        assert!((avg - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bfs_distances() {
        let g = line_graph();
        let dists = bfs_distances(&g, "A");
        assert_eq!(dists["A"], 0);
        assert_eq!(dists["B"], 1);
        assert_eq!(dists["C"], 2);
        assert_eq!(dists["D"], 3);
    }

    #[test]
    fn test_diameter_line() {
        let g = line_graph();
        assert_eq!(diameter(&g).unwrap(), Some(3));
    }

    #[test]
    fn test_diameter_triangle() {
        let g = triangle_graph();
        assert_eq!(diameter(&g).unwrap(), Some(1));
    }

    #[test]
    fn test_diameter_disconnected() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        // No edges — disconnected
        assert_eq!(diameter(&g).unwrap(), None);
    }

    #[test]
    fn test_connected_components() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        g.add_node("C", "C").unwrap();
        g.add_node("D", "D").unwrap();
        g.add_edge("A", "B", 1.0).unwrap();
        g.add_edge("C", "D", 1.0).unwrap();
        let comps = connected_components(&g);
        assert_eq!(comps.len(), 2);
    }

    #[test]
    fn test_shortest_path() {
        let g = line_graph();
        let path = shortest_path(&g, "A", "D").unwrap().unwrap();
        assert_eq!(path, vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn test_shortest_path_no_path() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        assert!(shortest_path(&g, "A", "B").unwrap().is_none());
    }

    #[test]
    fn test_shortest_path_same_node() {
        let g = line_graph();
        let path = shortest_path(&g, "A", "A").unwrap().unwrap();
        assert_eq!(path, vec!["A"]);
    }

    #[test]
    fn test_top_k() {
        let g = line_graph();
        let dc = degree_centrality(&g).unwrap();
        let top = dc.top_k(2);
        assert_eq!(top.len(), 2);
        // B and C should be in top 2
        let top_ids: Vec<&str> = top.iter().map(|(id, _)| *id).collect();
        assert!(top_ids.contains(&"B"));
        assert!(top_ids.contains(&"C"));
    }
}
