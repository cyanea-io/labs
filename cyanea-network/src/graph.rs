//! Core graph data structure for biological networks.
//!
//! Provides a flexible adjacency-list graph supporting directed and undirected
//! edges, weighted edges, and arbitrary string-keyed node/edge attributes.

use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

/// Whether a graph is directed or undirected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GraphType {
    Directed,
    Undirected,
}

/// A node in the graph.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Node {
    /// Unique node identifier.
    pub id: String,
    /// Human-readable label (e.g., gene symbol, protein name).
    pub label: String,
    /// Key-value attributes (e.g., "type" → "kinase", "organism" → "human").
    pub attributes: HashMap<String, String>,
}

impl Node {
    pub fn new(id: &str, label: &str) -> Self {
        Self {
            id: id.to_string(),
            label: label.to_string(),
            attributes: HashMap::new(),
        }
    }

    pub fn with_attr(mut self, key: &str, value: &str) -> Self {
        self.attributes.insert(key.to_string(), value.to_string());
        self
    }
}

/// An edge in the graph.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Edge {
    /// Source node ID.
    pub source: String,
    /// Target node ID.
    pub target: String,
    /// Edge weight (default 1.0).
    pub weight: f64,
    /// Edge type / interaction type (e.g., "activation", "inhibition", "binding").
    pub edge_type: Option<String>,
    /// Key-value attributes.
    pub attributes: HashMap<String, String>,
}

impl Edge {
    pub fn new(source: &str, target: &str) -> Self {
        Self {
            source: source.to_string(),
            target: target.to_string(),
            weight: 1.0,
            edge_type: None,
            attributes: HashMap::new(),
        }
    }

    pub fn weighted(source: &str, target: &str, weight: f64) -> Self {
        Self {
            source: source.to_string(),
            target: target.to_string(),
            weight,
            edge_type: None,
            attributes: HashMap::new(),
        }
    }

    pub fn with_type(mut self, edge_type: &str) -> Self {
        self.edge_type = Some(edge_type.to_string());
        self
    }

    pub fn with_attr(mut self, key: &str, value: &str) -> Self {
        self.attributes.insert(key.to_string(), value.to_string());
        self
    }
}

/// Adjacency-list graph for biological networks.
///
/// Supports directed and undirected graphs with weighted edges and attributes.
/// Nodes are identified by string IDs. For undirected graphs, edges are stored
/// in both directions internally.
///
/// # Example
///
/// ```
/// use cyanea_network::Graph;
/// use cyanea_network::graph::GraphType;
///
/// let mut g = Graph::new(GraphType::Undirected);
/// g.add_node("A", "Gene A").unwrap();
/// g.add_node("B", "Gene B").unwrap();
/// g.add_edge("A", "B", 0.95).unwrap();
///
/// assert_eq!(g.node_count(), 2);
/// assert_eq!(g.edge_count(), 1);
/// assert_eq!(g.neighbors("A").unwrap().len(), 1);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Graph {
    /// Directed or undirected.
    pub graph_type: GraphType,
    /// Node ID → Node.
    pub nodes: HashMap<String, Node>,
    /// Source ID → list of (target ID, edge index).
    adjacency: HashMap<String, Vec<(String, usize)>>,
    /// All edges (canonical storage).
    pub edges: Vec<Edge>,
    /// Graph-level attributes.
    pub attributes: HashMap<String, String>,
}

impl Graph {
    /// Create a new empty graph.
    pub fn new(graph_type: GraphType) -> Self {
        Self {
            graph_type,
            nodes: HashMap::new(),
            adjacency: HashMap::new(),
            edges: Vec::new(),
            attributes: HashMap::new(),
        }
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges (canonical count — undirected edges counted once).
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Add a node. Returns error if node ID already exists.
    pub fn add_node(&mut self, id: &str, label: &str) -> Result<()> {
        if self.nodes.contains_key(id) {
            return Err(CyaneaError::InvalidInput(
                format!("node '{}' already exists", id),
            ));
        }
        self.nodes.insert(id.to_string(), Node::new(id, label));
        self.adjacency.entry(id.to_string()).or_default();
        Ok(())
    }

    /// Add a node with a [`Node`] struct.
    pub fn add_node_struct(&mut self, node: Node) -> Result<()> {
        if self.nodes.contains_key(&node.id) {
            return Err(CyaneaError::InvalidInput(
                format!("node '{}' already exists", node.id),
            ));
        }
        let id = node.id.clone();
        self.nodes.insert(id.clone(), node);
        self.adjacency.entry(id).or_default();
        Ok(())
    }

    /// Add a weighted edge between two existing nodes.
    pub fn add_edge(&mut self, source: &str, target: &str, weight: f64) -> Result<()> {
        if !self.nodes.contains_key(source) {
            return Err(CyaneaError::InvalidInput(
                format!("source node '{}' not found", source),
            ));
        }
        if !self.nodes.contains_key(target) {
            return Err(CyaneaError::InvalidInput(
                format!("target node '{}' not found", target),
            ));
        }

        let edge_idx = self.edges.len();
        self.edges.push(Edge::weighted(source, target, weight));

        self.adjacency
            .entry(source.to_string())
            .or_default()
            .push((target.to_string(), edge_idx));

        if self.graph_type == GraphType::Undirected {
            self.adjacency
                .entry(target.to_string())
                .or_default()
                .push((source.to_string(), edge_idx));
        }

        Ok(())
    }

    /// Add an edge with an [`Edge`] struct.
    pub fn add_edge_struct(&mut self, edge: Edge) -> Result<()> {
        if !self.nodes.contains_key(&edge.source) {
            return Err(CyaneaError::InvalidInput(
                format!("source node '{}' not found", edge.source),
            ));
        }
        if !self.nodes.contains_key(&edge.target) {
            return Err(CyaneaError::InvalidInput(
                format!("target node '{}' not found", edge.target),
            ));
        }

        let edge_idx = self.edges.len();
        let source = edge.source.clone();
        let target = edge.target.clone();
        self.edges.push(edge);

        self.adjacency
            .entry(source.clone())
            .or_default()
            .push((target.clone(), edge_idx));

        if self.graph_type == GraphType::Undirected {
            self.adjacency
                .entry(target)
                .or_default()
                .push((source, edge_idx));
        }

        Ok(())
    }

    /// Get neighbors of a node as (neighbor_id, edge_weight) pairs.
    pub fn neighbors(&self, node_id: &str) -> Result<Vec<(&str, f64)>> {
        let adj = self.adjacency.get(node_id).ok_or_else(|| {
            CyaneaError::InvalidInput(format!("node '{}' not found", node_id))
        })?;

        Ok(adj
            .iter()
            .map(|(target, edge_idx)| (target.as_str(), self.edges[*edge_idx].weight))
            .collect())
    }

    /// Get the degree of a node (number of edges).
    pub fn degree(&self, node_id: &str) -> Result<usize> {
        let adj = self.adjacency.get(node_id).ok_or_else(|| {
            CyaneaError::InvalidInput(format!("node '{}' not found", node_id))
        })?;
        Ok(adj.len())
    }

    /// Get in-degree for directed graphs. Same as degree for undirected.
    pub fn in_degree(&self, node_id: &str) -> Result<usize> {
        if !self.nodes.contains_key(node_id) {
            return Err(CyaneaError::InvalidInput(
                format!("node '{}' not found", node_id),
            ));
        }
        if self.graph_type == GraphType::Undirected {
            return self.degree(node_id);
        }
        let count = self.edges.iter().filter(|e| e.target == node_id).count();
        Ok(count)
    }

    /// Get out-degree for directed graphs. Same as degree for undirected.
    pub fn out_degree(&self, node_id: &str) -> Result<usize> {
        if self.graph_type == GraphType::Undirected {
            return self.degree(node_id);
        }
        self.degree(node_id)
    }

    /// Check if an edge exists between two nodes.
    pub fn has_edge(&self, source: &str, target: &str) -> bool {
        self.adjacency
            .get(source)
            .map(|adj| adj.iter().any(|(t, _)| t == target))
            .unwrap_or(false)
    }

    /// Get edge weight between two nodes. Returns None if no edge exists.
    pub fn edge_weight(&self, source: &str, target: &str) -> Option<f64> {
        self.adjacency.get(source).and_then(|adj| {
            adj.iter()
                .find(|(t, _)| t == target)
                .map(|(_, idx)| self.edges[*idx].weight)
        })
    }

    /// Get all node IDs.
    pub fn node_ids(&self) -> Vec<&str> {
        self.nodes.keys().map(|s| s.as_str()).collect()
    }

    /// Get a node by ID.
    pub fn get_node(&self, id: &str) -> Option<&Node> {
        self.nodes.get(id)
    }

    /// Remove a node and all its edges.
    pub fn remove_node(&mut self, id: &str) -> Result<()> {
        if !self.nodes.contains_key(id) {
            return Err(CyaneaError::InvalidInput(
                format!("node '{}' not found", id),
            ));
        }

        // Remove edges involving this node
        self.edges.retain(|e| e.source != id && e.target != id);

        // Rebuild adjacency
        self.adjacency.clear();
        for (idx, edge) in self.edges.iter().enumerate() {
            self.adjacency
                .entry(edge.source.clone())
                .or_default()
                .push((edge.target.clone(), idx));
            if self.graph_type == GraphType::Undirected {
                self.adjacency
                    .entry(edge.target.clone())
                    .or_default()
                    .push((edge.source.clone(), idx));
            }
        }

        self.nodes.remove(id);
        self.adjacency.remove(id);

        Ok(())
    }

    /// Create a subgraph containing only the specified nodes and edges between them.
    pub fn subgraph(&self, node_ids: &[&str]) -> Result<Graph> {
        let mut sub = Graph::new(self.graph_type);
        let node_set: std::collections::HashSet<&str> = node_ids.iter().copied().collect();

        for id in node_ids {
            if let Some(node) = self.nodes.get(*id) {
                sub.add_node_struct(node.clone())?;
            }
        }

        for edge in &self.edges {
            if node_set.contains(edge.source.as_str())
                && node_set.contains(edge.target.as_str())
            {
                sub.add_edge_struct(edge.clone())?;
            }
        }

        sub.attributes = self.attributes.clone();
        Ok(sub)
    }

    /// Graph density: |E| / (|V| * (|V| - 1)) for directed, doubled for undirected.
    pub fn density(&self) -> f64 {
        let n = self.node_count() as f64;
        if n < 2.0 {
            return 0.0;
        }
        let e = self.edge_count() as f64;
        match self.graph_type {
            GraphType::Directed => e / (n * (n - 1.0)),
            GraphType::Undirected => 2.0 * e / (n * (n - 1.0)),
        }
    }

    /// Return the adjacency matrix as a Vec<Vec<f64>> with consistent node ordering.
    /// Returns (node_id_order, matrix).
    pub fn adjacency_matrix(&self) -> (Vec<String>, Vec<Vec<f64>>) {
        let mut ids: Vec<String> = self.nodes.keys().cloned().collect();
        ids.sort();
        let n = ids.len();
        let id_to_idx: HashMap<&str, usize> =
            ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

        let mut matrix = vec![vec![0.0; n]; n];
        for edge in &self.edges {
            if let (Some(&i), Some(&j)) =
                (id_to_idx.get(edge.source.as_str()), id_to_idx.get(edge.target.as_str()))
            {
                matrix[i][j] = edge.weight;
                if self.graph_type == GraphType::Undirected {
                    matrix[j][i] = edge.weight;
                }
            }
        }

        (ids, matrix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_undirected() -> Graph {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "Gene A").unwrap();
        g.add_node("B", "Gene B").unwrap();
        g.add_node("C", "Gene C").unwrap();
        g.add_node("D", "Gene D").unwrap();
        g.add_edge("A", "B", 1.0).unwrap();
        g.add_edge("B", "C", 0.8).unwrap();
        g.add_edge("C", "D", 0.6).unwrap();
        g.add_edge("A", "D", 0.5).unwrap();
        g
    }

    fn sample_directed() -> Graph {
        let mut g = Graph::new(GraphType::Directed);
        g.add_node("TF1", "Transcription Factor 1").unwrap();
        g.add_node("TF2", "Transcription Factor 2").unwrap();
        g.add_node("G1", "Gene 1").unwrap();
        g.add_node("G2", "Gene 2").unwrap();
        g.add_edge("TF1", "G1", 0.9).unwrap();
        g.add_edge("TF1", "G2", 0.7).unwrap();
        g.add_edge("TF2", "G1", 0.5).unwrap();
        g
    }

    #[test]
    fn test_new_graph() {
        let g = Graph::new(GraphType::Undirected);
        assert_eq!(g.node_count(), 0);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn test_add_nodes_and_edges() {
        let g = sample_undirected();
        assert_eq!(g.node_count(), 4);
        assert_eq!(g.edge_count(), 4);
    }

    #[test]
    fn test_duplicate_node() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        assert!(g.add_node("A", "A2").is_err());
    }

    #[test]
    fn test_edge_missing_node() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        assert!(g.add_edge("A", "Z", 1.0).is_err());
    }

    #[test]
    fn test_neighbors_undirected() {
        let g = sample_undirected();
        let neighbors = g.neighbors("A").unwrap();
        assert_eq!(neighbors.len(), 2);
        let neighbor_ids: Vec<&str> = neighbors.iter().map(|(id, _)| *id).collect();
        assert!(neighbor_ids.contains(&"B"));
        assert!(neighbor_ids.contains(&"D"));
    }

    #[test]
    fn test_neighbors_directed() {
        let g = sample_directed();
        let neighbors = g.neighbors("TF1").unwrap();
        assert_eq!(neighbors.len(), 2);
        let neighbors_g1 = g.neighbors("G1").unwrap();
        assert_eq!(neighbors_g1.len(), 0); // G1 has no outgoing edges
    }

    #[test]
    fn test_degree() {
        let g = sample_undirected();
        assert_eq!(g.degree("A").unwrap(), 2);
        assert_eq!(g.degree("B").unwrap(), 2);
    }

    #[test]
    fn test_in_degree_directed() {
        let g = sample_directed();
        assert_eq!(g.in_degree("G1").unwrap(), 2);
        assert_eq!(g.in_degree("TF1").unwrap(), 0);
    }

    #[test]
    fn test_has_edge() {
        let g = sample_undirected();
        assert!(g.has_edge("A", "B"));
        assert!(g.has_edge("B", "A")); // undirected
        assert!(!g.has_edge("A", "C"));
    }

    #[test]
    fn test_has_edge_directed() {
        let g = sample_directed();
        assert!(g.has_edge("TF1", "G1"));
        assert!(!g.has_edge("G1", "TF1"));
    }

    #[test]
    fn test_edge_weight() {
        let g = sample_undirected();
        assert_eq!(g.edge_weight("A", "B"), Some(1.0));
        assert_eq!(g.edge_weight("A", "C"), None);
    }

    #[test]
    fn test_density() {
        let g = sample_undirected();
        // 4 nodes, 4 edges, undirected: 2*4 / (4*3) = 8/12 ≈ 0.667
        assert!((g.density() - 0.6667).abs() < 0.01);
    }

    #[test]
    fn test_density_directed() {
        let g = sample_directed();
        // 4 nodes, 3 edges, directed: 3 / (4*3) = 0.25
        assert!((g.density() - 0.25).abs() < 0.01);
    }

    #[test]
    fn test_remove_node() {
        let mut g = sample_undirected();
        g.remove_node("B").unwrap();
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 2); // A-D and C-D remain
        assert!(!g.has_edge("A", "B"));
    }

    #[test]
    fn test_subgraph() {
        let g = sample_undirected();
        let sub = g.subgraph(&["A", "B", "C"]).unwrap();
        assert_eq!(sub.node_count(), 3);
        assert_eq!(sub.edge_count(), 2); // A-B and B-C
    }

    #[test]
    fn test_adjacency_matrix() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        g.add_edge("A", "B", 2.5).unwrap();
        let (ids, mat) = g.adjacency_matrix();
        let ai = ids.iter().position(|s| s == "A").unwrap();
        let bi = ids.iter().position(|s| s == "B").unwrap();
        assert!((mat[ai][bi] - 2.5).abs() < 1e-10);
        assert!((mat[bi][ai] - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_node_with_attributes() {
        let mut g = Graph::new(GraphType::Undirected);
        let node = Node::new("P53", "TP53").with_attr("type", "tumor_suppressor");
        g.add_node_struct(node).unwrap();
        let n = g.get_node("P53").unwrap();
        assert_eq!(n.attributes.get("type").unwrap(), "tumor_suppressor");
    }

    #[test]
    fn test_edge_with_type() {
        let mut g = Graph::new(GraphType::Directed);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        let edge = Edge::new("A", "B").with_type("activation").with_attr("source", "STRING");
        g.add_edge_struct(edge).unwrap();
        assert_eq!(g.edges[0].edge_type.as_deref(), Some("activation"));
        assert_eq!(g.edges[0].attributes.get("source").unwrap(), "STRING");
    }
}
