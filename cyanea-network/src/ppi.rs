//! Protein-protein interaction (PPI) network analysis.
//!
//! Tools for building, scoring, and analyzing PPI networks including
//! STRING-style combined scores, network propagation for guilt-by-association,
//! and hub/bottleneck identification.

use crate::graph::{Edge, Graph, GraphType};
use crate::topology;
use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

/// Evidence channels for PPI scoring (STRING-style).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PpiEvidence {
    /// Protein A identifier.
    pub protein_a: String,
    /// Protein B identifier.
    pub protein_b: String,
    /// Experimental evidence score (0.0–1.0).
    pub experimental: f64,
    /// Database/curated evidence score.
    pub database: f64,
    /// Text-mining evidence score.
    pub textmining: f64,
    /// Co-expression evidence score.
    pub coexpression: f64,
    /// Neighborhood (genomic context) evidence score.
    pub neighborhood: f64,
    /// Co-occurrence (phylogenetic profiles) evidence score.
    pub cooccurrence: f64,
}

impl PpiEvidence {
    pub fn new(protein_a: &str, protein_b: &str) -> Self {
        Self {
            protein_a: protein_a.to_string(),
            protein_b: protein_b.to_string(),
            experimental: 0.0,
            database: 0.0,
            textmining: 0.0,
            coexpression: 0.0,
            neighborhood: 0.0,
            cooccurrence: 0.0,
        }
    }

    /// Combine evidence channels using the STRING noisy-OR method.
    ///
    /// combined = 1 - Π(1 - s_i) for all non-zero channels.
    pub fn combined_score(&self) -> f64 {
        let channels = [
            self.experimental,
            self.database,
            self.textmining,
            self.coexpression,
            self.neighborhood,
            self.cooccurrence,
        ];

        let product: f64 = channels
            .iter()
            .filter(|&&s| s > 0.0)
            .map(|&s| 1.0 - s.clamp(0.0, 1.0))
            .product();

        1.0 - product
    }
}

/// Build a PPI graph from evidence data with a confidence threshold.
///
/// Only interactions with combined_score >= threshold are included.
pub fn build_ppi_network(
    evidence: &[PpiEvidence],
    threshold: f64,
) -> Result<Graph> {
    let mut graph = Graph::new(GraphType::Undirected);
    let mut seen_nodes = std::collections::HashSet::new();

    for ev in evidence {
        let score = ev.combined_score();
        if score < threshold {
            continue;
        }

        for protein in [&ev.protein_a, &ev.protein_b] {
            if seen_nodes.insert(protein.clone()) {
                graph.add_node(protein, protein)?;
            }
        }

        let edge = Edge::weighted(&ev.protein_a, &ev.protein_b, score)
            .with_type("ppi");
        graph.add_edge_struct(edge)?;
    }

    Ok(graph)
}

/// Network propagation using random walk with restart (RWR).
///
/// Propagates seed node scores through the network to find
/// functionally related proteins (guilt-by-association).
///
/// # Arguments
///
/// * `graph` - PPI network
/// * `seeds` - Node ID → seed score (e.g., disease-associated genes scored by p-value)
/// * `restart_prob` - Probability of returning to seed nodes (0.1–0.5 typical)
/// * `max_iter` - Maximum iterations
/// * `tolerance` - Convergence threshold
pub fn network_propagation(
    graph: &Graph,
    seeds: &HashMap<String, f64>,
    restart_prob: f64,
    max_iter: usize,
    tolerance: f64,
) -> Result<HashMap<String, f64>> {
    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let n = ids.len();

    if n == 0 {
        return Ok(HashMap::new());
    }

    // Normalize seed scores
    let seed_sum: f64 = seeds.values().sum();
    let seed_vec: HashMap<String, f64> = if seed_sum > 0.0 {
        seeds.iter().map(|(k, v)| (k.clone(), v / seed_sum)).collect()
    } else {
        HashMap::new()
    };

    // Initialize scores
    let mut scores: HashMap<String, f64> = ids
        .iter()
        .map(|id| (id.clone(), seed_vec.get(id).copied().unwrap_or(0.0)))
        .collect();

    // Build normalized adjacency
    let (sorted_ids, matrix) = graph.adjacency_matrix();
    let id_to_idx: HashMap<&str, usize> =
        sorted_ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

    // Row-normalize
    let mut norm_matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        let row_sum: f64 = matrix[i].iter().sum();
        if row_sum > 0.0 {
            for j in 0..n {
                norm_matrix[i][j] = matrix[i][j] / row_sum;
            }
        }
    }

    for _ in 0..max_iter {
        let mut new_scores = HashMap::new();

        for id in &sorted_ids {
            let i = id_to_idx[id.as_str()];
            // Random walk step
            let walk_score: f64 = (0..n)
                .map(|j| norm_matrix[j][i] * scores[&sorted_ids[j]])
                .sum();

            let seed_score = seed_vec.get(id).copied().unwrap_or(0.0);
            let new_score = (1.0 - restart_prob) * walk_score + restart_prob * seed_score;
            new_scores.insert(id.clone(), new_score);
        }

        // Check convergence
        let max_diff = sorted_ids
            .iter()
            .map(|id| (new_scores[id] - scores[id]).abs())
            .fold(0.0f64, f64::max);

        scores = new_scores;
        if max_diff < tolerance {
            break;
        }
    }

    Ok(scores)
}

/// Hub and bottleneck analysis for PPI networks.
///
/// Identifies hub proteins (high degree) and bottleneck proteins
/// (high betweenness centrality). Hub-bottleneck proteins are often
/// essential and functionally important.
#[derive(Debug, Clone)]
pub struct HubBottleneckResult {
    /// Node ID → degree.
    pub degrees: HashMap<String, usize>,
    /// Node ID → betweenness centrality.
    pub betweenness: HashMap<String, f64>,
    /// Nodes classified as hubs (degree > threshold).
    pub hubs: Vec<String>,
    /// Nodes classified as bottlenecks (betweenness > threshold).
    pub bottlenecks: Vec<String>,
    /// Nodes that are both hubs and bottlenecks.
    pub hub_bottlenecks: Vec<String>,
}

/// Identify hub and bottleneck proteins in a PPI network.
///
/// * `degree_percentile` - Nodes above this percentile of degree are hubs (e.g., 0.9 for top 10%)
/// * `betweenness_percentile` - Nodes above this percentile of betweenness are bottlenecks
pub fn hub_bottleneck_analysis(
    graph: &Graph,
    degree_percentile: f64,
    betweenness_percentile: f64,
) -> Result<HubBottleneckResult> {
    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();

    let degrees: HashMap<String, usize> = ids
        .iter()
        .map(|id| (id.clone(), graph.degree(id).unwrap_or(0)))
        .collect();

    let bc = topology::betweenness_centrality(graph)?;

    // Compute thresholds
    let mut deg_vals: Vec<usize> = degrees.values().copied().collect();
    deg_vals.sort();
    let deg_threshold = if deg_vals.is_empty() {
        0
    } else {
        let idx = ((deg_vals.len() as f64) * degree_percentile) as usize;
        deg_vals[idx.min(deg_vals.len() - 1)]
    };

    let mut bc_vals: Vec<f64> = bc.scores.values().copied().collect();
    bc_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let bc_threshold = if bc_vals.is_empty() {
        0.0
    } else {
        let idx = ((bc_vals.len() as f64) * betweenness_percentile) as usize;
        bc_vals[idx.min(bc_vals.len() - 1)]
    };

    let hubs: Vec<String> = degrees
        .iter()
        .filter(|(_, &d)| d >= deg_threshold && deg_threshold > 0)
        .map(|(id, _)| id.clone())
        .collect();

    let bottlenecks: Vec<String> = bc
        .scores
        .iter()
        .filter(|(_, &b)| b >= bc_threshold && bc_threshold > 0.0)
        .map(|(id, _)| id.clone())
        .collect();

    let hub_set: std::collections::HashSet<&str> = hubs.iter().map(|s| s.as_str()).collect();
    let hub_bottlenecks: Vec<String> = bottlenecks
        .iter()
        .filter(|id| hub_set.contains(id.as_str()))
        .cloned()
        .collect();

    Ok(HubBottleneckResult {
        degrees,
        betweenness: bc.scores,
        hubs,
        bottlenecks,
        hub_bottlenecks,
    })
}

/// Compute the Jaccard similarity between neighborhoods of two nodes.
///
/// Used to predict potential interactions: nodes with similar neighborhoods
/// are more likely to interact.
pub fn neighborhood_similarity(graph: &Graph, node_a: &str, node_b: &str) -> Result<f64> {
    let neighbors_a: std::collections::HashSet<String> = graph
        .neighbors(node_a)?
        .iter()
        .map(|(id, _)| id.to_string())
        .collect();

    let neighbors_b: std::collections::HashSet<String> = graph
        .neighbors(node_b)?
        .iter()
        .map(|(id, _)| id.to_string())
        .collect();

    let intersection = neighbors_a.intersection(&neighbors_b).count();
    let union = neighbors_a.union(&neighbors_b).count();

    if union == 0 {
        Ok(0.0)
    } else {
        Ok(intersection as f64 / union as f64)
    }
}

/// Link prediction scores for all non-existing edges.
///
/// Returns (node_a, node_b, Jaccard similarity) tuples sorted by score descending.
pub fn predict_interactions(graph: &Graph, top_k: usize) -> Result<Vec<(String, String, f64)>> {
    let ids: Vec<String> = graph.node_ids().into_iter().map(|s| s.to_string()).collect();
    let mut predictions = Vec::new();

    for i in 0..ids.len() {
        for j in (i + 1)..ids.len() {
            if !graph.has_edge(&ids[i], &ids[j]) {
                let sim = neighborhood_similarity(graph, &ids[i], &ids[j])?;
                if sim > 0.0 {
                    predictions.push((ids[i].clone(), ids[j].clone(), sim));
                }
            }
        }
    }

    predictions.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
    predictions.truncate(top_k);
    Ok(predictions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_combined_score_single_channel() {
        let mut ev = PpiEvidence::new("P1", "P2");
        ev.experimental = 0.9;
        assert!((ev.combined_score() - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_combined_score_multiple() {
        let mut ev = PpiEvidence::new("P1", "P2");
        ev.experimental = 0.7;
        ev.database = 0.8;
        // 1 - (1-0.7)*(1-0.8) = 1 - 0.3*0.2 = 1 - 0.06 = 0.94
        assert!((ev.combined_score() - 0.94).abs() < 1e-10);
    }

    #[test]
    fn test_combined_score_zero() {
        let ev = PpiEvidence::new("P1", "P2");
        assert!((ev.combined_score() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_ppi_network() {
        let mut ev1 = PpiEvidence::new("TP53", "MDM2");
        ev1.experimental = 0.95;
        let mut ev2 = PpiEvidence::new("TP53", "BAX");
        ev2.experimental = 0.8;
        let mut ev3 = PpiEvidence::new("LOW1", "LOW2");
        ev3.textmining = 0.1;

        let graph = build_ppi_network(&[ev1, ev2, ev3], 0.5).unwrap();
        assert_eq!(graph.node_count(), 3); // TP53, MDM2, BAX (LOW1/LOW2 filtered)
        assert_eq!(graph.edge_count(), 2);
    }

    #[test]
    fn test_build_ppi_empty_threshold() {
        let mut ev = PpiEvidence::new("A", "B");
        ev.experimental = 0.3;
        let graph = build_ppi_network(&[ev], 0.5).unwrap();
        assert_eq!(graph.node_count(), 0);
    }

    #[test]
    fn test_network_propagation() {
        let mut graph = Graph::new(GraphType::Undirected);
        for id in &["A", "B", "C", "D", "E"] {
            graph.add_node(id, id).unwrap();
        }
        graph.add_edge("A", "B", 1.0).unwrap();
        graph.add_edge("B", "C", 1.0).unwrap();
        graph.add_edge("C", "D", 1.0).unwrap();
        graph.add_edge("D", "E", 1.0).unwrap();

        let seeds: HashMap<String, f64> = vec![("A".into(), 1.0)].into_iter().collect();
        let scores = network_propagation(&graph, &seeds, 0.3, 100, 1e-8).unwrap();

        // A should have highest score, then B, C, D, E
        assert!(scores["A"] > scores["B"]);
        assert!(scores["B"] > scores["C"]);
        assert!(scores["C"] > scores["D"]);
    }

    #[test]
    fn test_network_propagation_two_seeds() {
        let mut graph = Graph::new(GraphType::Undirected);
        for id in &["A", "B", "C"] {
            graph.add_node(id, id).unwrap();
        }
        graph.add_edge("A", "B", 1.0).unwrap();
        graph.add_edge("B", "C", 1.0).unwrap();

        let seeds: HashMap<String, f64> =
            vec![("A".into(), 1.0), ("C".into(), 1.0)].into_iter().collect();
        let scores = network_propagation(&graph, &seeds, 0.3, 100, 1e-8).unwrap();

        // B is between both seeds, should receive propagated signal
        assert!(scores["B"] > 0.0);
    }

    #[test]
    fn test_hub_bottleneck() {
        // Star graph: center is hub
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("center", "center").unwrap();
        for i in 0..6 {
            let id = format!("leaf{}", i);
            g.add_node(&id, &id).unwrap();
            g.add_edge("center", &id, 1.0).unwrap();
        }

        let result = hub_bottleneck_analysis(&g, 0.8, 0.8).unwrap();
        assert!(result.hubs.contains(&"center".to_string()));
    }

    #[test]
    fn test_neighborhood_similarity() {
        let mut g = Graph::new(GraphType::Undirected);
        for id in &["A", "B", "C", "D", "E"] {
            g.add_node(id, id).unwrap();
        }
        g.add_edge("A", "C", 1.0).unwrap();
        g.add_edge("A", "D", 1.0).unwrap();
        g.add_edge("B", "C", 1.0).unwrap();
        g.add_edge("B", "D", 1.0).unwrap();
        g.add_edge("B", "E", 1.0).unwrap();

        let sim = neighborhood_similarity(&g, "A", "B").unwrap();
        // A neighbors: {C, D}, B neighbors: {C, D, E}
        // Jaccard = 2 / 3 = 0.667
        assert!((sim - 2.0 / 3.0).abs() < 0.01);
    }

    #[test]
    fn test_predict_interactions() {
        let mut g = Graph::new(GraphType::Undirected);
        for id in &["A", "B", "C", "D"] {
            g.add_node(id, id).unwrap();
        }
        g.add_edge("A", "C", 1.0).unwrap();
        g.add_edge("A", "D", 1.0).unwrap();
        g.add_edge("B", "C", 1.0).unwrap();
        g.add_edge("B", "D", 1.0).unwrap();

        let preds = predict_interactions(&g, 5).unwrap();
        // A-B and C-D should both be predicted (shared neighbors)
        assert!(!preds.is_empty());
        let pairs: Vec<(&str, &str)> = preds.iter().map(|(a, b, _)| (a.as_str(), b.as_str())).collect();
        let has_ab = pairs.iter().any(|&(a, b)| (a == "A" && b == "B") || (a == "B" && b == "A"));
        assert!(has_ab, "A-B should be predicted");
    }
}
