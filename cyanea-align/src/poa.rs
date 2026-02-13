//! Partial Order Alignment (POA) for multi-sequence consensus.
//!
//! Implements the POA algorithm described in:
//! C. Lee, C. Grasso, M. F. Sharlow, "Multiple sequence alignment using
//! partial order graphs", *Bioinformatics*, 18(3):452--464, 2002.
//!
//! Sequences are progressively aligned to a DAG and integrated, building a
//! partial order graph that represents the multiple alignment. Consensus is
//! extracted by finding the heaviest weighted path through the DAG.
//!
//! # Example
//!
//! ```
//! use cyanea_align::poa::{PoaGraph, PoaScoring};
//!
//! let mut graph = PoaGraph::from_sequence(b"ACGTACGT");
//! let scoring = PoaScoring::default();
//! graph.add_sequence(b"ACGTACGT", &scoring).unwrap();
//! graph.add_sequence(b"ACGAACGT", &scoring).unwrap();
//! let consensus = graph.consensus();
//! assert_eq!(consensus, b"ACGTACGT");
//! ```

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// A node in the POA graph.
#[derive(Debug, Clone)]
struct PoaNode {
    /// The base at this node.
    base: u8,
    /// Indices of successor nodes (outgoing edges).
    successors: Vec<usize>,
    /// Indices of predecessor nodes (incoming edges).
    predecessors: Vec<usize>,
    /// Edge weights to successors (parallel to `successors` vec).
    weights: Vec<usize>,
    /// Alignment info: which nodes from other sequences this aligns to.
    aligned_to: Vec<usize>,
}

impl PoaNode {
    fn new(base: u8) -> Self {
        Self {
            base,
            successors: Vec::new(),
            predecessors: Vec::new(),
            weights: Vec::new(),
            aligned_to: Vec::new(),
        }
    }
}

/// Partial Order Alignment graph for multi-sequence consensus.
///
/// Implements the Lee 2002 algorithm: sequences are progressively aligned
/// to the graph and integrated, building a DAG that represents the
/// multiple alignment.
#[derive(Debug, Clone)]
pub struct PoaGraph {
    nodes: Vec<PoaNode>,
    /// Start nodes (no predecessors) for topological traversal.
    start_nodes: Vec<usize>,
    /// Number of sequences that have been added to the graph.
    n_sequences: usize,
}

/// Scoring parameters for POA.
#[derive(Debug, Clone, Copy)]
pub struct PoaScoring {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_score: i32,
}

impl Default for PoaScoring {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_score: -2,
        }
    }
}

// ---------------------------------------------------------------------------
// Traceback direction constants
// ---------------------------------------------------------------------------

const TB_NONE: u8 = 0;
const TB_DIAG: u8 = 1;
const TB_UP: u8 = 2; // gap in sequence (consume graph node)
const TB_LEFT: u8 = 3; // gap in graph (consume sequence position)

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

impl PoaGraph {
    /// Create an empty POA graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            start_nodes: Vec::new(),
            n_sequences: 0,
        }
    }

    /// Initialize a POA graph from the first sequence (linear chain of nodes).
    pub fn from_sequence(seq: &[u8]) -> Self {
        let mut graph = Self::new();
        if seq.is_empty() {
            return graph;
        }

        for (i, &base) in seq.iter().enumerate() {
            let node = PoaNode::new(base);
            graph.nodes.push(node);

            if i > 0 {
                graph.add_edge(i - 1, i, 1);
            }
        }

        graph.start_nodes.push(0);
        graph.n_sequences = 1;
        graph
    }

    /// Add a directed edge from `from` to `to` with the given weight.
    /// If the edge already exists, increment its weight.
    fn add_edge(&mut self, from: usize, to: usize, weight: usize) {
        if let Some(pos) = self.nodes[from].successors.iter().position(|&s| s == to) {
            self.nodes[from].weights[pos] += weight;
        } else {
            self.nodes[from].successors.push(to);
            self.nodes[from].weights.push(weight);
            self.nodes[to].predecessors.push(from);
        }
    }

    /// Align a sequence to the graph using DP (Lee 2002), then integrate
    /// the aligned sequence into the graph.
    ///
    /// # Errors
    ///
    /// Returns an error if the sequence is empty or the graph is empty.
    pub fn add_sequence(&mut self, seq: &[u8], scoring: &PoaScoring) -> Result<()> {
        if seq.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "sequence must not be empty".into(),
            ));
        }
        if self.nodes.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "graph is empty; use from_sequence to initialize".into(),
            ));
        }

        let topo = self.topological_sort();
        let alignment = self.align_to_graph(seq, scoring, &topo);
        self.integrate_alignment(seq, &alignment);
        self.n_sequences += 1;

        Ok(())
    }

    /// Align a sequence to the graph via DP over the topological ordering.
    ///
    /// Returns a list of `(Option<graph_node>, Option<seq_pos>)` pairs
    /// representing the alignment from traceback.
    fn align_to_graph(
        &self,
        seq: &[u8],
        scoring: &PoaScoring,
        topo: &[usize],
    ) -> Vec<(Option<usize>, Option<usize>)> {
        let n_nodes = self.nodes.len();
        let seq_len = seq.len();

        // Map from node index to topological rank for fast lookup.
        let mut topo_rank = vec![0usize; n_nodes];
        for (rank, &node_idx) in topo.iter().enumerate() {
            topo_rank[node_idx] = rank;
        }

        // DP score matrix: score[topo_rank][seq_pos + 1]
        // Row 0 is the "before any graph node" row.
        // Col 0 is the "before any sequence base" column.
        let rows = topo.len() + 1;
        let cols = seq_len + 1;
        let mut score = vec![vec![i32::MIN / 2; cols]; rows];
        let mut trace = vec![vec![(TB_NONE, 0usize); cols]; rows];

        // Initialize: starting before any node or base is score 0.
        score[0][0] = 0;

        // Initialize first row: gaps in the graph (consuming sequence bases).
        for j in 1..cols {
            score[0][j] = scoring.gap_score * j as i32;
            trace[0][j] = (TB_LEFT, 0);
        }

        // For each graph node in topological order.
        for (rank_idx, &node_idx) in topo.iter().enumerate() {
            let rank = rank_idx + 1; // +1 because row 0 is the "no node" row
            let node_base = self.nodes[node_idx].base;

            // Find the best score arriving at this node from predecessors
            // (or from the start if this is a start node).
            let pred_ranks: Vec<usize> = if self.nodes[node_idx].predecessors.is_empty() {
                vec![0] // Use the "no node" row
            } else {
                self.nodes[node_idx]
                    .predecessors
                    .iter()
                    .map(|&p| topo_rank[p] + 1)
                    .collect()
            };

            // Column 0: gap in the sequence (consuming graph node).
            for &pr in &pred_ranks {
                let gap_val = score[pr][0] + scoring.gap_score;
                if gap_val > score[rank][0] {
                    score[rank][0] = gap_val;
                    trace[rank][0] = (TB_UP, pr);
                }
            }

            // Fill DP for each sequence position.
            for j in 1..cols {
                let seq_base = seq[j - 1];

                // Diagonal: match/mismatch from a predecessor row.
                let sub = if seq_base.to_ascii_uppercase() == node_base.to_ascii_uppercase() {
                    scoring.match_score
                } else {
                    scoring.mismatch_score
                };
                for &pr in &pred_ranks {
                    let diag = score[pr][j - 1] + sub;
                    if diag > score[rank][j] {
                        score[rank][j] = diag;
                        trace[rank][j] = (TB_DIAG, pr);
                    }
                }

                // Up: gap in sequence (consume graph node, skip seq base).
                for &pr in &pred_ranks {
                    let up = score[pr][j] + scoring.gap_score;
                    if up > score[rank][j] {
                        score[rank][j] = up;
                        trace[rank][j] = (TB_UP, pr);
                    }
                }

                // Left: gap in graph (consume seq base, stay at same graph position).
                let left = score[rank][j - 1] + scoring.gap_score;
                if left > score[rank][j] {
                    score[rank][j] = left;
                    trace[rank][j] = (TB_LEFT, rank);
                }
            }
        }

        // Find the best ending position: last column of any node row.
        let mut best_score = i32::MIN;
        let mut best_rank = 0;
        for rank in 1..rows {
            if score[rank][seq_len] > best_score {
                best_score = score[rank][seq_len];
                best_rank = rank;
            }
        }

        // Traceback.
        let mut alignment = Vec::new();
        let mut r = best_rank;
        let mut c = seq_len;

        while r > 0 || c > 0 {
            let (dir, prev_rank) = trace[r][c];

            match dir {
                TB_DIAG => {
                    let node_idx = topo[r - 1];
                    alignment.push((Some(node_idx), Some(c - 1)));
                    r = prev_rank;
                    c -= 1;
                }
                TB_UP => {
                    let node_idx = topo[r - 1];
                    alignment.push((Some(node_idx), None));
                    r = prev_rank;
                }
                TB_LEFT => {
                    alignment.push((None, Some(c - 1)));
                    c -= 1;
                }
                _ => {
                    // TB_NONE at (0, 0) — we are done.
                    break;
                }
            }
        }

        alignment.reverse();
        alignment
    }

    /// Integrate an aligned sequence into the graph.
    ///
    /// For matched positions where bases agree, we reuse the existing node and
    /// increase edge weight. For mismatches, we add a new node and record it
    /// in `aligned_to`. For insertions (sequence bases with no graph node),
    /// we add new nodes to the graph.
    fn integrate_alignment(
        &mut self,
        seq: &[u8],
        alignment: &[(Option<usize>, Option<usize>)],
    ) {
        let mut prev_node: Option<usize> = None;

        for &(graph_node, seq_pos) in alignment {
            match (graph_node, seq_pos) {
                (Some(gn), Some(sp)) => {
                    // Aligned pair: graph node gn aligns to sequence position sp.
                    let seq_base = seq[sp];
                    let node_base = self.nodes[gn].base;

                    if seq_base.to_ascii_uppercase() == node_base.to_ascii_uppercase() {
                        // Match: reuse the node, just add/increase edge from prev.
                        if let Some(prev) = prev_node {
                            self.add_edge(prev, gn, 1);
                        } else if !self.start_nodes.contains(&gn) {
                            // This node becomes a start node for this sequence.
                            self.start_nodes.push(gn);
                        }
                        prev_node = Some(gn);
                    } else {
                        // Mismatch: add a new node with the sequence base.
                        let new_idx = self.nodes.len();
                        self.nodes.push(PoaNode::new(seq_base));

                        // Record alignment relationship.
                        self.nodes[gn].aligned_to.push(new_idx);
                        self.nodes[new_idx].aligned_to.push(gn);

                        if let Some(prev) = prev_node {
                            self.add_edge(prev, new_idx, 1);
                        } else {
                            self.start_nodes.push(new_idx);
                        }
                        prev_node = Some(new_idx);
                    }
                }
                (None, Some(sp)) => {
                    // Insertion in sequence: add a new node.
                    let seq_base = seq[sp];
                    let new_idx = self.nodes.len();
                    self.nodes.push(PoaNode::new(seq_base));

                    if let Some(prev) = prev_node {
                        self.add_edge(prev, new_idx, 1);
                    } else {
                        self.start_nodes.push(new_idx);
                    }
                    prev_node = Some(new_idx);
                }
                (Some(_gn), None) => {
                    // Deletion in sequence: graph node not consumed.
                    // We skip this node — no contribution from this sequence.
                }
                (None, None) => {
                    // Should not happen, but ignore.
                }
            }
        }
    }

    /// Extract consensus by finding the heaviest weighted path through the DAG.
    ///
    /// Uses topological-order DP: for each node, the score is the maximum
    /// incoming edge weight plus the predecessor's best score. The consensus
    /// is the path from the highest-scoring node traced back through best
    /// predecessors.
    pub fn consensus(&self) -> Vec<u8> {
        if self.nodes.is_empty() {
            return Vec::new();
        }

        let topo = self.topological_sort();

        // best_score[node] = maximum total weight of any path ending at node.
        let mut best_score = vec![0i64; self.nodes.len()];
        // best_pred[node] = predecessor on the heaviest path (None for start).
        let mut best_pred: Vec<Option<usize>> = vec![None; self.nodes.len()];

        for &node_idx in &topo {
            // For start nodes (no predecessors), score stays 0.
            for &pred in &self.nodes[node_idx].predecessors {
                // Find the weight of this edge in the predecessor's successor list.
                let weight = self.edge_weight(pred, node_idx);
                let candidate = best_score[pred] + weight as i64;

                if candidate > best_score[node_idx] {
                    best_score[node_idx] = candidate;
                    best_pred[node_idx] = Some(pred);
                }
            }
        }

        // Find the node with the maximum score.
        let mut best_end = 0;
        let mut max_score = i64::MIN;
        for &node_idx in &topo {
            if best_score[node_idx] > max_score {
                max_score = best_score[node_idx];
                best_end = node_idx;
            }
        }

        // Traceback: follow best_pred from best_end to a start node.
        let mut path = Vec::new();
        let mut current = best_end;
        path.push(self.nodes[current].base);

        while let Some(pred) = best_pred[current] {
            path.push(self.nodes[pred].base);
            current = pred;
        }

        path.reverse();
        path
    }

    /// Look up the weight of the edge from `from` to `to`.
    fn edge_weight(&self, from: usize, to: usize) -> usize {
        if let Some(pos) = self.nodes[from].successors.iter().position(|&s| s == to) {
            self.nodes[from].weights[pos]
        } else {
            0
        }
    }

    /// Topological sort of the graph using Kahn's algorithm.
    pub fn topological_sort(&self) -> Vec<usize> {
        let n = self.nodes.len();
        let mut in_degree = vec![0usize; n];

        for node in &self.nodes {
            for &succ in &node.successors {
                in_degree[succ] += 1;
            }
        }

        let mut queue: Vec<usize> = (0..n).filter(|&i| in_degree[i] == 0).collect();
        let mut order = Vec::with_capacity(n);

        while let Some(node_idx) = queue.pop() {
            order.push(node_idx);
            for &succ in &self.nodes[node_idx].successors {
                in_degree[succ] -= 1;
                if in_degree[succ] == 0 {
                    queue.push(succ);
                }
            }
        }

        order
    }

    /// Produce a DOT-format string for visualization of the graph.
    ///
    /// Nodes are labeled with their base character and index. Edge labels
    /// show the weight.
    pub fn to_dot(&self) -> String {
        let mut dot = String::from("digraph POA {\n");
        dot.push_str("    rankdir=LR;\n");

        for (i, node) in self.nodes.iter().enumerate() {
            let base = node.base as char;
            dot.push_str(&format!("    n{} [label=\"{}({})\"];\n", i, base, i));
        }

        for (i, node) in self.nodes.iter().enumerate() {
            for (edge_idx, &succ) in node.successors.iter().enumerate() {
                let w = node.weights[edge_idx];
                dot.push_str(&format!(
                    "    n{} -> n{} [label=\"{}\"];\n",
                    i, succ, w
                ));
            }
        }

        // Show aligned_to relationships as dashed edges.
        for (i, node) in self.nodes.iter().enumerate() {
            for &aligned in &node.aligned_to {
                if i < aligned {
                    dot.push_str(&format!(
                        "    n{} -> n{} [style=dashed, color=gray, dir=none];\n",
                        i, aligned
                    ));
                }
            }
        }

        dot.push_str("}\n");
        dot
    }

    /// Number of nodes in the graph.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of sequences that have been added to (or used to initialize)
    /// the graph.
    pub fn num_sequences(&self) -> usize {
        self.n_sequences
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_sequence_consensus() {
        let seq = b"ACGTACGT";
        let graph = PoaGraph::from_sequence(seq);
        let consensus = graph.consensus();
        assert_eq!(consensus, seq.to_vec());
    }

    #[test]
    fn two_identical_sequences() {
        let seq = b"ACGTACGT";
        let mut graph = PoaGraph::from_sequence(seq);
        let scoring = PoaScoring::default();
        graph.add_sequence(seq, &scoring).unwrap();
        let consensus = graph.consensus();
        assert_eq!(consensus, seq.to_vec());
        assert_eq!(graph.num_sequences(), 2);
    }

    #[test]
    fn two_sequences_single_snp() {
        // Original sequence is the reference. A second sequence has one SNP.
        // With only two sequences, the consensus should be the original because
        // the original path has the same weight as the alternate path, and
        // topological order favors the first-added nodes.
        let mut graph = PoaGraph::from_sequence(b"ACGTACGT");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"ACGAACGT", &scoring).unwrap();
        let consensus = graph.consensus();
        // Both paths have weight 1 through the SNP. The consensus should
        // be one of the two valid sequences.
        assert!(
            consensus == b"ACGTACGT" || consensus == b"ACGAACGT",
            "unexpected consensus: {:?}",
            String::from_utf8_lossy(&consensus)
        );
    }

    #[test]
    fn three_sequences_majority_rules() {
        // Two sequences say T, one says A at position 3.
        // Majority should win.
        let mut graph = PoaGraph::from_sequence(b"ACGTAC");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"ACGTAC", &scoring).unwrap();
        graph.add_sequence(b"ACGAAC", &scoring).unwrap();
        let consensus = graph.consensus();
        assert_eq!(
            consensus,
            b"ACGTAC".to_vec(),
            "consensus should follow majority: {:?}",
            String::from_utf8_lossy(&consensus)
        );
    }

    #[test]
    fn three_sequences_with_insertion() {
        // Two sequences agree on original, one has an extra base inserted.
        // Majority rules: consensus should follow the two without insertion.
        let mut graph = PoaGraph::from_sequence(b"ACGTAC");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"ACGTAC", &scoring).unwrap();
        graph.add_sequence(b"ACGXTAC", &scoring).unwrap();
        let consensus = graph.consensus();
        assert_eq!(
            consensus,
            b"ACGTAC".to_vec(),
            "consensus should follow majority without insertion: {:?}",
            String::from_utf8_lossy(&consensus)
        );
    }

    #[test]
    fn two_sequences_with_deletion() {
        // Second sequence has one base deleted.
        let mut graph = PoaGraph::from_sequence(b"ACGTAC");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"ACGAC", &scoring).unwrap();
        let consensus = graph.consensus();
        // The original path has all edges weight >= 1, and the shortcut
        // edge (skipping T) has weight 1. Original path wins because it
        // has more edges contributing total weight.
        assert_eq!(
            consensus,
            b"ACGTAC".to_vec(),
            "consensus should be original with deletion: {:?}",
            String::from_utf8_lossy(&consensus)
        );
    }

    #[test]
    fn empty_sequence_error() {
        let mut graph = PoaGraph::from_sequence(b"ACGT");
        let scoring = PoaScoring::default();
        let result = graph.add_sequence(b"", &scoring);
        assert!(result.is_err());
    }

    #[test]
    fn add_to_empty_graph_error() {
        let mut graph = PoaGraph::new();
        let scoring = PoaScoring::default();
        let result = graph.add_sequence(b"ACGT", &scoring);
        assert!(result.is_err());
    }

    #[test]
    fn from_sequence_creates_linear_graph() {
        let seq = b"ACGT";
        let graph = PoaGraph::from_sequence(seq);
        assert_eq!(graph.num_nodes(), 4);
        assert_eq!(graph.num_sequences(), 1);
        assert_eq!(graph.start_nodes, vec![0]);

        // Check linear chain structure.
        for i in 0..3 {
            assert_eq!(graph.nodes[i].successors, vec![i + 1]);
            assert_eq!(graph.nodes[i].weights, vec![1]);
        }
        assert!(graph.nodes[3].successors.is_empty());
        assert!(graph.nodes[0].predecessors.is_empty());
        for i in 1..4 {
            assert_eq!(graph.nodes[i].predecessors, vec![i - 1]);
        }
    }

    #[test]
    fn topological_sort_linear() {
        let graph = PoaGraph::from_sequence(b"ACGT");
        let topo = graph.topological_sort();
        // For a linear graph, topological order should visit all 4 nodes.
        assert_eq!(topo.len(), 4);
        // Each node should appear before its successor.
        for i in 0..3 {
            let pos_i = topo.iter().position(|&x| x == i).unwrap();
            let pos_next = topo.iter().position(|&x| x == i + 1).unwrap();
            assert!(pos_i < pos_next);
        }
    }

    #[test]
    fn to_dot_valid_syntax() {
        let mut graph = PoaGraph::from_sequence(b"AC");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"AT", &scoring).unwrap();

        let dot = graph.to_dot();
        assert!(dot.starts_with("digraph POA {"));
        assert!(dot.ends_with("}\n"));
        assert!(dot.contains("->"));
        assert!(dot.contains("rankdir=LR"));
        // Should have node labels.
        assert!(dot.contains("[label="));
    }

    #[test]
    fn empty_graph_consensus() {
        let graph = PoaGraph::new();
        let consensus = graph.consensus();
        assert!(consensus.is_empty());
    }

    #[test]
    fn single_base_sequences() {
        let mut graph = PoaGraph::from_sequence(b"A");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"A", &scoring).unwrap();
        let consensus = graph.consensus();
        assert_eq!(consensus, b"A".to_vec());
    }

    #[test]
    fn num_sequences_tracks_correctly() {
        let mut graph = PoaGraph::from_sequence(b"ACGT");
        assert_eq!(graph.num_sequences(), 1);
        let scoring = PoaScoring::default();
        graph.add_sequence(b"ACGT", &scoring).unwrap();
        assert_eq!(graph.num_sequences(), 2);
        graph.add_sequence(b"ACGT", &scoring).unwrap();
        assert_eq!(graph.num_sequences(), 3);
    }

    #[test]
    fn five_sequences_strong_consensus() {
        let mut graph = PoaGraph::from_sequence(b"ACGTACGT");
        let scoring = PoaScoring::default();
        // Four copies of the same sequence.
        for _ in 0..4 {
            graph.add_sequence(b"ACGTACGT", &scoring).unwrap();
        }
        let consensus = graph.consensus();
        assert_eq!(consensus, b"ACGTACGT".to_vec());
        assert_eq!(graph.num_sequences(), 5);
    }

    #[test]
    fn from_sequence_empty() {
        let graph = PoaGraph::from_sequence(b"");
        assert_eq!(graph.num_nodes(), 0);
        assert_eq!(graph.num_sequences(), 0);
        assert!(graph.consensus().is_empty());
    }

    #[test]
    fn custom_scoring() {
        let scoring = PoaScoring {
            match_score: 5,
            mismatch_score: -4,
            gap_score: -8,
        };
        let mut graph = PoaGraph::from_sequence(b"ACGT");
        graph.add_sequence(b"ACGT", &scoring).unwrap();
        let consensus = graph.consensus();
        assert_eq!(consensus, b"ACGT".to_vec());
    }

    #[test]
    fn three_way_majority_with_snp() {
        // Position 2: two sequences have G, one has A.
        let mut graph = PoaGraph::from_sequence(b"AAGAA");
        let scoring = PoaScoring::default();
        graph.add_sequence(b"AAGAA", &scoring).unwrap();
        graph.add_sequence(b"AAAAA", &scoring).unwrap();
        let consensus = graph.consensus();
        assert_eq!(
            consensus,
            b"AAGAA".to_vec(),
            "majority G should win: {:?}",
            String::from_utf8_lossy(&consensus)
        );
    }
}
