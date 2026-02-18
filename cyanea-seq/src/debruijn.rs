//! De Bruijn graph construction and unitig extraction.
//!
//! Build a De Bruijn graph from sequences or raw k-mers, then extract
//! unitigs (maximal non-branching paths) with average coverage.

use std::collections::BTreeMap;

use cyanea_core::{CyaneaError, Result};

/// A node-centric De Bruijn graph built from k-mers.
///
/// Each k-mer is an edge from its (k-1)-prefix to its (k-1)-suffix.
/// Nodes are the (k-1)-mers; edges are the original k-mers.
#[derive(Debug, Clone)]
pub struct DeBruijnGraph {
    k: usize,
    /// Adjacency: prefix (k-1)-mer → list of suffix (k-1)-mers.
    edges: BTreeMap<Vec<u8>, Vec<Vec<u8>>>,
    /// Count of each k-mer seen during construction.
    kmer_counts: BTreeMap<Vec<u8>, usize>,
}

/// A unitig: a maximal non-branching path through the De Bruijn graph.
#[derive(Debug, Clone)]
pub struct Unitig {
    /// The assembled sequence of the unitig.
    pub sequence: Vec<u8>,
    /// Mean k-mer coverage along the unitig.
    pub coverage: f64,
}

impl DeBruijnGraph {
    /// Build a De Bruijn graph from a set of input sequences.
    ///
    /// Extracts all k-mers from each sequence and adds them as edges.
    /// Bases are uppercased; only A/C/G/T bases are accepted.
    ///
    /// # Errors
    ///
    /// Returns an error if `k < 2`, any sequence is shorter than `k`,
    /// or a sequence contains non-ACGT characters.
    pub fn from_sequences(sequences: &[&[u8]], k: usize) -> Result<Self> {
        if k < 2 {
            return Err(CyaneaError::InvalidInput(
                "k must be at least 2 for De Bruijn graph construction".into(),
            ));
        }
        let mut graph = Self {
            k,
            edges: BTreeMap::new(),
            kmer_counts: BTreeMap::new(),
        };
        for seq in sequences {
            if seq.len() < k {
                return Err(CyaneaError::InvalidInput(format!(
                    "sequence length {} is shorter than k={}",
                    seq.len(),
                    k
                )));
            }
            let upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            for b in &upper {
                if !matches!(b, b'A' | b'C' | b'G' | b'T') {
                    return Err(CyaneaError::InvalidInput(format!(
                        "invalid base '{}' for De Bruijn graph",
                        *b as char
                    )));
                }
            }
            for window in upper.windows(k) {
                graph.add_kmer(window);
            }
        }
        Ok(graph)
    }

    /// Build a De Bruijn graph from pre-extracted k-mers.
    ///
    /// All k-mers must have the same length (≥ 2).
    ///
    /// # Errors
    ///
    /// Returns an error if the k-mer slice is empty or k-mers differ in length.
    pub fn from_kmers(kmers: &[&[u8]]) -> Result<Self> {
        if kmers.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "at least one k-mer is required".into(),
            ));
        }
        let k = kmers[0].len();
        if k < 2 {
            return Err(CyaneaError::InvalidInput(
                "k-mers must have length at least 2".into(),
            ));
        }
        for kmer in kmers {
            if kmer.len() != k {
                return Err(CyaneaError::InvalidInput(format!(
                    "all k-mers must have the same length; expected {} but got {}",
                    k,
                    kmer.len()
                )));
            }
        }
        let mut graph = Self {
            k,
            edges: BTreeMap::new(),
            kmer_counts: BTreeMap::new(),
        };
        for kmer in kmers {
            let upper: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();
            graph.add_kmer(&upper);
        }
        Ok(graph)
    }

    fn add_kmer(&mut self, kmer: &[u8]) {
        let prefix = kmer[..self.k - 1].to_vec();
        let suffix = kmer[1..].to_vec();
        self.edges.entry(prefix).or_default().push(suffix.clone());
        // Ensure suffix node exists in the graph even if it has no outgoing edges.
        self.edges.entry(suffix).or_default();
        *self.kmer_counts.entry(kmer.to_vec()).or_insert(0) += 1;
    }

    /// Number of distinct (k-1)-mer nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.edges.len()
    }

    /// Number of k-mer edges in the graph (including duplicates).
    pub fn edge_count(&self) -> usize {
        self.kmer_counts.len()
    }

    /// Check whether a specific k-mer exists as an edge in the graph.
    pub fn contains_kmer(&self, kmer: &[u8]) -> bool {
        let upper: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();
        self.kmer_counts.contains_key(&upper)
    }

    /// Extract all unitigs (maximal non-branching paths).
    ///
    /// A unitig is a path where every internal node has exactly one incoming
    /// and one outgoing edge. Coverage is the mean count of constituent k-mers.
    pub fn unitigs(&self) -> Vec<Unitig> {
        // Compute in-degree for each node.
        let mut in_degree: BTreeMap<&Vec<u8>, usize> = BTreeMap::new();
        for (_, successors) in &self.edges {
            for s in successors {
                *in_degree.entry(s).or_insert(0) += 1;
            }
        }

        let out_degree = |node: &Vec<u8>| -> usize {
            self.edges.get(node).map_or(0, |v| v.len())
        };
        let in_deg = |node: &Vec<u8>| -> usize { in_degree.get(node).copied().unwrap_or(0) };

        // A node is a branching point if in-degree != 1 or out-degree != 1.
        let is_start = |node: &Vec<u8>| -> bool {
            in_deg(node) != 1 || out_degree(node) != 1
        };

        let mut visited: BTreeMap<(Vec<u8>, Vec<u8>), bool> = BTreeMap::new();
        let mut unitigs = Vec::new();

        // Start unitig walks from branching nodes.
        let nodes: Vec<Vec<u8>> = self.edges.keys().cloned().collect();
        for node in &nodes {
            if !is_start(node) {
                continue;
            }
            let successors = match self.edges.get(node) {
                Some(s) => s.clone(),
                None => continue,
            };
            for succ in &successors {
                let edge_key = (node.clone(), succ.clone());
                if visited.contains_key(&edge_key) {
                    continue;
                }
                visited.insert(edge_key, true);

                // Build the unitig sequence: start with the prefix node, then extend.
                let mut path_nodes: Vec<Vec<u8>> = vec![node.clone(), succ.clone()];
                let mut current = succ.clone();

                // Extend forward while the path is non-branching.
                while !is_start(&current) {
                    let next_successors = match self.edges.get(&current) {
                        Some(s) if s.len() == 1 => s,
                        _ => break,
                    };
                    let next = &next_successors[0];
                    let edge_key = (current.clone(), next.clone());
                    if visited.contains_key(&edge_key) {
                        break;
                    }
                    visited.insert(edge_key, true);
                    path_nodes.push(next.clone());
                    current = next.clone();
                }

                // Build sequence: first (k-1)-mer, then one base per subsequent node.
                let mut sequence = path_nodes[0].clone();
                for p in &path_nodes[1..] {
                    sequence.push(*p.last().unwrap());
                }

                // Compute mean coverage from constituent k-mers.
                let k = self.k;
                let mut total_count = 0usize;
                let mut n_kmers = 0usize;
                for window in sequence.windows(k) {
                    if let Some(&c) = self.kmer_counts.get(window) {
                        total_count += c;
                        n_kmers += 1;
                    }
                }
                let coverage = if n_kmers > 0 {
                    total_count as f64 / n_kmers as f64
                } else {
                    0.0
                };

                unitigs.push(Unitig { sequence, coverage });
            }
        }

        // Handle isolated cycles (all nodes have in=1, out=1).
        // Find any unvisited edges.
        for (node, successors) in &self.edges {
            for succ in successors {
                let edge_key = (node.clone(), succ.clone());
                if visited.contains_key(&edge_key) {
                    continue;
                }
                visited.insert(edge_key.clone(), true);

                let mut path_nodes: Vec<Vec<u8>> = vec![node.clone(), succ.clone()];
                let mut current = succ.clone();

                loop {
                    let next_successors = match self.edges.get(&current) {
                        Some(s) if s.len() == 1 => s,
                        _ => break,
                    };
                    let next = &next_successors[0];
                    let ek = (current.clone(), next.clone());
                    if visited.contains_key(&ek) {
                        break;
                    }
                    visited.insert(ek, true);
                    path_nodes.push(next.clone());
                    current = next.clone();
                }

                let mut sequence = path_nodes[0].clone();
                for p in &path_nodes[1..] {
                    sequence.push(*p.last().unwrap());
                }

                let k = self.k;
                let mut total_count = 0usize;
                let mut n_kmers = 0usize;
                for window in sequence.windows(k) {
                    if let Some(&c) = self.kmer_counts.get(window) {
                        total_count += c;
                        n_kmers += 1;
                    }
                }
                let coverage = if n_kmers > 0 {
                    total_count as f64 / n_kmers as f64
                } else {
                    0.0
                };

                unitigs.push(Unitig { sequence, coverage });
            }
        }

        unitigs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn debruijn_simple_sequence() {
        // ACGTACGT with k=3 should build a valid graph.
        let graph = DeBruijnGraph::from_sequences(&[b"ACGTACGT"], 3).unwrap();
        assert!(graph.node_count() > 0);
        assert!(graph.edge_count() > 0);
    }

    #[test]
    fn debruijn_node_edge_counts() {
        // ACGT with k=3: k-mers are ACG, CGT → 2 edges.
        // Nodes (2-mers): AC, CG, GT → 3 nodes.
        let graph = DeBruijnGraph::from_sequences(&[b"ACGT"], 3).unwrap();
        assert_eq!(graph.edge_count(), 2);
        assert_eq!(graph.node_count(), 3);
    }

    #[test]
    fn debruijn_contains_kmer() {
        let graph = DeBruijnGraph::from_sequences(&[b"ACGTACGT"], 3).unwrap();
        assert!(graph.contains_kmer(b"ACG"));
        assert!(graph.contains_kmer(b"CGT"));
        assert!(!graph.contains_kmer(b"AAA"));
    }

    #[test]
    fn unitig_extraction_simple() {
        // Linear sequence with no branching should yield unitig(s).
        let graph = DeBruijnGraph::from_sequences(&[b"ACGTACGT"], 3).unwrap();
        let unitigs = graph.unitigs();
        assert!(!unitigs.is_empty());
        // At least one unitig should contain our original k-mers.
        let total_len: usize = unitigs.iter().map(|u| u.sequence.len()).sum();
        assert!(total_len >= 3); // at least one k-mer length
    }

    #[test]
    fn unitig_coverage_correct() {
        // Feed the same sequence twice — coverage should be 2.0.
        let graph =
            DeBruijnGraph::from_sequences(&[b"ACGT", b"ACGT"], 3).unwrap();
        let unitigs = graph.unitigs();
        assert!(!unitigs.is_empty());
        for u in &unitigs {
            assert!((u.coverage - 2.0).abs() < 1e-10);
        }
    }
}
