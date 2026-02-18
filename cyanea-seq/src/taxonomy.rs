//! Taxonomic classification — taxonomy trees, LCA queries, k-mer classifiers.
//!
//! Build a taxonomy tree, compute lowest common ancestors, and classify
//! sequences using a Kraken-style k-mer approach.

use std::collections::BTreeMap;

use cyanea_core::{CyaneaError, Result};

/// Taxonomic rank in the NCBI hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TaxonRank {
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Unranked,
}

/// A node in the taxonomy tree.
#[derive(Debug, Clone)]
pub struct TaxonomyNode {
    /// Unique node identifier (index in the tree's node vector).
    pub id: usize,
    /// Taxon name.
    pub name: String,
    /// Taxonomic rank.
    pub rank: TaxonRank,
    /// Parent node id, or `None` for the root.
    pub parent: Option<usize>,
}

/// A rooted taxonomy tree.
///
/// Nodes are stored in a flat vector; parent links encode the tree structure.
#[derive(Debug, Clone)]
pub struct TaxonomyTree {
    nodes: Vec<TaxonomyNode>,
}

impl TaxonomyTree {
    /// Create a new empty taxonomy tree.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    /// Add a node to the tree. Returns the node's id.
    ///
    /// The root node should have `parent = None`. All other nodes must
    /// reference a valid parent id.
    pub fn add_node(&mut self, name: &str, rank: TaxonRank, parent: Option<usize>) -> usize {
        let id = self.nodes.len();
        self.nodes.push(TaxonomyNode {
            id,
            name: name.to_string(),
            rank,
            parent,
        });
        id
    }

    /// Compute the lowest common ancestor of a set of node ids.
    ///
    /// Returns the deepest node that is an ancestor of all input nodes.
    ///
    /// # Errors
    ///
    /// Returns an error if `ids` is empty or any id is out of range.
    pub fn lca(&self, ids: &[usize]) -> Result<usize> {
        if ids.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "at least one taxon id is required for LCA".into(),
            ));
        }
        for &id in ids {
            if id >= self.nodes.len() {
                return Err(CyaneaError::InvalidInput(format!(
                    "taxon id {} is out of range (tree has {} nodes)",
                    id,
                    self.nodes.len()
                )));
            }
        }

        // Get ancestor set for the first node.
        let mut common_ancestors = self.ancestor_set(ids[0]);

        // Intersect with ancestor sets of remaining nodes.
        for &id in &ids[1..] {
            let ancestors = self.ancestor_set(id);
            common_ancestors.retain(|a| ancestors.contains(a));
        }

        // The LCA is the deepest (greatest depth) common ancestor.
        common_ancestors
            .into_iter()
            .max_by_key(|&a| self.depth(a))
            .ok_or_else(|| {
                CyaneaError::InvalidInput("no common ancestor found".into())
            })
    }

    /// Get the full lineage (path to root) for a node, ordered from the node to root.
    pub fn lineage(&self, id: usize) -> Vec<usize> {
        let mut path = Vec::new();
        let mut current = id;
        loop {
            if current >= self.nodes.len() {
                break;
            }
            path.push(current);
            match self.nodes[current].parent {
                Some(p) => current = p,
                None => break,
            }
        }
        path
    }

    /// Depth of a node (root = 0).
    pub fn depth(&self, id: usize) -> usize {
        self.lineage(id).len().saturating_sub(1)
    }

    /// Get the set of all ancestors of a node (including itself).
    fn ancestor_set(&self, id: usize) -> Vec<usize> {
        self.lineage(id)
    }
}

impl Default for TaxonomyTree {
    fn default() -> Self {
        Self::new()
    }
}

/// K-mer based taxonomic classifier (Kraken-style).
///
/// Maps k-mer hashes to sets of taxon ids. Classification of a query
/// sequence takes the LCA of all k-mer hits.
#[derive(Debug, Clone)]
pub struct KmerClassifier {
    /// k-mer hash → taxon ids that contain this k-mer.
    db: BTreeMap<u64, Vec<usize>>,
    /// The taxonomy tree for LCA computation.
    taxonomy: TaxonomyTree,
    /// k-mer length.
    k: usize,
}

impl KmerClassifier {
    /// Create a new classifier with the given taxonomy and k-mer size.
    pub fn new(taxonomy: TaxonomyTree, k: usize) -> Self {
        Self {
            db: BTreeMap::new(),
            taxonomy,
            k,
        }
    }

    /// Index a reference sequence under the given taxon id.
    ///
    /// Extracts all k-mers, hashes them, and associates them with `taxon_id`.
    pub fn add_reference(&mut self, sequence: &[u8], taxon_id: usize) {
        if sequence.len() < self.k {
            return;
        }
        let upper: Vec<u8> = sequence.iter().map(|b| b.to_ascii_uppercase()).collect();
        for window in upper.windows(self.k) {
            let hash = hash_kmer(window);
            let entry = self.db.entry(hash).or_default();
            if !entry.contains(&taxon_id) {
                entry.push(taxon_id);
            }
        }
    }

    /// Classify a query sequence by k-mer matching.
    ///
    /// For each k-mer in the query, looks up matching taxa. Returns the LCA
    /// of all matched taxa, preferring deeper (more specific) assignments.
    /// Returns `None` if no k-mers match.
    pub fn classify(&self, sequence: &[u8]) -> Option<usize> {
        if sequence.len() < self.k {
            return None;
        }

        let upper: Vec<u8> = sequence.iter().map(|b| b.to_ascii_uppercase()).collect();
        let mut hit_taxa: Vec<usize> = Vec::new();

        for window in upper.windows(self.k) {
            let hash = hash_kmer(window);
            if let Some(taxa) = self.db.get(&hash) {
                // For each k-mer, compute LCA of its taxa and use that.
                if let Ok(kmer_lca) = self.taxonomy.lca(taxa) {
                    hit_taxa.push(kmer_lca);
                }
            }
        }

        if hit_taxa.is_empty() {
            return None;
        }

        // Final classification: LCA of all per-kmer assignments.
        self.taxonomy.lca(&hit_taxa).ok()
    }
}

/// Simple 2-bit hash for DNA k-mers.
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut h: u64 = 0;
    for &b in kmer {
        h = h.wrapping_mul(4);
        h = h.wrapping_add(match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        });
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_tree() -> TaxonomyTree {
        let mut tree = TaxonomyTree::new();
        // 0: root
        tree.add_node("root", TaxonRank::Unranked, None);
        // 1: Bacteria
        tree.add_node("Bacteria", TaxonRank::Domain, Some(0));
        // 2: Proteobacteria
        tree.add_node("Proteobacteria", TaxonRank::Phylum, Some(1));
        // 3: Firmicutes
        tree.add_node("Firmicutes", TaxonRank::Phylum, Some(1));
        // 4: E. coli (under Proteobacteria)
        tree.add_node("E. coli", TaxonRank::Species, Some(2));
        // 5: Salmonella (under Proteobacteria)
        tree.add_node("Salmonella", TaxonRank::Species, Some(2));
        // 6: B. subtilis (under Firmicutes)
        tree.add_node("B. subtilis", TaxonRank::Species, Some(3));
        tree
    }

    #[test]
    fn taxonomy_tree_construction() {
        let tree = sample_tree();
        assert_eq!(tree.nodes.len(), 7);
        assert_eq!(tree.depth(0), 0); // root
        assert_eq!(tree.depth(4), 3); // E. coli: root → Bacteria → Proteo → E. coli
    }

    #[test]
    fn lca_sibling_species() {
        let tree = sample_tree();
        // LCA(E. coli, Salmonella) = Proteobacteria
        let lca = tree.lca(&[4, 5]).unwrap();
        assert_eq!(lca, 2);
    }

    #[test]
    fn lca_same_node() {
        let tree = sample_tree();
        let lca = tree.lca(&[4, 4]).unwrap();
        assert_eq!(lca, 4);
    }

    #[test]
    fn classify_exact_match() {
        let tree = sample_tree();
        let mut classifier = KmerClassifier::new(tree, 4);
        // Add a reference for E. coli (taxon 4).
        classifier.add_reference(b"ACGTACGTACGT", 4);
        // Query with the same sequence should classify to E. coli.
        let result = classifier.classify(b"ACGTACGTACGT");
        assert_eq!(result, Some(4));
    }

    #[test]
    fn classify_ambiguous_lca() {
        let tree = sample_tree();
        let mut classifier = KmerClassifier::new(tree, 4);
        // Add overlapping references for two species.
        classifier.add_reference(b"ACGTACGT", 4); // E. coli
        classifier.add_reference(b"ACGTACGT", 5); // Salmonella
        // Query with shared sequence should resolve to LCA = Proteobacteria.
        let result = classifier.classify(b"ACGTACGT");
        assert_eq!(result, Some(2));
    }
}
