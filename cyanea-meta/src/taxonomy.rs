//! Taxonomic classification and management for metagenomics.
//!
//! Provides:
//! - [`TaxonomyDB`] — in-memory taxonomy tree with LCA support
//! - [`TaxonNode`] — taxonomy nodes with rank, name, parent links
//! - [`classify_sequence`] — k-mer based LCA classification (Kraken-style)
//! - [`get_lineage`] — full lineage from taxid to root
//! - [`lca`] — lowest common ancestor of multiple taxa

use std::collections::{BTreeMap, HashMap};
use crate::error::{MetaError, Result};

/// NCBI-style taxonomic rank.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum TaxonRank {
    /// Domain/Superkingdom
    Domain = 0,
    /// Phylum
    Phylum = 1,
    /// Class
    Class = 2,
    /// Order
    Order = 3,
    /// Family
    Family = 4,
    /// Genus
    Genus = 5,
    /// Species
    Species = 6,
    /// Strain or subspecies
    Strain = 7,
    /// Unranked
    Unranked = 8,
}

impl TaxonRank {
    /// Parse from string (e.g., "species" → Species).
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "domain" | "superkingdom" => TaxonRank::Domain,
            "phylum" => TaxonRank::Phylum,
            "class" => TaxonRank::Class,
            "order" => TaxonRank::Order,
            "family" => TaxonRank::Family,
            "genus" => TaxonRank::Genus,
            "species" => TaxonRank::Species,
            "strain" | "subspecies" => TaxonRank::Strain,
            _ => TaxonRank::Unranked,
        }
    }
}

/// A node in the taxonomy tree.
#[derive(Debug, Clone)]
pub struct TaxonNode {
    /// Unique taxon ID.
    pub taxid: u32,
    /// Parent taxon ID, or None for root.
    pub parent_id: Option<u32>,
    /// Taxonomic rank.
    pub rank: TaxonRank,
    /// Scientific name.
    pub name: String,
    /// Child taxon IDs.
    pub children: Vec<u32>,
}

/// In-memory taxonomy database (NCBI-style tree).
///
/// Maps taxon ID → node, supporting LCA queries and lineage retrieval.
#[derive(Debug, Clone)]
pub struct TaxonomyDB {
    nodes: HashMap<u32, TaxonNode>,
    #[allow(dead_code)]
    root_id: u32,
    /// k-mer hash → set of taxon IDs (for Kraken-style classification).
    kmer_db: BTreeMap<u64, Vec<u32>>,
    k: usize,
}

impl TaxonomyDB {
    /// Create a new empty taxonomy database with the given root taxon ID.
    pub fn new(root_id: u32) -> Self {
        let mut nodes = HashMap::new();
        nodes.insert(
            root_id,
            TaxonNode {
                taxid: root_id,
                parent_id: None,
                rank: TaxonRank::Domain,
                name: "root".to_string(),
                children: Vec::new(),
            },
        );
        Self {
            nodes,
            root_id,
            kmer_db: BTreeMap::new(),
            k: 31,
        }
    }

    /// Add a taxon node to the database.
    ///
    /// # Errors
    /// Returns an error if the parent taxon does not exist.
    pub fn add_node(&mut self, taxid: u32, parent_id: u32, rank: TaxonRank, name: &str) -> Result<()> {
        if !self.nodes.contains_key(&parent_id) {
            return Err(MetaError::Taxonomy(format!(
                "parent taxon {} not found",
                parent_id
            )));
        }

        let node = TaxonNode {
            taxid,
            parent_id: Some(parent_id),
            rank,
            name: name.to_string(),
            children: Vec::new(),
        };

        // Add to parent's children list
        if let Some(parent) = self.nodes.get_mut(&parent_id) {
            parent.children.push(taxid);
        }

        self.nodes.insert(taxid, node);
        Ok(())
    }

    /// Get the full lineage from a taxid to root (inclusive).
    ///
    /// # Errors
    /// Returns an error if the taxid does not exist.
    pub fn get_lineage(&self, taxid: u32) -> Result<Vec<u32>> {
        if !self.nodes.contains_key(&taxid) {
            return Err(MetaError::Taxonomy(format!("taxon {} not found", taxid)));
        }

        let mut lineage = Vec::new();
        let mut current = taxid;
        loop {
            lineage.push(current);
            match self.nodes.get(&current).and_then(|n| n.parent_id) {
                Some(parent) => current = parent,
                None => break,
            }
        }
        Ok(lineage)
    }

    /// Get the taxonomic rank for a taxon.
    ///
    /// # Errors
    /// Returns an error if the taxid does not exist.
    pub fn get_rank(&self, taxid: u32) -> Result<TaxonRank> {
        self.nodes
            .get(&taxid)
            .map(|n| n.rank)
            .ok_or_else(|| MetaError::Taxonomy(format!("taxon {} not found", taxid)))
    }

    /// Get the name for a taxon.
    ///
    /// # Errors
    /// Returns an error if the taxid does not exist.
    pub fn get_name(&self, taxid: u32) -> Result<String> {
        self.nodes
            .get(&taxid)
            .map(|n| n.name.clone())
            .ok_or_else(|| MetaError::Taxonomy(format!("taxon {} not found", taxid)))
    }

    /// Compute the lowest common ancestor (LCA) of multiple taxa.
    ///
    /// # Errors
    /// Returns an error if no taxa are provided or any taxid is invalid.
    pub fn lca(&self, taxa: &[u32]) -> Result<u32> {
        if taxa.is_empty() {
            return Err(MetaError::Taxonomy(
                "at least one taxon required for LCA".into(),
            ));
        }

        if taxa.len() == 1 {
            return Ok(taxa[0]);
        }

        // Get ancestor sets for all taxa
        let first_lineage = self.get_lineage(taxa[0])?;
        let mut common: std::collections::HashSet<u32> = first_lineage.iter().copied().collect();

        for &taxid in &taxa[1..] {
            let lineage = self.get_lineage(taxid)?;
            common.retain(|&t| lineage.contains(&t));
        }

        // Find the deepest (closest to leaf) common ancestor
        first_lineage
            .iter()
            .find(|&&t| common.contains(&t))
            .copied()
            .ok_or_else(|| MetaError::Taxonomy("no common ancestor found".into()))
    }

    /// Add a reference sequence for k-mer based classification.
    pub fn add_reference(&mut self, sequence: &[u8], taxon_id: u32) -> Result<()> {
        if !self.nodes.contains_key(&taxon_id) {
            return Err(MetaError::Taxonomy(format!("taxon {} not found", taxon_id)));
        }

        if sequence.len() < self.k {
            return Ok(());
        }

        let upper: Vec<u8> = sequence.iter().map(|b| b.to_ascii_uppercase()).collect();
        for window in upper.windows(self.k) {
            let hash = hash_kmer(window);
            let entry = self.kmer_db.entry(hash).or_default();
            if !entry.contains(&taxon_id) {
                entry.push(taxon_id);
            }
        }
        Ok(())
    }

    /// Classify a sequence using k-mer LCA method (Kraken-style).
    ///
    /// Returns the LCA of all k-mer matches, or None if no k-mers match.
    pub fn classify_sequence(&self, sequence: &[u8]) -> Result<Option<u32>> {
        if sequence.len() < self.k {
            return Ok(None);
        }

        let upper: Vec<u8> = sequence.iter().map(|b| b.to_ascii_uppercase()).collect();
        let mut hit_taxa: Vec<u32> = Vec::new();

        for window in upper.windows(self.k) {
            let hash = hash_kmer(window);
            if let Some(taxa) = self.kmer_db.get(&hash) {
                // Compute LCA of taxa for this k-mer
                if let Ok(kmer_lca) = self.lca(taxa) {
                    hit_taxa.push(kmer_lca);
                }
            }
        }

        if hit_taxa.is_empty() {
            return Ok(None);
        }

        // Final classification: LCA of all per-kmer LCAs
        self.lca(&hit_taxa).map(Some)
    }

    /// Get a node from the database.
    pub fn get_node(&self, taxid: u32) -> Option<&TaxonNode> {
        self.nodes.get(&taxid)
    }

    /// Set the k-mer size for classification.
    pub fn set_k(&mut self, k: usize) {
        self.k = k;
    }

    /// Number of taxa in the database.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Check if database is empty.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
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

    fn sample_db() -> TaxonomyDB {
        let mut db = TaxonomyDB::new(1);
        db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
        db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();
        db.add_node(202, 2, TaxonRank::Phylum, "Firmicutes").unwrap();
        db.add_node(2011, 201, TaxonRank::Species, "E. coli").unwrap();
        db.add_node(2012, 201, TaxonRank::Species, "Salmonella").unwrap();
        db.add_node(2021, 202, TaxonRank::Species, "B. subtilis").unwrap();
        db
    }

    #[test]
    fn lineage_to_root() {
        let db = sample_db();
        let lineage = db.get_lineage(2011).unwrap();
        assert_eq!(lineage, vec![2011, 201, 2, 1]);
    }

    #[test]
    fn lca_sibling_species() {
        let db = sample_db();
        // LCA(E. coli, Salmonella) = Proteobacteria
        let lca = db.lca(&[2011, 2012]).unwrap();
        assert_eq!(lca, 201);
    }

    #[test]
    fn lca_cross_phylum() {
        let db = sample_db();
        // LCA(E. coli, B. subtilis) = Bacteria
        let lca = db.lca(&[2011, 2021]).unwrap();
        assert_eq!(lca, 2);
    }

    #[test]
    fn lca_same_taxon() {
        let db = sample_db();
        let lca = db.lca(&[2011, 2011]).unwrap();
        assert_eq!(lca, 2011);
    }

    #[test]
    fn get_rank() {
        let db = sample_db();
        let rank = db.get_rank(2011).unwrap();
        assert_eq!(rank, TaxonRank::Species);
    }

    #[test]
    fn get_name() {
        let db = sample_db();
        let name = db.get_name(2011).unwrap();
        assert_eq!(name, "E. coli");
    }

    #[test]
    fn taxon_not_found() {
        let db = sample_db();
        assert!(db.get_rank(9999).is_err());
    }

    #[test]
    fn classify_exact_match() {
        let mut db = sample_db();
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        db.add_reference(seq, 2011).unwrap();
        let result = db.classify_sequence(seq).unwrap();
        assert_eq!(result, Some(2011));
    }

    #[test]
    fn classify_no_match() {
        let db = sample_db();
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let result = db.classify_sequence(seq).unwrap();
        assert_eq!(result, None);
    }
}
