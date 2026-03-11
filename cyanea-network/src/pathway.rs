//! Biological pathway representation and topology-based analysis.
//!
//! Provides pathway data structures compatible with KEGG and Reactome,
//! topology-based pathway scoring, and pathway crosstalk analysis.

use crate::graph::{Graph, GraphType};
use cyanea_core::{CyaneaError, Result};
use std::collections::{HashMap, HashSet};

/// A biological pathway.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Pathway {
    /// Pathway identifier (e.g., "hsa04110" for KEGG, "R-HSA-69278" for Reactome).
    pub id: String,
    /// Pathway name.
    pub name: String,
    /// Source database.
    pub source: PathwaySource,
    /// Gene/protein members.
    pub members: Vec<String>,
    /// Category / top-level classification.
    pub category: Option<String>,
    /// Organism.
    pub organism: Option<String>,
}

/// Source database for pathway definitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PathwaySource {
    KEGG,
    Reactome,
    WikiPathways,
    BioCyc,
    Custom,
}

impl Pathway {
    pub fn new(id: &str, name: &str, source: PathwaySource, members: Vec<String>) -> Self {
        Self {
            id: id.to_string(),
            name: name.to_string(),
            source,
            members,
            category: None,
            organism: None,
        }
    }

    /// Number of members.
    pub fn size(&self) -> usize {
        self.members.len()
    }

    /// Check if a gene is in this pathway.
    pub fn contains(&self, gene: &str) -> bool {
        self.members.iter().any(|m| m == gene)
    }
}

/// A collection of pathways for enrichment and crosstalk analysis.
#[derive(Debug, Clone)]
pub struct PathwayDatabase {
    /// Pathways indexed by ID.
    pub pathways: HashMap<String, Pathway>,
}

impl PathwayDatabase {
    pub fn new() -> Self {
        Self {
            pathways: HashMap::new(),
        }
    }

    pub fn add(&mut self, pathway: Pathway) {
        self.pathways.insert(pathway.id.clone(), pathway);
    }

    /// Find all pathways containing a gene.
    pub fn pathways_for_gene(&self, gene: &str) -> Vec<&Pathway> {
        self.pathways
            .values()
            .filter(|p| p.contains(gene))
            .collect()
    }

    /// Number of pathways.
    pub fn len(&self) -> usize {
        self.pathways.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pathways.is_empty()
    }
}

/// Result of topology-based pathway scoring.
#[derive(Debug, Clone)]
pub struct PathwayScore {
    /// Pathway ID.
    pub pathway_id: String,
    /// Pathway name.
    pub pathway_name: String,
    /// Number of pathway members in the network.
    pub members_in_network: usize,
    /// Total pathway members.
    pub total_members: usize,
    /// Average centrality of pathway members in the network.
    pub avg_centrality: f64,
    /// Fraction of internal edges vs. total possible.
    pub internal_density: f64,
    /// Combined topology score.
    pub topology_score: f64,
}

/// Score pathways by network topology.
///
/// For each pathway, measures how densely its members are connected in
/// the network and how central they are. High scores indicate pathways
/// whose members form coherent network modules.
pub fn score_pathways_by_topology(
    network: &Graph,
    db: &PathwayDatabase,
) -> Result<Vec<PathwayScore>> {
    let network_nodes: HashSet<&str> = network.node_ids().into_iter().collect();
    let mut results = Vec::new();

    for pathway in db.pathways.values() {
        let members_in_net: Vec<&str> = pathway
            .members
            .iter()
            .filter(|m| network_nodes.contains(m.as_str()))
            .map(|m| m.as_str())
            .collect();

        let n_in = members_in_net.len();
        if n_in < 2 {
            results.push(PathwayScore {
                pathway_id: pathway.id.clone(),
                pathway_name: pathway.name.clone(),
                members_in_network: n_in,
                total_members: pathway.size(),
                avg_centrality: 0.0,
                internal_density: 0.0,
                topology_score: 0.0,
            });
            continue;
        }

        // Internal density
        let member_set: HashSet<&str> = members_in_net.iter().copied().collect();
        let mut internal_edges = 0;
        for &m in &members_in_net {
            if let Ok(neighbors) = network.neighbors(m) {
                for (n, _) in &neighbors {
                    if member_set.contains(n) {
                        internal_edges += 1;
                    }
                }
            }
        }

        let max_internal = match network.graph_type {
            GraphType::Undirected => n_in * (n_in - 1), // double-counted
            GraphType::Directed => n_in * (n_in - 1),
        };
        let density = if max_internal > 0 {
            internal_edges as f64 / max_internal as f64
        } else {
            0.0
        };

        // Average degree centrality of members
        let avg_degree: f64 = members_in_net
            .iter()
            .map(|m| network.degree(m).unwrap_or(0) as f64)
            .sum::<f64>()
            / n_in as f64;

        let max_degree = (network.node_count() - 1).max(1) as f64;
        let avg_centrality = avg_degree / max_degree;

        // Combined score
        let topology_score = 0.5 * density + 0.5 * avg_centrality;

        results.push(PathwayScore {
            pathway_id: pathway.id.clone(),
            pathway_name: pathway.name.clone(),
            members_in_network: n_in,
            total_members: pathway.size(),
            avg_centrality,
            internal_density: density,
            topology_score,
        });
    }

    results.sort_by(|a, b| {
        b.topology_score
            .partial_cmp(&a.topology_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(results)
}

/// Pathway crosstalk: overlap and connectivity between pathway pairs.
#[derive(Debug, Clone)]
pub struct CrosstalkResult {
    /// Pathway A ID.
    pub pathway_a: String,
    /// Pathway B ID.
    pub pathway_b: String,
    /// Number of shared members.
    pub shared_members: usize,
    /// Jaccard similarity of member sets.
    pub jaccard: f64,
    /// Number of cross-edges between pathways in the network.
    pub cross_edges: usize,
}

/// Compute pathway crosstalk for all pathway pairs.
///
/// Identifies pathways that share members or are connected in the network.
pub fn pathway_crosstalk(
    network: &Graph,
    db: &PathwayDatabase,
    min_overlap: usize,
) -> Result<Vec<CrosstalkResult>> {
    let pathway_ids: Vec<String> = db.pathways.keys().cloned().collect();
    let mut results = Vec::new();

    for i in 0..pathway_ids.len() {
        for j in (i + 1)..pathway_ids.len() {
            let pa = &db.pathways[&pathway_ids[i]];
            let pb = &db.pathways[&pathway_ids[j]];

            let set_a: HashSet<&str> = pa.members.iter().map(|m| m.as_str()).collect();
            let set_b: HashSet<&str> = pb.members.iter().map(|m| m.as_str()).collect();

            let shared = set_a.intersection(&set_b).count();
            if shared < min_overlap {
                continue;
            }

            let union = set_a.union(&set_b).count();
            let jaccard = if union > 0 {
                shared as f64 / union as f64
            } else {
                0.0
            };

            // Count cross-edges in network
            let network_nodes: HashSet<&str> = network.node_ids().into_iter().collect();
            let mut cross_edges = 0;
            for &m in &set_a {
                if !network_nodes.contains(m) {
                    continue;
                }
                if let Ok(neighbors) = network.neighbors(m) {
                    for (n, _) in &neighbors {
                        if set_b.contains(n) && !set_a.contains(n) {
                            cross_edges += 1;
                        }
                    }
                }
            }

            results.push(CrosstalkResult {
                pathway_a: pathway_ids[i].clone(),
                pathway_b: pathway_ids[j].clone(),
                shared_members: shared,
                jaccard,
                cross_edges,
            });
        }
    }

    results.sort_by(|a, b| {
        b.jaccard
            .partial_cmp(&a.jaccard)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(results)
}

/// Parse a simple GMT (Gene Matrix Transposed) format.
///
/// Each line: pathway_id\tpathway_name\tgene1\tgene2\t...
pub fn parse_gmt(content: &str) -> Result<PathwayDatabase> {
    let mut db = PathwayDatabase::new();

    for line in content.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let id = fields[0];
        let name = fields[1];
        let members: Vec<String> = fields[2..].iter().map(|s| s.to_string()).collect();

        db.add(Pathway::new(id, name, PathwaySource::Custom, members));
    }

    if db.is_empty() {
        return Err(CyaneaError::Parse("no pathways found in GMT".into()));
    }

    Ok(db)
}

/// Write pathways in GMT format.
pub fn write_gmt(db: &PathwayDatabase) -> String {
    let mut lines = Vec::new();
    let mut ids: Vec<&String> = db.pathways.keys().collect();
    ids.sort();

    for id in ids {
        let pw = &db.pathways[id];
        let mut parts = vec![pw.id.as_str(), pw.name.as_str()];
        for member in &pw.members {
            parts.push(member);
        }
        lines.push(parts.join("\t"));
    }

    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_pathways() -> PathwayDatabase {
        let mut db = PathwayDatabase::new();
        db.add(Pathway::new(
            "PW001",
            "Cell Cycle",
            PathwaySource::KEGG,
            vec!["CDK1".into(), "CDK2".into(), "CCNB1".into(), "TP53".into()],
        ));
        db.add(Pathway::new(
            "PW002",
            "Apoptosis",
            PathwaySource::KEGG,
            vec!["TP53".into(), "BAX".into(), "BCL2".into(), "CASP3".into()],
        ));
        db.add(Pathway::new(
            "PW003",
            "DNA Repair",
            PathwaySource::Reactome,
            vec!["BRCA1".into(), "TP53".into(), "ATM".into()],
        ));
        db
    }

    fn sample_ppi_graph() -> Graph {
        let mut g = Graph::new(GraphType::Undirected);
        for id in &[
            "CDK1", "CDK2", "CCNB1", "TP53", "BAX", "BCL2", "CASP3", "BRCA1", "ATM",
        ] {
            g.add_node(id, id).unwrap();
        }
        g.add_edge("CDK1", "CDK2", 0.9).unwrap();
        g.add_edge("CDK1", "CCNB1", 0.95).unwrap();
        g.add_edge("CDK2", "TP53", 0.8).unwrap();
        g.add_edge("TP53", "BAX", 0.85).unwrap();
        g.add_edge("TP53", "BCL2", 0.7).unwrap();
        g.add_edge("BAX", "BCL2", 0.9).unwrap();
        g.add_edge("BAX", "CASP3", 0.8).unwrap();
        g.add_edge("TP53", "BRCA1", 0.85).unwrap();
        g.add_edge("BRCA1", "ATM", 0.9).unwrap();
        g
    }

    #[test]
    fn test_pathway_basics() {
        let pw = Pathway::new(
            "PW001",
            "Test",
            PathwaySource::KEGG,
            vec!["A".into(), "B".into()],
        );
        assert_eq!(pw.size(), 2);
        assert!(pw.contains("A"));
        assert!(!pw.contains("C"));
    }

    #[test]
    fn test_pathway_database() {
        let db = sample_pathways();
        assert_eq!(db.len(), 3);
        let tp53_pathways = db.pathways_for_gene("TP53");
        assert_eq!(tp53_pathways.len(), 3); // In all three pathways
    }

    #[test]
    fn test_score_pathways() {
        let g = sample_ppi_graph();
        let db = sample_pathways();
        let scores = score_pathways_by_topology(&g, &db).unwrap();
        assert_eq!(scores.len(), 3);
        // All pathways should have positive topology scores
        for s in &scores {
            assert!(s.members_in_network > 0);
        }
    }

    #[test]
    fn test_pathway_crosstalk() {
        let g = sample_ppi_graph();
        let db = sample_pathways();
        let crosstalk = pathway_crosstalk(&g, &db, 1).unwrap();
        // PW001 and PW002 share TP53, PW002 and PW003 share TP53, PW001 and PW003 share TP53
        assert_eq!(crosstalk.len(), 3);
        for ct in &crosstalk {
            assert!(ct.shared_members >= 1);
        }
    }

    #[test]
    fn test_crosstalk_min_overlap() {
        let g = sample_ppi_graph();
        let db = sample_pathways();
        let crosstalk = pathway_crosstalk(&g, &db, 2).unwrap();
        // Only pairs with 2+ shared members
        for ct in &crosstalk {
            assert!(ct.shared_members >= 2);
        }
    }

    #[test]
    fn test_parse_gmt() {
        let gmt = "PW001\tCell Cycle\tCDK1\tCDK2\tCCNB1\n\
                   PW002\tApoptosis\tTP53\tBAX\tBCL2\n";
        let db = parse_gmt(gmt).unwrap();
        assert_eq!(db.len(), 2);
        assert_eq!(db.pathways["PW001"].members.len(), 3);
        assert_eq!(db.pathways["PW002"].name, "Apoptosis");
    }

    #[test]
    fn test_write_gmt() {
        let mut db = PathwayDatabase::new();
        db.add(Pathway::new(
            "PW001",
            "Test",
            PathwaySource::Custom,
            vec!["A".into(), "B".into()],
        ));
        let gmt = write_gmt(&db);
        assert!(gmt.contains("PW001\tTest\tA\tB"));
    }

    #[test]
    fn test_gmt_roundtrip() {
        let gmt_str = "PW001\tPathway One\tGENE1\tGENE2\tGENE3\n\
                       PW002\tPathway Two\tGENE4\tGENE5\n";
        let db = parse_gmt(gmt_str).unwrap();
        let output = write_gmt(&db);
        let db2 = parse_gmt(&output).unwrap();
        assert_eq!(db.len(), db2.len());
        assert_eq!(
            db.pathways["PW001"].members.len(),
            db2.pathways["PW001"].members.len()
        );
    }

    #[test]
    fn test_parse_gmt_empty() {
        assert!(parse_gmt("").is_err());
    }
}
