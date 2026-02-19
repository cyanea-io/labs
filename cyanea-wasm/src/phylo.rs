//! Phylogenetics WASM bindings: Newick parsing, evolutionary distances,
//! tree construction (UPGMA, NJ), tree comparison (Robinson-Foulds),
//! NEXUS I/O, and simulation (sequence evolution, coalescent trees).

use serde::{Deserialize, Serialize};

use cyanea_ml::DistanceMatrix;
use cyanea_phylo::compare::{robinson_foulds, robinson_foulds_normalized};
use cyanea_phylo::distance::{jukes_cantor, kimura_2p, p_distance};
use cyanea_phylo::newick;
use cyanea_phylo::nexus;
use cyanea_phylo::simulation;
use cyanea_phylo::subst_model::{Hky85Model, Jc69Model};
use cyanea_phylo::tree::PhyloTree;
use cyanea_phylo::construct::{neighbor_joining, upgma};

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// Serializable tree info.
#[derive(Debug, Serialize)]
pub struct JsTreeInfo {
    pub leaf_count: usize,
    pub internal_count: usize,
    pub total_nodes: usize,
    pub leaf_names: Vec<String>,
    pub newick: String,
}

/// Serializable Robinson-Foulds distance.
#[derive(Debug, Serialize)]
pub struct JsRFDistance {
    pub distance: usize,
    pub normalized: f64,
}

/// Serializable NEXUS file representation.
#[derive(Debug, Serialize)]
pub struct JsNexusFile {
    pub taxa: Vec<String>,
    pub trees: Vec<JsNamedTree>,
}

/// Serializable named tree from a NEXUS TREES block.
#[derive(Debug, Serialize)]
pub struct JsNamedTree {
    pub name: String,
    pub newick: String,
    pub leaf_count: usize,
}

/// Serializable simulated alignment.
#[derive(Debug, Serialize)]
pub struct JsSimulatedAlignment {
    pub names: Vec<String>,
    pub sequences: Vec<String>,
    pub n_substitutions: usize,
}

/// Serializable coalescent tree result.
#[derive(Debug, Serialize)]
pub struct JsCoalescentTree {
    pub newick: String,
    pub n_leaves: usize,
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Parse a Newick string and return tree info as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn newick_info(newick: &str) -> String {
    let tree = match PhyloTree::from_newick(newick) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };
    let leaf_count = tree.leaf_count();
    let total = tree.node_count();
    let js = JsTreeInfo {
        leaf_count,
        internal_count: total - leaf_count,
        total_nodes: total,
        leaf_names: tree.leaf_names(),
        newick: tree.to_newick(),
    };
    wasm_ok(&js)
}

/// Compute evolutionary distance between two sequences.
///
/// `model` is one of `"p"`, `"jc"`, or `"k2p"`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn evolutionary_distance(seq1: &str, seq2: &str, model: &str) -> String {
    match model.to_ascii_lowercase().as_str() {
        "p" => match p_distance(seq1.as_bytes(), seq2.as_bytes()) {
            Ok(d) => wasm_ok(&d),
            Err(e) => wasm_err(e),
        },
        "jc" => {
            let p = match p_distance(seq1.as_bytes(), seq2.as_bytes()) {
                Ok(d) => d,
                Err(e) => return wasm_err(e),
            };
            match jukes_cantor(p) {
                Ok(d) => wasm_ok(&d),
                Err(e) => wasm_err(e),
            }
        }
        "k2p" => {
            // For K2P we need transition/transversion proportions.
            // Compute them from the sequences.
            let a = seq1.as_bytes();
            let b = seq2.as_bytes();
            if a.len() != b.len() || a.is_empty() {
                return wasm_err("sequences must be non-empty and equal length");
            }
            let total = a.len() as f64;
            let mut transitions = 0usize;
            let mut transversions = 0usize;
            for (&x, &y) in a.iter().zip(b.iter()) {
                let xu = x.to_ascii_uppercase();
                let yu = y.to_ascii_uppercase();
                if xu == yu {
                    continue;
                }
                match (xu, yu) {
                    (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => {
                        transitions += 1;
                    }
                    (b'A', b'C')
                    | (b'C', b'A')
                    | (b'A', b'T')
                    | (b'T', b'A')
                    | (b'G', b'C')
                    | (b'C', b'G')
                    | (b'G', b'T')
                    | (b'T', b'G') => {
                        transversions += 1;
                    }
                    _ => {}
                }
            }
            match kimura_2p(transitions as f64 / total, transversions as f64 / total) {
                Ok(d) => wasm_ok(&d),
                Err(e) => wasm_err(e),
            }
        }
        _ => wasm_err(format!(
            "unknown distance model: {model:?} (expected \"p\", \"jc\", or \"k2p\")"
        )),
    }
}

/// Build a UPGMA tree from a distance matrix.
///
/// `labels_json`: JSON array of strings. `matrix_json`: JSON 2D array of f64.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn build_upgma(labels_json: &str, matrix_json: &str) -> String {
    let labels: Vec<String> = match serde_json::from_str(labels_json) {
        Ok(l) => l,
        Err(e) => return wasm_err(format!("invalid labels JSON: {e}")),
    };
    let matrix: Vec<Vec<f64>> = match serde_json::from_str(matrix_json) {
        Ok(m) => m,
        Err(e) => return wasm_err(format!("invalid matrix JSON: {e}")),
    };
    let dm = match matrix_to_distance_matrix(&labels, &matrix) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match upgma(&dm, &labels) {
        Ok(tree) => wasm_ok(&tree.to_newick()),
        Err(e) => wasm_err(e),
    }
}

/// Build a Neighbor-Joining tree from a distance matrix.
///
/// `labels_json`: JSON array of strings. `matrix_json`: JSON 2D array of f64.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn build_nj(labels_json: &str, matrix_json: &str) -> String {
    let labels: Vec<String> = match serde_json::from_str(labels_json) {
        Ok(l) => l,
        Err(e) => return wasm_err(format!("invalid labels JSON: {e}")),
    };
    let matrix: Vec<Vec<f64>> = match serde_json::from_str(matrix_json) {
        Ok(m) => m,
        Err(e) => return wasm_err(format!("invalid matrix JSON: {e}")),
    };
    let dm = match matrix_to_distance_matrix(&labels, &matrix) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match neighbor_joining(&dm, &labels) {
        Ok(tree) => wasm_ok(&tree.to_newick()),
        Err(e) => wasm_err(e),
    }
}

/// Robinson-Foulds distance between two Newick trees.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn rf_distance(newick1: &str, newick2: &str) -> String {
    let t1 = match PhyloTree::from_newick(newick1) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };
    let t2 = match PhyloTree::from_newick(newick2) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };
    let dist = match robinson_foulds(&t1, &t2) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let norm = match robinson_foulds_normalized(&t1, &t2) {
        Ok(n) => n,
        Err(e) => return wasm_err(e),
    };
    let js = JsRFDistance {
        distance: dist,
        normalized: norm,
    };
    wasm_ok(&js)
}

/// Parse a NEXUS format string and return file contents as JSON.
///
/// Returns `JsNexusFile` with taxa list and named trees.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_nexus(text: &str) -> String {
    let nf = match nexus::parse(text) {
        Ok(n) => n,
        Err(e) => return wasm_err(e),
    };
    let trees: Vec<JsNamedTree> = nf
        .trees
        .iter()
        .map(|nt| JsNamedTree {
            name: nt.name.clone(),
            newick: nt.tree.to_newick(),
            leaf_count: nt.tree.leaf_count(),
        })
        .collect();
    let js = JsNexusFile {
        taxa: nf.taxa,
        trees,
    };
    wasm_ok(&js)
}

/// Write a NEXUS file from taxa and trees JSON.
///
/// `taxa_json`: JSON array of taxon name strings.
/// `trees_json`: JSON array of `{"name": "...", "newick": "..."}` objects.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn write_nexus(taxa_json: &str, trees_json: &str) -> String {
    let taxa: Vec<String> = match serde_json::from_str(taxa_json) {
        Ok(t) => t,
        Err(e) => return wasm_err(format!("invalid taxa JSON: {e}")),
    };

    #[derive(Deserialize)]
    struct TreeInput {
        name: String,
        newick: String,
    }

    let tree_inputs: Vec<TreeInput> = match serde_json::from_str(trees_json) {
        Ok(t) => t,
        Err(e) => return wasm_err(format!("invalid trees JSON: {e}")),
    };

    let mut parsed_trees: Vec<(String, PhyloTree)> = Vec::with_capacity(tree_inputs.len());
    for ti in &tree_inputs {
        let tree = match newick::parse(&ti.newick) {
            Ok(t) => t,
            Err(e) => return wasm_err(format!("failed to parse newick for tree '{}': {e}", ti.name)),
        };
        parsed_trees.push((ti.name.clone(), tree));
    }

    let pairs: Vec<(&str, &PhyloTree)> = parsed_trees
        .iter()
        .map(|(name, tree)| (name.as_str(), tree))
        .collect();

    let output = nexus::write(&taxa, &pairs);
    wasm_ok(&output)
}

/// Simulate sequence evolution along a phylogenetic tree.
///
/// `newick`: Newick string for the guide tree.
/// `seq_length`: length of sequences to simulate.
/// `model`: substitution model — `"jc69"` or `"k2p"` (defaults to JC69 for unknown).
/// `seed`: random seed for reproducibility.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn simulate_evolution(newick_str: &str, seq_length: usize, model: &str, seed: u64) -> String {
    let tree = match PhyloTree::from_newick(newick_str) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };

    // Build the substitution model based on the model name.
    // K2P is equivalent to HKY85 with equal base frequencies and kappa=2.
    let jc69 = Jc69Model::new();
    let k2p_result = Hky85Model::new(2.0, [0.25; 4]);
    let model_obj: &dyn cyanea_phylo::subst_model::SubstitutionModel =
        match model.to_ascii_lowercase().as_str() {
            "k2p" => match k2p_result {
                Ok(ref m) => m,
                Err(e) => return wasm_err(e),
            },
            // "jc69" and any unknown model default to JC69
            _ => &jc69,
        };

    let aln = match simulation::simulate_evolution(&tree, model_obj, seq_length, seed) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };

    let js = JsSimulatedAlignment {
        names: aln.names,
        sequences: aln
            .sequences
            .iter()
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect(),
        n_substitutions: aln.n_substitutions,
    };
    wasm_ok(&js)
}

/// Simulate a coalescent (Kingman) tree for a constant-size population.
///
/// `n_samples`: number of sampled lineages (>= 2).
/// `pop_size`: effective population size (> 0).
/// `seed`: random seed for reproducibility.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn simulate_coalescent(n_samples: usize, pop_size: f64, seed: u64) -> String {
    let tree = match simulation::simulate_coalescent(n_samples, pop_size, seed) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };
    let js = JsCoalescentTree {
        newick: tree.to_newick(),
        n_leaves: tree.leaf_count(),
    };
    wasm_ok(&js)
}

/// Simulate a coalescent tree with exponential population growth.
///
/// `n_samples`: number of sampled lineages (>= 2).
/// `pop_size`: current effective population size (> 0).
/// `growth_rate`: exponential growth rate (>= 0).
/// `seed`: random seed for reproducibility.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn simulate_coalescent_growth(
    n_samples: usize,
    pop_size: f64,
    growth_rate: f64,
    seed: u64,
) -> String {
    let tree = match simulation::simulate_coalescent_growth(n_samples, pop_size, growth_rate, seed) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };
    let js = JsCoalescentTree {
        newick: tree.to_newick(),
        n_leaves: tree.leaf_count(),
    };
    wasm_ok(&js)
}

// ── Helpers ──────────────────────────────────────────────────────────────

/// Convert a full (symmetric) distance matrix to a condensed DistanceMatrix.
fn matrix_to_distance_matrix(
    labels: &[String],
    matrix: &[Vec<f64>],
) -> Result<DistanceMatrix, String> {
    let n = labels.len();
    if n < 2 {
        return Err("need at least 2 taxa".into());
    }
    if matrix.len() != n {
        return Err(format!(
            "matrix rows ({}) don't match label count ({})",
            matrix.len(),
            n
        ));
    }
    for (i, row) in matrix.iter().enumerate() {
        if row.len() != n {
            return Err(format!(
                "matrix row {} has {} columns, expected {}",
                i,
                row.len(),
                n
            ));
        }
    }
    // Extract upper triangle into condensed form
    let mut condensed = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            condensed.push(matrix[i][j]);
        }
    }
    DistanceMatrix::from_condensed(condensed, n).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newick_info_basic() {
        let json = newick_info("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let info = &v["ok"];
        assert_eq!(info["leaf_count"], 4);
        assert_eq!(info["total_nodes"], 7);
        assert_eq!(info["internal_count"], 3);
        let names = info["leaf_names"].as_array().unwrap();
        assert_eq!(names.len(), 4);
        assert!(info["newick"].as_str().unwrap().contains(';'));
    }

    #[test]
    fn newick_info_invalid() {
        let json = newick_info("not a newick string");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn evolutionary_distance_p() {
        let json = evolutionary_distance("ACGT", "ACGT", "p");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap()).abs() < 1e-10);
    }

    #[test]
    fn evolutionary_distance_jc() {
        let json = evolutionary_distance("AAAAAAAAAA", "AAAAAAAAAT", "jc");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn evolutionary_distance_k2p() {
        let json = evolutionary_distance("AAAAAAAAAA", "AAAAAAAAGT", "k2p");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn evolutionary_distance_invalid_model() {
        let json = evolutionary_distance("ACGT", "ACGT", "foobar");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].as_str().unwrap().contains("unknown distance model"));
    }

    #[test]
    fn build_upgma_basic() {
        let labels = r#"["A","B","C"]"#;
        let matrix = r#"[[0,2,4],[2,0,4],[4,4,0]]"#;
        let json = build_upgma(labels, matrix);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let newick = v["ok"].as_str().unwrap();
        assert!(newick.contains(';'));
        assert!(newick.contains('A'));
        assert!(newick.contains('B'));
        assert!(newick.contains('C'));
    }

    #[test]
    fn build_nj_basic() {
        let labels = r#"["A","B","C"]"#;
        let matrix = r#"[[0,5,9],[5,0,10],[9,10,0]]"#;
        let json = build_nj(labels, matrix);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let newick = v["ok"].as_str().unwrap();
        assert!(newick.contains(';'));
        assert!(newick.contains('A'));
    }

    #[test]
    fn build_upgma_invalid_json() {
        let json = build_upgma("not json", "[[0]]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn rf_distance_identical() {
        let json = rf_distance("((A,B),(C,D));", "((A,B),(C,D));");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let rf = &v["ok"];
        assert_eq!(rf["distance"], 0);
        assert!((rf["normalized"].as_f64().unwrap()).abs() < 1e-10);
    }

    #[test]
    fn rf_distance_different() {
        let json = rf_distance("((A,B),(C,D));", "((A,C),(B,D));");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let rf = &v["ok"];
        assert!(rf["distance"].as_u64().unwrap() > 0);
        assert!(rf["normalized"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn rf_distance_different_leaves_error() {
        let json = rf_distance("((A,B),(C,D));", "((A,B),(C,E));");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    // ── NEXUS tests ──────────────────────────────────────────────────────

    #[test]
    fn parse_nexus_basic() {
        let input = "#NEXUS\nbegin taxa;\n dimensions ntax=3;\n taxlabels A B C;\nend;\nbegin trees;\n tree t1 = (A,(B,C));\nend;\n";
        let json = parse_nexus(input);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let nf = &v["ok"];
        assert_eq!(nf["trees"].as_array().unwrap().len(), 1);
        assert!(nf["trees"][0]["newick"].as_str().unwrap().contains(';'));
        assert!(nf["trees"][0]["leaf_count"].as_u64().unwrap() >= 3);
    }

    #[test]
    fn parse_nexus_taxa() {
        let input = "#NEXUS\nbegin taxa;\n dimensions ntax=3;\n taxlabels A B C;\nend;\nbegin trees;\n tree t1 = (A,(B,C));\nend;\n";
        let json = parse_nexus(input);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let taxa = v["ok"]["taxa"].as_array().unwrap();
        assert_eq!(taxa.len(), 3);
        let taxa_strs: Vec<&str> = taxa.iter().map(|t| t.as_str().unwrap()).collect();
        assert!(taxa_strs.contains(&"A"));
        assert!(taxa_strs.contains(&"B"));
        assert!(taxa_strs.contains(&"C"));
    }

    #[test]
    fn write_nexus_basic() {
        let taxa_json = r#"["A","B","C"]"#;
        let trees_json = r#"[{"name":"t1","newick":"(A,(B,C));"}]"#;
        let json = write_nexus(taxa_json, trees_json);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let output = v["ok"].as_str().unwrap();
        assert!(output.contains("#NEXUS"));
        assert!(output.contains("TREE t1"));
    }

    // ── Simulation tests ─────────────────────────────────────────────────

    #[test]
    fn simulate_evolution_jc69() {
        let json = simulate_evolution("(A:0.1,B:0.1);", 100, "jc69", 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].is_object(), "expected ok result, got: {json}");
        assert_eq!(v["ok"]["names"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn simulate_evolution_sequences() {
        let json = simulate_evolution("(A:0.1,B:0.1);", 100, "jc69", 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let seqs = v["ok"]["sequences"].as_array().unwrap();
        assert_eq!(seqs.len(), 2);
        for seq in seqs {
            assert_eq!(seq.as_str().unwrap().len(), 100);
        }
    }

    #[test]
    fn simulate_coalescent_basic() {
        let json = simulate_coalescent(5, 1000.0, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].is_object(), "expected ok result, got: {json}");
        assert!(v["ok"]["newick"].as_str().unwrap().contains(';'));
    }

    #[test]
    fn simulate_coalescent_tree() {
        let json = simulate_coalescent(5, 1000.0, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"]["n_leaves"].as_u64().unwrap(), 5);
    }

    #[test]
    fn simulate_coalescent_growth_basic() {
        let json = simulate_coalescent_growth(5, 1000.0, 0.01, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].is_object(), "expected ok result, got: {json}");
        assert!(v["ok"]["newick"].as_str().unwrap().contains(';'));
        assert_eq!(v["ok"]["n_leaves"].as_u64().unwrap(), 5);
    }
}
