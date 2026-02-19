//! UniFrac distance metrics and Faith's phylogenetic diversity.
//!
//! Provides phylogenetic beta-diversity (UniFrac) and alpha-diversity (Faith's PD)
//! metrics that incorporate evolutionary relationships via a phylogenetic tree.
//!
//! - **Faith's PD** — sum of branch lengths connecting observed taxa
//! - **Unweighted UniFrac** — fraction of branch length unique to one sample
//! - **Weighted UniFrac** — abundance-weighted phylogenetic dissimilarity
//! - **Generalized UniFrac** — parameterized family (Chen et al. 2012)

use std::collections::{HashMap, HashSet};

use cyanea_core::{CyaneaError, Result};

use crate::tree::{NodeId, PhyloTree};

/// UniFrac method variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnifracMethod {
    /// Unweighted UniFrac (presence/absence).
    Unweighted,
    /// Weighted UniFrac (abundance-weighted).
    Weighted,
    /// Generalized UniFrac with parameter alpha.
    Generalized(f64),
}

/// Result of a pairwise UniFrac distance matrix computation.
#[derive(Debug, Clone)]
pub struct UnifracResult {
    /// Pairwise distance matrix (n × n).
    pub distances: Vec<Vec<f64>>,
    /// Sample names in matrix order.
    pub sample_names: Vec<String>,
}

/// Faith's Phylogenetic Diversity: sum of branch lengths in the minimum spanning
/// subtree connecting the observed taxa to the root.
///
/// # Arguments
///
/// * `tree` — a rooted phylogenetic tree
/// * `taxa` — set of taxon names present in the sample
///
/// # Errors
///
/// Returns an error if `taxa` is empty or no taxa match tree leaf names.
pub fn faiths_pd(tree: &PhyloTree, taxa: &HashSet<String>) -> Result<f64> {
    if taxa.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "faiths_pd: taxa set must be non-empty".into(),
        ));
    }

    // Map leaf names to node ids
    let leaf_set: HashSet<NodeId> = tree
        .nodes()
        .iter()
        .filter(|n| n.is_leaf())
        .filter(|n| n.name.as_ref().map_or(false, |name| taxa.contains(name)))
        .map(|n| n.id)
        .collect();

    if leaf_set.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "faiths_pd: no taxa match tree leaf names".into(),
        ));
    }

    // Post-order: mark nodes that have a marked descendant (or are marked leaves)
    let n = tree.node_count();
    let mut has_marked = vec![false; n];

    for id in tree.iter_postorder() {
        if leaf_set.contains(&id) {
            has_marked[id] = true;
        } else {
            let node = tree.get_node(id).unwrap();
            has_marked[id] = node.children.iter().any(|&c| has_marked[c]);
        }
    }

    // Sum branch lengths of marked nodes (excluding root, whose branch goes nowhere)
    let mut pd = 0.0;
    for id in 0..n {
        if has_marked[id] && id != tree.root() {
            pd += tree.get_node(id).unwrap().branch_length.unwrap_or(0.0);
        }
    }

    Ok(pd)
}

/// Unweighted UniFrac distance between two samples.
///
/// `U = unique_branch_length / total_branch_length` where unique branches
/// are those leading to taxa present in one sample but not the other.
///
/// # Arguments
///
/// * `tree` — a rooted phylogenetic tree
/// * `sample_a`, `sample_b` — sparse abundance maps (taxon name → count)
///
/// # Errors
///
/// Returns an error if either sample is empty or no taxa match tree leaves.
pub fn unweighted_unifrac(
    tree: &PhyloTree,
    sample_a: &HashMap<String, f64>,
    sample_b: &HashMap<String, f64>,
) -> Result<f64> {
    validate_samples(sample_a, sample_b)?;

    let n = tree.node_count();
    let leaf_map = build_leaf_map(tree);

    // Post-order: propagate presence
    let mut in_a = vec![false; n];
    let mut in_b = vec![false; n];

    for id in tree.iter_postorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_leaf() {
            if let Some(name) = &node.name {
                in_a[id] = sample_a.get(name).map_or(false, |&v| v > 0.0);
                in_b[id] = sample_b.get(name).map_or(false, |&v| v > 0.0);
            }
        } else {
            in_a[id] = node.children.iter().any(|&c| in_a[c]);
            in_b[id] = node.children.iter().any(|&c| in_b[c]);
        }
    }

    let mut unique_bl = 0.0;
    let mut total_bl = 0.0;
    for id in 0..n {
        if id == tree.root() {
            continue;
        }
        let bl = tree.get_node(id).unwrap().branch_length.unwrap_or(0.0);
        if in_a[id] || in_b[id] {
            total_bl += bl;
            if in_a[id] != in_b[id] {
                unique_bl += bl;
            }
        }
    }

    if total_bl == 0.0 {
        return Ok(0.0);
    }
    let _ = leaf_map; // suppress unused warning
    Ok(unique_bl / total_bl)
}

/// Weighted UniFrac distance between two samples.
///
/// `W = Σ bl * |pA - pB| / Σ bl * (pA + pB)` where proportions are
/// propagated up from leaves.
///
/// # Errors
///
/// Returns an error if either sample is empty or no taxa match tree leaves.
pub fn weighted_unifrac(
    tree: &PhyloTree,
    sample_a: &HashMap<String, f64>,
    sample_b: &HashMap<String, f64>,
) -> Result<f64> {
    validate_samples(sample_a, sample_b)?;

    let n = tree.node_count();

    // Compute totals for normalization to proportions
    let total_a: f64 = sample_a.values().sum();
    let total_b: f64 = sample_b.values().sum();
    if total_a == 0.0 && total_b == 0.0 {
        return Ok(0.0);
    }

    // Post-order: propagate proportions
    let mut prop_a = vec![0.0; n];
    let mut prop_b = vec![0.0; n];

    for id in tree.iter_postorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_leaf() {
            if let Some(name) = &node.name {
                if total_a > 0.0 {
                    prop_a[id] = sample_a.get(name).copied().unwrap_or(0.0) / total_a;
                }
                if total_b > 0.0 {
                    prop_b[id] = sample_b.get(name).copied().unwrap_or(0.0) / total_b;
                }
            }
        } else {
            prop_a[id] = node.children.iter().map(|&c| prop_a[c]).sum();
            prop_b[id] = node.children.iter().map(|&c| prop_b[c]).sum();
        }
    }

    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for id in 0..n {
        if id == tree.root() {
            continue;
        }
        let bl = tree.get_node(id).unwrap().branch_length.unwrap_or(0.0);
        numerator += bl * (prop_a[id] - prop_b[id]).abs();
        denominator += bl * (prop_a[id] + prop_b[id]);
    }

    if denominator == 0.0 {
        return Ok(0.0);
    }
    Ok(numerator / denominator)
}

/// Generalized UniFrac distance (Chen et al. 2012).
///
/// Uses `(pA + pB)^alpha` weighting. Special cases:
/// - alpha = 0: approximates unweighted UniFrac
/// - alpha = 1: equivalent to weighted UniFrac
///
/// # Errors
///
/// Returns an error if either sample is empty or alpha is negative.
pub fn generalized_unifrac(
    tree: &PhyloTree,
    sample_a: &HashMap<String, f64>,
    sample_b: &HashMap<String, f64>,
    alpha: f64,
) -> Result<f64> {
    if alpha < 0.0 {
        return Err(CyaneaError::InvalidInput(
            "generalized_unifrac: alpha must be non-negative".into(),
        ));
    }
    validate_samples(sample_a, sample_b)?;

    let n = tree.node_count();

    let total_a: f64 = sample_a.values().sum();
    let total_b: f64 = sample_b.values().sum();
    if total_a == 0.0 && total_b == 0.0 {
        return Ok(0.0);
    }

    // Post-order: propagate proportions
    let mut prop_a = vec![0.0; n];
    let mut prop_b = vec![0.0; n];

    for id in tree.iter_postorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_leaf() {
            if let Some(name) = &node.name {
                if total_a > 0.0 {
                    prop_a[id] = sample_a.get(name).copied().unwrap_or(0.0) / total_a;
                }
                if total_b > 0.0 {
                    prop_b[id] = sample_b.get(name).copied().unwrap_or(0.0) / total_b;
                }
            }
        } else {
            prop_a[id] = node.children.iter().map(|&c| prop_a[c]).sum();
            prop_b[id] = node.children.iter().map(|&c| prop_b[c]).sum();
        }
    }

    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for id in 0..n {
        if id == tree.root() {
            continue;
        }
        let bl = tree.get_node(id).unwrap().branch_length.unwrap_or(0.0);
        let sum_prop = prop_a[id] + prop_b[id];
        if sum_prop > 0.0 {
            let weight = sum_prop.powf(alpha);
            let diff = (prop_a[id] - prop_b[id]).abs() / sum_prop;
            numerator += bl * weight * diff;
            denominator += bl * weight;
        }
    }

    if denominator == 0.0 {
        return Ok(0.0);
    }
    Ok(numerator / denominator)
}

/// Compute a pairwise UniFrac distance matrix for multiple samples.
///
/// # Arguments
///
/// * `tree` — a rooted phylogenetic tree
/// * `samples` — named samples with sparse abundance data
/// * `method` — UniFrac variant to use
///
/// # Errors
///
/// Returns an error if fewer than 2 samples or computation fails.
pub fn unifrac_matrix(
    tree: &PhyloTree,
    samples: &[(&str, HashMap<String, f64>)],
    method: UnifracMethod,
) -> Result<UnifracResult> {
    if samples.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "unifrac_matrix: at least 2 samples required".into(),
        ));
    }

    let n = samples.len();
    let mut distances = vec![vec![0.0; n]; n];
    let sample_names: Vec<String> = samples.iter().map(|(name, _)| name.to_string()).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            let d = match method {
                UnifracMethod::Unweighted => {
                    unweighted_unifrac(tree, &samples[i].1, &samples[j].1)?
                }
                UnifracMethod::Weighted => {
                    weighted_unifrac(tree, &samples[i].1, &samples[j].1)?
                }
                UnifracMethod::Generalized(alpha) => {
                    generalized_unifrac(tree, &samples[i].1, &samples[j].1, alpha)?
                }
            };
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }

    Ok(UnifracResult {
        distances,
        sample_names,
    })
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn validate_samples(a: &HashMap<String, f64>, b: &HashMap<String, f64>) -> Result<()> {
    if a.is_empty() && b.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "both samples are empty".into(),
        ));
    }
    Ok(())
}

fn build_leaf_map(tree: &PhyloTree) -> HashMap<String, NodeId> {
    let mut map = HashMap::new();
    for node in tree.nodes() {
        if node.is_leaf() {
            if let Some(ref name) = node.name {
                map.insert(name.clone(), node.id);
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_star_tree() -> PhyloTree {
        // Star tree: root → A:1.0, B:1.0, C:1.0
        let mut tree = PhyloTree::new();
        tree.add_child(0, Some("A".into()), Some(1.0)).unwrap();
        tree.add_child(0, Some("B".into()), Some(1.0)).unwrap();
        tree.add_child(0, Some("C".into()), Some(1.0)).unwrap();
        tree
    }

    fn make_balanced_tree() -> PhyloTree {
        // ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6)
        PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap()
    }

    fn sample(pairs: &[(&str, f64)]) -> HashMap<String, f64> {
        pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
    }

    // ── Faith's PD ───────────────────────────────────────────────────

    #[test]
    fn faiths_pd_all_taxa_equals_total_bl() {
        let tree = make_balanced_tree();
        let taxa: HashSet<String> = ["A", "B", "C", "D"].iter().map(|s| s.to_string()).collect();
        let pd = faiths_pd(&tree, &taxa).unwrap();
        let total = tree.total_branch_length();
        assert!(
            (pd - total).abs() < 1e-10,
            "PD={} total_bl={}",
            pd,
            total
        );
    }

    #[test]
    fn faiths_pd_single_taxon_root_to_leaf() {
        let tree = make_balanced_tree();
        let taxa: HashSet<String> = ["A"].iter().map(|s| s.to_string()).collect();
        let pd = faiths_pd(&tree, &taxa).unwrap();
        // Path from A to root: A:0.1 + AB:0.3 = 0.4
        assert!((pd - 0.4).abs() < 1e-10, "PD={}", pd);
    }

    #[test]
    fn faiths_pd_empty_taxa_error() {
        let tree = make_balanced_tree();
        let taxa: HashSet<String> = HashSet::new();
        assert!(faiths_pd(&tree, &taxa).is_err());
    }

    #[test]
    fn faiths_pd_no_match_error() {
        let tree = make_balanced_tree();
        let taxa: HashSet<String> = ["Z"].iter().map(|s| s.to_string()).collect();
        assert!(faiths_pd(&tree, &taxa).is_err());
    }

    // ── Unweighted UniFrac ───────────────────────────────────────────

    #[test]
    fn unweighted_identical_is_zero() {
        let tree = make_balanced_tree();
        let s = sample(&[("A", 1.0), ("B", 1.0), ("C", 1.0), ("D", 1.0)]);
        let d = unweighted_unifrac(&tree, &s, &s).unwrap();
        assert!(d.abs() < 1e-10, "d={}", d);
    }

    #[test]
    fn unweighted_disjoint_is_one_star() {
        let tree = make_star_tree();
        let sa = sample(&[("A", 1.0)]);
        let sb = sample(&[("B", 1.0)]);
        let d = unweighted_unifrac(&tree, &sa, &sb).unwrap();
        assert!((d - 1.0).abs() < 1e-10, "d={}", d);
    }

    #[test]
    fn unweighted_symmetric() {
        let tree = make_balanced_tree();
        let sa = sample(&[("A", 1.0), ("B", 1.0)]);
        let sb = sample(&[("C", 1.0), ("D", 1.0)]);
        let d1 = unweighted_unifrac(&tree, &sa, &sb).unwrap();
        let d2 = unweighted_unifrac(&tree, &sb, &sa).unwrap();
        assert!((d1 - d2).abs() < 1e-10);
    }

    // ── Weighted UniFrac ─────────────────────────────────────────────

    #[test]
    fn weighted_identical_is_zero() {
        let tree = make_balanced_tree();
        let s = sample(&[("A", 10.0), ("B", 20.0), ("C", 30.0), ("D", 40.0)]);
        let d = weighted_unifrac(&tree, &s, &s).unwrap();
        assert!(d.abs() < 1e-10, "d={}", d);
    }

    #[test]
    fn weighted_different_abundances() {
        let tree = make_star_tree();
        let sa = sample(&[("A", 100.0), ("B", 0.0), ("C", 0.0)]);
        let sb = sample(&[("A", 0.0), ("B", 100.0), ("C", 0.0)]);
        let d = weighted_unifrac(&tree, &sa, &sb).unwrap();
        // On star tree, fully disjoint → d = 1.0
        assert!((d - 1.0).abs() < 1e-10, "d={}", d);
    }

    // ── Generalized UniFrac ──────────────────────────────────────────

    #[test]
    fn generalized_alpha1_approx_weighted() {
        let tree = make_balanced_tree();
        let sa = sample(&[("A", 10.0), ("B", 20.0)]);
        let sb = sample(&[("C", 30.0), ("D", 40.0)]);
        let guf = generalized_unifrac(&tree, &sa, &sb, 1.0).unwrap();
        let wuf = weighted_unifrac(&tree, &sa, &sb).unwrap();
        assert!(
            (guf - wuf).abs() < 1e-10,
            "generalized(α=1)={} weighted={}",
            guf,
            wuf
        );
    }

    #[test]
    fn generalized_negative_alpha_error() {
        let tree = make_balanced_tree();
        let s = sample(&[("A", 1.0)]);
        assert!(generalized_unifrac(&tree, &s, &s, -1.0).is_err());
    }

    // ── UniFrac matrix ───────────────────────────────────────────────

    #[test]
    fn unifrac_matrix_symmetric() {
        let tree = make_star_tree();
        let samples = vec![
            ("s1", sample(&[("A", 1.0), ("B", 1.0)])),
            ("s2", sample(&[("B", 1.0), ("C", 1.0)])),
            ("s3", sample(&[("A", 1.0), ("C", 1.0)])),
        ];
        let result = unifrac_matrix(&tree, &samples, UnifracMethod::Unweighted).unwrap();
        assert_eq!(result.sample_names, vec!["s1", "s2", "s3"]);
        for i in 0..3 {
            assert!(result.distances[i][i].abs() < 1e-10);
            for j in 0..3 {
                assert!((result.distances[i][j] - result.distances[j][i]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn unifrac_matrix_too_few_samples_error() {
        let tree = make_star_tree();
        let samples = vec![("s1", sample(&[("A", 1.0)]))];
        assert!(unifrac_matrix(&tree, &samples, UnifracMethod::Unweighted).is_err());
    }
}
