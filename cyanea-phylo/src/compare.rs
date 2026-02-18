//! Tree comparison metrics.
//!
//! Provides Robinson-Foulds distance and branch score distance for comparing
//! phylogenetic tree topologies.

use crate::tree::PhyloTree;
use cyanea_core::{CyaneaError, Result};
use std::collections::BTreeSet;

/// A bipartition (split) of leaf names induced by an internal edge,
/// paired with the branch length on that edge.
type Bipartition = (BTreeSet<String>, Option<f64>);

/// Robinson-Foulds (symmetric difference) distance between two trees.
///
/// Counts the number of bipartitions present in one tree but not the other.
/// Both trees must have the same set of leaf names.
pub fn robinson_foulds(t1: &PhyloTree, t2: &PhyloTree) -> Result<usize> {
    validate_same_leaves(t1, t2)?;
    let bp1 = extract_bipartitions(t1);
    let bp2 = extract_bipartitions(t2);

    let set1: BTreeSet<_> = bp1.iter().map(|(s, _)| s).collect();
    let set2: BTreeSet<_> = bp2.iter().map(|(s, _)| s).collect();

    let symmetric_diff = set1.symmetric_difference(&set2).count();
    Ok(symmetric_diff)
}

/// Normalized Robinson-Foulds distance (0.0 to 1.0).
///
/// For rooted binary trees with `n` leaves, the maximum RF distance is `2*(n-2)`.
/// Returns 0.0 for trees with fewer than 3 leaves (no non-trivial bipartitions).
pub fn robinson_foulds_normalized(t1: &PhyloTree, t2: &PhyloTree) -> Result<f64> {
    let n = t1.leaf_count();
    if n < 3 {
        validate_same_leaves(t1, t2)?;
        return Ok(0.0);
    }
    let rf = robinson_foulds(t1, t2)?;
    let max_rf = 2 * (n - 2);
    if max_rf == 0 {
        return Ok(0.0);
    }
    Ok(rf as f64 / max_rf as f64)
}

/// Branch score distance between two trees.
///
/// For each bipartition, computes the squared difference in branch lengths.
/// Bipartitions unique to one tree contribute their squared branch length.
/// Returns the square root of the sum.
pub fn branch_score_distance(t1: &PhyloTree, t2: &PhyloTree) -> Result<f64> {
    validate_same_leaves(t1, t2)?;
    let bp1 = extract_bipartitions(t1);
    let bp2 = extract_bipartitions(t2);

    let map1: std::collections::BTreeMap<_, _> = bp1.into_iter().collect();
    let map2: std::collections::BTreeMap<_, _> = bp2.into_iter().collect();

    let mut sum_sq = 0.0;

    // Bipartitions in tree 1
    for (split, len1) in &map1 {
        let l1 = len1.unwrap_or(0.0);
        match map2.get(split) {
            Some(len2) => {
                let l2 = len2.unwrap_or(0.0);
                sum_sq += (l1 - l2).powi(2);
            }
            None => {
                sum_sq += l1.powi(2);
            }
        }
    }

    // Bipartitions unique to tree 2
    for (split, len2) in &map2 {
        if !map1.contains_key(split) {
            let l2 = len2.unwrap_or(0.0);
            sum_sq += l2.powi(2);
        }
    }

    Ok(sum_sq.sqrt())
}

/// Extract non-trivial bipartitions from a tree.
///
/// For each internal edge (excluding the root split), collects the sorted set
/// of leaf names in the subtree below that edge.
fn extract_bipartitions(tree: &PhyloTree) -> Vec<Bipartition> {
    let all_leaves: BTreeSet<String> = tree.leaf_names().into_iter().collect();
    let n_leaves = all_leaves.len();
    let mut bipartitions = Vec::new();

    for node_id in tree.iter_preorder() {
        let node = tree.get_node(node_id).unwrap();
        // Skip leaves and the root
        if node.is_leaf() || node.is_root() {
            continue;
        }
        let subtree_leaves = tree.subtree_leaf_names(node_id);
        // Skip trivial splits (single leaf or all leaves minus one)
        if subtree_leaves.len() <= 1 || subtree_leaves.len() >= n_leaves - 1 {
            continue;
        }
        bipartitions.push((subtree_leaves, node.branch_length));
    }

    bipartitions
}

fn validate_same_leaves(t1: &PhyloTree, t2: &PhyloTree) -> Result<()> {
    let leaves1 = t1.leaf_names();
    let leaves2 = t2.leaf_names();
    if leaves1 != leaves2 {
        return Err(CyaneaError::InvalidInput(format!(
            "trees have different leaf sets: {:?} vs {:?}",
            leaves1, leaves2
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree_from_newick(s: &str) -> PhyloTree {
        crate::newick::parse(s).unwrap()
    }

    #[test]
    fn rf_identical_trees() {
        let t1 = tree_from_newick("((A,B),(C,D));");
        let t2 = tree_from_newick("((A,B),(C,D));");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
    }

    #[test]
    fn rf_different_topologies() {
        let t1 = tree_from_newick("((A,B),(C,D));");
        let t2 = tree_from_newick("((A,C),(B,D));");
        let rf = robinson_foulds(&t1, &t2).unwrap();
        assert!(rf > 0);
    }

    #[test]
    fn rf_five_leaves() {
        let t1 = tree_from_newick("(((A,B),C),(D,E));");
        let t2 = tree_from_newick("(((A,C),B),(D,E));");
        let rf = robinson_foulds(&t1, &t2).unwrap();
        assert!(rf > 0);
    }

    #[test]
    fn rf_normalized_range() {
        let t1 = tree_from_newick("((A,B),(C,D));");
        let t2 = tree_from_newick("((A,C),(B,D));");
        let nrf = robinson_foulds_normalized(&t1, &t2).unwrap();
        assert!(nrf >= 0.0 && nrf <= 1.0, "nRF = {}", nrf);
    }

    #[test]
    fn rf_normalized_identical() {
        let t1 = tree_from_newick("((A,B),(C,D));");
        let t2 = tree_from_newick("((A,B),(C,D));");
        assert_eq!(robinson_foulds_normalized(&t1, &t2).unwrap(), 0.0);
    }

    #[test]
    fn rf_normalized_small_tree() {
        let t1 = tree_from_newick("(A,B,C);");
        let t2 = tree_from_newick("(A,B,C);");
        // < 4 leaves: always 0.0
        assert_eq!(robinson_foulds_normalized(&t1, &t2).unwrap(), 0.0);
    }

    #[test]
    fn branch_score_identical() {
        let t1 = tree_from_newick("((A:1,B:2):3,(C:4,D:5):6);");
        let t2 = tree_from_newick("((A:1,B:2):3,(C:4,D:5):6);");
        let bs = branch_score_distance(&t1, &t2).unwrap();
        assert!(bs.abs() < 1e-10, "bs = {}", bs);
    }

    #[test]
    fn branch_score_different_lengths() {
        let t1 = tree_from_newick("((A:1,B:2):3,(C:4,D:5):6);");
        let t2 = tree_from_newick("((A:1,B:2):5,(C:4,D:5):8);");
        let bs = branch_score_distance(&t1, &t2).unwrap();
        assert!(bs > 0.0, "bs should be positive");
    }

    #[test]
    fn branch_score_different_topology() {
        let t1 = tree_from_newick("((A:1,B:1):1,(C:1,D:1):1);");
        let t2 = tree_from_newick("((A:1,C:1):1,(B:1,D:1):1);");
        let bs = branch_score_distance(&t1, &t2).unwrap();
        assert!(bs > 0.0);
    }

    #[test]
    fn leaf_mismatch_error() {
        let t1 = tree_from_newick("((A,B),(C,D));");
        let t2 = tree_from_newick("((A,B),(C,E));");
        assert!(robinson_foulds(&t1, &t2).is_err());
        assert!(branch_score_distance(&t1, &t2).is_err());
    }
}
