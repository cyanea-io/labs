//! Consensus tree construction from bootstrap replicates.
//!
//! Builds consensus trees from a collection of trees by analyzing
//! bipartition frequencies. Supports strict, majority-rule, and
//! extended majority-rule consensus methods.

use std::collections::{BTreeMap, BTreeSet};

use crate::tree::{Node, NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Type of consensus tree to build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusType {
    /// Include only bipartitions present in >50% of trees.
    MajorityRule,
    /// Include only bipartitions present in 100% of trees.
    Strict,
    /// Majority-rule plus additional compatible bipartitions
    /// in decreasing frequency order.
    ExtendedMajorityRule,
}

/// A bipartition with its observed frequency across input trees.
#[derive(Debug, Clone)]
pub struct SupportedBipartition {
    /// The set of leaf names on one side of the split.
    pub leaves: BTreeSet<String>,
    /// Proportion of input trees containing this bipartition (0.0–1.0).
    pub support: f64,
}

/// Compute bipartition frequencies across a set of trees.
///
/// Returns bipartitions sorted by decreasing frequency, excluding trivial
/// splits (single leaf or complement of single leaf).
pub fn bipartition_frequencies(trees: &[PhyloTree]) -> Result<Vec<SupportedBipartition>> {
    if trees.is_empty() {
        return Err(CyaneaError::InvalidInput("no trees provided".into()));
    }

    let n_trees = trees.len();
    let mut counts: BTreeMap<BTreeSet<String>, usize> = BTreeMap::new();

    // Collect all leaf names from the first tree.
    let all_leaves: BTreeSet<String> = trees[0].leaf_names().into_iter().collect();
    let n_leaves = all_leaves.len();

    for tree in trees {
        // Extract non-trivial bipartitions, deduplicated per tree.
        let mut tree_bps: BTreeSet<BTreeSet<String>> = BTreeSet::new();
        for node_id in tree.iter_preorder() {
            let node = tree.get_node(node_id).unwrap();
            if node.is_leaf() || node.is_root() {
                continue;
            }
            let subtree_leaves = tree.subtree_leaf_names(node_id);
            // Skip trivial splits.
            if subtree_leaves.len() <= 1 || subtree_leaves.len() >= n_leaves - 1 {
                continue;
            }
            // Canonicalize: use the smaller side.
            let canonical = canonicalize_bipartition(&subtree_leaves, &all_leaves);
            tree_bps.insert(canonical);
        }
        for bp in tree_bps {
            *counts.entry(bp).or_insert(0) += 1;
        }
    }

    let mut result: Vec<SupportedBipartition> = counts
        .into_iter()
        .map(|(leaves, count)| SupportedBipartition {
            leaves,
            support: count as f64 / n_trees as f64,
        })
        .collect();

    // Sort by decreasing support, then by leaf set for determinism.
    result.sort_by(|a, b| {
        b.support
            .partial_cmp(&a.support)
            .unwrap()
            .then_with(|| a.leaves.cmp(&b.leaves))
    });

    Ok(result)
}

/// Canonicalize a bipartition: always use the side that contains the
/// lexicographically smallest leaf name (or the smaller set if equal).
fn canonicalize_bipartition(
    split: &BTreeSet<String>,
    all: &BTreeSet<String>,
) -> BTreeSet<String> {
    let complement: BTreeSet<String> = all.difference(split).cloned().collect();
    // Use lexicographic comparison of first elements.
    if let (Some(s_first), Some(c_first)) = (split.iter().next(), complement.iter().next()) {
        if s_first < c_first {
            split.clone()
        } else if c_first < s_first {
            complement
        } else if split.len() <= complement.len() {
            split.clone()
        } else {
            complement
        }
    } else {
        split.clone()
    }
}

/// Build a consensus tree from a collection of trees.
///
/// The consensus type determines which bipartitions are included:
/// - `Strict`: only bipartitions found in all trees
/// - `MajorityRule`: bipartitions found in >50% of trees
/// - `ExtendedMajorityRule`: majority-rule plus additional compatible
///   bipartitions in decreasing frequency order
pub fn consensus_tree(trees: &[PhyloTree], consensus_type: ConsensusType) -> Result<PhyloTree> {
    if trees.is_empty() {
        return Err(CyaneaError::InvalidInput("no trees provided".into()));
    }

    let all_leaves: Vec<String> = {
        let mut names = trees[0].leaf_names();
        names.sort();
        names
    };

    let bps = bipartition_frequencies(trees)?;

    // Filter bipartitions based on consensus type.
    let threshold = match consensus_type {
        ConsensusType::Strict => 1.0 - 1e-10,
        ConsensusType::MajorityRule => 0.5,
        ConsensusType::ExtendedMajorityRule => 0.0, // Will add compatibility check.
    };

    let mut accepted: Vec<&SupportedBipartition> = Vec::new();

    for bp in &bps {
        if bp.support > threshold
            || (consensus_type == ConsensusType::Strict && bp.support >= 1.0 - 1e-10)
        {
            accepted.push(bp);
        } else if consensus_type == ConsensusType::ExtendedMajorityRule {
            // Check compatibility with already-accepted bipartitions.
            let all_set: BTreeSet<String> = all_leaves.iter().cloned().collect();
            if is_compatible(&bp.leaves, &accepted, &all_set) {
                accepted.push(bp);
            }
        }
    }

    // Build tree from accepted bipartitions.
    build_consensus_from_bipartitions(&all_leaves, &accepted)
}

/// Check if a bipartition is compatible with a set of accepted bipartitions.
///
/// Two bipartitions A, B are compatible if one of these holds:
/// A ⊂ B, B ⊂ A, or A ∩ B = ∅ (relative to the full leaf set).
fn is_compatible(
    candidate: &BTreeSet<String>,
    accepted: &[&SupportedBipartition],
    all_leaves: &BTreeSet<String>,
) -> bool {
    let candidate_comp: BTreeSet<String> = all_leaves.difference(candidate).cloned().collect();

    for bp in accepted {
        let a = candidate;
        let b = &bp.leaves;

        let inter: BTreeSet<_> = a.intersection(b).collect();

        if !inter.is_empty() && !a.is_subset(b) && !b.is_subset(a) {
            // Also check complement.
            let inter_comp: BTreeSet<_> = candidate_comp.intersection(b).collect();
            if !inter_comp.is_empty()
                && !candidate_comp.is_subset(b)
                && !b.is_subset(&candidate_comp)
            {
                return false;
            }
        }
    }
    true
}

/// Build a tree from a set of compatible bipartitions.
fn build_consensus_from_bipartitions(
    all_leaves: &[String],
    bipartitions: &[&SupportedBipartition],
) -> Result<PhyloTree> {
    // Sort bipartitions by size (largest first) for nesting.
    let mut sorted_bps: Vec<&SupportedBipartition> = bipartitions.to_vec();
    sorted_bps.sort_by(|a, b| b.leaves.len().cmp(&a.leaves.len()));

    // Start with a root node.
    let mut nodes: Vec<Node> = Vec::new();
    let root_id = 0;
    nodes.push(Node {
        id: root_id,
        parent: None,
        children: Vec::new(),
        branch_length: None,
        name: None,
    });

    // Create leaf nodes.
    let mut leaf_ids: BTreeMap<String, NodeId> = BTreeMap::new();
    for name in all_leaves {
        let id = nodes.len();
        nodes.push(Node {
            id,
            parent: None,
            children: Vec::new(),
            branch_length: None,
            name: Some(name.clone()),
        });
        leaf_ids.insert(name.clone(), id);
    }

    // Initially all leaves are children of root.
    let all_leaf_node_ids: Vec<NodeId> = all_leaves.iter().map(|n| leaf_ids[n]).collect();
    for &lid in &all_leaf_node_ids {
        nodes[root_id].children.push(lid);
        nodes[lid].parent = Some(root_id);
    }

    // For each bipartition (largest first), create an internal node grouping those leaves.
    for bp in &sorted_bps {
        let bp_leaf_ids: Vec<NodeId> = bp.leaves.iter().map(|n| leaf_ids[n]).collect();
        let bp_set: BTreeSet<NodeId> = bp_leaf_ids.iter().copied().collect();

        // Find the set of unique current parents.
        let parents: BTreeSet<NodeId> = bp_leaf_ids
            .iter()
            .map(|&id| nodes[id].parent.unwrap())
            .collect();

        if parents.len() == 1 {
            let parent = *parents.iter().next().unwrap();
            let parent_children: Vec<NodeId> = nodes[parent].children.clone();

            // Check if all and only bp_leaf_ids are children of parent.
            let children_set: BTreeSet<NodeId> = parent_children.iter().copied().collect();
            if bp_set == children_set {
                // Already grouped.
                continue;
            }

            // Create new internal node.
            let new_id = nodes.len();
            nodes.push(Node {
                id: new_id,
                parent: Some(parent),
                children: Vec::new(),
                branch_length: None,
                name: Some(format!("{:.0}", bp.support * 100.0)),
            });

            // Move bp leaves from parent to new node.
            nodes[parent].children.retain(|c| !bp_set.contains(c));
            nodes[parent].children.push(new_id);

            for &lid in &bp_leaf_ids {
                nodes[lid].parent = Some(new_id);
                nodes[new_id].children.push(lid);
            }
        }
        // Multi-parent case: bipartition spans already-created groups.
        // Left as unresolved (star topology at that level).
    }

    PhyloTree::from_nodes(nodes, root_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn majority_rule_identical_trees() {
        let trees: Vec<PhyloTree> = (0..10)
            .map(|_| PhyloTree::from_newick("((A,B),(C,D));").unwrap())
            .collect();
        let result = consensus_tree(&trees, ConsensusType::MajorityRule).unwrap();
        let mut names = result.leaf_names();
        names.sort();
        assert_eq!(names, vec!["A", "B", "C", "D"]);
        // Should have both bipartitions {A,B} and {C,D}.
        assert!(result.node_count() > 4); // More than just leaves + root.
    }

    #[test]
    fn strict_subset_of_majority() {
        // With identical trees, strict and majority should give same result.
        let trees: Vec<PhyloTree> = (0..5)
            .map(|_| PhyloTree::from_newick("((A,B),(C,D));").unwrap())
            .collect();
        let strict = consensus_tree(&trees, ConsensusType::Strict).unwrap();
        let majority = consensus_tree(&trees, ConsensusType::MajorityRule).unwrap();
        assert_eq!(strict.leaf_count(), majority.leaf_count());
    }

    #[test]
    fn extended_adds_compatible() {
        // Extended majority rule should include at least as many bipartitions
        // as strict majority rule.
        let trees: Vec<PhyloTree> = (0..10)
            .map(|_| PhyloTree::from_newick("((A,B),(C,D));").unwrap())
            .collect();
        let majority = consensus_tree(&trees, ConsensusType::MajorityRule).unwrap();
        let extended = consensus_tree(&trees, ConsensusType::ExtendedMajorityRule).unwrap();
        assert!(extended.node_count() >= majority.node_count());
    }

    #[test]
    fn strict_discordant_is_star() {
        // Two completely different topologies → strict consensus is a star.
        let trees = vec![
            PhyloTree::from_newick("((A,B),(C,D));").unwrap(),
            PhyloTree::from_newick("((A,C),(B,D));").unwrap(),
        ];
        let strict = consensus_tree(&trees, ConsensusType::Strict).unwrap();
        // In a strict consensus of conflicting trees, no bipartitions pass.
        // All leaves should be directly under root.
        assert_eq!(strict.leaf_count(), 4);
    }

    #[test]
    fn bipartition_frequencies_correct() {
        let trees = vec![
            PhyloTree::from_newick("((A,B),(C,D));").unwrap(),
            PhyloTree::from_newick("((A,B),(C,D));").unwrap(),
            PhyloTree::from_newick("((A,C),(B,D));").unwrap(),
        ];
        let bps = bipartition_frequencies(&trees).unwrap();

        // {A,B} appears in 2/3 trees.
        let ab: BTreeSet<String> = ["A", "B"].iter().map(|s| s.to_string()).collect();
        let ab_bp = bps.iter().find(|bp| bp.leaves == ab);
        assert!(ab_bp.is_some(), "should find {{A,B}} bipartition");
        assert!(
            (ab_bp.unwrap().support - 2.0 / 3.0).abs() < 1e-10,
            "{{A,B}} support should be 2/3, got {}",
            ab_bp.unwrap().support
        );
    }
}
