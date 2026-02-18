//! Bootstrap support values for phylogenetic trees.
//!
//! Implements non-parametric bootstrapping: resample alignment columns with
//! replacement, rebuild trees from each replicate, and compute the proportion
//! of replicates that recover each bipartition (split) from the original tree.

use std::collections::BTreeSet;

use cyanea_core::{CyaneaError, Result};

use crate::tree::PhyloTree;

/// Compute bootstrap support values for each internal edge of the original tree.
///
/// For each replicate:
/// 1. Resample alignment columns with replacement.
/// 2. Rebuild a tree from the resampled alignment using `tree_builder`.
/// 3. Compare bipartitions in the replicate tree to the original tree.
///
/// Returns a vector of support values (0.0 to 1.0), one per internal edge
/// of the original tree. The ordering corresponds to the internal edges as
/// encountered in a pre-order traversal (skipping root and leaves).
///
/// # Arguments
///
/// * `sequences` - Aligned sequences (one per taxon). All must have the same length.
/// * `tree_builder` - A function that builds a tree from resampled sequences.
///   The sequences are provided as owned `Vec<u8>` slices.
/// * `n_replicates` - Number of bootstrap replicates to generate.
///
/// # Errors
///
/// Returns an error if sequences are empty, have mismatched lengths, or if
/// `tree_builder` fails on any replicate.
pub fn bootstrap_support(
    sequences: &[&[u8]],
    original_tree: &PhyloTree,
    tree_builder: impl Fn(&[Vec<u8>]) -> Result<PhyloTree>,
    n_replicates: usize,
) -> Result<Vec<f64>> {
    if sequences.is_empty() {
        return Err(CyaneaError::InvalidInput("no sequences provided".into()));
    }
    let seq_len = sequences[0].len();
    if seq_len == 0 {
        return Err(CyaneaError::InvalidInput("empty sequences".into()));
    }
    for (i, seq) in sequences.iter().enumerate() {
        if seq.len() != seq_len {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence {} has length {}, expected {}",
                i,
                seq.len(),
                seq_len
            )));
        }
    }

    // Extract bipartitions from the original tree.
    let original_bps = bipartitions(original_tree);
    if original_bps.is_empty() {
        return Ok(Vec::new());
    }

    // Count how many times each original bipartition appears in replicate trees.
    let mut counts = vec![0usize; original_bps.len()];

    // Simple deterministic PRNG (xorshift64) for reproducibility without
    // external dependencies.
    let mut rng_state: u64 = 42;
    let mut next_rand = move || -> u64 {
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        rng_state
    };

    for _ in 0..n_replicates {
        // Resample columns with replacement.
        let mut resampled: Vec<Vec<u8>> = vec![Vec::with_capacity(seq_len); sequences.len()];
        for _ in 0..seq_len {
            let col = (next_rand() as usize) % seq_len;
            for (i, seq) in sequences.iter().enumerate() {
                resampled[i].push(seq[col]);
            }
        }

        // Build a tree from the resampled alignment.
        let replicate_tree = match tree_builder(&resampled) {
            Ok(t) => t,
            Err(_) => continue, // Skip failed replicates.
        };

        // Extract bipartitions from the replicate tree.
        let replicate_bps = bipartitions(&replicate_tree);
        let replicate_set: BTreeSet<&BTreeSet<String>> =
            replicate_bps.iter().collect();

        // Check which original bipartitions are present in this replicate.
        for (i, bp) in original_bps.iter().enumerate() {
            if replicate_set.contains(bp) {
                counts[i] += 1;
            }
        }
    }

    // Convert counts to proportions.
    let support: Vec<f64> = if n_replicates == 0 {
        vec![0.0; original_bps.len()]
    } else {
        counts
            .iter()
            .map(|&c| c as f64 / n_replicates as f64)
            .collect()
    };

    Ok(support)
}

/// Extract non-trivial bipartitions (splits) from a tree.
///
/// For each internal edge (non-root, non-leaf), collects the sorted set
/// of leaf names in the subtree below that edge. Trivial splits (single
/// leaf or all-but-one leaf) are excluded.
///
/// The bipartitions are returned in pre-order traversal order.
pub fn bipartitions(tree: &PhyloTree) -> Vec<BTreeSet<String>> {
    let all_leaves: BTreeSet<String> = tree.leaf_names().into_iter().collect();
    let n_leaves = all_leaves.len();
    let mut result = Vec::new();

    for node_id in tree.iter_preorder() {
        let node = tree.get_node(node_id).unwrap();
        // Skip leaves and the root.
        if node.is_leaf() || node.is_root() {
            continue;
        }
        let subtree_leaves = tree.subtree_leaf_names(node_id);
        // Skip trivial splits.
        if subtree_leaves.len() <= 1 || subtree_leaves.len() >= n_leaves - 1 {
            continue;
        }
        result.push(subtree_leaves);
    }

    result
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bipartitions_four_leaf_tree() {
        let tree = PhyloTree::from_newick("((A,B),(C,D));").unwrap();
        let bps = bipartitions(&tree);
        // For ((A,B),(C,D)), the non-trivial bipartitions are {A,B} and {C,D}.
        assert_eq!(bps.len(), 2);

        let ab: BTreeSet<String> = ["A", "B"].iter().map(|s| s.to_string()).collect();
        let cd: BTreeSet<String> = ["C", "D"].iter().map(|s| s.to_string()).collect();
        assert!(bps.contains(&ab), "missing {{A,B}} bipartition");
        assert!(bps.contains(&cd), "missing {{C,D}} bipartition");
    }

    #[test]
    fn bipartitions_three_leaf_tree() {
        // ((A,B),C) — the only internal edge gives split {A,B}|{C},
        // which is trivial (one side has a single leaf), so 0 non-trivial bipartitions.
        let tree = PhyloTree::from_newick("((A,B),C);").unwrap();
        let bps = bipartitions(&tree);
        assert_eq!(bps.len(), 0);
    }

    #[test]
    fn bipartitions_five_leaf_tree() {
        let tree = PhyloTree::from_newick("(((A,B),C),(D,E));").unwrap();
        let bps = bipartitions(&tree);
        // Non-trivial splits: {A,B}, {A,B,C}, {D,E}
        assert_eq!(bps.len(), 3);
    }

    #[test]
    fn bipartitions_star_tree() {
        // Star tree (A,B,C,D) has no internal edges, no non-trivial bipartitions.
        let tree = PhyloTree::from_newick("(A,B,C,D);").unwrap();
        let bps = bipartitions(&tree);
        assert!(bps.is_empty());
    }

    #[test]
    fn bootstrap_values_in_range() {
        // Simple test: 4 sequences, 4 taxa.
        let seqs: Vec<Vec<u8>> = vec![
            b"AAACCCAAACCC".to_vec(),
            b"AAACCCAAACCC".to_vec(),
            b"CCCAAACCCAAA".to_vec(),
            b"CCCAAACCCAAA".to_vec(),
        ];
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        // Build an original tree.
        let original = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();

        // A trivial tree builder that always returns the same topology.
        let builder = |_seqs: &[Vec<u8>]| -> Result<PhyloTree> {
            PhyloTree::from_newick("((A,B),(C,D));")
        };

        let support = bootstrap_support(&seq_refs, &original, builder, 100).unwrap();

        for &s in &support {
            assert!(
                (0.0..=1.0).contains(&s),
                "bootstrap support {} out of range",
                s
            );
        }
    }

    #[test]
    fn bootstrap_perfect_support() {
        // If the tree builder always returns the same topology as the original,
        // all bipartitions should have 100% support.
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"TGCATGCA".to_vec(),
            b"TGCATGCA".to_vec(),
        ];
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let original = PhyloTree::from_newick("((A,B),(C,D));").unwrap();

        let builder = |_seqs: &[Vec<u8>]| -> Result<PhyloTree> {
            PhyloTree::from_newick("((A,B),(C,D));")
        };

        let support = bootstrap_support(&seq_refs, &original, builder, 50).unwrap();

        for &s in &support {
            assert!(
                (s - 1.0).abs() < 1e-10,
                "expected 100% support, got {}",
                s
            );
        }
    }

    #[test]
    fn bootstrap_empty_sequences_error() {
        let seqs: Vec<&[u8]> = vec![];
        let tree = PhyloTree::from_newick("(A,B);").unwrap();
        let builder = |_: &[Vec<u8>]| -> Result<PhyloTree> {
            PhyloTree::from_newick("(A,B);")
        };
        assert!(bootstrap_support(&seqs, &tree, builder, 10).is_err());
    }

    #[test]
    fn bootstrap_zero_replicates() {
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
            b"TGCA".to_vec(),
            b"TGCA".to_vec(),
        ];
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let original = PhyloTree::from_newick("((A,B),(C,D));").unwrap();
        let builder = |_: &[Vec<u8>]| -> Result<PhyloTree> {
            PhyloTree::from_newick("((A,B),(C,D));")
        };

        let support = bootstrap_support(&seq_refs, &original, builder, 0).unwrap();
        for &s in &support {
            assert_eq!(s, 0.0);
        }
    }

    #[test]
    fn bootstrap_star_tree_returns_empty() {
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
            b"TGCA".to_vec(),
            b"TGCA".to_vec(),
        ];
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let original = PhyloTree::from_newick("(A,B,C,D);").unwrap();
        let builder = |_: &[Vec<u8>]| -> Result<PhyloTree> {
            PhyloTree::from_newick("(A,B,C,D);")
        };

        let support = bootstrap_support(&seq_refs, &original, builder, 10).unwrap();
        assert!(support.is_empty());
    }
}
