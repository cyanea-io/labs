//! Maximum likelihood phylogenetics: Felsenstein pruning and NNI search.
//!
//! Provides tree likelihood computation via Felsenstein's pruning algorithm
//! and hill-climbing tree search via Nearest Neighbor Interchange (NNI).

use crate::models::{nucleotide_index, JC69_FREQ, NUM_STATES};
use crate::tree::{Node, NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Compute the log-likelihood of a phylogenetic tree given aligned sequences
/// under a substitution model.
///
/// Uses Felsenstein's pruning algorithm (post-order traversal) to compute
/// partial likelihoods at each internal node, then sums log-likelihoods
/// across all alignment columns.
///
/// # Arguments
///
/// * `tree` - A rooted phylogenetic tree with branch lengths.
/// * `sequences` - Aligned sequences, one per leaf, in the same order as
///   the tree's sorted leaf names. All sequences must be the same length.
/// * `model` - A substitution model function that takes a branch length
///   and returns a 4x4 transition probability matrix.
///
/// # Errors
///
/// Returns an error if the number of sequences doesn't match the leaf count,
/// sequences have different lengths, or sequences are empty.
pub fn tree_likelihood(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: fn(f64) -> [[f64; 4]; 4],
) -> Result<f64> {
    let leaves = tree.leaves();
    let n_leaves = leaves.len();

    if sequences.len() != n_leaves {
        return Err(CyaneaError::InvalidInput(format!(
            "expected {} sequences for {} leaves, got {}",
            n_leaves,
            n_leaves,
            sequences.len()
        )));
    }
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

    // Build a mapping from leaf NodeId to sequence index.
    // Leaves are sorted by name to match the sequence order.
    let mut leaf_ids_sorted: Vec<NodeId> = leaves.clone();
    leaf_ids_sorted.sort_by(|&a, &b| {
        let name_a = tree
            .get_node(a)
            .and_then(|n| n.name.as_ref())
            .cloned()
            .unwrap_or_default();
        let name_b = tree
            .get_node(b)
            .and_then(|n| n.name.as_ref())
            .cloned()
            .unwrap_or_default();
        name_a.cmp(&name_b)
    });

    let n_nodes = tree.node_count();
    let mut leaf_seq_index = vec![0usize; n_nodes];
    for (seq_idx, &node_id) in leaf_ids_sorted.iter().enumerate() {
        leaf_seq_index[node_id] = seq_idx;
    }

    // Precompute transition probability matrices for each node's branch length.
    let mut trans_probs: Vec<Option<[[f64; 4]; 4]>> = vec![None; n_nodes];
    for id in 0..n_nodes {
        if let Some(node) = tree.get_node(id) {
            if let Some(bl) = node.branch_length {
                trans_probs[id] = Some(model(bl));
            }
        }
    }

    // Compute log-likelihood over all sites.
    let mut total_log_likelihood = 0.0;

    // Partial likelihoods: partials[node][state]
    let mut partials = vec![[0.0f64; NUM_STATES]; n_nodes];

    for site in 0..seq_len {
        // Initialize partials for all nodes.
        for p in partials.iter_mut() {
            *p = [0.0; NUM_STATES];
        }

        // Set leaf partials.
        for &leaf_id in &leaves {
            let seq_idx = leaf_seq_index[leaf_id];
            let base = sequences[seq_idx][site];
            if let Some(state) = nucleotide_index(base) {
                partials[leaf_id][state] = 1.0;
            } else {
                // Ambiguous or unknown base: equal probability for all states.
                partials[leaf_id] = [1.0; NUM_STATES];
            }
        }

        // Post-order traversal: compute internal node partials.
        for id in tree.iter_postorder() {
            let node = tree.get_node(id).unwrap();
            if node.is_leaf() {
                continue;
            }

            let mut node_partial = [1.0f64; NUM_STATES];

            for &child_id in &node.children {
                // Ensure we have a transition probability matrix for this edge.
                if trans_probs[child_id].is_none() {
                    // No branch length: use a small default.
                    trans_probs[child_id] = Some(model(1e-6));
                }
                let prob_matrix = trans_probs[child_id].as_ref().unwrap();
                let child_partial = &partials[child_id];

                for s in 0..NUM_STATES {
                    let mut sum = 0.0;
                    for t in 0..NUM_STATES {
                        sum += prob_matrix[s][t] * child_partial[t];
                    }
                    node_partial[s] *= sum;
                }
            }

            partials[id] = node_partial;
        }

        // Site likelihood at the root with uniform prior.
        let root_id = tree.root();
        let root_partial = &partials[root_id];
        let mut site_likelihood = 0.0;
        for s in 0..NUM_STATES {
            site_likelihood += JC69_FREQ * root_partial[s];
        }

        if site_likelihood <= 0.0 {
            // This shouldn't happen with valid data, but guard against it.
            return Err(CyaneaError::InvalidInput(format!(
                "zero or negative site likelihood at site {}",
                site
            )));
        }

        total_log_likelihood += site_likelihood.ln();
    }

    Ok(total_log_likelihood)
}

/// Perform a hill-climbing search using Nearest Neighbor Interchange (NNI)
/// to find a tree with improved likelihood.
///
/// For each internal edge, tries the two possible NNI rearrangements. Keeps
/// the rearrangement that most improves the log-likelihood. Repeats until
/// no further improvement is found.
///
/// # Arguments
///
/// * `tree` - An initial rooted phylogenetic tree with branch lengths.
/// * `sequences` - Aligned sequences, one per leaf, sorted by leaf name.
/// * `model` - A substitution model function.
///
/// # Returns
///
/// The tree with the best log-likelihood found by the NNI search.
pub fn nni_search(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: fn(f64) -> [[f64; 4]; 4],
) -> Result<PhyloTree> {
    let mut best_tree = tree.clone();
    let mut best_ll = tree_likelihood(&best_tree, sequences, model)?;

    loop {
        let mut improved = false;

        // Collect internal edges: edges where the child is an internal node
        // (has children) and itself is not the root.
        let internal_edges: Vec<(NodeId, NodeId)> = {
            let mut edges = Vec::new();
            for id in 0..best_tree.node_count() {
                let node = best_tree.get_node(id).unwrap();
                if !node.is_leaf() && !node.is_root() && node.children.len() == 2 {
                    if let Some(parent_id) = node.parent {
                        let parent = best_tree.get_node(parent_id).unwrap();
                        if parent.children.len() == 2 {
                            edges.push((parent_id, id));
                        }
                    }
                }
            }
            edges
        };

        for &(parent_id, child_id) in &internal_edges {
            let parent = best_tree.get_node(parent_id).unwrap();
            let child = best_tree.get_node(child_id).unwrap();

            // Identify the four subtrees around this internal edge.
            // Parent has two children: child_id and sibling.
            // Child has two children: nephew_a and nephew_b.
            let sibling_id = if parent.children[0] == child_id {
                parent.children[1]
            } else {
                parent.children[0]
            };
            let nephew_a = child.children[0];
            let nephew_b = child.children[1];

            // NNI swap 1: swap sibling with nephew_a
            if let Ok(swapped) =
                nni_swap(&best_tree, parent_id, child_id, sibling_id, nephew_a)
            {
                if let Ok(ll) = tree_likelihood(&swapped, sequences, model) {
                    if ll > best_ll {
                        best_tree = swapped;
                        best_ll = ll;
                        improved = true;
                        break; // Restart from the beginning
                    }
                }
            }

            // NNI swap 2: swap sibling with nephew_b
            if let Ok(swapped) =
                nni_swap(&best_tree, parent_id, child_id, sibling_id, nephew_b)
            {
                if let Ok(ll) = tree_likelihood(&swapped, sequences, model) {
                    if ll > best_ll {
                        best_tree = swapped;
                        best_ll = ll;
                        improved = true;
                        break; // Restart from the beginning
                    }
                }
            }
        }

        if !improved {
            break;
        }
    }

    Ok(best_tree)
}

/// Perform a single NNI swap: exchange `sibling_id` and `nephew_id`
/// across the internal edge between `parent_id` and `child_id`.
///
/// Returns a new tree with the swap applied.
fn nni_swap(
    tree: &PhyloTree,
    parent_id: NodeId,
    child_id: NodeId,
    sibling_id: NodeId,
    nephew_id: NodeId,
) -> Result<PhyloTree> {
    // Clone all nodes.
    let mut nodes: Vec<Node> = (0..tree.node_count())
        .map(|id| tree.get_node(id).unwrap().clone())
        .collect();

    // Remove nephew from child's children, add sibling.
    nodes[child_id].children.retain(|&c| c != nephew_id);
    nodes[child_id].children.push(sibling_id);

    // Remove sibling from parent's children, add nephew.
    nodes[parent_id].children.retain(|&c| c != sibling_id);
    nodes[parent_id].children.push(nephew_id);

    // Update parent pointers.
    nodes[sibling_id].parent = Some(child_id);
    nodes[nephew_id].parent = Some(parent_id);

    // Swap branch lengths of sibling and nephew to maintain consistency.
    // Actually, in NNI the branch lengths of the moved subtrees stay with
    // the subtrees, so we don't swap branch lengths -- the subtree root
    // keeps its own branch length.

    PhyloTree::from_nodes(nodes, tree.root())
}

/// Compute the log-likelihood of a tree under a general substitution model
/// with optional discrete gamma rate heterogeneity.
///
/// This is a generalized version of [`tree_likelihood`] that accepts any
/// substitution model as a closure and custom root frequencies.
///
/// When `gamma` is provided, the site likelihood is averaged over the
/// discrete gamma rate categories.
pub fn tree_likelihood_gtr(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: &dyn Fn(f64) -> [[f64; 4]; 4],
    freqs: &[f64; 4],
    gamma: Option<&crate::models::GammaRates>,
) -> Result<f64> {
    let leaves = tree.leaves();
    let n_leaves = leaves.len();

    if sequences.len() != n_leaves {
        return Err(CyaneaError::InvalidInput(format!(
            "expected {} sequences for {} leaves, got {}",
            n_leaves, n_leaves, sequences.len()
        )));
    }
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
                i, seq.len(), seq_len
            )));
        }
    }

    // Build leaf → sequence index mapping (sorted by leaf name).
    let mut leaf_ids_sorted: Vec<NodeId> = leaves.clone();
    leaf_ids_sorted.sort_by(|&a, &b| {
        let name_a = tree
            .get_node(a)
            .and_then(|n| n.name.as_ref())
            .cloned()
            .unwrap_or_default();
        let name_b = tree
            .get_node(b)
            .and_then(|n| n.name.as_ref())
            .cloned()
            .unwrap_or_default();
        name_a.cmp(&name_b)
    });

    let n_nodes = tree.node_count();
    let mut leaf_seq_index = vec![0usize; n_nodes];
    for (seq_idx, &node_id) in leaf_ids_sorted.iter().enumerate() {
        leaf_seq_index[node_id] = seq_idx;
    }

    // Get rate categories.
    let rates = match gamma {
        Some(g) => g.category_rates(),
        None => vec![1.0],
    };
    let n_cats = rates.len();

    // Precompute transition probability matrices for each (node, rate_category).
    let mut trans_probs: Vec<Vec<Option<[[f64; 4]; 4]>>> = vec![vec![None; n_cats]; n_nodes];
    for id in 0..n_nodes {
        if let Some(node) = tree.get_node(id) {
            if let Some(bl) = node.branch_length {
                for (c, &rate) in rates.iter().enumerate() {
                    trans_probs[id][c] = Some(model(bl * rate));
                }
            }
        }
    }

    let mut total_log_likelihood = 0.0;
    let mut partials = vec![[0.0f64; NUM_STATES]; n_nodes];
    let postorder: Vec<NodeId> = tree.iter_postorder().collect();

    for site in 0..seq_len {
        let mut site_likelihood = 0.0;

        for (cat, &_rate) in rates.iter().enumerate() {
            // Initialize partials.
            for p in partials.iter_mut() {
                *p = [0.0; NUM_STATES];
            }

            // Set leaf partials.
            for &leaf_id in &leaves {
                let seq_idx = leaf_seq_index[leaf_id];
                let base = sequences[seq_idx][site];
                if let Some(state) = nucleotide_index(base) {
                    partials[leaf_id][state] = 1.0;
                } else {
                    partials[leaf_id] = [1.0; NUM_STATES];
                }
            }

            // Post-order traversal.
            for &id in &postorder {
                let node = tree.get_node(id).unwrap();
                if node.is_leaf() {
                    continue;
                }

                let mut node_partial = [1.0f64; NUM_STATES];

                for &child_id in &node.children {
                    let prob_matrix = match &trans_probs[child_id][cat] {
                        Some(m) => m,
                        None => {
                            // No branch length: use small default.
                            if trans_probs[child_id][cat].is_none() {
                                trans_probs[child_id][cat] = Some(model(1e-6 * rates[cat]));
                            }
                            trans_probs[child_id][cat].as_ref().unwrap()
                        }
                    };
                    let child_partial = &partials[child_id];

                    for s in 0..NUM_STATES {
                        let mut sum = 0.0;
                        for t in 0..NUM_STATES {
                            sum += prob_matrix[s][t] * child_partial[t];
                        }
                        node_partial[s] *= sum;
                    }
                }

                partials[id] = node_partial;
            }

            // Root likelihood for this category.
            let root_id = tree.root();
            let root_partial = &partials[root_id];
            let mut cat_likelihood = 0.0;
            for s in 0..NUM_STATES {
                cat_likelihood += freqs[s] * root_partial[s];
            }
            site_likelihood += cat_likelihood / n_cats as f64;
        }

        if site_likelihood <= 0.0 {
            return Err(CyaneaError::InvalidInput(format!(
                "zero or negative site likelihood at site {}",
                site
            )));
        }

        total_log_likelihood += site_likelihood.ln();
    }

    Ok(total_log_likelihood)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::jc69_probability;

    /// Helper: build a simple 4-leaf tree with branch lengths and
    /// return aligned sequences that match the tree topology.
    fn test_tree_and_sequences() -> (PhyloTree, Vec<Vec<u8>>) {
        // Tree: ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();

        // Sequences (sorted by leaf name: A, B, C, D).
        // A and B are similar, C and D are similar, but differ from A/B.
        let seqs = vec![
            b"ACGTACGTACGT".to_vec(), // A
            b"ACGTACGTACGT".to_vec(), // B
            b"TGCATGCATGCA".to_vec(), // C
            b"TGCATGCATGCA".to_vec(), // D
        ];

        (tree, seqs)
    }

    #[test]
    fn likelihood_is_negative() {
        let (tree, seqs) = test_tree_and_sequences();
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let ll = tree_likelihood(&tree, &seq_refs, jc69_probability).unwrap();
        assert!(
            ll < 0.0,
            "log-likelihood should be negative, got {}",
            ll
        );
    }

    #[test]
    fn likelihood_identical_sequences() {
        // When all sequences are identical, likelihood should be higher
        // (closer to zero) than when sequences differ.
        let tree =
            PhyloTree::from_newick("((A:0.01,B:0.01):0.01,(C:0.01,D:0.01):0.01);")
                .unwrap();

        let identical: Vec<Vec<u8>> = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
        ];
        let different: Vec<Vec<u8>> = vec![
            b"AAAAAAAA".to_vec(),
            b"CCCCCCCC".to_vec(),
            b"GGGGGGGG".to_vec(),
            b"TTTTTTTT".to_vec(),
        ];

        let refs_id: Vec<&[u8]> = identical.iter().map(|s| s.as_slice()).collect();
        let refs_diff: Vec<&[u8]> = different.iter().map(|s| s.as_slice()).collect();

        let ll_id = tree_likelihood(&tree, &refs_id, jc69_probability).unwrap();
        let ll_diff = tree_likelihood(&tree, &refs_diff, jc69_probability).unwrap();

        assert!(
            ll_id > ll_diff,
            "identical sequences should have higher likelihood: {} vs {}",
            ll_id,
            ll_diff
        );
    }

    #[test]
    fn likelihood_wrong_sequence_count() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT"]; // 2 seqs for 3 leaves
        assert!(tree_likelihood(&tree, &seqs, jc69_probability).is_err());
    }

    #[test]
    fn likelihood_empty_sequences() {
        let tree = PhyloTree::from_newick("(A:0.1,B:0.1);").unwrap();
        let seqs: Vec<&[u8]> = vec![b"", b""];
        assert!(tree_likelihood(&tree, &seqs, jc69_probability).is_err());
    }

    #[test]
    fn likelihood_length_mismatch() {
        let tree = PhyloTree::from_newick("(A:0.1,B:0.1);").unwrap();
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"AC"];
        assert!(tree_likelihood(&tree, &seqs, jc69_probability).is_err());
    }

    #[test]
    fn nni_preserves_or_improves_likelihood() {
        // Start with a deliberately suboptimal tree topology.
        // Sequences: A and C are similar, B and D are similar.
        // But the tree groups A with B and C with D -- suboptimal.
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();

        let seqs: Vec<Vec<u8>> = vec![
            b"AAACCCAAACCC".to_vec(), // A
            b"CCCAAACCCAAA".to_vec(), // B - differs from A
            b"AAACCCAAACCC".to_vec(), // C - same as A
            b"CCCAAACCCAAA".to_vec(), // D - same as B
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let ll_before = tree_likelihood(&tree, &refs, jc69_probability).unwrap();
        let improved = nni_search(&tree, &refs, jc69_probability).unwrap();
        let ll_after = tree_likelihood(&improved, &refs, jc69_probability).unwrap();

        assert!(
            ll_after >= ll_before - 1e-10,
            "NNI should not decrease likelihood: {} -> {}",
            ll_before,
            ll_after
        );
    }

    #[test]
    fn nni_on_optimal_tree_is_stable() {
        // If the tree already matches the data, NNI should not change it.
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();

        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGT".to_vec(), // A
            b"ACGTACGTACGT".to_vec(), // B - same as A
            b"TGCATGCATGCA".to_vec(), // C
            b"TGCATGCATGCA".to_vec(), // D - same as C
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let ll_before = tree_likelihood(&tree, &refs, jc69_probability).unwrap();
        let result = nni_search(&tree, &refs, jc69_probability).unwrap();
        let ll_after = tree_likelihood(&result, &refs, jc69_probability).unwrap();

        assert!(
            ll_after >= ll_before - 1e-10,
            "NNI on optimal tree should maintain likelihood: {} -> {}",
            ll_before,
            ll_after
        );
    }

    #[test]
    fn tree_likelihood_gtr_finite() {
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTACGT".to_vec(),
            b"TGCATGCATGCA".to_vec(),
            b"TGCATGCATGCA".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let params =
            crate::models::GtrParams::new([1.0, 2.0, 1.0, 1.0, 2.0, 1.0], [0.3, 0.2, 0.2, 0.3])
                .unwrap();
        let prob_fn = params.probability_fn();
        let gamma = crate::models::GammaRates::new(0.5, 4).unwrap();

        let ll = tree_likelihood_gtr(&tree, &refs, &prob_fn, &params.freqs, Some(&gamma)).unwrap();
        assert!(ll.is_finite(), "log-likelihood should be finite, got {}", ll);
        assert!(ll < 0.0, "log-likelihood should be negative, got {}", ll);
    }

    #[test]
    fn nni_five_leaves() {
        let tree = PhyloTree::from_newick(
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,(D:0.1,E:0.1):0.1);",
        )
        .unwrap();

        let seqs: Vec<Vec<u8>> = vec![
            b"AACCGGTTAACCGGTT".to_vec(), // A
            b"AACCGGTTAACCGGTT".to_vec(), // B
            b"AACCGGTTAACCGGTT".to_vec(), // C
            b"TTGGCCAATTGGCCAA".to_vec(), // D
            b"TTGGCCAATTGGCCAA".to_vec(), // E
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let ll_before = tree_likelihood(&tree, &refs, jc69_probability).unwrap();
        let result = nni_search(&tree, &refs, jc69_probability).unwrap();
        let ll_after = tree_likelihood(&result, &refs, jc69_probability).unwrap();

        assert!(ll_after >= ll_before - 1e-10);
        assert_eq!(result.leaf_count(), 5);
    }
}
