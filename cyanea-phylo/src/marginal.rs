//! Marginal ancestral sequence reconstruction via maximum likelihood.
//!
//! Implements a two-pass algorithm to compute posterior probabilities of
//! each nucleotide state at every internal node, given the tree topology,
//! branch lengths, aligned sequences, and a substitution model.
//!
//! The algorithm:
//! 1. **Upward (postorder):** Felsenstein pruning → partial likelihoods from below
//! 2. **Downward (preorder):** propagate information from above each node
//! 3. **Combine:** posterior ∝ upward × downward, normalized

use crate::models::{nucleotide_index, NUM_STATES};
use crate::tree::{NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Posterior probability distribution over 4 nucleotide states at one node/site.
#[derive(Debug, Clone)]
pub struct MarginalPosterior {
    /// Probabilities for A, C, G, T (indices 0–3).
    pub probs: [f64; 4],
}

/// Result of marginal ancestral reconstruction.
#[derive(Debug, Clone)]
pub struct MarginalReconstruction {
    /// Posterior distributions: `posteriors[node_id][site]`.
    pub posteriors: Vec<Vec<MarginalPosterior>>,
    /// MAP (maximum a posteriori) state at each node/site: `map_states[node_id][site]`.
    pub map_states: Vec<Vec<u8>>,
}

/// Perform marginal ancestral reconstruction.
///
/// Computes posterior probabilities of each nucleotide state at every node
/// (leaves and internals) for each alignment column.
///
/// # Arguments
///
/// * `tree` - A rooted phylogenetic tree with branch lengths.
/// * `sequences` - Aligned sequences, one per leaf, sorted by leaf name.
/// * `model` - Substitution model: closure mapping branch length → 4×4 transition matrix.
/// * `freqs` - Root equilibrium frequencies [π_A, π_C, π_G, π_T].
pub fn marginal_reconstruct(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: &dyn Fn(f64) -> [[f64; 4]; 4],
    freqs: &[f64; 4],
) -> Result<MarginalReconstruction> {
    let leaves = tree.leaves();
    let n_leaves = leaves.len();

    if sequences.len() != n_leaves {
        return Err(CyaneaError::InvalidInput(format!(
            "expected {} sequences, got {}",
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

    // Build leaf → sequence index mapping (sorted by name).
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

    // Precompute transition probability matrices.
    let mut trans_probs: Vec<[[f64; 4]; 4]> = vec![[[0.0; 4]; 4]; n_nodes];
    for id in 0..n_nodes {
        if let Some(node) = tree.get_node(id) {
            let bl = node.branch_length.unwrap_or(1e-6);
            trans_probs[id] = model(bl);
        }
    }

    let postorder: Vec<NodeId> = tree.iter_postorder().collect();
    let preorder: Vec<NodeId> = tree.iter_preorder().collect();
    let leaf_set: Vec<bool> = (0..n_nodes)
        .map(|id| tree.get_node(id).map_or(false, |n| n.is_leaf()))
        .collect();

    let mut all_posteriors: Vec<Vec<MarginalPosterior>> = (0..n_nodes)
        .map(|_| Vec::with_capacity(seq_len))
        .collect();
    let mut all_map_states: Vec<Vec<u8>> = (0..n_nodes)
        .map(|_| Vec::with_capacity(seq_len))
        .collect();

    // Process each site.
    let mut l_down = vec![[0.0f64; NUM_STATES]; n_nodes];
    let mut l_up = vec![[0.0f64; NUM_STATES]; n_nodes];

    for site in 0..seq_len {
        // === Pass 1: Upward (postorder) ===
        for p in l_down.iter_mut() {
            *p = [0.0; NUM_STATES];
        }

        // Set leaf partials.
        for &leaf_id in &leaves {
            let seq_idx = leaf_seq_index[leaf_id];
            let base = sequences[seq_idx][site];
            if let Some(state) = nucleotide_index(base) {
                l_down[leaf_id][state] = 1.0;
            } else {
                l_down[leaf_id] = [1.0; NUM_STATES];
            }
        }

        for &id in &postorder {
            if leaf_set[id] {
                continue;
            }
            let node = tree.get_node(id).unwrap();
            let mut partial = [1.0f64; NUM_STATES];
            for &child_id in &node.children {
                let pm = &trans_probs[child_id];
                for s in 0..NUM_STATES {
                    let mut sum = 0.0;
                    for t in 0..NUM_STATES {
                        sum += pm[s][t] * l_down[child_id][t];
                    }
                    partial[s] *= sum;
                }
            }
            l_down[id] = partial;
        }

        // === Pass 2: Downward (preorder) ===
        for p in l_up.iter_mut() {
            *p = [0.0; NUM_STATES];
        }

        let root_id = tree.root();
        for s in 0..NUM_STATES {
            l_up[root_id][s] = freqs[s];
        }

        for &id in &preorder {
            if id == root_id {
                continue;
            }
            let node = tree.get_node(id).unwrap();
            let parent_id = node.parent.unwrap();
            let parent_node = tree.get_node(parent_id).unwrap();
            let pm = &trans_probs[id]; // P matrix for edge id→parent

            // L_up[v][s] = Σ_t P_v(s|t) · L_up[parent][t] · Π_{sibling w} (Σ_r P_w(r|t) · L_down[w][r])
            for s in 0..NUM_STATES {
                let mut sum = 0.0;
                for t in 0..NUM_STATES {
                    let mut contrib = pm[t][s] * l_up[parent_id][t];
                    // Multiply by sibling contributions.
                    for &sib_id in &parent_node.children {
                        if sib_id == id {
                            continue;
                        }
                        let sib_pm = &trans_probs[sib_id];
                        let mut sib_sum = 0.0;
                        for r in 0..NUM_STATES {
                            sib_sum += sib_pm[t][r] * l_down[sib_id][r];
                        }
                        contrib *= sib_sum;
                    }
                    sum += contrib;
                }
                l_up[id][s] = sum;
            }
        }

        // === Pass 3: Combine ===
        for id in 0..n_nodes {
            let mut posterior = [0.0f64; NUM_STATES];
            for s in 0..NUM_STATES {
                posterior[s] = l_down[id][s] * l_up[id][s];
            }
            // Normalize.
            let total: f64 = posterior.iter().sum();
            if total > 0.0 {
                for s in 0..NUM_STATES {
                    posterior[s] /= total;
                }
            }

            // MAP state.
            let map = (0..NUM_STATES)
                .max_by(|&a, &b| posterior[a].partial_cmp(&posterior[b]).unwrap())
                .unwrap();
            let map_byte = match map {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };

            all_posteriors[id].push(MarginalPosterior { probs: posterior });
            all_map_states[id].push(map_byte);
        }
    }

    Ok(MarginalReconstruction {
        posteriors: all_posteriors,
        map_states: all_map_states,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::jc69_probability;

    fn get_leaf_id(tree: &PhyloTree, name: &str) -> NodeId {
        for id in tree.leaves() {
            if tree.get_node(id).unwrap().name.as_deref() == Some(name) {
                return id;
            }
        }
        panic!("leaf {} not found", name);
    }

    #[test]
    fn marginal_identical_sequences() {
        // If all leaves have the same base, all internal nodes should reconstruct it.
        let tree = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![b"AAAA".to_vec(), b"AAAA".to_vec(), b"AAAA".to_vec()];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let freqs = [0.25; 4];

        let result =
            marginal_reconstruct(&tree, &refs, &(jc69_probability as fn(f64) -> _), &freqs)
                .unwrap();

        // All MAP states should be 'A'.
        for id in 0..tree.node_count() {
            for site in 0..4 {
                assert_eq!(
                    result.map_states[id][site], b'A',
                    "node {} site {} should be A",
                    id, site
                );
            }
        }
    }

    #[test]
    fn marginal_posteriors_sum_to_one() {
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"TGCATGCA".to_vec(),
            b"TGCATGCA".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let freqs = [0.25; 4];

        let result =
            marginal_reconstruct(&tree, &refs, &(jc69_probability as fn(f64) -> _), &freqs)
                .unwrap();

        for id in 0..tree.node_count() {
            for site in 0..8 {
                let sum: f64 = result.posteriors[id][site].probs.iter().sum();
                assert!(
                    (sum - 1.0).abs() < 1e-8,
                    "posteriors at node {} site {} sum to {} (expected 1.0)",
                    id,
                    site,
                    sum
                );
            }
        }
    }

    #[test]
    fn marginal_map_agrees_with_fitch() {
        // For strongly-supported sites, marginal MAP should agree with Fitch parsimony.
        let tree = PhyloTree::from_newick("((A:0.01,B:0.01):0.01,C:0.01);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![b"AA".to_vec(), b"AA".to_vec(), b"AT".to_vec()];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let freqs = [0.25; 4];

        let result =
            marginal_reconstruct(&tree, &refs, &(jc69_probability as fn(f64) -> _), &freqs)
                .unwrap();

        // Site 0: all A → root should be A.
        let root_id = tree.root();
        assert_eq!(result.map_states[root_id][0], b'A');

        // Internal node (parent of A,B) at site 0 should also be A.
        let a_id = get_leaf_id(&tree, "A");
        let ab_parent = tree.get_node(a_id).unwrap().parent.unwrap();
        assert_eq!(result.map_states[ab_parent][0], b'A');
    }

    #[test]
    fn marginal_correct_dimensions() {
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let freqs = [0.25; 4];

        let result =
            marginal_reconstruct(&tree, &refs, &(jc69_probability as fn(f64) -> _), &freqs)
                .unwrap();

        assert_eq!(result.posteriors.len(), tree.node_count());
        assert_eq!(result.map_states.len(), tree.node_count());
        for id in 0..tree.node_count() {
            assert_eq!(result.posteriors[id].len(), 4);
            assert_eq!(result.map_states[id].len(), 4);
        }
    }
}
