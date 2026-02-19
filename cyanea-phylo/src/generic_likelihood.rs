//! Generic N-state Felsenstein pruning algorithm.
//!
//! Works with any [`SubstitutionModel`], supporting both nucleotide (4-state)
//! and protein (20-state) models via dynamic dispatch.

use crate::models::{nucleotide_index, GammaRates};
use crate::protein_models::amino_acid_index;
use crate::subst_model::SubstitutionModel;
use crate::tree::{NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Compute the log-likelihood of a phylogenetic tree under any substitution model.
///
/// Uses Felsenstein's pruning algorithm with N-state partials determined by
/// the model's `n_states()`. Automatically detects DNA (4) or protein (20)
/// sequences based on the model.
pub fn generic_tree_likelihood(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: &dyn SubstitutionModel,
    gamma: Option<&GammaRates>,
) -> Result<f64> {
    let site_lls = site_likelihoods(tree, sequences, model, gamma)?;
    Ok(site_lls.iter().sum())
}

/// Compute per-site log-likelihoods under any substitution model.
///
/// Returns a vector of log-likelihoods, one per alignment column.
pub fn site_likelihoods(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: &dyn SubstitutionModel,
    gamma: Option<&GammaRates>,
) -> Result<Vec<f64>> {
    let n_states = model.n_states();
    let freqs = model.frequencies();
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
                i,
                seq.len(),
                seq_len
            )));
        }
    }

    // Map leaves to sequence indices (sorted by name).
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

    // Determine state indexing function.
    let state_fn: fn(u8) -> Option<usize> = if n_states == 4 {
        nucleotide_index
    } else {
        amino_acid_index
    };

    // Rate categories.
    let rates = match gamma {
        Some(g) => g.category_rates(),
        None => vec![1.0],
    };
    let n_cats = rates.len();

    // Precompute transition probability matrices for each (node, rate_category).
    let mut trans_probs: Vec<Vec<Option<Vec<Vec<f64>>>>> =
        vec![vec![None; n_cats]; n_nodes];
    for id in 0..n_nodes {
        if let Some(node) = tree.get_node(id) {
            if let Some(bl) = node.branch_length {
                for (c, &rate) in rates.iter().enumerate() {
                    trans_probs[id][c] = Some(model.transition_probs(bl * rate));
                }
            }
        }
    }

    let postorder: Vec<NodeId> = tree.iter_postorder().collect();
    let mut site_lls = Vec::with_capacity(seq_len);

    // Reusable partial likelihood buffers.
    let mut partials = vec![vec![0.0f64; n_states]; n_nodes];

    for site in 0..seq_len {
        let mut site_likelihood = 0.0;

        for (cat, &_rate) in rates.iter().enumerate() {
            // Reset partials.
            for p in partials.iter_mut() {
                for v in p.iter_mut() {
                    *v = 0.0;
                }
            }

            // Set leaf partials.
            for &leaf_id in &leaves {
                let seq_idx = leaf_seq_index[leaf_id];
                let base = sequences[seq_idx][site];
                if let Some(state) = state_fn(base) {
                    if state < n_states {
                        partials[leaf_id][state] = 1.0;
                    } else {
                        // Unknown state: equal probability.
                        for v in partials[leaf_id].iter_mut() {
                            *v = 1.0;
                        }
                    }
                } else {
                    // Ambiguous: equal probability for all states.
                    for v in partials[leaf_id].iter_mut() {
                        *v = 1.0;
                    }
                }
            }

            // Post-order: compute internal node partials.
            for &id in &postorder {
                let node = tree.get_node(id).unwrap();
                if node.is_leaf() {
                    continue;
                }

                // Initialize to 1.0 (multiplicative identity).
                for v in partials[id].iter_mut() {
                    *v = 1.0;
                }

                for &child_id in &node.children {
                    let prob_matrix = match &trans_probs[child_id][cat] {
                        Some(m) => m,
                        None => {
                            // No branch length: compute with small default.
                            trans_probs[child_id][cat] =
                                Some(model.transition_probs(1e-6 * rates[cat]));
                            trans_probs[child_id][cat].as_ref().unwrap()
                        }
                    };

                    for s in 0..n_states {
                        let mut sum = 0.0;
                        for t in 0..n_states {
                            sum += prob_matrix[s][t] * partials[child_id][t];
                        }
                        partials[id][s] *= sum;
                    }
                }
            }

            // Root likelihood for this category.
            let root_id = tree.root();
            let mut cat_likelihood = 0.0;
            for s in 0..n_states {
                cat_likelihood += freqs[s] * partials[root_id][s];
            }
            site_likelihood += cat_likelihood / n_cats as f64;
        }

        if site_likelihood <= 0.0 {
            return Err(CyaneaError::InvalidInput(format!(
                "zero or negative site likelihood at site {}",
                site
            )));
        }

        site_lls.push(site_likelihood.ln());
    }

    Ok(site_lls)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::jc69_probability;
    use crate::subst_model::Jc69Model;

    fn test_tree() -> PhyloTree {
        PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap()
    }

    #[test]
    fn generic_jc69_matches_original() {
        let tree = test_tree();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTACGT".to_vec(),
            b"TGCATGCATGCA".to_vec(),
            b"TGCATGCATGCA".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let model = Jc69Model::new();
        let generic_ll =
            generic_tree_likelihood(&tree, &refs, &model, None).unwrap();
        let original_ll =
            crate::likelihood::tree_likelihood(&tree, &refs, jc69_probability).unwrap();

        assert!(
            (generic_ll - original_ll).abs() < 1e-6,
            "generic {} vs original {}",
            generic_ll,
            original_ll
        );
    }

    #[test]
    fn generic_gtr_gamma_matches_original() {
        let tree = test_tree();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTACGT".to_vec(),
            b"TGCATGCATGCA".to_vec(),
            b"TGCATGCATGCA".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let params = crate::models::GtrParams::new(
            [1.0, 2.0, 1.0, 1.0, 2.0, 1.0],
            [0.3, 0.2, 0.2, 0.3],
        )
        .unwrap();
        let gamma = GammaRates::new(0.5, 4).unwrap();

        let gtr_model = crate::subst_model::GtrModel::new(params.clone());
        let generic_ll =
            generic_tree_likelihood(&tree, &refs, &gtr_model, Some(&gamma)).unwrap();

        let prob_fn = params.probability_fn();
        let original_ll = crate::likelihood::tree_likelihood_gtr(
            &tree,
            &refs,
            &prob_fn,
            &params.freqs,
            Some(&gamma),
        )
        .unwrap();

        assert!(
            (generic_ll - original_ll).abs() < 1e-4,
            "generic {} vs original {}",
            generic_ll,
            original_ll
        );
    }

    #[test]
    fn protein_likelihood_computable() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.1,(C:0.3,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACDEFGHIKL".to_vec(),
            b"ACDEFGHIKL".to_vec(),
            b"MNPQRSTVWY".to_vec(),
            b"MNPQRSTVWY".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let model = crate::protein_models::LgModel::new();
        let ll = generic_tree_likelihood(&tree, &refs, &model, None).unwrap();
        assert!(ll.is_finite(), "protein likelihood should be finite");
        assert!(ll < 0.0, "protein likelihood should be negative");
    }

    #[test]
    fn site_likelihoods_sum_matches_total() {
        let tree = test_tree();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTAC".to_vec(),
            b"ACGTAC".to_vec(),
            b"TGCATG".to_vec(),
            b"TGCATG".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let model = Jc69Model::new();
        let total = generic_tree_likelihood(&tree, &refs, &model, None).unwrap();
        let per_site = site_likelihoods(&tree, &refs, &model, None).unwrap();
        let sum: f64 = per_site.iter().sum();

        assert!(
            (total - sum).abs() < 1e-10,
            "total {} vs site sum {}",
            total,
            sum
        );
        assert_eq!(per_site.len(), 6);
    }

    #[test]
    fn empty_sequences_error() {
        let tree = PhyloTree::from_newick("(A:0.1,B:0.1);").unwrap();
        let seqs: Vec<&[u8]> = vec![b"", b""];
        let model = Jc69Model::new();
        assert!(generic_tree_likelihood(&tree, &seqs, &model, None).is_err());
    }

    #[test]
    fn mismatched_leaf_count_error() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT"]; // 2 seqs for 3 leaves
        let model = Jc69Model::new();
        assert!(generic_tree_likelihood(&tree, &seqs, &model, None).is_err());
    }

    #[test]
    fn ambiguous_bases_handled() {
        let tree = PhyloTree::from_newick("(A:0.1,B:0.1);").unwrap();
        let seqs: Vec<&[u8]> = vec![b"ACNGT", b"ACNGT"];
        let model = Jc69Model::new();
        let ll = generic_tree_likelihood(&tree, &seqs, &model, None).unwrap();
        assert!(ll.is_finite());
    }
}
