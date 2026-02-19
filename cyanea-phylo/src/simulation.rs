//! Sequence evolution and coalescent tree simulation.
//!
//! Provides functions for simulating DNA sequence evolution along a phylogenetic
//! tree under a given substitution model, and for generating random coalescent
//! trees (constant population size or exponential growth).

use cyanea_core::{CyaneaError, Result};

use crate::subst_model::SubstitutionModel;
use crate::tree::{Node, NodeId, PhyloTree};

/// Nucleotide bytes indexed by state (A=0, C=1, G=2, T=3).
const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Result of simulating sequence evolution along a phylogenetic tree.
#[derive(Debug, Clone)]
pub struct SimulatedAlignment {
    /// Leaf taxon names, in the order they appear in the tree.
    pub names: Vec<String>,
    /// Evolved sequences for each leaf (same order as `names`).
    pub sequences: Vec<Vec<u8>>,
    /// Total number of substitution events across all branches and sites.
    pub n_substitutions: usize,
}

/// Simple xorshift64 pseudo-random number generator.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Ensure state is never zero (xorshift requires nonzero state).
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        self.next_u64() as f64 / u64::MAX as f64
    }
}

/// Sample a discrete state from cumulative frequency distribution.
///
/// `cumulative` is a slice where `cumulative[i]` is the cumulative probability
/// up to and including state `i`. The function returns the first index where
/// the uniform draw `u` falls below the cumulative probability.
fn sample_state(cumulative: &[f64], u: f64) -> usize {
    for (i, &c) in cumulative.iter().enumerate() {
        if u <= c {
            return i;
        }
    }
    // Fallback to last state (handles floating-point edge cases).
    cumulative.len() - 1
}

/// Build a cumulative distribution from a probability row.
fn cumulative_from_row(row: &[f64]) -> Vec<f64> {
    let mut cum = Vec::with_capacity(row.len());
    let mut acc = 0.0;
    for &p in row {
        acc += p;
        cum.push(acc);
    }
    // Ensure the last entry is exactly 1.0 to avoid rounding issues.
    if let Some(last) = cum.last_mut() {
        *last = 1.0;
    }
    cum
}

/// Simulate sequence evolution along a phylogenetic tree.
///
/// Starting from a root sequence drawn from the model's equilibrium
/// frequencies, this function evolves sequences down every branch of `tree`
/// according to the transition probability matrices of `model`. The returned
/// alignment contains one sequence per leaf, labelled with the leaf's name
/// (or `"leaf_{id}"` when unnamed).
///
/// # Errors
///
/// Returns an error if `seq_length` is zero or if the tree has no leaves.
///
/// # Examples
///
/// ```
/// use cyanea_phylo::tree::PhyloTree;
/// use cyanea_phylo::subst_model::Jc69Model;
/// use cyanea_phylo::simulation::simulate_evolution;
///
/// let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
/// let model = Jc69Model::new();
/// let aln = simulate_evolution(&tree, &model, 100, 42).unwrap();
/// assert_eq!(aln.names.len(), 3);
/// assert!(aln.sequences.iter().all(|s| s.len() == 100));
/// ```
pub fn simulate_evolution(
    tree: &PhyloTree,
    model: &dyn SubstitutionModel,
    seq_length: usize,
    seed: u64,
) -> Result<SimulatedAlignment> {
    // --- validation ---
    if seq_length == 0 {
        return Err(CyaneaError::InvalidInput(
            "seq_length must be > 0".into(),
        ));
    }
    let leaves = tree.leaves();
    if leaves.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "tree has no leaf nodes".into(),
        ));
    }

    let n_states = model.n_states();
    let freqs = model.frequencies();

    let mut rng = Xorshift64::new(seed);

    // --- build cumulative equilibrium frequencies ---
    let eq_cum = cumulative_from_row(freqs);

    // --- generate root sequence from equilibrium distribution ---
    let n_nodes = tree.node_count();
    let mut node_seqs: Vec<Option<Vec<usize>>> = vec![None; n_nodes];

    let root_id = tree.root();
    let root_seq: Vec<usize> = (0..seq_length)
        .map(|_| sample_state(&eq_cum, rng.next_f64()))
        .collect();
    node_seqs[root_id] = Some(root_seq);

    // --- evolve along tree in pre-order ---
    let mut n_substitutions: usize = 0;

    for node_id in tree.iter_preorder() {
        let node = tree.get_node(node_id).unwrap();
        if node.is_root() && node_id == root_id {
            // Root already has its sequence.
            continue;
        }

        // This node must have a parent whose sequence is already assigned
        // (pre-order guarantees parent is visited first).
        let parent_id = match node.parent {
            Some(p) => p,
            None => continue, // Should not happen in a well-formed tree.
        };

        let branch_len = node.branch_length.unwrap_or(0.0);
        let trans_prob = model.transition_probs(branch_len);

        // Precompute cumulative distributions for each parent state.
        let cum_rows: Vec<Vec<f64>> = (0..n_states)
            .map(|i| cumulative_from_row(&trans_prob[i]))
            .collect();

        // Clone the parent sequence (we need immutable access while building child).
        let parent_seq = node_seqs[parent_id]
            .as_ref()
            .expect("parent sequence must exist in preorder traversal")
            .clone();

        let mut child_seq = Vec::with_capacity(seq_length);
        for &parent_state in &parent_seq {
            let child_state = sample_state(&cum_rows[parent_state], rng.next_f64());
            if child_state != parent_state {
                n_substitutions += 1;
            }
            child_seq.push(child_state);
        }

        node_seqs[node_id] = Some(child_seq);
    }

    // --- collect leaf sequences ---
    let mut names = Vec::with_capacity(leaves.len());
    let mut sequences = Vec::with_capacity(leaves.len());

    for &leaf_id in &leaves {
        let node = tree.get_node(leaf_id).unwrap();
        let name = node
            .name
            .clone()
            .unwrap_or_else(|| format!("leaf_{}", leaf_id));
        names.push(name);

        let state_seq = node_seqs[leaf_id]
            .as_ref()
            .expect("leaf must have a sequence after evolution");
        let nuc_seq: Vec<u8> = state_seq
            .iter()
            .map(|&s| {
                if s < NUCLEOTIDES.len() {
                    NUCLEOTIDES[s]
                } else {
                    b'N' // Fallback for non-nucleotide models.
                }
            })
            .collect();
        sequences.push(nuc_seq);
    }

    Ok(SimulatedAlignment {
        names,
        sequences,
        n_substitutions,
    })
}

/// Simulate a coalescent (Kingman) tree for a constant-size population.
///
/// Generates a random genealogy for `n_samples` lineages under the standard
/// neutral coalescent with effective population size `pop_size`. Coalescence
/// times are exponentially distributed with rate `k*(k-1) / (4*pop_size)` when
/// `k` lineages remain.
///
/// # Errors
///
/// Returns an error if `n_samples < 2` or `pop_size <= 0`.
///
/// # Examples
///
/// ```
/// use cyanea_phylo::simulation::simulate_coalescent;
///
/// let tree = simulate_coalescent(10, 1000.0, 42).unwrap();
/// assert_eq!(tree.leaf_count(), 10);
/// ```
pub fn simulate_coalescent(
    n_samples: usize,
    pop_size: f64,
    seed: u64,
) -> Result<PhyloTree> {
    if n_samples < 2 {
        return Err(CyaneaError::InvalidInput(
            "n_samples must be >= 2".into(),
        ));
    }
    if pop_size <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "pop_size must be > 0".into(),
        ));
    }

    let mut rng = Xorshift64::new(seed);

    // Create leaf nodes. Each leaf starts at time 0.
    let mut nodes: Vec<Node> = Vec::new();
    let mut lineages: Vec<(NodeId, f64)> = Vec::new(); // (node_id, current_time)

    for i in 0..n_samples {
        nodes.push(Node {
            id: i,
            parent: None,
            children: Vec::new(),
            branch_length: None,
            name: Some(format!("tip_{}", i)),
        });
        lineages.push((i, 0.0));
    }

    let mut current_time = 0.0;

    while lineages.len() > 1 {
        let k = lineages.len() as f64;
        let rate = k * (k - 1.0) / (4.0 * pop_size);

        // Sample exponential waiting time: -ln(U) / rate.
        let u = rng.next_f64();
        // Guard against u == 0 (ln(0) = -inf).
        let u_safe = if u < 1e-300 { 1e-300 } else { u };
        let wait = -u_safe.ln() / rate;
        current_time += wait;

        // Pick two random lineages to coalesce.
        let n_lin = lineages.len();
        let idx1 = (rng.next_u64() as usize) % n_lin;
        let mut idx2 = (rng.next_u64() as usize) % (n_lin - 1);
        if idx2 >= idx1 {
            idx2 += 1;
        }

        let (child1_id, child1_time) = lineages[idx1];
        let (child2_id, child2_time) = lineages[idx2];

        // Create internal (coalescent) node.
        let new_id = nodes.len();
        nodes.push(Node {
            id: new_id,
            parent: None,
            children: vec![child1_id, child2_id],
            branch_length: None,
            name: None,
        });

        // Set parent and branch lengths for children.
        nodes[child1_id].parent = Some(new_id);
        nodes[child1_id].branch_length = Some(current_time - child1_time);

        nodes[child2_id].parent = Some(new_id);
        nodes[child2_id].branch_length = Some(current_time - child2_time);

        // Remove the two coalesced lineages and add the new one.
        // Remove higher index first to avoid shifting issues.
        let (remove_first, remove_second) = if idx1 > idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };
        lineages.swap_remove(remove_first);
        lineages.swap_remove(remove_second);
        lineages.push((new_id, current_time));
    }

    let root_id = lineages[0].0;
    PhyloTree::from_nodes(nodes, root_id)
}

/// Simulate a coalescent tree with exponential population growth.
///
/// The population size looking backward in time is
/// `N(t) = current_pop * exp(-growth_rate * t)`, where `t` is time into the
/// past. When `growth_rate > 0` the population was smaller in the past,
/// leading to faster coalescence and shorter trees than under a constant-size
/// model.
///
/// The waiting time for the next coalescence event is computed exactly by
/// inverting the CDF of the integrated coalescence rate:
///
/// ```text
/// t = ln(1 + 4 * growth_rate * current_pop * (-ln(U)) / (k*(k-1))) / growth_rate
/// ```
///
/// where `U ~ Uniform(0, 1)` and `k` is the number of active lineages.
///
/// # Errors
///
/// Returns an error if `n_samples < 2`, `current_pop <= 0`, or
/// `growth_rate < 0`.
///
/// # Examples
///
/// ```
/// use cyanea_phylo::simulation::simulate_coalescent_growth;
///
/// let tree = simulate_coalescent_growth(10, 1000.0, 0.01, 42).unwrap();
/// assert_eq!(tree.leaf_count(), 10);
/// ```
pub fn simulate_coalescent_growth(
    n_samples: usize,
    current_pop: f64,
    growth_rate: f64,
    seed: u64,
) -> Result<PhyloTree> {
    if n_samples < 2 {
        return Err(CyaneaError::InvalidInput(
            "n_samples must be >= 2".into(),
        ));
    }
    if current_pop <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "current_pop must be > 0".into(),
        ));
    }
    if growth_rate < 0.0 {
        return Err(CyaneaError::InvalidInput(
            "growth_rate must be >= 0".into(),
        ));
    }

    // If growth_rate is effectively zero, delegate to constant-size coalescent.
    if growth_rate < 1e-15 {
        return simulate_coalescent(n_samples, current_pop, seed);
    }

    let mut rng = Xorshift64::new(seed);

    // Create leaf nodes.
    let mut nodes: Vec<Node> = Vec::new();
    let mut lineages: Vec<(NodeId, f64)> = Vec::new();

    for i in 0..n_samples {
        nodes.push(Node {
            id: i,
            parent: None,
            children: Vec::new(),
            branch_length: None,
            name: Some(format!("tip_{}", i)),
        });
        lineages.push((i, 0.0));
    }

    let mut current_time = 0.0;

    while lineages.len() > 1 {
        let k = lineages.len() as f64;
        let coal_choose = k * (k - 1.0);

        // For exponential growth backward in time: N(t) = current_pop * exp(-growth_rate * t)
        // Rate at time s (measured from current_time): lambda(s) = coal_choose / (4 * N(current_time + s))
        //   = coal_choose / (4 * current_pop * exp(-growth_rate * (current_time + s)))
        //   = coal_choose * exp(growth_rate * (current_time + s)) / (4 * current_pop)
        //
        // Integrated rate from 0 to t:
        //   Lambda(t) = coal_choose / (4 * current_pop) * integral_0^t exp(growth_rate * (current_time + s)) ds
        //             = coal_choose * exp(growth_rate * current_time) / (4 * current_pop * growth_rate)
        //               * (exp(growth_rate * t) - 1)
        //
        // CDF: P(T <= t) = 1 - exp(-Lambda(t))
        // Inversion: -ln(U) = Lambda(t)
        //   => exp(growth_rate * t) - 1 = 4 * current_pop * growth_rate * (-ln(U))
        //                                  / (coal_choose * exp(growth_rate * current_time))
        //   => t = ln(1 + 4 * current_pop * growth_rate * (-ln(U))
        //            / (coal_choose * exp(growth_rate * current_time))) / growth_rate

        let u = rng.next_f64();
        let u_safe = if u < 1e-300 { 1e-300 } else { u };
        let neg_ln_u = -u_safe.ln();

        let exp_gt = (growth_rate * current_time).exp();
        let arg = 1.0 + 4.0 * current_pop * growth_rate * neg_ln_u / (coal_choose * exp_gt);

        // arg should always be > 1, but protect against numerical issues.
        let wait = if arg <= 0.0 {
            // Extremely rare edge case; fall back to a tiny step.
            1e-15
        } else {
            arg.ln() / growth_rate
        };

        current_time += wait;

        // Pick two random lineages.
        let n_lin = lineages.len();
        let idx1 = (rng.next_u64() as usize) % n_lin;
        let mut idx2 = (rng.next_u64() as usize) % (n_lin - 1);
        if idx2 >= idx1 {
            idx2 += 1;
        }

        let (child1_id, child1_time) = lineages[idx1];
        let (child2_id, child2_time) = lineages[idx2];

        let new_id = nodes.len();
        nodes.push(Node {
            id: new_id,
            parent: None,
            children: vec![child1_id, child2_id],
            branch_length: None,
            name: None,
        });

        nodes[child1_id].parent = Some(new_id);
        nodes[child1_id].branch_length = Some(current_time - child1_time);

        nodes[child2_id].parent = Some(new_id);
        nodes[child2_id].branch_length = Some(current_time - child2_time);

        let (remove_first, remove_second) = if idx1 > idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };
        lineages.swap_remove(remove_first);
        lineages.swap_remove(remove_second);
        lineages.push((new_id, current_time));
    }

    let root_id = lineages[0].0;
    PhyloTree::from_nodes(nodes, root_id)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::subst_model::Jc69Model;

    /// Helper: build a simple balanced 4-leaf tree.
    fn balanced_tree() -> PhyloTree {
        // ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);
        let mut tree = PhyloTree::new();
        let ab = tree
            .add_child(0, Some("AB".into()), Some(0.1))
            .unwrap();
        let cd = tree
            .add_child(0, Some("CD".into()), Some(0.1))
            .unwrap();
        tree.add_child(ab, Some("A".into()), Some(0.1)).unwrap();
        tree.add_child(ab, Some("B".into()), Some(0.1)).unwrap();
        tree.add_child(cd, Some("C".into()), Some(0.1)).unwrap();
        tree.add_child(cd, Some("D".into()), Some(0.1)).unwrap();
        tree
    }

    /// Helper: build a tree where all branch lengths are zero.
    fn zero_branch_tree() -> PhyloTree {
        let mut tree = PhyloTree::new();
        let ab = tree.add_child(0, None, Some(0.0)).unwrap();
        let cd = tree.add_child(0, None, Some(0.0)).unwrap();
        tree.add_child(ab, Some("A".into()), Some(0.0)).unwrap();
        tree.add_child(ab, Some("B".into()), Some(0.0)).unwrap();
        tree.add_child(cd, Some("C".into()), Some(0.0)).unwrap();
        tree.add_child(cd, Some("D".into()), Some(0.0)).unwrap();
        tree
    }

    /// Helper: build a tree with very long branches to force divergence.
    fn long_branch_tree() -> PhyloTree {
        let mut tree = PhyloTree::new();
        tree.add_child(0, Some("A".into()), Some(50.0)).unwrap();
        tree.add_child(0, Some("B".into()), Some(50.0)).unwrap();
        tree
    }

    // ---------------------------------------------------------------
    // Test 1: Correct leaf count matches tree leaves.
    // ---------------------------------------------------------------
    #[test]
    fn evolution_leaf_count_matches_tree() {
        let tree = balanced_tree();
        let model = Jc69Model::new();
        let aln = simulate_evolution(&tree, &model, 200, 12345).unwrap();
        assert_eq!(aln.names.len(), tree.leaf_count());
        assert_eq!(aln.sequences.len(), tree.leaf_count());
    }

    // ---------------------------------------------------------------
    // Test 2: Zero branch lengths produce identical sequences.
    // ---------------------------------------------------------------
    #[test]
    fn zero_branches_identical_sequences() {
        let tree = zero_branch_tree();
        let model = Jc69Model::new();
        let aln = simulate_evolution(&tree, &model, 500, 99).unwrap();

        // Under zero branch lengths, P(t=0) is the identity matrix, so every
        // leaf sequence must be identical to the root sequence (and therefore
        // identical to each other).
        let first = &aln.sequences[0];
        for seq in &aln.sequences[1..] {
            assert_eq!(first, seq, "all leaf sequences must be identical when branches are zero");
        }
        assert_eq!(aln.n_substitutions, 0);
    }

    // ---------------------------------------------------------------
    // Test 3: Correct sequence length.
    // ---------------------------------------------------------------
    #[test]
    fn evolution_correct_sequence_length() {
        let tree = balanced_tree();
        let model = Jc69Model::new();
        for &len in &[1, 10, 100, 1000] {
            let aln = simulate_evolution(&tree, &model, len, 42).unwrap();
            for seq in &aln.sequences {
                assert_eq!(seq.len(), len);
            }
        }
    }

    // ---------------------------------------------------------------
    // Test 4: Long branches cause divergence.
    // ---------------------------------------------------------------
    #[test]
    fn long_branches_cause_divergence() {
        let tree = long_branch_tree();
        let model = Jc69Model::new();
        let aln = simulate_evolution(&tree, &model, 1000, 7).unwrap();

        // With branch length 50, sequences should have converged to equilibrium
        // frequencies and differ at roughly 75% of sites. Check that they are
        // not identical.
        let diffs: usize = aln.sequences[0]
            .iter()
            .zip(aln.sequences[1].iter())
            .filter(|(a, b)| a != b)
            .count();
        assert!(
            diffs > 100,
            "expected significant divergence with long branches, got {} diffs out of 1000",
            diffs
        );
    }

    // ---------------------------------------------------------------
    // Test 5: Coalescent produces correct leaf count.
    // ---------------------------------------------------------------
    #[test]
    fn coalescent_correct_leaf_count() {
        for &n in &[2, 5, 10, 20] {
            let tree = simulate_coalescent(n, 1000.0, 42).unwrap();
            assert_eq!(
                tree.leaf_count(),
                n,
                "coalescent tree should have {} leaves",
                n
            );
        }
    }

    // ---------------------------------------------------------------
    // Test 6: Coalescent tree height scales with pop_size.
    // ---------------------------------------------------------------
    #[test]
    fn coalescent_height_scales_with_pop_size() {
        let n = 10;
        let seed = 42;

        // Average over several replicates for robustness.
        let mut total_height_small = 0.0;
        let mut total_height_large = 0.0;
        let replicates = 20;

        for rep in 0..replicates {
            let s = seed + rep as u64;
            let tree_small = simulate_coalescent(n, 100.0, s).unwrap();
            let tree_large = simulate_coalescent(n, 10000.0, s + 1000).unwrap();
            total_height_small += tree_small.total_branch_length();
            total_height_large += tree_large.total_branch_length();
        }

        // Larger population size should produce a taller (longer total branch
        // length) tree on average. The ratio of expected total branch lengths
        // is proportional to the ratio of population sizes.
        assert!(
            total_height_large > total_height_small,
            "larger pop_size should produce taller trees: small={}, large={}",
            total_height_small,
            total_height_large
        );
    }

    // ---------------------------------------------------------------
    // Test 7: Growth coalescent produces shorter trees than constant-size.
    // ---------------------------------------------------------------
    #[test]
    fn growth_coalescent_shorter_than_constant() {
        let n = 10;
        let pop = 5000.0;
        let replicates = 20;

        let mut total_constant = 0.0;
        let mut total_growth = 0.0;

        for rep in 0..replicates {
            let s = 100 + rep as u64;
            let tree_const = simulate_coalescent(n, pop, s).unwrap();
            let tree_growth =
                simulate_coalescent_growth(n, pop, 0.05, s).unwrap();
            total_constant += tree_const.total_branch_length();
            total_growth += tree_growth.total_branch_length();
        }

        // With positive growth rate, the ancestral population was smaller,
        // causing faster coalescence and shorter trees.
        assert!(
            total_growth < total_constant,
            "growth coalescent should produce shorter trees: growth={}, constant={}",
            total_growth,
            total_constant
        );
    }

    // ---------------------------------------------------------------
    // Test 8: All branch lengths are non-negative.
    // ---------------------------------------------------------------
    #[test]
    fn coalescent_branch_lengths_nonnegative() {
        for &n in &[3, 10, 25] {
            let tree = simulate_coalescent(n, 500.0, 77).unwrap();
            for node in tree.nodes() {
                if let Some(bl) = node.branch_length {
                    assert!(
                        bl >= 0.0,
                        "branch length must be >= 0, got {} for node {}",
                        bl,
                        node.id
                    );
                }
            }

            let tree_growth =
                simulate_coalescent_growth(n, 500.0, 0.01, 77).unwrap();
            for node in tree_growth.nodes() {
                if let Some(bl) = node.branch_length {
                    assert!(
                        bl >= 0.0,
                        "branch length must be >= 0, got {} for node {}",
                        bl,
                        node.id
                    );
                }
            }
        }
    }
}
