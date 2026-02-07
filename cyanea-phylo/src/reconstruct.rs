//! Ancestral sequence reconstruction via parsimony methods.
//!
//! Provides algorithms for inferring ancestral character states at internal
//! nodes of a phylogenetic tree:
//!
//! - **Fitch parsimony** — unweighted maximum parsimony. Assigns states by
//!   a two-pass algorithm: bottom-up intersection/union, then top-down selection.
//!   Optimal for minimizing the total number of changes on the tree.
//!
//! - **Sankoff parsimony** — weighted (step-matrix) parsimony. Each state
//!   transition can have a different cost. The bottom-up pass computes the
//!   minimum cost for each state at each node; the top-down pass selects
//!   optimal states.

use std::collections::HashSet;

use cyanea_core::{CyaneaError, Result};

use crate::tree::{NodeId, PhyloTree};

/// Result of ancestral reconstruction for a single character (column).
#[derive(Debug, Clone)]
pub struct AncestralStates {
    /// Reconstructed state for each node (index = NodeId).
    /// Leaf nodes have the observed state; internal nodes have the inferred state.
    pub states: Vec<u8>,
    /// Number of changes (substitutions) on the tree.
    pub n_changes: usize,
}

/// Fitch parsimony reconstruction for a single character.
///
/// `leaf_states` maps leaf `NodeId` → observed character state.
///
/// # Errors
///
/// Returns an error if leaf states are incomplete or the tree is empty.
pub fn fitch(tree: &PhyloTree, leaf_states: &[(NodeId, u8)]) -> Result<AncestralStates> {
    if tree.node_count() == 0 {
        return Err(CyaneaError::InvalidInput("empty tree".into()));
    }

    let n = tree.node_count();
    let mut state_sets: Vec<HashSet<u8>> = vec![HashSet::new(); n];

    // Initialize leaf states
    for &(id, state) in leaf_states {
        if id >= n {
            return Err(CyaneaError::InvalidInput(format!(
                "node id {} out of range ({})",
                id, n
            )));
        }
        state_sets[id].insert(state);
    }

    // Verify all leaves have states
    for &leaf_id in &tree.leaves() {
        if state_sets[leaf_id].is_empty() {
            return Err(CyaneaError::InvalidInput(format!(
                "no state for leaf node {}",
                leaf_id
            )));
        }
    }

    // Bottom-up pass (post-order)
    let mut n_changes = 0;
    for id in tree.iter_postorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_leaf() {
            continue;
        }

        let child_sets: Vec<HashSet<u8>> = node
            .children
            .iter()
            .map(|&c| state_sets[c].clone())
            .collect();

        // Intersection of all children
        let mut intersection = child_sets[0].clone();
        for cs in &child_sets[1..] {
            intersection = intersection.intersection(cs).copied().collect();
        }

        if intersection.is_empty() {
            // Union if no common states
            let mut union = HashSet::new();
            for cs in &child_sets {
                union = union.union(cs).copied().collect();
            }
            state_sets[id] = union;
            n_changes += 1;
        } else {
            state_sets[id] = intersection;
        }
    }

    // Top-down pass (pre-order)
    let mut states = vec![0u8; n];

    for id in tree.iter_preorder() {
        let node = tree.get_node(id).unwrap();

        if node.is_root() || node.parent.is_none() {
            // Pick the smallest state from the set (deterministic choice)
            states[id] = *state_sets[id].iter().min().unwrap_or(&0);
        } else {
            let parent_state = states[node.parent.unwrap()];
            if state_sets[id].contains(&parent_state) {
                states[id] = parent_state;
            } else {
                states[id] = *state_sets[id].iter().min().unwrap_or(&0);
            }
        }
    }

    Ok(AncestralStates { states, n_changes })
}

/// Cost matrix for Sankoff parsimony.
///
/// A symmetric matrix giving the cost of transitioning between each pair
/// of character states. States are represented as indices 0..n_states.
#[derive(Debug, Clone)]
pub struct CostMatrix {
    /// Flattened n_states × n_states cost matrix.
    costs: Vec<f64>,
    n_states: usize,
}

impl CostMatrix {
    /// Create a uniform cost matrix where all transitions cost 1.0 and
    /// staying in the same state costs 0.0.
    pub fn uniform(n_states: usize) -> Self {
        let mut costs = vec![1.0; n_states * n_states];
        for i in 0..n_states {
            costs[i * n_states + i] = 0.0;
        }
        Self { costs, n_states }
    }

    /// Create from a flat cost vector.
    pub fn from_flat(costs: Vec<f64>, n_states: usize) -> Result<Self> {
        if costs.len() != n_states * n_states {
            return Err(CyaneaError::InvalidInput(format!(
                "cost matrix size {} does not match n_states²={}",
                costs.len(),
                n_states * n_states
            )));
        }
        Ok(Self { costs, n_states })
    }

    /// Cost of transitioning from state `from` to state `to`.
    pub fn cost(&self, from: usize, to: usize) -> f64 {
        self.costs[from * self.n_states + to]
    }

    /// Number of states.
    pub fn n_states(&self) -> usize {
        self.n_states
    }
}

/// Sankoff parsimony reconstruction for a single character.
///
/// `leaf_states` maps leaf `NodeId` → state index (0-based, must be < cost_matrix.n_states).
///
/// # Errors
///
/// Returns an error if the tree is empty or leaf states reference invalid state indices.
pub fn sankoff(
    tree: &PhyloTree,
    leaf_states: &[(NodeId, usize)],
    cost_matrix: &CostMatrix,
) -> Result<AncestralStates> {
    if tree.node_count() == 0 {
        return Err(CyaneaError::InvalidInput("empty tree".into()));
    }

    let n = tree.node_count();
    let k = cost_matrix.n_states();
    let inf = f64::MAX / 2.0;

    // Cost table: cost[node][state] = minimum cost of subtree if node has that state
    let mut cost_table = vec![vec![inf; k]; n];

    // Initialize leaves
    for &(id, state) in leaf_states {
        if id >= n {
            return Err(CyaneaError::InvalidInput(format!(
                "node id {} out of range ({})",
                id, n
            )));
        }
        if state >= k {
            return Err(CyaneaError::InvalidInput(format!(
                "state {} >= n_states {}",
                state, k
            )));
        }
        for s in 0..k {
            cost_table[id][s] = cost_matrix.cost(state, s);
        }
    }

    // Bottom-up (post-order)
    for id in tree.iter_postorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_leaf() {
            continue;
        }

        for s in 0..k {
            let mut total = 0.0;
            for &child in &node.children {
                let mut min_child_cost = inf;
                for t in 0..k {
                    let c = cost_matrix.cost(s, t) + cost_table[child][t];
                    if c < min_child_cost {
                        min_child_cost = c;
                    }
                }
                total += min_child_cost;
            }
            cost_table[id][s] = total;
        }
    }

    // Build a lookup of leaf observed states
    let mut leaf_observed: Vec<Option<usize>> = vec![None; n];
    for &(id, state) in leaf_states {
        leaf_observed[id] = Some(state);
    }

    // Top-down (pre-order) to select optimal states
    let mut states = vec![0usize; n];

    for id in tree.iter_preorder() {
        let node = tree.get_node(id).unwrap();

        // Leaves keep their observed state
        if let Some(observed) = leaf_observed[id] {
            states[id] = observed;
            continue;
        }

        if node.is_root() || node.parent.is_none() {
            // Pick state with minimum cost at root
            let mut min_cost = inf;
            let mut min_state = 0;
            for s in 0..k {
                if cost_table[id][s] < min_cost {
                    min_cost = cost_table[id][s];
                    min_state = s;
                }
            }
            states[id] = min_state;
        } else {
            let parent_state = states[node.parent.unwrap()];
            let mut min_cost = inf;
            let mut min_state = 0;
            for s in 0..k {
                let c = cost_matrix.cost(parent_state, s) + cost_table[id][s];
                if c < min_cost {
                    min_cost = c;
                    min_state = s;
                }
            }
            states[id] = min_state;
        }
    }

    // Count changes
    let mut n_changes = 0;
    for id in tree.iter_preorder() {
        let node = tree.get_node(id).unwrap();
        if let Some(parent) = node.parent {
            if states[id] != states[parent] {
                n_changes += 1;
            }
        }
    }

    // Convert to byte states
    let byte_states: Vec<u8> = states.iter().map(|&s| s as u8).collect();

    Ok(AncestralStates {
        states: byte_states,
        n_changes,
    })
}

/// Reconstruct ancestral sequences for multiple alignment columns.
///
/// `alignment` is a matrix of leaf sequences: each entry maps `NodeId` → sequence byte.
/// Returns a state vector for each column.
pub fn reconstruct_sequences(
    tree: &PhyloTree,
    alignment: &[(NodeId, &[u8])],
) -> Result<Vec<AncestralStates>> {
    if alignment.is_empty() {
        return Err(CyaneaError::InvalidInput("empty alignment".into()));
    }
    let seq_len = alignment[0].1.len();
    for &(id, seq) in alignment {
        if seq.len() != seq_len {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence at node {} has length {}, expected {}",
                id,
                seq.len(),
                seq_len
            )));
        }
    }

    let mut results = Vec::with_capacity(seq_len);
    for col in 0..seq_len {
        let leaf_states: Vec<(NodeId, u8)> = alignment
            .iter()
            .map(|&(id, seq)| (id, seq[col]))
            .collect();
        results.push(fitch(tree, &leaf_states)?);
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_tree() -> PhyloTree {
        // ((A,B),C)
        PhyloTree::from_newick("((A,B),C);").unwrap()
    }

    fn get_leaf_id(tree: &PhyloTree, name: &str) -> NodeId {
        for id in tree.leaves() {
            if tree.get_node(id).unwrap().name.as_deref() == Some(name) {
                return id;
            }
        }
        panic!("leaf {} not found", name);
    }

    #[test]
    fn fitch_all_same() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let b = get_leaf_id(&tree, "B");
        let c = get_leaf_id(&tree, "C");
        let result = fitch(&tree, &[(a, b'A'), (b, b'A'), (c, b'A')]).unwrap();
        assert_eq!(result.n_changes, 0);
        // All nodes should be 'A'
        for &s in &result.states {
            assert_eq!(s, b'A');
        }
    }

    #[test]
    fn fitch_one_change() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let b = get_leaf_id(&tree, "B");
        let c = get_leaf_id(&tree, "C");
        let result = fitch(&tree, &[(a, b'A'), (b, b'A'), (c, b'T')]).unwrap();
        // A and B agree, C differs → 1 change
        assert_eq!(result.n_changes, 1);
    }

    #[test]
    fn fitch_two_changes() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let b = get_leaf_id(&tree, "B");
        let c = get_leaf_id(&tree, "C");
        let result = fitch(&tree, &[(a, b'A'), (b, b'T'), (c, b'G')]).unwrap();
        // All different → at least 2 changes
        assert!(result.n_changes >= 2);
    }

    #[test]
    fn fitch_empty_tree_error() {
        let tree = PhyloTree::from_nodes(vec![], 0);
        assert!(tree.is_err());
    }

    #[test]
    fn fitch_missing_leaf_state() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        // Missing B and C
        let result = fitch(&tree, &[(a, b'A')]);
        assert!(result.is_err());
    }

    #[test]
    fn sankoff_uniform_all_same() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let b = get_leaf_id(&tree, "B");
        let c = get_leaf_id(&tree, "C");
        let cost = CostMatrix::uniform(4);
        let result = sankoff(&tree, &[(a, 0), (b, 0), (c, 0)], &cost).unwrap();
        assert_eq!(result.n_changes, 0);
    }

    #[test]
    fn sankoff_one_different() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let b = get_leaf_id(&tree, "B");
        let c = get_leaf_id(&tree, "C");
        let cost = CostMatrix::uniform(4);
        let result = sankoff(&tree, &[(a, 0), (b, 0), (c, 1)], &cost).unwrap();
        assert_eq!(result.n_changes, 1);
    }

    #[test]
    fn sankoff_invalid_state() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let cost = CostMatrix::uniform(4);
        assert!(sankoff(&tree, &[(a, 10)], &cost).is_err());
    }

    #[test]
    fn reconstruct_sequences_basic() {
        let tree = simple_tree();
        let a = get_leaf_id(&tree, "A");
        let b = get_leaf_id(&tree, "B");
        let c = get_leaf_id(&tree, "C");
        let alignment = vec![
            (a, b"ACG" as &[u8]),
            (b, b"ACG" as &[u8]),
            (c, b"ACG" as &[u8]),
        ];
        let results = reconstruct_sequences(&tree, &alignment).unwrap();
        assert_eq!(results.len(), 3);
        for r in &results {
            assert_eq!(r.n_changes, 0);
        }
    }

    #[test]
    fn cost_matrix_from_flat() {
        let costs = vec![0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0];
        let cm = CostMatrix::from_flat(costs, 3).unwrap();
        assert_eq!(cm.cost(0, 1), 1.0);
        assert_eq!(cm.cost(0, 2), 2.0);
    }
}
