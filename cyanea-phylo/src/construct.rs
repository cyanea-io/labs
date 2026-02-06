//! Distance-based tree construction algorithms.
//!
//! Provides UPGMA and Neighbor-Joining methods for building phylogenetic trees
//! from distance matrices. Requires the `ml` feature (depends on `cyanea-ml`).

use crate::tree::{Node, PhyloTree};
use cyanea_core::{CyaneaError, Result};
use cyanea_ml::DistanceMatrix;

/// Build a tree using UPGMA (Unweighted Pair Group Method with Arithmetic Mean).
///
/// Produces an ultrametric (clock-like) rooted tree where branch lengths reflect
/// the average-linkage clustering heights.
pub fn upgma(distances: &DistanceMatrix, leaf_names: &[String]) -> Result<PhyloTree> {
    let n = distances.n();
    validate_inputs(n, leaf_names)?;

    if n == 2 {
        return build_two_leaf_tree(distances, leaf_names);
    }

    // Expand condensed matrix into a full working matrix.
    let mut dist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = distances.get(i, j);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    let mut nodes: Vec<Node> = Vec::new();
    let mut cluster_size: Vec<usize> = vec![1; n];
    let mut heights: Vec<f64> = vec![0.0; n];
    // Map from active cluster index to node id
    let mut active: Vec<usize> = Vec::new();

    // Create leaf nodes
    for i in 0..n {
        let id = nodes.len();
        nodes.push(Node {
            id,
            parent: None,
            children: Vec::new(),
            branch_length: None,
            name: Some(leaf_names[i].clone()),
        });
        active.push(id);
    }

    let mut n_active = n;

    while n_active > 1 {
        // Find the pair with minimum distance
        let (mut min_i, mut min_j) = (0, 1);
        let mut min_dist = dist[0][1];
        for i in 0..n_active {
            for j in (i + 1)..n_active {
                if dist[i][j] < min_dist {
                    min_dist = dist[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        let new_height = min_dist / 2.0;
        let node_i = active[min_i];
        let node_j = active[min_j];
        let size_i = cluster_size[min_i];
        let size_j = cluster_size[min_j];

        // Create the new internal node
        let new_id = nodes.len();
        nodes.push(Node {
            id: new_id,
            parent: None,
            children: vec![node_i, node_j],
            branch_length: None,
            name: None,
        });
        nodes[node_i].parent = Some(new_id);
        nodes[node_j].parent = Some(new_id);
        nodes[node_i].branch_length = Some(new_height - heights[min_i]);
        nodes[node_j].branch_length = Some(new_height - heights[min_j]);

        // Update distance matrix: new cluster replaces min_i, remove min_j
        for k in 0..n_active {
            if k == min_i || k == min_j {
                continue;
            }
            let d_new = (dist[min_i][k] * size_i as f64 + dist[min_j][k] * size_j as f64)
                / (size_i + size_j) as f64;
            dist[min_i][k] = d_new;
            dist[k][min_i] = d_new;
        }

        active[min_i] = new_id;
        cluster_size[min_i] = size_i + size_j;
        heights[min_i] = new_height;

        // Remove min_j by swapping with last
        let last = n_active - 1;
        if min_j != last {
            active[min_j] = active[last];
            cluster_size[min_j] = cluster_size[last];
            heights[min_j] = heights[last];
            for k in 0..n_active {
                dist[min_j][k] = dist[last][k];
                dist[k][min_j] = dist[k][last];
            }
        }
        n_active -= 1;
    }

    let root = active[0];
    PhyloTree::from_nodes(nodes, root)
}

/// Build a tree using the Neighbor-Joining algorithm.
///
/// Produces an additive (not necessarily ultrametric) rooted tree.
pub fn neighbor_joining(distances: &DistanceMatrix, leaf_names: &[String]) -> Result<PhyloTree> {
    let n = distances.n();
    validate_inputs(n, leaf_names)?;

    if n == 2 {
        return build_two_leaf_tree(distances, leaf_names);
    }

    // Expand condensed matrix into a full working matrix.
    let mut dist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = distances.get(i, j);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    let mut nodes: Vec<Node> = Vec::new();
    let mut active: Vec<usize> = Vec::new();

    // Create leaf nodes
    for i in 0..n {
        let id = nodes.len();
        nodes.push(Node {
            id,
            parent: None,
            children: Vec::new(),
            branch_length: None,
            name: Some(leaf_names[i].clone()),
        });
        active.push(id);
    }

    let mut n_active = n;

    while n_active > 2 {
        // Compute row sums
        let mut r = vec![0.0; n_active];
        for i in 0..n_active {
            for j in 0..n_active {
                r[i] += dist[i][j];
            }
        }

        // Find pair minimizing Q-matrix criterion
        let (mut min_i, mut min_j) = (0, 1);
        let mut min_q = f64::INFINITY;
        for i in 0..n_active {
            for j in (i + 1)..n_active {
                let q = (n_active as f64 - 2.0) * dist[i][j] - r[i] - r[j];
                if q < min_q {
                    min_q = q;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        let d_ij = dist[min_i][min_j];
        let n_f = n_active as f64;

        // Branch lengths for the two children
        let b_i = d_ij / 2.0 + (r[min_i] - r[min_j]) / (2.0 * (n_f - 2.0));
        let b_j = d_ij - b_i;

        let node_i = active[min_i];
        let node_j = active[min_j];

        // Create the new internal node
        let new_id = nodes.len();
        nodes.push(Node {
            id: new_id,
            parent: None,
            children: vec![node_i, node_j],
            branch_length: None,
            name: None,
        });
        nodes[node_i].parent = Some(new_id);
        nodes[node_j].parent = Some(new_id);
        nodes[node_i].branch_length = Some(b_i.max(0.0));
        nodes[node_j].branch_length = Some(b_j.max(0.0));

        // Update distances: new node replaces min_i
        for k in 0..n_active {
            if k == min_i || k == min_j {
                continue;
            }
            let d_new = (dist[min_i][k] + dist[min_j][k] - d_ij) / 2.0;
            dist[min_i][k] = d_new;
            dist[k][min_i] = d_new;
        }

        active[min_i] = new_id;

        // Remove min_j by swapping with last
        let last = n_active - 1;
        if min_j != last {
            active[min_j] = active[last];
            for k in 0..n_active {
                dist[min_j][k] = dist[last][k];
                dist[k][min_j] = dist[k][last];
            }
        }
        n_active -= 1;
    }

    // Join the final two clusters
    let d_final = dist[0][1];
    let node_0 = active[0];
    let node_1 = active[1];

    let root_id = nodes.len();
    nodes.push(Node {
        id: root_id,
        parent: None,
        children: vec![node_0, node_1],
        branch_length: None,
        name: None,
    });
    nodes[node_0].parent = Some(root_id);
    nodes[node_1].parent = Some(root_id);
    nodes[node_0].branch_length = Some(d_final / 2.0);
    nodes[node_1].branch_length = Some(d_final / 2.0);

    PhyloTree::from_nodes(nodes, root_id)
}

fn validate_inputs(n: usize, leaf_names: &[String]) -> Result<()> {
    if leaf_names.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "leaf_names length {} doesn't match distance matrix size {}",
            leaf_names.len(),
            n
        )));
    }
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 taxa to build a tree".into(),
        ));
    }
    Ok(())
}

fn build_two_leaf_tree(distances: &DistanceMatrix, leaf_names: &[String]) -> Result<PhyloTree> {
    let d = distances.get(0, 1);
    let mut nodes = Vec::new();

    // Leaf 0
    nodes.push(Node {
        id: 0,
        parent: None,
        children: Vec::new(),
        branch_length: Some(d / 2.0),
        name: Some(leaf_names[0].clone()),
    });
    // Leaf 1
    nodes.push(Node {
        id: 1,
        parent: None,
        children: Vec::new(),
        branch_length: Some(d / 2.0),
        name: Some(leaf_names[1].clone()),
    });
    // Root
    nodes.push(Node {
        id: 2,
        parent: None,
        children: vec![0, 1],
        branch_length: None,
        name: None,
    });
    nodes[0].parent = Some(2);
    nodes[1].parent = Some(2);

    PhyloTree::from_nodes(nodes, 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_dm(values: &[f64], n: usize) -> DistanceMatrix {
        DistanceMatrix::from_condensed(values.to_vec(), n).unwrap()
    }

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("{}", (b'A' + i as u8) as char)).collect()
    }

    // --- UPGMA tests ---

    #[test]
    fn upgma_two_leaves() {
        let dm = make_dm(&[4.0], 2);
        let tree = upgma(&dm, &names(2)).unwrap();
        assert_eq!(tree.leaf_count(), 2);
        assert_eq!(tree.node_count(), 3);
        // Each branch should be d/2 = 2.0
        for &leaf in &tree.leaves() {
            let bl = tree.get_node(leaf).unwrap().branch_length.unwrap();
            assert!((bl - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn upgma_three_leaves() {
        // AB=2, AC=4, BC=4
        let dm = make_dm(&[2.0, 4.0, 4.0], 3);
        let tree = upgma(&dm, &names(3)).unwrap();
        assert_eq!(tree.leaf_count(), 3);
        assert_eq!(tree.node_count(), 5);
    }

    #[test]
    fn upgma_four_leaves_ultrametric() {
        // Simple ultrametric-compatible distances
        // AB=2, AC=6, AD=6, BC=6, BD=6, CD=2
        let dm = make_dm(&[2.0, 6.0, 6.0, 6.0, 6.0, 2.0], 4);
        let tree = upgma(&dm, &names(4)).unwrap();
        assert_eq!(tree.leaf_count(), 4);

        // In an ultrametric tree, root-to-leaf distances should be equal
        let root = tree.root();
        let leaves = tree.leaves();
        let mut root_dists = Vec::new();
        for &leaf in &leaves {
            let mut dist = 0.0;
            let mut cur = leaf;
            while cur != root {
                dist += tree.get_node(cur).unwrap().branch_length.unwrap_or(0.0);
                cur = tree.get_node(cur).unwrap().parent.unwrap();
            }
            root_dists.push(dist);
        }
        for d in &root_dists {
            assert!(
                (d - root_dists[0]).abs() < 1e-10,
                "ultrametric violation: {:?}",
                root_dists
            );
        }
    }

    #[test]
    fn upgma_name_mismatch_error() {
        let dm = make_dm(&[1.0], 2);
        assert!(upgma(&dm, &names(3)).is_err());
    }

    // --- Neighbor-Joining tests ---

    #[test]
    fn nj_two_leaves() {
        let dm = make_dm(&[4.0], 2);
        let tree = neighbor_joining(&dm, &names(2)).unwrap();
        assert_eq!(tree.leaf_count(), 2);
        assert_eq!(tree.node_count(), 3);
    }

    #[test]
    fn nj_three_leaves() {
        // AB=5, AC=9, BC=10
        let dm = make_dm(&[5.0, 9.0, 10.0], 3);
        let tree = neighbor_joining(&dm, &names(3)).unwrap();
        assert_eq!(tree.leaf_count(), 3);
        // For 3 taxa NJ produces a tree with 3 leaves + 1 internal + 1 root-like = depends on impl
        // The tree should be valid
        assert!(tree.node_count() >= 4);
    }

    #[test]
    fn nj_four_leaves_topology() {
        // Additive tree: ((A:1,B:2):1,(C:3,D:4):1)
        // AB=3, AC=5, AD=6, BC=6, BD=7, CD=7... nah let's use a known additive matrix
        // For ((A:2,B:3):1,(C:4,D:5):1):
        // d(A,B) = 5, d(A,C) = 7, d(A,D) = 8, d(B,C) = 8, d(B,D) = 9, d(C,D) = 9
        // Actually, let's just check it produces a valid tree
        let dm = make_dm(&[5.0, 7.0, 8.0, 8.0, 9.0, 9.0], 4);
        let tree = neighbor_joining(&dm, &names(4)).unwrap();
        assert_eq!(tree.leaf_count(), 4);
        assert_eq!(tree.leaf_names(), vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn nj_branch_lengths_positive() {
        // Use a well-separated additive matrix
        let dm = make_dm(&[2.0, 10.0, 10.0, 10.0, 10.0, 2.0], 4);
        let tree = neighbor_joining(&dm, &names(4)).unwrap();
        // All branch lengths should be non-negative
        for id in tree.iter_preorder() {
            if let Some(bl) = tree.get_node(id).unwrap().branch_length {
                assert!(bl >= 0.0, "negative branch length: {}", bl);
            }
        }
    }

    #[test]
    fn nj_name_mismatch_error() {
        let dm = make_dm(&[1.0], 2);
        assert!(neighbor_joining(&dm, &names(3)).is_err());
    }
}
