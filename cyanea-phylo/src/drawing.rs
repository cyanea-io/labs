//! Tree drawing coordinates for visualization.
//!
//! Computes 2D layout coordinates for phylogenetic trees in three styles:
//! rectangular (phylogram), cladogram, and radial.

use crate::tree::{NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Coordinates of a single node in the layout.
#[derive(Debug, Clone)]
pub struct NodeCoord {
    /// Node index in the tree.
    pub node_id: NodeId,
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}

/// An edge connecting two nodes.
#[derive(Debug, Clone)]
pub struct Edge {
    /// Source node.
    pub from: NodeId,
    /// Target node.
    pub to: NodeId,
}

/// Complete layout result.
#[derive(Debug, Clone)]
pub struct TreeLayout {
    /// Node coordinates.
    pub coords: Vec<NodeCoord>,
    /// Edge list.
    pub edges: Vec<Edge>,
}

/// Layout style for tree drawing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LayoutStyle {
    /// Phylogram: x = cumulative branch length from root.
    Rectangular,
    /// Cladogram: x = node depth (number of edges from root).
    Cladogram,
    /// Radial: polar transform of rectangular layout.
    Radial,
}

/// Compute a 2D layout for a phylogenetic tree.
///
/// # Arguments
///
/// * `tree` - The tree to lay out.
/// * `style` - Layout style (Rectangular, Cladogram, or Radial).
/// * `width` - Optional target width (default 1.0).
/// * `height` - Optional target height (default 1.0).
pub fn tree_layout(
    tree: &PhyloTree,
    style: LayoutStyle,
    width: Option<f64>,
    height: Option<f64>,
) -> Result<TreeLayout> {
    if tree.node_count() == 0 {
        return Err(CyaneaError::InvalidInput("empty tree".into()));
    }

    let w = width.unwrap_or(1.0);
    let h = height.unwrap_or(1.0);
    let n_nodes = tree.node_count();

    // Step 1: Compute raw x coordinates.
    let mut x_raw = vec![0.0f64; n_nodes];
    match style {
        LayoutStyle::Rectangular | LayoutStyle::Radial => {
            // x = cumulative branch length from root (preorder).
            for id in tree.iter_preorder() {
                let node = tree.get_node(id).unwrap();
                if node.is_root() {
                    x_raw[id] = 0.0;
                } else {
                    let parent_x = x_raw[node.parent.unwrap()];
                    x_raw[id] = parent_x + node.branch_length.unwrap_or(0.0);
                }
            }
        }
        LayoutStyle::Cladogram => {
            // First compute max depth.
            let max_depth = tree.depth() as f64;
            // x = depth, with leaves at max depth.
            let mut depth = vec![0usize; n_nodes];
            for id in tree.iter_preorder() {
                let node = tree.get_node(id).unwrap();
                if node.is_root() {
                    depth[id] = 0;
                } else {
                    depth[id] = depth[node.parent.unwrap()] + 1;
                }
            }
            // Leaves get max_depth, internals get their actual depth.
            for id in 0..n_nodes {
                let node = tree.get_node(id).unwrap();
                if node.is_leaf() {
                    x_raw[id] = max_depth;
                } else {
                    x_raw[id] = depth[id] as f64;
                }
            }
        }
    }

    // Step 2: Compute y coordinates.
    // Leaves: evenly spaced 0, 1, 2, ... in postorder.
    // Internals: midpoint of children's y range.
    let mut y_raw = vec![0.0f64; n_nodes];
    let mut leaf_counter = 0.0;
    for id in tree.iter_postorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_leaf() {
            y_raw[id] = leaf_counter;
            leaf_counter += 1.0;
        } else {
            let child_ys: Vec<f64> = node.children.iter().map(|&c| y_raw[c]).collect();
            let min_y = child_ys.iter().cloned().fold(f64::MAX, f64::min);
            let max_y = child_ys.iter().cloned().fold(f64::MIN, f64::max);
            y_raw[id] = (min_y + max_y) / 2.0;
        }
    }

    // Step 3: Scale to bounding box.
    let x_max = x_raw.iter().cloned().fold(0.0f64, f64::max);
    let y_max = y_raw.iter().cloned().fold(0.0f64, f64::max);

    let x_scale = if x_max > 0.0 { w / x_max } else { 1.0 };
    let y_scale = if y_max > 0.0 { h / y_max } else { 1.0 };

    // Step 4: Build coordinates.
    let coords = match style {
        LayoutStyle::Rectangular | LayoutStyle::Cladogram => {
            (0..n_nodes)
                .map(|id| NodeCoord {
                    node_id: id,
                    x: x_raw[id] * x_scale,
                    y: y_raw[id] * y_scale,
                })
                .collect()
        }
        LayoutStyle::Radial => {
            let n_leaves = tree.leaf_count() as f64;
            let angle_scale = if n_leaves > 1.0 {
                2.0 * std::f64::consts::PI / n_leaves
            } else {
                0.0
            };
            (0..n_nodes)
                .map(|id| {
                    let radius = x_raw[id] * x_scale;
                    let angle = y_raw[id] * angle_scale;
                    NodeCoord {
                        node_id: id,
                        x: radius * angle.cos(),
                        y: radius * angle.sin(),
                    }
                })
                .collect()
        }
    };

    // Step 5: Build edges.
    let mut edges = Vec::new();
    for id in 0..n_nodes {
        let node = tree.get_node(id).unwrap();
        if let Some(parent) = node.parent {
            edges.push(Edge {
                from: parent,
                to: id,
            });
        }
    }

    Ok(TreeLayout { coords, edges })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangular_leaf_x_equals_distance() {
        // Leaf x should be proportional to root-to-leaf distance.
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        let layout = tree_layout(&tree, LayoutStyle::Rectangular, None, None).unwrap();

        // Find leaf coordinates.
        let leaf_ids = tree.leaves();
        for &id in &leaf_ids {
            let coord = layout.coords.iter().find(|c| c.node_id == id).unwrap();
            assert!(coord.x > 0.0, "leaf {} should have x > 0", id);
        }
        // Verify relative ordering: leaf with longer path should have larger x.
        // A has path 0.1+0.3=0.4, B has 0.2+0.3=0.5, C has 0.4+0.6=1.0, D has 0.5+0.6=1.1
        let get_x = |name: &str| {
            let id = tree
                .leaves()
                .into_iter()
                .find(|&id| tree.get_node(id).unwrap().name.as_deref() == Some(name))
                .unwrap();
            layout.coords.iter().find(|c| c.node_id == id).unwrap().x
        };
        assert!(get_x("A") < get_x("D"));
        assert!(get_x("B") < get_x("C"));
    }

    #[test]
    fn rectangular_leaves_evenly_spaced() {
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let layout = tree_layout(&tree, LayoutStyle::Rectangular, None, None).unwrap();

        let mut leaf_ys: Vec<f64> = tree
            .leaves()
            .iter()
            .map(|&id| layout.coords.iter().find(|c| c.node_id == id).unwrap().y)
            .collect();
        leaf_ys.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Check that spacing is uniform.
        if leaf_ys.len() >= 2 {
            let spacing = leaf_ys[1] - leaf_ys[0];
            for i in 2..leaf_ys.len() {
                let s = leaf_ys[i] - leaf_ys[i - 1];
                assert!(
                    (s - spacing).abs() < 1e-10,
                    "non-uniform leaf spacing: {} vs {}",
                    s,
                    spacing
                );
            }
        }
    }

    #[test]
    fn cladogram_leaves_same_x() {
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        let layout = tree_layout(&tree, LayoutStyle::Cladogram, None, None).unwrap();

        let leaf_xs: Vec<f64> = tree
            .leaves()
            .iter()
            .map(|&id| layout.coords.iter().find(|c| c.node_id == id).unwrap().x)
            .collect();

        let first_x = leaf_xs[0];
        for &x in &leaf_xs[1..] {
            assert!(
                (x - first_x).abs() < 1e-10,
                "cladogram leaves should have same x: {} vs {}",
                x,
                first_x
            );
        }
    }

    #[test]
    fn radial_leaves_same_radius() {
        let tree =
            PhyloTree::from_newick("((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);").unwrap();
        let layout = tree_layout(&tree, LayoutStyle::Radial, None, None).unwrap();

        let leaf_radii: Vec<f64> = tree
            .leaves()
            .iter()
            .map(|&id| {
                let c = layout.coords.iter().find(|c| c.node_id == id).unwrap();
                (c.x * c.x + c.y * c.y).sqrt()
            })
            .collect();

        let first_r = leaf_radii[0];
        for &r in &leaf_radii[1..] {
            assert!(
                (r - first_r).abs() < 1e-10,
                "radial leaves should have same radius: {} vs {}",
                r,
                first_r
            );
        }
    }
}
