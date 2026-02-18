//! Molecular clock dating for phylogenetic trees.
//!
//! Provides two methods for estimating divergence times:
//! - **Strict clock:** assumes a constant substitution rate across the tree,
//!   calibrated from a single node with a known age.
//! - **Root-to-tip regression:** linear regression of root-to-tip genetic
//!   distance against sampling dates, useful for serially-sampled data.

use crate::tree::{NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// A calibration point: a node with a known age.
#[derive(Debug, Clone)]
pub struct Calibration {
    /// The calibrated node.
    pub node_id: NodeId,
    /// The age of this node (time units, e.g. Mya).
    pub age: f64,
}

/// Result of molecular clock dating.
#[derive(Debug, Clone)]
pub struct DatingResult {
    /// Estimated age of each node (index = NodeId).
    pub node_ages: Vec<f64>,
    /// Estimated substitution rate.
    pub rate: f64,
}

/// Estimate divergence times under a strict molecular clock.
///
/// Assumes a constant substitution rate across the tree. The rate is
/// calibrated from a single node with a known age:
///
///   rate = root_to_calibration_distance / calibration_age
///
/// Node ages are then computed as:
///
///   age(v) = (max_root_to_tip - root_to_v_distance) / rate
pub fn strict_clock(tree: &PhyloTree, calibration: &Calibration) -> Result<DatingResult> {
    let n = tree.node_count();
    if calibration.node_id >= n {
        return Err(CyaneaError::InvalidInput("calibration node out of range".into()));
    }
    if calibration.age <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "calibration age must be positive".into(),
        ));
    }

    // Compute root-to-node distances for all nodes.
    let mut root_dist = vec![0.0f64; n];
    for id in tree.iter_preorder() {
        let node = tree.get_node(id).unwrap();
        if node.is_root() {
            root_dist[id] = 0.0;
        } else {
            root_dist[id] =
                root_dist[node.parent.unwrap()] + node.branch_length.unwrap_or(0.0);
        }
    }

    // Rate = distance_to_calibration / age.
    let cal_dist = root_dist[calibration.node_id];
    if cal_dist <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "calibration node has zero distance from root".into(),
        ));
    }
    let rate = cal_dist / calibration.age;

    // Maximum root-to-tip distance (needed to compute ages as time-from-present).
    let max_dist = tree
        .leaves()
        .iter()
        .map(|&id| root_dist[id])
        .fold(0.0f64, f64::max);

    // Node ages: age = (max_dist - root_dist) / rate.
    let node_ages: Vec<f64> = (0..n)
        .map(|id| {
            let age = (max_dist - root_dist[id]) / rate;
            if age < 0.0 { 0.0 } else { age }
        })
        .collect();

    Ok(DatingResult { node_ages, rate })
}

/// Estimate divergence times via root-to-tip regression.
///
/// Performs linear regression of root-to-tip genetic distance against
/// sampling dates for tip nodes. The substitution rate is estimated as
/// the slope, and the root age as -intercept/slope.
///
/// `tip_dates` maps leaf NodeIds to their sampling dates (in the same
/// time units as desired output ages).
pub fn root_to_tip_regression(
    tree: &PhyloTree,
    tip_dates: &[(NodeId, f64)],
) -> Result<DatingResult> {
    let n = tree.node_count();
    if tip_dates.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 tip dates for regression".into(),
        ));
    }
    for &(id, _) in tip_dates {
        if id >= n {
            return Err(CyaneaError::InvalidInput(format!(
                "tip node {} out of range",
                id
            )));
        }
    }

    // Compute root-to-node distances.
    let mut root_dist = vec![0.0f64; n];
    for id in tree.iter_preorder() {
        let node = tree.get_node(id).unwrap();
        if !node.is_root() {
            root_dist[id] =
                root_dist[node.parent.unwrap()] + node.branch_length.unwrap_or(0.0);
        }
    }

    // Linear regression: distance = rate * date + intercept.
    let m = tip_dates.len() as f64;
    let sum_x: f64 = tip_dates.iter().map(|&(_, d)| d).sum();
    let sum_y: f64 = tip_dates.iter().map(|&(id, _)| root_dist[id]).sum();
    let sum_xx: f64 = tip_dates.iter().map(|&(_, d)| d * d).sum();
    let sum_xy: f64 = tip_dates.iter().map(|&(id, d)| d * root_dist[id]).sum();

    let denom = m * sum_xx - sum_x * sum_x;
    if denom.abs() < 1e-30 {
        return Err(CyaneaError::InvalidInput(
            "degenerate regression: all dates are the same".into(),
        ));
    }

    let rate = (m * sum_xy - sum_x * sum_y) / denom;
    let intercept = (sum_y - rate * sum_x) / m;

    if rate <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "negative or zero substitution rate from regression".into(),
        ));
    }

    // Root age = -intercept / rate (time when distance = 0).
    let root_age = -intercept / rate;

    // Node ages: proportional to distance from root.
    // Max root-to-tip distance corresponds to most recent tip.
    let max_date: f64 = tip_dates
        .iter()
        .map(|&(_, d)| d)
        .fold(f64::MIN, f64::max);

    let node_ages: Vec<f64> = (0..n)
        .map(|id| {
            // date(node) = root_age + root_dist[node] / rate
            // age = max_date - date(node)
            let node_date = root_age + root_dist[id] / rate;
            let age = max_date - node_date;
            if age < 0.0 { 0.0 } else { age }
        })
        .collect();

    Ok(DatingResult { node_ages, rate })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strict_clock_ultrametric() {
        // Ultrametric tree: all root-to-tip distances equal.
        let tree =
            PhyloTree::from_newick("((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);").unwrap();
        let root_id = tree.root();
        // Calibrate root at age 1.0.
        // We need a non-root node for calibration. Use internal node.
        let internal = tree.get_node(root_id).unwrap().children[0];
        let cal = Calibration {
            node_id: internal,
            age: 0.5, // distance to internal = 0.5, so rate = 1.0
        };
        let result = strict_clock(&tree, &cal).unwrap();

        // All tips should have age ~0 (present).
        for &leaf_id in &tree.leaves() {
            assert!(
                result.node_ages[leaf_id].abs() < 1e-10,
                "leaf {} should have age 0, got {}",
                leaf_id,
                result.node_ages[leaf_id]
            );
        }
        // Root should have age ~1.0.
        assert!(
            (result.node_ages[root_id] - 1.0).abs() < 1e-10,
            "root age should be 1.0, got {}",
            result.node_ages[root_id]
        );
    }

    #[test]
    fn strict_clock_rate() {
        let tree = PhyloTree::from_newick("((A:0.2,B:0.2):0.3,C:0.5);").unwrap();
        let a_id = tree
            .leaves()
            .into_iter()
            .find(|&id| tree.get_node(id).unwrap().name.as_deref() == Some("A"))
            .unwrap();
        // Distance to A = 0.2 + 0.3 = 0.5, calibrate at age 10.0
        let cal = Calibration {
            node_id: a_id,
            age: 10.0,
        };
        let result = strict_clock(&tree, &cal).unwrap();
        assert!(
            (result.rate - 0.05).abs() < 1e-10,
            "rate should be 0.05, got {}",
            result.rate
        );
    }

    #[test]
    fn root_to_tip_linear() {
        // Tree with tips at different distances matching linear relationship.
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.2):0.1,C:0.3);").unwrap();
        let a_id = tree
            .leaves()
            .into_iter()
            .find(|&id| tree.get_node(id).unwrap().name.as_deref() == Some("A"))
            .unwrap();
        let b_id = tree
            .leaves()
            .into_iter()
            .find(|&id| tree.get_node(id).unwrap().name.as_deref() == Some("B"))
            .unwrap();
        let c_id = tree
            .leaves()
            .into_iter()
            .find(|&id| tree.get_node(id).unwrap().name.as_deref() == Some("C"))
            .unwrap();

        // Distances: A=0.2, B=0.3, C=0.3
        // Dates: choose so that distance ~ rate * date + intercept
        let tip_dates = vec![(a_id, 2.0), (b_id, 3.0), (c_id, 3.0)];
        let result = root_to_tip_regression(&tree, &tip_dates).unwrap();
        assert!(
            result.rate > 0.0,
            "rate should be positive, got {}",
            result.rate
        );
        assert!(result.rate.is_finite());
    }

    #[test]
    fn dating_validates_input() {
        let tree = PhyloTree::from_newick("(A:0.1,B:0.2);").unwrap();
        // Invalid calibration node.
        let cal = Calibration {
            node_id: 999,
            age: 1.0,
        };
        assert!(strict_clock(&tree, &cal).is_err());

        // Negative age.
        let cal = Calibration {
            node_id: 0,
            age: -1.0,
        };
        assert!(strict_clock(&tree, &cal).is_err());

        // Too few tip dates.
        let tip_dates = vec![(0, 1.0)];
        assert!(root_to_tip_regression(&tree, &tip_dates).is_err());
    }
}
