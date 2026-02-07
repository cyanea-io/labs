//! Structural superposition via the Kabsch algorithm.
//!
//! Finds the optimal rigid-body rotation (and translation) that minimizes RMSD
//! between two sets of corresponding points.

use cyanea_core::{CyaneaError, Result, Scored};

use crate::geometry::center_of_mass_points;
use crate::linalg::{svd_3x3, Matrix3x3};
use crate::types::{Atom, Point3D};

use alloc::format;
use alloc::vec::Vec;

/// Result of a Kabsch superposition.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SuperpositionResult {
    /// RMSD after optimal superposition.
    pub rmsd: f64,
    /// 3x3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Translation vector applied after rotation.
    pub translation: Point3D,
    /// Transformed coordinates of the mobile set after superposition.
    pub transformed_coords: Vec<Point3D>,
}

impl Scored for SuperpositionResult {
    fn score(&self) -> f64 {
        -self.rmsd
    }
}

/// Kabsch superposition on atom coordinates.
///
/// `atoms1` is the reference (fixed) set, `atoms2` is the mobile set that gets
/// rotated and translated to minimize RMSD.
pub fn kabsch(atoms1: &[&Atom], atoms2: &[&Atom]) -> Result<SuperpositionResult> {
    let p1: Vec<Point3D> = atoms1.iter().map(|a| a.coords).collect();
    let p2: Vec<Point3D> = atoms2.iter().map(|a| a.coords).collect();
    kabsch_points(&p1, &p2)
}

/// Kabsch superposition on point coordinates.
pub fn kabsch_points(
    points1: &[Point3D],
    points2: &[Point3D],
) -> Result<SuperpositionResult> {
    if points1.len() != points2.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "point set sizes differ: {} vs {}",
            points1.len(),
            points2.len()
        )));
    }
    if points1.len() < 3 {
        return Err(CyaneaError::InvalidInput(
            "need at least 3 points for Kabsch superposition".into(),
        ));
    }

    let n = points1.len();

    // Step 1: center both sets
    let com1 = center_of_mass_points(points1);
    let com2 = center_of_mass_points(points2);

    let centered1: Vec<Point3D> = points1.iter().map(|p| p.sub(&com1)).collect();
    let centered2: Vec<Point3D> = points2.iter().map(|p| p.sub(&com2)).collect();

    // Step 2: compute cross-covariance matrix H = P2^T * P1
    let mut h = Matrix3x3::zeros();
    for i in 0..n {
        let p = &centered2[i];
        let q = &centered1[i];
        h.data[0][0] += p.x * q.x;
        h.data[0][1] += p.x * q.y;
        h.data[0][2] += p.x * q.z;
        h.data[1][0] += p.y * q.x;
        h.data[1][1] += p.y * q.y;
        h.data[1][2] += p.y * q.z;
        h.data[2][0] += p.z * q.x;
        h.data[2][1] += p.z * q.y;
        h.data[2][2] += p.z * q.z;
    }

    // Step 3: SVD of H
    let svd = svd_3x3(&h);

    // Step 4: R = V * U^T, with reflection correction
    let v = svd.vt.transpose();
    let ut = svd.u.transpose();
    let mut r = v.multiply(&ut);

    // Fix reflection: if det(R) < 0, negate the column of V corresponding
    // to the smallest singular value (always column 2 after sorting)
    if r.determinant() < 0.0 {
        let mut v_fixed = v;
        for row in 0..3 {
            v_fixed.data[row][2] = -v_fixed.data[row][2];
        }
        r = v_fixed.multiply(&ut);
    }

    // Step 5: apply rotation and compute translation + RMSD
    let mut transformed = Vec::with_capacity(n);
    let mut sum_sq = 0.0;
    for i in 0..n {
        let rotated = r.apply(&centered2[i]);
        let final_point = rotated.add(&com1);
        let diff = final_point.sub(&points1[i]);
        sum_sq += diff.dot(&diff);
        transformed.push(final_point);
    }

    let rmsd = (sum_sq / n as f64).sqrt();

    Ok(SuperpositionResult {
        rmsd,
        rotation: r.data,
        translation: com1.sub(&r.apply(&com2)),
        transformed_coords: transformed,
    })
}

/// Align two structures using only alpha-carbon atoms.
pub fn align_structures_by_ca(
    atoms1: &[&Atom],
    atoms2: &[&Atom],
) -> Result<SuperpositionResult> {
    let ca1: Vec<&Atom> = atoms1.iter().copied().filter(|a| a.is_alpha_carbon()).collect();
    let ca2: Vec<&Atom> = atoms2.iter().copied().filter(|a| a.is_alpha_carbon()).collect();

    if ca1.len() != ca2.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "different number of CA atoms: {} vs {}",
            ca1.len(),
            ca2.len()
        )));
    }

    kabsch(&ca1, &ca2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use crate::types::Atom;

    fn make_atom(name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial: 1,
            name: name.into(),
            alt_loc: None,
            coords: Point3D::new(x, y, z),
            occupancy: 1.0,
            temp_factor: 0.0,
            element: None,
            charge: None,
            is_hetatm: false,
        }
    }

    #[test]
    fn identical_points_rmsd_zero() {
        let points = vec![
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(0.0, 0.0, 1.0),
        ];
        let result = kabsch_points(&points, &points).unwrap();
        assert!(result.rmsd < 1e-6, "RMSD should be ~0, got {}", result.rmsd);
    }

    #[test]
    fn translated_points() {
        let p1 = vec![
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(0.0, 0.0, 1.0),
        ];
        let p2: Vec<Point3D> = p1.iter().map(|p: &Point3D| p.add(&Point3D::new(10.0, 20.0, 30.0))).collect();
        let result = kabsch_points(&p1, &p2).unwrap();
        assert!(result.rmsd < 1e-6, "RMSD should be ~0 for translated set, got {}", result.rmsd);
    }

    #[test]
    fn rotated_points() {
        // 90-degree rotation around Z axis
        let p1 = vec![
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(-1.0, 0.0, 0.0),
            Point3D::new(0.0, -1.0, 0.0),
        ];
        let p2 = vec![
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(-1.0, 0.0, 0.0),
            Point3D::new(0.0, -1.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
        ];
        let result = kabsch_points(&p1, &p2).unwrap();
        assert!(result.rmsd < 1e-6, "RMSD should be ~0 for rotated set, got {}", result.rmsd);
    }

    #[test]
    fn mismatched_lengths_error() {
        let p1 = vec![Point3D::new(0.0, 0.0, 0.0); 3];
        let p2 = vec![Point3D::new(0.0, 0.0, 0.0); 4];
        assert!(kabsch_points(&p1, &p2).is_err());
    }

    #[test]
    fn align_by_ca() {
        let atoms1 = vec![
            make_atom("N", 0.0, 0.0, 0.0),
            make_atom("CA", 1.0, 0.0, 0.0),
            make_atom("C", 2.0, 0.0, 0.0),
            make_atom("N", 3.0, 0.0, 0.0),
            make_atom("CA", 4.0, 0.0, 0.0),
            make_atom("C", 5.0, 0.0, 0.0),
            make_atom("N", 6.0, 0.0, 0.0),
            make_atom("CA", 7.0, 0.0, 0.0),
            make_atom("C", 8.0, 0.0, 0.0),
        ];
        let atoms2 = vec![
            make_atom("N", 0.0, 0.0, 5.0),
            make_atom("CA", 1.0, 0.0, 5.0),
            make_atom("C", 2.0, 0.0, 5.0),
            make_atom("N", 3.0, 0.0, 5.0),
            make_atom("CA", 4.0, 0.0, 5.0),
            make_atom("C", 5.0, 0.0, 5.0),
            make_atom("N", 6.0, 0.0, 5.0),
            make_atom("CA", 7.0, 0.0, 5.0),
            make_atom("C", 8.0, 0.0, 5.0),
        ];
        let refs1: Vec<&Atom> = atoms1.iter().collect();
        let refs2: Vec<&Atom> = atoms2.iter().collect();
        let result = align_structures_by_ca(&refs1, &refs2).unwrap();
        assert!(result.rmsd < 1e-6);
    }
}
