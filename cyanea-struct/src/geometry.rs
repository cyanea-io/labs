//! Coordinate geometry: distances, angles, dihedrals, RMSD.

use cyanea_core::{CyaneaError, Result};

use crate::types::{Atom, Point3D};

use alloc::format;

/// Euclidean distance between two atoms.
pub fn distance(a1: &Atom, a2: &Atom) -> f64 {
    distance_points(&a1.coords, &a2.coords)
}

/// Euclidean distance between two points.
pub fn distance_points(p1: &Point3D, p2: &Point3D) -> f64 {
    p1.distance_to(p2)
}

/// Bond angle in degrees at the central atom `a2`.
pub fn angle(a1: &Atom, a2: &Atom, a3: &Atom) -> f64 {
    angle_points(&a1.coords, &a2.coords, &a3.coords)
}

/// Bond angle in degrees at the central point `p2`.
pub fn angle_points(p1: &Point3D, p2: &Point3D, p3: &Point3D) -> f64 {
    let v1 = p1.sub(p2);
    let v2 = p3.sub(p2);
    let cos_angle = v1.dot(&v2) / (v1.norm() * v2.norm());
    // Clamp for numerical safety
    cos_angle.clamp(-1.0, 1.0).acos().to_degrees()
}

/// Dihedral (torsion) angle in degrees defined by four atoms.
pub fn dihedral(a1: &Atom, a2: &Atom, a3: &Atom, a4: &Atom) -> f64 {
    dihedral_points(&a1.coords, &a2.coords, &a3.coords, &a4.coords)
}

/// Dihedral (torsion) angle in degrees defined by four points.
pub fn dihedral_points(p1: &Point3D, p2: &Point3D, p3: &Point3D, p4: &Point3D) -> f64 {
    let b1 = p2.sub(p1);
    let b2 = p3.sub(p2);
    let b3 = p4.sub(p3);

    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);

    let m1 = n1.cross(&b2.normalize());

    let x = n1.dot(&n2);
    let y = m1.dot(&n2);

    (-y).atan2(x).to_degrees()
}

/// Geometric center of mass (unweighted) of a slice of atoms.
pub fn center_of_mass(atoms: &[&Atom]) -> Point3D {
    let points: alloc::vec::Vec<Point3D> = atoms.iter().map(|a| a.coords).collect();
    center_of_mass_points(&points)
}

/// Geometric center of mass (unweighted) of a slice of points.
pub fn center_of_mass_points(points: &[Point3D]) -> Point3D {
    if points.is_empty() {
        return Point3D::zero();
    }
    let mut sum = Point3D::zero();
    for p in points {
        sum = sum.add(p);
    }
    sum.scale(1.0 / points.len() as f64)
}

/// RMSD between two equal-length sets of atoms (no alignment, direct comparison).
pub fn rmsd(atoms1: &[&Atom], atoms2: &[&Atom]) -> Result<f64> {
    if atoms1.len() != atoms2.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "atom set sizes differ: {} vs {}",
            atoms1.len(),
            atoms2.len()
        )));
    }
    if atoms1.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "cannot compute RMSD of empty atom sets".into(),
        ));
    }
    let p1: alloc::vec::Vec<Point3D> = atoms1.iter().map(|a| a.coords).collect();
    let p2: alloc::vec::Vec<Point3D> = atoms2.iter().map(|a| a.coords).collect();
    rmsd_points(&p1, &p2)
}

/// RMSD between two equal-length slices of points.
pub fn rmsd_points(points1: &[Point3D], points2: &[Point3D]) -> Result<f64> {
    if points1.len() != points2.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "point set sizes differ: {} vs {}",
            points1.len(),
            points2.len()
        )));
    }
    if points1.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "cannot compute RMSD of empty point sets".into(),
        ));
    }
    let sum: f64 = points1
        .iter()
        .zip(points2.iter())
        .map(|(a, b)| {
            let d = a.sub(b);
            d.dot(&d)
        })
        .sum();
    Ok((sum / points1.len() as f64).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use crate::types::Atom;

    fn make_atom(x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial: 1,
            name: "CA".into(),
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
    fn test_distance() {
        let a = make_atom(0.0, 0.0, 0.0);
        let b = make_atom(3.0, 4.0, 0.0);
        assert!((distance(&a, &b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_90() {
        let p1 = Point3D::new(1.0, 0.0, 0.0);
        let p2 = Point3D::new(0.0, 0.0, 0.0);
        let p3 = Point3D::new(0.0, 1.0, 0.0);
        assert!((angle_points(&p1, &p2, &p3) - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_180() {
        let p1 = Point3D::new(-1.0, 0.0, 0.0);
        let p2 = Point3D::new(0.0, 0.0, 0.0);
        let p3 = Point3D::new(1.0, 0.0, 0.0);
        assert!((angle_points(&p1, &p2, &p3) - 180.0).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral() {
        // Trans conformation: ~180 degrees
        let p1 = Point3D::new(1.0, 0.0, 0.0);
        let p2 = Point3D::new(0.0, 0.0, 0.0);
        let p3 = Point3D::new(0.0, 1.0, 0.0);
        let p4 = Point3D::new(-1.0, 1.0, 0.0);
        let d = dihedral_points(&p1, &p2, &p3, &p4);
        assert!((d.abs() - 180.0).abs() < 1e-10);
    }

    #[test]
    fn test_center_of_mass() {
        let points = vec![
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(2.0, 0.0, 0.0),
            Point3D::new(0.0, 2.0, 0.0),
        ];
        let com = center_of_mass_points(&points);
        assert!((com.x - 2.0 / 3.0).abs() < 1e-10);
        assert!((com.y - 2.0 / 3.0).abs() < 1e-10);
        assert!((com.z).abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_identical() {
        let points = vec![
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(0.0, 0.0, 1.0),
        ];
        assert!((rmsd_points(&points, &points).unwrap()).abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_mismatch_error() {
        let p1 = vec![Point3D::new(0.0, 0.0, 0.0)];
        let p2 = vec![
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
        ];
        assert!(rmsd_points(&p1, &p2).is_err());
    }
}
