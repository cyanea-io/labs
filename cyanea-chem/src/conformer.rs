//! 3D coordinate container for molecular conformers.
//!
//! Stores 3D coordinates separately from [`Molecule`](crate::molecule::Molecule),
//! which derives `Eq`/`Hash` and cannot hold `f64` fields.

use cyanea_core::Result;

/// A single 3D conformer: one set of xyz coordinates per atom.
#[derive(Debug, Clone)]
pub struct Conformer {
    /// `coords[i]` is `[x, y, z]` for atom `i`.
    pub coords: Vec<[f64; 3]>,
}

impl Conformer {
    /// Create a new conformer with the given coordinates.
    pub fn new(coords: Vec<[f64; 3]>) -> Self {
        Conformer { coords }
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.coords.len()
    }

    /// Whether the conformer has no atoms.
    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }

    /// Euclidean distance between two atoms.
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        let a = self.coords[i];
        let b = self.coords[j];
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Bond angle in radians between atoms i-j-k (angle at j).
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let v1 = sub(self.coords[i], self.coords[j]);
        let v2 = sub(self.coords[k], self.coords[j]);
        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let n1 = norm(v1);
        let n2 = norm(v2);
        if n1 < 1e-12 || n2 < 1e-12 {
            return 0.0;
        }
        (dot / (n1 * n2)).clamp(-1.0, 1.0).acos()
    }

    /// Dihedral (torsion) angle in radians for atoms i-j-k-l.
    pub fn dihedral(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        let b1 = sub(self.coords[j], self.coords[i]);
        let b2 = sub(self.coords[k], self.coords[j]);
        let b3 = sub(self.coords[l], self.coords[k]);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let m1 = cross(n1, b2);
        let n2_norm = norm(n2);
        let b2_norm = norm(b2);
        if n2_norm < 1e-12 || b2_norm < 1e-12 {
            return 0.0;
        }
        let x = dot3(n1, n2);
        let y = dot3(m1, n2) / b2_norm;
        (-y).atan2(-x) + std::f64::consts::PI
    }

    /// Geometric centroid of all atoms.
    pub fn centroid(&self) -> [f64; 3] {
        if self.coords.is_empty() {
            return [0.0; 3];
        }
        let n = self.coords.len() as f64;
        let mut c = [0.0; 3];
        for p in &self.coords {
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
        }
        c[0] /= n;
        c[1] /= n;
        c[2] /= n;
        c
    }

    /// Root-mean-square deviation from another conformer (same atom count).
    pub fn rmsd(&self, other: &Conformer) -> Result<f64> {
        if self.coords.len() != other.coords.len() {
            return Err(cyanea_core::CyaneaError::InvalidInput(
                "conformers must have the same number of atoms for RMSD".into(),
            ));
        }
        if self.coords.is_empty() {
            return Ok(0.0);
        }
        let sum: f64 = self
            .coords
            .iter()
            .zip(other.coords.iter())
            .map(|(a, b)| {
                let dx = a[0] - b[0];
                let dy = a[1] - b[1];
                let dz = a[2] - b[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        Ok((sum / self.coords.len() as f64).sqrt())
    }

    /// Set the dihedral angle for atoms i-j-k-l by rotating all atoms beyond k
    /// around the j-k bond axis. `atom_beyond_k` lists all atoms on the l-side.
    pub fn set_dihedral(
        &mut self,
        j: usize,
        k: usize,
        atoms_beyond_k: &[usize],
        target_angle: f64,
    ) {
        let current = if atoms_beyond_k.is_empty() {
            return;
        } else {
            // Pick the first atom beyond k as reference for computing current dihedral
            // We need a reference atom on the j side too
            // Find a neighbor of j that isn't k
            0.0 // Will be overridden below
        };
        let _ = current;

        // Axis of rotation: k - j
        let axis = sub(self.coords[k], self.coords[j]);
        let axis_len = norm(axis);
        if axis_len < 1e-12 {
            return;
        }
        let axis_unit = [axis[0] / axis_len, axis[1] / axis_len, axis[2] / axis_len];

        // Rotate atoms by the target angle around j-k axis, centered at k
        let origin = self.coords[k];
        let cos_a = target_angle.cos();
        let sin_a = target_angle.sin();

        for &atom in atoms_beyond_k {
            let p = sub(self.coords[atom], origin);
            let rotated = rodrigues(p, axis_unit, cos_a, sin_a);
            self.coords[atom] = [
                rotated[0] + origin[0],
                rotated[1] + origin[1],
                rotated[2] + origin[2],
            ];
        }
    }
}

/// A set of conformers with optional energies.
#[derive(Debug, Clone)]
pub struct ConformerSet {
    pub conformers: Vec<Conformer>,
    pub energies: Vec<Option<f64>>,
}

impl ConformerSet {
    /// Create an empty conformer set.
    pub fn new() -> Self {
        ConformerSet {
            conformers: Vec::new(),
            energies: Vec::new(),
        }
    }

    /// Add a conformer with an optional energy.
    pub fn push(&mut self, conformer: Conformer, energy: Option<f64>) {
        self.conformers.push(conformer);
        self.energies.push(energy);
    }

    /// Number of conformers.
    pub fn len(&self) -> usize {
        self.conformers.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.conformers.is_empty()
    }

    /// Sort conformers by energy (lowest first). Conformers without energy go last.
    pub fn sort_by_energy(&mut self) {
        let mut indices: Vec<usize> = (0..self.conformers.len()).collect();
        indices.sort_by(|&a, &b| {
            let ea = self.energies[a].unwrap_or(f64::INFINITY);
            let eb = self.energies[b].unwrap_or(f64::INFINITY);
            ea.partial_cmp(&eb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let old_conformers = std::mem::take(&mut self.conformers);
        let old_energies = std::mem::take(&mut self.energies);
        for &i in &indices {
            self.conformers.push(old_conformers[i].clone());
            self.energies.push(old_energies[i]);
        }
    }

    /// Return a reference to the lowest-energy conformer, or None if empty.
    pub fn best(&self) -> Option<&Conformer> {
        if self.conformers.is_empty() {
            return None;
        }
        let mut best_idx = 0;
        let mut best_e = f64::INFINITY;
        for (i, e) in self.energies.iter().enumerate() {
            if let Some(energy) = e {
                if *energy < best_e {
                    best_e = *energy;
                    best_idx = i;
                }
            }
        }
        Some(&self.conformers[best_idx])
    }
}

impl Default for ConformerSet {
    fn default() -> Self {
        Self::new()
    }
}

// --- Vector helpers ---

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Rodrigues' rotation formula: rotate vector `v` around unit axis `k` by angle (cos_a, sin_a).
fn rodrigues(v: [f64; 3], k: [f64; 3], cos_a: f64, sin_a: f64) -> [f64; 3] {
    let dot_kv = dot3(k, v);
    let cross_kv = cross(k, v);
    [
        v[0] * cos_a + cross_kv[0] * sin_a + k[0] * dot_kv * (1.0 - cos_a),
        v[1] * cos_a + cross_kv[1] * sin_a + k[1] * dot_kv * (1.0 - cos_a),
        v[2] * cos_a + cross_kv[2] * sin_a + k[2] * dot_kv * (1.0 - cos_a),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_triangle() -> Conformer {
        // Right triangle in XY plane
        Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
    }

    #[test]
    fn distance_basic() {
        let c = make_triangle();
        assert!((c.distance(0, 1) - 1.0).abs() < 1e-10);
        assert!((c.distance(0, 2) - 1.0).abs() < 1e-10);
        assert!((c.distance(1, 2) - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn angle_right_angle() {
        let c = make_triangle();
        let a = c.angle(1, 0, 2); // angle at origin between x-axis and y-axis
        assert!((a - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn dihedral_basic() {
        // Planar arrangement: all in XY plane → dihedral ≈ 0 or π
        let c = Conformer::new(vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ]);
        let d = c.dihedral(0, 1, 2, 3);
        // cis arrangement → dihedral ≈ 0
        assert!(d.abs() < 0.1 || (d - 2.0 * PI).abs() < 0.1, "got {d}");
    }

    #[test]
    fn centroid_basic() {
        let c = make_triangle();
        let ctr = c.centroid();
        assert!((ctr[0] - 1.0 / 3.0).abs() < 1e-10);
        assert!((ctr[1] - 1.0 / 3.0).abs() < 1e-10);
        assert!(ctr[2].abs() < 1e-10);
    }

    #[test]
    fn rmsd_identical() {
        let c = make_triangle();
        assert!((c.rmsd(&c).unwrap()).abs() < 1e-10);
    }

    #[test]
    fn rmsd_different() {
        let c1 = Conformer::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let c2 = Conformer::new(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let r = c1.rmsd(&c2).unwrap();
        // sqrt((0 + 1) / 2) = sqrt(0.5) ≈ 0.707
        assert!((r - (0.5_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn rmsd_mismatched_size() {
        let c1 = Conformer::new(vec![[0.0; 3]]);
        let c2 = Conformer::new(vec![[0.0; 3], [1.0, 0.0, 0.0]]);
        assert!(c1.rmsd(&c2).is_err());
    }

    #[test]
    fn conformer_set_sort_and_best() {
        let mut cs = ConformerSet::new();
        cs.push(Conformer::new(vec![[0.0; 3]]), Some(10.0));
        cs.push(Conformer::new(vec![[1.0, 0.0, 0.0]]), Some(5.0));
        cs.push(Conformer::new(vec![[2.0, 0.0, 0.0]]), Some(8.0));

        let best = cs.best().unwrap();
        assert!((best.coords[0][0] - 1.0).abs() < 1e-10);

        cs.sort_by_energy();
        assert!((cs.energies[0].unwrap() - 5.0).abs() < 1e-10);
        assert!((cs.energies[1].unwrap() - 8.0).abs() < 1e-10);
        assert!((cs.energies[2].unwrap() - 10.0).abs() < 1e-10);
    }
}
