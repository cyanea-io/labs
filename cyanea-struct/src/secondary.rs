//! Simplified secondary structure assignment.
//!
//! Uses hydrogen-bond geometry heuristics (simplified DSSP) to assign helix,
//! sheet, turn, and coil to each residue in a polypeptide chain.

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::geometry::dihedral_points;
use crate::types::{Chain, Point3D, Residue};

use alloc::format;
use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

/// Secondary structure classification for a single residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    Turn,
    Coil,
}

impl SecondaryStructure {
    /// Single-character DSSP code.
    pub fn code(&self) -> char {
        match self {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'E',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Coil => 'C',
        }
    }
}

/// Secondary structure assignment for a chain.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SecondaryStructureAssignment {
    /// Chain identifier.
    pub chain_id: char,
    /// One assignment per residue.
    pub assignments: Vec<SecondaryStructure>,
}

impl SecondaryStructureAssignment {
    /// Count of each secondary structure type.
    pub fn counts(&self) -> (usize, usize, usize, usize) {
        let mut h = 0;
        let mut e = 0;
        let mut t = 0;
        let mut c = 0;
        for ss in &self.assignments {
            match ss {
                SecondaryStructure::Helix => h += 1,
                SecondaryStructure::Sheet => e += 1,
                SecondaryStructure::Turn => t += 1,
                SecondaryStructure::Coil => c += 1,
            }
        }
        (h, e, t, c)
    }

    /// Fraction of residues in helix.
    pub fn helix_fraction(&self) -> f64 {
        if self.assignments.is_empty() {
            return 0.0;
        }
        let (h, _, _, _) = self.counts();
        h as f64 / self.assignments.len() as f64
    }

    /// Fraction of residues in sheet.
    pub fn sheet_fraction(&self) -> f64 {
        if self.assignments.is_empty() {
            return 0.0;
        }
        let (_, e, _, _) = self.counts();
        e as f64 / self.assignments.len() as f64
    }
}

impl Summarizable for SecondaryStructureAssignment {
    fn summary(&self) -> String {
        let (h, e, t, c) = self.counts();
        format!(
            "Chain {} SS: {} residue(s) — H:{} E:{} T:{} C:{}",
            self.chain_id,
            self.assignments.len(),
            h,
            e,
            t,
            c,
        )
    }
}

/// Assign secondary structure to each residue in a chain.
///
/// Uses a simplified DSSP-like algorithm:
/// - **Helix**: backbone O(i) to N(i+4) distance < 3.5 Å for at least 3 consecutive residues
/// - **Turn**: backbone O(i) to N(i+3) distance < 3.5 Å
/// - **Sheet**: detected via antiparallel beta-bridge patterns (O(i)→N(j) and O(j)→N(i) both < 3.5 Å with |i-j| > 4)
/// - **Coil**: everything else
pub fn assign_secondary_structure(chain: &Chain) -> Result<SecondaryStructureAssignment> {
    let n = chain.residues.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "cannot assign SS to empty chain".into(),
        ));
    }

    let mut assignments = vec![SecondaryStructure::Coil; n];

    // Precompute backbone atom positions (N, O) for each residue
    let backbone: Vec<Option<(Point3D, Point3D)>> = chain
        .residues
        .iter()
        .map(|r| {
            let n_atom = r.get_atom("N")?;
            let o_atom = r.get_atom("O")?;
            Some((n_atom.coords, o_atom.coords))
        })
        .collect();

    let hbond_cutoff = 3.5;

    // Check i→i+4 H-bonds (helix pattern)
    let mut helix_hbond = vec![false; n];
    for i in 0..n.saturating_sub(4) {
        if let (Some((_, o_i)), Some((n_i4, _))) = (&backbone[i], &backbone[i + 4]) {
            if o_i.distance_to(n_i4) < hbond_cutoff {
                helix_hbond[i] = true;
            }
        }
    }

    // Helix: 3+ consecutive i→i+4 H-bonds
    for i in 0..n.saturating_sub(6) {
        if helix_hbond[i] && helix_hbond[i + 1] && helix_hbond[i + 2] {
            // Mark residues i through i+6 (the full helix span) as Helix
            for j in i..=(i + 5).min(n - 1) {
                assignments[j] = SecondaryStructure::Helix;
            }
        }
    }

    // Check antiparallel beta-bridges (sheet pattern)
    for i in 0..n {
        for j in (i + 5)..n {
            if let (Some((n_i, o_i)), Some((n_j, o_j))) = (&backbone[i], &backbone[j]) {
                let oi_nj = o_i.distance_to(n_j);
                let oj_ni = o_j.distance_to(n_i);
                if oi_nj < hbond_cutoff && oj_ni < hbond_cutoff {
                    if assignments[i] == SecondaryStructure::Coil {
                        assignments[i] = SecondaryStructure::Sheet;
                    }
                    if assignments[j] == SecondaryStructure::Coil {
                        assignments[j] = SecondaryStructure::Sheet;
                    }
                }
            }
        }
    }

    // Check i→i+3 (turn pattern)
    for i in 0..n.saturating_sub(3) {
        if let (Some((_, o_i)), Some((n_i3, _))) = (&backbone[i], &backbone[i + 3]) {
            if o_i.distance_to(n_i3) < hbond_cutoff {
                // Only assign turn if not already helix/sheet
                for j in i..=(i + 3).min(n - 1) {
                    if assignments[j] == SecondaryStructure::Coil {
                        assignments[j] = SecondaryStructure::Turn;
                    }
                }
            }
        }
    }

    Ok(SecondaryStructureAssignment {
        chain_id: chain.id,
        assignments,
    })
}

/// Compute backbone dihedral angles (phi, psi) for a residue given its neighbors.
///
/// - `prev`: previous residue (needed for phi: C(i-1)–N(i)–CA(i)–C(i))
/// - `curr`: current residue
/// - `next`: next residue (needed for psi: N(i)–CA(i)–C(i)–N(i+1))
///
/// Returns `Some((phi, psi))` in degrees if all required atoms are present.
pub fn backbone_dihedrals(
    prev: Option<&Residue>,
    curr: &Residue,
    next: Option<&Residue>,
) -> Option<(f64, f64)> {
    let n_i = curr.get_atom("N")?.coords;
    let ca_i = curr.get_atom("CA")?.coords;
    let c_i = curr.get_atom("C")?.coords;

    let phi = if let Some(prev_res) = prev {
        let c_prev = prev_res.get_atom("C")?.coords;
        dihedral_points(&c_prev, &n_i, &ca_i, &c_i)
    } else {
        return None;
    };

    let psi = if let Some(next_res) = next {
        let n_next = next_res.get_atom("N")?.coords;
        dihedral_points(&n_i, &ca_i, &c_i, &n_next)
    } else {
        return None;
    };

    Some((phi, psi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Atom, Chain, Residue};

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

    fn make_helix_chain() -> Chain {
        // Build a synthetic alpha-helix with proper H-bond geometry.
        // Alpha helix: 3.6 residues/turn, rise 1.5Å/residue, radius 2.3Å.
        // Key: O(i) must be < 3.5Å from N(i+4).
        //
        // Strategy: first place all N atoms on the helix, then place O atoms
        // offset toward N(i+4) so the H-bond distance criterion is met.
        let n_residues = 12;
        let rise = 1.5_f64; // Å per residue
        let r = 2.3_f64; // helix radius
        let angle_per_res = 100.0_f64.to_radians(); // ~100° per residue

        // Compute N positions first
        let n_positions: Vec<Point3D> = (0..n_residues)
            .map(|i| {
                let angle = i as f64 * angle_per_res;
                Point3D::new(r * angle.cos(), r * angle.sin(), i as f64 * rise)
            })
            .collect();

        let mut residues = Vec::new();
        for i in 0..n_residues {
            let n_pos = n_positions[i];
            let angle = i as f64 * angle_per_res;
            let z = i as f64 * rise;

            let ca_pos = Point3D::new(
                r * (angle + 0.3).cos(),
                r * (angle + 0.3).sin(),
                z + 0.4,
            );
            let c_pos = Point3D::new(
                r * (angle + 0.6).cos(),
                r * (angle + 0.6).sin(),
                z + 0.8,
            );

            // Place O near N(i+4) so that dist(O_i, N_{i+4}) < 3.5Å.
            // We place O 2.5Å from C along the C→N(i+4) direction, giving
            // an O→N(i+4) distance of ~2.7Å (well within the 3.5Å cutoff).
            let o_pos = if i + 4 < n_residues {
                let target = n_positions[i + 4];
                let dir = target.sub(&c_pos);
                let d = dir.norm();
                if d > 1e-10 {
                    c_pos.add(&dir.scale(2.5 / d))
                } else {
                    Point3D::new(c_pos.x, c_pos.y + 1.24, c_pos.z)
                }
            } else {
                Point3D::new(c_pos.x, c_pos.y + 1.24, c_pos.z)
            };

            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i as i32 + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", n_pos.x, n_pos.y, n_pos.z),
                    make_atom("CA", ca_pos.x, ca_pos.y, ca_pos.z),
                    make_atom("C", c_pos.x, c_pos.y, c_pos.z),
                    make_atom("O", o_pos.x, o_pos.y, o_pos.z),
                ],
            });
        }
        Chain::new('A', residues)
    }

    #[test]
    fn short_chain_all_coil() {
        let chain = Chain::new(
            'A',
            vec![
                Residue {
                    name: "ALA".into(),
                    seq_num: 1,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 0.0, 0.0, 0.0),
                        make_atom("CA", 1.5, 0.0, 0.0),
                        make_atom("C", 3.0, 0.0, 0.0),
                        make_atom("O", 3.5, 1.0, 0.0),
                    ],
                },
                Residue {
                    name: "GLY".into(),
                    seq_num: 2,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 4.0, 0.0, 0.0),
                        make_atom("CA", 5.5, 0.0, 0.0),
                        make_atom("C", 7.0, 0.0, 0.0),
                        make_atom("O", 7.5, 1.0, 0.0),
                    ],
                },
            ],
        );
        let ss = assign_secondary_structure(&chain).unwrap();
        assert_eq!(ss.assignments.len(), 2);
        // Too short for helix/sheet, should be coil or turn
        for a in &ss.assignments {
            assert!(
                *a == SecondaryStructure::Coil || *a == SecondaryStructure::Turn,
                "short chain should be coil/turn"
            );
        }
    }

    #[test]
    fn helix_detection() {
        let chain = make_helix_chain();
        let ss = assign_secondary_structure(&chain).unwrap();
        // At least some residues should be helix
        let (h, _, _, _) = ss.counts();
        assert!(
            h > 0,
            "helix chain should have helical residues, got: {:?}",
            ss.assignments
        );
    }

    #[test]
    fn sheet_detection() {
        // Build two strands running antiparallel with close O-N contacts.
        let mut residues = Vec::new();
        // Strand 1: residues 1-5, running along +x
        for i in 0..5 {
            let x = i as f64 * 3.5;
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", x, 0.0, 0.0),
                    make_atom("CA", x + 1.0, 0.0, 0.0),
                    make_atom("C", x + 2.0, 0.0, 0.0),
                    make_atom("O", x + 2.0, 1.0, 0.0),
                ],
            });
        }
        // Spacer residues 6-10 (far away)
        for i in 5..10 {
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", 50.0, 50.0, i as f64 * 10.0),
                    make_atom("CA", 51.0, 50.0, i as f64 * 10.0),
                    make_atom("C", 52.0, 50.0, i as f64 * 10.0),
                    make_atom("O", 52.0, 51.0, i as f64 * 10.0),
                ],
            });
        }
        // Strand 2: residues 11-15, antiparallel (running along -x)
        // O of strand1[i] close to N of strand2[4-i] and vice versa
        for i in 0..5 {
            let x = (4 - i) as f64 * 3.5;
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: (i + 11) as i32,
                i_code: None,
                atoms: vec![
                    make_atom("N", x + 2.0, 2.5, 0.0), // Close to O of strand 1
                    make_atom("CA", x + 1.0, 2.5, 0.0),
                    make_atom("C", x, 2.5, 0.0),
                    make_atom("O", x, 1.5, 0.0), // Close to N of strand 1
                ],
            });
        }

        let chain = Chain::new('A', residues);
        let ss = assign_secondary_structure(&chain).unwrap();
        let (_, e, _, _) = ss.counts();
        assert!(e > 0, "should detect sheet residues, got: {:?}", ss.assignments);
    }

    #[test]
    fn backbone_dihedral_computation() {
        let prev = Residue {
            name: "ALA".into(),
            seq_num: 1,
            i_code: None,
            atoms: vec![
                make_atom("N", -1.0, 0.0, 0.0),
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("C", 1.0, 0.0, 0.0),
                make_atom("O", 1.0, 1.0, 0.0),
            ],
        };
        let curr = Residue {
            name: "GLY".into(),
            seq_num: 2,
            i_code: None,
            atoms: vec![
                make_atom("N", 2.0, 0.0, 0.0),
                make_atom("CA", 3.0, 0.0, 0.0),
                make_atom("C", 4.0, 0.0, 0.0),
                make_atom("O", 4.0, 1.0, 0.0),
            ],
        };
        let next = Residue {
            name: "ALA".into(),
            seq_num: 3,
            i_code: None,
            atoms: vec![
                make_atom("N", 5.0, 0.0, 0.0),
                make_atom("CA", 6.0, 0.0, 0.0),
                make_atom("C", 7.0, 0.0, 0.0),
                make_atom("O", 7.0, 1.0, 0.0),
            ],
        };
        let result = backbone_dihedrals(Some(&prev), &curr, Some(&next));
        assert!(result.is_some(), "should compute phi/psi");
    }
}
