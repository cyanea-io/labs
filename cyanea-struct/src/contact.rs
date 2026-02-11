//! Residue-residue contact maps.
//!
//! Provides CA-CA distance matrices and all-atom (minimum heavy-atom distance)
//! contact maps for protein chains.

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::types::{Chain, Point3D};

use alloc::format;
use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

/// A symmetric distance matrix for residue-residue contacts.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ContactMap {
    /// Chain identifier.
    pub chain_id: char,
    /// Number of residues (matrix is size × size).
    pub size: usize,
    /// Row-major n×n distance matrix.
    pub distances: Vec<f64>,
}

impl ContactMap {
    /// Get the distance between residues i and j.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.distances[i * self.size + j]
    }

    /// Count residue pairs with distance below the cutoff (excluding diagonal).
    pub fn count_contacts(&self, cutoff: f64) -> usize {
        let mut count = 0;
        for i in 0..self.size {
            for j in (i + 1)..self.size {
                if self.get(i, j) < cutoff {
                    count += 1;
                }
            }
        }
        count
    }

    /// Return (i, j) pairs of residues in contact below the cutoff.
    pub fn contacts_below(&self, cutoff: f64) -> Vec<(usize, usize)> {
        let mut contacts = Vec::new();
        for i in 0..self.size {
            for j in (i + 1)..self.size {
                if self.get(i, j) < cutoff {
                    contacts.push((i, j));
                }
            }
        }
        contacts
    }

    /// Contact density: fraction of possible pairs below the cutoff.
    pub fn contact_density(&self, cutoff: f64) -> f64 {
        if self.size < 2 {
            return 0.0;
        }
        let total_pairs = self.size * (self.size - 1) / 2;
        self.count_contacts(cutoff) as f64 / total_pairs as f64
    }
}

impl Summarizable for ContactMap {
    fn summary(&self) -> String {
        let contacts_8 = self.count_contacts(8.0);
        format!(
            "ContactMap chain {} — {} residues, {} contacts (<8Å)",
            self.chain_id, self.size, contacts_8,
        )
    }
}

/// Compute a CA-CA distance contact map for a chain.
///
/// Each entry (i, j) is the Euclidean distance between the alpha-carbon atoms
/// of residues i and j. Residues without a CA atom are assigned `f64::INFINITY`.
pub fn compute_contact_map(chain: &Chain) -> Result<ContactMap> {
    let n = chain.residues.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "cannot compute contact map for empty chain".into(),
        ));
    }

    // Extract CA positions
    let ca_positions: Vec<Option<Point3D>> = chain
        .residues
        .iter()
        .map(|r| r.get_alpha_carbon().map(|a| a.coords))
        .collect();

    let mut distances = vec![0.0f64; n * n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dist = match (&ca_positions[i], &ca_positions[j]) {
                (Some(pi), Some(pj)) => pi.distance_to(pj),
                _ => f64::INFINITY,
            };
            distances[i * n + j] = dist;
            distances[j * n + i] = dist;
        }
    }

    Ok(ContactMap {
        chain_id: chain.id,
        size: n,
        distances,
    })
}

/// Compute an all-atom contact map for a chain.
///
/// Each entry (i, j) is the minimum distance between any pair of non-hydrogen
/// atoms in residues i and j. This gives a tighter contact definition than CA-CA.
pub fn compute_contact_map_allatom(chain: &Chain) -> Result<ContactMap> {
    let n = chain.residues.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "cannot compute contact map for empty chain".into(),
        ));
    }

    #[cfg(feature = "parallel")]
    let distances = {
        use rayon::prelude::*;
        let upper: Vec<Vec<(usize, f64)>> = (0..n)
            .into_par_iter()
            .map(|i| {
                ((i + 1)..n)
                    .map(|j| {
                        let mut min_dist = f64::INFINITY;
                        for a1 in &chain.residues[i].atoms {
                            if a1.element.as_deref() == Some("H") {
                                continue;
                            }
                            for a2 in &chain.residues[j].atoms {
                                if a2.element.as_deref() == Some("H") {
                                    continue;
                                }
                                let d = a1.coords.distance_to(&a2.coords);
                                if d < min_dist {
                                    min_dist = d;
                                }
                            }
                        }
                        (j, min_dist)
                    })
                    .collect()
            })
            .collect();
        let mut distances = vec![0.0f64; n * n];
        for (i, row) in upper.into_iter().enumerate() {
            for (j, d) in row {
                distances[i * n + j] = d;
                distances[j * n + i] = d;
            }
        }
        distances
    };

    #[cfg(not(feature = "parallel"))]
    let distances = {
        let mut distances = vec![0.0f64; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let mut min_dist = f64::INFINITY;
                for a1 in &chain.residues[i].atoms {
                    if a1.element.as_deref() == Some("H") {
                        continue;
                    }
                    for a2 in &chain.residues[j].atoms {
                        if a2.element.as_deref() == Some("H") {
                            continue;
                        }
                        let d = a1.coords.distance_to(&a2.coords);
                        if d < min_dist {
                            min_dist = d;
                        }
                    }
                }
                distances[i * n + j] = min_dist;
                distances[j * n + i] = min_dist;
            }
        }
        distances
    };

    Ok(ContactMap {
        chain_id: chain.id,
        size: n,
        distances,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Atom, Chain, Point3D, Residue};

    fn make_atom(name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial: 1,
            name: name.into(),
            alt_loc: None,
            coords: Point3D::new(x, y, z),
            occupancy: 1.0,
            temp_factor: 0.0,
            element: Some("C".into()),
            charge: None,
            is_hetatm: false,
        }
    }

    fn make_test_chain() -> Chain {
        Chain::new(
            'A',
            vec![
                Residue {
                    name: "ALA".into(),
                    seq_num: 1,
                    i_code: None,
                    atoms: vec![make_atom("CA", 0.0, 0.0, 0.0)],
                },
                Residue {
                    name: "GLY".into(),
                    seq_num: 2,
                    i_code: None,
                    atoms: vec![make_atom("CA", 3.8, 0.0, 0.0)],
                },
                Residue {
                    name: "VAL".into(),
                    seq_num: 3,
                    i_code: None,
                    atoms: vec![make_atom("CA", 7.6, 0.0, 0.0)],
                },
                Residue {
                    name: "LEU".into(),
                    seq_num: 4,
                    i_code: None,
                    atoms: vec![make_atom("CA", 11.4, 0.0, 0.0)],
                },
            ],
        )
    }

    #[test]
    fn ca_contact_map() {
        let chain = make_test_chain();
        let cm = compute_contact_map(&chain).unwrap();
        assert_eq!(cm.size, 4);
        // Diagonal should be zero
        assert!((cm.get(0, 0)).abs() < 1e-10);
        // Distance between residues 0 and 1 should be 3.8
        assert!((cm.get(0, 1) - 3.8).abs() < 1e-10);
        // Symmetric
        assert!((cm.get(0, 1) - cm.get(1, 0)).abs() < 1e-10);
    }

    #[test]
    fn get_distance() {
        let chain = make_test_chain();
        let cm = compute_contact_map(&chain).unwrap();
        assert!((cm.get(0, 2) - 7.6).abs() < 1e-10);
    }

    #[test]
    fn count_contacts_cutoff() {
        let chain = make_test_chain();
        let cm = compute_contact_map(&chain).unwrap();
        // With cutoff 8.0: (0,1)=3.8, (0,2)=7.6, (1,2)=3.8, (1,3)=7.6, (2,3)=3.8 → 5 contacts
        assert_eq!(cm.count_contacts(8.0), 5);
        // With cutoff 4.0: only adjacent pairs 3.8 < 4.0
        assert_eq!(cm.count_contacts(4.0), 3);
    }

    #[test]
    fn contact_density() {
        let chain = make_test_chain();
        let cm = compute_contact_map(&chain).unwrap();
        // 4 residues = 6 pairs. 5 contacts at 8Å cutoff → density = 5/6
        let density = cm.contact_density(8.0);
        assert!((density - 5.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn allatom_contact_map() {
        let chain = Chain::new(
            'A',
            vec![
                Residue {
                    name: "ALA".into(),
                    seq_num: 1,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 0.0, 0.0, 0.0),
                        make_atom("CA", 1.0, 0.0, 0.0),
                        make_atom("C", 2.0, 0.0, 0.0),
                    ],
                },
                Residue {
                    name: "GLY".into(),
                    seq_num: 2,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 3.0, 0.0, 0.0),
                        make_atom("CA", 4.0, 0.0, 0.0),
                        make_atom("C", 5.0, 0.0, 0.0),
                    ],
                },
            ],
        );
        let cm = compute_contact_map_allatom(&chain).unwrap();
        // Min distance: C of res1 (2.0) to N of res2 (3.0) = 1.0
        assert!((cm.get(0, 1) - 1.0).abs() < 1e-10);
    }
}
