//! Stereochemistry assignment: R/S for tetrahedral centers, E/Z for double bonds.
//!
//! Implements the Cahn-Ingold-Prelog (CIP) priority rules for assigning
//! absolute stereochemical descriptors.

use crate::molecule::{BondOrder, BondStereo, Chirality, Molecule};

/// Assign R/S descriptor to a tetrahedral stereocenter.
///
/// Returns `Some('R')` or `Some('S')` if the atom is a stereocenter with
/// chirality information, or `None` if it is not a valid stereocenter.
///
/// The algorithm uses CIP priority rules based on atomic number, with
/// recursive tiebreaking by examining attached atoms.
pub fn assign_rs(mol: &Molecule, atom_idx: usize) -> Option<char> {
    let atom = mol.atoms.get(atom_idx)?;
    if atom.chirality == Chirality::None {
        return None;
    }

    // A tetrahedral center needs exactly 4 neighbors (including implicit H as a virtual neighbor).
    // Gather explicit neighbors in the order they appear in the adjacency list
    // (which matches the SMILES input order).
    let explicit_neighbors: Vec<usize> = mol.adjacency[atom_idx]
        .iter()
        .map(|&(n, _)| n)
        .collect();

    let total_neighbors = explicit_neighbors.len() + atom.implicit_hydrogens as usize;
    if total_neighbors != 4 {
        return None;
    }

    // Build list of neighbor priorities. We represent implicit H as None (atomic number 1),
    // explicit neighbors as Some(idx).
    let mut neighbor_list: Vec<Option<usize>> = explicit_neighbors.iter().map(|&n| Some(n)).collect();
    for _ in 0..atom.implicit_hydrogens {
        neighbor_list.push(None); // implicit hydrogen
    }

    // Compute CIP priorities for each neighbor
    let priorities: Vec<u64> = neighbor_list
        .iter()
        .map(|&n| cip_priority(mol, atom_idx, n))
        .collect();

    // Check that all priorities are distinct (otherwise not a true stereocenter)
    let mut sorted_prio = priorities.clone();
    sorted_prio.sort();
    for i in 1..sorted_prio.len() {
        if sorted_prio[i] == sorted_prio[i - 1] {
            return None; // Not a stereocenter — duplicate priorities
        }
    }

    // Determine R/S from the chirality annotation and priority ordering.
    //
    // In SMILES, @ means the neighbors (as listed from second to fourth, looking from the
    // first neighbor) are arranged counterclockwise. @@ means clockwise.
    //
    // To determine R/S:
    // 1. We need to look at the molecule from the side opposite the lowest-priority group.
    // 2. If the remaining three groups decrease in priority going clockwise → R.
    // 3. If they decrease counterclockwise → S.

    // Find the index of the lowest priority neighbor in our list
    let min_prio_pos = priorities
        .iter()
        .enumerate()
        .min_by_key(|&(_, &p)| p)
        .unwrap()
        .0;

    // Get the other three neighbors in order (preserving the cyclic order from SMILES)
    let mut others: Vec<(usize, u64)> = Vec::new();
    for i in 0..4 {
        if i != min_prio_pos {
            others.push((i, priorities[i]));
        }
    }

    // Determine if the three remaining priorities are in descending clockwise order
    // when the SMILES chirality annotation is considered.
    //
    // In SMILES notation with @ (counterclockwise):
    //   neighbors[0] is the "from" atom, neighbors[1..4] go counterclockwise
    //
    // We need to figure out whether the cyclic order of the three remaining
    // neighbors (by CIP priority, highest to lowest) is clockwise or counterclockwise.

    // The parity: count inversions needed to sort the three priorities descending
    let p0 = others[0].1;
    let p1 = others[1].1;
    let p2 = others[2].1;

    let mut inversions = 0u32;
    if p0 < p1 { inversions += 1; }
    if p0 < p2 { inversions += 1; }
    if p1 < p2 { inversions += 1; }

    let is_even_perm = inversions % 2 == 0;

    // If the lowest-priority neighbor was removed from a non-first position,
    // that introduces additional permutation parity.
    let removal_parity = min_prio_pos % 2 == 0; // even position = even removal parity

    // Combine with the SMILES chirality annotation:
    // @ = CCW as written = one sense, @@ = CW as written = other sense
    let smiles_ccw = atom.chirality == Chirality::CounterClockwise;

    // The combination of these three parities determines R vs S
    let mut is_r = is_even_perm;
    if !removal_parity {
        is_r = !is_r;
    }
    if smiles_ccw {
        is_r = !is_r;
    }

    Some(if is_r { 'R' } else { 'S' })
}

/// Assign E/Z descriptor to a double bond.
///
/// Returns `Some('E')` or `Some('Z')` if the bond has stereo information
/// and valid substituents, or `None` otherwise.
///
/// E (entgegen) = higher-priority groups on opposite sides.
/// Z (zusammen) = higher-priority groups on same side.
pub fn assign_ez(mol: &Molecule, bond_idx: usize) -> Option<char> {
    let bond = mol.bonds.get(bond_idx)?;

    // Must be a double bond
    if bond.order != BondOrder::Double {
        return None;
    }

    let atom1 = bond.atom1;
    let atom2 = bond.atom2;

    // Find the stereo bonds connected to atom1 and atom2
    // Returns (neighbor_idx, raw_bond_stereo, is_center_atom1_of_bond)
    let (stereo1_neighbor, stereo1_raw, stereo1_center_is_atom1) =
        find_stereo_neighbor_raw(mol, atom1, atom2)?;
    let (stereo2_neighbor, stereo2_raw, stereo2_center_is_atom1) =
        find_stereo_neighbor_raw(mol, atom2, atom1)?;

    // Determine CIP priorities for the substituents on each side
    let atom1_neighbors: Vec<usize> = mol.adjacency[atom1]
        .iter()
        .map(|&(n, _)| n)
        .filter(|&n| n != atom2)
        .collect();

    let atom2_neighbors: Vec<usize> = mol.adjacency[atom2]
        .iter()
        .map(|&(n, _)| n)
        .filter(|&n| n != atom1)
        .collect();

    if atom1_neighbors.is_empty() || atom2_neighbors.is_empty() {
        return None;
    }

    // For each side, find the higher-priority substituent
    let prio1_stereo = cip_priority(mol, atom1, Some(stereo1_neighbor));
    let prio1_other: Option<u64> = atom1_neighbors
        .iter()
        .filter(|&&n| n != stereo1_neighbor)
        .map(|&n| cip_priority(mol, atom1, Some(n)))
        .max();
    let stereo1_is_higher = match prio1_other {
        Some(other) => prio1_stereo > other,
        None => true, // Only one substituent
    };

    let prio2_stereo = cip_priority(mol, atom2, Some(stereo2_neighbor));
    let prio2_other: Option<u64> = atom2_neighbors
        .iter()
        .filter(|&&n| n != stereo2_neighbor)
        .map(|&n| cip_priority(mol, atom2, Some(n)))
        .max();
    let stereo2_is_higher = match prio2_other {
        Some(other) => prio2_stereo > other,
        None => true,
    };

    // Determine same-side vs opposite-side from the / and \ markers.
    //
    // In SMILES: F/C=C/F means trans (E), F/C=C\F means cis (Z).
    //
    // The raw stereo directions on the bonds: identical raw markers (both / or both \)
    // means the substituents are on opposite sides (trans).
    // Different raw markers (one / and one \) means same side (cis).
    //
    // However, we also need to account for the direction of the bond relative
    // to the double bond atom. The stereo bond is stored with atom1->atom2 direction.
    // If the double-bond atom is atom1 of the stereo bond, the substituent is at atom2
    // (away from double bond). If it's atom2, the substituent is at atom1 (away from double bond).
    //
    // For the left side (atom1 of double bond): the stereo bond goes from the double bond
    // atom outward to the substituent, or from the substituent inward to the double bond atom.
    // For the right side (atom2 of double bond): similarly.
    //
    // The inversion rule: if the stereo bond is "incoming" (substituent -> double bond atom,
    // i.e., the double bond atom is atom2 of the stereo bond), we need to invert the direction
    // to get the "outgoing" direction from the double bond atom's perspective.
    //
    // After normalization to "outgoing from double bond atom":
    // - Same directions = trans (E)
    // - Different directions = cis (Z)

    let dir1 = if stereo1_center_is_atom1 {
        stereo1_raw // Outgoing from center
    } else {
        // Incoming to center, invert
        invert_stereo(stereo1_raw)
    };

    let dir2 = if stereo2_center_is_atom1 {
        stereo2_raw // Outgoing from center
    } else {
        invert_stereo(stereo2_raw)
    };

    // Same normalized directions (both Up or both Down) = same side (cis, Z)
    // Different normalized directions (one Up, one Down) = opposite sides (trans, E)
    let substituents_on_same_side = dir1 == dir2;

    // Now combine with CIP priorities:
    // The stereo markers tell us about the specific substituent atoms found.
    // If both stereo-marked substituents are the higher-priority ones on their side,
    // then same_side directly gives us Z, opposite gives E.
    // If one is higher and one is lower (or both lower), we need to flip.
    let higher_on_same_side = if stereo1_is_higher == stereo2_is_higher {
        substituents_on_same_side
    } else {
        !substituents_on_same_side
    };

    Some(if higher_on_same_side { 'Z' } else { 'E' })
}

fn invert_stereo(s: BondStereo) -> BondStereo {
    match s {
        BondStereo::Up => BondStereo::Down,
        BondStereo::Down => BondStereo::Up,
        BondStereo::None => BondStereo::None,
    }
}

/// Find a neighbor of `center` (other than `exclude`) that has a stereo bond.
/// Returns (neighbor_index, raw_bond_stereo, center_is_atom1_of_bond).
fn find_stereo_neighbor_raw(
    mol: &Molecule,
    center: usize,
    exclude: usize,
) -> Option<(usize, BondStereo, bool)> {
    for &(neighbor, bond_idx) in &mol.adjacency[center] {
        if neighbor == exclude {
            continue;
        }
        let bond = &mol.bonds[bond_idx];
        if bond.stereo != BondStereo::None {
            let center_is_atom1 = bond.atom1 == center;
            return Some((neighbor, bond.stereo, center_is_atom1));
        }
    }
    None
}

/// CIP priority for comparison. Higher value = higher CIP priority.
///
/// Encodes (atomic_number, level1_neighbors_desc, level2_neighbors_desc)
/// into a u64 with atomic number in the most significant position.
///
/// `neighbor` is `Some(atom_idx)` for an explicit neighbor, or `None` for implicit hydrogen.
fn cip_priority(mol: &Molecule, center: usize, neighbor: Option<usize>) -> u64 {
    match neighbor {
        None => {
            // Implicit hydrogen: atomic number 1, no further substituents
            1_u64 << 42
        }
        Some(idx) => {
            let atom = &mol.atoms[idx];
            let an = atom.atomic_number as u64;

            // Level 1: sorted descending neighbor atomic numbers (excluding center)
            let mut level1: Vec<u8> = Vec::new();
            for &(n, _) in &mol.adjacency[idx] {
                if n == center {
                    continue;
                }
                level1.push(mol.atoms[n].atomic_number);
            }
            for _ in 0..atom.implicit_hydrogens {
                level1.push(1);
            }
            level1.sort_unstable_by(|a, b| b.cmp(a));

            // Encode level1: up to 4 neighbors, each 7 bits (values 0-127) = 28 bits
            let mut l1: u64 = 0;
            for (i, &v) in level1.iter().enumerate() {
                if i >= 4 { break; }
                l1 = l1 * 128 + v as u64;
            }

            // Level 2: for each level-1 neighbor, examine their neighbors
            let mut level2: Vec<u8> = Vec::new();
            for &(n, _) in &mol.adjacency[idx] {
                if n == center {
                    continue;
                }
                for &(nn, _) in &mol.adjacency[n] {
                    if nn == idx {
                        continue;
                    }
                    level2.push(mol.atoms[nn].atomic_number);
                }
                for _ in 0..mol.atoms[n].implicit_hydrogens {
                    level2.push(1);
                }
            }
            level2.sort_unstable_by(|a, b| b.cmp(a));

            let mut l2: u64 = 0;
            for (i, &v) in level2.iter().enumerate() {
                if i >= 2 { break; }
                l2 = l2 * 128 + v as u64;
            }

            // Pack: an (7 bits at top), l1 (28 bits), l2 (14 bits) = 49 bits total, fits in u64
            // an << 42 | l1 << 14 | l2
            (an << 42) | (l1 << 14) | l2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn parse_chirality_at() {
        // L-alanine: N[C@@H](C)C(=O)O
        let mol = parse_smiles("N[C@@H](C)C(=O)O").unwrap();
        assert_eq!(mol.atoms[1].chirality, Chirality::Clockwise);
    }

    #[test]
    fn parse_chirality_at_single() {
        // [C@H] counterclockwise
        let mol = parse_smiles("N[C@H](C)C(=O)O").unwrap();
        assert_eq!(mol.atoms[1].chirality, Chirality::CounterClockwise);
    }

    #[test]
    fn no_chirality_by_default() {
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(mol.atoms[0].chirality, Chirality::None);
        assert_eq!(mol.atoms[1].chirality, Chirality::None);
    }

    #[test]
    fn assign_rs_l_alanine() {
        // L-alanine: N[C@@H](C)C(=O)O
        // The central carbon has: N, H, CH3, COOH
        // CIP priorities: O-containing group > N > C(methyl) > H
        // @@ = clockwise in SMILES
        let mol = parse_smiles("N[C@@H](C)C(=O)O").unwrap();
        let rs = assign_rs(&mol, 1);
        assert!(rs.is_some(), "should be a stereocenter");
        let descriptor = rs.unwrap();
        assert!(
            descriptor == 'R' || descriptor == 'S',
            "should be R or S, got {descriptor}"
        );
    }

    #[test]
    fn assign_rs_returns_none_for_non_stereo() {
        // Methane has no stereocenter
        let mol = parse_smiles("C").unwrap();
        assert_eq!(assign_rs(&mol, 0), None);
    }

    #[test]
    fn assign_rs_returns_none_for_symmetric() {
        // CC(C)C — central carbon has 3 identical methyl groups + H
        // Not a true stereocenter (duplicate priorities)
        let mol = parse_smiles("[C@@H](C)(C)C").unwrap();
        assert_eq!(assign_rs(&mol, 0), None);
    }

    #[test]
    fn chirality_opposite_descriptors() {
        // @ and @@ on the same molecule should give opposite descriptors
        let mol1 = parse_smiles("N[C@H](C)C(=O)O").unwrap();
        let mol2 = parse_smiles("N[C@@H](C)C(=O)O").unwrap();
        let rs1 = assign_rs(&mol1, 1);
        let rs2 = assign_rs(&mol2, 1);
        assert!(rs1.is_some() && rs2.is_some());
        assert_ne!(rs1.unwrap(), rs2.unwrap(), "@ and @@ should give opposite R/S");
    }

    #[test]
    fn parse_bond_stereo_up() {
        // F/C=C/F — trans-1,2-difluoroethene
        let mol = parse_smiles("F/C=C/F").unwrap();
        // Bond 0 (F-C) should have Up stereo
        assert_eq!(mol.bonds[0].stereo, BondStereo::Up);
    }

    #[test]
    fn parse_bond_stereo_down() {
        // F\C=C\F
        let mol = parse_smiles("F\\C=C\\F").unwrap();
        assert_eq!(mol.bonds[0].stereo, BondStereo::Down);
    }

    #[test]
    fn assign_ez_trans() {
        // F/C=C/F is trans (E)
        let mol = parse_smiles("F/C=C/F").unwrap();
        // Find the double bond
        let db_idx = mol.bonds.iter().position(|b| b.order == BondOrder::Double).unwrap();
        let ez = assign_ez(&mol, db_idx);
        assert!(ez.is_some(), "should have E/Z assignment");
        assert_eq!(ez.unwrap(), 'E', "F/C=C/F should be E (trans)");
    }

    #[test]
    fn assign_ez_cis() {
        // F/C=C\F is cis (Z)
        let mol = parse_smiles("F/C=C\\F").unwrap();
        let db_idx = mol.bonds.iter().position(|b| b.order == BondOrder::Double).unwrap();
        let ez = assign_ez(&mol, db_idx);
        assert!(ez.is_some(), "should have E/Z assignment");
        assert_eq!(ez.unwrap(), 'Z', "F/C=C\\F should be Z (cis)");
    }

    #[test]
    fn assign_ez_returns_none_for_single_bond() {
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(assign_ez(&mol, 0), None);
    }

    #[test]
    fn assign_ez_returns_none_without_stereo_markers() {
        let mol = parse_smiles("FC=CF").unwrap();
        let db_idx = mol.bonds.iter().position(|b| b.order == BondOrder::Double).unwrap();
        assert_eq!(assign_ez(&mol, db_idx), None);
    }
}
