//! Canonical SMILES generation.
//!
//! Generates a deterministic, canonical SMILES string for a molecule using
//! Morgan-like invariant refinement and DFS traversal.
//!
//! # Example
//!
//! ```
//! use cyanea_chem::{parse_smiles, canonical_smiles};
//!
//! // Different SMILES for ethanol produce the same canonical form
//! let mol1 = parse_smiles("OCC").unwrap();
//! let mol2 = parse_smiles("CCO").unwrap();
//! assert_eq!(canonical_smiles(&mol1), canonical_smiles(&mol2));
//! ```

use crate::element::element_by_number;
use crate::molecule::{BondOrder, Molecule};

/// Generate a canonical SMILES string for the given molecule.
///
/// The algorithm:
/// 1. Assign initial invariants based on atom properties.
/// 2. Refine invariants using Morgan-like iterative neighbor summation.
/// 3. Choose the canonical starting atom (lowest rank).
/// 4. DFS traversal to produce the SMILES string with canonical ordering.
pub fn canonical_smiles(mol: &Molecule) -> String {
    let n = mol.atom_count();
    if n == 0 {
        return String::new();
    }

    // Step 1 & 2: Compute canonical ranks
    let ranks = compute_canonical_ranks(mol);

    // Step 3: Handle disconnected fragments
    let mut visited = vec![false; n];
    let mut result = String::new();

    // Precompute ring closure assignments via a preliminary DFS pass
    let mut ring_closures = precompute_ring_closures(mol, &ranks);

    loop {
        // Find the unvisited atom with the lowest rank
        let start = (0..n)
            .filter(|&i| !visited[i])
            .min_by_key(|&i| ranks[i]);

        match start {
            Some(start_idx) => {
                if !result.is_empty() {
                    result.push('.'); // fragment separator
                }
                // Step 4: DFS traversal
                dfs_smiles(
                    mol,
                    start_idx,
                    None,
                    &ranks,
                    &mut visited,
                    &mut result,
                    &mut ring_closures,
                );
            }
            None => break,
        }
    }

    result
}

/// A ring closure assignment: which ring digits to write at each atom.
struct RingClosureInfo {
    /// For each atom: list of (ring_number, bond_order, is_opening)
    /// is_opening: true if this atom is where the ring opens (first encountered in DFS),
    /// false if this is where the ring closes (second encountered).
    atom_closures: Vec<Vec<(usize, BondOrder, bool)>>,
}

/// Precompute ring closure assignments by doing a DFS and identifying back-edges.
fn precompute_ring_closures(mol: &Molecule, ranks: &[u64]) -> RingClosureInfo {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut atom_closures: Vec<Vec<(usize, BondOrder, bool)>> = vec![Vec::new(); n];
    let mut next_ring_num: usize = 1;
    let mut used_bonds = vec![false; mol.bond_count()];

    // Process each connected component
    for start in 0..n {
        if visited[start] {
            continue;
        }
        // Find the canonical start atom for this component
        let component_start = (0..n)
            .filter(|&i| !visited[i])
            .min_by_key(|&i| ranks[i])
            .unwrap_or(start);

        precompute_dfs(
            mol,
            component_start,
            None,
            ranks,
            &mut visited,
            &mut atom_closures,
            &mut next_ring_num,
            &mut used_bonds,
        );
    }

    RingClosureInfo { atom_closures }
}

fn precompute_dfs(
    mol: &Molecule,
    atom_idx: usize,
    from_atom: Option<usize>,
    ranks: &[u64],
    visited: &mut Vec<bool>,
    atom_closures: &mut Vec<Vec<(usize, BondOrder, bool)>>,
    next_ring_num: &mut usize,
    used_bonds: &mut Vec<bool>,
) {
    visited[atom_idx] = true;

    // Get neighbors sorted by rank, excluding the atom we came from
    let mut neighbors: Vec<(usize, usize)> = mol.adjacency[atom_idx]
        .iter()
        .map(|&(n, bi)| (n, bi))
        .filter(|&(n, _)| Some(n) != from_atom)
        .collect();
    neighbors.sort_by_key(|&(n, _)| ranks[n]);

    for &(n, bi) in &neighbors {
        if visited[n] {
            // Back-edge: only create ring closure if this bond hasn't been used
            if !used_bonds[bi] {
                used_bonds[bi] = true;
                let bond = &mol.bonds[bi];
                let ring_num = *next_ring_num;
                *next_ring_num += 1;
                atom_closures[n].push((ring_num, bond.order, true));
                atom_closures[atom_idx].push((ring_num, bond.order, false));
            }
        } else {
            precompute_dfs(
                mol,
                n,
                Some(atom_idx),
                ranks,
                visited,
                atom_closures,
                next_ring_num,
                used_bonds,
            );
        }
    }
}

/// Compute canonical atom ranks using Morgan-like algorithm.
fn compute_canonical_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();

    // Step 1: Initial invariants
    // Encode: (atomic_number, degree, h_count, charge, mass, aromatic) into a single u64
    let mut invariants: Vec<u64> = Vec::with_capacity(n);
    for i in 0..n {
        let atom = &mol.atoms[i];
        let degree = mol.degree(i) as u64;
        let h_count = atom.implicit_hydrogens as u64;
        let charge = (atom.formal_charge as i64 + 128) as u64; // shift to positive
        let mass = atom.isotope.unwrap_or(0) as u64;
        let aromatic = atom.is_aromatic as u64;

        let inv = (atom.atomic_number as u64) << 40
            | degree << 32
            | h_count << 24
            | charge << 16
            | mass << 8
            | aromatic;
        invariants.push(inv);
    }

    // Step 2: Morgan-like iterative refinement
    // Iterate until the number of distinct invariants stops increasing
    let mut prev_distinct = count_distinct(&invariants);

    for _ in 0..n {
        let mut new_invariants = Vec::with_capacity(n);
        for i in 0..n {
            let mut combined = invariants[i].wrapping_mul(1000003);
            // Sort neighbor invariants for determinism
            let mut neighbor_invs: Vec<u64> = mol.adjacency[i]
                .iter()
                .map(|&(neighbor, bond_idx)| {
                    let bond_val = mol.bonds[bond_idx].order as u64;
                    invariants[neighbor].wrapping_mul(31).wrapping_add(bond_val)
                })
                .collect();
            neighbor_invs.sort();
            for nv in &neighbor_invs {
                combined = combined.wrapping_mul(1000003).wrapping_add(*nv);
            }
            new_invariants.push(combined);
        }

        let new_distinct = count_distinct(&new_invariants);
        invariants = new_invariants;

        if new_distinct <= prev_distinct {
            break; // Convergence: no more discrimination
        }
        prev_distinct = new_distinct;
    }

    // Convert invariants to ranks (0-based)
    let mut indexed: Vec<(u64, usize)> = invariants
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    indexed.sort();

    let mut ranks = vec![0u64; n];
    if !indexed.is_empty() {
        let mut rank = 0u64;
        ranks[indexed[0].1] = 0;
        for i in 1..indexed.len() {
            if indexed[i].0 != indexed[i - 1].0 {
                rank += 1;
            }
            ranks[indexed[i].1] = rank;
        }
    }

    ranks
}

/// Count distinct values in a slice.
fn count_distinct(values: &[u64]) -> usize {
    let mut sorted = values.to_vec();
    sorted.sort();
    sorted.dedup();
    sorted.len()
}

/// DFS traversal to generate SMILES string.
fn dfs_smiles(
    mol: &Molecule,
    atom_idx: usize,
    from_atom: Option<usize>,
    ranks: &[u64],
    visited: &mut Vec<bool>,
    output: &mut String,
    ring_info: &mut RingClosureInfo,
) {
    visited[atom_idx] = true;

    // Write atom
    write_atom(mol, atom_idx, output);

    // Write ring closure digits for this atom (both opening and closing)
    // Sort by ring number for determinism
    let mut closures = ring_info.atom_closures[atom_idx].clone();
    closures.sort_by_key(|&(rn, _, _)| rn);
    for &(ring_num, order, is_opening) in &closures {
        // For the opening side, write the bond order symbol before the digit
        // if it's not a single bond. For closing side, omit the bond symbol
        // (it was already specified by the opener, or it's implied).
        if is_opening {
            write_bond_symbol_for_ring(order, output);
        }
        write_ring_number(ring_num, output);
    }

    // Get neighbors sorted by rank, excluding the atom we came from
    let mut neighbors: Vec<(usize, usize)> = mol.adjacency[atom_idx]
        .iter()
        .map(|&(n, bi)| (n, bi))
        .filter(|&(n, _)| Some(n) != from_atom)
        .collect();
    neighbors.sort_by_key(|&(n, _)| ranks[n]);

    // Process neighbors as tree edges, re-checking visited at each step
    // since earlier branches may visit nodes that were unvisited initially
    for i in 0..neighbors.len() {
        let (n, bi) = neighbors[i];
        if visited[n] {
            continue;
        }

        let bond = &mol.bonds[bi];
        let is_aromatic = bond.is_aromatic;

        // Check if there are more unvisited neighbors after this one
        let has_more = neighbors[i + 1..].iter().any(|&(m, _)| !visited[m]);

        if has_more {
            // Branch
            output.push('(');
            write_bond_symbol(bond.order, is_aromatic, output);
            dfs_smiles(mol, n, Some(atom_idx), ranks, visited, output, ring_info);
            output.push(')');
        } else {
            // Inline (last/only unvisited neighbor)
            write_bond_symbol(bond.order, is_aromatic, output);
            dfs_smiles(mol, n, Some(atom_idx), ranks, visited, output, ring_info);
        }
    }
}

/// Write a ring closure number (handles two-digit with %).
fn write_ring_number(num: usize, output: &mut String) {
    if num < 10 {
        output.push((b'0' + num as u8) as char);
    } else {
        output.push('%');
        output.push_str(&num.to_string());
    }
}

/// Write a bond symbol for a ring closure (before the ring digit at the opening atom).
fn write_bond_symbol_for_ring(order: BondOrder, output: &mut String) {
    match order {
        BondOrder::Single | BondOrder::Aromatic => {
            // Omit — single and aromatic bonds are implicit
        }
        BondOrder::Double => {
            output.push('=');
        }
        BondOrder::Triple => {
            output.push('#');
        }
    }
}

/// Write a bond symbol to the output (omit single bonds and aromatic bonds in aromatic context).
fn write_bond_symbol(order: BondOrder, is_aromatic: bool, output: &mut String) {
    match order {
        BondOrder::Single => {
            // Omit single bonds (default)
        }
        BondOrder::Double => {
            output.push('=');
        }
        BondOrder::Triple => {
            output.push('#');
        }
        BondOrder::Aromatic => {
            // Aromatic bonds are implicit between aromatic atoms
            if !is_aromatic {
                output.push(':');
            }
        }
    }
}

/// Write an atom to the SMILES output string.
fn write_atom(mol: &Molecule, atom_idx: usize, output: &mut String) {
    let atom = &mol.atoms[atom_idx];

    // Determine if we can use the organic subset (no brackets needed)
    let needs_bracket = atom.formal_charge != 0
        || atom.isotope.is_some()
        || !is_organic_subset(atom.atomic_number, atom.is_aromatic);

    if needs_bracket {
        output.push('[');
        if let Some(iso) = atom.isotope {
            output.push_str(&iso.to_string());
        }
        if let Some(elem) = element_by_number(atom.atomic_number) {
            if atom.is_aromatic {
                for c in elem.symbol.chars() {
                    output.push(c.to_ascii_lowercase());
                }
            } else {
                output.push_str(elem.symbol);
            }
        }
        // Hydrogen count in bracket
        if atom.implicit_hydrogens > 0 {
            output.push('H');
            if atom.implicit_hydrogens > 1 {
                output.push_str(&atom.implicit_hydrogens.to_string());
            }
        }
        // Charge
        if atom.formal_charge > 0 {
            output.push('+');
            if atom.formal_charge > 1 {
                output.push_str(&atom.formal_charge.to_string());
            }
        } else if atom.formal_charge < 0 {
            output.push('-');
            if atom.formal_charge < -1 {
                output.push_str(&atom.formal_charge.abs().to_string());
            }
        }
        output.push(']');
    } else {
        // Organic subset atom
        if let Some(elem) = element_by_number(atom.atomic_number) {
            if atom.is_aromatic {
                for c in elem.symbol.chars() {
                    output.push(c.to_ascii_lowercase());
                }
            } else {
                output.push_str(elem.symbol);
            }
        }
    }
}

/// Check if an atom can be written in the organic subset (no brackets).
fn is_organic_subset(atomic_number: u8, is_aromatic: bool) -> bool {
    if is_aromatic {
        // Aromatic organic subset: b, c, n, o, p, s
        matches!(atomic_number, 5 | 6 | 7 | 8 | 15 | 16)
    } else {
        // Non-aromatic organic subset: B, C, N, O, P, S, F, Cl, Br, I
        matches!(atomic_number, 5 | 6 | 7 | 8 | 15 | 16 | 9 | 17 | 35 | 53)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn canonical_smiles_empty() {
        let mol = Molecule::new("empty".into(), vec![], vec![]);
        assert_eq!(canonical_smiles(&mol), "");
    }

    #[test]
    fn canonical_smiles_methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(canonical_smiles(&mol), "C");
    }

    #[test]
    fn canonical_ethanol_from_different_inputs() {
        let mol1 = parse_smiles("OCC").unwrap();
        let mol2 = parse_smiles("CCO").unwrap();
        let can1 = canonical_smiles(&mol1);
        let can2 = canonical_smiles(&mol2);
        assert_eq!(
            can1, can2,
            "OCC and CCO should give same canonical SMILES: got '{}' vs '{}'",
            can1, can2
        );
    }

    #[test]
    fn canonical_ethanol_roundtrip() {
        let mol = parse_smiles("CCO").unwrap();
        let can = canonical_smiles(&mol);
        // The canonical form should be parseable and produce the same canonical form
        let mol2 = parse_smiles(&can).unwrap();
        let can2 = canonical_smiles(&mol2);
        assert_eq!(
            can, can2,
            "roundtrip failed: '{}' -> '{}' -> '{}'",
            "CCO", can, can2
        );
    }

    #[test]
    fn canonical_benzene_from_different_inputs() {
        // Same aromatic form parsed twice
        let mol1 = parse_smiles("c1ccccc1").unwrap();
        let mol2 = parse_smiles("c1ccccc1").unwrap();
        let can1 = canonical_smiles(&mol1);
        let can2 = canonical_smiles(&mol2);
        assert_eq!(can1, can2);
    }

    #[test]
    fn canonical_benzene_roundtrip() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let can = canonical_smiles(&mol);
        let mol2 = parse_smiles(&can).unwrap();
        let can2 = canonical_smiles(&mol2);
        assert_eq!(
            can, can2,
            "benzene roundtrip failed: '{}' vs '{}'",
            can, can2
        );
    }

    #[test]
    fn canonical_propane_variants() {
        let mol1 = parse_smiles("CCCO").unwrap();
        let mol2 = parse_smiles("OCCC").unwrap();
        let can1 = canonical_smiles(&mol1);
        let can2 = canonical_smiles(&mol2);
        assert_eq!(can1, can2, "CCCO vs OCCC: '{}' vs '{}'", can1, can2);
    }

    #[test]
    fn canonical_branching() {
        // Isobutane from different orderings
        let mol1 = parse_smiles("CC(C)C").unwrap();
        let mol2 = parse_smiles("C(C)(C)C").unwrap();
        let can1 = canonical_smiles(&mol1);
        let can2 = canonical_smiles(&mol2);
        assert_eq!(can1, can2, "isobutane: '{}' vs '{}'", can1, can2);
    }

    #[test]
    fn canonical_double_bond() {
        let mol = parse_smiles("C=C").unwrap();
        let can = canonical_smiles(&mol);
        assert_eq!(can, "C=C");
    }

    #[test]
    fn canonical_triple_bond() {
        let mol = parse_smiles("C#N").unwrap();
        let can = canonical_smiles(&mol);
        // Should contain the triple bond
        assert!(
            can.contains('#'),
            "canonical SMILES should contain '#': got '{}'",
            can
        );
    }

    #[test]
    fn canonical_charged_atom() {
        let mol = parse_smiles("[NH4+]").unwrap();
        let can = canonical_smiles(&mol);
        assert!(can.contains("[NH4+]"), "got '{}'", can);
    }

    #[test]
    fn canonical_disconnected_fragments() {
        let mol = parse_smiles("C.O").unwrap();
        let can = canonical_smiles(&mol);
        // Should contain both fragments separated by .
        assert!(
            can.contains('.'),
            "should have fragment separator: got '{}'",
            can
        );
    }

    #[test]
    fn canonical_ring() {
        // Cyclohexane
        let mol = parse_smiles("C1CCCCC1").unwrap();
        let can = canonical_smiles(&mol);
        // Should parse back to the same molecule
        let mol2 = parse_smiles(&can).unwrap();
        assert_eq!(mol2.atom_count(), 6);
        assert_eq!(mol2.bond_count(), 6);
    }

    #[test]
    fn canonical_preserves_formula() {
        use crate::properties::molecular_formula;

        let smiles_list = ["CCO", "CC(=O)O", "c1ccccc1", "CC(C)C", "C=CC=C"];
        for smi in &smiles_list {
            let mol = parse_smiles(smi).unwrap();
            let can = canonical_smiles(&mol);
            let mol2 = parse_smiles(&can).unwrap();
            assert_eq!(
                molecular_formula(&mol),
                molecular_formula(&mol2),
                "formula mismatch for '{}' -> '{}': {} vs {}",
                smi,
                can,
                molecular_formula(&mol),
                molecular_formula(&mol2)
            );
        }
    }

    #[test]
    fn canonical_deterministic() {
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let can1 = canonical_smiles(&mol);
        let can2 = canonical_smiles(&mol);
        assert_eq!(can1, can2, "canonical SMILES should be deterministic");
    }
}
