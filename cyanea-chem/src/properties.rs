//! Molecular property calculation.

use std::collections::BTreeMap;

use cyanea_core::Summarizable;

use crate::element::element_by_number;
use crate::molecule::{BondOrder, Molecule};
use crate::ring;

/// Computed molecular properties.
#[derive(Debug, Clone)]
pub struct MolecularProperties {
    pub molecular_weight: f64,
    pub exact_mass: f64,
    pub formula: String,
    pub heavy_atom_count: usize,
    pub hydrogen_bond_donors: usize,
    pub hydrogen_bond_acceptors: usize,
    pub rotatable_bonds: usize,
    pub ring_count: usize,
    pub aromatic_ring_count: usize,
}

impl Summarizable for MolecularProperties {
    fn summary(&self) -> String {
        format!(
            "MW={:.2} Formula={} HBD={} HBA={} RotBonds={} Rings={}",
            self.molecular_weight,
            self.formula,
            self.hydrogen_bond_donors,
            self.hydrogen_bond_acceptors,
            self.rotatable_bonds,
            self.ring_count,
        )
    }
}

/// Compute all molecular properties at once.
pub fn compute_properties(mol: &Molecule) -> MolecularProperties {
    let rings = ring::find_sssr(mol);
    let aromatic_ring_count = rings
        .iter()
        .filter(|r| r.iter().all(|&i| mol.atoms[i].is_aromatic))
        .count();

    MolecularProperties {
        molecular_weight: molecular_weight(mol),
        exact_mass: molecular_weight(mol), // simplified: same as MW for now
        formula: molecular_formula(mol),
        heavy_atom_count: mol.heavy_atom_count(),
        hydrogen_bond_donors: hbd_count(mol),
        hydrogen_bond_acceptors: hba_count(mol),
        rotatable_bonds: rotatable_bond_count(mol, &rings),
        ring_count: rings.len(),
        aromatic_ring_count,
    }
}

/// Calculate the molecular weight (sum of atomic weights including implicit H).
pub fn molecular_weight(mol: &Molecule) -> f64 {
    let h_weight = 1.008;
    let mut mw = 0.0;
    for atom in &mol.atoms {
        if let Some(elem) = element_by_number(atom.atomic_number) {
            mw += elem.atomic_weight;
        }
        mw += atom.implicit_hydrogens as f64 * h_weight;
    }
    mw
}

/// Generate the molecular formula in Hill system order (C first, then H, then alphabetical).
pub fn molecular_formula(mol: &Molecule) -> String {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();

    for atom in &mol.atoms {
        if let Some(elem) = element_by_number(atom.atomic_number) {
            *counts.entry(elem.symbol).or_insert(0) += 1;
        }
        if atom.implicit_hydrogens > 0 {
            *counts.entry("H").or_insert(0) += atom.implicit_hydrogens as usize;
        }
    }

    let mut formula = String::new();

    // Hill system: C first, then H, then alphabetical
    if let Some(&c_count) = counts.get("C") {
        formula.push('C');
        if c_count > 1 {
            formula.push_str(&c_count.to_string());
        }
        counts.remove("C");

        if let Some(&h_count) = counts.get("H") {
            formula.push('H');
            if h_count > 1 {
                formula.push_str(&h_count.to_string());
            }
            counts.remove("H");
        }
    }

    // Remaining elements in alphabetical order (BTreeMap already sorted)
    for (symbol, count) in &counts {
        formula.push_str(symbol);
        if *count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    formula
}

/// Count hydrogen bond donors (N or O atoms with at least one attached H).
pub fn hbd_count(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| {
            (a.atomic_number == 7 || a.atomic_number == 8) && a.implicit_hydrogens > 0
        })
        .count()
}

/// Count hydrogen bond acceptors (N or O atoms).
pub fn hba_count(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| a.atomic_number == 7 || a.atomic_number == 8)
        .count()
}

/// Count rotatable bonds: single bonds not in a ring, not to terminal atoms, not amide C-N.
pub fn rotatable_bond_count(mol: &Molecule, rings: &[Vec<usize>]) -> usize {
    // Build a set of ring bonds
    let ring_bond_set = ring_bond_set(mol, rings);

    mol.bonds
        .iter()
        .enumerate()
        .filter(|&(bi, bond)| {
            // Must be single bond
            if bond.order != BondOrder::Single {
                return false;
            }
            // Not in a ring
            if ring_bond_set.contains(&bi) {
                return false;
            }
            // Not to terminal atoms
            if mol.degree(bond.atom1) <= 1 || mol.degree(bond.atom2) <= 1 {
                return false;
            }
            // Not an amide C-N bond
            if is_amide_bond(mol, bond.atom1, bond.atom2) {
                return false;
            }
            true
        })
        .count()
}

/// Build a set of bond indices that are part of rings.
fn ring_bond_set(mol: &Molecule, rings: &[Vec<usize>]) -> Vec<usize> {
    let mut bond_indices = Vec::new();
    for ring in rings {
        for i in 0..ring.len() {
            let a1 = ring[i];
            let a2 = ring[(i + 1) % ring.len()];
            for &(neighbor, bond_idx) in &mol.adjacency[a1] {
                if neighbor == a2 {
                    bond_indices.push(bond_idx);
                }
            }
        }
    }
    bond_indices.sort();
    bond_indices.dedup();
    bond_indices
}

/// Check if the bond between a1 and a2 is an amide C-N bond (C(=O)-N).
fn is_amide_bond(mol: &Molecule, a1: usize, a2: usize) -> bool {
    let (c_idx, n_idx) = if mol.atoms[a1].atomic_number == 6 && mol.atoms[a2].atomic_number == 7 {
        (a1, a2)
    } else if mol.atoms[a1].atomic_number == 7 && mol.atoms[a2].atomic_number == 6 {
        (a2, a1)
    } else {
        return false;
    };

    // Check if C has a double bond to O
    for &(neighbor, bond_idx) in &mol.adjacency[c_idx] {
        if neighbor == n_idx {
            continue;
        }
        if mol.atoms[neighbor].atomic_number == 8 && mol.bonds[bond_idx].order == BondOrder::Double {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn mw_of_water() {
        // Water: [OH2] — bracket atom with 2 explicit H
        let mol = parse_smiles("[OH2]").unwrap();
        let mw = molecular_weight(&mol);
        // O=15.999, 2H=2.016 → 18.015
        assert!((mw - 18.015).abs() < 0.01, "got {mw}");
    }

    #[test]
    fn formula_of_glucose() {
        // Glucose (open chain): OC(CO)C(O)C(O)C(O)C=O
        // C6H12O6
        let mol = parse_smiles("OC(CO)C(O)C(O)C(O)C=O").unwrap();
        let formula = molecular_formula(&mol);
        assert_eq!(formula, "C6H12O6");
    }

    #[test]
    fn hbd_hba_of_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(hbd_count(&mol), 1); // OH
        assert_eq!(hba_count(&mol), 1); // O
    }

    #[test]
    fn rotatable_bonds_of_butane() {
        let mol = parse_smiles("CCCC").unwrap();
        let rings = ring::find_sssr(&mol);
        assert_eq!(rotatable_bond_count(&mol, &rings), 1);
    }

    #[test]
    fn properties_of_aspirin() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let props = compute_properties(&mol);
        // C9H8O4, MW ≈ 180.16
        assert!((props.molecular_weight - 180.16).abs() < 0.1, "MW={}", props.molecular_weight);
        assert_eq!(props.formula, "C9H8O4");
        assert_eq!(props.ring_count, 1);
        assert_eq!(props.aromatic_ring_count, 1);
        assert!(props.hydrogen_bond_donors >= 1); // COOH
        assert!(props.hydrogen_bond_acceptors >= 2); // multiple O atoms
    }
}
