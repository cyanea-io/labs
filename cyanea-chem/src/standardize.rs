//! Molecule standardization pipeline.
//!
//! Preprocessing for consistent molecular representation: salt stripping,
//! fragment selection, charge neutralization, and tautomer canonicalization.

use std::collections::VecDeque;

use cyanea_core::Result;

use crate::molecule::{Bond, BondOrder, MolAtom, Molecule};
use crate::properties::molecular_weight;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A step in the standardization pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StandardizeStep {
    LargestFragment,
    StripSalts,
    Neutralize,
    CanonicalTautomer,
}

/// Configuration for the standardization pipeline.
#[derive(Debug, Clone)]
pub struct StandardizeConfig {
    pub steps: Vec<StandardizeStep>,
}

impl Default for StandardizeConfig {
    fn default() -> Self {
        StandardizeConfig {
            steps: vec![
                StandardizeStep::StripSalts,
                StandardizeStep::LargestFragment,
                StandardizeStep::Neutralize,
                StandardizeStep::CanonicalTautomer,
            ],
        }
    }
}

// ---------------------------------------------------------------------------
// Connected components helper
// ---------------------------------------------------------------------------

fn connected_components(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(curr) = queue.pop_front() {
            component.push(curr);
            for &(neighbor, _) in &mol.adjacency[curr] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }

    components
}

/// Build a sub-molecule from a subset of atom indices.
fn extract_fragment(mol: &Molecule, component: &[usize]) -> Molecule {
    let mut index_map = vec![usize::MAX; mol.atom_count()];
    let mut atoms = Vec::new();

    for (new_idx, &old_idx) in component.iter().enumerate() {
        index_map[old_idx] = new_idx;
        atoms.push(mol.atoms[old_idx].clone());
    }

    let mut bonds = Vec::new();
    for bond in &mol.bonds {
        let a1 = index_map[bond.atom1];
        let a2 = index_map[bond.atom2];
        if a1 != usize::MAX && a2 != usize::MAX {
            bonds.push(Bond {
                atom1: a1,
                atom2: a2,
                order: bond.order,
                is_aromatic: bond.is_aromatic,
                stereo: bond.stereo,
            });
        }
    }

    Molecule::new(mol.name.clone(), atoms, bonds)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Keep only the largest connected fragment (by heavy atom count, then MW).
pub fn largest_fragment(mol: &Molecule) -> Molecule {
    let components = connected_components(mol);
    if components.len() <= 1 {
        return mol.clone();
    }

    let best = components
        .iter()
        .max_by(|a, b| {
            let ha = a.iter().filter(|&&i| mol.atoms[i].atomic_number != 1).count();
            let hb = b.iter().filter(|&&i| mol.atoms[i].atomic_number != 1).count();
            ha.cmp(&hb).then_with(|| {
                let fa = extract_fragment(mol, a);
                let fb = extract_fragment(mol, b);
                molecular_weight(&fa)
                    .partial_cmp(&molecular_weight(&fb))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        })
        .unwrap();

    extract_fragment(mol, best)
}

/// Common salt/counterion atomic numbers and small fragment signatures.
const SALT_ATOMS: &[(u8, i8)] = &[
    (11, 1),  // Na+
    (19, 1),  // K+
    (3, 1),   // Li+
    (17, -1), // Cl-
    (35, -1), // Br-
    (53, -1), // I-
    (9, -1),  // F-
    (20, 2),  // Ca2+
    (12, 2),  // Mg2+
];

/// Remove common salt fragments (counterions, small inorganic fragments).
pub fn strip_salts(mol: &Molecule) -> Molecule {
    let components = connected_components(mol);
    if components.len() <= 1 {
        return mol.clone();
    }

    // Filter out salt-like fragments
    let organic: Vec<&Vec<usize>> = components
        .iter()
        .filter(|comp| !is_salt_fragment(mol, comp))
        .collect();

    if organic.is_empty() {
        // If everything looks like a salt, return the largest fragment
        return largest_fragment(mol);
    }

    if organic.len() == 1 {
        return extract_fragment(mol, organic[0]);
    }

    // Multiple organic fragments remain — return the largest
    let best = organic
        .iter()
        .max_by_key(|comp| comp.iter().filter(|&&i| mol.atoms[i].atomic_number != 1).count())
        .unwrap();
    extract_fragment(mol, best)
}

fn is_salt_fragment(mol: &Molecule, component: &[usize]) -> bool {
    // Single atom salts
    if component.len() == 1 {
        let atom = &mol.atoms[component[0]];
        return SALT_ATOMS.iter().any(|&(an, ch)| atom.atomic_number == an && atom.formal_charge == ch);
    }

    // Small fragments (≤ 3 atoms) that are all inorganic
    if component.len() <= 3 {
        let all_inorganic = component.iter().all(|&i| {
            let an = mol.atoms[i].atomic_number;
            // Not C or Si (organic backbone elements)
            an != 6 && an != 14
        });
        if all_inorganic {
            return true;
        }
    }

    // Acetate, TFA, and other small counterion fragments
    if component.len() <= 4 {
        let heavy_count = component.iter().filter(|&&i| mol.atoms[i].atomic_number != 1).count();
        if heavy_count <= 3 {
            let has_negative = component.iter().any(|&i| mol.atoms[i].formal_charge < 0);
            let no_carbon = component.iter().all(|&i| mol.atoms[i].atomic_number != 6);
            if has_negative && no_carbon {
                return true;
            }
        }
    }

    false
}

/// Neutralize common charged groups.
///
/// N+ with H → remove charge and one H; O-/S- → add H; keeps isolated carboxylates.
pub fn neutralize(mol: &Molecule) -> Molecule {
    let mut atoms: Vec<MolAtom> = mol.atoms.clone();
    let bonds: Vec<Bond> = mol.bonds.clone();

    for (i, atom) in atoms.iter_mut().enumerate() {
        // Positively charged nitrogen with H: remove charge and H
        if atom.atomic_number == 7 && atom.formal_charge > 0 && atom.implicit_hydrogens > 0 {
            atom.formal_charge -= 1;
            atom.implicit_hydrogens -= 1;
        }

        // Negatively charged oxygen: protonate (unless it's a carboxylate)
        if atom.atomic_number == 8 && atom.formal_charge == -1 {
            // Check if this is a carboxylate (O- attached to C=O)
            let is_carboxylate = mol.adjacency[i].iter().any(|&(neighbor, _)| {
                mol.atoms[neighbor].atomic_number == 6
                    && mol.adjacency[neighbor].iter().any(|&(n2, bi2)| {
                        n2 != i
                            && mol.atoms[n2].atomic_number == 8
                            && mol.bonds[bi2].order == BondOrder::Double
                    })
            });

            if !is_carboxylate {
                atom.formal_charge = 0;
                atom.implicit_hydrogens += 1;
            }
        }

        // Negatively charged sulfur: protonate
        if atom.atomic_number == 16 && atom.formal_charge == -1 {
            atom.formal_charge = 0;
            atom.implicit_hydrogens += 1;
        }
    }

    Molecule::new(mol.name.clone(), atoms, bonds)
}

/// Canonical tautomer selection.
///
/// Applies common tautomeric transforms and selects the canonical form
/// using a scoring function that prefers fewer charges, more aromatic atoms,
/// and fewer proton rearrangements.
pub fn canonical_tautomer(mol: &Molecule) -> Molecule {
    let mut best = mol.clone();
    let mut best_score = tautomer_score(&best);

    // Apply tautomeric transforms and keep the best
    let transforms = get_tautomer_transforms();
    let mut changed = true;
    let mut iterations = 0;

    while changed && iterations < 10 {
        changed = false;
        iterations += 1;

        for transform in &transforms {
            if let Some(result) = apply_tautomer_transform(&best, transform) {
                let score = tautomer_score(&result);
                if score > best_score {
                    best = result;
                    best_score = score;
                    changed = true;
                }
            }
        }
    }

    best
}

struct TautomerTransform {
    /// (atomic_number, min_h) for donor
    donor: (u8, u8),
    /// (atomic_number) for acceptor
    acceptor: u8,
    /// Bond order between donor/acceptor before transform
    bond_before: BondOrder,
    /// Bond order after transform
    bond_after: BondOrder,
}

fn get_tautomer_transforms() -> Vec<TautomerTransform> {
    vec![
        // Keto-enol: C(=O)-C-H ↔ C(-OH)=C
        TautomerTransform { donor: (6, 1), acceptor: 8, bond_before: BondOrder::Single, bond_after: BondOrder::Double },
        // Amide-imidic acid: C(=O)-N-H ↔ C(-OH)=N
        TautomerTransform { donor: (7, 1), acceptor: 8, bond_before: BondOrder::Single, bond_after: BondOrder::Double },
        // 1,3-proton shift N: N-C=N-H ↔ N=C-N (acceptor N)
        TautomerTransform { donor: (7, 1), acceptor: 7, bond_before: BondOrder::Single, bond_after: BondOrder::Double },
    ]
}

fn apply_tautomer_transform(mol: &Molecule, transform: &TautomerTransform) -> Option<Molecule> {
    // Find donor-acceptor pairs connected through a shared atom
    for (i, atom) in mol.atoms.iter().enumerate() {
        if atom.atomic_number != transform.donor.0 || atom.implicit_hydrogens < transform.donor.1 {
            continue;
        }

        for &(middle, bi) in &mol.adjacency[i] {
            if mol.bonds[bi].order != transform.bond_before {
                continue;
            }

            for &(acceptor, bj) in &mol.adjacency[middle] {
                if acceptor == i { continue; }
                if mol.atoms[acceptor].atomic_number != transform.acceptor { continue; }
                if mol.bonds[bj].order != transform.bond_after { continue; }

                // Can do the reverse transform: move H from donor to get new tautomer
                let mut new_atoms = mol.atoms.clone();
                let mut new_bonds = mol.bonds.clone();

                new_atoms[i].implicit_hydrogens -= 1;
                new_atoms[acceptor].implicit_hydrogens += 1;
                new_bonds[bi].order = transform.bond_after;
                new_bonds[bj].order = transform.bond_before;

                return Some(Molecule::new(mol.name.clone(), new_atoms, new_bonds));
            }
        }
    }
    None
}

fn tautomer_score(mol: &Molecule) -> i64 {
    let mut score: i64 = 0;

    // Prefer fewer charges
    let total_charge: i64 = mol.atoms.iter().map(|a| a.formal_charge.unsigned_abs() as i64).sum();
    score -= total_charge * 10;

    // Prefer more aromatic atoms
    let aromatic_count = mol.atoms.iter().filter(|a| a.is_aromatic).count() as i64;
    score += aromatic_count * 5;

    // Prefer carbonyl over enol form (more double bonds to O)
    let carbonyl_count = mol
        .bonds
        .iter()
        .filter(|b| {
            b.order == BondOrder::Double
                && (mol.atoms[b.atom1].atomic_number == 8
                    || mol.atoms[b.atom2].atomic_number == 8)
        })
        .count() as i64;
    score += carbonyl_count * 3;

    score
}

/// Apply the full standardization pipeline.
pub fn standardize(mol: &Molecule, config: &StandardizeConfig) -> Result<Molecule> {
    let mut result = mol.clone();

    for step in &config.steps {
        result = match step {
            StandardizeStep::StripSalts => strip_salts(&result),
            StandardizeStep::LargestFragment => largest_fragment(&result),
            StandardizeStep::Neutralize => neutralize(&result),
            StandardizeStep::CanonicalTautomer => canonical_tautomer(&result),
        };
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn largest_fragment_keeps_biggest() {
        // CC.CCCC → butane (4 atoms)
        let mol = parse_smiles("CC.CCCC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 4);
    }

    #[test]
    fn single_fragment_unchanged() {
        let mol = parse_smiles("CCCC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), mol.atom_count());
    }

    #[test]
    fn strip_salts_removes_sodium() {
        // [Na+].CC(=O)[O-] → acetic acid fragment
        let mol = parse_smiles("[Na+].CC(=O)[O-]").unwrap();
        let result = strip_salts(&mol);
        // Should have removed Na+
        assert!(
            result.atom_count() < mol.atom_count(),
            "expected fewer atoms after salt strip: got {} vs {}",
            result.atom_count(),
            mol.atom_count()
        );
        // No sodium in result
        assert!(!result.atoms.iter().any(|a| a.atomic_number == 11));
    }

    #[test]
    fn strip_salts_removes_chloride() {
        // [Cl-].CCN → should remove Cl-
        let mol = parse_smiles("[Cl-].CCN").unwrap();
        let result = strip_salts(&mol);
        assert!(!result.atoms.iter().any(|a| a.atomic_number == 17));
    }

    #[test]
    fn neutralize_amine() {
        // [NH3+]CC → neutralize to NH2-CC
        let mol = parse_smiles("[NH3+]CC").unwrap();
        let result = neutralize(&mol);
        // The nitrogen should have reduced charge
        let n_atom = result.atoms.iter().find(|a| a.atomic_number == 7).unwrap();
        assert!(
            n_atom.formal_charge < 1,
            "charge={}",
            n_atom.formal_charge
        );
    }

    #[test]
    fn neutralize_thiolate() {
        // CC[S-] → CCSH
        let mol = parse_smiles("CC[S-]").unwrap();
        let result = neutralize(&mol);
        let s_atom = result.atoms.iter().find(|a| a.atomic_number == 16).unwrap();
        assert_eq!(s_atom.formal_charge, 0);
        assert_eq!(s_atom.implicit_hydrogens, 1);
    }

    #[test]
    fn pipeline_idempotent() {
        let mol = parse_smiles("CCO").unwrap();
        let config = StandardizeConfig::default();
        let r1 = standardize(&mol, &config).unwrap();
        let r2 = standardize(&r1, &config).unwrap();
        assert_eq!(r1.atom_count(), r2.atom_count());
        assert_eq!(r1.bond_count(), r2.bond_count());
    }

    #[test]
    fn standardize_sodium_acetate() {
        let mol = parse_smiles("[Na+].CC(=O)[O-]").unwrap();
        let config = StandardizeConfig::default();
        let result = standardize(&mol, &config).unwrap();
        // No sodium
        assert!(!result.atoms.iter().any(|a| a.atomic_number == 11));
    }

    #[test]
    fn empty_molecule() {
        let mol = Molecule::new(String::new(), Vec::new(), Vec::new());
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn canonical_tautomer_returns_molecule() {
        let mol = parse_smiles("CC(=O)CC").unwrap();
        let result = canonical_tautomer(&mol);
        assert_eq!(result.atom_count(), mol.atom_count());
    }
}
