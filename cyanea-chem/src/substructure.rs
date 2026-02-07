//! Substructure search via VF2 subgraph isomorphism.

use crate::molecule::{BondOrder, Molecule};

/// A mapping from pattern atoms to target atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct SubstructureMatch {
    /// Pairs of (pattern_atom_idx, target_atom_idx).
    pub atom_mapping: Vec<(usize, usize)>,
}

/// Check if `target` contains `pattern` as a substructure.
pub fn has_substructure(target: &Molecule, pattern: &Molecule) -> bool {
    let mut state = Vf2State::new(target, pattern);
    state.search(true);
    !state.matches.is_empty()
}

/// Find all substructure matches of `pattern` in `target`.
pub fn find_substructure_matches(target: &Molecule, pattern: &Molecule) -> Vec<SubstructureMatch> {
    let mut state = Vf2State::new(target, pattern);
    state.search(false);
    state.matches
}

struct Vf2State<'a> {
    target: &'a Molecule,
    pattern: &'a Molecule,
    // core_target[t] = Some(p) means target atom t is mapped to pattern atom p
    core_target: Vec<Option<usize>>,
    // core_pattern[p] = Some(t) means pattern atom p is mapped to target atom t
    core_pattern: Vec<Option<usize>>,
    matches: Vec<SubstructureMatch>,
}

impl<'a> Vf2State<'a> {
    fn new(target: &'a Molecule, pattern: &'a Molecule) -> Self {
        Vf2State {
            target,
            pattern,
            core_target: vec![None; target.atom_count()],
            core_pattern: vec![None; pattern.atom_count()],
            matches: Vec::new(),
        }
    }

    fn search(&mut self, early_exit: bool) {
        // Pre-filter: pattern can't be larger than target
        if self.pattern.atom_count() > self.target.atom_count() {
            return;
        }
        if self.pattern.bond_count() > self.target.bond_count() {
            return;
        }

        // Element count pre-filter
        if !self.element_count_compatible() {
            return;
        }

        self.match_recursive(0, early_exit);
    }

    fn element_count_compatible(&self) -> bool {
        let mut pattern_counts = [0u16; 55];
        let mut target_counts = [0u16; 55];

        for atom in &self.pattern.atoms {
            pattern_counts[atom.atomic_number as usize] += 1;
        }
        for atom in &self.target.atoms {
            target_counts[atom.atomic_number as usize] += 1;
        }

        for i in 0..55 {
            if pattern_counts[i] > target_counts[i] {
                return false;
            }
        }
        true
    }

    fn match_recursive(&mut self, depth: usize, early_exit: bool) {
        if early_exit && !self.matches.is_empty() {
            return;
        }

        if depth == self.pattern.atom_count() {
            // Complete match found
            let mapping: Vec<(usize, usize)> = self
                .core_pattern
                .iter()
                .enumerate()
                .map(|(p, t)| (p, t.unwrap()))
                .collect();
            self.matches.push(SubstructureMatch {
                atom_mapping: mapping,
            });
            return;
        }

        // Next unmapped pattern atom (in order)
        let pattern_atom = depth;

        // Find candidate target atoms
        let candidates = self.find_candidates(pattern_atom);

        for target_atom in candidates {
            if self.core_target[target_atom].is_some() {
                continue;
            }

            if self.is_feasible(pattern_atom, target_atom) {
                // Add to mapping
                self.core_pattern[pattern_atom] = Some(target_atom);
                self.core_target[target_atom] = Some(pattern_atom);

                self.match_recursive(depth + 1, early_exit);

                // Remove from mapping
                self.core_pattern[pattern_atom] = None;
                self.core_target[target_atom] = None;

                if early_exit && !self.matches.is_empty() {
                    return;
                }
            }
        }
    }

    fn find_candidates(&self, pattern_atom: usize) -> Vec<usize> {
        // If the pattern atom has already-mapped neighbors, restrict candidates
        // to neighbors of those mapped target atoms
        let mut candidates: Option<Vec<usize>> = None;

        for &(p_neighbor, _) in &self.pattern.adjacency[pattern_atom] {
            if let Some(t_mapped) = self.core_pattern[p_neighbor] {
                let t_neighbors: Vec<usize> = self
                    .target
                    .adjacency[t_mapped]
                    .iter()
                    .map(|&(n, _)| n)
                    .filter(|&n| self.core_target[n].is_none())
                    .collect();

                candidates = Some(match candidates {
                    None => t_neighbors,
                    Some(existing) => existing
                        .into_iter()
                        .filter(|n| t_neighbors.contains(n))
                        .collect(),
                });
            }
        }

        match candidates {
            Some(c) => c,
            None => (0..self.target.atom_count())
                .filter(|&i| self.core_target[i].is_none())
                .collect(),
        }
    }

    fn is_feasible(&self, pattern_atom: usize, target_atom: usize) -> bool {
        let p_atom = &self.pattern.atoms[pattern_atom];
        let t_atom = &self.target.atoms[target_atom];

        // Atom compatibility: same atomic number
        if p_atom.atomic_number != t_atom.atomic_number {
            return false;
        }

        // Check bonds to already-mapped neighbors
        for &(p_neighbor, p_bond_idx) in &self.pattern.adjacency[pattern_atom] {
            if let Some(t_mapped) = self.core_pattern[p_neighbor] {
                // There must be a bond in target between target_atom and t_mapped
                let t_bond = self.target.get_bond(target_atom, t_mapped);
                match t_bond {
                    None => return false,
                    Some(tb) => {
                        let pb = &self.pattern.bonds[p_bond_idx];
                        if !bond_compatible(pb.order, tb.order) {
                            return false;
                        }
                    }
                }
            }
        }

        true
    }
}

/// Check bond compatibility: pattern bond must match target bond.
/// Aromatic matches aromatic; otherwise exact match.
fn bond_compatible(pattern: BondOrder, target: BondOrder) -> bool {
    match (pattern, target) {
        (BondOrder::Aromatic, BondOrder::Aromatic) => true,
        (BondOrder::Single, BondOrder::Single) => true,
        (BondOrder::Double, BondOrder::Double) => true,
        (BondOrder::Triple, BondOrder::Triple) => true,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn benzene_in_phenol() {
        let phenol = parse_smiles("Oc1ccccc1").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert!(has_substructure(&phenol, &benzene));
        let matches = find_substructure_matches(&phenol, &benzene);
        assert!(!matches.is_empty());
    }

    #[test]
    fn no_match_benzene_in_cyclohexane() {
        let cyclohexane = parse_smiles("C1CCCCC1").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        // Benzene is aromatic, cyclohexane is not — no match
        assert!(!has_substructure(&cyclohexane, &benzene));
    }

    #[test]
    fn multiple_matches() {
        // Naphthalene has 2 benzene substructures
        let naphthalene = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let matches = find_substructure_matches(&naphthalene, &benzene);
        assert!(matches.len() >= 2, "got {} matches", matches.len());
    }
}
