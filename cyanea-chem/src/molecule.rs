//! Molecular graph representation.

use cyanea_core::{Annotated, ContentAddressable, Summarizable};
use sha2::{Digest, Sha256};

/// Tetrahedral chirality at a stereocenter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Chirality {
    /// No chirality specified.
    #[default]
    None,
    /// Counterclockwise (`@` in SMILES).
    CounterClockwise,
    /// Clockwise (`@@` in SMILES).
    Clockwise,
}

/// Cis-trans stereo bond direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondStereo {
    /// No stereo bond.
    #[default]
    None,
    /// Up bond (`/` in SMILES).
    Up,
    /// Down bond (`\` in SMILES).
    Down,
}

/// Bond order classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    /// Numeric bond order for valence calculations.
    pub fn as_f64(self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
        }
    }
}

/// An atom in a molecular graph.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MolAtom {
    pub atomic_number: u8,
    pub formal_charge: i8,
    pub isotope: Option<u16>,
    pub is_aromatic: bool,
    pub implicit_hydrogens: u8,
    pub chirality: Chirality,
}

/// A bond between two atoms.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
    pub is_aromatic: bool,
    pub stereo: BondStereo,
}

/// A molecular graph with atoms, bonds, and adjacency information.
#[derive(Debug, Clone)]
pub struct Molecule {
    pub name: String,
    pub atoms: Vec<MolAtom>,
    pub bonds: Vec<Bond>,
    /// adjacency[atom_idx] = Vec<(neighbor_atom_idx, bond_idx)>
    pub adjacency: Vec<Vec<(usize, usize)>>,
}

impl Molecule {
    /// Create a new molecule, building the adjacency list from atoms and bonds.
    pub fn new(name: String, atoms: Vec<MolAtom>, bonds: Vec<Bond>) -> Self {
        let mut adjacency = vec![Vec::new(); atoms.len()];
        for (bi, bond) in bonds.iter().enumerate() {
            adjacency[bond.atom1].push((bond.atom2, bi));
            adjacency[bond.atom2].push((bond.atom1, bi));
        }
        Molecule { name, atoms, bonds, adjacency }
    }

    /// Number of atoms (including implicit hydrogens conceptually, but counting graph nodes).
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Number of bonds.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Number of non-hydrogen atoms.
    pub fn heavy_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| a.atomic_number != 1).count()
    }

    /// Neighbor atom indices for a given atom.
    pub fn neighbors(&self, atom_idx: usize) -> Vec<usize> {
        self.adjacency[atom_idx].iter().map(|&(n, _)| n).collect()
    }

    /// Graph degree of an atom (number of explicit bonds).
    pub fn degree(&self, atom_idx: usize) -> usize {
        self.adjacency[atom_idx].len()
    }

    /// Find the bond between two atoms, if any.
    pub fn get_bond(&self, a1: usize, a2: usize) -> Option<&Bond> {
        self.adjacency[a1]
            .iter()
            .find(|&&(n, _)| n == a2)
            .map(|&(_, bi)| &self.bonds[bi])
    }

    /// Total hydrogen count (implicit + explicit H atoms).
    pub fn total_hydrogen_count(&self) -> usize {
        let explicit: usize = self.atoms.iter().filter(|a| a.atomic_number == 1).count();
        let implicit: usize = self.atoms.iter().map(|a| a.implicit_hydrogens as usize).sum();
        explicit + implicit
    }
}

impl Annotated for Molecule {
    fn name(&self) -> &str {
        &self.name
    }
}

impl Summarizable for Molecule {
    fn summary(&self) -> String {
        format!(
            "{}: {} atoms, {} bonds",
            if self.name.is_empty() { "Molecule" } else { &self.name },
            self.atom_count(),
            self.bond_count()
        )
    }
}

impl ContentAddressable for Molecule {
    fn content_hash(&self) -> String {
        let mut hasher = Sha256::new();
        // Sort atoms by (atomic_number, charge, isotope, aromatic, implicit_h)
        let mut sorted_atoms: Vec<_> = self.atoms.iter().enumerate().collect();
        sorted_atoms.sort_by_key(|(_, a)| {
            (a.atomic_number, a.formal_charge, a.isotope, a.is_aromatic, a.implicit_hydrogens)
        });
        for (_, atom) in &sorted_atoms {
            hasher.update([atom.atomic_number]);
            hasher.update(atom.formal_charge.to_le_bytes());
            hasher.update(atom.implicit_hydrogens.to_le_bytes());
            if let Some(iso) = atom.isotope {
                hasher.update(iso.to_le_bytes());
            }
            hasher.update([atom.is_aromatic as u8]);
        }
        // Sort bonds by (min_atom, max_atom, order)
        let mut sorted_bonds: Vec<_> = self.bonds.iter().collect();
        sorted_bonds.sort_by_key(|b| {
            let (a, c) = if b.atom1 <= b.atom2 { (b.atom1, b.atom2) } else { (b.atom2, b.atom1) };
            (a, c, b.order as u8)
        });
        for bond in &sorted_bonds {
            hasher.update(bond.atom1.to_le_bytes());
            hasher.update(bond.atom2.to_le_bytes());
            hasher.update([bond.order as u8]);
        }
        hex::encode(hasher.finalize())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_water() -> Molecule {
        let atoms = vec![
            MolAtom { atomic_number: 8, formal_charge: 0, isotope: None, is_aromatic: false, implicit_hydrogens: 2, chirality: Chirality::None },
        ];
        Molecule::new("water".into(), atoms, vec![])
    }

    fn make_ethane() -> Molecule {
        let atoms = vec![
            MolAtom { atomic_number: 6, formal_charge: 0, isotope: None, is_aromatic: false, implicit_hydrogens: 3, chirality: Chirality::None },
            MolAtom { atomic_number: 6, formal_charge: 0, isotope: None, is_aromatic: false, implicit_hydrogens: 3, chirality: Chirality::None },
        ];
        let bonds = vec![
            Bond { atom1: 0, atom2: 1, order: BondOrder::Single, is_aromatic: false, stereo: BondStereo::None },
        ];
        Molecule::new("ethane".into(), atoms, bonds)
    }

    #[test]
    fn construction_and_adjacency() {
        let mol = make_ethane();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.adjacency[0].len(), 1);
        assert_eq!(mol.adjacency[1].len(), 1);
    }

    #[test]
    fn neighbors_and_degree() {
        let mol = make_ethane();
        assert_eq!(mol.neighbors(0), vec![1]);
        assert_eq!(mol.degree(0), 1);
        assert_eq!(mol.degree(1), 1);
    }

    #[test]
    fn heavy_atom_count() {
        let mol = make_water();
        assert_eq!(mol.heavy_atom_count(), 1);
        assert_eq!(mol.total_hydrogen_count(), 2);
    }

    #[test]
    fn summarizable_and_content_addressable() {
        let mol = make_ethane();
        assert!(mol.summary().contains("2 atoms"));
        let hash = mol.content_hash();
        assert_eq!(hash.len(), 64);
        // Deterministic
        assert_eq!(hash, mol.content_hash());
    }
}
