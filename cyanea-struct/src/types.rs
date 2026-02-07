//! Core types for macromolecular 3D structure representation.

use cyanea_core::{Annotated, ContentAddressable, Summarizable};
use sha2::{Digest, Sha256};

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

/// A point in 3D Cartesian space.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    /// Create a new point.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The origin.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Euclidean distance to another point.
    pub fn distance_to(&self, other: &Point3D) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Dot product.
    pub fn dot(&self, other: &Point3D) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    pub fn cross(&self, other: &Point3D) -> Point3D {
        Point3D {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Vector magnitude.
    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Unit vector in the same direction, or zero if magnitude is zero.
    pub fn normalize(&self) -> Point3D {
        let n = self.norm();
        if n < 1e-15 {
            Point3D::zero()
        } else {
            Point3D {
                x: self.x / n,
                y: self.y / n,
                z: self.z / n,
            }
        }
    }

    /// Vector addition.
    pub fn add(&self, other: &Point3D) -> Point3D {
        Point3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    /// Vector subtraction.
    pub fn sub(&self, other: &Point3D) -> Point3D {
        Point3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    /// Scalar multiplication.
    pub fn scale(&self, s: f64) -> Point3D {
        Point3D {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

/// A single atom in a macromolecular structure.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Atom {
    /// Atom serial number.
    pub serial: u32,
    /// Atom name (e.g. "CA", "N", "CB").
    pub name: String,
    /// Alternate location indicator.
    pub alt_loc: Option<char>,
    /// 3D coordinates in Angstroms.
    pub coords: Point3D,
    /// Occupancy factor.
    pub occupancy: f64,
    /// Temperature factor (B-factor).
    pub temp_factor: f64,
    /// Element symbol.
    pub element: Option<String>,
    /// Formal charge.
    pub charge: Option<i8>,
    /// Whether this is a HETATM record.
    pub is_hetatm: bool,
}

impl Atom {
    /// Whether this atom is a backbone atom (N, CA, C, O).
    pub fn is_backbone(&self) -> bool {
        let trimmed = self.name.trim();
        matches!(trimmed, "N" | "CA" | "C" | "O")
    }

    /// Whether this is an alpha carbon.
    pub fn is_alpha_carbon(&self) -> bool {
        self.name.trim() == "CA"
    }
}

/// A residue (amino acid or nucleotide) in a chain.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Residue {
    /// Three-letter residue name (e.g. "ALA", "GLY").
    pub name: String,
    /// Sequence number from the PDB file.
    pub seq_num: i32,
    /// Insertion code.
    pub i_code: Option<char>,
    /// Atoms belonging to this residue.
    pub atoms: Vec<Atom>,
}

impl Residue {
    /// Get an atom by name, returning the first match.
    pub fn get_atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.name.trim() == name)
    }

    /// Get the alpha carbon atom.
    pub fn get_alpha_carbon(&self) -> Option<&Atom> {
        self.get_atom("CA")
    }

    /// Return all backbone atoms (N, CA, C, O).
    pub fn backbone_atoms(&self) -> Vec<&Atom> {
        self.atoms.iter().filter(|a| a.is_backbone()).collect()
    }

    /// Geometric center of mass (unweighted) of all atoms.
    pub fn center_of_mass(&self) -> Point3D {
        if self.atoms.is_empty() {
            return Point3D::zero();
        }
        let mut sum = Point3D::zero();
        for atom in &self.atoms {
            sum = sum.add(&atom.coords);
        }
        sum.scale(1.0 / self.atoms.len() as f64)
    }
}

impl Annotated for Residue {
    fn name(&self) -> &str {
        &self.name
    }
}

/// A polypeptide or polynucleotide chain.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Chain {
    /// Single-character chain identifier.
    pub id: char,
    /// Residues in this chain, in sequence order.
    pub residues: Vec<Residue>,
    /// String form of chain ID for trait impl.
    chain_id_str: String,
}

impl Chain {
    /// Create a new chain.
    pub fn new(id: char, residues: Vec<Residue>) -> Self {
        Self {
            id,
            residues,
            chain_id_str: format!("Chain {}", id),
        }
    }

    /// Number of residues.
    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    /// Total number of atoms across all residues.
    pub fn atom_count(&self) -> usize {
        self.residues.iter().map(|r| r.atoms.len()).sum()
    }

    /// Collect all alpha carbon atoms in this chain.
    pub fn alpha_carbons(&self) -> Vec<&Atom> {
        self.residues
            .iter()
            .filter_map(|r| r.get_alpha_carbon())
            .collect()
    }
}

impl Annotated for Chain {
    fn name(&self) -> &str {
        &self.chain_id_str
    }
}

/// A complete macromolecular structure (one or more chains).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Structure {
    /// PDB identifier or user-supplied name.
    pub id: String,
    /// Chains in this structure.
    pub chains: Vec<Chain>,
}

impl Structure {
    /// Number of chains.
    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    /// Total residues across all chains.
    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|c| c.residue_count()).sum()
    }

    /// Total atoms across all chains.
    pub fn atom_count(&self) -> usize {
        self.chains.iter().map(|c| c.atom_count()).sum()
    }

    /// Get a chain by its single-character ID.
    pub fn get_chain(&self, id: char) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }

    /// Collect all atoms across all chains.
    pub fn all_atoms(&self) -> Vec<&Atom> {
        self.chains
            .iter()
            .flat_map(|c| c.residues.iter().flat_map(|r| r.atoms.iter()))
            .collect()
    }

    /// Collect all alpha carbon atoms across all chains.
    pub fn alpha_carbons(&self) -> Vec<&Atom> {
        self.chains.iter().flat_map(|c| c.alpha_carbons()).collect()
    }

    /// Geometric center of mass of all atoms.
    pub fn center_of_mass(&self) -> Point3D {
        let atoms = self.all_atoms();
        if atoms.is_empty() {
            return Point3D::zero();
        }
        let mut sum = Point3D::zero();
        for atom in &atoms {
            sum = sum.add(&atom.coords);
        }
        sum.scale(1.0 / atoms.len() as f64)
    }
}

impl Annotated for Structure {
    fn name(&self) -> &str {
        &self.id
    }
}

impl Summarizable for Structure {
    fn summary(&self) -> String {
        format!(
            "Structure {} — {} chain(s), {} residue(s), {} atom(s)",
            self.id,
            self.chain_count(),
            self.residue_count(),
            self.atom_count(),
        )
    }
}

impl ContentAddressable for Structure {
    fn content_hash(&self) -> String {
        let mut hasher = Sha256::new();
        hasher.update(self.id.as_bytes());
        for chain in &self.chains {
            hasher.update(&[chain.id as u8]);
            for residue in &chain.residues {
                hasher.update(residue.name.as_bytes());
                hasher.update(&residue.seq_num.to_le_bytes());
                for atom in &residue.atoms {
                    hasher.update(atom.name.as_bytes());
                    hasher.update(&atom.coords.x.to_le_bytes());
                    hasher.update(&atom.coords.y.to_le_bytes());
                    hasher.update(&atom.coords.z.to_le_bytes());
                }
            }
        }
        hex::encode(hasher.finalize())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

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

    #[test]
    fn point3d_arithmetic() {
        let a = Point3D::new(1.0, 2.0, 3.0);
        let b = Point3D::new(4.0, 5.0, 6.0);
        assert_eq!(a.add(&b), Point3D::new(5.0, 7.0, 9.0));
        assert_eq!(a.sub(&b), Point3D::new(-3.0, -3.0, -3.0));
        assert!((a.dot(&b) - 32.0).abs() < 1e-10);
        assert!((a.scale(2.0).x - 2.0).abs() < 1e-10);
        assert!((a.distance_to(&b) - (27.0_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn point3d_cross_product() {
        let x = Point3D::new(1.0, 0.0, 0.0);
        let y = Point3D::new(0.0, 1.0, 0.0);
        let z = x.cross(&y);
        assert!((z.x).abs() < 1e-10);
        assert!((z.y).abs() < 1e-10);
        assert!((z.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn atom_backbone_detection() {
        let ca = make_atom("CA", 0.0, 0.0, 0.0);
        let cb = make_atom("CB", 0.0, 0.0, 0.0);
        let n = make_atom("N", 0.0, 0.0, 0.0);
        assert!(ca.is_backbone());
        assert!(ca.is_alpha_carbon());
        assert!(!cb.is_backbone());
        assert!(!cb.is_alpha_carbon());
        assert!(n.is_backbone());
    }

    #[test]
    fn residue_get_alpha_carbon() {
        let r = Residue {
            name: "ALA".into(),
            seq_num: 1,
            i_code: None,
            atoms: vec![
                make_atom("N", 0.0, 0.0, 0.0),
                make_atom("CA", 1.0, 0.0, 0.0),
                make_atom("C", 2.0, 0.0, 0.0),
            ],
        };
        assert!(r.get_alpha_carbon().is_some());
        assert_eq!(r.backbone_atoms().len(), 3);
    }

    #[test]
    fn structure_summary_and_hash() {
        let chain = Chain::new(
            'A',
            vec![Residue {
                name: "GLY".into(),
                seq_num: 1,
                i_code: None,
                atoms: vec![make_atom("CA", 1.0, 2.0, 3.0)],
            }],
        );
        let s = Structure {
            id: "1ABC".into(),
            chains: vec![chain],
        };
        assert!(s.summary().contains("1ABC"));
        assert!(s.summary().contains("1 chain"));
        assert!(s.summary().contains("1 residue"));
        assert!(s.summary().contains("1 atom"));

        let hash = s.content_hash();
        assert_eq!(hash.len(), 64); // SHA-256 hex
        // Deterministic
        assert_eq!(hash, s.content_hash());
    }
}
