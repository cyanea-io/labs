//! Protein and nucleic acid 3D structures for the Cyanea bioinformatics ecosystem.
//!
//! - **PDB parsing** — Read macromolecular structure files with [`pdb::parse_pdb`]
//! - **Coordinate geometry** — Distance, angle, dihedral, RMSD in [`geometry`]
//! - **Secondary structure** — Simplified DSSP assignment in [`secondary`]
//! - **Superposition** — Kabsch structural alignment in [`superposition`]
//! - **Contact maps** — Residue-residue contact analysis in [`contact`]
//!
//! # Quick start
//!
//! ```
//! use cyanea_struct::pdb::parse_pdb;
//! use cyanea_struct::types::Structure;
//! use cyanea_core::Summarizable;
//!
//! let pdb_text = "\
//! HEADER                                                        1TST
//! ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
//! ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
//! ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
//! ATOM      4  O   ALA A   1       3.000   3.000   3.000  1.00  0.00           O
//! TER
//! END
//! ";
//!
//! let structure = parse_pdb(pdb_text).unwrap();
//! assert_eq!(structure.chain_count(), 1);
//! assert!(structure.summary().contains("1TST"));
//! ```

#![no_std]

extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

pub mod contact;
pub mod geometry;
mod linalg;
pub mod pdb;
pub mod secondary;
pub mod superposition;
pub mod types;

pub use contact::{compute_contact_map, compute_contact_map_allatom, ContactMap};
pub use geometry::{angle, angle_points, center_of_mass, dihedral, dihedral_points, distance};
pub use pdb::parse_pdb;
pub use secondary::{
    assign_secondary_structure, backbone_dihedrals, SecondaryStructure,
    SecondaryStructureAssignment,
};
pub use superposition::{kabsch, kabsch_points, SuperpositionResult};
pub use types::{Atom, Chain, Point3D, Residue, Structure};

#[cfg(test)]
mod tests {
    use super::*;
    use cyanea_core::Summarizable;

    #[test]
    fn integration_parse_and_analyze() {
        let pdb_text = "\
HEADER                                                        1INT\n\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n\
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n\
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n\
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O\n\
ATOM      5  N   GLY A   2       3.325   1.506   0.000  1.00  0.00           N\n\
ATOM      6  CA  GLY A   2       3.988   2.802   0.000  1.00  0.00           C\n\
ATOM      7  C   GLY A   2       5.504   2.714   0.000  1.00  0.00           C\n\
ATOM      8  O   GLY A   2       6.092   1.635   0.000  1.00  0.00           O\n\
ATOM      9  N   VAL A   3       6.120   3.898   0.000  1.00  0.00           N\n\
ATOM     10  CA  VAL A   3       7.574   3.984   0.000  1.00  0.00           C\n\
ATOM     11  C   VAL A   3       8.173   2.578   0.000  1.00  0.00           C\n\
ATOM     12  O   VAL A   3       9.398   2.445   0.000  1.00  0.00           O\n\
TER\n\
END\n";

        let s = parse_pdb(pdb_text).unwrap();
        assert_eq!(s.id, "1INT");
        assert_eq!(s.chain_count(), 1);
        assert_eq!(s.residue_count(), 3);
        assert!(s.summary().contains("3 residue"));

        // Contact map
        let chain = s.get_chain('A').unwrap();
        let cm = compute_contact_map(chain).unwrap();
        assert_eq!(cm.size, 3);
        assert!(cm.get(0, 0).abs() < 1e-10); // diagonal

        // SS assignment
        let ss = assign_secondary_structure(chain).unwrap();
        assert_eq!(ss.assignments.len(), 3);
    }
}
