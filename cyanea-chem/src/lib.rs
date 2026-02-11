//! Chemistry and small-molecule analysis for the Cyanea bioinformatics ecosystem.
//!
//! Provides molecular graph representation, SMILES/SDF parsing, property calculation,
//! Morgan fingerprints with Tanimoto similarity, substructure search, stereochemistry
//! assignment (R/S, E/Z), and canonical SMILES generation.
//!
//! # Example
//!
//! ```
//! use cyanea_chem::{parse_smiles, compute_properties, morgan_fingerprint, tanimoto_similarity, canonical_smiles};
//!
//! // Parse ethanol from SMILES
//! let ethanol = parse_smiles("CCO").unwrap();
//! assert_eq!(ethanol.atom_count(), 3);
//!
//! // Compute properties
//! let props = compute_properties(&ethanol);
//! assert_eq!(props.formula, "C2H6O");
//!
//! // Generate fingerprint and compare
//! let fp1 = morgan_fingerprint(&ethanol, 2, 2048);
//! let fp2 = morgan_fingerprint(&ethanol, 2, 2048);
//! assert!((tanimoto_similarity(&fp1, &fp2) - 1.0).abs() < 1e-10);
//!
//! // Canonical SMILES: same molecule from different inputs gives same output
//! let mol1 = parse_smiles("OCC").unwrap();
//! let mol2 = parse_smiles("CCO").unwrap();
//! assert_eq!(canonical_smiles(&mol1), canonical_smiles(&mol2));
//! ```

pub mod canon;
pub mod element;
pub mod fingerprint;
pub mod molecule;
pub mod properties;
pub mod sdf;
pub mod smiles;
pub mod stereo;
pub mod substructure;

mod ring;

pub use canon::canonical_smiles;
pub use element::{element_by_number, element_by_symbol, Element};
pub use fingerprint::{morgan_fingerprint, tanimoto_bulk, tanimoto_similarity, Fingerprint};
pub use molecule::{Bond, BondOrder, BondStereo, Chirality, MolAtom, Molecule};
pub use properties::{compute_properties, molecular_formula, molecular_weight, MolecularProperties};
pub use sdf::{parse_mol_v2000, parse_sdf};
pub use smiles::{parse_smiles, parse_smiles_named};
pub use stereo::{assign_ez, assign_rs};
pub use substructure::{find_substructure_matches, has_substructure, SubstructureMatch};

#[cfg(feature = "std")]
pub use sdf::parse_sdf_file;
