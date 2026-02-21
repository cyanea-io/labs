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
pub mod conformer;
pub mod descriptors;
pub mod druglikeness;
pub mod element;
pub mod embed;
pub mod fingerprint;
pub mod forcefield;
pub mod gasteiger;
pub mod maccs;
pub mod molecule;
pub mod properties;
pub mod reaction;
pub mod scaffold;
pub mod sdf;
pub mod smarts;
pub mod smiles;
pub mod standardize;
pub mod stereo;
pub mod substructure;

mod ring;

pub use canon::canonical_smiles;
pub use element::{element_by_number, element_by_symbol, Element};
pub use fingerprint::{morgan_fingerprint, tanimoto_bulk, tanimoto_similarity, Fingerprint};
pub use maccs::maccs_fingerprint;
pub use molecule::{Bond, BondOrder, BondStereo, Chirality, MolAtom, Molecule};
pub use properties::{compute_properties, molecular_formula, molecular_weight, MolecularProperties};
pub use sdf::{parse_mol_v2000, parse_mol_v3000, parse_sdf};
pub use smiles::{parse_smiles, parse_smiles_named};
pub use stereo::{assign_ez, assign_rs};
pub use substructure::{find_substructure_matches, has_substructure, SubstructureMatch};

pub use smarts::{parse_smarts, smarts_find_all, smarts_match, SmartsPattern};
pub use descriptors::{
    autocorrelation_descriptors, balaban_j, bertz_ct, chi_connectivity,
    compute_all_descriptors, estate_indices, fraction_sp3, kappa_shape_indices,
    ring_count_details, tpsa, wiener_index, wildman_crippen_logp, zagreb_indices,
    AutocorrelationResult, DescriptorSet, RingDetails,
};
pub use standardize::{
    canonical_tautomer, largest_fragment, neutralize, standardize, strip_salts,
    StandardizeConfig, StandardizeStep,
};
pub use druglikeness::{
    brenk_filter, drug_likeness_report, lead_likeness, lipinski, pains_filter, qed, veber,
    AlertFilterResult, AlertResult, DrugLikenessReport, LipinskiResult, QedResult, VeberResult,
};
pub use scaffold::{
    generic_scaffold, maximum_common_substructure, murcko_scaffold, r_group_decomposition,
    McsResult, MurckoResult, RGroupResult,
};

pub use conformer::{Conformer, ConformerSet};
pub use embed::{embed_molecule, embed_multiple, EmbedConfig, ForceFieldType};
pub use forcefield::{
    assign_mmff94_types, assign_uff_types, mmff94_energy, mmff94_minimize, uff_energy,
    uff_minimize, EnergyComponents, MinimizeConfig, MinimizeMethod, MinimizeResult,
    Mmff94AtomType, UffAtomType,
};
pub use gasteiger::gasteiger_charges;
pub use reaction::{
    apply_reaction, atom_atom_map, enumerate_reactions, parse_reaction, retrosynthetic_disconnections,
    AtomAtomMapping, Disconnection, Reaction, ReactionProduct,
};

#[cfg(feature = "std")]
pub use sdf::parse_sdf_file;
