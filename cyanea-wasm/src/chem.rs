//! Chemistry WASM bindings: SMILES parsing, molecular properties, fingerprints,
//! similarity, and substructure search.

use serde::Serialize;

use cyanea_chem::{
    canonical_smiles, compute_properties, find_substructure_matches, has_substructure,
    morgan_fingerprint, parse_smiles, tanimoto_similarity,
};
use cyanea_chem::sdf;
use cyanea_chem::maccs;

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// Serializable molecular properties.
#[derive(Debug, Serialize)]
pub struct JsMolecularProperties {
    pub formula: String,
    pub molecular_weight: f64,
    pub atom_count: usize,
    pub bond_count: usize,
    pub ring_count: usize,
    pub rotatable_bonds: usize,
    pub hbd: usize,
    pub hba: usize,
}

/// Serializable fingerprint summary.
#[derive(Debug, Serialize)]
pub struct JsFingerprint {
    pub bits: Vec<usize>,
    pub n_bits: usize,
    pub n_on_bits: usize,
}

/// Serializable substructure result.
#[derive(Debug, Serialize)]
pub struct JsSubstructureResult {
    pub has_match: bool,
    pub match_count: usize,
}

/// Serializable SDF molecule summary.
#[derive(Debug, Serialize)]
pub struct JsSdfMolecule {
    pub name: String,
    pub formula: String,
    pub n_atoms: usize,
    pub n_bonds: usize,
    pub molecular_weight: f64,
}

/// Serializable MACCS fingerprint result.
#[derive(Debug, Serialize)]
pub struct JsMaccsFingerprint {
    pub bits: Vec<usize>,
    pub n_bits_set: usize,
    pub n_total_bits: usize,
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Parse SMILES and return molecular properties as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn smiles_properties(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let props = compute_properties(&mol);
    let js = JsMolecularProperties {
        formula: props.formula,
        molecular_weight: props.molecular_weight,
        atom_count: mol.atom_count(),
        bond_count: mol.bond_count(),
        ring_count: props.ring_count,
        rotatable_bonds: props.rotatable_bonds,
        hbd: props.hydrogen_bond_donors,
        hba: props.hydrogen_bond_acceptors,
    };
    wasm_ok(&js)
}

/// Generate canonical SMILES from an input SMILES string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn canonical(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    wasm_ok(&canonical_smiles(&mol))
}

/// Compute a Morgan fingerprint and return set bits as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn smiles_fingerprint(smiles: &str, radius: usize, n_bits: usize) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let fp = morgan_fingerprint(&mol, radius, n_bits);
    let bits: Vec<usize> = (0..fp.nbits()).filter(|&i| fp.get_bit(i)).collect();
    let js = JsFingerprint {
        n_on_bits: bits.len(),
        bits,
        n_bits: fp.nbits(),
    };
    wasm_ok(&js)
}

/// Tanimoto similarity between two SMILES strings.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn tanimoto(smiles1: &str, smiles2: &str) -> String {
    let mol1 = match parse_smiles(smiles1) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let mol2 = match parse_smiles(smiles2) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let fp1 = morgan_fingerprint(&mol1, 2, 2048);
    let fp2 = morgan_fingerprint(&mol2, 2, 2048);
    let sim = tanimoto_similarity(&fp1, &fp2);
    wasm_ok(&sim)
}

/// Check for substructure match between molecule and pattern SMILES.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn smiles_substructure(molecule: &str, pattern: &str) -> String {
    let mol = match parse_smiles(molecule) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let pat = match parse_smiles(pattern) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let matches = find_substructure_matches(&mol, &pat);
    let js = JsSubstructureResult {
        has_match: has_substructure(&mol, &pat),
        match_count: matches.len(),
    };
    wasm_ok(&js)
}

/// Parse an SDF string and return molecule summaries as JSON.
///
/// Each successfully parsed molecule is returned with its name, formula,
/// atom/bond counts, and molecular weight.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_sdf(sdf_text: &str) -> String {
    let results = sdf::parse_sdf(sdf_text);
    let molecules: Vec<JsSdfMolecule> = results
        .into_iter()
        .filter_map(|r| r.ok())
        .map(|mol| {
            let props = compute_properties(&mol);
            JsSdfMolecule {
                name: mol.name.clone(),
                formula: props.formula,
                n_atoms: mol.atom_count(),
                n_bonds: mol.bond_count(),
                molecular_weight: props.molecular_weight,
            }
        })
        .collect();
    wasm_ok(&molecules)
}

/// Compute a MACCS 166-key structural fingerprint from a SMILES string.
///
/// Returns the set bit positions, count of set bits, and total bit count.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn maccs_fingerprint(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let fp = maccs::maccs_fingerprint(&mol);
    let bits: Vec<usize> = (0..fp.nbits()).filter(|&i| fp.get_bit(i)).collect();
    let js = JsMaccsFingerprint {
        n_bits_set: bits.len(),
        bits,
        n_total_bits: fp.nbits(),
    };
    wasm_ok(&js)
}

// ── 3D Coordinate Generation ─────────────────────────────────────────────

/// Serializable conformer result.
#[derive(Debug, Serialize)]
pub struct JsConformer {
    pub coords: Vec<[f64; 3]>,
    pub n_atoms: usize,
}

/// Serializable conformer set result.
#[derive(Debug, Serialize)]
pub struct JsConformerSet {
    pub conformers: Vec<JsConformer>,
    pub energies: Vec<Option<f64>>,
    pub n_conformers: usize,
}

/// Serializable energy result.
#[derive(Debug, Serialize)]
pub struct JsEnergyResult {
    pub bond_stretch: f64,
    pub angle_bend: f64,
    pub torsion: f64,
    pub van_der_waals: f64,
    pub electrostatic: f64,
    pub out_of_plane: f64,
    pub total: f64,
}

/// Serializable minimize result.
#[derive(Debug, Serialize)]
pub struct JsMinimizeResult {
    pub initial_energy: f64,
    pub final_energy: f64,
    pub n_steps: usize,
    pub converged: bool,
    pub coords: Vec<[f64; 3]>,
}

/// Serializable reaction product.
#[derive(Debug, Serialize)]
pub struct JsReactionProduct {
    pub smiles: String,
}

/// Serializable disconnection.
#[derive(Debug, Serialize)]
pub struct JsDisconnection {
    pub transform_name: String,
    pub smirks: String,
    pub precursors: Vec<String>,
}

/// Serializable atom-atom mapping.
#[derive(Debug, Serialize)]
pub struct JsAtomMapping {
    pub mapping: Vec<[usize; 2]>,
    pub unmapped_reactant: Vec<usize>,
    pub unmapped_product: Vec<usize>,
}

/// Embed a single 3D conformer from SMILES.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn embed_conformer(smiles: &str, seed: u64, use_torsion_prefs: bool, max_minimize_steps: usize) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_chem::EmbedConfig {
        random_seed: seed,
        use_torsion_prefs,
        max_minimize_steps,
        ..Default::default()
    };
    match cyanea_chem::embed_molecule(&mol, &config) {
        Ok(conf) => {
            let js = JsConformer {
                n_atoms: conf.len(),
                coords: conf.coords,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Embed multiple 3D conformers from SMILES with RMSD pruning.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn embed_conformers(smiles: &str, max_conformers: usize, rmsd_threshold: f64, seed: u64) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_chem::EmbedConfig {
        max_conformers,
        rmsd_threshold,
        random_seed: seed,
        ..Default::default()
    };
    match cyanea_chem::embed_multiple(&mol, &config) {
        Ok(cset) => {
            let conformers: Vec<JsConformer> = cset.conformers.iter().map(|c| JsConformer {
                n_atoms: c.len(),
                coords: c.coords.clone(),
            }).collect();
            let js = JsConformerSet {
                n_conformers: conformers.len(),
                energies: cset.energies.clone(),
                conformers,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Compute UFF energy for a SMILES string (auto-embeds 3D coordinates).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn uff_energy_js(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = match cyanea_chem::embed_molecule(&mol, &config) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    match cyanea_chem::uff_energy(&mol, &conf) {
        Ok(e) => wasm_ok(&JsEnergyResult {
            bond_stretch: e.bond_stretch,
            angle_bend: e.angle_bend,
            torsion: e.torsion,
            van_der_waals: e.van_der_waals,
            electrostatic: e.electrostatic,
            out_of_plane: e.out_of_plane,
            total: e.total,
        }),
        Err(e) => wasm_err(e),
    }
}

/// Compute MMFF94 energy for a SMILES string (auto-embeds 3D coordinates).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn mmff94_energy_js(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = match cyanea_chem::embed_molecule(&mol, &config) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    match cyanea_chem::mmff94_energy(&mol, &conf) {
        Ok(e) => wasm_ok(&JsEnergyResult {
            bond_stretch: e.bond_stretch,
            angle_bend: e.angle_bend,
            torsion: e.torsion,
            van_der_waals: e.van_der_waals,
            electrostatic: e.electrostatic,
            out_of_plane: e.out_of_plane,
            total: e.total,
        }),
        Err(e) => wasm_err(e),
    }
}

/// Minimize UFF energy for a SMILES string (auto-embeds, then minimizes).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn minimize_uff(smiles: &str, max_steps: usize, gradient_threshold: f64) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let embed_config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = match cyanea_chem::embed_molecule(&mol, &embed_config) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    let min_config = cyanea_chem::MinimizeConfig {
        max_steps,
        gradient_threshold,
        method: cyanea_chem::MinimizeMethod::SteepestDescent,
    };
    match cyanea_chem::uff_minimize(&mol, &conf, &min_config) {
        Ok(r) => wasm_ok(&JsMinimizeResult {
            initial_energy: r.initial_energy,
            final_energy: r.final_energy,
            n_steps: r.n_steps,
            converged: r.converged,
            coords: r.conformer.coords,
        }),
        Err(e) => wasm_err(e),
    }
}

/// Minimize MMFF94 energy for a SMILES string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn minimize_mmff94(smiles: &str, max_steps: usize, gradient_threshold: f64) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let embed_config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = match cyanea_chem::embed_molecule(&mol, &embed_config) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    let min_config = cyanea_chem::MinimizeConfig {
        max_steps,
        gradient_threshold,
        method: cyanea_chem::MinimizeMethod::SteepestDescent,
    };
    match cyanea_chem::mmff94_minimize(&mol, &conf, &min_config) {
        Ok(r) => wasm_ok(&JsMinimizeResult {
            initial_energy: r.initial_energy,
            final_energy: r.final_energy,
            n_steps: r.n_steps,
            converged: r.converged,
            coords: r.conformer.coords,
        }),
        Err(e) => wasm_err(e),
    }
}

// ── Chemical Reactions ──────────────────────────────────────────────────

/// Apply a SMIRKS reaction to a molecule.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn apply_reaction_js(smiles: &str, smirks: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let rxn = match cyanea_chem::parse_reaction(smirks) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_chem::apply_reaction(&mol, &rxn) {
        Ok(products) => {
            let js: Vec<JsReactionProduct> = products.iter().map(|p| JsReactionProduct {
                smiles: p.smiles.clone(),
            }).collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Find retrosynthetic disconnections for a target SMILES.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn retrosynthetic_disconnect(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let disconnections = cyanea_chem::retrosynthetic_disconnections(&mol);
    let js: Vec<JsDisconnection> = disconnections.iter().map(|d| JsDisconnection {
        transform_name: d.transform_name.clone(),
        smirks: d.smirks.clone(),
        precursors: d.precursors.clone(),
    }).collect();
    wasm_ok(&js)
}

/// Compute atom-atom mapping between two SMILES strings.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn atom_atom_mapping(smiles1: &str, smiles2: &str) -> String {
    let mol1 = match parse_smiles(smiles1) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let mol2 = match parse_smiles(smiles2) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let mapping = cyanea_chem::atom_atom_map(&mol1, &mol2);
    let js = JsAtomMapping {
        mapping: mapping.mapping.iter().map(|&(a, b)| [a, b]).collect(),
        unmapped_reactant: mapping.unmapped_reactant,
        unmapped_product: mapping.unmapped_product,
    };
    wasm_ok(&js)
}

/// Compute Gasteiger-Marsili partial charges from SMILES.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn gasteiger_charges_js(smiles: &str) -> String {
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    match cyanea_chem::gasteiger_charges(&mol) {
        Ok(charges) => wasm_ok(&charges),
        Err(e) => wasm_err(e),
    }
}

/// Compute Tanimoto similarity between two SMILES strings using MACCS fingerprints.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn tanimoto_maccs(smiles1: &str, smiles2: &str) -> String {
    let mol1 = match parse_smiles(smiles1) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let mol2 = match parse_smiles(smiles2) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let fp1 = maccs::maccs_fingerprint(&mol1);
    let fp2 = maccs::maccs_fingerprint(&mol2);
    let sim = tanimoto_similarity(&fp1, &fp2);
    wasm_ok(&sim)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smiles_properties_ethanol() {
        let json = smiles_properties("CCO");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let props = &v["ok"];
        assert_eq!(props["formula"].as_str().unwrap(), "C2H6O");
        assert!(props["molecular_weight"].as_f64().unwrap() > 40.0);
        assert_eq!(props["atom_count"], 3);
        assert_eq!(props["bond_count"], 2);
        assert_eq!(props["hbd"], 1); // OH
        assert_eq!(props["hba"], 1); // O
    }

    #[test]
    fn smiles_properties_invalid() {
        let json = smiles_properties("INVALID_SMILES!!!");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn canonical_smiles_deterministic() {
        let json1 = canonical("OCC");
        let json2 = canonical("CCO");
        let v1: serde_json::Value = serde_json::from_str(&json1).unwrap();
        let v2: serde_json::Value = serde_json::from_str(&json2).unwrap();
        assert_eq!(v1["ok"].as_str().unwrap(), v2["ok"].as_str().unwrap());
    }

    #[test]
    fn canonical_invalid() {
        let json = canonical("INVALID!!!");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn fingerprint_nonzero() {
        let json = smiles_fingerprint("c1ccccc1", 2, 2048);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let fp = &v["ok"];
        assert_eq!(fp["n_bits"], 2048);
        assert!(fp["n_on_bits"].as_u64().unwrap() > 0);
        assert!(!fp["bits"].as_array().unwrap().is_empty());
    }

    #[test]
    fn tanimoto_identical() {
        let json = tanimoto("CCO", "CCO");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn tanimoto_different() {
        let json = tanimoto("CCO", "c1ccccc1");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let sim = v["ok"].as_f64().unwrap();
        assert!(sim >= 0.0 && sim < 1.0);
    }

    #[test]
    fn substructure_benzene_in_phenol() {
        let json = smiles_substructure("Oc1ccccc1", "c1ccccc1");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["has_match"], true);
        assert!(result["match_count"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn substructure_no_match() {
        let json = smiles_substructure("CCCC", "c1ccccc1");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["has_match"], false);
        assert_eq!(result["match_count"], 0);
    }

    // ── SDF parsing tests ────────────────────────────────────────────────

    #[test]
    fn parse_sdf_basic() {
        let sdf_text = "\
\x20\x20Molecule
     CDK    \n
  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$
";
        let json = parse_sdf(sdf_text);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let mols = v["ok"].as_array().unwrap();
        assert_eq!(mols.len(), 1);
        assert_eq!(mols[0]["n_atoms"], 2);
        assert_eq!(mols[0]["n_bonds"], 1);
    }

    #[test]
    fn parse_sdf_empty() {
        let json = parse_sdf("");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let mols = v["ok"].as_array().unwrap();
        assert!(mols.is_empty());
    }

    // ── MACCS fingerprint tests ──────────────────────────────────────────

    #[test]
    fn maccs_fingerprint_basic() {
        let json = maccs_fingerprint("CCO");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let fp = &v["ok"];
        assert!(fp["n_bits_set"].as_u64().unwrap() > 0);
        assert_eq!(fp["n_total_bits"], 166);
    }

    #[test]
    fn maccs_fingerprint_benzene() {
        let json = maccs_fingerprint("c1ccccc1");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let fp = &v["ok"];
        assert!(fp["n_bits_set"].as_u64().unwrap() > 0);
        assert!(!fp["bits"].as_array().unwrap().is_empty());
    }

    // ── Tanimoto MACCS tests ─────────────────────────────────────────────

    #[test]
    fn tanimoto_maccs_identical() {
        let json = tanimoto_maccs("CCO", "CCO");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 1.0).abs() < 1e-10);
    }

    // ── 3D embedding tests ──────────────────────────────────────────────

    #[test]
    fn embed_conformer_basic() {
        let json = embed_conformer("CCO", 42, true, 100);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let conf = &v["ok"];
        assert_eq!(conf["n_atoms"], 3);
        assert!(!conf["coords"].as_array().unwrap().is_empty());
    }

    #[test]
    fn embed_conformers_multiple() {
        let json = embed_conformers("CCCC", 3, 0.1, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let set = &v["ok"];
        assert!(set["n_conformers"].as_u64().unwrap() >= 1);
    }

    // ── Energy tests ────────────────────────────────────────────────────

    #[test]
    fn uff_energy_basic() {
        let json = uff_energy_js("CC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let e = &v["ok"];
        assert!(e["total"].as_f64().unwrap().is_finite());
    }

    #[test]
    fn mmff94_energy_basic() {
        let json = mmff94_energy_js("CC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let e = &v["ok"];
        assert!(e["total"].as_f64().unwrap().is_finite());
    }

    // ── Minimization tests ──────────────────────────────────────────────

    #[test]
    fn minimize_uff_basic() {
        let json = minimize_uff("CC", 50, 1.0);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let r = &v["ok"];
        assert!(r["final_energy"].as_f64().is_some());
        assert!(r["n_steps"].as_u64().is_some());
    }

    #[test]
    fn minimize_mmff94_basic() {
        let json = minimize_mmff94("CC", 50, 1.0);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let r = &v["ok"];
        assert!(r["final_energy"].as_f64().is_some());
    }

    // ── Reaction tests ──────────────────────────────────────────────────

    #[test]
    fn apply_reaction_basic() {
        let json = apply_reaction_js("CCO", "[C:1][OH:2]>>[C:1]=[O:2]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let products = v["ok"].as_array().unwrap();
        assert!(!products.is_empty());
    }

    #[test]
    fn retrosynthetic_basic() {
        let json = retrosynthetic_disconnect("CC(=O)NC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let disconnections = v["ok"].as_array().unwrap();
        assert!(!disconnections.is_empty());
    }

    #[test]
    fn atom_mapping_basic() {
        let json = atom_atom_mapping("CCO", "CC=O");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let m = &v["ok"];
        assert!(!m["mapping"].as_array().unwrap().is_empty());
    }

    // ── Gasteiger charges test ──────────────────────────────────────────

    #[test]
    fn gasteiger_charges_basic() {
        let json = gasteiger_charges_js("CCO");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let charges = v["ok"].as_array().unwrap();
        assert_eq!(charges.len(), 3);
    }
}
