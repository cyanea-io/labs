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
}
