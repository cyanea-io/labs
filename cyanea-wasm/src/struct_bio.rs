//! Structural biology WASM bindings: PDB parsing, secondary structure, RMSD.

use serde::Serialize;

use cyanea_struct::geometry::rmsd_points;
use cyanea_struct::pdb::parse_pdb;
use cyanea_struct::secondary::assign_secondary_structure;
use cyanea_struct::types::Point3D;

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// Serializable structure info.
#[derive(Debug, Serialize)]
pub struct JsStructureInfo {
    pub id: String,
    pub chain_count: usize,
    pub residue_count: usize,
    pub atom_count: usize,
    pub chains: Vec<JsChainInfo>,
}

/// Serializable chain info.
#[derive(Debug, Serialize)]
pub struct JsChainInfo {
    pub id: String,
    pub residue_count: usize,
    pub atom_count: usize,
}

/// Serializable secondary structure assignment.
#[derive(Debug, Serialize)]
pub struct JsSecondaryStructure {
    pub chain_id: String,
    pub assignments: Vec<JsSSAssignment>,
}

/// Serializable per-residue assignment.
#[derive(Debug, Serialize)]
pub struct JsSSAssignment {
    pub residue_num: i32,
    pub residue_name: String,
    pub structure: String,
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Parse PDB text and return structure info as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pdb_info(pdb_text: &str) -> String {
    let structure = match parse_pdb(pdb_text) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let chains: Vec<JsChainInfo> = structure
        .chains
        .iter()
        .map(|c| JsChainInfo {
            id: c.id.to_string(),
            residue_count: c.residue_count(),
            atom_count: c.atom_count(),
        })
        .collect();
    let js = JsStructureInfo {
        id: structure.id.clone(),
        chain_count: structure.chain_count(),
        residue_count: structure.residue_count(),
        atom_count: structure.atom_count(),
        chains,
    };
    wasm_ok(&js)
}

/// Assign secondary structure from PDB text and return JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pdb_secondary_structure(pdb_text: &str) -> String {
    let structure = match parse_pdb(pdb_text) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    if structure.chains.is_empty() {
        return wasm_err("no chains found in PDB");
    }
    let chain = &structure.chains[0];
    let ss = match assign_secondary_structure(chain) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let assignments: Vec<JsSSAssignment> = chain
        .residues
        .iter()
        .zip(ss.assignments.iter())
        .map(|(res, assignment)| {
            let structure_name = match assignment {
                cyanea_struct::secondary::SecondaryStructure::Helix => "Helix",
                cyanea_struct::secondary::SecondaryStructure::Sheet => "Sheet",
                cyanea_struct::secondary::SecondaryStructure::Turn => "Turn",
                cyanea_struct::secondary::SecondaryStructure::Coil => "Coil",
            };
            JsSSAssignment {
                residue_num: res.seq_num,
                residue_name: res.name.clone(),
                structure: structure_name.to_string(),
            }
        })
        .collect();
    let js = JsSecondaryStructure {
        chain_id: chain.id.to_string(),
        assignments,
    };
    wasm_ok(&js)
}

/// Compute RMSD between two coordinate sets (JSON arrays of [x,y,z]).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn rmsd(coords1_json: &str, coords2_json: &str) -> String {
    let coords1: Vec<[f64; 3]> = match serde_json::from_str(coords1_json) {
        Ok(c) => c,
        Err(e) => return wasm_err(format!("invalid JSON coords1: {e}")),
    };
    let coords2: Vec<[f64; 3]> = match serde_json::from_str(coords2_json) {
        Ok(c) => c,
        Err(e) => return wasm_err(format!("invalid JSON coords2: {e}")),
    };
    let points1: Vec<Point3D> = coords1
        .iter()
        .map(|c| Point3D::new(c[0], c[1], c[2]))
        .collect();
    let points2: Vec<Point3D> = coords2
        .iter()
        .map(|c| Point3D::new(c[0], c[1], c[2]))
        .collect();
    match rmsd_points(&points1, &points2) {
        Ok(val) => wasm_ok(&val),
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_PDB: &str = "\
HEADER                                                        1TST
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.000   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   2.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   2.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   2.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.000   3.000   3.000  1.00  0.00           O
TER
END
";

    #[test]
    fn pdb_info_basic() {
        let json = pdb_info(SAMPLE_PDB);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let info = &v["ok"];
        assert_eq!(info["id"].as_str().unwrap(), "1TST");
        assert_eq!(info["chain_count"], 1);
        assert_eq!(info["residue_count"], 2);
        assert_eq!(info["atom_count"], 8);
        let chains = info["chains"].as_array().unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0]["id"].as_str().unwrap(), "A");
    }

    #[test]
    fn pdb_info_invalid() {
        let json = pdb_info("not a pdb file");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn secondary_structure_basic() {
        let json = pdb_secondary_structure(SAMPLE_PDB);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let ss = &v["ok"];
        assert_eq!(ss["chain_id"].as_str().unwrap(), "A");
        let assignments = ss["assignments"].as_array().unwrap();
        assert_eq!(assignments.len(), 2);
        // Both should be Coil for such a short chain
        for a in assignments {
            let structure = a["structure"].as_str().unwrap();
            assert!(
                structure == "Coil" || structure == "Turn" || structure == "Helix" || structure == "Sheet",
                "unexpected structure: {structure}"
            );
        }
    }

    #[test]
    fn rmsd_identical() {
        let coords = "[[1,0,0],[0,1,0],[0,0,1]]";
        let json = rmsd(coords, coords);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap()).abs() < 1e-10);
    }

    #[test]
    fn rmsd_different() {
        let c1 = "[[0,0,0],[1,0,0],[0,1,0]]";
        let c2 = "[[1,0,0],[2,0,0],[1,1,0]]";
        let json = rmsd(c1, c2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let val = v["ok"].as_f64().unwrap();
        assert!(val > 0.0);
        // Shifted by (1,0,0): RMSD = 1.0
        assert!((val - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rmsd_mismatched_lengths() {
        let c1 = "[[0,0,0]]";
        let c2 = "[[0,0,0],[1,0,0]]";
        let json = rmsd(c1, c2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn rmsd_invalid_json() {
        let json = rmsd("not json", "[[0,0,0]]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }
}
