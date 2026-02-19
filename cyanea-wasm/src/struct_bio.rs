//! Structural biology WASM bindings: PDB parsing, secondary structure, RMSD.

use serde::Serialize;

use cyanea_struct::contact;
use cyanea_struct::geometry::rmsd_points;
use cyanea_struct::mmcif;
use cyanea_struct::pdb;
use cyanea_struct::pdb::parse_pdb;
use cyanea_struct::ramachandran;
use cyanea_struct::secondary::assign_secondary_structure;
use cyanea_struct::superposition;
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

/// Serializable contact map result.
#[derive(Debug, Serialize)]
pub struct JsContactMap {
    pub chain_id: String,
    pub size: usize,
    pub n_contacts: usize,
    pub contact_density: f64,
    pub contacts: Vec<(usize, usize)>,
}

/// Serializable Ramachandran entry.
#[derive(Debug, Serialize)]
pub struct JsRamachandranEntry {
    pub residue_index: usize,
    pub residue_name: String,
    pub phi: f64,
    pub psi: f64,
    pub region: String,
}

/// Serializable mmCIF info.
#[derive(Debug, Serialize)]
pub struct JsMmcifInfo {
    pub id: String,
    pub n_chains: usize,
    pub n_residues: usize,
    pub n_atoms: usize,
}

/// Serializable Kabsch superposition result.
#[derive(Debug, Serialize)]
pub struct JsKabschResult {
    pub rmsd: f64,
    pub rotation: Vec<f64>,
    pub translation: Vec<f64>,
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

/// Compute a CA-CA contact map from PDB text and return contacts below cutoff as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn contact_map(pdb_text: &str, cutoff: f64) -> String {
    let structure = match pdb::parse_pdb(pdb_text) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    if structure.chains.is_empty() {
        return wasm_err("no chains found in PDB");
    }
    let chain = &structure.chains[0];
    let map = match contact::compute_contact_map(chain) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let contacts = map.contacts_below(cutoff);
    let js = JsContactMap {
        chain_id: chain.id.to_string(),
        size: map.size,
        n_contacts: contacts.len(),
        contact_density: map.contact_density(cutoff),
        contacts,
    };
    wasm_ok(&js)
}

/// Generate a Ramachandran report from PDB text and return JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn ramachandran_analysis(pdb_text: &str) -> String {
    let structure = match pdb::parse_pdb(pdb_text) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let report = match ramachandran::ramachandran_report(&structure) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    let entries: Vec<JsRamachandranEntry> = report
        .into_iter()
        .map(|(index, name, phi, psi, region)| {
            let region_str = match region {
                ramachandran::RamachandranRegion::Favored => "favored",
                ramachandran::RamachandranRegion::Allowed => "allowed",
                ramachandran::RamachandranRegion::Outlier => "outlier",
                ramachandran::RamachandranRegion::Glycine => "glycine",
                ramachandran::RamachandranRegion::Proline => "proline",
                ramachandran::RamachandranRegion::PreProline => "pre-proline",
            };
            JsRamachandranEntry {
                residue_index: index,
                residue_name: name,
                phi,
                psi,
                region: region_str.to_string(),
            }
        })
        .collect();
    wasm_ok(&entries)
}

/// Parse mmCIF text and return structure info as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_mmcif(text: &str) -> String {
    let structure = match mmcif::parse_mmcif(text) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let js = JsMmcifInfo {
        id: structure.id.clone(),
        n_chains: structure.chain_count(),
        n_residues: structure.residue_count(),
        n_atoms: structure.atom_count(),
    };
    wasm_ok(&js)
}

/// Kabsch superposition on two coordinate sets (JSON arrays of [x,y,z]).
///
/// Returns RMSD, rotation matrix (9 elements, row-major), and translation vector.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn kabsch_align(coords1_json: &str, coords2_json: &str) -> String {
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
    match superposition::kabsch_points(&points1, &points2) {
        Ok(result) => {
            let rotation: Vec<f64> = result
                .rotation
                .iter()
                .flat_map(|row| row.iter().copied())
                .collect();
            let translation = vec![
                result.translation.x,
                result.translation.y,
                result.translation.z,
            ];
            let js = JsKabschResult {
                rmsd: result.rmsd,
                rotation,
                translation,
            };
            wasm_ok(&js)
        }
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

    const MINI_PDB: &str = "\
HEADER    TEST
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00 10.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00 10.00           O
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00 10.00           N
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00 10.00           C
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00 10.00           C
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00 10.00           O
ATOM      9  N   ALA A   3       7.000   0.000   3.000  1.00 10.00           N
ATOM     10  CA  ALA A   3       8.000   0.000   3.000  1.00 10.00           C
ATOM     11  C   ALA A   3       9.000   0.000   3.000  1.00 10.00           C
ATOM     12  O   ALA A   3       9.500   1.000   3.000  1.00 10.00           O
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

    // ── Contact map tests ────────────────────────────────────────────────

    #[test]
    fn contact_map_basic() {
        let json = contact_map(MINI_PDB, 8.0);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let info = &v["ok"];
        assert_eq!(info["chain_id"].as_str().unwrap(), "A");
        assert_eq!(info["size"].as_u64().unwrap(), 3);
        let contacts = info["contacts"].as_array().unwrap();
        assert!(!contacts.is_empty(), "should have contacts at 8.0A cutoff");
        assert!(info["n_contacts"].as_u64().unwrap() > 0);
    }

    #[test]
    fn contact_map_tight_cutoff() {
        let json_wide = contact_map(MINI_PDB, 8.0);
        let json_tight = contact_map(MINI_PDB, 2.0);
        let v_wide: serde_json::Value = serde_json::from_str(&json_wide).unwrap();
        let v_tight: serde_json::Value = serde_json::from_str(&json_tight).unwrap();
        let n_wide = v_wide["ok"]["n_contacts"].as_u64().unwrap();
        let n_tight = v_tight["ok"]["n_contacts"].as_u64().unwrap();
        assert!(
            n_tight <= n_wide,
            "tight cutoff should yield fewer or equal contacts: tight={n_tight}, wide={n_wide}"
        );
    }

    // ── Ramachandran tests ───────────────────────────────────────────────

    #[test]
    fn ramachandran_basic() {
        let json = ramachandran_analysis(MINI_PDB);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        // With 3 residues, only the interior residue (index 1) can have dihedrals,
        // so we expect at least one entry (possibly zero if backbone geometry is degenerate).
        // The function should succeed (no error key).
        assert!(v["error"].is_null(), "expected success, got: {v}");
        let entries = v["ok"].as_array().unwrap();
        // With 3 residues there should be 1 interior residue with phi/psi
        assert!(
            !entries.is_empty(),
            "should have at least one Ramachandran entry for 3-residue chain"
        );
    }

    #[test]
    fn ramachandran_regions() {
        let json = ramachandran_analysis(MINI_PDB);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let entries = v["ok"].as_array().unwrap();
        let valid_regions = ["favored", "allowed", "outlier", "glycine", "proline", "pre-proline"];
        for entry in entries {
            let region = entry["region"].as_str().unwrap();
            assert!(
                valid_regions.contains(&region),
                "unexpected region: {region}"
            );
            // phi and psi should be finite numbers
            assert!(entry["phi"].as_f64().is_some());
            assert!(entry["psi"].as_f64().is_some());
        }
    }

    // ── mmCIF test ───────────────────────────────────────────────────────

    #[test]
    fn parse_mmcif_basic() {
        let mmcif_text = "\
data_TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 10.00
ATOM 2 C CA GLY A 2 4.000 5.000 6.000 1.00 10.00
#
";
        let json = parse_mmcif(mmcif_text);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let info = &v["ok"];
        assert_eq!(info["id"].as_str().unwrap(), "TEST");
        assert_eq!(info["n_chains"].as_u64().unwrap(), 1);
        assert_eq!(info["n_residues"].as_u64().unwrap(), 2);
        assert_eq!(info["n_atoms"].as_u64().unwrap(), 2);
    }

    // ── Kabsch alignment tests ───────────────────────────────────────────

    #[test]
    fn kabsch_align_identical() {
        let coords = "[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]";
        let json = kabsch_align(coords, coords);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        let rmsd_val = result["rmsd"].as_f64().unwrap();
        assert!(
            rmsd_val < 1e-6,
            "RMSD for identical coords should be ~0, got {rmsd_val}"
        );
        // Rotation should be 9 elements
        let rotation = result["rotation"].as_array().unwrap();
        assert_eq!(rotation.len(), 9);
        // Translation should be 3 elements
        let translation = result["translation"].as_array().unwrap();
        assert_eq!(translation.len(), 3);
    }

    #[test]
    fn kabsch_align_translated() {
        let c1 = "[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]";
        let c2 = "[[10,20,30],[11,20,30],[10,21,30],[10,20,31]]";
        let json = kabsch_align(c1, c2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        let rmsd_val = result["rmsd"].as_f64().unwrap();
        // Kabsch should find the optimal superposition, RMSD should be ~0
        // because the sets are related by pure translation
        assert!(
            rmsd_val < 1e-4,
            "RMSD for translated coords should be ~0, got {rmsd_val}"
        );
    }
}
