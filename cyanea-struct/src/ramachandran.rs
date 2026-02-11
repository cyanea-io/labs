//! Ramachandran plot validation.
//!
//! Classifies residue backbone dihedral angles (phi, psi) into favored, allowed,
//! or outlier regions of the Ramachandran plot, with special handling for glycine,
//! proline, and pre-proline residues.

use cyanea_core::{CyaneaError, Result};

use crate::secondary::backbone_dihedrals;
use crate::types::{Chain, Structure};

use alloc::string::String;
use alloc::vec::Vec;

/// Classification of a residue's position in the Ramachandran plot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RamachandranRegion {
    /// In a highly populated (favored) region of the plot.
    Favored,
    /// In an additionally allowed region.
    Allowed,
    /// Outside all allowed regions.
    Outlier,
    /// Glycine — uses the expanded symmetric allowed region.
    Glycine,
    /// Proline — uses the restricted proline-specific regions.
    Proline,
    /// Pre-proline — the residue immediately before a proline.
    PreProline,
}

impl RamachandranRegion {
    /// Single-character code for the classification.
    pub fn code(&self) -> char {
        match self {
            RamachandranRegion::Favored => 'F',
            RamachandranRegion::Allowed => 'A',
            RamachandranRegion::Outlier => 'O',
            RamachandranRegion::Glycine => 'G',
            RamachandranRegion::Proline => 'P',
            RamachandranRegion::PreProline => 'p',
        }
    }
}

/// Validate a single residue's phi/psi angles against the Ramachandran plot.
///
/// Uses rectangular approximations of the standard Ramachandran regions:
/// - **General residues**: alpha-helix (~phi=-60, psi=-45), beta-sheet (~phi=-120, psi=130)
/// - **Glycine**: expanded symmetric allowed regions
/// - **Proline**: restricted regions
/// - **Pre-proline**: shifted regions
///
/// `residue_type` should be the three-letter amino acid code (e.g., "ALA", "GLY", "PRO").
/// For pre-proline classification, pass "PRE-PRO" or use the report function which
/// handles this automatically.
pub fn validate_ramachandran(phi: f64, psi: f64, residue_type: &str) -> RamachandranRegion {
    match residue_type {
        "GLY" => validate_glycine(phi, psi),
        "PRO" => validate_proline(phi, psi),
        "PRE-PRO" => validate_pre_proline(phi, psi),
        _ => validate_general(phi, psi),
    }
}

/// Generate a Ramachandran report for all residues in a structure.
///
/// Returns a tuple of (residue_number, residue_name, phi, psi, region) for each
/// residue that has computable backbone dihedrals (requires previous and next
/// residues, so terminal residues are excluded).
///
/// # Errors
///
/// Returns an error if the structure has no chains.
pub fn ramachandran_report(
    structure: &Structure,
) -> Result<Vec<(usize, String, f64, f64, RamachandranRegion)>> {
    if structure.chains.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "cannot generate Ramachandran report for empty structure".into(),
        ));
    }

    let mut report = Vec::new();

    for chain in &structure.chains {
        ramachandran_report_chain(chain, &mut report);
    }

    Ok(report)
}

/// Generate Ramachandran report for a single chain, appending to `report`.
fn ramachandran_report_chain(
    chain: &Chain,
    report: &mut Vec<(usize, String, f64, f64, RamachandranRegion)>,
) {
    let n = chain.residues.len();
    if n < 3 {
        return;
    }

    for i in 1..n - 1 {
        let prev = &chain.residues[i - 1];
        let curr = &chain.residues[i];
        let next = &chain.residues[i + 1];

        if let Some((phi, psi)) = backbone_dihedrals(Some(prev), curr, Some(next)) {
            // Determine if this is a pre-proline residue
            let res_type = if i + 1 < n && chain.residues[i + 1].name == "PRO" {
                "PRE-PRO"
            } else {
                curr.name.as_str()
            };

            let region = validate_ramachandran(phi, psi, res_type);
            report.push((
                curr.seq_num as usize,
                curr.name.clone(),
                phi,
                psi,
                region,
            ));
        }
    }
}

// ---- Region validators using rectangular approximations ----

/// General residues (non-Gly, non-Pro).
fn validate_general(phi: f64, psi: f64) -> RamachandranRegion {
    // Favored regions (rectangular approximations):
    // Alpha-helix region: phi in [-100, -30], psi in [-67, -7]
    // Beta-sheet region: phi in [-170, -70], psi in [90, 180]
    // Left-handed helix (rare): phi in [30, 90], psi in [-10, 60]
    if in_rect(phi, psi, -100.0, -30.0, -67.0, -7.0)
        || in_rect(phi, psi, -170.0, -70.0, 90.0, 180.0)
        || in_rect(phi, psi, 30.0, 90.0, -10.0, 60.0)
    {
        return RamachandranRegion::Favored;
    }

    // Allowed regions (wider envelopes around favored):
    // Extended alpha: phi in [-130, -10], psi in [-80, 10]
    // Extended beta: phi in [-180.0, -40.0], psi in [70.0, 180.0]
    // Extended left-handed: phi in [20, 110], psi in [-30, 80]
    // Bridge region: phi in [-130, -50], psi in [-10, 90]
    if in_rect(phi, psi, -130.0, -10.0, -80.0, 10.0)
        || in_rect(phi, psi, -180.0, -40.0, 70.0, 180.0)
        || in_rect(phi, psi, 20.0, 110.0, -30.0, 80.0)
        || in_rect(phi, psi, -130.0, -50.0, -10.0, 90.0)
    {
        return RamachandranRegion::Allowed;
    }

    RamachandranRegion::Outlier
}

/// Glycine — expanded symmetric regions (no CB clash).
fn validate_glycine(phi: f64, psi: f64) -> RamachandranRegion {
    // Glycine has no side chain, so both positive and negative phi are well-populated.
    // Favored: standard alpha + mirror, standard beta + mirror, and the wide center
    if in_rect(phi, psi, -120.0, -20.0, -70.0, 0.0)
        || in_rect(phi, psi, 20.0, 120.0, 0.0, 70.0)
        || in_rect(phi, psi, -180.0, -50.0, 80.0, 180.0)
        || in_rect(phi, psi, 50.0, 180.0, -180.0, -80.0)
        || in_rect(phi, psi, 40.0, 120.0, -70.0, 0.0)
        || in_rect(phi, psi, -120.0, -40.0, 0.0, 70.0)
    {
        return RamachandranRegion::Glycine;
    }

    // Allowed: almost everything for glycine except extreme outliers
    if in_rect(phi, psi, -180.0, -10.0, -100.0, 30.0)
        || in_rect(phi, psi, 10.0, 180.0, -30.0, 100.0)
        || in_rect(phi, psi, -180.0, -20.0, 50.0, 180.0)
        || in_rect(phi, psi, 20.0, 180.0, -180.0, -50.0)
        || in_rect(phi, psi, -180.0, 180.0, -180.0, 180.0)
    {
        // Glycine is permissive, return Glycine for any recognized region
        return RamachandranRegion::Glycine;
    }

    RamachandranRegion::Glycine
}

/// Proline — restricted regions.
fn validate_proline(phi: f64, psi: f64) -> RamachandranRegion {
    // Proline phi is restricted to ~-75 to -55 due to ring constraint
    // Favored: phi in [-85, -45], psi in [-55, -15] (alpha)
    //          phi in [-85, -45], psi in [110, 170] (beta-like polyproline II)
    if in_rect(phi, psi, -85.0, -45.0, -55.0, -15.0)
        || in_rect(phi, psi, -85.0, -45.0, 110.0, 170.0)
    {
        return RamachandranRegion::Proline;
    }

    // Allowed (wider): phi in [-100, -30], psi in [-70, 0]
    //                  phi in [-100, -30], psi in [90, 180]
    if in_rect(phi, psi, -100.0, -30.0, -70.0, 0.0)
        || in_rect(phi, psi, -100.0, -30.0, 90.0, 180.0)
    {
        return RamachandranRegion::Proline;
    }

    RamachandranRegion::Outlier
}

/// Pre-proline — shifted regions due to steric effects.
fn validate_pre_proline(phi: f64, psi: f64) -> RamachandranRegion {
    // Pre-proline residues have psi shifted toward more negative values
    // Favored: phi in [-100, -40], psi in [-60, -20] (alpha)
    //          phi in [-170, -80], psi in [100, 160] (beta)
    if in_rect(phi, psi, -100.0, -40.0, -60.0, -20.0)
        || in_rect(phi, psi, -170.0, -80.0, 100.0, 160.0)
    {
        return RamachandranRegion::PreProline;
    }

    // Allowed (wider)
    if in_rect(phi, psi, -130.0, -20.0, -80.0, 10.0)
        || in_rect(phi, psi, -180.0, -50.0, 70.0, 180.0)
    {
        return RamachandranRegion::PreProline;
    }

    RamachandranRegion::Outlier
}

/// Check if (phi, psi) is within a rectangular region.
#[inline]
fn in_rect(phi: f64, psi: f64, phi_min: f64, phi_max: f64, psi_min: f64, psi_max: f64) -> bool {
    phi >= phi_min && phi <= phi_max && psi >= psi_min && psi <= psi_max
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Atom, Chain, Point3D, Residue, Structure};
    use alloc::vec;

    #[test]
    fn alpha_helix_favored() {
        let region = validate_ramachandran(-60.0, -45.0, "ALA");
        assert_eq!(region, RamachandranRegion::Favored);
    }

    #[test]
    fn beta_sheet_favored() {
        let region = validate_ramachandran(-120.0, 130.0, "ALA");
        assert_eq!(region, RamachandranRegion::Favored);
    }

    #[test]
    fn left_handed_helix_favored() {
        let region = validate_ramachandran(60.0, 30.0, "ALA");
        assert_eq!(region, RamachandranRegion::Favored);
    }

    #[test]
    fn general_outlier() {
        let region = validate_ramachandran(0.0, 0.0, "ALA");
        assert_eq!(region, RamachandranRegion::Outlier);
    }

    #[test]
    fn glycine_expanded() {
        // Glycine mirror-beta region (positive phi, negative psi) is allowed
        // for glycine but is an outlier for general residues
        let gly_region = validate_ramachandran(150.0, -150.0, "GLY");
        assert_eq!(gly_region, RamachandranRegion::Glycine);

        let ala_region = validate_ramachandran(150.0, -150.0, "ALA");
        assert_eq!(ala_region, RamachandranRegion::Outlier);
    }

    #[test]
    fn proline_restricted() {
        // Standard proline alpha region
        let region = validate_ramachandran(-65.0, -35.0, "PRO");
        assert_eq!(region, RamachandranRegion::Proline);

        // Out of proline range (positive phi)
        let region2 = validate_ramachandran(60.0, 30.0, "PRO");
        assert_eq!(region2, RamachandranRegion::Outlier);
    }

    #[test]
    fn pre_proline_shifted() {
        let region = validate_ramachandran(-70.0, -40.0, "PRE-PRO");
        assert_eq!(region, RamachandranRegion::PreProline);
    }

    #[test]
    fn region_codes() {
        assert_eq!(RamachandranRegion::Favored.code(), 'F');
        assert_eq!(RamachandranRegion::Allowed.code(), 'A');
        assert_eq!(RamachandranRegion::Outlier.code(), 'O');
        assert_eq!(RamachandranRegion::Glycine.code(), 'G');
        assert_eq!(RamachandranRegion::Proline.code(), 'P');
        assert_eq!(RamachandranRegion::PreProline.code(), 'p');
    }

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
    fn ramachandran_report_empty_structure() {
        let s = Structure {
            id: "EMPTY".into(),
            chains: vec![],
        };
        assert!(ramachandran_report(&s).is_err());
    }

    #[test]
    fn ramachandran_report_basic() {
        // Build a 4-residue chain with known backbone geometry.
        // We place atoms to give a well-defined alpha-helix-like phi/psi.
        // Using idealized positions where the dihedral angles are calculable.
        let residues = vec![
            Residue {
                name: "ALA".into(),
                seq_num: 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", 0.0, 0.0, 0.0),
                    make_atom("CA", 1.458, 0.0, 0.0),
                    make_atom("C", 2.009, 1.420, 0.0),
                    make_atom("O", 1.246, 2.390, 0.0),
                ],
            },
            Residue {
                name: "ALA".into(),
                seq_num: 2,
                i_code: None,
                atoms: vec![
                    make_atom("N", 3.325, 1.506, 0.0),
                    make_atom("CA", 3.988, 2.802, 0.0),
                    make_atom("C", 5.504, 2.714, 0.0),
                    make_atom("O", 6.092, 1.635, 0.0),
                ],
            },
            Residue {
                name: "ALA".into(),
                seq_num: 3,
                i_code: None,
                atoms: vec![
                    make_atom("N", 6.120, 3.898, 0.0),
                    make_atom("CA", 7.574, 3.984, 0.0),
                    make_atom("C", 8.173, 2.578, 0.0),
                    make_atom("O", 9.398, 2.445, 0.0),
                ],
            },
            Residue {
                name: "ALA".into(),
                seq_num: 4,
                i_code: None,
                atoms: vec![
                    make_atom("N", 7.340, 1.540, 0.0),
                    make_atom("CA", 7.810, 0.160, 0.0),
                    make_atom("C", 9.330, 0.070, 0.0),
                    make_atom("O", 9.920, -1.010, 0.0),
                ],
            },
        ];

        let chain = Chain::new('A', residues);
        let structure = Structure {
            id: "TEST".into(),
            chains: vec![chain],
        };

        let report = ramachandran_report(&structure).unwrap();
        // Interior residues (2 and 3) should have entries
        assert_eq!(report.len(), 2);
        // Each entry should have the residue name "ALA"
        for (_, name, _, _, _) in &report {
            assert_eq!(name, "ALA");
        }
    }

    #[test]
    fn allowed_region() {
        // phi=-50, psi=5 is in the extended alpha allowed region
        let region = validate_ramachandran(-50.0, 5.0, "ALA");
        assert_eq!(region, RamachandranRegion::Allowed);
    }
}
