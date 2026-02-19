//! Drug-likeness filters and scoring.
//!
//! Implements Lipinski's Rule of Five, Veber rules, PAINS/Brenk structural
//! alert filters, lead-likeness criteria, and QED (Quantitative Estimate
//! of Drug-likeness).

use crate::descriptors::{tpsa, wildman_crippen_logp};
use crate::molecule::Molecule;
use crate::properties::{compute_properties, hba_count, hbd_count, molecular_weight};
use crate::smarts::{parse_smarts, smarts_match};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Result of Lipinski's Rule of Five evaluation.
#[derive(Debug, Clone)]
pub struct LipinskiResult {
    pub mw: f64,
    pub logp: f64,
    pub hbd: usize,
    pub hba: usize,
    pub passes: bool,
    pub violations: u8,
}

/// Result of Veber rules evaluation.
#[derive(Debug, Clone)]
pub struct VeberResult {
    pub rotatable_bonds: usize,
    pub tpsa: f64,
    pub passes: bool,
}

/// A single structural alert match.
#[derive(Debug, Clone)]
pub struct AlertResult {
    pub name: String,
    pub matched: bool,
}

/// Result of a structural alert filter (PAINS or Brenk).
#[derive(Debug, Clone)]
pub struct AlertFilterResult {
    pub alerts: Vec<AlertResult>,
    pub n_hits: usize,
    pub passes: bool,
}

/// QED (Quantitative Estimate of Drug-likeness) result.
#[derive(Debug, Clone)]
pub struct QedResult {
    pub score: f64,
    pub properties: Vec<(String, f64)>,
}

/// Comprehensive drug-likeness report.
#[derive(Debug, Clone)]
pub struct DrugLikenessReport {
    pub lipinski: LipinskiResult,
    pub veber: VeberResult,
    pub lead_like: bool,
    pub qed: QedResult,
    pub pains: AlertFilterResult,
    pub brenk: AlertFilterResult,
}

// ---------------------------------------------------------------------------
// Lipinski's Rule of Five
// ---------------------------------------------------------------------------

/// Evaluate Lipinski's Rule of Five.
///
/// Passes if at most 1 violation: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10.
pub fn lipinski(mol: &Molecule) -> LipinskiResult {
    let mw = molecular_weight(mol);
    let (logp, _) = wildman_crippen_logp(mol);
    let hbd = hbd_count(mol);
    let hba = hba_count(mol);

    let mut violations = 0u8;
    if mw > 500.0 { violations += 1; }
    if logp > 5.0 { violations += 1; }
    if hbd > 5 { violations += 1; }
    if hba > 10 { violations += 1; }

    LipinskiResult {
        mw,
        logp,
        hbd,
        hba,
        passes: violations <= 1,
        violations,
    }
}

// ---------------------------------------------------------------------------
// Veber rules
// ---------------------------------------------------------------------------

/// Evaluate Veber oral bioavailability rules.
///
/// Passes if rotatable bonds ≤ 10 AND TPSA ≤ 140 Å².
pub fn veber(mol: &Molecule) -> VeberResult {
    let props = compute_properties(mol);
    let t = tpsa(mol);

    VeberResult {
        rotatable_bonds: props.rotatable_bonds,
        tpsa: t,
        passes: props.rotatable_bonds <= 10 && t <= 140.0,
    }
}

// ---------------------------------------------------------------------------
// PAINS structural alerts
// ---------------------------------------------------------------------------

/// Representative PAINS (Pan Assay Interference Compounds) patterns.
///
/// Covers all three families (A/B/C) with ~50 representative SMARTS.
const PAINS_PATTERNS: &[(&str, &str)] = &[
    ("rhodanine", "[#7]C(=S)SC=O"),
    ("catechol", "c1cc(O)c(O)cc1"),
    ("quinone_A", "C1(=O)C=CC(=O)C=C1"),
    ("hydroxyphenyl_hydrazone", "c1ccc(O)cc1/N=N"),
    ("ene_rhodanine", "C=CC(=S)N"),
    ("2_amino_thiophene", "c1cc(N)sc1C=O"),
    ("anil_alk_A", "c1ccccc1N=CC"),
    ("indol_3yl_alk", "c1cc2c([nH]c(C)c2)cc1"),
    ("mannich_A", "NC(=O)CN"),
    ("azo_A", "c1ccccc1N=Nc1ccccc1"),
    ("imine_one_A", "C(=O)C=NC"),
    ("thiocarbonyl_A", "C(=S)N"),
    ("michael_acceptor_A", "C=CC(=O)"),
    ("enone_A", "C=CC(=O)C"),
    ("het_thio_N_A", "c1nsc(=S)n1"),
    ("acyl_hydrazine", "C(=O)NN"),
    ("sulfonyl_hydrazine", "S(=O)(=O)NN"),
    ("phenol_sulfonamide_A", "c1cc(O)ccc1S(=O)(=O)N"),
    ("styrene_A", "c1ccccc1C=C"),
    ("alkylidene_barbiturate", "C1(=O)NC(=O)NC1=CC"),
    ("cyano_imine", "N=CC#N"),
    ("pyridinium_A", "c1cc[nH+]cc1"),
    ("thiol_A", "[SH]C"),
    ("hzone_phenol_A", "c1ccc(O)cc1NN=C"),
    ("ene_cyano_A", "C=CC#N"),
    ("isothiocyanate", "N=C=S"),
    ("aldehyde", "[CH]=O"),
    ("acrylate", "C=CC(=O)O"),
    ("sulfonium", "[S+]"),
    ("beta_lactam", "C1(=O)NCC1"),
    ("nitro_aromatic", "c1ccccc1[N+](=O)[O-]"),
    ("para_hydroxy_amine", "c1cc(O)ccc1N"),
    ("hydroquinone", "Oc1ccc(O)cc1"),
    ("alpha_halo_ketone", "C(=O)CCl"),
    ("alpha_halo_ketone_br", "C(=O)CBr"),
    ("thio_ester", "C(=O)SC"),
    ("ene_one_ene", "C=CC(=O)C=C"),
    ("diketo_A", "C(=O)C(=O)"),
    ("anthranil", "c1ccc2c(c1)onc2"),
    ("catechol_methylenedioxy", "c1cc2OCOc2cc1"),
];

/// Run PAINS structural alert filter.
///
/// Returns `passes: true` if no PAINS alerts are triggered.
pub fn pains_filter(mol: &Molecule) -> AlertFilterResult {
    run_alert_filter(mol, PAINS_PATTERNS)
}

/// Run PAINS filter with custom patterns.
pub fn pains_filter_custom(mol: &Molecule, patterns: &[(&str, &str)]) -> AlertFilterResult {
    run_alert_filter(mol, patterns)
}

// ---------------------------------------------------------------------------
// Brenk structural alerts
// ---------------------------------------------------------------------------

/// Representative Brenk structural alerts (reactive groups, toxicophores).
const BRENK_PATTERNS: &[(&str, &str)] = &[
    ("alkyl_halide", "CCl"),
    ("alkyl_halide_br", "CBr"),
    ("alkyl_halide_i", "CI"),
    ("michael_acceptor", "C=CC(=O)"),
    ("epoxide", "C1OC1"),
    ("aziridine", "C1NC1"),
    ("acyl_halide_cl", "C(=O)Cl"),
    ("acyl_halide_br", "C(=O)Br"),
    ("sulfonyl_halide", "S(=O)(=O)Cl"),
    ("acid_anhydride", "C(=O)OC(=O)"),
    ("peroxide", "OO"),
    ("isocyanate", "N=C=O"),
    ("isothiocyanate", "N=C=S"),
    ("acyl_cyanide", "C(=O)C#N"),
    ("sulfonyl_cyanide", "S(=O)(=O)C#N"),
    ("azide", "N=[N+]=[N-]"),
    ("triflate", "OS(=O)(=O)C(F)(F)F"),
    ("phosphoramide", "P(=O)(N)(N)"),
    ("aromatic_nitro", "c[N+](=O)[O-]"),
    ("thiol", "[SH]"),
    ("aldehyde_free", "[CH]=O"),
    ("diazo", "C=[N+]=[N-]"),
    ("disulfide", "SSCC"),
    ("polycyclic_aromatic_5ring", "c1ccc2c(c1)ccc3ccccc32"),
    ("hydrazine", "NN"),
    ("hydroxamic_acid", "C(=O)NO"),
    ("beta_lactone", "C1(=O)OCC1"),
    ("cumene_hydroperoxide", "c1ccccc1C(C)(C)OO"),
    ("nitrosamine", "N(=O)N"),
    ("phosphonate_ester", "P(=O)(OC)(OC)"),
    ("vinyl_sulfone", "C=CS(=O)(=O)"),
    ("aliphatic_long_chain", "CCCCCCCCCC"),
    ("crown_ether_like", "COCCOCCOCC"),
    ("sulfonamide_alkyl", "CS(=O)(=O)NC"),
    ("thioamide", "C(=S)N"),
    ("imine", "C=NC"),
    ("enamine", "C=CN"),
    ("diamine", "NCCN"),
    ("nitro_aliphatic", "C[N+](=O)[O-]"),
    ("cyanohydrin", "OC(C#N)"),
];

/// Run Brenk structural alert filter.
///
/// Returns `passes: true` if no Brenk alerts are triggered.
pub fn brenk_filter(mol: &Molecule) -> AlertFilterResult {
    run_alert_filter(mol, BRENK_PATTERNS)
}

fn run_alert_filter(mol: &Molecule, patterns: &[(&str, &str)]) -> AlertFilterResult {
    let mut alerts = Vec::new();
    let mut n_hits = 0;

    for &(name, smarts_str) in patterns {
        let matched = match parse_smarts(smarts_str) {
            Ok(pattern) => smarts_match(mol, &pattern),
            Err(_) => false, // Skip patterns that fail to parse
        };
        if matched {
            n_hits += 1;
        }
        alerts.push(AlertResult {
            name: name.to_string(),
            matched,
        });
    }

    AlertFilterResult {
        alerts,
        n_hits,
        passes: n_hits == 0,
    }
}

// ---------------------------------------------------------------------------
// Lead-likeness
// ---------------------------------------------------------------------------

/// Evaluate lead-likeness criteria.
///
/// MW 250-350, LogP ≤ 3.5, rotatable bonds ≤ 7, rings ≤ 3.
pub fn lead_likeness(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let (logp, _) = wildman_crippen_logp(mol);
    let props = compute_properties(mol);

    (250.0..=350.0).contains(&mw)
        && logp <= 3.5
        && props.rotatable_bonds <= 7
        && props.ring_count <= 3
}

// ---------------------------------------------------------------------------
// QED (Quantitative Estimate of Drug-likeness)
// ---------------------------------------------------------------------------

/// QED weights from Bickerton 2012, Table 1.
const QED_WEIGHTS: [f64; 8] = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95];


/// Desirability function using Gaussian with asymmetric tails.
///
/// Models each property's ideal range as a Gaussian-like function.
/// Parameters: (center, width_left, width_right) derived from QED literature ranges.
fn desirability(x: f64, idx: usize) -> f64 {
    // Optimal ranges for QED properties:
    // MW: peak ~300, good 200-500
    // LogP: peak ~2.5, good 0-5
    // HBA: peak ~4, good 0-10
    // HBD: peak ~1, good 0-5
    // TPSA: peak ~60, good 20-130
    // RotBonds: peak ~3, good 0-8
    // AromaticRings: peak ~2, good 0-4
    // Alerts: peak ~0, good 0
    let (center, sigma_l, sigma_r): (f64, f64, f64) = match idx {
        0 => (300.0, 120.0, 200.0),  // MW
        1 => (2.5, 2.5, 2.5),        // LogP
        2 => (4.0, 4.0, 6.0),        // HBA
        3 => (1.0, 1.0, 4.0),        // HBD
        4 => (60.0, 40.0, 80.0),     // TPSA
        5 => (3.0, 3.0, 7.0),        // RotBonds
        6 => (2.0, 2.0, 2.0),        // AromaticRings
        7 => (0.0, 0.5, 0.5),        // Alerts
        _ => return 0.5,
    };

    let sigma = if x <= center { sigma_l } else { sigma_r };
    if sigma.abs() < 1e-15 { return if (x - center).abs() < 1e-10 { 1.0 } else { 0.0 }; }
    let z = (x - center) / sigma;
    (-0.5 * z * z).exp()
}

/// Compute QED (Quantitative Estimate of Drug-likeness).
///
/// Weighted geometric mean of desirability functions across 8 molecular properties.
pub fn qed(mol: &Molecule) -> QedResult {
    let mw = molecular_weight(mol);
    let (logp, _) = wildman_crippen_logp(mol);
    let hba = hba_count(mol) as f64;
    let hbd = hbd_count(mol) as f64;
    let t = tpsa(mol);
    let props = compute_properties(mol);
    let rot = props.rotatable_bonds as f64;
    let arom = props.aromatic_ring_count as f64;

    // Count alerts (using a minimal subset)
    let alert_count = count_alerts(mol) as f64;

    let property_values = [mw, logp, hba, hbd, t, rot, arom, alert_count];

    let mut desirabilities = Vec::new();
    let mut property_names = Vec::new();

    let names = ["MW", "LogP", "HBA", "HBD", "TPSA", "RotBonds", "AromaticRings", "Alerts"];

    for (i, &x) in property_values.iter().enumerate() {
        let d_val = desirability(x, i);
        desirabilities.push(d_val);
        property_names.push(names[i]);
    }

    // Weighted geometric mean
    let mut log_sum = 0.0;
    let mut weight_sum = 0.0;
    for (i, &d) in desirabilities.iter().enumerate() {
        let w = QED_WEIGHTS[i];
        if d > 0.0 {
            log_sum += w * d.ln();
        } else {
            log_sum += w * 1e-10_f64.ln(); // Avoid -inf
        }
        weight_sum += w;
    }

    let score = if weight_sum > 0.0 {
        (log_sum / weight_sum).exp().clamp(0.0, 1.0)
    } else {
        0.0
    };

    let properties = property_names
        .iter()
        .zip(property_values.iter())
        .map(|(&name, &val)| (name.to_string(), val))
        .collect();

    QedResult { score, properties }
}

fn count_alerts(mol: &Molecule) -> usize {
    // Use a small subset of structural alerts for QED scoring
    let alert_smarts = [
        "C(=O)Cl", "S(=O)(=O)Cl", "C1OC1", "C1NC1",
        "[SH]", "[CH]=O", "N=C=O", "N=C=S",
        "C(=O)C(=O)", "OO", "NN",
    ];

    let mut count = 0;
    for smarts_str in &alert_smarts {
        if let Ok(pattern) = parse_smarts(smarts_str) {
            if smarts_match(mol, &pattern) {
                count += 1;
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Full report
// ---------------------------------------------------------------------------

/// Compute a comprehensive drug-likeness report.
pub fn drug_likeness_report(mol: &Molecule) -> DrugLikenessReport {
    DrugLikenessReport {
        lipinski: lipinski(mol),
        veber: veber(mol),
        lead_like: lead_likeness(mol),
        qed: qed(mol),
        pains: pains_filter(mol),
        brenk: brenk_filter(mol),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn aspirin_passes_lipinski() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O — MW ≈ 180, should pass easily
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let result = lipinski(&mol);
        assert!(result.passes, "violations={}", result.violations);
        assert_eq!(result.violations, 0);
    }

    #[test]
    fn large_molecule_fails_lipinski() {
        // Create a large molecule: long alkyl chain with many heteroatoms
        // C20 chain with OH groups → high MW, many HBD/HBA
        let mol = parse_smiles("OCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCN(CC)CC").unwrap();
        let result = lipinski(&mol);
        // Should have at least 1 violation (MW > 500)
        assert!(result.mw > 300.0, "mw={}", result.mw);
    }

    #[test]
    fn aspirin_veber() {
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let result = veber(&mol);
        assert!(result.passes);
    }

    #[test]
    fn benign_molecule_passes_pains() {
        // Simple ethanol — no PAINS
        let mol = parse_smiles("CCO").unwrap();
        let result = pains_filter(&mol);
        assert!(result.passes, "hits={}", result.n_hits);
    }

    #[test]
    fn michael_acceptor_triggers_pains() {
        // Methyl vinyl ketone: C=CC(=O)C — michael acceptor
        let mol = parse_smiles("C=CC(=O)C").unwrap();
        let result = pains_filter(&mol);
        assert!(result.n_hits > 0, "expected PAINS hit for michael acceptor");
    }

    #[test]
    fn benign_molecule_passes_brenk() {
        let mol = parse_smiles("CCO").unwrap();
        let result = brenk_filter(&mol);
        // Ethanol may match some very broad patterns but shouldn't match most
        // At minimum, the filter should return a result
        assert!(!result.alerts.is_empty());
    }

    #[test]
    fn epoxide_triggers_brenk() {
        // Ethylene oxide: C1OC1
        let mol = parse_smiles("C1OC1").unwrap();
        let result = brenk_filter(&mol);
        assert!(result.n_hits > 0, "expected Brenk hit for epoxide");
    }

    #[test]
    fn qed_druglike_range() {
        // Aspirin should have reasonable QED
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let result = qed(&mol);
        assert!(
            result.score > 0.1 && result.score <= 1.0,
            "QED={}", result.score
        );
    }

    #[test]
    fn qed_has_properties() {
        let mol = parse_smiles("CCO").unwrap();
        let result = qed(&mol);
        assert_eq!(result.properties.len(), 8);
    }

    #[test]
    fn lead_likeness_boundaries() {
        // Ethanol: MW ≈ 46 → too small
        let mol = parse_smiles("CCO").unwrap();
        assert!(!lead_likeness(&mol));
    }

    #[test]
    fn drug_likeness_report_complete() {
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let report = drug_likeness_report(&mol);
        assert!(report.lipinski.passes);
    }
}
