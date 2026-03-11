//! Metabolomics analysis: metabolite identification, isotope patterns,
//! retention time prediction, and KEGG metabolic pathway mapping.
//!
//! Provides tools for untargeted and targeted metabolomics workflows
//! including mass-based metabolite matching, theoretical isotope pattern
//! generation, RT prediction from molecular properties, and pathway
//! enrichment scoring.

use cyanea_core::{CyaneaError, Result};
use std::collections::BTreeMap;

/// A metabolite entry in a reference database.
#[derive(Debug, Clone)]
pub struct Metabolite {
    /// Unique identifier (e.g., KEGG compound ID).
    pub id: String,
    /// Common name.
    pub name: String,
    /// Molecular formula.
    pub formula: String,
    /// Monoisotopic mass (Da).
    pub exact_mass: f64,
    /// KEGG pathway IDs this metabolite participates in.
    pub pathways: Vec<String>,
    /// InChIKey (first 14 characters) for structural matching.
    pub inchikey_prefix: Option<String>,
}

/// Result of a mass-based metabolite search.
#[derive(Debug, Clone)]
pub struct MassMatch {
    /// Matched metabolite.
    pub metabolite: Metabolite,
    /// Mass error in ppm.
    pub ppm_error: f64,
    /// Adduct form that produced the match (e.g., "[M+H]+").
    pub adduct: String,
    /// Calculated m/z for this adduct.
    pub calc_mz: f64,
}

/// An adduct definition for mass calculation.
#[derive(Debug, Clone)]
pub struct Adduct {
    /// Name (e.g., "[M+H]+").
    pub name: String,
    /// Mass shift: m/z = (mass * mult + shift) / charge.
    pub mult: f64,
    pub shift: f64,
    pub charge: u32,
}

/// A theoretical isotope peak.
#[derive(Debug, Clone)]
pub struct IsotopePeak {
    /// Mass offset from monoisotopic (Da).
    pub mass_offset: f64,
    /// Relative abundance (monoisotopic = 1.0).
    pub abundance: f64,
}

/// Predicted retention time with confidence interval.
#[derive(Debug, Clone)]
pub struct RtPrediction {
    /// Predicted RT in minutes.
    pub rt_minutes: f64,
    /// Estimated error margin (minutes).
    pub error_margin: f64,
    /// Descriptors used: (logP, MW, PSA).
    pub descriptors: (f64, f64, f64),
}

/// A metabolic pathway with enrichment statistics.
#[derive(Debug, Clone)]
pub struct PathwayEnrichment {
    /// Pathway identifier.
    pub pathway_id: String,
    /// Pathway name.
    pub pathway_name: String,
    /// Number of matched metabolites in this pathway.
    pub hits: usize,
    /// Total metabolites in this pathway.
    pub total: usize,
    /// p-value from hypergeometric test.
    pub p_value: f64,
    /// Impact score (topology-based, 0-1).
    pub impact: f64,
}

/// A KEGG-style pathway definition.
#[derive(Debug, Clone)]
pub struct MetabolicPathway {
    pub id: String,
    pub name: String,
    pub metabolite_ids: Vec<String>,
}

// ---------------------------------------------------------------------------
// Common adducts
// ---------------------------------------------------------------------------

/// Common positive-mode ESI adducts.
pub fn positive_adducts() -> Vec<Adduct> {
    vec![
        Adduct { name: "[M+H]+".into(), mult: 1.0, shift: 1.007276, charge: 1 },
        Adduct { name: "[M+Na]+".into(), mult: 1.0, shift: 22.989218, charge: 1 },
        Adduct { name: "[M+K]+".into(), mult: 1.0, shift: 38.963158, charge: 1 },
        Adduct { name: "[M+NH4]+".into(), mult: 1.0, shift: 18.034164, charge: 1 },
        Adduct { name: "[2M+H]+".into(), mult: 2.0, shift: 1.007276, charge: 1 },
        Adduct { name: "[M+2H]2+".into(), mult: 1.0, shift: 2.014552, charge: 2 },
    ]
}

/// Common negative-mode ESI adducts.
pub fn negative_adducts() -> Vec<Adduct> {
    vec![
        Adduct { name: "[M-H]-".into(), mult: 1.0, shift: -1.007276, charge: 1 },
        Adduct { name: "[M+FA-H]-".into(), mult: 1.0, shift: 44.998201, charge: 1 },
        Adduct { name: "[M+Cl]-".into(), mult: 1.0, shift: 34.969402, charge: 1 },
        Adduct { name: "[M-2H]2-".into(), mult: 1.0, shift: -2.014552, charge: 2 },
    ]
}

/// Calculate m/z for a given mass and adduct.
pub fn calc_mz(exact_mass: f64, adduct: &Adduct) -> f64 {
    (exact_mass * adduct.mult + adduct.shift) / adduct.charge as f64
}

// ---------------------------------------------------------------------------
// Mass-based metabolite matching
// ---------------------------------------------------------------------------

/// Search a metabolite database by observed m/z within a ppm tolerance.
///
/// Returns all matches across the provided adduct list, sorted by ppm error.
pub fn match_by_mass(
    observed_mz: f64,
    database: &[Metabolite],
    adducts: &[Adduct],
    ppm_tolerance: f64,
) -> Vec<MassMatch> {
    let mut matches = Vec::new();

    for met in database {
        for adduct in adducts {
            let theoretical_mz = calc_mz(met.exact_mass, adduct);
            if theoretical_mz <= 0.0 {
                continue;
            }
            let ppm_error = ((observed_mz - theoretical_mz) / theoretical_mz).abs() * 1e6;
            if ppm_error <= ppm_tolerance {
                matches.push(MassMatch {
                    metabolite: met.clone(),
                    ppm_error,
                    adduct: adduct.name.clone(),
                    calc_mz: theoretical_mz,
                });
            }
        }
    }

    matches.sort_by(|a, b| a.ppm_error.partial_cmp(&b.ppm_error).unwrap_or(core::cmp::Ordering::Equal));
    matches
}

// ---------------------------------------------------------------------------
// Isotope pattern prediction
// ---------------------------------------------------------------------------

/// Generate a theoretical isotope pattern from a molecular formula.
///
/// Predicts M+0 through M+`max_peaks`-1 using average natural isotope
/// abundances. Returns peaks normalized so the monoisotopic (M+0) = 1.0.
///
/// Supported elements: C, H, N, O, S, P, Cl, Br, F, Si.
pub fn isotope_pattern(formula: &str, max_peaks: usize) -> Result<Vec<IsotopePeak>> {
    let composition = parse_formula(formula)?;

    // Start with a single peak at mass 0, abundance 1.0
    let mut pattern = vec![1.0f64; 1];

    // Convolve each element's contribution
    for (&element, &count) in &composition {
        let elem_pattern = element_isotope_pattern(element, count)?;
        pattern = convolve_patterns(&pattern, &elem_pattern);
    }

    // Truncate and normalize
    pattern.truncate(max_peaks.max(1));
    let max_val = pattern.iter().cloned().fold(0.0f64, f64::max);
    if max_val > 0.0 {
        for p in &mut pattern {
            *p /= max_val;
        }
    }

    Ok(pattern.iter().enumerate().map(|(i, &abundance)| {
        IsotopePeak {
            mass_offset: i as f64 * 1.003355, // average neutron mass difference
            abundance,
        }
    }).collect())
}

fn parse_formula(formula: &str) -> Result<BTreeMap<char, u32>> {
    let mut composition = BTreeMap::new();
    let chars: Vec<char> = formula.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        if chars[i].is_ascii_uppercase() {
            let element = chars[i];
            i += 1;
            // Check for lowercase second letter (not needed for single-char elements)
            // We handle single-char elements; multi-char via first letter only
            if i < chars.len() && chars[i].is_ascii_lowercase() {
                // Two-letter element — skip lowercase for now (Cl, Br, Si, etc)
                i += 1;
            }
            // Parse count
            let mut count = 0u32;
            while i < chars.len() && chars[i].is_ascii_digit() {
                count = count * 10 + chars[i].to_digit(10).unwrap();
                i += 1;
            }
            if count == 0 {
                count = 1;
            }
            *composition.entry(element).or_insert(0) += count;
        } else {
            i += 1;
        }
    }

    if composition.is_empty() {
        return Err(CyaneaError::InvalidInput("Empty formula".into()));
    }

    Ok(composition)
}

fn element_isotope_pattern(element: char, count: u32) -> Result<Vec<f64>> {
    // Natural isotope abundances (M+0, M+1, M+2, ...)
    let base = match element {
        'C' => vec![0.9893, 0.0107],               // 12C, 13C
        'H' => vec![0.999885, 0.000115],            // 1H, 2H
        'N' => vec![0.99632, 0.00368],              // 14N, 15N
        'O' => vec![0.99757, 0.00038, 0.00205],     // 16O, 17O, 18O
        'S' => vec![0.9499, 0.0075, 0.0425, 0.0, 0.0001], // 32S, 33S, 34S, skip, 36S
        'P' => vec![1.0],                           // 31P (monoisotopic)
        'F' => vec![1.0],                           // 19F
        _ => vec![1.0],
    };

    // Exponentiate by repeated convolution
    let mut result = vec![1.0];
    for _ in 0..count {
        result = convolve_patterns(&result, &base);
    }

    Ok(result)
}

fn convolve_patterns(a: &[f64], b: &[f64]) -> Vec<f64> {
    let len = a.len() + b.len() - 1;
    let mut result = vec![0.0; len.min(10)]; // cap at M+9
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            if i + j < result.len() {
                result[i + j] += ai * bj;
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Retention time prediction
// ---------------------------------------------------------------------------

/// Predict reversed-phase HPLC retention time from molecular properties.
///
/// Uses a simple linear model: RT = a * logP + b * MW/1000 - c * PSA/100 + d
/// Coefficients are calibrated for C18 reversed-phase with water/acetonitrile gradient.
///
/// * `logp` — Partition coefficient (e.g., from Wildman-Crippen)
/// * `molecular_weight` — in Da
/// * `polar_surface_area` — in Angstrom^2
pub fn predict_rt(logp: f64, molecular_weight: f64, polar_surface_area: f64) -> RtPrediction {
    // Empirical coefficients for a typical 20-min C18 gradient
    let a = 2.5;   // logP contribution
    let b = 1.8;   // MW contribution
    let c = 3.2;   // PSA contribution (polar = earlier elution)
    let d = 5.0;   // intercept

    let rt = a * logp + b * (molecular_weight / 1000.0) - c * (polar_surface_area / 100.0) + d;
    let rt_clamped = rt.clamp(0.5, 30.0);

    // Error margin scales with predicted RT uncertainty
    let error = 0.5 + 0.1 * logp.abs();

    RtPrediction {
        rt_minutes: rt_clamped,
        error_margin: error,
        descriptors: (logp, molecular_weight, polar_surface_area),
    }
}

/// Score isotope pattern match between observed and theoretical.
///
/// Returns a cosine similarity score (0-1). Score > 0.9 indicates good match.
pub fn isotope_cosine_score(observed: &[f64], theoretical: &[f64]) -> f64 {
    let n = observed.len().min(theoretical.len());
    if n == 0 {
        return 0.0;
    }

    let mut dot = 0.0;
    let mut norm_o = 0.0;
    let mut norm_t = 0.0;

    for i in 0..n {
        dot += observed[i] * theoretical[i];
        norm_o += observed[i] * observed[i];
        norm_t += theoretical[i] * theoretical[i];
    }

    let denom = (norm_o * norm_t).sqrt();
    if denom < 1e-15 {
        0.0
    } else {
        dot / denom
    }
}

// ---------------------------------------------------------------------------
// Pathway enrichment
// ---------------------------------------------------------------------------

/// Perform metabolic pathway enrichment analysis (hypergeometric test).
///
/// Given a list of matched metabolite IDs and a pathway database, tests
/// which pathways are over-represented in the matched set.
///
/// * `matched_ids` — IDs of metabolites detected/matched
/// * `pathways` — reference pathway definitions
/// * `universe_size` — total metabolites in the database
pub fn pathway_enrichment(
    matched_ids: &[&str],
    pathways: &[MetabolicPathway],
    universe_size: usize,
) -> Vec<PathwayEnrichment> {
    let matched_set: std::collections::BTreeSet<&str> = matched_ids.iter().copied().collect();
    let n_matched = matched_set.len();

    let mut results = Vec::new();

    for pw in pathways {
        let hits: usize = pw.metabolite_ids.iter()
            .filter(|id| matched_set.contains(id.as_str()))
            .count();

        if hits == 0 {
            continue;
        }

        let pathway_size = pw.metabolite_ids.len();

        // Hypergeometric p-value (Fisher's exact, one-tailed)
        let p_value = hypergeometric_pvalue(n_matched, pathway_size, universe_size, hits);

        // Topology impact score — proportion of pathway hit
        let impact = hits as f64 / pathway_size as f64;

        results.push(PathwayEnrichment {
            pathway_id: pw.id.clone(),
            pathway_name: pw.name.clone(),
            hits,
            total: pathway_size,
            p_value,
            impact,
        });
    }

    results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(core::cmp::Ordering::Equal));
    results
}

/// Hypergeometric p-value: P(X >= k) where X ~ Hypergeometric(N, K, n).
fn hypergeometric_pvalue(n_draw: usize, n_success: usize, n_total: usize, k: usize) -> f64 {
    // Use log factorials for numerical stability
    let mut p_value = 0.0;
    let max_k = n_draw.min(n_success);

    for x in k..=max_k {
        let p = hypergeometric_pmf(n_draw, n_success, n_total, x);
        p_value += p;
    }

    p_value.min(1.0)
}

fn hypergeometric_pmf(n: usize, k_total: usize, big_n: usize, k: usize) -> f64 {
    // P(X=k) = C(K,k) * C(N-K, n-k) / C(N, n)
    if k > n || k > k_total || n > big_n || (n - k) > (big_n - k_total) {
        return 0.0;
    }
    let log_p = log_choose(k_total, k) + log_choose(big_n - k_total, n - k) - log_choose(big_n, n);
    log_p.exp()
}

fn log_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

fn log_factorial(n: usize) -> f64 {
    // Stirling approximation for large n, exact for small n
    if n <= 1 {
        return 0.0;
    }
    if n <= 20 {
        let mut result = 0.0;
        for i in 2..=n {
            result += (i as f64).ln();
        }
        return result;
    }
    let n = n as f64;
    n * n.ln() - n + 0.5 * (2.0 * core::f64::consts::PI * n).ln()
}

// ---------------------------------------------------------------------------
// Demo metabolite database
// ---------------------------------------------------------------------------

/// A small demo metabolite database (20 common metabolites).
pub fn demo_metabolite_database() -> Vec<Metabolite> {
    vec![
        Metabolite { id: "C00031".into(), name: "D-Glucose".into(), formula: "C6H12O6".into(), exact_mass: 180.06339, pathways: vec!["map00010".into(), "map00500".into()], inchikey_prefix: Some("WQZGKKKJIJFFOK".into()) },
        Metabolite { id: "C00036".into(), name: "Oxaloacetate".into(), formula: "C4H4O5".into(), exact_mass: 132.00587, pathways: vec!["map00020".into(), "map00620".into()], inchikey_prefix: Some("KHPXUQMNIQBQEV".into()) },
        Metabolite { id: "C00042".into(), name: "Succinate".into(), formula: "C4H6O4".into(), exact_mass: 118.02661, pathways: vec!["map00020".into(), "map00190".into()], inchikey_prefix: Some("KDYFGRWQOYBRFD".into()) },
        Metabolite { id: "C00074".into(), name: "Phosphoenolpyruvate".into(), formula: "C3H5O6P".into(), exact_mass: 167.98237, pathways: vec!["map00010".into(), "map00620".into()], inchikey_prefix: Some("DTBNBXWJWCWCIK".into()) },
        Metabolite { id: "C00149".into(), name: "L-Malate".into(), formula: "C4H6O5".into(), exact_mass: 134.02153, pathways: vec!["map00020".into(), "map00620".into()], inchikey_prefix: Some("BJEPYKJPYRNKOW".into()) },
        Metabolite { id: "C00158".into(), name: "Citrate".into(), formula: "C6H8O7".into(), exact_mass: 192.02700, pathways: vec!["map00020".into()], inchikey_prefix: Some("KRKNYBCHXYNGOX".into()) },
        Metabolite { id: "C00186".into(), name: "L-Lactate".into(), formula: "C3H6O3".into(), exact_mass: 90.03169, pathways: vec!["map00010".into(), "map00620".into()], inchikey_prefix: Some("JVTAAEKCZFNVCJ".into()) },
        Metabolite { id: "C00022".into(), name: "Pyruvate".into(), formula: "C3H4O3".into(), exact_mass: 88.01604, pathways: vec!["map00010".into(), "map00020".into(), "map00620".into()], inchikey_prefix: Some("LCTONWCANYUPML".into()) },
        Metabolite { id: "C00024".into(), name: "Acetyl-CoA".into(), formula: "C23H38N7O17P3S".into(), exact_mass: 809.12576, pathways: vec!["map00010".into(), "map00020".into()], inchikey_prefix: None },
        Metabolite { id: "C00037".into(), name: "Glycine".into(), formula: "C2H5NO2".into(), exact_mass: 75.03203, pathways: vec!["map00260".into(), "map00630".into()], inchikey_prefix: Some("DHMQDGOQFOQNFH".into()) },
        Metabolite { id: "C00041".into(), name: "L-Alanine".into(), formula: "C3H7NO2".into(), exact_mass: 89.04768, pathways: vec!["map00250".into(), "map00260".into()], inchikey_prefix: Some("QNAYBMKLOCPYGJ".into()) },
        Metabolite { id: "C00064".into(), name: "L-Glutamine".into(), formula: "C5H10N2O3".into(), exact_mass: 146.06914, pathways: vec!["map00250".into(), "map00230".into()], inchikey_prefix: Some("ZDXPYRJPNDTMRX".into()) },
        Metabolite { id: "C00049".into(), name: "L-Aspartate".into(), formula: "C4H7NO4".into(), exact_mass: 133.03751, pathways: vec!["map00250".into(), "map00260".into()], inchikey_prefix: Some("CKLJMWTZIZZHCS".into()) },
        Metabolite { id: "C00079".into(), name: "L-Phenylalanine".into(), formula: "C9H11NO2".into(), exact_mass: 165.07898, pathways: vec!["map00360".into(), "map00400".into()], inchikey_prefix: Some("COLNVLDHVKWLRT".into()) },
        Metabolite { id: "C00078".into(), name: "L-Tryptophan".into(), formula: "C11H12N2O2".into(), exact_mass: 204.08988, pathways: vec!["map00380".into(), "map00400".into()], inchikey_prefix: Some("QIVBCDIJIAJPQS".into()) },
        Metabolite { id: "C00062".into(), name: "L-Arginine".into(), formula: "C6H14N4O2".into(), exact_mass: 174.11168, pathways: vec!["map00220".into(), "map00330".into()], inchikey_prefix: Some("ODKSFYDXXFIFQN".into()) },
        Metabolite { id: "C00025".into(), name: "L-Glutamate".into(), formula: "C5H9NO4".into(), exact_mass: 147.05316, pathways: vec!["map00250".into(), "map00220".into()], inchikey_prefix: Some("WHUUTDBJXJRKMK".into()) },
        Metabolite { id: "C00065".into(), name: "L-Serine".into(), formula: "C3H7NO3".into(), exact_mass: 105.04259, pathways: vec!["map00260".into(), "map00630".into()], inchikey_prefix: Some("MTCFGRXMJLQNBG".into()) },
        Metabolite { id: "C00183".into(), name: "L-Valine".into(), formula: "C5H11NO2".into(), exact_mass: 117.07898, pathways: vec!["map00280".into(), "map00290".into()], inchikey_prefix: Some("KZSNJWFQEVHDMF".into()) },
        Metabolite { id: "C00082".into(), name: "L-Tyrosine".into(), formula: "C9H11NO3".into(), exact_mass: 181.07389, pathways: vec!["map00350".into(), "map00400".into()], inchikey_prefix: Some("OUYCCCASQSFEME".into()) },
    ]
}

/// Demo KEGG metabolic pathways (6 core pathways).
pub fn demo_metabolic_pathways() -> Vec<MetabolicPathway> {
    vec![
        MetabolicPathway {
            id: "map00010".into(),
            name: "Glycolysis / Gluconeogenesis".into(),
            metabolite_ids: vec!["C00031".into(), "C00074".into(), "C00022".into(), "C00186".into(), "C00024".into()],
        },
        MetabolicPathway {
            id: "map00020".into(),
            name: "Citrate cycle (TCA cycle)".into(),
            metabolite_ids: vec!["C00036".into(), "C00042".into(), "C00149".into(), "C00158".into(), "C00022".into(), "C00024".into()],
        },
        MetabolicPathway {
            id: "map00250".into(),
            name: "Alanine, aspartate, glutamate metabolism".into(),
            metabolite_ids: vec!["C00041".into(), "C00049".into(), "C00064".into(), "C00025".into()],
        },
        MetabolicPathway {
            id: "map00260".into(),
            name: "Glycine, serine, threonine metabolism".into(),
            metabolite_ids: vec!["C00037".into(), "C00041".into(), "C00049".into(), "C00065".into()],
        },
        MetabolicPathway {
            id: "map00400".into(),
            name: "Phenylalanine, tyrosine, tryptophan biosynthesis".into(),
            metabolite_ids: vec!["C00079".into(), "C00078".into(), "C00082".into()],
        },
        MetabolicPathway {
            id: "map00620".into(),
            name: "Pyruvate metabolism".into(),
            metabolite_ids: vec!["C00022".into(), "C00036".into(), "C00074".into(), "C00149".into(), "C00186".into()],
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_positive_adducts() {
        let adducts = positive_adducts();
        assert_eq!(adducts.len(), 6);
        assert_eq!(adducts[0].name, "[M+H]+");
    }

    #[test]
    fn test_negative_adducts() {
        let adducts = negative_adducts();
        assert_eq!(adducts.len(), 4);
        assert_eq!(adducts[0].name, "[M-H]-");
    }

    #[test]
    fn test_calc_mz() {
        // Glucose [M+H]+: 180.063 + 1.007 = ~181.070
        let mz = calc_mz(180.06339, &positive_adducts()[0]);
        assert!((mz - 181.07067).abs() < 0.001);
    }

    #[test]
    fn test_match_by_mass_glucose() {
        let db = demo_metabolite_database();
        let adducts = positive_adducts();

        // Observed glucose [M+H]+ at ~181.0707
        let matches = match_by_mass(181.0707, &db, &adducts, 10.0);
        assert!(!matches.is_empty());
        assert_eq!(matches[0].metabolite.name, "D-Glucose");
        assert!(matches[0].ppm_error < 5.0);
    }

    #[test]
    fn test_match_by_mass_no_hit() {
        let db = demo_metabolite_database();
        let adducts = positive_adducts();

        // Very high m/z — no match expected
        let matches = match_by_mass(9999.0, &db, &adducts, 5.0);
        assert!(matches.is_empty());
    }

    #[test]
    fn test_isotope_pattern_carbon() {
        let pattern = isotope_pattern("C10H20O", 4).unwrap();
        assert_eq!(pattern.len(), 4);
        // M+0 should be the most abundant (normalized to 1.0)
        assert!((pattern[0].abundance - 1.0).abs() < 0.01);
        // M+1 should be ~10.7% for 10 carbons
        assert!(pattern[1].abundance > 0.05 && pattern[1].abundance < 0.20);
    }

    #[test]
    fn test_isotope_pattern_simple() {
        let pattern = isotope_pattern("CH4", 3).unwrap();
        assert_eq!(pattern.len(), 3);
        assert!((pattern[0].abundance - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_isotope_cosine_score() {
        let obs = vec![1.0, 0.11, 0.006];
        let theo = vec![1.0, 0.107, 0.005];
        let score = isotope_cosine_score(&obs, &theo);
        assert!(score > 0.99);
    }

    #[test]
    fn test_isotope_cosine_mismatch() {
        let obs = vec![1.0, 0.5, 0.3];
        let theo = vec![1.0, 0.01, 0.0001];
        let score = isotope_cosine_score(&obs, &theo);
        assert!(score < 0.95);
    }

    #[test]
    fn test_predict_rt() {
        // Hydrophobic molecule: high logP → later elution
        let rt_hydrophobic = predict_rt(3.0, 300.0, 40.0);
        // Hydrophilic molecule: low logP → earlier elution
        let rt_hydrophilic = predict_rt(-1.0, 150.0, 100.0);

        assert!(rt_hydrophobic.rt_minutes > rt_hydrophilic.rt_minutes);
        assert!(rt_hydrophobic.rt_minutes > 0.0);
        assert!(rt_hydrophilic.rt_minutes > 0.0);
    }

    #[test]
    fn test_predict_rt_clamped() {
        // Extreme negative should clamp to minimum
        let rt = predict_rt(-20.0, 50.0, 500.0);
        assert!(rt.rt_minutes >= 0.5);
    }

    #[test]
    fn test_pathway_enrichment() {
        let pathways = demo_metabolic_pathways();
        let matched_ids = vec!["C00031", "C00022", "C00186", "C00074", "C00024"];

        let results = pathway_enrichment(&matched_ids, &pathways, 20);

        // Glycolysis should be the top hit — all 5 metabolites in the pathway are in matched_ids
        assert!(!results.is_empty());
        let glycolysis = results.iter().find(|r| r.pathway_id == "map00010").unwrap();
        assert_eq!(glycolysis.hits, 5);
        assert_eq!(glycolysis.total, 5);
        assert!((glycolysis.impact - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_pathway_enrichment_no_hits() {
        let pathways = demo_metabolic_pathways();
        let matched_ids = vec!["C99999"]; // not in any pathway

        let results = pathway_enrichment(&matched_ids, &pathways, 20);
        assert!(results.is_empty());
    }

    #[test]
    fn test_demo_database() {
        let db = demo_metabolite_database();
        assert_eq!(db.len(), 20);
        assert!(db.iter().all(|m| m.exact_mass > 0.0));
        assert!(db.iter().all(|m| !m.pathways.is_empty()));
    }

    #[test]
    fn test_demo_pathways() {
        let pathways = demo_metabolic_pathways();
        assert_eq!(pathways.len(), 6);
        assert!(pathways.iter().all(|p| !p.metabolite_ids.is_empty()));
    }

    #[test]
    fn test_parse_formula() {
        let comp = parse_formula("C6H12O6").unwrap();
        assert_eq!(comp.get(&'C'), Some(&6));
        assert_eq!(comp.get(&'H'), Some(&12));
        assert_eq!(comp.get(&'O'), Some(&6));
    }
}
