//! Protein sequence property analysis.
//!
//! Computes physicochemical and structural properties from amino acid sequences:
//!
//! - **Amino acid composition** — residue counts and fractions
//! - **Hydrophobicity profiles** — Kyte-Doolittle and Hopp-Woods sliding window
//! - **GRAVY** — grand average of hydropathicity
//! - **Isoelectric point** — pI via Henderson-Hasselbalch bisection
//! - **Extinction coefficient** — molar absorptivity at 280 nm
//! - **Secondary structure prediction** — Chou-Fasman nucleation/extension, GOR
//! - **Intrinsic disorder prediction** — windowed propensity with logistic smoothing

use cyanea_core::{CyaneaError, Result};

// ── Amino acid indexing ─────────────────────────────────────────

/// Map amino acid byte to index 0–19. Returns None for non-standard residues.
fn aa_index(aa: u8) -> Option<usize> {
    match aa {
        b'A' => Some(0),
        b'C' => Some(1),
        b'D' => Some(2),
        b'E' => Some(3),
        b'F' => Some(4),
        b'G' => Some(5),
        b'H' => Some(6),
        b'I' => Some(7),
        b'K' => Some(8),
        b'L' => Some(9),
        b'M' => Some(10),
        b'N' => Some(11),
        b'P' => Some(12),
        b'Q' => Some(13),
        b'R' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'V' => Some(17),
        b'W' => Some(18),
        b'Y' => Some(19),
        _ => None,
    }
}

/// Normalize a protein sequence: uppercase and validate standard 20 AAs.
fn normalize_protein(seq: &[u8]) -> Result<Vec<u8>> {
    if seq.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "empty protein sequence".to_string(),
        ));
    }
    seq.iter()
        .map(|&b| {
            let upper = b.to_ascii_uppercase();
            if aa_index(upper).is_some() {
                Ok(upper)
            } else {
                Err(CyaneaError::InvalidInput(format!(
                    "invalid amino acid '{}' in protein sequence",
                    b as char
                )))
            }
        })
        .collect()
}

// ── Hydrophobicity scales ───────────────────────────────────────

/// Kyte-Doolittle (1982) hydropathy values, indexed by aa_index.
const KYTE_DOOLITTLE: [f64; 20] = [
    1.8,  // A
    2.5,  // C
    -3.5, // D
    -3.5, // E
    2.8,  // F
    -0.4, // G
    -3.2, // H
    4.5,  // I
    -3.9, // K
    3.8,  // L
    1.9,  // M
    -3.5, // N
    -1.6, // P
    -3.5, // Q
    -4.5, // R
    -0.8, // S
    -0.7, // T
    4.2,  // V
    -0.9, // W
    -1.3, // Y
];

/// Hopp-Woods (1981) hydrophilicity values, indexed by aa_index.
const HOPP_WOODS: [f64; 20] = [
    -0.5, // A
    -1.0, // C
    3.0,  // D
    3.0,  // E
    -2.5, // F
    0.0,  // G
    -0.5, // H
    -1.8, // I
    3.0,  // K
    -1.8, // L
    -1.3, // M
    0.2,  // N
    0.0,  // P
    0.2,  // Q
    3.0,  // R
    0.3,  // S
    -0.4, // T
    -1.5, // V
    -3.4, // W
    -2.3, // Y
];

// ── pKa values (EMBOSS) ────────────────────────────────────────

const PKA_NTERM: f64 = 9.69;
const PKA_CTERM: f64 = 2.34;
const PKA_D: f64 = 3.65;
const PKA_E: f64 = 4.25;
const PKA_C: f64 = 8.18;
const PKA_Y: f64 = 10.07;
const PKA_H: f64 = 6.00;
const PKA_K: f64 = 10.53;
const PKA_R: f64 = 12.48;

// ── Extinction coefficient constants (280 nm) ──────────────────

const EXT_TRP: f64 = 5500.0;
const EXT_TYR: f64 = 1490.0;
const EXT_CYSTINE: f64 = 125.0;

// ── Chou-Fasman propensities (1978) ────────────────────────────

/// (P_alpha, P_beta, P_turn) for each amino acid, indexed by aa_index.
const CHOU_FASMAN: [(f64, f64, f64); 20] = [
    (1.42, 0.83, 0.66), // A — strong helix former
    (0.70, 1.19, 1.19), // C
    (1.01, 0.54, 1.46), // D
    (1.51, 0.37, 0.74), // E — strong helix former
    (1.13, 1.38, 0.60), // F
    (0.57, 0.75, 1.56), // G — strong turn former
    (1.00, 0.87, 0.95), // H
    (1.08, 1.60, 0.47), // I
    (1.16, 0.74, 1.01), // K
    (1.21, 1.30, 0.59), // L
    (1.45, 1.05, 0.60), // M — strong helix former
    (0.67, 0.89, 1.56), // N — strong turn former
    (0.57, 0.55, 1.52), // P — strong turn former, helix breaker
    (1.11, 1.10, 0.98), // Q
    (0.98, 0.93, 0.95), // R
    (0.77, 0.75, 1.43), // S
    (0.83, 1.19, 0.96), // T
    (1.06, 1.70, 0.50), // V — strong strand former
    (1.08, 1.37, 0.96), // W
    (0.69, 1.47, 1.14), // Y
];

// ── GOR information values ─────────────────────────────────────

/// Simplified GOR single-residue information values (I_helix, I_strand, I_coil),
/// indexed by aa_index. Derived from published GOR I statistics.
const GOR_INFO: [(f64, f64, f64); 20] = [
    (0.36, -0.23, -0.13),  // A
    (-0.20, 0.17, 0.03),   // C
    (0.07, -0.42, 0.35),   // D
    (0.42, -0.37, -0.05),  // E
    (-0.09, 0.32, -0.23),  // F
    (-0.43, -0.18, 0.61),  // G
    (0.04, -0.09, 0.05),   // H
    (-0.06, 0.42, -0.36),  // I
    (0.13, -0.25, 0.12),   // K
    (0.21, 0.22, -0.43),   // L
    (0.36, 0.03, -0.39),   // M
    (-0.29, -0.18, 0.47),  // N
    (-0.42, -0.37, 0.79),  // P
    (0.18, -0.10, -0.08),  // Q
    (-0.01, -0.15, 0.16),  // R
    (-0.15, -0.07, 0.22),  // S
    (-0.11, 0.16, -0.05),  // T
    (-0.06, 0.52, -0.46),  // V
    (-0.02, 0.27, -0.25),  // W
    (-0.17, 0.31, -0.14),  // Y
];

/// GOR window half-width.
const GOR_HALF_WIDTH: usize = 8;

// ── Disorder propensity scale ──────────────────────────────────

/// Per-residue disorder propensity (order→disorder frequency ratio, mean-centered).
/// Positive values indicate disorder tendency. Indexed by aa_index.
const DISORDER_PROPENSITY: [f64; 20] = [
    -0.26, // A — slightly ordered
    -0.20, // C — ordered (disulfides)
    0.70,  // D — disordered (charged)
    0.55,  // E — disordered (charged)
    -0.60, // F — ordered (hydrophobic)
    0.16,  // G — slightly disordered (flexible)
    0.06,  // H — neutral
    -0.70, // I — ordered (hydrophobic)
    0.60,  // K — disordered (charged)
    -0.50, // L — ordered (hydrophobic)
    -0.14, // M — slightly ordered
    0.38,  // N — disordered (polar)
    0.55,  // P — disordered (rigid kink)
    0.45,  // Q — disordered (polar)
    0.30,  // R — disordered (charged)
    0.20,  // S — slightly disordered (polar)
    0.12,  // T — slightly disordered (polar)
    -0.60, // V — ordered (hydrophobic)
    -0.55, // W — ordered (hydrophobic, bulky)
    -0.45, // Y — ordered (hydrophobic)
];

// ── Public types ────────────────────────────────────────────────

/// Choice of hydrophobicity scale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HydrophobicityScale {
    /// Kyte-Doolittle (1982) — positive = hydrophobic.
    KyteDoolittle,
    /// Hopp-Woods (1981) — positive = hydrophilic.
    HoppWoods,
}

/// Amino acid composition of a protein sequence.
#[derive(Debug, Clone)]
pub struct AminoAcidComposition {
    /// Absolute count for each of the 20 standard amino acids (indexed by `aa_index`).
    pub counts: [usize; 20],
    /// Fraction (0.0–1.0) for each amino acid.
    pub fractions: [f64; 20],
    /// Sequence length.
    pub length: usize,
}

/// Molar extinction coefficient at 280 nm.
#[derive(Debug, Clone)]
pub struct ExtinctionCoefficient {
    /// Extinction coefficient assuming all cysteines are reduced (M⁻¹ cm⁻¹).
    pub reduced: f64,
    /// Extinction coefficient assuming all cysteines form cystines (M⁻¹ cm⁻¹).
    pub cystine: f64,
    /// Absorbance at 0.1% (1 mg/mL), reduced cysteines.
    pub abs_01_reduced: f64,
    /// Absorbance at 0.1% (1 mg/mL), cystine bridges.
    pub abs_01_cystine: f64,
}

/// Predicted secondary structure state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecondaryStructure {
    Helix,
    Strand,
    Coil,
}

/// Result of secondary structure prediction.
#[derive(Debug, Clone)]
pub struct SecondaryStructurePrediction {
    /// Per-residue predicted state.
    pub states: Vec<SecondaryStructure>,
    /// Per-residue helix score.
    pub helix_scores: Vec<f64>,
    /// Per-residue strand score.
    pub strand_scores: Vec<f64>,
    /// Per-residue coil score.
    pub coil_scores: Vec<f64>,
    /// Fraction of residues predicted as helix.
    pub helix_fraction: f64,
    /// Fraction of residues predicted as strand.
    pub strand_fraction: f64,
    /// Fraction of residues predicted as coil.
    pub coil_fraction: f64,
}

/// Result of intrinsic disorder prediction.
#[derive(Debug, Clone)]
pub struct DisorderPrediction {
    /// Per-residue disorder score (0.0–1.0).
    pub scores: Vec<f64>,
    /// Per-residue disorder call (true if score > 0.5).
    pub disordered: Vec<bool>,
    /// Fraction of residues predicted as disordered.
    pub disorder_fraction: f64,
}

// ── Amino acid composition ──────────────────────────────────────

/// Compute amino acid composition of a protein sequence.
///
/// # Errors
///
/// Returns an error if the sequence is empty or contains non-standard residues.
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::amino_acid_composition;
///
/// let comp = amino_acid_composition(b"ACDEFGHIKLMNPQRSTVWY").unwrap();
/// assert_eq!(comp.length, 20);
/// // Each AA appears once
/// for &c in comp.counts.iter() {
///     assert_eq!(c, 1);
/// }
/// ```
pub fn amino_acid_composition(seq: &[u8]) -> Result<AminoAcidComposition> {
    let norm = normalize_protein(seq)?;
    let mut counts = [0usize; 20];
    for &aa in &norm {
        counts[aa_index(aa).unwrap()] += 1;
    }
    let len = norm.len() as f64;
    let mut fractions = [0.0f64; 20];
    for i in 0..20 {
        fractions[i] = counts[i] as f64 / len;
    }
    Ok(AminoAcidComposition {
        counts,
        fractions,
        length: norm.len(),
    })
}

// ── Hydrophobicity ──────────────────────────────────────────────

/// Compute a hydrophobicity profile using a sliding window average.
///
/// # Arguments
///
/// * `seq` — Protein sequence (standard 20 AAs)
/// * `window` — Window size (must be odd and ≥ 1)
/// * `scale` — Hydrophobicity scale to use
///
/// # Errors
///
/// Returns an error if the sequence is empty, contains invalid residues,
/// or window is even or larger than the sequence.
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::{hydrophobicity_profile, HydrophobicityScale};
///
/// let profile = hydrophobicity_profile(b"IIIIIIIII", 3, HydrophobicityScale::KyteDoolittle).unwrap();
/// assert!((profile[1] - 4.5).abs() < 1e-10); // I = 4.5 on KD scale
/// ```
pub fn hydrophobicity_profile(
    seq: &[u8],
    window: usize,
    scale: HydrophobicityScale,
) -> Result<Vec<f64>> {
    let norm = normalize_protein(seq)?;
    if window == 0 || window % 2 == 0 {
        return Err(CyaneaError::InvalidInput(
            "window size must be odd and >= 1".to_string(),
        ));
    }
    if window > norm.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "window size {} exceeds sequence length {}",
            window,
            norm.len()
        )));
    }

    let table = match scale {
        HydrophobicityScale::KyteDoolittle => &KYTE_DOOLITTLE,
        HydrophobicityScale::HoppWoods => &HOPP_WOODS,
    };

    let values: Vec<f64> = norm
        .iter()
        .map(|&aa| table[aa_index(aa).unwrap()])
        .collect();

    let n = norm.len();
    let mut profile = Vec::with_capacity(n - window + 1);

    // Initial window sum
    let mut sum: f64 = values[..window].iter().sum();
    profile.push(sum / window as f64);

    // Slide
    for i in 1..=(n - window) {
        sum += values[i + window - 1] - values[i - 1];
        profile.push(sum / window as f64);
    }

    Ok(profile)
}

/// Compute the GRAVY (grand average of hydropathicity) score.
///
/// GRAVY is the mean Kyte-Doolittle hydropathy over the entire sequence.
/// Positive values indicate overall hydrophobic character.
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::gravy;
///
/// let g = gravy(b"IIIII").unwrap();
/// assert!((g - 4.5).abs() < 1e-10);
/// ```
pub fn gravy(seq: &[u8]) -> Result<f64> {
    let norm = normalize_protein(seq)?;
    let sum: f64 = norm
        .iter()
        .map(|&aa| KYTE_DOOLITTLE[aa_index(aa).unwrap()])
        .sum();
    Ok(sum / norm.len() as f64)
}

// ── Isoelectric point ───────────────────────────────────────────

/// Net charge at a given pH via Henderson-Hasselbalch.
fn net_charge(seq: &[u8], ph: f64) -> f64 {
    let mut charge = 0.0;

    // N-terminus (positive)
    charge += 1.0 / (1.0 + 10_f64.powf(ph - PKA_NTERM));
    // C-terminus (negative)
    charge -= 1.0 / (1.0 + 10_f64.powf(PKA_CTERM - ph));

    for &aa in seq {
        match aa {
            b'D' => charge -= 1.0 / (1.0 + 10_f64.powf(PKA_D - ph)),
            b'E' => charge -= 1.0 / (1.0 + 10_f64.powf(PKA_E - ph)),
            b'C' => charge -= 1.0 / (1.0 + 10_f64.powf(PKA_C - ph)),
            b'Y' => charge -= 1.0 / (1.0 + 10_f64.powf(PKA_Y - ph)),
            b'H' => charge += 1.0 / (1.0 + 10_f64.powf(ph - PKA_H)),
            b'K' => charge += 1.0 / (1.0 + 10_f64.powf(ph - PKA_K)),
            b'R' => charge += 1.0 / (1.0 + 10_f64.powf(ph - PKA_R)),
            _ => {}
        }
    }
    charge
}

/// Compute the isoelectric point (pI) of a protein sequence.
///
/// Uses bisection on the Henderson-Hasselbalch charge equation with
/// EMBOSS pKa values. Converges to |charge| < 0.001 (~47 iterations).
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::isoelectric_point;
///
/// let pi = isoelectric_point(b"DDDDD").unwrap();
/// assert!(pi < 4.0); // highly acidic
/// ```
pub fn isoelectric_point(seq: &[u8]) -> Result<f64> {
    let norm = normalize_protein(seq)?;
    let mut lo = 0.0_f64;
    let mut hi = 14.0_f64;

    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        let charge = net_charge(&norm, mid);
        if charge.abs() < 0.001 {
            return Ok(mid);
        }
        if charge > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok((lo + hi) / 2.0)
}

// ── Extinction coefficient ──────────────────────────────────────

/// Estimate the molar extinction coefficient at 280 nm.
///
/// Uses the Pace et al. method based on Trp, Tyr, and Cys content.
/// Returns both reduced (all Cys free) and cystine (all Cys paired) estimates.
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::extinction_coefficient;
///
/// let ec = extinction_coefficient(b"W").unwrap();
/// assert!((ec.reduced - 5500.0).abs() < 1e-10);
/// ```
pub fn extinction_coefficient(seq: &[u8]) -> Result<ExtinctionCoefficient> {
    let norm = normalize_protein(seq)?;
    let comp = amino_acid_composition(&norm)?;

    let n_trp = comp.counts[aa_index(b'W').unwrap()] as f64;
    let n_tyr = comp.counts[aa_index(b'Y').unwrap()] as f64;
    let n_cys = comp.counts[aa_index(b'C').unwrap()] as f64;

    let reduced = n_trp * EXT_TRP + n_tyr * EXT_TYR;
    let n_cystine = (n_cys / 2.0).floor();
    let cystine = reduced + n_cystine * EXT_CYSTINE;

    // Molecular weight for A280 0.1% calculation
    let mw = molecular_weight_from_seq(&norm);

    let abs_01_reduced = if mw > 0.0 { reduced / mw } else { 0.0 };
    let abs_01_cystine = if mw > 0.0 { cystine / mw } else { 0.0 };

    Ok(ExtinctionCoefficient {
        reduced,
        cystine,
        abs_01_reduced,
        abs_01_cystine,
    })
}

/// Molecular weight from raw normalized sequence bytes.
fn molecular_weight_from_seq(seq: &[u8]) -> f64 {
    // Average monoisotopic weights
    const WEIGHTS: [f64; 20] = [
        89.09,  // A
        121.16, // C
        133.10, // D
        147.13, // E
        165.19, // F
        75.03,  // G
        155.16, // H
        131.17, // I
        146.19, // K
        131.17, // L
        149.21, // M
        132.12, // N
        115.13, // P
        146.15, // Q
        174.20, // R
        105.09, // S
        119.12, // T
        117.15, // V
        204.23, // W
        181.19, // Y
    ];
    if seq.is_empty() {
        return 0.0;
    }
    let sum: f64 = seq
        .iter()
        .map(|&aa| WEIGHTS[aa_index(aa).unwrap()])
        .sum();
    sum - (seq.len() as f64 - 1.0) * 18.015
}

// ── Secondary structure: Chou-Fasman ────────────────────────────

/// Predict secondary structure using the Chou-Fasman algorithm (1978).
///
/// Nucleation/extension method:
/// 1. Find helix nuclei: 6-residue windows with ≥4 having P_alpha ≥ 1.03
/// 2. Find strand nuclei: 5-residue windows with ≥3 having P_beta ≥ 1.05
/// 3. Extend nuclei while tetrapeptide average ≥ 1.00
/// 4. Resolve overlaps: higher average propensity wins
///
/// # Errors
///
/// Returns an error if the sequence is shorter than 6 residues.
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::{chou_fasman, SecondaryStructure};
///
/// // Poly-alanine is a strong helix former
/// let pred = chou_fasman(b"AAAAAAAAAAAAAAAAAAAA").unwrap();
/// assert!(pred.helix_fraction > 0.5);
/// ```
pub fn chou_fasman(seq: &[u8]) -> Result<SecondaryStructurePrediction> {
    let norm = normalize_protein(seq)?;
    let n = norm.len();
    if n < 6 {
        return Err(CyaneaError::InvalidInput(
            "sequence must be at least 6 residues for Chou-Fasman prediction".to_string(),
        ));
    }

    let props: Vec<(f64, f64, f64)> = norm
        .iter()
        .map(|&aa| CHOU_FASMAN[aa_index(aa).unwrap()])
        .collect();

    // State assignment: 0 = unassigned, 1 = helix, 2 = strand
    let mut assignment = vec![0u8; n];

    // Step 1: Find helix nuclei (6-residue windows, ≥4 with P_alpha ≥ 1.03)
    let mut helix_regions: Vec<(usize, usize)> = Vec::new();
    for i in 0..=n.saturating_sub(6) {
        let count = (0..6).filter(|&j| props[i + j].0 >= 1.03).count();
        if count >= 4 {
            helix_regions.push((i, i + 5));
        }
    }

    // Merge overlapping helix nuclei
    let helix_regions = merge_regions(helix_regions);

    // Extend helix nuclei while tetrapeptide average P_alpha ≥ 1.00
    let helix_regions: Vec<(usize, usize)> = helix_regions
        .into_iter()
        .map(|(start, end)| extend_region(&props, start, end, n, |p| p.0, 1.00))
        .collect();

    // Step 2: Find strand nuclei (5-residue windows, ≥3 with P_beta ≥ 1.05)
    let mut strand_regions: Vec<(usize, usize)> = Vec::new();
    for i in 0..=n.saturating_sub(5) {
        let count = (0..5).filter(|&j| props[i + j].1 >= 1.05).count();
        if count >= 3 {
            strand_regions.push((i, i + 4));
        }
    }

    let strand_regions = merge_regions(strand_regions);
    let strand_regions: Vec<(usize, usize)> = strand_regions
        .into_iter()
        .map(|(start, end)| extend_region(&props, start, end, n, |p| p.1, 1.00))
        .collect();

    // Assign helices
    for &(start, end) in &helix_regions {
        for i in start..=end {
            assignment[i] = 1;
        }
    }

    // Assign strands, resolving overlaps
    for &(start, end) in &strand_regions {
        for i in start..=end {
            if assignment[i] == 1 {
                // Overlap: compare average propensities
                let alpha_avg: f64 = props[i].0;
                let beta_avg: f64 = props[i].1;
                if beta_avg > alpha_avg {
                    assignment[i] = 2;
                }
                // else keep helix
            } else {
                assignment[i] = 2;
            }
        }
    }

    // Build scores and states
    let mut helix_scores = Vec::with_capacity(n);
    let mut strand_scores = Vec::with_capacity(n);
    let mut coil_scores = Vec::with_capacity(n);
    let mut states = Vec::with_capacity(n);

    for i in 0..n {
        let (pa, pb, pt) = props[i];
        helix_scores.push(pa);
        strand_scores.push(pb);
        coil_scores.push(pt);
        states.push(match assignment[i] {
            1 => SecondaryStructure::Helix,
            2 => SecondaryStructure::Strand,
            _ => SecondaryStructure::Coil,
        });
    }

    let helix_count = states.iter().filter(|&&s| s == SecondaryStructure::Helix).count();
    let strand_count = states.iter().filter(|&&s| s == SecondaryStructure::Strand).count();
    let coil_count = n - helix_count - strand_count;

    Ok(SecondaryStructurePrediction {
        states,
        helix_scores,
        strand_scores,
        coil_scores,
        helix_fraction: helix_count as f64 / n as f64,
        strand_fraction: strand_count as f64 / n as f64,
        coil_fraction: coil_count as f64 / n as f64,
    })
}

/// Merge overlapping (start, end) regions into non-overlapping regions.
fn merge_regions(mut regions: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if regions.is_empty() {
        return regions;
    }
    regions.sort_by_key(|r| r.0);
    let mut merged = vec![regions[0]];
    for &(start, end) in &regions[1..] {
        let last = merged.last_mut().unwrap();
        if start <= last.1 + 1 {
            last.1 = last.1.max(end);
        } else {
            merged.push((start, end));
        }
    }
    merged
}

/// Extend a nucleated region while tetrapeptide average propensity ≥ threshold.
fn extend_region(
    props: &[(f64, f64, f64)],
    start: usize,
    end: usize,
    n: usize,
    accessor: fn(&(f64, f64, f64)) -> f64,
    threshold: f64,
) -> (usize, usize) {
    let mut s = start;
    let mut e = end;

    // Extend left
    while s >= 4 {
        let avg: f64 = (0..4).map(|j| accessor(&props[s - 4 + j])).sum::<f64>() / 4.0;
        if avg >= threshold {
            s -= 1;
        } else {
            break;
        }
    }

    // Extend right
    while e + 4 < n {
        let avg: f64 = (0..4).map(|j| accessor(&props[e + 1 + j])).sum::<f64>() / 4.0;
        if avg >= threshold {
            e += 1;
        } else {
            break;
        }
    }

    (s, e)
}

// ── Secondary structure: GOR ────────────────────────────────────

/// Predict secondary structure using the GOR (Garnier-Osguthorpe-Robson) method.
///
/// Applies single-residue GOR information values with a triangular window
/// (half-width = 8). Isolated single-residue predictions surrounded by
/// a different state are smoothed to the surrounding state.
///
/// # Errors
///
/// Returns an error if the sequence is empty or contains invalid residues.
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::gor;
///
/// let pred = gor(b"EEEEEEEEEEEEEEEEEEEE").unwrap();
/// // E (glutamate) has high helix propensity in GOR
/// assert!(pred.helix_fraction > pred.strand_fraction);
/// ```
pub fn gor(seq: &[u8]) -> Result<SecondaryStructurePrediction> {
    let norm = normalize_protein(seq)?;
    let n = norm.len();

    let indices: Vec<usize> = norm
        .iter()
        .map(|&aa| aa_index(aa).unwrap())
        .collect();

    let mut helix_scores = vec![0.0f64; n];
    let mut strand_scores = vec![0.0f64; n];
    let mut coil_scores = vec![0.0f64; n];

    // Compute windowed GOR scores with triangular weighting
    for i in 0..n {
        let mut h = 0.0;
        let mut e = 0.0;
        let mut c = 0.0;
        let mut weight_sum = 0.0;

        let w_start = if i >= GOR_HALF_WIDTH { i - GOR_HALF_WIDTH } else { 0 };
        let w_end = (i + GOR_HALF_WIDTH).min(n - 1);

        for j in w_start..=w_end {
            let dist = if j >= i { j - i } else { i - j };
            let weight = 1.0 - (dist as f64 / (GOR_HALF_WIDTH as f64 + 1.0));
            let (ih, ie, ic) = GOR_INFO[indices[j]];
            h += ih * weight;
            e += ie * weight;
            c += ic * weight;
            weight_sum += weight;
        }

        if weight_sum > 0.0 {
            helix_scores[i] = h / weight_sum;
            strand_scores[i] = e / weight_sum;
            coil_scores[i] = c / weight_sum;
        }
    }

    // Assign states by argmax
    let mut states: Vec<SecondaryStructure> = (0..n)
        .map(|i| {
            let h = helix_scores[i];
            let e = strand_scores[i];
            let c = coil_scores[i];
            if h >= e && h >= c {
                SecondaryStructure::Helix
            } else if e >= h && e >= c {
                SecondaryStructure::Strand
            } else {
                SecondaryStructure::Coil
            }
        })
        .collect();

    // Smooth: flip isolated single-residue H/E surrounded by a different state
    if n >= 3 {
        for i in 1..n - 1 {
            let prev = states[i - 1];
            let next = states[i + 1];
            if prev == next && states[i] != prev {
                states[i] = prev;
            }
        }
    }

    let helix_count = states.iter().filter(|&&s| s == SecondaryStructure::Helix).count();
    let strand_count = states.iter().filter(|&&s| s == SecondaryStructure::Strand).count();
    let coil_count = n - helix_count - strand_count;

    Ok(SecondaryStructurePrediction {
        states,
        helix_scores,
        strand_scores,
        coil_scores,
        helix_fraction: helix_count as f64 / n as f64,
        strand_fraction: strand_count as f64 / n as f64,
        coil_fraction: coil_count as f64 / n as f64,
    })
}

// ── Intrinsic disorder ─────────────────────────────────────────

/// Predict intrinsic disorder using windowed amino acid propensities.
///
/// Averages per-residue disorder propensities in a centered window and
/// applies a logistic sigmoid (k=5.0) to produce a 0–1 score.
/// Residues with score > 0.5 are predicted disordered.
///
/// # Arguments
///
/// * `seq` — Protein sequence
/// * `window` — Window size (must be odd and ≥ 1)
///
/// # Example
///
/// ```
/// use cyanea_seq::protein_properties::predict_disorder;
///
/// // Charged residues → high disorder
/// let pred = predict_disorder(b"KKKKKKKKDDDDDDDDD", 7).unwrap();
/// assert!(pred.disorder_fraction > 0.5);
/// ```
pub fn predict_disorder(seq: &[u8], window: usize) -> Result<DisorderPrediction> {
    let norm = normalize_protein(seq)?;
    let n = norm.len();

    if window == 0 || window % 2 == 0 {
        return Err(CyaneaError::InvalidInput(
            "window size must be odd and >= 1".to_string(),
        ));
    }
    if window > n {
        return Err(CyaneaError::InvalidInput(format!(
            "window size {} exceeds sequence length {}",
            window, n
        )));
    }

    let raw_props: Vec<f64> = norm
        .iter()
        .map(|&aa| DISORDER_PROPENSITY[aa_index(aa).unwrap()])
        .collect();

    let half = window / 2;
    let mut scores = Vec::with_capacity(n);

    for i in 0..n {
        let w_start = if i >= half { i - half } else { 0 };
        let w_end = (i + half).min(n - 1);
        let avg: f64 = raw_props[w_start..=w_end].iter().sum::<f64>()
            / (w_end - w_start + 1) as f64;

        // Logistic sigmoid: 1 / (1 + exp(-k * x))
        let score = 1.0 / (1.0 + (-5.0 * avg).exp());
        scores.push(score);
    }

    let disordered: Vec<bool> = scores.iter().map(|&s| s > 0.5).collect();
    let n_disordered = disordered.iter().filter(|&&d| d).count();

    Ok(DisorderPrediction {
        scores,
        disordered,
        disorder_fraction: n_disordered as f64 / n as f64,
    })
}

// ── Tests ───────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── normalize_protein ──

    #[test]
    fn normalize_protein_uppercase() {
        let r = normalize_protein(b"ACDEFGHIKLMNPQRSTVWY").unwrap();
        assert_eq!(r, b"ACDEFGHIKLMNPQRSTVWY");
    }

    #[test]
    fn normalize_protein_lowercase() {
        let r = normalize_protein(b"acdefghiklmnpqrstvwy").unwrap();
        assert_eq!(r, b"ACDEFGHIKLMNPQRSTVWY");
    }

    #[test]
    fn normalize_protein_invalid() {
        assert!(normalize_protein(b"ABCDE").is_err()); // B is not standard
    }

    // ── amino_acid_composition ──

    #[test]
    fn composition_all_alanine() {
        let comp = amino_acid_composition(b"AAAAAAAAAA").unwrap();
        assert_eq!(comp.counts[0], 10); // A at index 0
        assert!((comp.fractions[0] - 1.0).abs() < 1e-10);
        assert_eq!(comp.length, 10);
        // All others are 0
        for i in 1..20 {
            assert_eq!(comp.counts[i], 0);
        }
    }

    #[test]
    fn composition_each_aa_once() {
        let comp = amino_acid_composition(b"ACDEFGHIKLMNPQRSTVWY").unwrap();
        assert_eq!(comp.length, 20);
        for i in 0..20 {
            assert_eq!(comp.counts[i], 1);
            assert!((comp.fractions[i] - 0.05).abs() < 1e-10);
        }
    }

    #[test]
    fn composition_empty_error() {
        assert!(amino_acid_composition(b"").is_err());
    }

    // ── hydrophobicity_profile ──

    #[test]
    fn hydro_poly_i_kd() {
        // I = 4.5 on Kyte-Doolittle
        let profile =
            hydrophobicity_profile(b"IIIIIIIII", 3, HydrophobicityScale::KyteDoolittle).unwrap();
        for &v in &profile {
            assert!((v - 4.5).abs() < 1e-10);
        }
    }

    #[test]
    fn hydro_poly_r_kd() {
        // R = -4.5 on Kyte-Doolittle
        let profile =
            hydrophobicity_profile(b"RRRRRRRRR", 3, HydrophobicityScale::KyteDoolittle).unwrap();
        for &v in &profile {
            assert!((v - (-4.5)).abs() < 1e-10);
        }
    }

    #[test]
    fn hydro_hopp_woods_sign() {
        // D = 3.0 on Hopp-Woods (hydrophilic, positive)
        let profile =
            hydrophobicity_profile(b"DDDDDDDDD", 3, HydrophobicityScale::HoppWoods).unwrap();
        for &v in &profile {
            assert!(v > 0.0);
        }
    }

    #[test]
    fn hydro_even_window_error() {
        assert!(
            hydrophobicity_profile(b"AAAAAA", 4, HydrophobicityScale::KyteDoolittle).is_err()
        );
    }

    #[test]
    fn hydro_window_1_raw() {
        // window=1 should return raw per-residue values
        let profile =
            hydrophobicity_profile(b"AIV", 1, HydrophobicityScale::KyteDoolittle).unwrap();
        assert_eq!(profile.len(), 3);
        assert!((profile[0] - 1.8).abs() < 1e-10); // A
        assert!((profile[1] - 4.5).abs() < 1e-10); // I
        assert!((profile[2] - 4.2).abs() < 1e-10); // V
    }

    // ── gravy ──

    #[test]
    fn gravy_poly_i() {
        let g = gravy(b"IIIII").unwrap();
        assert!((g - 4.5).abs() < 1e-10);
    }

    #[test]
    fn gravy_mixed() {
        // A=1.8, R=-4.5 → mean = (1.8 + -4.5) / 2 = -1.35
        let g = gravy(b"AR").unwrap();
        assert!((g - (-1.35)).abs() < 1e-10);
    }

    #[test]
    fn gravy_empty_error() {
        assert!(gravy(b"").is_err());
    }

    // ── isoelectric_point ──

    #[test]
    fn pi_poly_d_acidic() {
        let pi = isoelectric_point(b"DDDDD").unwrap();
        assert!(pi < 3.5, "poly-D pI should be < 3.5, got {}", pi);
    }

    #[test]
    fn pi_poly_k_basic() {
        let pi = isoelectric_point(b"KKKKK").unwrap();
        assert!(pi > 10.0, "poly-K pI should be > 10.0, got {}", pi);
    }

    #[test]
    fn pi_poly_g_neutral() {
        let pi = isoelectric_point(b"GGGGG").unwrap();
        // Glycine has no charged side chains; pI near (9.69 + 2.34)/2 ≈ 6.0
        assert!(pi > 5.0 && pi < 7.0, "poly-G pI should be ~6.0, got {}", pi);
    }

    #[test]
    fn pi_single_residue() {
        // Should not panic
        let pi = isoelectric_point(b"A").unwrap();
        assert!(pi > 0.0 && pi < 14.0);
    }

    #[test]
    fn pi_known_protein() {
        // Insulin B chain: FVNQHLCGSHLVEALYLVCGERGFFYTPKT
        // Expected pI ≈ 6.9 (literature value varies 6.5-7.2)
        let pi = isoelectric_point(b"FVNQHLCGSHLVEALYLVCGERGFFYTPKT").unwrap();
        assert!(
            pi > 5.5 && pi < 8.5,
            "insulin B chain pI should be ~6.9, got {}",
            pi
        );
    }

    // ── extinction_coefficient ──

    #[test]
    fn ext_trp_only() {
        let ec = extinction_coefficient(b"W").unwrap();
        assert!((ec.reduced - 5500.0).abs() < 1e-10);
        assert!((ec.cystine - 5500.0).abs() < 1e-10); // no Cys
    }

    #[test]
    fn ext_tyr_only() {
        let ec = extinction_coefficient(b"Y").unwrap();
        assert!((ec.reduced - 1490.0).abs() < 1e-10);
    }

    #[test]
    fn ext_cystine_contribution() {
        // Two cysteines can form one cystine bond
        let ec = extinction_coefficient(b"CC").unwrap();
        assert!((ec.reduced - 0.0).abs() < 1e-10); // no Trp/Tyr
        assert!((ec.cystine - 125.0).abs() < 1e-10); // one cystine
    }

    #[test]
    fn ext_no_absorbers() {
        // No Trp, Tyr, or Cys
        let ec = extinction_coefficient(b"AAAAAA").unwrap();
        assert!((ec.reduced - 0.0).abs() < 1e-10);
        assert!((ec.cystine - 0.0).abs() < 1e-10);
    }

    // ── chou_fasman ──

    #[test]
    fn cf_poly_a_helix() {
        // Alanine has high P_alpha (1.42) — strong helix former
        let pred = chou_fasman(b"AAAAAAAAAAAAAAAAAAA").unwrap();
        assert!(
            pred.helix_fraction > 0.5,
            "poly-A should be mostly helix, got helix_fraction={}",
            pred.helix_fraction
        );
    }

    #[test]
    fn cf_poly_v_strand() {
        // Valine has high P_beta (1.70) — strong strand former
        let pred = chou_fasman(b"VVVVVVVVVVVVVVVVVVV").unwrap();
        assert!(
            pred.strand_fraction > 0.5,
            "poly-V should be mostly strand, got strand_fraction={}",
            pred.strand_fraction
        );
    }

    #[test]
    fn cf_poly_g_coil() {
        // Glycine has low P_alpha (0.57) and P_beta (0.75), high P_turn (1.56)
        let pred = chou_fasman(b"GGGGGGGGGGGGGGGGGGG").unwrap();
        assert!(
            pred.coil_fraction > 0.5,
            "poly-G should be mostly coil, got coil_fraction={}",
            pred.coil_fraction
        );
    }

    #[test]
    fn cf_short_error() {
        assert!(chou_fasman(b"AAAA").is_err());
    }

    // ── gor ──

    #[test]
    fn gor_helix_rich() {
        // E (glutamate) has highest GOR helix info (0.42)
        let pred = gor(b"EEEEEEEEEEEEEEEEEEEE").unwrap();
        assert!(
            pred.helix_fraction > pred.strand_fraction,
            "E-rich should be mostly helix: H={}, E={}",
            pred.helix_fraction,
            pred.strand_fraction
        );
    }

    #[test]
    fn gor_strand_rich() {
        // V (valine) has highest GOR strand info (0.52)
        let pred = gor(b"VVVVVVVVVVVVVVVVVVVV").unwrap();
        assert!(
            pred.strand_fraction > pred.helix_fraction,
            "V-rich should be mostly strand: E={}, H={}",
            pred.strand_fraction,
            pred.helix_fraction
        );
    }

    #[test]
    fn gor_output_length() {
        let pred = gor(b"ACDEFGHIKLMNPQRSTVWY").unwrap();
        assert_eq!(pred.states.len(), 20);
        assert_eq!(pred.helix_scores.len(), 20);
        assert_eq!(pred.strand_scores.len(), 20);
        assert_eq!(pred.coil_scores.len(), 20);
    }

    #[test]
    fn gor_fractions_sum() {
        let pred = gor(b"ACDEFGHIKLMNPQRSTVWY").unwrap();
        let sum = pred.helix_fraction + pred.strand_fraction + pred.coil_fraction;
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "fractions should sum to 1.0, got {}",
            sum
        );
    }

    // ── predict_disorder ──

    #[test]
    fn disorder_charged_stretch() {
        // Highly charged (K, D, E) → disordered
        let pred = predict_disorder(b"KKKDDDEEEKKKDDDEEEK", 7).unwrap();
        assert!(
            pred.disorder_fraction > 0.5,
            "charged stretch should be mostly disordered, got {}",
            pred.disorder_fraction
        );
    }

    #[test]
    fn disorder_hydrophobic_stretch() {
        // Hydrophobic (I, V, L, F, W) → ordered
        let pred = predict_disorder(b"IVLFWIVLFWIVLFWIVLFW", 7).unwrap();
        assert!(
            pred.disorder_fraction < 0.5,
            "hydrophobic stretch should be mostly ordered, got {}",
            pred.disorder_fraction
        );
    }

    #[test]
    fn disorder_scores_bounded() {
        let pred = predict_disorder(b"ACDEFGHIKLMNPQRSTVWY", 5).unwrap();
        for &s in &pred.scores {
            assert!(s >= 0.0 && s <= 1.0, "score {} out of [0,1]", s);
        }
    }

    #[test]
    fn disorder_window_zero_error() {
        assert!(predict_disorder(b"AAAA", 0).is_err());
    }
}
