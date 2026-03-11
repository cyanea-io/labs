//! Peptide and protein digestion, amino acid masses, and fragment ion generation.
//!
//! Provides in-silico protein digestion (trypsin, chymotrypsin, etc.),
//! monoisotopic amino acid masses, and theoretical b/y ion series generation.

use crate::error::{ProteomicsError, Result};

/// Monoisotopic masses of the 20 standard amino acids (in Daltons).
pub fn amino_acid_mass(aa: u8) -> Option<f64> {
    match aa {
        b'G' => Some(57.02146),
        b'A' => Some(71.03711),
        b'V' => Some(99.06841),
        b'L' => Some(113.08406),
        b'I' => Some(113.08406),
        b'P' => Some(97.05276),
        b'F' => Some(147.06841),
        b'W' => Some(186.07931),
        b'M' => Some(131.04049),
        b'S' => Some(87.03203),
        b'T' => Some(101.04768),
        b'C' => Some(103.00919),
        b'Y' => Some(163.06333),
        b'H' => Some(137.05891),
        b'D' => Some(115.02694),
        b'E' => Some(129.04259),
        b'N' => Some(114.04293),
        b'Q' => Some(128.05858),
        b'K' => Some(128.09496),
        b'R' => Some(156.10111),
        _ => None,
    }
}

/// Proton mass (H+) used for m/z calculation.
pub const PROTON_MASS: f64 = 1.007276;
/// Water mass (H2O) added for peptide molecular weight.
pub const WATER_MASS: f64 = 18.010565;
/// Ammonia mass (NH3) for neutral loss calculation.
pub const AMMONIA_MASS: f64 = 17.026549;

/// Common post-translational modifications (mass shifts in Da).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Modification {
    /// Carbamidomethylation of cysteine (+57.021 Da).
    Carbamidomethyl,
    /// Oxidation of methionine (+15.995 Da).
    Oxidation,
    /// Phosphorylation of S/T/Y (+79.966 Da).
    Phospho,
    /// Acetylation of N-terminus or K (+42.011 Da).
    Acetyl,
    /// Deamidation of N/Q (+0.984 Da).
    Deamidation,
    /// TMT6/10/11/16plex tag (+229.163 Da).
    TMT,
    /// iTRAQ 4-plex (+144.102 Da).
    ITRAQ4,
    /// Custom modification.
    Custom(f64),
}

impl Modification {
    /// Get the mass shift for this modification.
    pub fn mass_shift(&self) -> f64 {
        match self {
            Modification::Carbamidomethyl => 57.02146,
            Modification::Oxidation => 15.99491,
            Modification::Phospho => 79.96633,
            Modification::Acetyl => 42.01057,
            Modification::Deamidation => 0.98402,
            Modification::TMT => 229.16293,
            Modification::ITRAQ4 => 144.10207,
            Modification::Custom(m) => *m,
        }
    }
}

/// A modified amino acid residue.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ModifiedResidue {
    /// Position in the peptide (0-based).
    pub position: usize,
    /// The modification applied.
    pub modification: Modification,
}

/// A peptide with optional modifications.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Peptide {
    /// Amino acid sequence (uppercase single-letter).
    pub sequence: Vec<u8>,
    /// Modifications on specific residues.
    pub modifications: Vec<ModifiedResidue>,
    /// Number of missed cleavages that produced this peptide.
    pub missed_cleavages: usize,
}

impl Peptide {
    /// Create a new unmodified peptide.
    pub fn new(sequence: &[u8]) -> Result<Self> {
        for &aa in sequence {
            if amino_acid_mass(aa).is_none() {
                return Err(ProteomicsError::Peptide(
                    format!("unknown amino acid: {}", aa as char),
                ));
            }
        }
        Ok(Self {
            sequence: sequence.to_vec(),
            modifications: Vec::new(),
            missed_cleavages: 0,
        })
    }

    /// Add a modification at a specific position.
    pub fn add_modification(&mut self, position: usize, modification: Modification) -> Result<()> {
        if position >= self.sequence.len() {
            return Err(ProteomicsError::Peptide("modification position out of range".into()));
        }
        self.modifications.push(ModifiedResidue { position, modification });
        Ok(())
    }

    /// Calculate the monoisotopic molecular weight (neutral mass).
    pub fn molecular_weight(&self) -> f64 {
        let base: f64 = self.sequence.iter()
            .map(|&aa| amino_acid_mass(aa).unwrap_or(0.0))
            .sum();
        let mods: f64 = self.modifications.iter()
            .map(|m| m.modification.mass_shift())
            .sum();
        base + mods + WATER_MASS
    }

    /// Calculate m/z for a given charge state.
    pub fn mz(&self, charge: i32) -> f64 {
        (self.molecular_weight() + charge as f64 * PROTON_MASS) / charge as f64
    }

    /// Length of the peptide.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the peptide is empty.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Sequence as a string.
    pub fn sequence_str(&self) -> String {
        String::from_utf8_lossy(&self.sequence).to_string()
    }
}

/// Fragment ion type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum IonType {
    /// b-ion (N-terminal fragment).
    B,
    /// y-ion (C-terminal fragment).
    Y,
    /// a-ion (b - CO).
    A,
    /// c-ion (for ETD fragmentation).
    C,
    /// z-ion (for ETD fragmentation).
    Z,
}

/// A theoretical fragment ion.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FragmentIon {
    /// Ion type (b, y, a, etc.).
    pub ion_type: IonType,
    /// Ion number (1-based, e.g., b3 = 3).
    pub number: usize,
    /// Charge state.
    pub charge: i32,
    /// Calculated m/z.
    pub mz: f64,
    /// Whether this includes a neutral loss.
    pub neutral_loss: Option<f64>,
}

/// Generate theoretical b and y ions for a peptide.
///
/// Returns fragment ions for charge states 1 through `max_charge`.
pub fn fragment_ions(peptide: &Peptide, max_charge: i32) -> Vec<FragmentIon> {
    let n = peptide.len();
    if n < 2 {
        return Vec::new();
    }

    // Precompute cumulative masses from N-terminus
    let mut cumulative = Vec::with_capacity(n);
    let mut running = 0.0;
    for (i, &aa) in peptide.sequence.iter().enumerate() {
        running += amino_acid_mass(aa).unwrap_or(0.0);
        // Add modification masses at this position
        for m in &peptide.modifications {
            if m.position == i {
                running += m.modification.mass_shift();
            }
        }
        cumulative.push(running);
    }

    let total_mass = cumulative[n - 1] + WATER_MASS;
    let mut ions = Vec::new();

    for i in 1..n {
        let b_neutral = cumulative[i - 1]; // b-ion neutral mass
        let y_neutral = total_mass - cumulative[i - 1]; // y-ion neutral mass

        for charge in 1..=max_charge {
            let c = charge as f64;

            // b-ion
            ions.push(FragmentIon {
                ion_type: IonType::B,
                number: i,
                charge,
                mz: (b_neutral + c * PROTON_MASS) / c,
                neutral_loss: None,
            });

            // y-ion
            ions.push(FragmentIon {
                ion_type: IonType::Y,
                number: n - i,
                charge,
                mz: (y_neutral + c * PROTON_MASS) / c,
                neutral_loss: None,
            });

            // a-ion (b - CO, 27.995 Da)
            ions.push(FragmentIon {
                ion_type: IonType::A,
                number: i,
                charge,
                mz: (b_neutral - 27.99491 + c * PROTON_MASS) / c,
                neutral_loss: None,
            });

            // b-ion with water loss
            ions.push(FragmentIon {
                ion_type: IonType::B,
                number: i,
                charge,
                mz: (b_neutral - WATER_MASS + c * PROTON_MASS) / c,
                neutral_loss: Some(WATER_MASS),
            });

            // y-ion with ammonia loss
            ions.push(FragmentIon {
                ion_type: IonType::Y,
                number: n - i,
                charge,
                mz: (y_neutral - AMMONIA_MASS + c * PROTON_MASS) / c,
                neutral_loss: Some(AMMONIA_MASS),
            });
        }
    }

    ions
}

/// Protease for in-silico digestion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Protease {
    /// Cleaves C-terminal to K/R (not before P).
    Trypsin,
    /// Cleaves C-terminal to K/R (allows before P).
    TrypsinP,
    /// Cleaves C-terminal to K only.
    LysC,
    /// Cleaves C-terminal to F/W/Y/L (not before P).
    Chymotrypsin,
    /// Cleaves C-terminal to D/E.
    AspN,
    /// Cleaves C-terminal to D.
    GluC,
}

impl Protease {
    /// Check if a cleavage occurs between residues at positions i and i+1.
    fn cleaves(&self, protein: &[u8], i: usize) -> bool {
        if i >= protein.len() - 1 {
            return false;
        }
        let current = protein[i];
        let next = protein[i + 1];

        match self {
            Protease::Trypsin => (current == b'K' || current == b'R') && next != b'P',
            Protease::TrypsinP => current == b'K' || current == b'R',
            Protease::LysC => current == b'K',
            Protease::Chymotrypsin => {
                (current == b'F' || current == b'W' || current == b'Y' || current == b'L')
                    && next != b'P'
            }
            Protease::AspN => next == b'D',
            Protease::GluC => current == b'D' || current == b'E',
        }
    }
}

/// Configuration for in-silico digestion.
#[derive(Debug, Clone)]
pub struct DigestConfig {
    /// Protease to use.
    pub protease: Protease,
    /// Maximum number of missed cleavages.
    pub max_missed_cleavages: usize,
    /// Minimum peptide length.
    pub min_length: usize,
    /// Maximum peptide length.
    pub max_length: usize,
    /// Minimum peptide mass (Da).
    pub min_mass: f64,
    /// Maximum peptide mass (Da).
    pub max_mass: f64,
}

impl Default for DigestConfig {
    fn default() -> Self {
        Self {
            protease: Protease::Trypsin,
            max_missed_cleavages: 2,
            min_length: 6,
            max_length: 50,
            min_mass: 400.0,
            max_mass: 5000.0,
        }
    }
}

/// Perform in-silico digestion of a protein sequence.
///
/// Returns all peptides satisfying the digest configuration constraints.
pub fn digest(protein: &[u8], config: &DigestConfig) -> Result<Vec<Peptide>> {
    if protein.is_empty() {
        return Err(ProteomicsError::Peptide("empty protein sequence".into()));
    }

    // Validate sequence
    for &aa in protein {
        if amino_acid_mass(aa).is_none() {
            return Err(ProteomicsError::Peptide(
                format!("unknown amino acid in protein: {}", aa as char),
            ));
        }
    }

    // Find cleavage sites
    let mut sites = Vec::new();
    for i in 0..protein.len() - 1 {
        if config.protease.cleaves(protein, i) {
            sites.push(i + 1); // cleavage after position i
        }
    }

    // Generate peptides with up to max_missed_cleavages
    let mut peptides = Vec::new();
    let mut starts: Vec<usize> = vec![0];
    starts.extend_from_slice(&sites);

    let mut ends: Vec<usize> = sites.clone();
    ends.push(protein.len());

    for i in 0..starts.len() {
        for mc in 0..=config.max_missed_cleavages {
            let end_idx = i + mc;
            if end_idx >= ends.len() {
                break;
            }
            let start = starts[i];
            let end = ends[end_idx];
            let seq = &protein[start..end];
            let len = seq.len();

            if len < config.min_length || len > config.max_length {
                continue;
            }

            let mut pep = Peptide::new(seq)?;
            pep.missed_cleavages = mc;
            let mass = pep.molecular_weight();

            if mass >= config.min_mass && mass <= config.max_mass {
                peptides.push(pep);
            }
        }
    }

    Ok(peptides)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_amino_acid_masses() {
        assert!((amino_acid_mass(b'G').unwrap() - 57.02146).abs() < 1e-5);
        assert!((amino_acid_mass(b'A').unwrap() - 71.03711).abs() < 1e-5);
        assert!(amino_acid_mass(b'X').is_none());
    }

    #[test]
    fn test_peptide_molecular_weight() {
        let pep = Peptide::new(b"PEPTIDE").unwrap();
        // P(97.05) + E(129.04) + P(97.05) + T(101.05) + I(113.08) + D(115.03) + E(129.04) + H2O
        let expected = 97.05276 + 129.04259 + 97.05276 + 101.04768 + 113.08406 + 115.02694 + 129.04259 + WATER_MASS;
        assert!((pep.molecular_weight() - expected).abs() < 0.01);
    }

    #[test]
    fn test_peptide_mz() {
        let pep = Peptide::new(b"PEPTIDE").unwrap();
        let mw = pep.molecular_weight();
        let mz2 = (mw + 2.0 * PROTON_MASS) / 2.0;
        assert!((pep.mz(2) - mz2).abs() < 1e-5);
    }

    #[test]
    fn test_modification() {
        let mut pep = Peptide::new(b"PEPTCIDE").unwrap();
        let mw_before = pep.molecular_weight();
        pep.add_modification(4, Modification::Carbamidomethyl).unwrap();
        let mw_after = pep.molecular_weight();
        assert!((mw_after - mw_before - 57.02146).abs() < 1e-4);
    }

    #[test]
    fn test_modification_out_of_range() {
        let mut pep = Peptide::new(b"PEP").unwrap();
        assert!(pep.add_modification(5, Modification::Oxidation).is_err());
    }

    #[test]
    fn test_fragment_ions_count() {
        let pep = Peptide::new(b"PEPTIDE").unwrap();
        let ions = fragment_ions(&pep, 1);
        // For a 7-residue peptide: 6 cleavage sites, each producing b, y, a, b-H2O, y-NH3
        // = 6 * 5 = 30 ions for charge 1
        assert_eq!(ions.len(), 30);
    }

    #[test]
    fn test_fragment_ions_b1_y6() {
        let pep = Peptide::new(b"PEPTIDE").unwrap();
        let ions = fragment_ions(&pep, 1);
        // b1 ion should be mass of P + proton
        let b1: Vec<&FragmentIon> = ions.iter().filter(|i| i.ion_type == IonType::B && i.number == 1 && i.neutral_loss.is_none()).collect();
        assert_eq!(b1.len(), 1);
        let expected_b1 = (amino_acid_mass(b'P').unwrap() + PROTON_MASS) / 1.0;
        assert!((b1[0].mz - expected_b1).abs() < 0.01);
    }

    #[test]
    fn test_trypsin_digest() {
        let config = DigestConfig {
            protease: Protease::Trypsin,
            max_missed_cleavages: 0,
            min_length: 1,
            max_length: 100,
            min_mass: 0.0,
            max_mass: 10000.0,
        };
        // PEPTIDEKSEQENCER => [PEPTIDEK, SEQENCER]
        let peptides = digest(b"PEPTIDEKSEQENCER", &config).unwrap();
        assert_eq!(peptides.len(), 2);
        assert_eq!(peptides[0].sequence_str(), "PEPTIDEK");
        assert_eq!(peptides[1].sequence_str(), "SEQENCER");
    }

    #[test]
    fn test_trypsin_proline_rule() {
        let config = DigestConfig {
            protease: Protease::Trypsin,
            max_missed_cleavages: 0,
            min_length: 1,
            max_length: 100,
            min_mass: 0.0,
            max_mass: 10000.0,
        };
        // Trypsin does NOT cleave before P: PEPTIDEKPANTHER should not split at K-P
        let peptides = digest(b"PEPTIDEKPANTHER", &config).unwrap();
        assert_eq!(peptides.len(), 1); // No cleavage at K-P
    }

    #[test]
    fn test_trypsinp_digest() {
        let config = DigestConfig {
            protease: Protease::TrypsinP,
            max_missed_cleavages: 0,
            min_length: 1,
            max_length: 100,
            min_mass: 0.0,
            max_mass: 10000.0,
        };
        // TrypsinP ignores proline rule
        let peptides = digest(b"PEPTIDEKPANTHER", &config).unwrap();
        assert_eq!(peptides.len(), 2);
    }

    #[test]
    fn test_missed_cleavages() {
        let config = DigestConfig {
            protease: Protease::Trypsin,
            max_missed_cleavages: 1,
            min_length: 1,
            max_length: 100,
            min_mass: 0.0,
            max_mass: 10000.0,
        };
        // AAKDDKCC => [AAK, DDK, CC] + missed: [AAKDDK, DDKCC]
        let peptides = digest(b"AAKDDKCC", &config).unwrap();
        assert_eq!(peptides.len(), 5);
    }

    #[test]
    fn test_lysc_digest() {
        let config = DigestConfig {
            protease: Protease::LysC,
            max_missed_cleavages: 0,
            min_length: 1,
            max_length: 100,
            min_mass: 0.0,
            max_mass: 10000.0,
        };
        // LysC cleaves after K only, not R
        let peptides = digest(b"AAKDDRCCAK", &config).unwrap();
        assert_eq!(peptides.len(), 2);
        assert_eq!(peptides[0].sequence_str(), "AAK");
        assert_eq!(peptides[1].sequence_str(), "DDRCCAK");
    }

    #[test]
    fn test_digest_mass_filter() {
        let config = DigestConfig {
            protease: Protease::Trypsin,
            max_missed_cleavages: 0,
            min_length: 1,
            max_length: 100,
            min_mass: 500.0, // Will filter short peptides
            max_mass: 10000.0,
        };
        let peptides = digest(b"AAKDDKCCDDEER", &config).unwrap();
        // Short peptides like AAK (~300 Da) should be filtered
        for p in &peptides {
            assert!(p.molecular_weight() >= 500.0);
        }
    }

    #[test]
    fn test_empty_protein() {
        let config = DigestConfig::default();
        assert!(digest(b"", &config).is_err());
    }

    #[test]
    fn test_invalid_protein() {
        let config = DigestConfig::default();
        assert!(digest(b"PEPTX1DE", &config).is_err());
    }

    #[test]
    fn test_modification_mass_shifts() {
        assert!((Modification::Carbamidomethyl.mass_shift() - 57.02146).abs() < 0.001);
        assert!((Modification::Oxidation.mass_shift() - 15.99491).abs() < 0.001);
        assert!((Modification::Phospho.mass_shift() - 79.96633).abs() < 0.001);
        assert!((Modification::TMT.mass_shift() - 229.16293).abs() < 0.001);
        assert!((Modification::Custom(42.0).mass_shift() - 42.0).abs() < 1e-10);
    }
}
