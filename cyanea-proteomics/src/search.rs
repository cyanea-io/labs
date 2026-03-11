//! Peptide-spectrum matching (database search).
//!
//! Implements scoring functions for matching experimental MS/MS spectra
//! against theoretical fragment ions from candidate peptides.

use crate::error::Result;
use crate::peptide::{fragment_ions, Peptide, IonType};
use crate::spectrum::MassSpectrum;

/// A peptide-spectrum match (PSM).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Psm {
    /// Spectrum identifier.
    pub spectrum_id: String,
    /// Matched peptide sequence.
    pub peptide_sequence: String,
    /// XCorr-like cross-correlation score.
    pub xcorr: f64,
    /// Hyperscore (sum of matched ion intensities).
    pub hyperscore: f64,
    /// Number of matched b-ions.
    pub matched_b: usize,
    /// Number of matched y-ions.
    pub matched_y: usize,
    /// Total possible b-ions.
    pub total_b: usize,
    /// Total possible y-ions.
    pub total_y: usize,
    /// Delta mass (precursor observed - peptide theoretical, in Da).
    pub delta_mass: f64,
    /// Charge state used.
    pub charge: i32,
    /// Whether this is a decoy hit.
    pub is_decoy: bool,
}

impl Psm {
    /// Fraction of matched ions out of total possible.
    pub fn ion_fraction(&self) -> f64 {
        let total = self.total_b + self.total_y;
        if total == 0 {
            return 0.0;
        }
        (self.matched_b + self.matched_y) as f64 / total as f64
    }
}

/// Configuration for database search.
#[derive(Debug, Clone)]
pub struct SearchConfig {
    /// Fragment ion mass tolerance in Da.
    pub fragment_tolerance: f64,
    /// Precursor mass tolerance in Da.
    pub precursor_tolerance: f64,
    /// Maximum fragment charge to consider.
    pub max_fragment_charge: i32,
    /// Minimum number of matched ions for a valid PSM.
    pub min_matched_ions: usize,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            fragment_tolerance: 0.02,
            precursor_tolerance: 10.0, // Da (use 0.01 for high-res)
            max_fragment_charge: 2,
            min_matched_ions: 4,
        }
    }
}

/// Score a single spectrum against a candidate peptide.
///
/// Returns a PSM with XCorr and hyperscore if enough ions match.
pub fn score_peptide(
    spectrum: &MassSpectrum,
    peptide: &Peptide,
    charge: i32,
    config: &SearchConfig,
) -> Result<Option<Psm>> {
    // Check precursor mass tolerance
    if let Some(ref pre) = spectrum.precursor {
        let theo_mz = peptide.mz(charge);
        let delta = (pre.mz - theo_mz).abs();
        if delta > config.precursor_tolerance {
            return Ok(None);
        }
    }

    let ions = fragment_ions(peptide, config.max_fragment_charge);

    // Count matched ions and compute scores
    let mut matched_b = 0usize;
    let mut matched_y = 0usize;
    let mut total_b = 0usize;
    let mut total_y = 0usize;
    let mut matched_intensity_sum = 0.0f64;
    let mut xcorr_sum = 0.0f64;

    for ion in &ions {
        if ion.neutral_loss.is_some() {
            continue; // Only count primary ions for PSM scoring
        }

        match ion.ion_type {
            IonType::B => total_b += 1,
            IonType::Y => total_y += 1,
            _ => continue,
        }

        if let Some(peak) = spectrum.find_peak(ion.mz, config.fragment_tolerance) {
            match ion.ion_type {
                IonType::B => matched_b += 1,
                IonType::Y => matched_y += 1,
                _ => {}
            }
            matched_intensity_sum += peak.intensity;
            // XCorr contribution: intensity normalized by number of peaks
            xcorr_sum += peak.intensity / spectrum.tic.max(1.0);
        }
    }

    let total_matched = matched_b + matched_y;
    if total_matched < config.min_matched_ions {
        return Ok(None);
    }

    // Compute hyperscore: log of product of factorials of matched ions * intensity
    // Simplified: log(matched_b! * matched_y! * intensity_sum)
    let log_b_fact = log_factorial(matched_b);
    let log_y_fact = log_factorial(matched_y);
    let hyperscore = log_b_fact + log_y_fact + matched_intensity_sum.ln().max(0.0);

    let delta_mass = if let Some(ref pre) = spectrum.precursor {
        pre.mz * charge as f64 - (charge as f64 * crate::peptide::PROTON_MASS) - peptide.molecular_weight()
    } else {
        0.0
    };

    Ok(Some(Psm {
        spectrum_id: spectrum.id.clone(),
        peptide_sequence: peptide.sequence_str(),
        xcorr: xcorr_sum,
        hyperscore,
        matched_b,
        matched_y,
        total_b,
        total_y,
        delta_mass,
        charge,
        is_decoy: false,
    }))
}

/// Search a spectrum against a database of peptides.
///
/// Returns the best-scoring PSM (if any passes filters).
pub fn search_spectrum(
    spectrum: &MassSpectrum,
    database: &[Peptide],
    config: &SearchConfig,
) -> Result<Option<Psm>> {
    let charge = spectrum.precursor.as_ref()
        .and_then(|p| p.charge)
        .unwrap_or(2);

    let mut best: Option<Psm> = None;

    for peptide in database {
        if let Some(psm) = score_peptide(spectrum, peptide, charge, config)? {
            if best.as_ref().map_or(true, |b| psm.hyperscore > b.hyperscore) {
                best = Some(psm);
            }
        }
    }

    Ok(best)
}

/// Search all spectra against a peptide database.
///
/// Returns the best PSM per spectrum.
pub fn search_all(
    spectra: &[MassSpectrum],
    database: &[Peptide],
    config: &SearchConfig,
) -> Result<Vec<Psm>> {
    let mut psms = Vec::new();
    for spectrum in spectra {
        if spectrum.ms_level != crate::spectrum::MsLevel::Ms2 {
            continue;
        }
        if let Some(psm) = search_spectrum(spectrum, database, config)? {
            psms.push(psm);
        }
    }
    Ok(psms)
}

/// Generate a reversed-sequence decoy database from target peptides.
pub fn generate_decoys(targets: &[Peptide]) -> Vec<Peptide> {
    targets.iter().map(|p| {
        let mut rev = p.sequence.clone();
        // Reverse all but last AA (keep tryptic C-terminus)
        if rev.len() > 1 {
            let last = rev.len() - 1;
            rev[..last].reverse();
        }
        Peptide {
            sequence: rev,
            modifications: Vec::new(),
            missed_cleavages: p.missed_cleavages,
        }
    }).collect()
}

fn log_factorial(n: usize) -> f64 {
    (1..=n).map(|i| (i as f64).ln()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectrum::{Peak, MsLevel, Precursor, FragmentationMethod};

    fn make_test_spectrum(peptide: &Peptide, charge: i32) -> MassSpectrum {
        // Generate theoretical ions and create a spectrum from some of them
        let ions = fragment_ions(peptide, 1);
        let peaks: Vec<Peak> = ions.iter()
            .filter(|i| i.neutral_loss.is_none() && (i.ion_type == IonType::B || i.ion_type == IonType::Y))
            .map(|i| Peak { mz: i.mz, intensity: 1000.0 })
            .collect();

        let pre_mz = peptide.mz(charge);
        MassSpectrum::new("test_scan", MsLevel::Ms2, 100.0, peaks)
            .unwrap()
            .with_precursor(Precursor {
                mz: pre_mz,
                charge: Some(charge),
                intensity: Some(5000.0),
                isolation_width: None,
                fragmentation: FragmentationMethod::HCD,
            })
    }

    #[test]
    fn test_score_perfect_match() {
        let pep = Peptide::new(b"PEPTIDE").unwrap();
        let spec = make_test_spectrum(&pep, 2);
        let config = SearchConfig::default();

        let psm = score_peptide(&spec, &pep, 2, &config).unwrap();
        assert!(psm.is_some());
        let psm = psm.unwrap();
        assert_eq!(psm.peptide_sequence, "PEPTIDE");
        assert!(psm.matched_b > 0);
        assert!(psm.matched_y > 0);
        assert!(psm.hyperscore > 0.0);
        assert!(psm.ion_fraction() > 0.0);
    }

    #[test]
    fn test_no_match_wrong_precursor() {
        let pep1 = Peptide::new(b"PEPTIDE").unwrap();
        let pep2 = Peptide::new(b"ACDEFGHIK").unwrap();
        let spec = make_test_spectrum(&pep1, 2);
        let config = SearchConfig {
            precursor_tolerance: 1.0, // tight tolerance
            ..SearchConfig::default()
        };

        let psm = score_peptide(&spec, &pep2, 2, &config).unwrap();
        assert!(psm.is_none());
    }

    #[test]
    fn test_search_spectrum() {
        let target = Peptide::new(b"PEPTIDE").unwrap();
        let decoy = Peptide::new(b"ACDEFGHIK").unwrap();
        let spec = make_test_spectrum(&target, 2);
        let config = SearchConfig {
            precursor_tolerance: 50.0,
            ..SearchConfig::default()
        };

        let database = vec![target.clone(), decoy];
        let best = search_spectrum(&spec, &database, &config).unwrap();
        assert!(best.is_some());
        assert_eq!(best.unwrap().peptide_sequence, "PEPTIDE");
    }

    #[test]
    fn test_generate_decoys() {
        let targets = vec![
            Peptide::new(b"PEPTIDEK").unwrap(),
            Peptide::new(b"SEQENCER").unwrap(),
        ];
        let decoys = generate_decoys(&targets);
        assert_eq!(decoys.len(), 2);
        // Last AA preserved, rest reversed
        assert_eq!(decoys[0].sequence_str(), "EDITPEPK"); // PEPTIDE reversed + K
        assert_eq!(decoys[1].sequence_str(), "ECNEQESR"); // SEQENCE reversed + R
    }

    #[test]
    fn test_search_all_skips_ms1() {
        let pep = Peptide::new(b"PEPTIDE").unwrap();
        let ms1 = MassSpectrum::new("ms1", MsLevel::Ms1, 10.0, vec![Peak { mz: 100.0, intensity: 1000.0 }]).unwrap();
        let ms2 = make_test_spectrum(&pep, 2);
        let config = SearchConfig::default();

        let psms = search_all(&[ms1, ms2], &[pep], &config).unwrap();
        assert_eq!(psms.len(), 1);
    }

    #[test]
    fn test_psm_ion_fraction() {
        let psm = Psm {
            spectrum_id: "s1".into(),
            peptide_sequence: "TEST".into(),
            xcorr: 1.0,
            hyperscore: 10.0,
            matched_b: 3,
            matched_y: 4,
            total_b: 6,
            total_y: 6,
            delta_mass: 0.0,
            charge: 2,
            is_decoy: false,
        };
        assert!((psm.ion_fraction() - 7.0 / 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_log_factorial() {
        assert!((log_factorial(0) - 0.0).abs() < 1e-10);
        assert!((log_factorial(1) - 0.0).abs() < 1e-10);
        assert!((log_factorial(5) - (120.0f64).ln()).abs() < 1e-10);
    }
}
