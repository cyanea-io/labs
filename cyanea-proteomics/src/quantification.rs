//! Protein quantification: label-free and labelled approaches.
//!
//! - **Spectral counting**: Count PSMs per protein
//! - **Intensity-based**: Sum precursor intensities (MS1-level)
//! - **TMT/iTRAQ**: Reporter ion quantification from MS2/MS3 spectra

use crate::search::Psm;
use crate::spectrum::MassSpectrum;
use crate::protein::ProteinGroup;

/// Quantification method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum QuantMethod {
    /// Spectral counting (number of PSMs).
    SpectralCount,
    /// Sum of precursor intensities.
    IntensitySum,
    /// Top-3 peptide intensity average.
    Top3,
}

/// Quantification result for a single protein group.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProteinQuant {
    /// Protein accession(s).
    pub accessions: Vec<String>,
    /// Raw quantification value.
    pub raw_value: f64,
    /// Normalized value (after normalization).
    pub normalized_value: f64,
    /// Number of quantified peptides.
    pub peptide_count: usize,
    /// Number of PSMs used.
    pub psm_count: usize,
}

/// Quantify proteins using spectral counting.
pub fn spectral_counting(
    groups: &[ProteinGroup],
    psms: &[Psm],
) -> Vec<ProteinQuant> {
    groups.iter().map(|g| {
        let all_peptides: Vec<&str> = g.unique_peptides.iter()
            .chain(g.shared_peptides.iter())
            .map(|s| s.as_str())
            .collect();

        let count = psms.iter()
            .filter(|p| all_peptides.contains(&p.peptide_sequence.as_str()))
            .count();

        ProteinQuant {
            accessions: g.accessions.clone(),
            raw_value: count as f64,
            normalized_value: count as f64,
            peptide_count: all_peptides.len(),
            psm_count: count,
        }
    }).collect()
}

/// Quantify proteins using precursor intensity (sum or top-3).
pub fn intensity_quantification(
    groups: &[ProteinGroup],
    psms: &[Psm],
    spectra: &[MassSpectrum],
    method: QuantMethod,
) -> Vec<ProteinQuant> {
    // Build spectrum lookup
    let spec_map: std::collections::HashMap<&str, &MassSpectrum> = spectra.iter()
        .map(|s| (s.id.as_str(), s))
        .collect();

    groups.iter().map(|g| {
        let all_peptides: Vec<&str> = g.unique_peptides.iter()
            .chain(g.shared_peptides.iter())
            .map(|s| s.as_str())
            .collect();

        // Collect intensities per peptide
        let mut peptide_intensities: std::collections::HashMap<String, f64> =
            std::collections::HashMap::new();

        for psm in psms {
            if !all_peptides.contains(&psm.peptide_sequence.as_str()) {
                continue;
            }
            if let Some(spec) = spec_map.get(psm.spectrum_id.as_str()) {
                let intensity = spec.precursor.as_ref()
                    .and_then(|p| p.intensity)
                    .unwrap_or(spec.tic);

                let entry = peptide_intensities.entry(psm.peptide_sequence.clone())
                    .or_insert(0.0);
                // Take max intensity per peptide across PSMs
                if intensity > *entry {
                    *entry = intensity;
                }
            }
        }

        let mut intensities: Vec<f64> = peptide_intensities.values().cloned().collect();
        intensities.sort_by(|a, b| b.partial_cmp(a).unwrap());

        let raw_value = match method {
            QuantMethod::Top3 => {
                let top = &intensities[..intensities.len().min(3)];
                if top.is_empty() { 0.0 } else { top.iter().sum::<f64>() / top.len() as f64 }
            }
            _ => intensities.iter().sum(),
        };

        let psm_count = psms.iter()
            .filter(|p| all_peptides.contains(&p.peptide_sequence.as_str()))
            .count();

        ProteinQuant {
            accessions: g.accessions.clone(),
            raw_value,
            normalized_value: raw_value,
            peptide_count: peptide_intensities.len(),
            psm_count,
        }
    }).collect()
}

/// TMT reporter ion masses (6-plex, 10-plex, 11-plex).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TmtPlex {
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
}

impl TmtPlex {
    /// Get reporter ion m/z values for this TMT plex.
    pub fn reporter_mzs(&self) -> Vec<f64> {
        match self {
            TmtPlex::Tmt6 => vec![126.1277, 127.1311, 128.1344, 129.1378, 130.1411, 131.1382],
            TmtPlex::Tmt10 => vec![
                126.1277, 127.1248, 127.1311, 128.1281, 128.1344,
                129.1315, 129.1378, 130.1348, 130.1411, 131.1382,
            ],
            TmtPlex::Tmt11 => {
                let mut mzs = TmtPlex::Tmt10.reporter_mzs();
                mzs.push(131.1445);
                mzs
            }
            TmtPlex::Tmt16 => vec![
                126.1277, 127.1248, 127.1311, 128.1281, 128.1344,
                129.1315, 129.1378, 130.1348, 130.1411, 131.1382,
                131.1445, 132.1416, 132.1479, 133.1449, 133.1512,
                134.1483,
            ],
        }
    }

    /// Number of channels.
    pub fn channels(&self) -> usize {
        self.reporter_mzs().len()
    }
}

/// TMT quantification result for a single spectrum.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TmtQuant {
    /// Spectrum identifier.
    pub spectrum_id: String,
    /// Reporter ion intensities (one per channel).
    pub channel_intensities: Vec<f64>,
    /// Isolation purity (fraction of signal from target ion).
    pub purity: f64,
}

/// Extract TMT reporter ion intensities from MS2/MS3 spectra.
pub fn quantify_tmt(
    spectra: &[MassSpectrum],
    plex: TmtPlex,
    tolerance: f64,
) -> Vec<TmtQuant> {
    let reporters = plex.reporter_mzs();

    spectra.iter().map(|spec| {
        let mut intensities = Vec::with_capacity(reporters.len());
        let mut total_reporter = 0.0f64;

        for &reporter_mz in &reporters {
            let intensity = spec.find_peak(reporter_mz, tolerance)
                .map(|p| p.intensity)
                .unwrap_or(0.0);
            intensities.push(intensity);
            total_reporter += intensity;
        }

        // Estimate purity: reporter signal / total signal in reporter region
        let reporter_region_signal: f64 = spec.peaks.iter()
            .filter(|p| p.mz >= 126.0 && p.mz <= 135.0)
            .map(|p| p.intensity)
            .sum();

        let purity = if reporter_region_signal > 0.0 {
            total_reporter / reporter_region_signal
        } else {
            0.0
        };

        TmtQuant {
            spectrum_id: spec.id.clone(),
            channel_intensities: intensities,
            purity,
        }
    }).collect()
}

/// Normalize quantification values across samples using median normalization.
pub fn median_normalize(quants: &mut [ProteinQuant]) {
    let mut values: Vec<f64> = quants.iter()
        .filter(|q| q.raw_value > 0.0)
        .map(|q| q.raw_value)
        .collect();

    if values.is_empty() {
        return;
    }

    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if values.len() % 2 == 0 {
        (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
    } else {
        values[values.len() / 2]
    };

    if median == 0.0 {
        return;
    }

    for q in quants.iter_mut() {
        q.normalized_value = q.raw_value / median;
    }
}

/// Log2-transform quantification values.
pub fn log2_transform(quants: &mut [ProteinQuant]) {
    for q in quants.iter_mut() {
        q.normalized_value = if q.raw_value > 0.0 {
            q.raw_value.log2()
        } else {
            0.0
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectrum::{MsLevel, Peak};

    fn make_group(accessions: &[&str], peptides: &[&str]) -> ProteinGroup {
        ProteinGroup {
            accessions: accessions.iter().map(|s| s.to_string()).collect(),
            unique_peptides: peptides.iter().map(|s| s.to_string()).collect(),
            shared_peptides: Vec::new(),
            psm_count: 0,
            best_score: 0.0,
            coverage: 0.0,
            is_decoy: false,
        }
    }

    fn make_psm(seq: &str, spec_id: &str, score: f64) -> Psm {
        Psm {
            spectrum_id: spec_id.to_string(),
            peptide_sequence: seq.to_string(),
            xcorr: 1.0,
            hyperscore: score,
            matched_b: 5,
            matched_y: 5,
            total_b: 6,
            total_y: 6,
            delta_mass: 0.0,
            charge: 2,
            is_decoy: false,
        }
    }

    #[test]
    fn test_spectral_counting() {
        let groups = vec![make_group(&["P1"], &["AAAK", "BBBK"])];
        let psms = vec![
            make_psm("AAAK", "s1", 20.0),
            make_psm("AAAK", "s2", 18.0),
            make_psm("BBBK", "s3", 15.0),
            make_psm("CCCK", "s4", 10.0), // not in group
        ];

        let quants = spectral_counting(&groups, &psms);
        assert_eq!(quants.len(), 1);
        assert_eq!(quants[0].raw_value, 3.0); // 3 PSMs match
        assert_eq!(quants[0].psm_count, 3);
    }

    #[test]
    fn test_tmt_reporter_channels() {
        assert_eq!(TmtPlex::Tmt6.channels(), 6);
        assert_eq!(TmtPlex::Tmt10.channels(), 10);
        assert_eq!(TmtPlex::Tmt11.channels(), 11);
        assert_eq!(TmtPlex::Tmt16.channels(), 16);
    }

    #[test]
    fn test_quantify_tmt() {
        let peaks: Vec<Peak> = TmtPlex::Tmt6.reporter_mzs().iter()
            .enumerate()
            .map(|(i, &mz)| Peak { mz, intensity: (i + 1) as f64 * 100.0 })
            .collect();

        let spec = MassSpectrum::new("s1", MsLevel::Ms2, 0.0, peaks).unwrap();
        let quants = quantify_tmt(&[spec], TmtPlex::Tmt6, 0.01);

        assert_eq!(quants.len(), 1);
        assert_eq!(quants[0].channel_intensities.len(), 6);
        assert!((quants[0].channel_intensities[0] - 100.0).abs() < 1e-5);
        assert!((quants[0].channel_intensities[5] - 600.0).abs() < 1e-5);
        assert!(quants[0].purity > 0.9);
    }

    #[test]
    fn test_median_normalize() {
        let mut quants = vec![
            ProteinQuant {
                accessions: vec!["P1".into()],
                raw_value: 100.0,
                normalized_value: 100.0,
                peptide_count: 2,
                psm_count: 3,
            },
            ProteinQuant {
                accessions: vec!["P2".into()],
                raw_value: 200.0,
                normalized_value: 200.0,
                peptide_count: 3,
                psm_count: 5,
            },
            ProteinQuant {
                accessions: vec!["P3".into()],
                raw_value: 300.0,
                normalized_value: 300.0,
                peptide_count: 1,
                psm_count: 2,
            },
        ];

        median_normalize(&mut quants);
        // Median is 200. Values should be divided by 200.
        assert!((quants[0].normalized_value - 0.5).abs() < 1e-10);
        assert!((quants[1].normalized_value - 1.0).abs() < 1e-10);
        assert!((quants[2].normalized_value - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_log2_transform() {
        let mut quants = vec![
            ProteinQuant {
                accessions: vec!["P1".into()],
                raw_value: 8.0,
                normalized_value: 0.0,
                peptide_count: 1,
                psm_count: 1,
            },
        ];
        log2_transform(&mut quants);
        assert!((quants[0].normalized_value - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_log2_zero() {
        let mut quants = vec![
            ProteinQuant {
                accessions: vec!["P1".into()],
                raw_value: 0.0,
                normalized_value: 0.0,
                peptide_count: 0,
                psm_count: 0,
            },
        ];
        log2_transform(&mut quants);
        assert_eq!(quants[0].normalized_value, 0.0);
    }
}
