//! MGF (Mascot Generic Format) parser.
//!
//! Parses `.mgf` files containing MS/MS spectra. The format is text-based
//! with BEGIN IONS / END IONS blocks.

use crate::error::{ProteomicsError, Result};
use crate::spectrum::{FragmentationMethod, MassSpectrum, MsLevel, Peak, Precursor};

/// Parse an MGF-format string into a vector of MS2 spectra.
///
/// Each `BEGIN IONS` ... `END IONS` block becomes one [`MassSpectrum`].
///
/// # Format
///
/// ```text
/// BEGIN IONS
/// TITLE=spectrum_1
/// PEPMASS=523.25 1500.0
/// CHARGE=2+
/// RTINSECONDS=120.5
/// 100.05 1000
/// 200.10 500
/// 300.15 750
/// END IONS
/// ```
pub fn parse_mgf(text: &str) -> Result<Vec<MassSpectrum>> {
    let mut spectra = Vec::new();
    let mut in_block = false;
    let mut title = String::new();
    let mut pepmass = 0.0f64;
    let mut pep_intensity: Option<f64> = None;
    let mut charge: Option<i32> = None;
    let mut rt = 0.0f64;
    let mut peaks = Vec::new();
    let mut scan_counter = 0usize;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line == "BEGIN IONS" {
            in_block = true;
            title.clear();
            pepmass = 0.0;
            pep_intensity = None;
            charge = None;
            rt = 0.0;
            peaks.clear();
            continue;
        }

        if line == "END IONS" {
            if !in_block {
                return Err(ProteomicsError::Parse("END IONS without BEGIN IONS".into()));
            }
            in_block = false;
            scan_counter += 1;

            if peaks.is_empty() {
                continue; // skip empty spectra
            }

            let id = if title.is_empty() {
                format!("scan_{}", scan_counter)
            } else {
                title.clone()
            };

            let mut spec = MassSpectrum::new(id, MsLevel::Ms2, rt, peaks.clone())?;

            if pepmass > 0.0 {
                spec = spec.with_precursor(Precursor {
                    mz: pepmass,
                    charge,
                    intensity: pep_intensity,
                    isolation_width: None,
                    fragmentation: FragmentationMethod::Unknown,
                });
            }

            spectra.push(spec);
            continue;
        }

        if !in_block {
            continue;
        }

        // Parse header fields
        if let Some(val) = line.strip_prefix("TITLE=") {
            title = val.to_string();
        } else if let Some(val) = line.strip_prefix("PEPMASS=") {
            let parts: Vec<&str> = val.split_whitespace().collect();
            pepmass = parts[0].parse::<f64>()
                .map_err(|_| ProteomicsError::Parse(format!("invalid PEPMASS: {}", val)))?;
            if parts.len() > 1 {
                pep_intensity = parts[1].parse::<f64>().ok();
            }
        } else if let Some(val) = line.strip_prefix("CHARGE=") {
            let val = val.trim_end_matches('+').trim_end_matches('-');
            charge = val.parse::<i32>().ok();
        } else if let Some(val) = line.strip_prefix("RTINSECONDS=") {
            rt = val.parse::<f64>()
                .map_err(|_| ProteomicsError::Parse(format!("invalid RTINSECONDS: {}", val)))?;
        } else if line.contains('=') {
            // Unknown header field — skip
        } else {
            // Peak line: mz intensity
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                let mz = parts[0].parse::<f64>()
                    .map_err(|_| ProteomicsError::Parse(format!("invalid m/z: {}", parts[0])))?;
                let intensity = parts[1].parse::<f64>()
                    .map_err(|_| ProteomicsError::Parse(format!("invalid intensity: {}", parts[1])))?;
                peaks.push(Peak { mz, intensity });
            }
        }
    }

    if in_block {
        return Err(ProteomicsError::Parse("unterminated BEGIN IONS block".into()));
    }

    Ok(spectra)
}

/// Write spectra to MGF format.
pub fn write_mgf(spectra: &[MassSpectrum]) -> String {
    let mut out = String::new();

    for spec in spectra {
        out.push_str("BEGIN IONS\n");
        out.push_str(&format!("TITLE={}\n", spec.id));
        out.push_str(&format!("RTINSECONDS={}\n", spec.retention_time));

        if let Some(ref pre) = spec.precursor {
            if let Some(intensity) = pre.intensity {
                out.push_str(&format!("PEPMASS={} {}\n", pre.mz, intensity));
            } else {
                out.push_str(&format!("PEPMASS={}\n", pre.mz));
            }
            if let Some(charge) = pre.charge {
                out.push_str(&format!("CHARGE={}+\n", charge));
            }
        }

        for peak in &spec.peaks {
            out.push_str(&format!("{} {}\n", peak.mz, peak.intensity));
        }

        out.push_str("END IONS\n\n");
    }

    out
}

/// Stats from an MGF file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MgfStats {
    /// Total number of spectra.
    pub num_spectra: usize,
    /// Spectra with precursor charge annotation.
    pub with_charge: usize,
    /// Min precursor m/z.
    pub min_precursor_mz: f64,
    /// Max precursor m/z.
    pub max_precursor_mz: f64,
    /// Min retention time (seconds).
    pub min_rt: f64,
    /// Max retention time (seconds).
    pub max_rt: f64,
    /// Mean peaks per spectrum.
    pub mean_peaks: f64,
}

/// Compute summary statistics for parsed MGF spectra.
pub fn mgf_stats(spectra: &[MassSpectrum]) -> Result<MgfStats> {
    if spectra.is_empty() {
        return Err(ProteomicsError::Parse("no spectra".into()));
    }

    let mut with_charge = 0usize;
    let mut min_pre = f64::MAX;
    let mut max_pre = f64::MIN;
    let mut min_rt = f64::MAX;
    let mut max_rt = f64::MIN;
    let mut total_peaks = 0usize;

    for s in spectra {
        if let Some(ref pre) = s.precursor {
            if pre.charge.is_some() {
                with_charge += 1;
            }
            if pre.mz < min_pre {
                min_pre = pre.mz;
            }
            if pre.mz > max_pre {
                max_pre = pre.mz;
            }
        }
        if s.retention_time < min_rt {
            min_rt = s.retention_time;
        }
        if s.retention_time > max_rt {
            max_rt = s.retention_time;
        }
        total_peaks += s.num_peaks();
    }

    Ok(MgfStats {
        num_spectra: spectra.len(),
        with_charge,
        min_precursor_mz: if min_pre == f64::MAX { 0.0 } else { min_pre },
        max_precursor_mz: if max_pre == f64::MIN { 0.0 } else { max_pre },
        min_rt: if min_rt == f64::MAX { 0.0 } else { min_rt },
        max_rt: if max_rt == f64::MIN { 0.0 } else { max_rt },
        mean_peaks: total_peaks as f64 / spectra.len() as f64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_MGF: &str = "\
BEGIN IONS
TITLE=spectrum_1
PEPMASS=523.25 1500.0
CHARGE=2+
RTINSECONDS=120.5
100.05 1000
200.10 500
300.15 750
END IONS

BEGIN IONS
TITLE=spectrum_2
PEPMASS=789.42
CHARGE=3+
RTINSECONDS=180.0
150.00 200
250.00 800
350.00 400
450.00 100
END IONS
";

    #[test]
    fn test_parse_mgf() {
        let spectra = parse_mgf(SAMPLE_MGF).unwrap();
        assert_eq!(spectra.len(), 2);

        let s1 = &spectra[0];
        assert_eq!(s1.id, "spectrum_1");
        assert_eq!(s1.num_peaks(), 3);
        assert!((s1.retention_time - 120.5).abs() < 1e-10);
        let pre1 = s1.precursor.as_ref().unwrap();
        assert!((pre1.mz - 523.25).abs() < 1e-10);
        assert_eq!(pre1.charge, Some(2));
        assert!((pre1.intensity.unwrap() - 1500.0).abs() < 1e-10);

        let s2 = &spectra[1];
        assert_eq!(s2.id, "spectrum_2");
        assert_eq!(s2.num_peaks(), 4);
        let pre2 = s2.precursor.as_ref().unwrap();
        assert!((pre2.mz - 789.42).abs() < 1e-10);
        assert_eq!(pre2.charge, Some(3));
        assert!(pre2.intensity.is_none());
    }

    #[test]
    fn test_mgf_roundtrip() {
        let spectra = parse_mgf(SAMPLE_MGF).unwrap();
        let written = write_mgf(&spectra);
        let reparsed = parse_mgf(&written).unwrap();
        assert_eq!(reparsed.len(), spectra.len());
        assert_eq!(reparsed[0].id, spectra[0].id);
        assert_eq!(reparsed[0].num_peaks(), spectra[0].num_peaks());
    }

    #[test]
    fn test_mgf_empty_spectrum_skipped() {
        let mgf = "BEGIN IONS\nTITLE=empty\nEND IONS\n";
        let spectra = parse_mgf(mgf).unwrap();
        assert_eq!(spectra.len(), 0);
    }

    #[test]
    fn test_mgf_unterminated() {
        let mgf = "BEGIN IONS\nTITLE=bad\n100 500\n";
        assert!(parse_mgf(mgf).is_err());
    }

    #[test]
    fn test_mgf_no_title() {
        let mgf = "BEGIN IONS\nPEPMASS=500.0\n100 200\nEND IONS\n";
        let spectra = parse_mgf(mgf).unwrap();
        assert_eq!(spectra[0].id, "scan_1");
    }

    #[test]
    fn test_mgf_stats() {
        let spectra = parse_mgf(SAMPLE_MGF).unwrap();
        let stats = mgf_stats(&spectra).unwrap();
        assert_eq!(stats.num_spectra, 2);
        assert_eq!(stats.with_charge, 2);
        assert!((stats.min_precursor_mz - 523.25).abs() < 1e-10);
        assert!((stats.max_precursor_mz - 789.42).abs() < 1e-10);
        assert!((stats.mean_peaks - 3.5).abs() < 1e-10);
    }
}
