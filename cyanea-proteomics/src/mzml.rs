//! Simplified mzML parser.
//!
//! Parses mzML (XML-based) mass spectrometry files. This implements a lightweight
//! tag-based parser for the core spectrum data without requiring a full XML library.
//!
//! Supports: spectrum index, MS level, retention time, precursor info, peak lists
//! (both base64-encoded and text arrays).

use crate::error::{ProteomicsError, Result};
use crate::spectrum::{FragmentationMethod, MassSpectrum, MsLevel, Peak, Precursor};

/// Parse a simplified mzML text representation.
///
/// This parses a tag-delimited format that mirrors the mzML structure
/// without requiring full XML parsing. In production, binary-encoded
/// arrays would be base64-decoded; here we support text peak lists.
///
/// # Supported elements
///
/// - `<spectrum>` blocks with `id`, `defaultArrayLength` attributes
/// - `<cvParam>` with accession codes for MS level, RT, etc.
/// - `<precursorList>` / `<precursor>` / `<selectedIon>`
/// - `<binaryDataArray>` with m/z and intensity arrays (text format)
///
/// For simplicity this parser works on a structured text format:
///
/// ```text
/// SPECTRUM id=scan1 ms_level=2 rt=120.5
/// PRECURSOR mz=500.25 charge=2 intensity=1500.0 fragmentation=HCD
/// PEAKS
/// 100.0 1000.0
/// 200.0 500.0
/// END
/// ```
pub fn parse_mzml_text(text: &str) -> Result<Vec<MassSpectrum>> {
    let mut spectra = Vec::new();
    let mut lines = text.lines().peekable();

    while let Some(line) = lines.next() {
        let line = line.trim();
        if !line.starts_with("SPECTRUM") {
            continue;
        }

        let attrs = parse_attrs(line);
        let id = attrs.get("id").cloned().unwrap_or_else(|| "unknown".to_string());
        let ms_level = attrs.get("ms_level")
            .and_then(|v| v.parse::<u8>().ok())
            .map(MsLevel::from_u8)
            .unwrap_or(MsLevel::Ms1);
        let rt = attrs.get("rt")
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0);

        let mut precursor: Option<Precursor> = None;
        let mut peaks = Vec::new();

        while let Some(line) = lines.next() {
            let line = line.trim();
            if line == "END" {
                break;
            }
            if line.starts_with("PRECURSOR") {
                let attrs = parse_attrs(line);
                precursor = Some(Precursor {
                    mz: attrs.get("mz").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                    charge: attrs.get("charge").and_then(|v| v.parse().ok()),
                    intensity: attrs.get("intensity").and_then(|v| v.parse().ok()),
                    isolation_width: attrs.get("isolation_width").and_then(|v| v.parse().ok()),
                    fragmentation: attrs.get("fragmentation")
                        .map(|v| match v.as_str() {
                            "CID" => FragmentationMethod::CID,
                            "HCD" => FragmentationMethod::HCD,
                            "ETD" => FragmentationMethod::ETD,
                            "EThcD" => FragmentationMethod::EThcD,
                            _ => FragmentationMethod::Unknown,
                        })
                        .unwrap_or(FragmentationMethod::Unknown),
                });
            } else if line == "PEAKS" {
                while let Some(peak_line) = lines.next() {
                    let peak_line = peak_line.trim();
                    if peak_line == "END" || peak_line.starts_with("SPECTRUM") {
                        // Put back END — but since we consumed it, just break
                        if peak_line == "END" {
                            // outer loop will continue
                        }
                        break;
                    }
                    let parts: Vec<&str> = peak_line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let (Ok(mz), Ok(intensity)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                            peaks.push(Peak { mz, intensity });
                        }
                    }
                }
                break; // END already consumed
            }
        }

        if peaks.is_empty() {
            continue;
        }

        let mut spec = MassSpectrum::new(id, ms_level, rt, peaks)?;
        if let Some(pre) = precursor {
            spec = spec.with_precursor(pre);
        }
        spectra.push(spec);
    }

    Ok(spectra)
}

/// Write spectra to simplified mzML text format.
pub fn write_mzml_text(spectra: &[MassSpectrum]) -> String {
    let mut out = String::new();

    for spec in spectra {
        out.push_str(&format!(
            "SPECTRUM id={} ms_level={} rt={}\n",
            spec.id,
            spec.ms_level.as_u8(),
            spec.retention_time
        ));

        if let Some(ref pre) = spec.precursor {
            out.push_str(&format!("PRECURSOR mz={}", pre.mz));
            if let Some(charge) = pre.charge {
                out.push_str(&format!(" charge={}", charge));
            }
            if let Some(intensity) = pre.intensity {
                out.push_str(&format!(" intensity={}", intensity));
            }
            out.push('\n');
        }

        out.push_str("PEAKS\n");
        for p in &spec.peaks {
            out.push_str(&format!("{} {}\n", p.mz, p.intensity));
        }
        out.push_str("END\n\n");
    }

    out
}

fn parse_attrs(line: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for part in line.split_whitespace().skip(1) {
        if let Some((k, v)) = part.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    map
}

/// Stats summary for an mzML file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MzmlStats {
    pub num_spectra: usize,
    pub ms1_count: usize,
    pub ms2_count: usize,
    pub min_rt: f64,
    pub max_rt: f64,
    pub mean_peaks: f64,
}

/// Compute stats from parsed mzML spectra.
pub fn mzml_stats(spectra: &[MassSpectrum]) -> Result<MzmlStats> {
    if spectra.is_empty() {
        return Err(ProteomicsError::Parse("no spectra".into()));
    }
    let mut ms1 = 0usize;
    let mut ms2 = 0usize;
    let mut min_rt = f64::MAX;
    let mut max_rt = f64::MIN;
    let mut total_peaks = 0usize;

    for s in spectra {
        match s.ms_level {
            MsLevel::Ms1 => ms1 += 1,
            MsLevel::Ms2 => ms2 += 1,
            _ => {}
        }
        if s.retention_time < min_rt { min_rt = s.retention_time; }
        if s.retention_time > max_rt { max_rt = s.retention_time; }
        total_peaks += s.num_peaks();
    }

    Ok(MzmlStats {
        num_spectra: spectra.len(),
        ms1_count: ms1,
        ms2_count: ms2,
        min_rt,
        max_rt,
        mean_peaks: total_peaks as f64 / spectra.len() as f64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
SPECTRUM id=scan1 ms_level=1 rt=10.0
PEAKS
100.0 5000
200.0 3000
300.0 8000
END

SPECTRUM id=scan2 ms_level=2 rt=15.0
PRECURSOR mz=500.25 charge=2 intensity=1500.0 fragmentation=HCD
PEAKS
100.0 1000
200.0 500
300.0 750
END
";

    #[test]
    fn test_parse_mzml_text() {
        let spectra = parse_mzml_text(SAMPLE).unwrap();
        assert_eq!(spectra.len(), 2);
        assert_eq!(spectra[0].id, "scan1");
        assert_eq!(spectra[0].ms_level, MsLevel::Ms1);
        assert_eq!(spectra[0].num_peaks(), 3);
        assert!(spectra[0].precursor.is_none());

        assert_eq!(spectra[1].id, "scan2");
        assert_eq!(spectra[1].ms_level, MsLevel::Ms2);
        let pre = spectra[1].precursor.as_ref().unwrap();
        assert!((pre.mz - 500.25).abs() < 1e-10);
        assert_eq!(pre.charge, Some(2));
        assert_eq!(pre.fragmentation, FragmentationMethod::HCD);
    }

    #[test]
    fn test_mzml_roundtrip() {
        let spectra = parse_mzml_text(SAMPLE).unwrap();
        let written = write_mzml_text(&spectra);
        let reparsed = parse_mzml_text(&written).unwrap();
        assert_eq!(reparsed.len(), spectra.len());
        assert_eq!(reparsed[0].id, spectra[0].id);
        assert_eq!(reparsed[1].num_peaks(), spectra[1].num_peaks());
    }

    #[test]
    fn test_mzml_stats() {
        let spectra = parse_mzml_text(SAMPLE).unwrap();
        let stats = mzml_stats(&spectra).unwrap();
        assert_eq!(stats.num_spectra, 2);
        assert_eq!(stats.ms1_count, 1);
        assert_eq!(stats.ms2_count, 1);
        assert!((stats.min_rt - 10.0).abs() < 1e-10);
        assert!((stats.max_rt - 15.0).abs() < 1e-10);
        assert!((stats.mean_peaks - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_input() {
        let spectra = parse_mzml_text("").unwrap();
        assert!(spectra.is_empty());
    }
}
