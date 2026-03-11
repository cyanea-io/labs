//! FCS (Flow Cytometry Standard) file parser.
//!
//! Supports FCS 2.0, 3.0, and 3.1 formats. Parses the HEADER, TEXT,
//! and DATA segments to extract parameter metadata and event data.
//!
//! # Format overview
//!
//! An FCS file has three segments:
//! - **HEADER** (bytes 0-57): version string, offsets to TEXT and DATA
//! - **TEXT**: key-value pairs delimited by a separator character
//! - **DATA**: binary or ASCII event data (one value per parameter per event)

use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

/// Parsed FCS file.
#[derive(Debug, Clone)]
pub struct FcsFile {
    /// FCS version string (e.g., "FCS3.1").
    pub version: String,
    /// TEXT segment key-value pairs (all keys uppercased).
    pub text: HashMap<String, String>,
    /// Parameter definitions.
    pub parameters: Vec<FcsParameter>,
    /// Event data: `events[event_idx][param_idx]`.
    pub events: Vec<Vec<f64>>,
}

/// A single FCS parameter (channel).
#[derive(Debug, Clone)]
pub struct FcsParameter {
    /// Short name ($PnN).
    pub name: String,
    /// Long name / stain ($PnS), if present.
    pub stain: Option<String>,
    /// Bit width ($PnB).
    pub bits: u32,
    /// Range ($PnR).
    pub range: f64,
    /// Amplification type ($PnE): (decades, log0). (0,0) = linear.
    pub amplification: (f64, f64),
}

/// Summary statistics for an FCS file.
#[derive(Debug, Clone)]
pub struct FcsStats {
    pub version: String,
    pub n_events: usize,
    pub n_parameters: usize,
    pub parameter_names: Vec<String>,
}

impl FcsFile {
    /// Number of events.
    pub fn n_events(&self) -> usize {
        self.events.len()
    }

    /// Number of parameters.
    pub fn n_parameters(&self) -> usize {
        self.parameters.len()
    }

    /// Get summary statistics.
    pub fn stats(&self) -> FcsStats {
        FcsStats {
            version: self.version.clone(),
            n_events: self.n_events(),
            n_parameters: self.n_parameters(),
            parameter_names: self.parameters.iter().map(|p| p.name.clone()).collect(),
        }
    }

    /// Get a TEXT keyword value.
    pub fn keyword(&self, key: &str) -> Option<&str> {
        self.text.get(&key.to_uppercase()).map(|s| s.as_str())
    }

    /// Get all event values for a parameter by index.
    pub fn parameter_data(&self, param_idx: usize) -> Option<Vec<f64>> {
        if param_idx >= self.n_parameters() {
            return None;
        }
        Some(self.events.iter().map(|e| e[param_idx]).collect())
    }

    /// Get all event values for a parameter by name.
    pub fn parameter_data_by_name(&self, name: &str) -> Option<Vec<f64>> {
        let idx = self.parameters.iter().position(|p| p.name == name)?;
        self.parameter_data(idx)
    }
}

/// Parse an FCS file from binary data.
pub fn parse_fcs(data: &[u8]) -> Result<FcsFile> {
    if data.len() < 58 {
        return Err(CyaneaError::Parse("FCS file too short (< 58 bytes)".into()));
    }

    // HEADER: bytes 0-5 = version, 6-9 = spaces
    let version = std::str::from_utf8(&data[0..6])
        .map_err(|_| CyaneaError::Parse("Invalid FCS version string".into()))?
        .trim()
        .to_string();

    if !version.starts_with("FCS") {
        return Err(CyaneaError::Parse(format!("Not an FCS file: '{}'", version)));
    }

    // Offsets (ASCII integers, right-justified in 8-byte fields)
    let text_start = parse_offset(&data[10..18])? as usize;
    let text_end = parse_offset(&data[18..26])? as usize;
    let data_start_hdr = parse_offset(&data[26..34])? as usize;
    let data_end_hdr = parse_offset(&data[34..42])? as usize;

    if text_start >= data.len() || text_end >= data.len() || text_start > text_end {
        return Err(CyaneaError::Parse("Invalid TEXT segment offsets".into()));
    }

    // Parse TEXT segment
    let text_bytes = &data[text_start..=text_end];
    let text_str = std::str::from_utf8(text_bytes)
        .map_err(|_| CyaneaError::Parse("Invalid UTF-8 in TEXT segment".into()))?;

    let text = parse_text_segment(text_str)?;

    // Resolve DATA segment offsets (TEXT keywords override header for large files)
    let data_begin = if data_start_hdr == 0 {
        text.get("$BEGINDATA")
            .and_then(|v| v.trim().parse::<usize>().ok())
            .ok_or_else(|| CyaneaError::Parse("Missing DATA offset".into()))?
    } else {
        data_start_hdr
    };
    let data_end = if data_end_hdr == 0 {
        text.get("$ENDDATA")
            .and_then(|v| v.trim().parse::<usize>().ok())
            .ok_or_else(|| CyaneaError::Parse("Missing DATA end offset".into()))?
    } else {
        data_end_hdr
    };

    // Parse parameters
    let n_params = text.get("$PAR")
        .and_then(|v| v.trim().parse::<usize>().ok())
        .ok_or_else(|| CyaneaError::Parse("Missing $PAR keyword".into()))?;

    let n_events = text.get("$TOT")
        .and_then(|v| v.trim().parse::<usize>().ok())
        .ok_or_else(|| CyaneaError::Parse("Missing $TOT keyword".into()))?;

    let mut parameters = Vec::with_capacity(n_params);
    for i in 1..=n_params {
        let name = text.get(&format!("$P{}N", i))
            .cloned()
            .unwrap_or_else(|| format!("P{}", i));
        let stain = text.get(&format!("$P{}S", i)).cloned();
        let bits = text.get(&format!("$P{}B", i))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(32);
        let range = text.get(&format!("$P{}R", i))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(262144.0);
        let amplification = text.get(&format!("$P{}E", i))
            .map(|v| parse_amplification(v))
            .unwrap_or((0.0, 0.0));

        parameters.push(FcsParameter { name, stain, bits, range, amplification });
    }

    // Parse DATA segment
    let datatype = text.get("$DATATYPE")
        .map(|s| s.trim().to_uppercase())
        .unwrap_or_else(|| "F".into());

    let byteord = text.get("$BYTEORD")
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "1,2,3,4".into());

    let little_endian = byteord.starts_with('1');

    let events = if data_begin < data.len() && data_end < data.len() && data_begin <= data_end {
        let data_bytes = &data[data_begin..=data_end.min(data.len() - 1)];
        parse_data_segment(data_bytes, &datatype, &parameters, n_events, little_endian)?
    } else {
        Vec::new()
    };

    Ok(FcsFile {
        version,
        text,
        parameters,
        events,
    })
}

/// Write a minimal FCS 3.1 file from event data.
pub fn write_fcs(fcs: &FcsFile) -> Result<Vec<u8>> {
    let n_params = fcs.parameters.len();
    let n_events = fcs.events.len();

    // Build TEXT segment
    let mut text_pairs: Vec<(String, String)> = Vec::new();
    text_pairs.push(("$BEGINANALYSIS".into(), "0".into()));
    text_pairs.push(("$ENDANALYSIS".into(), "0".into()));
    text_pairs.push(("$DATATYPE".into(), "F".into()));
    text_pairs.push(("$MODE".into(), "L".into()));
    text_pairs.push(("$BYTEORD".into(), "1,2,3,4".into()));
    text_pairs.push(("$PAR".into(), n_params.to_string()));
    text_pairs.push(("$TOT".into(), n_events.to_string()));

    for (i, p) in fcs.parameters.iter().enumerate() {
        let idx = i + 1;
        text_pairs.push((format!("$P{}N", idx), p.name.clone()));
        if let Some(ref s) = p.stain {
            text_pairs.push((format!("$P{}S", idx), s.clone()));
        }
        text_pairs.push((format!("$P{}B", idx), p.bits.to_string()));
        text_pairs.push((format!("$P{}R", idx), format!("{}", p.range as u64)));
        text_pairs.push((format!("$P{}E", idx), format!("{},{}", p.amplification.0, p.amplification.1)));
    }

    let sep = '/';
    let mut text_str = String::new();
    text_str.push(sep);
    for (k, v) in &text_pairs {
        text_str.push_str(k);
        text_str.push(sep);
        text_str.push_str(v);
        text_str.push(sep);
    }

    let text_bytes = text_str.as_bytes();

    // Calculate offsets
    let text_start: usize = 256; // standard offset for TEXT
    let text_end = text_start + text_bytes.len() - 1;
    let data_start = text_end + 1;
    let data_size = n_events * n_params * 4; // f32
    let data_end = data_start + data_size - 1;

    // Add data offset keywords
    let mut full_text = String::new();
    full_text.push(sep);
    text_pairs.push(("$BEGINDATA".into(), data_start.to_string()));
    text_pairs.push(("$ENDDATA".into(), data_end.to_string()));
    for (k, v) in &text_pairs {
        full_text.push_str(k);
        full_text.push(sep);
        full_text.push_str(v);
        full_text.push(sep);
    }

    let text_bytes = full_text.as_bytes();
    let text_end = text_start + text_bytes.len() - 1;
    let data_start = text_end + 1;
    let data_end = if data_size > 0 { data_start + data_size - 1 } else { data_start };

    // Build output
    let total_size = data_end + 1;
    let mut out = vec![0u8; total_size.max(256)];

    // HEADER
    out[0..6].copy_from_slice(b"FCS3.1");
    out[6..10].copy_from_slice(b"    ");
    write_offset(&mut out[10..18], text_start as u64);
    write_offset(&mut out[18..26], text_end as u64);
    write_offset(&mut out[26..34], data_start as u64);
    write_offset(&mut out[34..42], data_end as u64);
    write_offset(&mut out[42..50], 0); // ANALYSIS start
    write_offset(&mut out[50..58], 0); // ANALYSIS end

    // TEXT
    if text_start + text_bytes.len() <= out.len() {
        out[text_start..text_start + text_bytes.len()].copy_from_slice(text_bytes);
    }

    // DATA (float32, little-endian)
    let mut offset = data_start;
    for event in &fcs.events {
        for &val in event {
            if offset + 4 <= out.len() {
                let bytes = (val as f32).to_le_bytes();
                out[offset..offset + 4].copy_from_slice(&bytes);
            }
            offset += 4;
        }
    }

    Ok(out)
}

/// Compute basic statistics for FCS data.
pub fn fcs_stats(fcs: &FcsFile) -> FcsStats {
    fcs.stats()
}

// --- Internal helpers ---

fn parse_offset(field: &[u8]) -> Result<u64> {
    let s = std::str::from_utf8(field)
        .map_err(|_| CyaneaError::Parse("Invalid offset field".into()))?;
    s.trim().parse::<u64>()
        .map_err(|_| CyaneaError::Parse(format!("Invalid offset: '{}'", s.trim())))
}

fn write_offset(buf: &mut [u8], val: u64) {
    let s = format!("{:>8}", val);
    buf.copy_from_slice(s.as_bytes());
}

fn parse_text_segment(text: &str) -> Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    if text.is_empty() {
        return Ok(map);
    }

    // First character is the delimiter
    let sep = text.chars().next().unwrap();
    let content = &text[sep.len_utf8()..];

    // Split by delimiter, pair up as key/value
    let parts: Vec<&str> = content.split(sep).collect();
    let mut i = 0;
    while i + 1 < parts.len() {
        let key = parts[i].trim();
        let val = parts[i + 1].trim();
        if !key.is_empty() {
            map.insert(key.to_uppercase(), val.to_string());
        }
        i += 2;
    }

    Ok(map)
}

fn parse_amplification(s: &str) -> (f64, f64) {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() == 2 {
        let a = parts[0].trim().parse().unwrap_or(0.0);
        let b = parts[1].trim().parse().unwrap_or(0.0);
        (a, b)
    } else {
        (0.0, 0.0)
    }
}

fn parse_data_segment(
    data: &[u8],
    datatype: &str,
    parameters: &[FcsParameter],
    n_events: usize,
    little_endian: bool,
) -> Result<Vec<Vec<f64>>> {
    let n_params = parameters.len();
    let mut events = Vec::with_capacity(n_events);

    match datatype {
        "F" => {
            // 32-bit float
            let bytes_per_event = n_params * 4;
            for i in 0..n_events {
                let offset = i * bytes_per_event;
                if offset + bytes_per_event > data.len() {
                    break;
                }
                let mut event = Vec::with_capacity(n_params);
                for j in 0..n_params {
                    let o = offset + j * 4;
                    let val = if little_endian {
                        f32::from_le_bytes([data[o], data[o + 1], data[o + 2], data[o + 3]])
                    } else {
                        f32::from_be_bytes([data[o], data[o + 1], data[o + 2], data[o + 3]])
                    };
                    event.push(val as f64);
                }
                events.push(event);
            }
        }
        "D" => {
            // 64-bit double
            let bytes_per_event = n_params * 8;
            for i in 0..n_events {
                let offset = i * bytes_per_event;
                if offset + bytes_per_event > data.len() {
                    break;
                }
                let mut event = Vec::with_capacity(n_params);
                for j in 0..n_params {
                    let o = offset + j * 8;
                    let val = if little_endian {
                        f64::from_le_bytes([
                            data[o], data[o + 1], data[o + 2], data[o + 3],
                            data[o + 4], data[o + 5], data[o + 6], data[o + 7],
                        ])
                    } else {
                        f64::from_be_bytes([
                            data[o], data[o + 1], data[o + 2], data[o + 3],
                            data[o + 4], data[o + 5], data[o + 6], data[o + 7],
                        ])
                    };
                    event.push(val);
                }
                events.push(event);
            }
        }
        "I" => {
            // Integer — use per-parameter bit widths
            let total_bits: u32 = parameters.iter().map(|p| p.bits).sum();
            let bytes_per_event = (total_bits / 8) as usize;
            for i in 0..n_events {
                let offset = i * bytes_per_event;
                if offset + bytes_per_event > data.len() {
                    break;
                }
                let mut event = Vec::with_capacity(n_params);
                let mut bit_offset = 0usize;
                for p in parameters {
                    let byte_off = offset + bit_offset / 8;
                    let val = match p.bits {
                        8 => data[byte_off] as f64,
                        16 => {
                            if little_endian {
                                u16::from_le_bytes([data[byte_off], data[byte_off + 1]]) as f64
                            } else {
                                u16::from_be_bytes([data[byte_off], data[byte_off + 1]]) as f64
                            }
                        }
                        32 => {
                            if little_endian {
                                u32::from_le_bytes([
                                    data[byte_off], data[byte_off + 1],
                                    data[byte_off + 2], data[byte_off + 3],
                                ]) as f64
                            } else {
                                u32::from_be_bytes([
                                    data[byte_off], data[byte_off + 1],
                                    data[byte_off + 2], data[byte_off + 3],
                                ]) as f64
                            }
                        }
                        _ => 0.0,
                    };
                    event.push(val);
                    bit_offset += p.bits as usize;
                }
                events.push(event);
            }
        }
        "A" => {
            // ASCII — not commonly used, basic support
            let text = std::str::from_utf8(data)
                .map_err(|_| CyaneaError::Parse("Invalid ASCII DATA segment".into()))?;
            for line in text.lines().take(n_events) {
                let event: Vec<f64> = line.split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if event.len() == n_params {
                    events.push(event);
                }
            }
        }
        _ => {
            return Err(CyaneaError::Parse(format!("Unknown $DATATYPE: '{}'", datatype)));
        }
    }

    Ok(events)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal FCS 3.1 binary blob for testing.
    fn make_test_fcs(n_events: usize, n_params: usize, values: &[f32]) -> Vec<u8> {
        // Build TEXT segment
        let sep = '/';
        let mut text = String::new();
        text.push(sep);

        let params: Vec<String> = (0..n_params).map(|i| format!("P{}", i + 1)).collect();

        let pairs = vec![
            ("$DATATYPE", "F"),
            ("$MODE", "L"),
            ("$BYTEORD", "1,2,3,4"),
        ];
        for (k, v) in &pairs {
            text.push_str(k);
            text.push(sep);
            text.push_str(v);
            text.push(sep);
        }
        text.push_str("$PAR");
        text.push(sep);
        text.push_str(&n_params.to_string());
        text.push(sep);
        text.push_str("$TOT");
        text.push(sep);
        text.push_str(&n_events.to_string());
        text.push(sep);

        for (i, name) in params.iter().enumerate() {
            text.push_str(&format!("$P{}N", i + 1));
            text.push(sep);
            text.push_str(name);
            text.push(sep);
            text.push_str(&format!("$P{}B", i + 1));
            text.push(sep);
            text.push_str("32");
            text.push(sep);
            text.push_str(&format!("$P{}R", i + 1));
            text.push(sep);
            text.push_str("262144");
            text.push(sep);
            text.push_str(&format!("$P{}E", i + 1));
            text.push(sep);
            text.push_str("0,0");
            text.push(sep);
        }

        let text_bytes = text.as_bytes();
        let text_start: usize = 256;
        let text_end = text_start + text_bytes.len() - 1;
        let data_start = text_end + 1;
        let data_size = values.len() * 4;
        let data_end = if data_size > 0 { data_start + data_size - 1 } else { data_start };

        // Add data offsets to text
        let mut full_text = text.clone();
        full_text.push_str("$BEGINDATA");
        full_text.push(sep);
        full_text.push_str(&data_start.to_string());
        full_text.push(sep);
        full_text.push_str("$ENDDATA");
        full_text.push(sep);
        full_text.push_str(&data_end.to_string());
        full_text.push(sep);

        let text_bytes = full_text.as_bytes();
        let text_end = text_start + text_bytes.len() - 1;
        let data_start = text_end + 1;
        let data_end = if data_size > 0 { data_start + data_size - 1 } else { data_start };

        let total_size = data_end + 1;
        let mut buf = vec![0u8; total_size.max(256)];

        // HEADER
        buf[0..6].copy_from_slice(b"FCS3.1");
        buf[6..10].copy_from_slice(b"    ");
        write_offset(&mut buf[10..18], text_start as u64);
        write_offset(&mut buf[18..26], text_end as u64);
        write_offset(&mut buf[26..34], data_start as u64);
        write_offset(&mut buf[34..42], data_end as u64);
        write_offset(&mut buf[42..50], 0);
        write_offset(&mut buf[50..58], 0);

        buf[text_start..text_start + text_bytes.len()].copy_from_slice(text_bytes);

        let mut offset = data_start;
        for &val in values {
            let bytes = val.to_le_bytes();
            buf[offset..offset + 4].copy_from_slice(&bytes);
            offset += 4;
        }

        buf
    }

    #[test]
    fn test_parse_basic_fcs() {
        let values: Vec<f32> = vec![
            100.0, 200.0, 300.0,  // event 1
            150.0, 250.0, 350.0,  // event 2
            120.0, 220.0, 320.0,  // event 3
        ];
        let data = make_test_fcs(3, 3, &values);
        let fcs = parse_fcs(&data).unwrap();

        assert_eq!(fcs.version, "FCS3.1");
        assert_eq!(fcs.n_parameters(), 3);
        assert_eq!(fcs.n_events(), 3);
        assert_eq!(fcs.parameters[0].name, "P1");
        assert_eq!(fcs.parameters[1].name, "P2");
        assert_eq!(fcs.parameters[2].name, "P3");
    }

    #[test]
    fn test_event_values() {
        let values: Vec<f32> = vec![
            10.0, 20.0,
            30.0, 40.0,
        ];
        let data = make_test_fcs(2, 2, &values);
        let fcs = parse_fcs(&data).unwrap();

        assert_eq!(fcs.events.len(), 2);
        assert!((fcs.events[0][0] - 10.0).abs() < 0.01);
        assert!((fcs.events[0][1] - 20.0).abs() < 0.01);
        assert!((fcs.events[1][0] - 30.0).abs() < 0.01);
        assert!((fcs.events[1][1] - 40.0).abs() < 0.01);
    }

    #[test]
    fn test_parameter_data() {
        let values: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let data = make_test_fcs(3, 2, &values);
        let fcs = parse_fcs(&data).unwrap();

        let col0 = fcs.parameter_data(0).unwrap();
        assert_eq!(col0.len(), 3);
        assert!((col0[0] - 1.0).abs() < 0.01);
        assert!((col0[1] - 3.0).abs() < 0.01);
        assert!((col0[2] - 5.0).abs() < 0.01);

        let col1 = fcs.parameter_data_by_name("P2").unwrap();
        assert!((col1[0] - 2.0).abs() < 0.01);
        assert!((col1[1] - 4.0).abs() < 0.01);
    }

    #[test]
    fn test_keyword_access() {
        let data = make_test_fcs(1, 2, &[1.0, 2.0]);
        let fcs = parse_fcs(&data).unwrap();

        assert_eq!(fcs.keyword("$PAR"), Some("2"));
        assert_eq!(fcs.keyword("$TOT"), Some("1"));
        assert_eq!(fcs.keyword("$DATATYPE"), Some("F"));
        assert!(fcs.keyword("$NONEXISTENT").is_none());
    }

    #[test]
    fn test_fcs_stats() {
        let data = make_test_fcs(5, 4, &vec![0.0; 20]);
        let fcs = parse_fcs(&data).unwrap();
        let stats = fcs.stats();

        assert_eq!(stats.n_events, 5);
        assert_eq!(stats.n_parameters, 4);
        assert_eq!(stats.parameter_names, vec!["P1", "P2", "P3", "P4"]);
    }

    #[test]
    fn test_write_roundtrip() {
        let values: Vec<f32> = vec![100.0, 200.0, 300.0, 400.0, 500.0, 600.0];
        let data = make_test_fcs(3, 2, &values);
        let fcs = parse_fcs(&data).unwrap();

        let written = write_fcs(&fcs).unwrap();
        let fcs2 = parse_fcs(&written).unwrap();

        assert_eq!(fcs2.n_events(), 3);
        assert_eq!(fcs2.n_parameters(), 2);
        for i in 0..3 {
            for j in 0..2 {
                assert!((fcs2.events[i][j] - fcs.events[i][j]).abs() < 0.1);
            }
        }
    }

    #[test]
    fn test_empty_events() {
        let data = make_test_fcs(0, 3, &[]);
        let fcs = parse_fcs(&data).unwrap();
        assert_eq!(fcs.n_events(), 0);
        assert_eq!(fcs.n_parameters(), 3);
    }

    #[test]
    fn test_invalid_magic() {
        let mut data = make_test_fcs(1, 1, &[1.0]);
        data[0..6].copy_from_slice(b"NOTFCS");
        assert!(parse_fcs(&data).is_err());
    }

    #[test]
    fn test_too_short() {
        assert!(parse_fcs(&[0u8; 10]).is_err());
    }

    #[test]
    fn test_parse_text_segment() {
        let text = "/KEY1/val1/KEY2/val2/";
        let map = parse_text_segment(text).unwrap();
        assert_eq!(map.get("KEY1").unwrap(), "val1");
        assert_eq!(map.get("KEY2").unwrap(), "val2");
    }

    #[test]
    fn test_parameter_bits_and_range() {
        let data = make_test_fcs(1, 2, &[1.0, 2.0]);
        let fcs = parse_fcs(&data).unwrap();

        assert_eq!(fcs.parameters[0].bits, 32);
        assert!((fcs.parameters[0].range - 262144.0).abs() < 0.1);
        assert_eq!(fcs.parameters[0].amplification, (0.0, 0.0));
    }

    #[test]
    fn test_stain_field() {
        // The basic test FCS doesn't include stains, so they should be None
        let data = make_test_fcs(1, 2, &[1.0, 2.0]);
        let fcs = parse_fcs(&data).unwrap();
        assert!(fcs.parameters[0].stain.is_none());
    }

    #[test]
    fn test_many_events() {
        let n = 1000;
        let values: Vec<f32> = (0..n * 3).map(|i| i as f32).collect();
        let data = make_test_fcs(n, 3, &values);
        let fcs = parse_fcs(&data).unwrap();

        assert_eq!(fcs.n_events(), n);
        assert!((fcs.events[999][2] - 2999.0).abs() < 0.01);
    }
}
