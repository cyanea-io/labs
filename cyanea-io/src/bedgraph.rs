//! bedGraph and Wiggle format parser/writer.
//!
//! bedGraph is a BED-like format for continuous-valued data (e.g. signal
//! tracks, coverage). Wiggle (WIG) is a related format with two sub-formats:
//! `variableStep` and `fixedStep`, both of which are converted to
//! [`BedGraphRecord`]s on parsing.

use cyanea_core::{CyaneaError, Result};

/// A single bedGraph record: a genomic interval with an associated value.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BedGraphRecord {
    /// Chromosome name.
    pub chrom: String,
    /// 0-based start coordinate (inclusive).
    pub start: u64,
    /// 0-based end coordinate (exclusive).
    pub end: u64,
    /// Data value for this interval.
    pub value: f64,
}

/// Parse bedGraph format from a string.
///
/// Skips lines starting with `track`, `browser`, or `#`. Each data line must
/// have exactly 4 whitespace-separated fields: `chrom start end value`.
///
/// # Errors
///
/// Returns an error if a data line is malformed (wrong field count or
/// unparseable coordinates/values), with line number context.
///
/// # Examples
///
/// ```
/// # use cyanea_io::bedgraph::parse_bedgraph_str;
/// let data = "chr1\t0\t100\t1.5\nchr1\t100\t200\t2.3\n";
/// let records = parse_bedgraph_str(data).unwrap();
/// assert_eq!(records.len(), 2);
/// assert!((records[0].value - 1.5).abs() < f64::EPSILON);
/// ```
pub fn parse_bedgraph_str(data: &str) -> Result<Vec<BedGraphRecord>> {
    let mut records = Vec::new();

    for (line_idx, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with("track")
            || line.starts_with("browser")
            || line.starts_with('#')
        {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(CyaneaError::Parse(format!(
                "line {}: expected 4 fields, found {}",
                line_idx + 1,
                fields.len()
            )));
        }

        let chrom = fields[0].to_string();

        let start: u64 = fields[1].parse().map_err(|_| {
            CyaneaError::Parse(format!(
                "line {}: invalid start coordinate '{}'",
                line_idx + 1,
                fields[1]
            ))
        })?;

        let end: u64 = fields[2].parse().map_err(|_| {
            CyaneaError::Parse(format!(
                "line {}: invalid end coordinate '{}'",
                line_idx + 1,
                fields[2]
            ))
        })?;

        let value: f64 = fields[3].parse().map_err(|_| {
            CyaneaError::Parse(format!(
                "line {}: invalid value '{}'",
                line_idx + 1,
                fields[3]
            ))
        })?;

        records.push(BedGraphRecord {
            chrom,
            start,
            end,
            value,
        });
    }

    Ok(records)
}

/// Write bedGraph records to a string.
///
/// If `track_name` is `Some`, a track header line is emitted first.
/// Each record is written as `{chrom}\t{start}\t{end}\t{value}`.
///
/// # Examples
///
/// ```
/// # use cyanea_io::bedgraph::{BedGraphRecord, write_bedgraph_string};
/// let records = vec![
///     BedGraphRecord { chrom: "chr1".into(), start: 0, end: 100, value: 1.5 },
/// ];
/// let out = write_bedgraph_string(&records, Some("coverage"));
/// assert!(out.starts_with("track type=bedGraph"));
/// ```
pub fn write_bedgraph_string(records: &[BedGraphRecord], track_name: Option<&str>) -> String {
    let mut out = String::new();

    if let Some(name) = track_name {
        out.push_str(&format!("track type=bedGraph name=\"{}\"\n", name));
    }

    for rec in records {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            rec.chrom, rec.start, rec.end, rec.value
        ));
    }

    out
}

/// Parse a key=value pair from a wiggle declaration header.
fn parse_wig_kv<'a>(fields: &[&'a str], key: &str) -> Option<&'a str> {
    for field in fields {
        if let Some(val) = field.strip_prefix(&format!("{}=", key)) {
            return Some(val);
        }
    }
    None
}

/// Parse Wiggle (WIG) format from a string, converting to bedGraph records.
///
/// Supports both `variableStep` and `fixedStep` sub-formats. Wiggle uses
/// 1-based coordinates; these are converted to 0-based half-open intervals
/// in the returned [`BedGraphRecord`]s.
///
/// # variableStep
///
/// ```text
/// variableStep chrom=chr1 span=10
/// 100  1.5
/// 200  2.3
/// ```
///
/// Each data line contains `position value`. The interval is
/// `[position-1, position-1+span)` (default span=1).
///
/// # fixedStep
///
/// ```text
/// fixedStep chrom=chr1 start=100 step=10 span=5
/// 1.5
/// 2.3
/// ```
///
/// Data lines contain only `value`. Positions are computed from `start`
/// (1-based), incrementing by `step` each line.
///
/// # Errors
///
/// Returns an error for missing required fields in declaration lines or
/// malformed data lines.
pub fn parse_wiggle_str(data: &str) -> Result<Vec<BedGraphRecord>> {
    let mut records = Vec::new();

    // Current section state
    enum WigSection {
        None,
        VariableStep { chrom: String, span: u64 },
        FixedStep { chrom: String, _start: u64, step: u64, span: u64, current_pos: u64 },
    }

    let mut section = WigSection::None;

    for (line_idx, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with("track")
            || line.starts_with("browser")
            || line.starts_with('#')
        {
            continue;
        }

        if line.starts_with("variableStep") {
            let fields: Vec<&str> = line.split_whitespace().collect();
            let chrom = parse_wig_kv(&fields, "chrom").ok_or_else(|| {
                CyaneaError::Parse(format!(
                    "line {}: variableStep missing chrom",
                    line_idx + 1
                ))
            })?;
            let span: u64 = parse_wig_kv(&fields, "span")
                .map(|s| s.parse::<u64>())
                .transpose()
                .map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid span value",
                        line_idx + 1
                    ))
                })?
                .unwrap_or(1);

            section = WigSection::VariableStep {
                chrom: chrom.to_string(),
                span,
            };
            continue;
        }

        if line.starts_with("fixedStep") {
            let fields: Vec<&str> = line.split_whitespace().collect();
            let chrom = parse_wig_kv(&fields, "chrom").ok_or_else(|| {
                CyaneaError::Parse(format!(
                    "line {}: fixedStep missing chrom",
                    line_idx + 1
                ))
            })?;
            let start: u64 = parse_wig_kv(&fields, "start")
                .ok_or_else(|| {
                    CyaneaError::Parse(format!(
                        "line {}: fixedStep missing start",
                        line_idx + 1
                    ))
                })?
                .parse()
                .map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid start value",
                        line_idx + 1
                    ))
                })?;
            let step: u64 = parse_wig_kv(&fields, "step")
                .ok_or_else(|| {
                    CyaneaError::Parse(format!(
                        "line {}: fixedStep missing step",
                        line_idx + 1
                    ))
                })?
                .parse()
                .map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid step value",
                        line_idx + 1
                    ))
                })?;
            let span: u64 = parse_wig_kv(&fields, "span")
                .map(|s| s.parse::<u64>())
                .transpose()
                .map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid span value",
                        line_idx + 1
                    ))
                })?
                .unwrap_or(1);

            section = WigSection::FixedStep {
                chrom: chrom.to_string(),
                _start: start,
                step,
                span,
                current_pos: start,
            };
            continue;
        }

        // Data line — parse according to current section.
        match &mut section {
            WigSection::None => {
                return Err(CyaneaError::Parse(format!(
                    "line {}: data line before any declaration header",
                    line_idx + 1
                )));
            }
            WigSection::VariableStep { chrom, span } => {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() < 2 {
                    return Err(CyaneaError::Parse(format!(
                        "line {}: variableStep data line expected 2 fields, found {}",
                        line_idx + 1,
                        fields.len()
                    )));
                }
                let position: u64 = fields[0].parse().map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid position '{}'",
                        line_idx + 1,
                        fields[0]
                    ))
                })?;
                let value: f64 = fields[1].parse().map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid value '{}'",
                        line_idx + 1,
                        fields[1]
                    ))
                })?;

                // Wiggle is 1-based; convert to 0-based half-open.
                let start = position - 1;
                records.push(BedGraphRecord {
                    chrom: chrom.clone(),
                    start,
                    end: start + *span,
                    value,
                });
            }
            WigSection::FixedStep {
                chrom,
                step,
                span,
                current_pos,
                ..
            } => {
                let value: f64 = line.trim().parse().map_err(|_| {
                    CyaneaError::Parse(format!(
                        "line {}: invalid value '{}'",
                        line_idx + 1,
                        line.trim()
                    ))
                })?;

                // Wiggle is 1-based; convert to 0-based half-open.
                let start = *current_pos - 1;
                records.push(BedGraphRecord {
                    chrom: chrom.clone(),
                    start,
                    end: start + *span,
                    value,
                });
                *current_pos += *step;
            }
        }
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_bedgraph_with_track_header() {
        let data = "\
track type=bedGraph name=\"coverage\"
chr1\t0\t100\t1.5
chr1\t100\t200\t2.3
chr2\t0\t50\t0.8
";
        let records = parse_bedgraph_str(data).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].start, 0);
        assert_eq!(records[0].end, 100);
        assert!((records[0].value - 1.5).abs() < f64::EPSILON);
        assert_eq!(records[2].chrom, "chr2");
        assert!((records[2].value - 0.8).abs() < f64::EPSILON);
    }

    #[test]
    fn bedgraph_roundtrip() {
        let original = vec![
            BedGraphRecord { chrom: "chr1".into(), start: 0, end: 100, value: 1.5 },
            BedGraphRecord { chrom: "chr1".into(), start: 100, end: 200, value: 2.0 },
            BedGraphRecord { chrom: "chr2".into(), start: 50, end: 150, value: 3.7 },
        ];

        let written = write_bedgraph_string(&original, Some("test"));
        let parsed = parse_bedgraph_str(&written).unwrap();

        assert_eq!(parsed.len(), original.len());
        for (a, b) in original.iter().zip(parsed.iter()) {
            assert_eq!(a.chrom, b.chrom);
            assert_eq!(a.start, b.start);
            assert_eq!(a.end, b.end);
            assert!((a.value - b.value).abs() < f64::EPSILON);
        }
    }

    #[test]
    fn parse_wiggle_variable_step() {
        let data = "\
variableStep chrom=chr1 span=10
100\t1.5
200\t2.3
300\t3.1
";
        let records = parse_wiggle_str(data).unwrap();
        assert_eq!(records.len(), 3);

        // Position 100 (1-based) -> start=99, end=99+10=109
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].start, 99);
        assert_eq!(records[0].end, 109);
        assert!((records[0].value - 1.5).abs() < f64::EPSILON);

        assert_eq!(records[1].start, 199);
        assert_eq!(records[1].end, 209);

        assert_eq!(records[2].start, 299);
        assert_eq!(records[2].end, 309);
    }

    #[test]
    fn parse_wiggle_fixed_step() {
        let data = "\
fixedStep chrom=chr2 start=100 step=50 span=25
1.0
2.0
3.0
";
        let records = parse_wiggle_str(data).unwrap();
        assert_eq!(records.len(), 3);

        // start=100 (1-based) -> 0-based 99, span=25
        assert_eq!(records[0].chrom, "chr2");
        assert_eq!(records[0].start, 99);
        assert_eq!(records[0].end, 124);
        assert!((records[0].value - 1.0).abs() < f64::EPSILON);

        // next: 100+50=150 (1-based) -> 0-based 149
        assert_eq!(records[1].start, 149);
        assert_eq!(records[1].end, 174);

        // next: 150+50=200 (1-based) -> 0-based 199
        assert_eq!(records[2].start, 199);
        assert_eq!(records[2].end, 224);
    }

    #[test]
    fn malformed_bedgraph_returns_error() {
        // Missing value field
        let data = "chr1\t0\t100\n";
        let result = parse_bedgraph_str(data);
        assert!(result.is_err());

        // Non-numeric start
        let data2 = "chr1\tabc\t100\t1.5\n";
        let result2 = parse_bedgraph_str(data2);
        assert!(result2.is_err());
    }
}
