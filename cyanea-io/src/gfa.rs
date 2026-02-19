//! GFA v1 (Graphical Fragment Assembly) format parser/writer.
//!
//! GFA is a tab-separated format for representing sequence graphs (pangenomes,
//! assembly graphs). This module supports the three main record types of GFA v1:
//!
//! - **S** (Segment) — a named sequence
//! - **L** (Link) — an edge between two oriented segments
//! - **P** (Path) — a walk through the graph
//!
//! Optional header lines (`H`) and comments (`#`) are also handled.

use cyanea_core::{CyaneaError, Result};

/// A GFA segment: a named sequence in the graph.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GfaSegment {
    /// Segment name (unique identifier).
    pub name: String,
    /// Nucleotide sequence. Empty if the sequence was `*` (unknown).
    pub sequence: Vec<u8>,
    /// Sequence length. Taken from `LN:i:` tag if sequence is `*`, otherwise
    /// from the sequence itself.
    pub length: Option<u64>,
}

/// A GFA link: a directed edge between two oriented segments.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GfaLink {
    /// Name of the source segment.
    pub from_segment: String,
    /// Orientation of the source segment (`+` or `-`).
    pub from_orient: char,
    /// Name of the target segment.
    pub to_segment: String,
    /// Orientation of the target segment (`+` or `-`).
    pub to_orient: char,
    /// Overlap CIGAR string.
    pub overlap: String,
}

/// A GFA path: an ordered walk through the graph.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GfaPath {
    /// Path name.
    pub name: String,
    /// Ordered list of oriented segment names (e.g. `["seg1+", "seg2-"]`).
    pub segments: Vec<String>,
    /// Overlap CIGARs between consecutive segments, or `["*"]` if absent.
    pub overlaps: Vec<String>,
}

/// A complete GFA graph.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GfaGraph {
    /// All segments in the graph.
    pub segments: Vec<GfaSegment>,
    /// All links (edges) in the graph.
    pub links: Vec<GfaLink>,
    /// All paths through the graph.
    pub paths: Vec<GfaPath>,
    /// Optional header string (contents after `H\t`).
    pub header: Option<String>,
}

/// Summary statistics for a GFA graph.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GfaStats {
    /// Number of segments.
    pub n_segments: usize,
    /// Number of links.
    pub n_links: usize,
    /// Number of paths.
    pub n_paths: usize,
    /// Total sequence length across all segments.
    pub total_sequence_length: u64,
    /// Mean segment length (0.0 if no segments).
    pub mean_segment_length: f64,
}

/// Parse a GFA v1 graph from a string.
///
/// Each line is identified by its first tab-separated field:
/// - `H` — header
/// - `S` — segment
/// - `L` — link
/// - `P` — path
/// - `#` — comment (skipped)
///
/// # Errors
///
/// Returns an error if a record line has too few fields or contains
/// invalid data, with line number context.
///
/// # Examples
///
/// ```
/// # use cyanea_io::gfa::parse_gfa_str;
/// let data = "H\tVN:Z:1.0\nS\ts1\tACGT\nS\ts2\tTGCA\nL\ts1\t+\ts2\t+\t4M\n";
/// let graph = parse_gfa_str(data).unwrap();
/// assert_eq!(graph.segments.len(), 2);
/// assert_eq!(graph.links.len(), 1);
/// ```
pub fn parse_gfa_str(data: &str) -> Result<GfaGraph> {
    let mut segments = Vec::new();
    let mut links = Vec::new();
    let mut paths = Vec::new();
    let mut header: Option<String> = None;

    for (line_idx, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }

        match fields[0] {
            "H" => {
                // Header line: store everything after H\t
                if fields.len() > 1 {
                    header = Some(fields[1..].join("\t"));
                }
            }
            "S" => {
                // Segment: S\tname\tsequence[\tLN:i:123]
                if fields.len() < 3 {
                    return Err(CyaneaError::Parse(format!(
                        "line {}: segment (S) requires at least 3 fields, found {}",
                        line_idx + 1,
                        fields.len()
                    )));
                }
                let name = fields[1].to_string();
                let seq_str = fields[2];

                let (sequence, length) = if seq_str == "*" {
                    // Unknown sequence — look for LN:i: tag in optional fields.
                    let ln = fields[3..]
                        .iter()
                        .find_map(|f| f.strip_prefix("LN:i:"))
                        .map(|v| {
                            v.parse::<u64>().map_err(|_| {
                                CyaneaError::Parse(format!(
                                    "line {}: invalid LN tag value '{}'",
                                    line_idx + 1,
                                    v
                                ))
                            })
                        })
                        .transpose()?;
                    (Vec::new(), ln)
                } else {
                    let seq = seq_str.as_bytes().to_vec();
                    let len = seq.len() as u64;
                    (seq, Some(len))
                };

                segments.push(GfaSegment {
                    name,
                    sequence,
                    length,
                });
            }
            "L" => {
                // Link: L\tfrom\tfrom_orient\tto\tto_orient\toverlap
                if fields.len() < 6 {
                    return Err(CyaneaError::Parse(format!(
                        "line {}: link (L) requires at least 6 fields, found {}",
                        line_idx + 1,
                        fields.len()
                    )));
                }
                let from_orient = fields[2].chars().next().ok_or_else(|| {
                    CyaneaError::Parse(format!(
                        "line {}: empty from_orient field",
                        line_idx + 1
                    ))
                })?;
                let to_orient = fields[4].chars().next().ok_or_else(|| {
                    CyaneaError::Parse(format!(
                        "line {}: empty to_orient field",
                        line_idx + 1
                    ))
                })?;

                links.push(GfaLink {
                    from_segment: fields[1].to_string(),
                    from_orient,
                    to_segment: fields[3].to_string(),
                    to_orient,
                    overlap: fields[5].to_string(),
                });
            }
            "P" => {
                // Path: P\tname\tsegments\toverlaps
                if fields.len() < 4 {
                    return Err(CyaneaError::Parse(format!(
                        "line {}: path (P) requires at least 4 fields, found {}",
                        line_idx + 1,
                        fields.len()
                    )));
                }
                let name = fields[1].to_string();
                let segs: Vec<String> = fields[2]
                    .split(',')
                    .map(|s| s.to_string())
                    .collect();
                let overlaps: Vec<String> = fields[3]
                    .split(',')
                    .map(|s| s.to_string())
                    .collect();

                paths.push(GfaPath {
                    name,
                    segments: segs,
                    overlaps,
                });
            }
            _ => {
                // Unknown record type — skip silently (GFA allows extensions).
            }
        }
    }

    Ok(GfaGraph {
        segments,
        links,
        paths,
        header,
    })
}

/// Write a GFA graph to a string in GFA v1 format.
///
/// Output order: header, segments, links, paths. Segments with an empty
/// sequence are written as `*`.
///
/// # Examples
///
/// ```
/// # use cyanea_io::gfa::{GfaGraph, GfaSegment, write_gfa_string};
/// let graph = GfaGraph {
///     segments: vec![GfaSegment { name: "s1".into(), sequence: b"ACGT".to_vec(), length: Some(4) }],
///     links: vec![],
///     paths: vec![],
///     header: Some("VN:Z:1.0".into()),
/// };
/// let out = write_gfa_string(&graph);
/// assert!(out.starts_with("H\tVN:Z:1.0"));
/// ```
pub fn write_gfa_string(graph: &GfaGraph) -> String {
    let mut out = String::new();

    // Header
    if let Some(ref hdr) = graph.header {
        out.push_str(&format!("H\t{}\n", hdr));
    }

    // Segments
    for seg in &graph.segments {
        if seg.sequence.is_empty() {
            out.push_str(&format!("S\t{}\t*", seg.name));
            if let Some(ln) = seg.length {
                out.push_str(&format!("\tLN:i:{}", ln));
            }
            out.push('\n');
        } else {
            let seq_str = std::str::from_utf8(&seg.sequence).unwrap_or("*");
            out.push_str(&format!("S\t{}\t{}\n", seg.name, seq_str));
        }
    }

    // Links
    for link in &graph.links {
        out.push_str(&format!(
            "L\t{}\t{}\t{}\t{}\t{}\n",
            link.from_segment,
            link.from_orient,
            link.to_segment,
            link.to_orient,
            link.overlap
        ));
    }

    // Paths
    for path in &graph.paths {
        let segs = path.segments.join(",");
        let ovls = path.overlaps.join(",");
        out.push_str(&format!("P\t{}\t{}\t{}\n", path.name, segs, ovls));
    }

    out
}

/// Compute summary statistics for a GFA graph.
///
/// # Examples
///
/// ```
/// # use cyanea_io::gfa::{GfaGraph, GfaSegment, gfa_stats};
/// let graph = GfaGraph {
///     segments: vec![
///         GfaSegment { name: "s1".into(), sequence: b"ACGT".to_vec(), length: Some(4) },
///         GfaSegment { name: "s2".into(), sequence: b"TG".to_vec(), length: Some(2) },
///     ],
///     links: vec![],
///     paths: vec![],
///     header: None,
/// };
/// let stats = gfa_stats(&graph);
/// assert_eq!(stats.total_sequence_length, 6);
/// assert!((stats.mean_segment_length - 3.0).abs() < f64::EPSILON);
/// ```
pub fn gfa_stats(graph: &GfaGraph) -> GfaStats {
    let n_segments = graph.segments.len();
    let n_links = graph.links.len();
    let n_paths = graph.paths.len();

    let total_sequence_length: u64 = graph
        .segments
        .iter()
        .map(|seg| seg.length.unwrap_or(seg.sequence.len() as u64))
        .sum();

    let mean_segment_length = if n_segments > 0 {
        total_sequence_length as f64 / n_segments as f64
    } else {
        0.0
    };

    GfaStats {
        n_segments,
        n_links,
        n_paths,
        total_sequence_length,
        mean_segment_length,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_segments_and_links() {
        let data = "\
H\tVN:Z:1.0
S\ts1\tACGT
S\ts2\tTGCA
S\ts3\t*\tLN:i:100
L\ts1\t+\ts2\t+\t4M
L\ts2\t+\ts3\t-\t0M
";
        let graph = parse_gfa_str(data).unwrap();

        assert_eq!(graph.header, Some("VN:Z:1.0".to_string()));
        assert_eq!(graph.segments.len(), 3);
        assert_eq!(graph.links.len(), 2);

        // Segment s1
        assert_eq!(graph.segments[0].name, "s1");
        assert_eq!(graph.segments[0].sequence, b"ACGT");
        assert_eq!(graph.segments[0].length, Some(4));

        // Segment s3 (unknown sequence with LN tag)
        assert_eq!(graph.segments[2].name, "s3");
        assert!(graph.segments[2].sequence.is_empty());
        assert_eq!(graph.segments[2].length, Some(100));

        // Link s1+ -> s2+
        assert_eq!(graph.links[0].from_segment, "s1");
        assert_eq!(graph.links[0].from_orient, '+');
        assert_eq!(graph.links[0].to_segment, "s2");
        assert_eq!(graph.links[0].to_orient, '+');
        assert_eq!(graph.links[0].overlap, "4M");

        // Link s2+ -> s3-
        assert_eq!(graph.links[1].to_orient, '-');
    }

    #[test]
    fn parse_paths_with_orientations() {
        let data = "\
S\ts1\tACGT
S\ts2\tTGCA
S\ts3\tAAAA
P\tpath1\ts1+,s2+,s3-\t4M,0M
P\tpath2\ts3+,s1-\t*
";
        let graph = parse_gfa_str(data).unwrap();

        assert_eq!(graph.paths.len(), 2);

        assert_eq!(graph.paths[0].name, "path1");
        assert_eq!(
            graph.paths[0].segments,
            vec!["s1+", "s2+", "s3-"]
        );
        assert_eq!(graph.paths[0].overlaps, vec!["4M", "0M"]);

        assert_eq!(graph.paths[1].name, "path2");
        assert_eq!(graph.paths[1].segments, vec!["s3+", "s1-"]);
        assert_eq!(graph.paths[1].overlaps, vec!["*"]);
    }

    #[test]
    fn gfa_stats_correct() {
        let data = "\
S\ts1\tACGTACGT
S\ts2\tTGCA
S\ts3\t*\tLN:i:200
L\ts1\t+\ts2\t+\t4M
L\ts2\t+\ts3\t-\t0M
P\tref\ts1+,s2+,s3+\t4M,0M
";
        let graph = parse_gfa_str(data).unwrap();
        let stats = gfa_stats(&graph);

        assert_eq!(stats.n_segments, 3);
        assert_eq!(stats.n_links, 2);
        assert_eq!(stats.n_paths, 1);
        // 8 + 4 + 200 = 212
        assert_eq!(stats.total_sequence_length, 212);
        assert!((stats.mean_segment_length - 212.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn gfa_roundtrip() {
        let original = GfaGraph {
            header: Some("VN:Z:1.0".to_string()),
            segments: vec![
                GfaSegment {
                    name: "s1".to_string(),
                    sequence: b"ACGT".to_vec(),
                    length: Some(4),
                },
                GfaSegment {
                    name: "s2".to_string(),
                    sequence: b"TGCA".to_vec(),
                    length: Some(4),
                },
                GfaSegment {
                    name: "s3".to_string(),
                    sequence: Vec::new(),
                    length: Some(100),
                },
            ],
            links: vec![GfaLink {
                from_segment: "s1".to_string(),
                from_orient: '+',
                to_segment: "s2".to_string(),
                to_orient: '+',
                overlap: "4M".to_string(),
            }],
            paths: vec![GfaPath {
                name: "ref".to_string(),
                segments: vec!["s1+".to_string(), "s2+".to_string()],
                overlaps: vec!["4M".to_string()],
            }],
        };

        let written = write_gfa_string(&original);
        let parsed = parse_gfa_str(&written).unwrap();

        assert_eq!(parsed.header, original.header);
        assert_eq!(parsed.segments.len(), original.segments.len());
        assert_eq!(parsed.links.len(), original.links.len());
        assert_eq!(parsed.paths.len(), original.paths.len());

        for (a, b) in original.segments.iter().zip(parsed.segments.iter()) {
            assert_eq!(a.name, b.name);
            assert_eq!(a.sequence, b.sequence);
        }

        assert_eq!(parsed.links[0].from_segment, "s1");
        assert_eq!(parsed.links[0].overlap, "4M");
        assert_eq!(parsed.paths[0].name, "ref");
        assert_eq!(parsed.paths[0].segments, vec!["s1+", "s2+"]);
    }

    #[test]
    fn malformed_gfa_returns_error() {
        // Segment with too few fields
        let data = "S\ts1\n";
        let result = parse_gfa_str(data);
        assert!(result.is_err());

        // Link with too few fields
        let data2 = "L\ts1\t+\ts2\n";
        let result2 = parse_gfa_str(data2);
        assert!(result2.is_err());

        // Path with too few fields
        let data3 = "P\tpath1\ts1+,s2+\n";
        let result3 = parse_gfa_str(data3);
        assert!(result3.is_err());
    }
}
