//! Stockholm format parser and writer.
//!
//! Parses the Stockholm multiple sequence alignment format (version 1.0),
//! used by Pfam, Rfam, and HMMER. A single file can contain multiple
//! alignment blocks, each beginning with `# STOCKHOLM 1.0` and ending
//! with `//`.
//!
//! Supports all annotation types:
//! - `#=GF` — per-file (ignored in struct, kept as comments)
//! - `#=GC` — per-column (e.g. consensus structure)
//! - `#=GS` — per-sequence metadata
//! - `#=GR` — per-sequence per-column (e.g. posterior probability)

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

/// A single Stockholm alignment block.
///
/// Each block in a Stockholm file starts with `# STOCKHOLM 1.0` and
/// ends with `//`. Sequences may be interleaved (same name appears
/// multiple times), in which case their fragments are concatenated.
#[derive(Debug, Clone)]
pub struct StockholmAlignment {
    /// Aligned sequences as `(name, aligned_sequence)` pairs.
    pub sequences: Vec<(String, String)>,
    /// Per-column annotations from `#=GC` lines.
    pub gc_annotations: HashMap<String, String>,
    /// Per-sequence annotations from `#=GS` lines: `name -> (tag -> value)`.
    pub gs_annotations: HashMap<String, HashMap<String, String>>,
    /// Per-sequence per-column annotations from `#=GR` lines: `name -> (tag -> value)`.
    pub gr_annotations: HashMap<String, HashMap<String, String>>,
}

/// Parse one or more Stockholm alignment blocks from a string.
///
/// Returns an empty `Vec` if the input is empty or contains no blocks.
/// Returns an error if a block header is found but the block is malformed
/// (e.g. missing `//` terminator).
///
/// # Examples
///
/// ```
/// # use cyanea_io::stockholm::parse_stockholm;
/// let input = "# STOCKHOLM 1.0\nseq1 ACGT\nseq2 ACGT\n//\n";
/// let alns = parse_stockholm(input).unwrap();
/// assert_eq!(alns.len(), 1);
/// assert_eq!(alns[0].sequences.len(), 2);
/// ```
pub fn parse_stockholm(input: &str) -> Result<Vec<StockholmAlignment>> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    let mut alignments = Vec::new();
    let mut in_block = false;
    let mut seq_order: Vec<String> = Vec::new();
    let mut seq_data: HashMap<String, String> = HashMap::new();
    let mut gc_annotations: HashMap<String, String> = HashMap::new();
    let mut gs_annotations: HashMap<String, HashMap<String, String>> = HashMap::new();
    let mut gr_annotations: HashMap<String, HashMap<String, String>> = HashMap::new();
    let mut block_started = false;

    for (line_num, line) in input.lines().enumerate() {
        let line = line.trim_end();

        // Block header
        if line.starts_with("# STOCKHOLM 1.0") {
            if in_block {
                return Err(CyaneaError::Parse(format!(
                    "line {}: nested STOCKHOLM header without preceding //",
                    line_num + 1
                )));
            }
            in_block = true;
            block_started = true;
            seq_order.clear();
            seq_data.clear();
            gc_annotations.clear();
            gs_annotations.clear();
            gr_annotations.clear();
            continue;
        }

        if !in_block {
            continue;
        }

        // Block terminator
        if line == "//" {
            let sequences: Vec<(String, String)> = seq_order
                .iter()
                .map(|name| {
                    let seq = seq_data.get(name).cloned().unwrap_or_default();
                    (name.clone(), seq)
                })
                .collect();

            alignments.push(StockholmAlignment {
                sequences,
                gc_annotations: gc_annotations.clone(),
                gs_annotations: gs_annotations.clone(),
                gr_annotations: gr_annotations.clone(),
            });
            in_block = false;
            continue;
        }

        // Skip blank lines
        if line.trim().is_empty() {
            continue;
        }

        // #=GC tag value
        if let Some(rest) = line.strip_prefix("#=GC ") {
            let rest = rest.trim_start();
            if let Some((tag, value)) = rest.split_once(char::is_whitespace) {
                let tag = tag.trim();
                let value = value.trim();
                gc_annotations
                    .entry(tag.to_string())
                    .and_modify(|v| v.push_str(value))
                    .or_insert_with(|| value.to_string());
            }
            continue;
        }

        // #=GS seqname tag value
        if let Some(rest) = line.strip_prefix("#=GS ") {
            let rest = rest.trim_start();
            let parts: Vec<&str> = rest.splitn(3, char::is_whitespace).collect();
            if parts.len() >= 3 {
                let name = parts[0].trim();
                let tag = parts[1].trim();
                let value = parts[2].trim();
                gs_annotations
                    .entry(name.to_string())
                    .or_default()
                    .insert(tag.to_string(), value.to_string());
            }
            continue;
        }

        // #=GR seqname tag value
        if let Some(rest) = line.strip_prefix("#=GR ") {
            let rest = rest.trim_start();
            let parts: Vec<&str> = rest.splitn(3, char::is_whitespace).collect();
            if parts.len() >= 3 {
                let name = parts[0].trim();
                let tag = parts[1].trim();
                let value = parts[2].trim();
                gr_annotations
                    .entry(name.to_string())
                    .or_default()
                    .entry(tag.to_string())
                    .and_modify(|v| v.push_str(value))
                    .or_insert_with(|| value.to_string());
            }
            continue;
        }

        // #=GF lines and other comments — skip
        if line.starts_with('#') {
            continue;
        }

        // Sequence line: name sequence
        let parts: Vec<&str> = line.splitn(2, char::is_whitespace).collect();
        if parts.len() == 2 {
            let name = parts[0].trim();
            let seq = parts[1].trim();
            if !seq_data.contains_key(name) {
                seq_order.push(name.to_string());
            }
            seq_data
                .entry(name.to_string())
                .and_modify(|v| v.push_str(seq))
                .or_insert_with(|| seq.to_string());
        }
    }

    if in_block && block_started {
        return Err(CyaneaError::Parse(
            "unterminated Stockholm block (missing //)".to_string(),
        ));
    }

    Ok(alignments)
}

/// Write a Stockholm alignment block as a string.
///
/// Produces a single block with `# STOCKHOLM 1.0` header and `//` terminator.
///
/// # Examples
///
/// ```
/// # use cyanea_io::stockholm::{StockholmAlignment, write_stockholm};
/// # use std::collections::HashMap;
/// let aln = StockholmAlignment {
///     sequences: vec![("seq1".to_string(), "ACGT".to_string())],
///     gc_annotations: HashMap::new(),
///     gs_annotations: HashMap::new(),
///     gr_annotations: HashMap::new(),
/// };
/// let output = write_stockholm(&aln);
/// assert!(output.starts_with("# STOCKHOLM 1.0"));
/// ```
pub fn write_stockholm(aln: &StockholmAlignment) -> String {
    let mut out = String::new();
    out.push_str("# STOCKHOLM 1.0\n");

    // #=GS annotations
    for (name, tags) in &aln.gs_annotations {
        for (tag, value) in tags {
            out.push_str(&format!("#=GS {} {} {}\n", name, tag, value));
        }
    }

    // Sequence lines
    for (name, seq) in &aln.sequences {
        out.push_str(&format!("{} {}\n", name, seq));
    }

    // #=GR annotations
    for (name, tags) in &aln.gr_annotations {
        for (tag, value) in tags {
            out.push_str(&format!("#=GR {} {} {}\n", name, tag, value));
        }
    }

    // #=GC annotations
    for (tag, value) in &aln.gc_annotations {
        out.push_str(&format!("#=GC {} {}\n", tag, value));
    }

    out.push_str("//\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stockholm_round_trip() {
        let aln = StockholmAlignment {
            sequences: vec![
                ("seq1".to_string(), "ACGU--AACCU".to_string()),
                ("seq2".to_string(), "ACGU--AACCU".to_string()),
            ],
            gc_annotations: HashMap::new(),
            gs_annotations: HashMap::new(),
            gr_annotations: HashMap::new(),
        };

        let written = write_stockholm(&aln);
        let parsed = parse_stockholm(&written).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].sequences.len(), 2);
        assert_eq!(parsed[0].sequences[0].0, "seq1");
        assert_eq!(parsed[0].sequences[0].1, "ACGU--AACCU");
        assert_eq!(parsed[0].sequences[1].0, "seq2");
        assert_eq!(parsed[0].sequences[1].1, "ACGU--AACCU");
    }

    #[test]
    fn stockholm_multi_block() {
        let input = "\
# STOCKHOLM 1.0
seq1 ACGT
seq2 TGCA
//
# STOCKHOLM 1.0
seqA AAAA
seqB CCCC
seqC GGGG
//
";
        let alns = parse_stockholm(input).unwrap();
        assert_eq!(alns.len(), 2);
        assert_eq!(alns[0].sequences.len(), 2);
        assert_eq!(alns[1].sequences.len(), 3);
        assert_eq!(alns[1].sequences[2].0, "seqC");
        assert_eq!(alns[1].sequences[2].1, "GGGG");
    }

    #[test]
    fn stockholm_annotations_preserved() {
        let input = "\
# STOCKHOLM 1.0
#=GF AU Eddy SR
#=GS seq1 AC PF00001
#=GS seq2 AC PF00002
seq1 ACGU..AACCU
seq2 ACGU..AACCU
#=GR seq1 PP 99998877665
#=GC SS_cons ..<<..>>...
//
";
        let alns = parse_stockholm(input).unwrap();
        assert_eq!(alns.len(), 1);
        let aln = &alns[0];
        assert_eq!(aln.sequences.len(), 2);

        // GS annotations
        assert_eq!(
            aln.gs_annotations.get("seq1").unwrap().get("AC").unwrap(),
            "PF00001"
        );
        assert_eq!(
            aln.gs_annotations.get("seq2").unwrap().get("AC").unwrap(),
            "PF00002"
        );

        // GR annotations
        assert_eq!(
            aln.gr_annotations.get("seq1").unwrap().get("PP").unwrap(),
            "99998877665"
        );

        // GC annotations
        assert_eq!(
            aln.gc_annotations.get("SS_cons").unwrap(),
            "..<<..>>..."
        );
    }

    #[test]
    fn stockholm_empty_input() {
        let alns = parse_stockholm("").unwrap();
        assert!(alns.is_empty());

        let alns = parse_stockholm("   \n  \n").unwrap();
        assert!(alns.is_empty());
    }

    #[test]
    fn stockholm_unterminated_error() {
        let input = "# STOCKHOLM 1.0\nseq1 ACGT\nseq2 TGCA\n";
        let result = parse_stockholm(input);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("unterminated"));
    }

    #[test]
    fn stockholm_interleaved_sequences() {
        let input = "\
# STOCKHOLM 1.0
seq1 ACGT
seq2 TGCA

seq1 AAAA
seq2 CCCC
//
";
        let alns = parse_stockholm(input).unwrap();
        assert_eq!(alns.len(), 1);
        assert_eq!(alns[0].sequences[0].1, "ACGTAAAA");
        assert_eq!(alns[0].sequences[1].1, "TGCACCCC");
    }
}
