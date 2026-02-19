//! Clustal format parser and writer.
//!
//! Parses ClustalW / Clustal Omega multiple sequence alignment output.
//! The format starts with a `CLUSTAL` header line, followed by blocks
//! of interleaved sequence data. Each block contains one line per
//! sequence (name + aligned fragment) and optionally a conservation
//! line using `*`, `:`, and `.` characters.

use cyanea_core::{CyaneaError, Result};

/// A parsed Clustal multiple sequence alignment.
#[derive(Debug, Clone)]
pub struct ClustalAlignment {
    /// Aligned sequences as `(name, aligned_sequence)` pairs.
    pub sequences: Vec<(String, String)>,
    /// Conservation string assembled from all blocks.
    /// Uses `*` (fully conserved), `:` (strong), `.` (weak).
    pub conservation: Option<String>,
}

/// Parse a Clustal-format alignment from a string.
///
/// Expects the first non-blank line to start with `CLUSTAL`.
/// Returns an error if the header is missing.
///
/// # Examples
///
/// ```
/// # use cyanea_io::clustal::parse_clustal;
/// let input = "CLUSTAL W (1.83) multiple sequence alignment\n\n\
///              seq1    ACGT\n\
///              seq2    ACGT\n\
///              \x20       ****\n";
/// let aln = parse_clustal(input).unwrap();
/// assert_eq!(aln.sequences.len(), 2);
/// ```
pub fn parse_clustal(input: &str) -> Result<ClustalAlignment> {
    let mut lines = input.lines().peekable();

    // Find and validate the CLUSTAL header
    let mut found_header = false;
    while let Some(line) = lines.peek() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            lines.next();
            continue;
        }
        if trimmed.starts_with("CLUSTAL") {
            found_header = true;
            lines.next();
            break;
        } else {
            break;
        }
    }

    if !found_header {
        return Err(CyaneaError::Parse(
            "missing CLUSTAL header line".to_string(),
        ));
    }

    let mut seq_order: Vec<String> = Vec::new();
    let mut seq_data: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();
    let mut conservation = String::new();
    let mut block_seq_count = 0usize;
    let mut expected_seqs_per_block: Option<usize> = None;

    for line in lines {
        let line = line.trim_end();

        // Blank line = end of block
        if line.trim().is_empty() {
            if block_seq_count > 0 {
                if expected_seqs_per_block.is_none() {
                    expected_seqs_per_block = Some(block_seq_count);
                }
                block_seq_count = 0;
            }
            continue;
        }

        // Try to determine if this is a conservation line.
        // Conservation lines contain only spaces, '*', ':', '.', and
        // are typically indented to align under the sequence data.
        // They don't start with an alphanumeric name character.
        let trimmed = line.trim();
        if !trimmed.is_empty() && trimmed.chars().all(|c| c == '*' || c == ':' || c == '.' || c == ' ') {
            conservation.push_str(trimmed);
            continue;
        }

        // Sequence line: name sequence [optional_number]
        // The name is the first whitespace-delimited token
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }

        // If the line has at least 2 tokens and the first looks like a name
        if parts.len() >= 2 {
            let name = parts[0];
            // The sequence is the second token (ignore optional trailing count)
            let seq = parts[1];

            if !seq_data.contains_key(name) {
                seq_order.push(name.to_string());
            }
            seq_data
                .entry(name.to_string())
                .and_modify(|v| v.push_str(seq))
                .or_insert_with(|| seq.to_string());
            block_seq_count += 1;
        }
    }

    let sequences: Vec<(String, String)> = seq_order
        .into_iter()
        .map(|name| {
            let seq = seq_data.remove(&name).unwrap_or_default();
            (name, seq)
        })
        .collect();

    let conservation = if conservation.is_empty() {
        None
    } else {
        Some(conservation)
    };

    Ok(ClustalAlignment {
        sequences,
        conservation,
    })
}

/// Write a Clustal alignment block as a string.
///
/// Produces output with a `CLUSTAL` header, followed by interleaved
/// blocks of up to 60 residues per line.
///
/// # Examples
///
/// ```
/// # use cyanea_io::clustal::{ClustalAlignment, write_clustal};
/// let aln = ClustalAlignment {
///     sequences: vec![("seq1".to_string(), "ACGT".to_string())],
///     conservation: None,
/// };
/// let output = write_clustal(&aln);
/// assert!(output.starts_with("CLUSTAL"));
/// ```
pub fn write_clustal(aln: &ClustalAlignment) -> String {
    const BLOCK_SIZE: usize = 60;

    let mut out = String::new();
    out.push_str("CLUSTAL W multiple sequence alignment\n\n");

    if aln.sequences.is_empty() {
        return out;
    }

    // Find the longest name for padding
    let max_name_len = aln.sequences.iter().map(|(n, _)| n.len()).max().unwrap_or(0);
    let pad = max_name_len + 4;

    let total_len = aln.sequences.iter().map(|(_, s)| s.len()).max().unwrap_or(0);
    let conservation = aln.conservation.as_deref().unwrap_or("");

    let mut offset = 0;
    while offset < total_len {
        let end = std::cmp::min(offset + BLOCK_SIZE, total_len);

        for (name, seq) in &aln.sequences {
            let fragment = if offset < seq.len() {
                &seq[offset..std::cmp::min(end, seq.len())]
            } else {
                ""
            };
            out.push_str(&format!("{:<width$}{}\n", name, fragment, width = pad));
        }

        // Conservation line
        if !conservation.is_empty() && offset < conservation.len() {
            let cons_fragment = &conservation[offset..std::cmp::min(end, conservation.len())];
            out.push_str(&format!("{:width$}{}\n", "", cons_fragment, width = pad));
        }

        out.push('\n');
        offset = end;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clustal_round_trip() {
        let aln = ClustalAlignment {
            sequences: vec![
                ("seq1".to_string(), "ACGT--AACCU".to_string()),
                ("seq2".to_string(), "ACGU--AACCU".to_string()),
            ],
            conservation: Some("****  *****".to_string()),
        };

        let written = write_clustal(&aln);
        let parsed = parse_clustal(&written).unwrap();
        assert_eq!(parsed.sequences.len(), 2);
        assert_eq!(parsed.sequences[0].0, "seq1");
        assert_eq!(parsed.sequences[0].1, "ACGT--AACCU");
        assert_eq!(parsed.sequences[1].0, "seq2");
        assert_eq!(parsed.sequences[1].1, "ACGU--AACCU");
        assert!(parsed.conservation.is_some());
        assert_eq!(parsed.conservation.unwrap(), "****  *****");
    }

    #[test]
    fn clustal_multi_block_interleaved() {
        let input = "\
CLUSTAL W (1.83) multiple sequence alignment

seq1    ACGT
seq2    TGCA
        *..*

seq1    AAAA
seq2    CCCC
";
        let aln = parse_clustal(input).unwrap();
        assert_eq!(aln.sequences.len(), 2);
        assert_eq!(aln.sequences[0].1, "ACGTAAAA");
        assert_eq!(aln.sequences[1].1, "TGCACCCC");
    }

    #[test]
    fn clustal_conservation_preserved() {
        let input = "\
CLUSTAL W (1.83) multiple sequence alignment

seq1    ACGT
seq2    ACGT
        ****
";
        let aln = parse_clustal(input).unwrap();
        assert_eq!(aln.conservation, Some("****".to_string()));
    }

    #[test]
    fn clustal_missing_header_error() {
        let input = "seq1    ACGT\nseq2    ACGT\n";
        let result = parse_clustal(input);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("CLUSTAL"));
    }

    #[test]
    fn clustal_empty_sequences() {
        let input = "CLUSTAL W (1.83) multiple sequence alignment\n\n";
        let aln = parse_clustal(input).unwrap();
        assert!(aln.sequences.is_empty());
        assert!(aln.conservation.is_none());
    }
}
