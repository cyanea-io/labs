//! In-memory FASTA parsing for WASM environments.
//!
//! File-based parsing (needletail) is unavailable in the browser, so this
//! module provides byte-slice alternatives that produce the same [`FastaStats`]
//! type used by `cyanea-seq`.

use cyanea_core::{CyaneaError, Result};
use cyanea_seq::FastaStats;

use crate::error::{wasm_ok, wasm_result};

/// Parse FASTA-formatted bytes and compute summary statistics.
///
/// Recognises `>` header lines; everything else is treated as sequence data.
/// Returns an error if no sequences are found.
pub fn parse_fasta_bytes(data: &[u8]) -> Result<FastaStats> {
    let mut sequence_count: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut gc_count: u64 = 0;
    let mut in_sequence = false;

    for line in data.split(|&b| b == b'\n') {
        let line = if line.last() == Some(&b'\r') {
            &line[..line.len() - 1]
        } else {
            line
        };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            sequence_count += 1;
            in_sequence = true;
            continue;
        }
        if !in_sequence {
            continue;
        }
        for &base in line {
            total_bases += 1;
            match base {
                b'G' | b'g' | b'C' | b'c' => gc_count += 1,
                _ => {}
            }
        }
    }

    if sequence_count == 0 {
        return Err(CyaneaError::Parse("no FASTA sequences found".into()));
    }

    let gc_content = if total_bases > 0 {
        (gc_count as f64 / total_bases as f64) * 100.0
    } else {
        0.0
    };

    let avg_length = if sequence_count > 0 {
        total_bases as f64 / sequence_count as f64
    } else {
        0.0
    };

    Ok(FastaStats {
        sequence_count,
        total_bases,
        gc_content,
        avg_length,
    })
}

/// GC content (as a percentage) of a raw nucleotide byte slice.
pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
        .count();
    (gc as f64 / seq.len() as f64) * 100.0
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Parse FASTA from a string and return JSON.
pub fn parse_fasta(data: &str) -> String {
    wasm_result(parse_fasta_bytes(data.as_bytes()))
}

/// GC content of a raw nucleotide string, returned as JSON.
pub fn gc_content_json(seq: &str) -> String {
    wasm_ok(&gc_content(seq.as_bytes()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_sequence() {
        let data = b">seq1\nATCGATCG\n";
        let stats = parse_fasta_bytes(data).unwrap();
        assert_eq!(stats.sequence_count, 1);
        assert_eq!(stats.total_bases, 8);
        assert!((stats.gc_content - 50.0).abs() < 0.01);
        assert!((stats.avg_length - 8.0).abs() < 0.01);
    }

    #[test]
    fn multi_sequence() {
        let data = b">seq1\nATCGATCG\n>seq2\nGCGCGCGC\n";
        let stats = parse_fasta_bytes(data).unwrap();
        assert_eq!(stats.sequence_count, 2);
        assert_eq!(stats.total_bases, 16);
        assert!((stats.gc_content - 75.0).abs() < 0.01);
    }

    #[test]
    fn mixed_case() {
        let data = b">s\nacgtACGT\n";
        let stats = parse_fasta_bytes(data).unwrap();
        assert_eq!(stats.total_bases, 8);
        assert!((stats.gc_content - 50.0).abs() < 0.01);
    }

    #[test]
    fn empty_input_error() {
        assert!(parse_fasta_bytes(b"").is_err());
    }

    #[test]
    fn no_header_error() {
        assert!(parse_fasta_bytes(b"ACGTACGT\n").is_err());
    }

    #[test]
    fn gc_content_known() {
        assert!((gc_content(b"GCGC") - 100.0).abs() < 0.01);
        assert!((gc_content(b"ATAT") - 0.0).abs() < 0.01);
        assert!((gc_content(b"ACGT") - 50.0).abs() < 0.01);
    }

    #[test]
    fn gc_content_empty() {
        assert!((gc_content(b"") - 0.0).abs() < 0.01);
    }

    #[test]
    fn json_parse_fasta() {
        let json = parse_fasta(">s\nACGT\n");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"]["sequence_count"], 1);
    }

    #[test]
    fn json_gc_content() {
        let json = gc_content_json("GCGC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 100.0).abs() < 0.01);
    }
}
