//! In-memory FASTA parsing for WASM environments.
//!
//! File-based parsing (needletail) is unavailable in the browser, so this
//! module provides byte-slice alternatives that produce the same [`FastaStats`]
//! type used by `cyanea-seq`.

use cyanea_core::{CyaneaError, Result, Sequence};
use cyanea_seq::FastaStats;

use crate::error::{wasm_err, wasm_ok, wasm_result};

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

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

/// Parse FASTA from a string and return JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_fasta(data: &str) -> String {
    wasm_result(parse_fasta_bytes(data.as_bytes()))
}

/// GC content of a raw nucleotide string, returned as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn gc_content_json(seq: &str) -> String {
    wasm_ok(&gc_content(seq.as_bytes()))
}

/// Reverse complement of a DNA sequence string, returned as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn reverse_complement(seq: &str) -> String {
    match cyanea_seq::DnaSequence::new(seq.as_bytes()) {
        Ok(dna) => {
            let rc = dna.reverse_complement();
            wasm_ok(&String::from_utf8_lossy(rc.as_bytes()).into_owned())
        }
        Err(e) => wasm_err(e),
    }
}

/// Transcribe DNA to RNA, returned as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn transcribe(seq: &str) -> String {
    match cyanea_seq::DnaSequence::new(seq.as_bytes()) {
        Ok(dna) => {
            let rna = dna.transcribe();
            wasm_ok(&String::from_utf8_lossy(rna.as_bytes()).into_owned())
        }
        Err(e) => wasm_err(e),
    }
}

/// Translate DNA to protein (standard codon table), returned as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn translate(seq: &str) -> String {
    match cyanea_seq::DnaSequence::new(seq.as_bytes()) {
        Ok(dna) => match dna.translate() {
            Ok(protein) => wasm_ok(&String::from_utf8_lossy(protein.as_bytes()).into_owned()),
            Err(e) => wasm_err(e),
        },
        Err(e) => wasm_err(e),
    }
}

/// Validate a sequence against an alphabet ("dna", "rna", or "protein").
///
/// Returns JSON `{"ok": true}` if valid, or `{"error": "..."}` if invalid.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn validate(seq: &str, alphabet: &str) -> String {
    let result = match alphabet.to_ascii_lowercase().as_str() {
        "dna" => cyanea_seq::DnaSequence::new(seq.as_bytes()).map(|_| ()),
        "rna" => cyanea_seq::RnaSequence::new(seq.as_bytes()).map(|_| ()),
        "protein" => cyanea_seq::ProteinSequence::new(seq.as_bytes()).map(|_| ()),
        _ => Err(cyanea_core::CyaneaError::InvalidInput(format!(
            "unknown alphabet: {alphabet:?} (expected \"dna\", \"rna\", or \"protein\")"
        ))),
    };
    match result {
        Ok(()) => wasm_ok(&true),
        Err(e) => wasm_err(e),
    }
}

/// Parse FASTQ from a string and return JSON array of records.
///
/// Each record has `name`, `sequence`, and `quality` fields.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_fastq(data: &str) -> String {
    match parse_fastq_records(data.as_bytes()) {
        Ok(records) => wasm_ok(&records),
        Err(e) => wasm_err(e),
    }
}

/// A FASTQ record for JSON serialization.
#[derive(Debug, serde::Serialize)]
struct JsFastqRecord {
    name: String,
    sequence: String,
    quality: String,
}

fn parse_fastq_records(data: &[u8]) -> cyanea_core::Result<Vec<JsFastqRecord>> {
    let text = std::str::from_utf8(data)
        .map_err(|e| cyanea_core::CyaneaError::Parse(format!("invalid UTF-8: {e}")))?;
    let lines: Vec<&str> = text.lines().collect();

    if lines.is_empty() {
        return Err(cyanea_core::CyaneaError::Parse(
            "no FASTQ records found".into(),
        ));
    }

    let mut records = Vec::new();
    let mut i = 0;
    while i + 3 < lines.len() {
        let header = lines[i];
        if !header.starts_with('@') {
            return Err(cyanea_core::CyaneaError::Parse(format!(
                "expected '@' header at line {}, got: {}",
                i + 1,
                header
            )));
        }
        let name = header[1..].to_string();
        let sequence = lines[i + 1].to_string();
        // lines[i + 2] is the '+' separator
        let quality = lines[i + 3].to_string();
        records.push(JsFastqRecord {
            name,
            sequence,
            quality,
        });
        i += 4;
    }

    if records.is_empty() {
        return Err(cyanea_core::CyaneaError::Parse(
            "no FASTQ records found".into(),
        ));
    }

    Ok(records)
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

    #[test]
    fn reverse_complement_known() {
        let json = reverse_complement("ACGT");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"].as_str().unwrap(), "ACGT"); // palindrome
    }

    #[test]
    fn reverse_complement_asymmetric() {
        let json = reverse_complement("AACG");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"].as_str().unwrap(), "CGTT");
    }

    #[test]
    fn transcribe_known() {
        let json = transcribe("ACGT");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"].as_str().unwrap(), "ACGU");
    }

    #[test]
    fn translate_known() {
        let json = translate("ATGAAA");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"].as_str().unwrap(), "MK");
    }

    #[test]
    fn validate_dna_ok() {
        let json = validate("ACGTACGT", "dna");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], true);
    }

    #[test]
    fn validate_dna_invalid() {
        let json = validate("ACGTXYZ", "dna");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn validate_protein_ok() {
        let json = validate("HEAGAWGHEE", "protein");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], true);
    }

    #[test]
    fn validate_unknown_alphabet() {
        let json = validate("ACGT", "foobar");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].as_str().unwrap().contains("unknown alphabet"));
    }

    #[test]
    fn parse_fastq_single_record() {
        let json = parse_fastq("@seq1\nACGT\n+\nIIII\n");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let records = v["ok"].as_array().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0]["name"], "seq1");
        assert_eq!(records[0]["sequence"], "ACGT");
        assert_eq!(records[0]["quality"], "IIII");
    }

    #[test]
    fn parse_fastq_empty_error() {
        let json = parse_fastq("");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }
}
