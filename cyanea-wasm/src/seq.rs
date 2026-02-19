//! In-memory FASTA parsing for WASM environments.
//!
//! File-based parsing (needletail) is unavailable in the browser, so this
//! module provides byte-slice alternatives that produce the same [`FastaStats`]
//! type used by `cyanea-seq`.

use cyanea_core::{Annotated, CyaneaError, Result, Sequence};
use cyanea_seq::FastaStats;
use cyanea_seq::rna_structure;
use cyanea_seq::protein_properties;
use cyanea_seq::read_sim::{self, ReadSimConfig};
use cyanea_seq::codon::CodonUsage;
use cyanea_seq::assembly;

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

// ── Paired-end FASTQ ──────────────────────────────────────────────────────

/// Strip `/1`, `/2`, `_1`, `_2` mate suffixes from a read name.
///
/// Mirrors `cyanea_seq::strip_read_suffix` which is behind `#[cfg(feature = "std")]`.
fn strip_read_suffix(name: &str) -> &str {
    if name.len() >= 2 {
        let suffix = &name[name.len() - 2..];
        if suffix == "/1" || suffix == "/2" || suffix == "_1" || suffix == "_2" {
            return &name[..name.len() - 2];
        }
    }
    name
}

/// A paired FASTQ record for JSON serialization.
#[derive(Debug, serde::Serialize)]
struct JsPairedFastqRecord {
    r1: JsFastqRecord,
    r2: JsFastqRecord,
}

/// Parse paired FASTQ data from two separate strings.
///
/// `validation`: `"strict"`, `"relaxed"`, or `"none"`.
/// Returns JSON array of paired records.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_paired_fastq(r1_data: &str, r2_data: &str, validation: &str) -> String {
    match parse_paired_fastq_impl(r1_data, r2_data, validation) {
        Ok(pairs) => wasm_ok(&pairs),
        Err(e) => wasm_err(e),
    }
}

fn parse_paired_fastq_impl(
    r1_data: &str,
    r2_data: &str,
    validation: &str,
) -> cyanea_core::Result<Vec<JsPairedFastqRecord>> {
    let r1_records = parse_fastq_records(r1_data.as_bytes())?;
    let r2_records = parse_fastq_records(r2_data.as_bytes())?;

    if r1_records.len() != r2_records.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "R1 has {} records but R2 has {}",
            r1_records.len(),
            r2_records.len()
        )));
    }

    let pairs: cyanea_core::Result<Vec<_>> = r1_records
        .into_iter()
        .zip(r2_records)
        .map(|(r1, r2)| {
            match validation {
                "strict" => {
                    let n1 = r1.name.split_whitespace().next().unwrap_or(&r1.name);
                    let n2 = r2.name.split_whitespace().next().unwrap_or(&r2.name);
                    let s1 = strip_read_suffix(n1);
                    let s2 = strip_read_suffix(n2);
                    if s1 != s2 {
                        return Err(CyaneaError::InvalidInput(format!(
                            "mate pair name mismatch: '{}' vs '{}'",
                            n1, n2
                        )));
                    }
                    if !n1.ends_with("/1") {
                        return Err(CyaneaError::InvalidInput(format!(
                            "R1 read name '{}' does not end with /1",
                            n1
                        )));
                    }
                    if !n2.ends_with("/2") {
                        return Err(CyaneaError::InvalidInput(format!(
                            "R2 read name '{}' does not end with /2",
                            n2
                        )));
                    }
                }
                "relaxed" => {
                    let n1 = r1.name.split_whitespace().next().unwrap_or(&r1.name);
                    let n2 = r2.name.split_whitespace().next().unwrap_or(&r2.name);
                    let s1 = strip_read_suffix(n1);
                    let s2 = strip_read_suffix(n2);
                    if s1 != s2 {
                        return Err(CyaneaError::InvalidInput(format!(
                            "mate pair name mismatch: '{}' vs '{}'",
                            n1, n2
                        )));
                    }
                }
                "none" | _ => {}
            }
            Ok(JsPairedFastqRecord { r1, r2 })
        })
        .collect();
    pairs
}

/// Parse interleaved FASTQ data (alternating R1/R2 records).
///
/// `validation`: `"strict"`, `"relaxed"`, or `"none"`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_interleaved_fastq(data: &str, validation: &str) -> String {
    match parse_interleaved_fastq_impl(data, validation) {
        Ok(pairs) => wasm_ok(&pairs),
        Err(e) => wasm_err(e),
    }
}

fn parse_interleaved_fastq_impl(
    data: &str,
    validation: &str,
) -> cyanea_core::Result<Vec<JsPairedFastqRecord>> {
    let records = parse_fastq_records(data.as_bytes())?;
    if records.len() % 2 != 0 {
        return Err(CyaneaError::InvalidInput(format!(
            "odd number of records ({}) in interleaved FASTQ",
            records.len()
        )));
    }

    let mut pairs = Vec::with_capacity(records.len() / 2);
    let mut iter = records.into_iter();
    while let Some(r1) = iter.next() {
        let r2 = iter.next().unwrap(); // safe: we checked even count

        match validation {
            "strict" | "relaxed" => {
                let n1 = r1.name.split_whitespace().next().unwrap_or(&r1.name);
                let n2 = r2.name.split_whitespace().next().unwrap_or(&r2.name);
                let s1 = strip_read_suffix(n1);
                let s2 = strip_read_suffix(n2);
                if s1 != s2 {
                    return Err(CyaneaError::InvalidInput(format!(
                        "mate pair name mismatch: '{}' vs '{}'",
                        n1, n2
                    )));
                }
            }
            _ => {}
        }

        pairs.push(JsPairedFastqRecord { r1, r2 });
    }
    Ok(pairs)
}

// ── Single-end trimming ───────────────────────────────────────────────────

/// Trim configuration deserialized from JSON.
#[derive(Debug, serde::Deserialize)]
struct JsTrimConfig {
    min_quality: Option<f64>,
    window_size: Option<usize>,
    min_length: Option<usize>,
    max_length: Option<usize>,
    adapters: Option<Vec<String>>,
}

fn build_pipeline(config: &JsTrimConfig) -> cyanea_seq::trim::TrimPipeline {
    let mut p = cyanea_seq::trim::TrimPipeline::new();
    if let (Some(ws), Some(mq)) = (config.window_size, config.min_quality) {
        p = p.sliding_window(ws, mq);
    } else if let Some(mq) = config.min_quality {
        p = p.sliding_window(4, mq);
    }
    if let Some(min) = config.min_length {
        p = p.min_length(min);
    }
    if let Some(max) = config.max_length {
        p = p.max_length(max);
    }
    if let Some(ref adapters) = config.adapters {
        for a in adapters {
            p = p.adapter(a.as_bytes());
        }
    }
    p
}

/// Build a FastqRecord from raw WASM strings for pipeline processing.
fn make_fastq_record(
    name: &str,
    seq_str: &str,
    qual_str: &str,
) -> cyanea_core::Result<cyanea_seq::FastqRecord> {
    let sequence = cyanea_seq::DnaSequence::new(seq_str.as_bytes())?;
    let qual_bytes: Vec<u8> = qual_str.as_bytes().iter().map(|&q| q.saturating_sub(33)).collect();
    let quality = cyanea_seq::QualityScores::from_raw(qual_bytes);
    cyanea_seq::FastqRecord::new(name.to_string(), None, sequence, quality)
}

/// Trim single-end FASTQ records.
///
/// `config_json`: JSON object with optional fields: `min_quality`, `window_size`,
/// `min_length`, `max_length`, `adapters`.
/// Returns JSON array of trimmed records.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn trim_fastq(data: &str, config_json: &str) -> String {
    match trim_fastq_impl(data, config_json) {
        Ok(records) => wasm_ok(&records),
        Err(e) => wasm_err(e),
    }
}

fn trim_fastq_impl(data: &str, config_json: &str) -> cyanea_core::Result<Vec<JsFastqRecord>> {
    let config: JsTrimConfig = serde_json::from_str(config_json)
        .map_err(|e| CyaneaError::InvalidInput(format!("invalid trim config: {e}")))?;
    let pipeline = build_pipeline(&config);
    let records = parse_fastq_records(data.as_bytes())?;

    let mut results = Vec::new();
    for rec in &records {
        let fq = make_fastq_record(&rec.name, &rec.sequence, &rec.quality)?;
        if let Some(trimmed) = pipeline.process(&fq) {
            let quality: Vec<u8> = trimmed.quality().as_slice().iter().map(|&q| q + 33).collect();
            results.push(JsFastqRecord {
                name: trimmed.name().to_string(),
                sequence: String::from_utf8_lossy(trimmed.sequence().as_bytes()).into_owned(),
                quality: String::from_utf8_lossy(&quality).into_owned(),
            });
        }
    }
    Ok(results)
}

// ── Paired-end trimming ───────────────────────────────────────────────────

/// Paired-end trim statistics for JSON serialization.
#[derive(Debug, serde::Serialize)]
struct JsPairedTrimStats {
    total_input: usize,
    both_passed: usize,
    r1_only_passed: usize,
    r2_only_passed: usize,
    both_failed: usize,
    survival_rate: f64,
}

/// Paired-end trim result for JSON serialization.
#[derive(Debug, serde::Serialize)]
struct JsPairedTrimResult {
    pairs: Vec<JsPairedFastqRecord>,
    stats: JsPairedTrimStats,
}

/// Trim paired FASTQ records.
///
/// `config_json`: JSON trim config (see `trim_fastq`).
/// `orphan_policy`: `"drop_both"`, `"keep_first"`, or `"keep_second"`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn trim_paired_fastq(
    r1_data: &str,
    r2_data: &str,
    config_json: &str,
    orphan_policy: &str,
) -> String {
    match trim_paired_fastq_impl(r1_data, r2_data, config_json, orphan_policy) {
        Ok(result) => wasm_ok(&result),
        Err(e) => wasm_err(e),
    }
}

fn trim_paired_fastq_impl(
    r1_data: &str,
    r2_data: &str,
    config_json: &str,
    orphan_policy: &str,
) -> cyanea_core::Result<JsPairedTrimResult> {
    let config: JsTrimConfig = serde_json::from_str(config_json)
        .map_err(|e| CyaneaError::InvalidInput(format!("invalid trim config: {e}")))?;
    let pipeline = build_pipeline(&config);

    // Validate orphan_policy string early.
    match orphan_policy {
        "drop_both" | "keep_first" | "keep_second" => {}
        _ => {
            return Err(CyaneaError::InvalidInput(format!(
                "invalid orphan_policy: '{}' (expected 'drop_both', 'keep_first', or 'keep_second')",
                orphan_policy
            )))
        }
    }

    let r1_records = parse_fastq_records(r1_data.as_bytes())?;
    let r2_records = parse_fastq_records(r2_data.as_bytes())?;

    if r1_records.len() != r2_records.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "R1 has {} records but R2 has {}",
            r1_records.len(),
            r2_records.len()
        )));
    }

    let total_input = r1_records.len();
    let mut both_passed: usize = 0;
    let mut r1_only_passed: usize = 0;
    let mut r2_only_passed: usize = 0;
    let mut both_failed: usize = 0;
    let mut pairs = Vec::new();

    for (r1_raw, r2_raw) in r1_records.iter().zip(r2_records.iter()) {
        let r1_fq = make_fastq_record(&r1_raw.name, &r1_raw.sequence, &r1_raw.quality)?;
        let r2_fq = make_fastq_record(&r2_raw.name, &r2_raw.sequence, &r2_raw.quality)?;

        // Inline paired trim logic (process_paired is gated behind std in cyanea-seq).
        let r1_result = pipeline.process(&r1_fq);
        let r2_result = pipeline.process(&r2_fq);

        match (r1_result, r2_result) {
            (Some(t1), Some(t2)) => {
                both_passed += 1;
                let q1: Vec<u8> = t1.quality().as_slice().iter().map(|&q| q + 33).collect();
                let q2: Vec<u8> = t2.quality().as_slice().iter().map(|&q| q + 33).collect();
                pairs.push(JsPairedFastqRecord {
                    r1: JsFastqRecord {
                        name: t1.name().to_string(),
                        sequence: String::from_utf8_lossy(t1.sequence().as_bytes()).into_owned(),
                        quality: String::from_utf8_lossy(&q1).into_owned(),
                    },
                    r2: JsFastqRecord {
                        name: t2.name().to_string(),
                        sequence: String::from_utf8_lossy(t2.sequence().as_bytes()).into_owned(),
                        quality: String::from_utf8_lossy(&q2).into_owned(),
                    },
                });
            }
            (Some(_), None) => {
                if orphan_policy == "keep_first" {
                    r1_only_passed += 1;
                } else {
                    both_failed += 1;
                }
            }
            (None, Some(_)) => {
                if orphan_policy == "keep_second" {
                    r2_only_passed += 1;
                } else {
                    both_failed += 1;
                }
            }
            (None, None) => both_failed += 1,
        }
    }

    let survival_rate = if total_input > 0 {
        both_passed as f64 / total_input as f64
    } else {
        0.0
    };

    Ok(JsPairedTrimResult {
        pairs,
        stats: JsPairedTrimStats {
            total_input,
            both_passed,
            r1_only_passed,
            r2_only_passed,
            both_failed,
            survival_rate,
        },
    })
}

// ── MinHash ─────────────────────────────────────────────────────────────

/// Serializable MinHash sketch result.
#[derive(Debug, serde::Serialize)]
pub struct JsMinHashSketch {
    pub k: usize,
    pub sketch_size: usize,
    pub num_hashes: usize,
    pub hashes: Vec<u64>,
}

/// Serializable MinHash comparison result.
#[derive(Debug, serde::Serialize)]
pub struct JsMinHashComparison {
    pub jaccard: f64,
    pub containment_a_in_b: f64,
    pub containment_b_in_a: f64,
    pub ani: f64,
}

#[derive(Debug, serde::Serialize)]
pub struct JsRnaStructure {
    pub structure: String,
    pub pairs: Vec<(usize, usize)>,
    pub free_energy: Option<f64>,
    pub n_pairs: usize,
}

#[derive(Debug, serde::Serialize)]
pub struct JsProteinProperties {
    pub molecular_weight: f64,
    pub isoelectric_point: f64,
    pub gravy: f64,
    pub length: usize,
}

#[derive(Debug, serde::Serialize)]
pub struct JsSimulatedRead {
    pub name: String,
    pub sequence: String,
    pub quality: String,
    pub position: u64,
}

#[derive(Debug, serde::Serialize)]
pub struct JsCodonUsage {
    pub codons: std::collections::HashMap<String, usize>,
    pub total: usize,
}

#[derive(Debug, serde::Serialize)]
pub struct JsAssemblyStats {
    pub n_contigs: usize,
    pub total_length: usize,
    pub largest_contig: usize,
    pub smallest_contig: usize,
    pub gc_content: f64,
    pub n50: usize,
    pub l50: usize,
    pub n90: usize,
    pub l90: usize,
}

/// Create a MinHash sketch of a nucleotide sequence.
///
/// Returns a JSON object with k, sketch_size, num_hashes, and the hash values.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn minhash_sketch(seq: &str, k: usize, sketch_size: usize) -> String {
    match cyanea_seq::minhash::MinHash::from_sequence(seq.as_bytes(), k, sketch_size) {
        Ok(sketch) => {
            let js = JsMinHashSketch {
                k: sketch.k(),
                sketch_size: sketch.sketch_size(),
                num_hashes: sketch.len(),
                hashes: sketch.hashes().to_vec(),
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Compare two sequences using MinHash and return similarity metrics.
///
/// Returns Jaccard similarity, containment (both directions), and ANI estimate.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn minhash_compare(seq_a: &str, seq_b: &str, k: usize, sketch_size: usize) -> String {
    let sketch_a = match cyanea_seq::minhash::MinHash::from_sequence(
        seq_a.as_bytes(), k, sketch_size,
    ) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let sketch_b = match cyanea_seq::minhash::MinHash::from_sequence(
        seq_b.as_bytes(), k, sketch_size,
    ) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let jaccard = match sketch_a.jaccard(&sketch_b) {
        Ok(j) => j,
        Err(e) => return wasm_err(e),
    };
    let containment_a_in_b = match sketch_a.containment(&sketch_b) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    let containment_b_in_a = match sketch_b.containment(&sketch_a) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    let ani = match sketch_a.ani(&sketch_b) {
        Ok(a) => a,
        Err(_) => 0.0, // ANI undefined when Jaccard is zero
    };
    let js = JsMinHashComparison {
        jaccard,
        containment_a_in_b,
        containment_b_in_a,
        ani,
    };
    wasm_ok(&js)
}

// ── RNA structure prediction ────────────────────────────────────────────

/// Predict RNA secondary structure using the Nussinov algorithm (maximize base pairs).
///
/// Returns JSON with dot-bracket structure, base pairs, and pair count.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn rna_fold_nussinov(seq: &str) -> String {
    match rna_structure::nussinov(seq.as_bytes(), 3) {
        Ok(result) => {
            let pairs = result.structure.base_pairs();
            let n_pairs = result.structure.num_pairs();
            let js = JsRnaStructure {
                structure: result.structure.to_dot_bracket(),
                pairs,
                free_energy: None,
                n_pairs,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Predict RNA secondary structure using the Zuker MFE algorithm.
///
/// Returns JSON with dot-bracket structure, base pairs, free energy, and pair count.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn rna_fold_zuker(seq: &str) -> String {
    match rna_structure::zuker_mfe(seq.as_bytes()) {
        Ok(result) => {
            let pairs = result.structure.base_pairs();
            let n_pairs = result.structure.num_pairs();
            let js = JsRnaStructure {
                structure: result.structure.to_dot_bracket(),
                pairs,
                free_energy: Some(result.energy),
                n_pairs,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Protein properties ──────────────────────────────────────────────────

/// Compute basic protein sequence properties.
///
/// Returns JSON with molecular weight (estimated), isoelectric point, GRAVY, and length.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn protein_props(seq: &str) -> String {
    let pi = match protein_properties::isoelectric_point(seq.as_bytes()) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let gravy_val = match protein_properties::gravy(seq.as_bytes()) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    // Simple MW estimate: average amino acid MW ~128.16 Da minus water (18.02) per peptide bond
    let mw = seq.len() as f64 * 128.16 - (seq.len().saturating_sub(1)) as f64 * 18.02;
    let js = JsProteinProperties {
        molecular_weight: mw,
        isoelectric_point: pi,
        gravy: gravy_val,
        length: seq.len(),
    };
    wasm_ok(&js)
}

// ── Read simulation ─────────────────────────────────────────────────────

/// Simulate sequencing reads from a reference sequence.
///
/// `config_json`: JSON object with optional fields: `read_length`, `coverage`,
/// `error_rate`, `seed`. Returns JSON array of simulated reads.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn simulate_reads(ref_seq: &str, config_json: &str) -> String {
    match simulate_reads_impl(ref_seq, config_json) {
        Ok(reads) => wasm_ok(&reads),
        Err(e) => wasm_err(e),
    }
}

/// JSON config for read simulation.
#[derive(Debug, serde::Deserialize)]
struct JsReadSimConfig {
    read_length: Option<usize>,
    coverage: Option<f64>,
    error_rate: Option<f64>,
    seed: Option<u64>,
}

fn simulate_reads_impl(
    ref_seq: &str,
    config_json: &str,
) -> cyanea_core::Result<Vec<JsSimulatedRead>> {
    let js_config: JsReadSimConfig = serde_json::from_str(config_json)
        .map_err(|e| CyaneaError::InvalidInput(format!("invalid read sim config: {e}")))?;

    let mut config = ReadSimConfig::default();
    if let Some(rl) = js_config.read_length {
        config.read_length = rl;
    }
    if let Some(cov) = js_config.coverage {
        config.coverage = cov;
    }
    if let Some(er) = js_config.error_rate {
        config.error_rate = er;
    }
    if let Some(s) = js_config.seed {
        config.seed = s;
    }
    config.paired = false; // single-end for WASM simplicity

    let reads = read_sim::simulate_reads(ref_seq.as_bytes(), "wasm_ref", &config);
    let js_reads: Vec<JsSimulatedRead> = reads
        .into_iter()
        .map(|r| JsSimulatedRead {
            name: r.name,
            sequence: String::from_utf8_lossy(&r.sequence).into_owned(),
            quality: String::from_utf8_lossy(&r.quality).into_owned(),
            position: r.true_position,
        })
        .collect();
    Ok(js_reads)
}

// ── Codon usage ─────────────────────────────────────────────────────────

/// Compute codon usage from a coding DNA sequence.
///
/// The sequence should be in-frame coding DNA. Returns JSON with codon counts and total.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn codon_usage(seq: &str) -> String {
    let usage = CodonUsage::from_sequence(seq.as_bytes());
    let bases: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut codons = std::collections::HashMap::new();
    let mut total: usize = 0;

    for &b1 in &bases {
        for &b2 in &bases {
            for &b3 in &bases {
                let codon = [b1, b2, b3];
                let count = usage.count(&codon) as usize;
                if count > 0 {
                    let key = String::from_utf8_lossy(&codon).into_owned();
                    codons.insert(key, count);
                    total += count;
                }
            }
        }
    }

    let js = JsCodonUsage { codons, total };
    wasm_ok(&js)
}

// ── Assembly statistics ─────────────────────────────────────────────────

/// Compute assembly statistics from a set of contigs.
///
/// `contigs_json`: JSON array of contig sequences (strings).
/// Returns JSON with n_contigs, total_length, N50, L50, N90, L90, GC content, etc.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn assembly_stats_json(contigs_json: &str) -> String {
    match assembly_stats_impl(contigs_json) {
        Ok(stats) => wasm_ok(&stats),
        Err(e) => wasm_err(e),
    }
}

fn assembly_stats_impl(contigs_json: &str) -> cyanea_core::Result<JsAssemblyStats> {
    let contigs: Vec<String> = serde_json::from_str(contigs_json)
        .map_err(|e| CyaneaError::InvalidInput(format!("invalid contigs JSON: {e}")))?;

    let contig_bytes: Vec<&[u8]> = contigs.iter().map(|s| s.as_bytes()).collect();
    let stats = assembly::assembly_stats(&contig_bytes)?;

    Ok(JsAssemblyStats {
        n_contigs: stats.n_contigs,
        total_length: stats.total_length,
        largest_contig: stats.largest_contig,
        smallest_contig: stats.smallest_contig,
        gc_content: stats.gc_content,
        n50: stats.n50,
        l50: stats.l50,
        n90: stats.n90,
        l90: stats.l90,
    })
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

    // --- Paired FASTQ parsing ---

    #[test]
    fn parse_paired_fastq_matching() {
        let r1 = "@read1/1\nACGT\n+\nIIII\n";
        let r2 = "@read1/2\nTGCA\n+\nIIII\n";
        let json = parse_paired_fastq(r1, r2, "relaxed");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let pairs = v["ok"].as_array().unwrap();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0]["r1"]["sequence"], "ACGT");
        assert_eq!(pairs[0]["r2"]["sequence"], "TGCA");
    }

    #[test]
    fn parse_paired_fastq_mismatched_count() {
        let r1 = "@r1/1\nACGT\n+\nIIII\n@r2/1\nTGCA\n+\nIIII\n";
        let r2 = "@r1/2\nACGT\n+\nIIII\n";
        let json = parse_paired_fastq(r1, r2, "none");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn parse_paired_fastq_strict_valid() {
        let r1 = "@read1/1\nACGT\n+\nIIII\n";
        let r2 = "@read1/2\nTGCA\n+\nIIII\n";
        let json = parse_paired_fastq(r1, r2, "strict");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].is_array());
    }

    #[test]
    fn parse_paired_fastq_strict_missing_suffix() {
        let r1 = "@read1\nACGT\n+\nIIII\n";
        let r2 = "@read1/2\nTGCA\n+\nIIII\n";
        let json = parse_paired_fastq(r1, r2, "strict");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn parse_paired_fastq_relaxed_mismatch() {
        let r1 = "@read1/1\nACGT\n+\nIIII\n";
        let r2 = "@read2/2\nTGCA\n+\nIIII\n";
        let json = parse_paired_fastq(r1, r2, "relaxed");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn parse_paired_fastq_none_no_validation() {
        let r1 = "@foo\nACGT\n+\nIIII\n";
        let r2 = "@bar\nTGCA\n+\nIIII\n";
        let json = parse_paired_fastq(r1, r2, "none");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"].is_array());
    }

    // --- Interleaved FASTQ ---

    #[test]
    fn parse_interleaved_fastq_valid() {
        let data = "@read1/1\nACGT\n+\nIIII\n@read1/2\nTGCA\n+\nIIII\n";
        let json = parse_interleaved_fastq(data, "relaxed");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let pairs = v["ok"].as_array().unwrap();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0]["r1"]["sequence"], "ACGT");
        assert_eq!(pairs[0]["r2"]["sequence"], "TGCA");
    }

    #[test]
    fn parse_interleaved_fastq_odd_count() {
        let data = "@r1/1\nACGT\n+\nIIII\n@r1/2\nTGCA\n+\nIIII\n@r2/1\nAAAA\n+\nIIII\n";
        let json = parse_interleaved_fastq(data, "none");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
        assert!(v["error"].as_str().unwrap().contains("odd"));
    }

    // --- Trim single-end ---

    #[test]
    fn trim_fastq_basic() {
        // All high quality, should pass through
        let data = "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
        let config = r#"{"min_quality": 20, "window_size": 4, "min_length": 4}"#;
        let json = trim_fastq(data, config);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let records = v["ok"].as_array().unwrap();
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn trim_fastq_filters_short() {
        // After trimming low quality, record should be too short
        // Quality: !!!! = Phred 0, which at min_quality 20 gets trimmed to nothing
        let data = "@r1\nACGT\n+\n!!!!\n";
        let config = r#"{"min_quality": 20, "window_size": 4, "min_length": 4}"#;
        let json = trim_fastq(data, config);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let records = v["ok"].as_array().unwrap();
        assert_eq!(records.len(), 0);
    }

    // --- Trim paired-end ---

    #[test]
    fn trim_paired_fastq_both_pass() {
        let r1 = "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
        let r2 = "@r1\nTGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIII\n";
        let config = r#"{"min_quality": 20, "window_size": 4, "min_length": 4}"#;
        let json = trim_paired_fastq(r1, r2, config, "drop_both");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["stats"]["total_input"], 1);
        assert_eq!(result["stats"]["both_passed"], 1);
        assert_eq!(result["pairs"].as_array().unwrap().len(), 1);
    }

    #[test]
    fn trim_paired_fastq_both_fail() {
        let r1 = "@r1\nACGT\n+\n!!!!\n";
        let r2 = "@r1\nTGCA\n+\n!!!!\n";
        let config = r#"{"min_quality": 20, "window_size": 4, "min_length": 4}"#;
        let json = trim_paired_fastq(r1, r2, config, "drop_both");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["stats"]["both_failed"], 1);
        assert_eq!(result["pairs"].as_array().unwrap().len(), 0);
    }

    #[test]
    fn trim_paired_fastq_one_fails_drop_both() {
        let r1 = "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
        let r2 = "@r1\nTGCA\n+\n!!!!\n";
        let config = r#"{"min_quality": 20, "window_size": 4, "min_length": 4}"#;
        let json = trim_paired_fastq(r1, r2, config, "drop_both");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        // R1 passes but R2 fails, so with drop_both both get dropped
        assert_eq!(result["pairs"].as_array().unwrap().len(), 0);
    }

    #[test]
    fn trim_paired_fastq_invalid_policy() {
        let r1 = "@r1\nACGT\n+\nIIII\n";
        let r2 = "@r1\nTGCA\n+\nIIII\n";
        let config = r#"{"min_length": 1}"#;
        let json = trim_paired_fastq(r1, r2, config, "invalid");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn trim_paired_fastq_survival_rate() {
        let r1 = "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
                  @r2\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
        let r2 = "@r1\nTGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIII\n\
                  @r2\nTGCA\n+\n!!!!\n";
        let config = r#"{"min_quality": 20, "window_size": 4, "min_length": 4}"#;
        let json = trim_paired_fastq(r1, r2, config, "drop_both");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let rate = v["ok"]["stats"]["survival_rate"].as_f64().unwrap();
        assert!((rate - 0.5).abs() < 1e-10);
    }

    // --- MinHash ---

    #[test]
    fn minhash_sketch_basic() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let json = minhash_sketch(seq, 4, 10);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["k"], 4);
        assert_eq!(obj["sketch_size"], 10);
        assert!(obj["num_hashes"].as_u64().unwrap() > 0);
        assert!(obj["hashes"].as_array().unwrap().len() > 0);
    }

    #[test]
    fn minhash_sketch_invalid_k() {
        let json = minhash_sketch("ACGT", 0, 10);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn minhash_compare_identical() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let json = minhash_compare(seq, seq, 4, 10);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!((obj["jaccard"].as_f64().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn minhash_compare_different() {
        let seq_a = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq_b = "TGCATGCATGCATGCATGCATGCATGCATGCA";
        let json = minhash_compare(seq_a, seq_b, 4, 10);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!(obj["jaccard"].as_f64().unwrap() < 1.0);
    }

    // --- RNA structure prediction ---

    #[test]
    fn rna_fold_nussinov_basic() {
        let json = rna_fold_nussinov("GGGAAACCC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!(obj["n_pairs"].as_u64().unwrap() > 0);
        assert!(!obj["structure"].as_str().unwrap().is_empty());
        assert!(obj["free_energy"].is_null());
    }

    #[test]
    fn rna_fold_nussinov_no_pairs() {
        let json = rna_fold_nussinov("AAAA");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_pairs"].as_u64().unwrap(), 0);
        assert_eq!(obj["pairs"].as_array().unwrap().len(), 0);
    }

    #[test]
    fn rna_fold_zuker_basic() {
        let json = rna_fold_zuker("GGGAAACCC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!(!obj["structure"].as_str().unwrap().is_empty());
        assert!(obj["free_energy"].is_number());
    }

    #[test]
    fn rna_fold_zuker_pairs() {
        let json = rna_fold_zuker("GGGAAACCC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        let pairs = obj["pairs"].as_array().unwrap();
        // Each pair should be [i, j] with i < j
        for pair in pairs {
            let arr = pair.as_array().unwrap();
            assert!(arr[0].as_u64().unwrap() < arr[1].as_u64().unwrap());
        }
    }

    // --- Protein properties ---

    #[test]
    fn protein_properties_basic() {
        let json = protein_props("MKWVTFISLLFLFSSAYS");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!(obj["isoelectric_point"].as_f64().unwrap() > 0.0);
        assert!(obj["isoelectric_point"].as_f64().unwrap() < 14.0);
        assert!(obj["molecular_weight"].as_f64().unwrap() > 0.0);
        assert_eq!(obj["length"].as_u64().unwrap(), 18);
    }

    #[test]
    fn protein_properties_short() {
        let json = protein_props("MK");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["length"].as_u64().unwrap(), 2);
        assert!(obj["molecular_weight"].as_f64().unwrap() > 0.0);
        assert!(obj["gravy"].is_number());
    }

    // --- Read simulation ---

    #[test]
    fn simulate_reads_basic() {
        // Reference needs to be longer than read_length
        let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let config = r#"{"read_length": 10, "coverage": 5.0, "error_rate": 0.0, "seed": 1}"#;
        let json = simulate_reads(ref_seq, config);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let reads = v["ok"].as_array().unwrap();
        assert!(!reads.is_empty());
        // Each read should have name, sequence, quality, and position
        for read in reads {
            assert!(read["name"].is_string());
            assert!(read["sequence"].as_str().unwrap().len() > 0);
            assert!(read["quality"].as_str().unwrap().len() > 0);
        }
    }

    // --- Codon usage ---

    #[test]
    fn codon_usage_basic() {
        let json = codon_usage("ATGAAAGGG");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["total"].as_u64().unwrap(), 3);
        let codons = obj["codons"].as_object().unwrap();
        assert_eq!(codons["ATG"].as_u64().unwrap(), 1);
        assert_eq!(codons["AAA"].as_u64().unwrap(), 1);
        assert_eq!(codons["GGG"].as_u64().unwrap(), 1);
    }

    // --- Assembly statistics ---

    #[test]
    fn assembly_stats_basic() {
        let contigs = r#"["ACGTACGTACGTACGT", "ACGTACGT", "ACGT"]"#;
        let json = assembly_stats_json(contigs);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_contigs"].as_u64().unwrap(), 3);
        assert_eq!(obj["total_length"].as_u64().unwrap(), 28);
        assert_eq!(obj["largest_contig"].as_u64().unwrap(), 16);
        assert_eq!(obj["smallest_contig"].as_u64().unwrap(), 4);
        assert!(obj["n50"].as_u64().unwrap() > 0);
    }

    #[test]
    fn assembly_stats_single() {
        let contigs = r#"["ACGTACGTACGTACGT"]"#;
        let json = assembly_stats_json(contigs);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_contigs"].as_u64().unwrap(), 1);
        assert_eq!(obj["total_length"].as_u64().unwrap(), 16);
        assert_eq!(obj["n50"].as_u64().unwrap(), 16);
        assert_eq!(obj["l50"].as_u64().unwrap(), 1);
    }
}
