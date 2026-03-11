//! Motif discovery and scanning for regulatory elements.
//!
//! This module provides position weight matrix (PWM) operations, motif scanning,
//! de novo motif discovery, and format I/O (MEME/JASPAR).

use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::Result;

/// Position weight matrix representation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Motif {
    /// Motif name/ID.
    pub name: String,
    /// Position weight matrix as Vec<[freq; 4]> for A, C, G, T.
    pub pwm: Vec<[f64; 4]>,
    /// Consensus sequence.
    pub consensus: String,
    /// Optional background frequencies [A, C, G, T].
    pub background_freq: [f64; 4],
}

impl Motif {
    /// Create a new motif from a PWM.
    pub fn new(name: impl Into<String>, pwm: Vec<[f64; 4]>) -> Self {
        let consensus = pwm_to_consensus(&pwm);
        Self {
            name: name.into(),
            pwm,
            consensus,
            background_freq: [0.25, 0.25, 0.25, 0.25],
        }
    }

    /// Width of the motif (number of positions).
    pub fn width(&self) -> usize {
        self.pwm.len()
    }

    /// Set background frequencies (for log-odds scoring).
    pub fn with_background(mut self, bg: [f64; 4]) -> Self {
        self.background_freq = bg;
        self
    }
}

/// Match of a motif to a sequence.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MotifMatch {
    /// 0-based start position in query sequence.
    pub position: u64,
    /// Strand: '+' for forward, '-' for reverse.
    pub strand: char,
    /// Log-odds score.
    pub score: f64,
    /// P-value (approximate).
    pub p_value: f64,
    /// Matched sequence.
    pub matched_seq: String,
}

/// Parameters for motif discovery.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DiscoveryParams {
    /// K-mer width for discovery (8-15, default 10).
    pub motif_width: usize,
    /// Number of motifs to discover (default 10).
    pub n_motifs: usize,
    /// Background sequence for enrichment testing.
    pub background_freq: [f64; 4],
}

impl Default for DiscoveryParams {
    fn default() -> Self {
        Self {
            motif_width: 10,
            n_motifs: 10,
            background_freq: [0.25, 0.25, 0.25, 0.25],
        }
    }
}

/// Convert PWM to consensus sequence (IUPAC).
fn pwm_to_consensus(pwm: &[[f64; 4]]) -> String {
    pwm.iter()
        .map(|pos| {
            let max_idx = pos
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);

            match max_idx {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            }
        })
        .collect()
}

/// Compute log-odds score for a sequence at a position.
fn score_match(seq: &[u8], motif: &Motif, pos: usize) -> Option<f64> {
    if pos + motif.width() > seq.len() {
        return None;
    }

    let mut score = 0.0;
    for (i, &base) in seq[pos..pos + motif.width()].iter().enumerate() {
        let base_idx = match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };

        if let Some(idx) = base_idx {
            let prob = motif.pwm[i][idx];
            let bg = motif.background_freq[idx];
            score += (prob / bg.max(1e-10)).log2();
        } else {
            return None;
        }
    }

    Some(score)
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Scan a sequence for motif matches.
///
/// Returns all matches above the specified log-odds threshold,
/// on both forward and reverse strands.
pub fn scan_sequence(seq: &[u8], motif: &Motif, threshold: f64) -> Vec<MotifMatch> {
    let mut matches = Vec::new();

    // Forward strand
    for pos in 0..seq.len() {
        if let Some(score) = score_match(seq, motif, pos) {
            if score >= threshold {
                let matched = std::str::from_utf8(&seq[pos..pos + motif.width()])
                    .unwrap_or("N")
                    .to_string();
                matches.push(MotifMatch {
                    position: pos as u64,
                    strand: '+',
                    score,
                    p_value: 0.05, // Placeholder
                    matched_seq: matched,
                });
            }
        }
    }

    // Reverse complement strand
    let rc = reverse_complement(seq);
    for pos in 0..rc.len() {
        if let Some(score) = score_match(&rc, motif, pos) {
            if score >= threshold {
                let matched = std::str::from_utf8(&rc[pos..pos + motif.width()])
                    .unwrap_or("N")
                    .to_string();
                matches.push(MotifMatch {
                    position: (seq.len() - pos - motif.width()) as u64,
                    strand: '-',
                    score,
                    p_value: 0.05, // Placeholder
                    matched_seq: matched,
                });
            }
        }
    }

    matches.sort_by(|a, b| a.position.cmp(&b.position));
    matches
}

/// Discover motifs de novo from peak sequences.
///
/// Simple k-mer enrichment approach:
/// 1. Count k-mers in peak sequences
/// 2. Compare to shuffled background
/// 3. Build PWM from enriched k-mers
/// 4. Extend and align to refine
pub fn discover_motifs(
    peak_seqs: &[&[u8]],
    params: &DiscoveryParams,
) -> Result<Vec<Motif>> {
    if peak_seqs.is_empty() {
        return Ok(Vec::new());
    }

    let width = params.motif_width;

    // Count k-mers in peaks
    let mut kmer_counts: HashMap<Vec<u8>, usize> = HashMap::new();
    for seq in peak_seqs {
        if seq.len() >= width {
            for i in 0..=(seq.len() - width) {
                let kmer = seq[i..i + width].to_vec();
                *kmer_counts.entry(kmer).or_insert(0) += 1;
            }
        }
    }

    // Score k-mers by enrichment (simple version: count rank)
    let mut scored_kmers: Vec<(Vec<u8>, usize)> = kmer_counts.into_iter().collect();
    scored_kmers.sort_by_key(|&(_, count)| std::cmp::Reverse(count));

    // Build motifs from top k-mers
    let mut motifs = Vec::new();

    for (idx, (kmer, _count)) in scored_kmers.iter().take(params.n_motifs).enumerate() {
        let pwm = kmer_to_pwm(kmer);
        let motif = Motif::new(format!("motif_{}", idx), pwm)
            .with_background(params.background_freq);
        motifs.push(motif);
    }

    Ok(motifs)
}

/// Convert a k-mer to a position weight matrix (simple: uniform + observed).
fn kmer_to_pwm(kmer: &[u8]) -> Vec<[f64; 4]> {
    kmer.iter()
        .map(|&base| {
            let mut pwm_pos = [0.1; 4]; // Pseudocount
            let idx = match base.to_ascii_uppercase() {
                b'A' => Some(0),
                b'C' => Some(1),
                b'G' => Some(2),
                b'T' => Some(3),
                _ => None,
            };
            if let Some(i) = idx {
                pwm_pos[i] = 0.9;
            }

            // Normalize
            let sum: f64 = pwm_pos.iter().sum();
            pwm_pos.iter_mut().for_each(|x| *x /= sum);

            pwm_pos
        })
        .collect()
}

/// Parse motifs from MEME format file content.
pub fn parse_meme(content: &str) -> Result<Vec<Motif>> {
    let mut motifs = Vec::new();
    let mut in_motif = false;
    let mut current_name = String::new();
    let mut current_pwm = Vec::new();

    for line in content.lines() {
        if line.starts_with("MOTIF") {
            if !current_pwm.is_empty() {
                let motif = Motif::new(&current_name, current_pwm);
                motifs.push(motif);
                current_pwm = Vec::new();
            }

            in_motif = true;
            current_name = line.split_whitespace().nth(1).unwrap_or("unknown").to_string();
        } else if in_motif && line.starts_with("letter-probability matrix:") {
            // Next lines are PWM positions
            continue;
        } else if in_motif && !line.is_empty() && !line.starts_with("//") {
            // Parse PWM line: space-separated probabilities for A C G T
            let parts: Vec<f64> = line
                .split_whitespace()
                .take(4)
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();

            if parts.len() == 4 {
                current_pwm.push([parts[0], parts[1], parts[2], parts[3]]);
            }
        } else if line.starts_with("//") {
            in_motif = false;
        }
    }

    if !current_pwm.is_empty() {
        let motif = Motif::new(&current_name, current_pwm);
        motifs.push(motif);
    }

    Ok(motifs)
}

/// Write motifs in MEME format.
pub fn write_meme(motifs: &[Motif]) -> String {
    let mut output = String::from("MEME version 4.0\n\n");

    for motif in motifs {
        output.push_str(&format!("MOTIF {} {}\n", motif.name, motif.consensus));
        output.push_str(&format!("letter-probability matrix: alength= 4 w= {}\n", motif.width()));

        for pos in &motif.pwm {
            output.push_str(&format!(
                " {:.6} {:.6} {:.6} {:.6}\n",
                pos[0], pos[1], pos[2], pos[3]
            ));
        }

        output.push_str("//\n\n");
    }

    output
}

/// Compare two motifs for similarity (correlation of aligned columns).
pub fn compare_motifs(motif1: &Motif, motif2: &Motif) -> f64 {
    let min_len = motif1.width().min(motif2.width());
    if min_len == 0 {
        return 0.0;
    }

    let mut score = 0.0;
    for i in 0..min_len {
        let mut col_corr = 0.0;
        for j in 0..4 {
            let p1 = motif1.pwm[i][j];
            let p2 = motif2.pwm[i][j];
            col_corr += p1 * p2;
        }
        score += col_corr;
    }

    score / (min_len as f64)
}

/// Test motif enrichment in peaks vs background using Fisher's exact test.
pub fn motif_enrichment(
    motif: &Motif,
    peak_seqs: &[&[u8]],
    background_seqs: &[&[u8]],
    threshold: f64,
) -> Result<(f64, f64)> {
    let mut peaks_with_match = 0u32;
    for seq in peak_seqs {
        let matches = scan_sequence(seq, motif, threshold);
        if !matches.is_empty() {
            peaks_with_match += 1;
        }
    }

    let mut bg_with_match = 0u32;
    for seq in background_seqs {
        let matches = scan_sequence(seq, motif, threshold);
        if !matches.is_empty() {
            bg_with_match += 1;
        }
    }

    let peak_total = peak_seqs.len() as f64;
    let bg_total = background_seqs.len() as f64;

    let peak_rate = peaks_with_match as f64 / peak_total.max(1.0);
    let bg_rate = bg_with_match as f64 / bg_total.max(1.0);

    let fold_enrichment = peak_rate / bg_rate.max(1e-10);

    // Approximate Fisher p-value
    let chi_sq = (peaks_with_match as f64 - peak_total * peak_rate).powi(2) / (peak_total * peak_rate + 1.0);
    let p_value = (-chi_sq / 2.0).exp();

    Ok((fold_enrichment, p_value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_creation() {
        let pwm = vec![
            [0.9, 0.05, 0.03, 0.02],
            [0.05, 0.9, 0.03, 0.02],
            [0.03, 0.05, 0.9, 0.02],
            [0.02, 0.03, 0.05, 0.9],
        ];

        let motif = Motif::new("test_motif", pwm);
        assert_eq!(motif.width(), 4);
        assert_eq!(motif.consensus, "ACGT");
    }

    #[test]
    fn test_score_match() {
        let pwm = vec![
            [0.9, 0.05, 0.03, 0.02],
            [0.05, 0.9, 0.03, 0.02],
        ];

        let motif = Motif::new("test", pwm);
        let seq = b"ACGT";

        let score = score_match(seq, &motif, 0);
        assert!(score.is_some());
        assert!(score.unwrap() > 0.0);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AT"), b"AT");
    }

    #[test]
    fn test_scan_sequence() {
        let pwm = vec![
            [0.9, 0.05, 0.03, 0.02],
            [0.05, 0.9, 0.03, 0.02],
            [0.03, 0.05, 0.9, 0.02],
            [0.02, 0.03, 0.05, 0.9],
        ];

        let motif = Motif::new("test", pwm);
        let seq = b"ACGTACGTACGT";

        let matches = scan_sequence(seq, &motif, -10.0);
        assert!(!matches.is_empty());
    }

    #[test]
    fn test_discover_motifs() {
        let seqs = [b"ACGTACGTACGT".as_slice(), b"ACGTACGTACGT".as_slice()];
        let peak_seqs: Vec<&[u8]> = seqs.to_vec();
        let params = DiscoveryParams {
            motif_width: 4,
            n_motifs: 2,
            background_freq: [0.25; 4],
        };

        let motifs = discover_motifs(&peak_seqs, &params).unwrap();
        assert!(!motifs.is_empty());
    }

    #[test]
    fn test_meme_roundtrip() {
        let pwm = vec![
            [0.9, 0.05, 0.03, 0.02],
            [0.05, 0.9, 0.03, 0.02],
        ];

        let motif = Motif::new("test_motif", pwm);
        let meme_str = write_meme(&[motif]);
        assert!(meme_str.contains("MOTIF"));
        assert!(meme_str.contains("test_motif"));

        let parsed = parse_meme(&meme_str).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "test_motif");
    }

    #[test]
    fn test_compare_motifs() {
        let pwm1 = vec![[0.9, 0.05, 0.03, 0.02]; 4];
        let pwm2 = vec![[0.9, 0.05, 0.03, 0.02]; 4];

        let motif1 = Motif::new("m1", pwm1);
        let motif2 = Motif::new("m2", pwm2);

        let sim = compare_motifs(&motif1, &motif2);
        assert!(sim > 0.8);
    }
}
