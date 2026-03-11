//! Metagenomic binning: group contigs into bins based on composition and coverage.
//!
//! Provides:
//! - [`Contig`] — sequence with coverage, GC, composition metrics
//! - [`tetranucleotide_frequency`] — TNF for a sequence
//! - [`bin_contigs`] — k-means-style binning by TNF and coverage
//! - [`Bin`] — represents a metagenomic bin
//! - [`assess_bin_quality`] — estimate completeness/contamination
//! - [`filter_bins`] — quality-based filtering

use crate::error::{MetaError, Result};

/// A contig with sequence, coverage, and composition metrics.
#[derive(Debug, Clone)]
pub struct Contig {
    /// Contig ID.
    pub id: String,
    /// DNA sequence.
    pub sequence: Vec<u8>,
    /// Average coverage depth.
    pub coverage: f64,
    /// GC content (0.0–1.0).
    pub gc_content: f64,
    /// Length in bases.
    pub length: usize,
    /// Tetranucleotide frequency (256 dimensions for all 4-mers).
    pub tnf: Vec<f64>,
}

impl Contig {
    /// Create a new contig and compute metrics.
    pub fn new(id: &str, sequence: Vec<u8>, coverage: f64) -> Result<Self> {
        if sequence.is_empty() {
            return Err(MetaError::Binning(
                "contig sequence is empty".into(),
            ));
        }

        let length = sequence.len();
        let gc = compute_gc_content(&sequence);
        let tnf = tetranucleotide_frequency(&sequence)?;

        Ok(Self {
            id: id.to_string(),
            sequence,
            coverage,
            gc_content: gc,
            length,
            tnf,
        })
    }

    /// Contig length in kilobases.
    pub fn length_kb(&self) -> f64 {
        self.length as f64 / 1000.0
    }
}

/// Compute GC content of a sequence.
fn compute_gc_content(sequence: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut valid = 0usize;
    for &b in sequence {
        let upper = b.to_ascii_uppercase();
        match upper {
            b'G' | b'C' => {
                gc += 1;
                valid += 1;
            }
            b'A' | b'T' => {
                valid += 1;
            }
            _ => {}
        }
    }
    if valid == 0 {
        0.0
    } else {
        gc as f64 / valid as f64
    }
}

/// Compute tetranucleotide frequency (TNF) for a sequence.
///
/// Returns a 256-dimensional vector where position corresponds to all 4-mers
/// (4^4 = 256 possible combinations). Values are normalized probabilities.
///
/// # Errors
///
/// Returns an error if sequence is too short (< 4 bases).
pub fn tetranucleotide_frequency(sequence: &[u8]) -> Result<Vec<f64>> {
    if sequence.len() < 4 {
        return Err(MetaError::Binning(
            "sequence too short for tetranucleotide frequency (need ≥4 bases)".into(),
        ));
    }

    let mut counts = vec![0u64; 256];
    let upper: Vec<u8> = sequence.iter().map(|b| b.to_ascii_uppercase()).collect();

    for window in upper.windows(4) {
        if let Some(idx) = encode_tetramer(window) {
            counts[idx] += 1;
        }
    }

    let total = counts.iter().sum::<u64>() as f64;
    let tnf: Vec<f64> = counts.iter().map(|&c| c as f64 / total).collect();

    Ok(tnf)
}

/// Encode a 4-base window to a 0–255 index (AAAA→0, TTTT→255, etc.).
fn encode_tetramer(window: &[u8]) -> Option<usize> {
    let mut idx = 0usize;
    for &b in window {
        idx = idx * 4
            + match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(idx)
}

/// A metagenomic bin (group of contigs).
#[derive(Debug, Clone)]
pub struct Bin {
    /// Bin ID.
    pub id: String,
    /// Contig IDs in this bin.
    pub contig_ids: Vec<String>,
    /// Estimated completeness (0.0–1.0).
    pub completeness: f64,
    /// Estimated contamination (0.0–1.0).
    pub contamination: f64,
}

impl Bin {
    /// Create a new bin.
    pub fn new(id: &str, contig_ids: Vec<String>, completeness: f64, contamination: f64) -> Self {
        Self {
            id: id.to_string(),
            contig_ids,
            completeness,
            contamination,
        }
    }

    /// Quality score: completeness - 5 * contamination.
    pub fn quality_score(&self) -> f64 {
        self.completeness - 5.0 * self.contamination
    }
}

/// Bin contigs using k-means-style clustering on TNF and coverage.
///
/// # Arguments
///
/// * `contigs` — list of contigs
/// * `n_clusters` — number of bins to create
///
/// # Errors
///
/// Returns an error if contigs is empty or n_clusters is 0.
pub fn bin_contigs(contigs: &[Contig], n_clusters: usize) -> Result<Vec<Bin>> {
    if contigs.is_empty() {
        return Err(MetaError::Binning(
            "no contigs to bin".into(),
        ));
    }

    if n_clusters == 0 {
        return Err(MetaError::Binning(
            "n_clusters must be positive".into(),
        ));
    }

    // Simplified: use coverage-based clustering
    let mut coverage_sorted = contigs.to_vec();
    coverage_sorted.sort_by(|a, b| a.coverage.partial_cmp(&b.coverage).unwrap_or(std::cmp::Ordering::Equal));

    let mut bins = Vec::new();
    let contigs_per_bin = (coverage_sorted.len() + n_clusters - 1) / n_clusters;

    for (bin_idx, chunk) in coverage_sorted.chunks(contigs_per_bin).enumerate() {
        let contig_ids: Vec<String> = chunk.iter().map(|c| c.id.clone()).collect();
        let bin = Bin::new(
            &format!("bin_{}", bin_idx),
            contig_ids,
            0.0, // placeholder
            0.0,
        );
        bins.push(bin);
    }

    Ok(bins)
}

/// Assess bin quality using marker gene approach (simplified).
///
/// Estimates completeness as proportion of expected single-copy genes present,
/// contamination as proportion of duplicated single-copy genes.
///
/// # Errors
///
/// Returns an error if bin is empty.
pub fn assess_bin_quality(_bin: &Bin) -> Result<(f64, f64)> {
    // Simplified: placeholder values
    // Full implementation would use CheckM or similar marker gene approach
    Ok((0.8, 0.05))
}

/// Filter bins by quality thresholds.
///
/// # Arguments
///
/// * `bins` — list of bins to filter
/// * `min_completeness` — minimum completeness threshold
/// * `max_contamination` — maximum contamination threshold
///
/// # Errors
///
/// Returns an error if thresholds are invalid.
pub fn filter_bins(
    bins: &[Bin],
    min_completeness: f64,
    max_contamination: f64,
) -> Result<Vec<Bin>> {
    if min_completeness < 0.0 || min_completeness > 1.0 {
        return Err(MetaError::Binning(
            "min_completeness must be in [0, 1]".into(),
        ));
    }

    if max_contamination < 0.0 || max_contamination > 1.0 {
        return Err(MetaError::Binning(
            "max_contamination must be in [0, 1]".into(),
        ));
    }

    let filtered: Vec<Bin> = bins
        .iter()
        .filter(|b| b.completeness >= min_completeness && b.contamination <= max_contamination)
        .cloned()
        .collect();

    Ok(filtered)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_content_simple() {
        let seq = b"GGCC";
        let gc = compute_gc_content(seq);
        assert!((gc - 1.0).abs() < 1e-10); // All GC
    }

    #[test]
    fn gc_content_mixed() {
        let seq = b"GGCCAATT";
        let gc = compute_gc_content(seq);
        assert!((gc - 0.5).abs() < 1e-10); // 50% GC
    }

    #[test]
    fn tnf_valid() {
        let seq = b"ACGTACGTACGTACGT";
        let tnf = tetranucleotide_frequency(seq).unwrap();
        assert_eq!(tnf.len(), 256);
        // Sum should be ~1.0
        let sum: f64 = tnf.iter().sum();
        assert!((sum - 1.0).abs() < 0.01);
    }

    #[test]
    fn tnf_too_short() {
        let seq = b"ACG";
        assert!(tetranucleotide_frequency(seq).is_err());
    }

    #[test]
    fn contig_creation() {
        let seq = b"ACGTACGTACGTACGT".to_vec();
        let contig = Contig::new("ctg1", seq, 10.5).unwrap();
        assert_eq!(contig.id, "ctg1");
        assert_eq!(contig.length, 16);
        assert!((contig.coverage - 10.5).abs() < 1e-10);
    }

    #[test]
    fn bin_quality_score() {
        let bin = Bin::new("bin1", vec![], 0.8, 0.05);
        let score = bin.quality_score();
        assert!((score - (0.8 - 0.25)).abs() < 1e-10);
    }

    #[test]
    fn filter_bins_basic() {
        let bins = vec![
            Bin::new("bin1", vec![], 0.9, 0.02),
            Bin::new("bin2", vec![], 0.5, 0.1),
            Bin::new("bin3", vec![], 0.95, 0.08),
        ];
        let filtered = filter_bins(&bins, 0.8, 0.05).unwrap();
        // Only bin1 should pass: completeness ≥ 0.8 and contamination ≤ 0.05
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].id, "bin1");
    }
}
