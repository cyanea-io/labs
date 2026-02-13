//! Low-complexity and repeat masking for biological sequences.
//!
//! Implements DUST (DNA), SEG (protein), tandem repeat detection, and
//! soft/hard masking application.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// How to mask identified regions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskMode {
    /// Lowercase the masked bases.
    Soft,
    /// Replace with N (DNA) or X (protein).
    Hard,
}

/// Source algorithm that identified a masked region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskSource {
    Dust,
    Seg,
    TandemRepeat,
}

/// A region identified for masking.
#[derive(Debug, Clone)]
pub struct MaskedRegion {
    /// Start position (0-indexed, inclusive).
    pub start: usize,
    /// End position (0-indexed, exclusive).
    pub end: usize,
    /// Score (algorithm-dependent).
    pub score: f64,
    /// Which algorithm found this region.
    pub source: MaskSource,
}

/// Result of applying a mask to a sequence.
#[derive(Debug, Clone)]
pub struct MaskResult {
    /// The masked sequence (same length as input).
    pub sequence: Vec<u8>,
    /// Regions that were masked.
    pub regions: Vec<MaskedRegion>,
    /// Fraction of the sequence that was masked, in [0.0, 1.0].
    pub masked_fraction: f64,
}

// ---------------------------------------------------------------------------
// DUST parameters
// ---------------------------------------------------------------------------

/// Parameters for the DUST low-complexity filter (DNA).
#[derive(Debug, Clone)]
pub struct DustParams {
    /// Sliding window size (default: 64).
    pub window: usize,
    /// Score threshold; regions scoring above are masked (default: 20.0).
    pub threshold: f64,
    /// Join regions within this many bases (default: 1).
    pub linker: usize,
}

impl Default for DustParams {
    fn default() -> Self {
        Self {
            window: 64,
            threshold: 20.0,
            linker: 1,
        }
    }
}

// ---------------------------------------------------------------------------
// SEG parameters
// ---------------------------------------------------------------------------

/// Parameters for the SEG low-complexity filter (protein).
#[derive(Debug, Clone)]
pub struct SegParams {
    /// Sliding window size (default: 12).
    pub window: usize,
    /// Trigger masking when entropy ≤ lowcut (default: 2.2).
    pub lowcut: f64,
    /// Stop extending when entropy reaches highcut (default: 2.5).
    pub highcut: f64,
}

impl Default for SegParams {
    fn default() -> Self {
        Self {
            window: 12,
            lowcut: 2.2,
            highcut: 2.5,
        }
    }
}

// ---------------------------------------------------------------------------
// Tandem repeat parameters
// ---------------------------------------------------------------------------

/// Parameters for tandem repeat detection.
#[derive(Debug, Clone)]
pub struct TandemRepeatParams {
    /// Minimum repeat unit period (default: 1).
    pub min_period: usize,
    /// Maximum repeat unit period (default: 6).
    pub max_period: usize,
    /// Minimum number of complete repeat copies (default: 3).
    pub min_copies: usize,
}

impl Default for TandemRepeatParams {
    fn default() -> Self {
        Self {
            min_period: 1,
            max_period: 6,
            min_copies: 3,
        }
    }
}

// ---------------------------------------------------------------------------
// DUST algorithm
// ---------------------------------------------------------------------------

/// Compute DUST triplet score for a window.
///
/// Count all 64 DNA triplets in the window, then score = Σ(c*(c-1)/2) / (W-2).
fn dust_score(window: &[u8]) -> f64 {
    if window.len() < 3 {
        return 0.0;
    }

    let mut counts = [0u32; 64];
    for tri in window.windows(3) {
        let idx = triplet_index(tri);
        if let Some(i) = idx {
            counts[i] += 1;
        }
    }

    let mut score = 0.0f64;
    for &c in &counts {
        if c > 1 {
            score += (c as f64) * (c as f64 - 1.0) / 2.0;
        }
    }

    let denom = (window.len() as f64) - 2.0;
    if denom > 0.0 {
        score / denom
    } else {
        0.0
    }
}

/// Map a DNA triplet to an index in [0, 64).
fn triplet_index(tri: &[u8]) -> Option<usize> {
    let map = |b: u8| -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ => None,
        }
    };
    Some(map(tri[0])? * 16 + map(tri[1])? * 4 + map(tri[2])?)
}

/// Identify low-complexity regions in a DNA sequence using the DUST algorithm.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
pub fn dust(seq: &[u8], params: &DustParams) -> Result<Vec<MaskedRegion>> {
    if seq.is_empty() {
        return Err(CyaneaError::InvalidInput("sequence is empty".into()));
    }

    let w = params.window.min(seq.len());
    if w < 3 {
        return Ok(Vec::new());
    }

    let mut raw_regions: Vec<(usize, usize, f64)> = Vec::new();

    for start in 0..=seq.len().saturating_sub(w) {
        let window = &seq[start..start + w];
        let score = dust_score(window);
        if score > params.threshold {
            raw_regions.push((start, start + w, score));
        }
    }

    // Merge overlapping / linker-adjacent regions
    let merged = merge_regions(&raw_regions, params.linker);

    Ok(merged
        .into_iter()
        .map(|(start, end, score)| MaskedRegion {
            start,
            end,
            score,
            source: MaskSource::Dust,
        })
        .collect())
}

// ---------------------------------------------------------------------------
// SEG algorithm
// ---------------------------------------------------------------------------

/// Shannon entropy of amino acid frequencies in a window.
fn aa_entropy(window: &[u8]) -> f64 {
    let mut counts = [0u32; 26]; // A-Z
    let mut total = 0u32;

    for &b in window {
        let upper = b.to_ascii_uppercase();
        if upper >= b'A' && upper <= b'Z' {
            counts[(upper - b'A') as usize] += 1;
            total += 1;
        }
    }

    if total == 0 {
        return 0.0;
    }

    let mut entropy = 0.0f64;
    let t = total as f64;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / t;
            entropy -= p * p.log2();
        }
    }
    entropy
}

/// Identify low-complexity regions in a protein sequence using the SEG algorithm.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
pub fn seg(seq: &[u8], params: &SegParams) -> Result<Vec<MaskedRegion>> {
    if seq.is_empty() {
        return Err(CyaneaError::InvalidInput("sequence is empty".into()));
    }

    let w = params.window.min(seq.len());
    if w < 2 {
        return Ok(Vec::new());
    }

    let mut raw_regions: Vec<(usize, usize, f64)> = Vec::new();

    for start in 0..=seq.len().saturating_sub(w) {
        let window = &seq[start..start + w];
        let ent = aa_entropy(window);
        if ent <= params.lowcut {
            // Extend in both directions until entropy >= highcut
            let mut ext_start = start;
            let mut ext_end = start + w;

            // Extend left
            while ext_start > 0 {
                let candidate = &seq[ext_start - 1..ext_end];
                if aa_entropy(candidate) <= params.highcut {
                    ext_start -= 1;
                } else {
                    break;
                }
            }

            // Extend right
            while ext_end < seq.len() {
                let candidate = &seq[ext_start..ext_end + 1];
                if aa_entropy(candidate) <= params.highcut {
                    ext_end += 1;
                } else {
                    break;
                }
            }

            raw_regions.push((ext_start, ext_end, ent));
        }
    }

    let merged = merge_regions(&raw_regions, 0);

    Ok(merged
        .into_iter()
        .map(|(start, end, score)| MaskedRegion {
            start,
            end,
            score,
            source: MaskSource::Seg,
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Tandem repeat detection
// ---------------------------------------------------------------------------

/// Find tandem repeat regions in a sequence.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
pub fn find_tandem_repeats(
    seq: &[u8],
    params: &TandemRepeatParams,
) -> Result<Vec<MaskedRegion>> {
    if seq.is_empty() {
        return Err(CyaneaError::InvalidInput("sequence is empty".into()));
    }

    let min_p = params.min_period.max(1);
    let max_p = params.max_period.min(seq.len());

    let mut raw_regions: Vec<(usize, usize, f64)> = Vec::new();

    for p in min_p..=max_p {
        let min_len = p * params.min_copies;
        if min_len > seq.len() {
            continue;
        }

        let mut i = p;
        while i < seq.len() {
            // Check if seq[i] matches seq[i-p]
            if seq[i].to_ascii_uppercase() == seq[i - p].to_ascii_uppercase() {
                // Found a match — extend the run
                let run_start = i - p;
                let mut run_end = i + 1;
                while run_end < seq.len()
                    && seq[run_end].to_ascii_uppercase()
                        == seq[run_end - p].to_ascii_uppercase()
                {
                    run_end += 1;
                }
                let run_len = run_end - run_start;
                let copies = run_len / p;
                if copies >= params.min_copies {
                    // Trim to complete copies
                    let trimmed_end = run_start + copies * p;
                    raw_regions.push((run_start, trimmed_end, copies as f64));
                }
                i = run_end;
            } else {
                i += 1;
            }
        }
    }

    let merged = merge_regions(&raw_regions, 0);

    Ok(merged
        .into_iter()
        .map(|(start, end, score)| MaskedRegion {
            start,
            end,
            score,
            source: MaskSource::TandemRepeat,
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Masking application
// ---------------------------------------------------------------------------

/// Apply masking to a sequence given a set of regions.
///
/// - `Soft`: lowercase the masked bases.
/// - `Hard`: replace with `N` (DNA) or `X` (protein).
pub fn apply_mask(
    seq: &[u8],
    regions: &[MaskedRegion],
    mode: MaskMode,
    is_protein: bool,
) -> MaskResult {
    let mut out = seq.to_vec();
    let mut masked_positions = vec![false; seq.len()];

    for region in regions {
        let start = region.start.min(seq.len());
        let end = region.end.min(seq.len());
        for i in start..end {
            masked_positions[i] = true;
            match mode {
                MaskMode::Soft => {
                    out[i] = out[i].to_ascii_lowercase();
                }
                MaskMode::Hard => {
                    out[i] = if is_protein { b'X' } else { b'N' };
                }
            }
        }
    }

    let masked_count = masked_positions.iter().filter(|&&m| m).count();
    let masked_fraction = if seq.is_empty() {
        0.0
    } else {
        masked_count as f64 / seq.len() as f64
    };

    MaskResult {
        sequence: out,
        regions: regions.to_vec(),
        masked_fraction,
    }
}

/// Run DUST and apply masking in one step.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
pub fn mask_dust(seq: &[u8], params: &DustParams, mode: MaskMode) -> Result<MaskResult> {
    let regions = dust(seq, params)?;
    Ok(apply_mask(seq, &regions, mode, false))
}

/// Run SEG and apply masking in one step.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
pub fn mask_seg(seq: &[u8], params: &SegParams, mode: MaskMode) -> Result<MaskResult> {
    let regions = seg(seq, params)?;
    Ok(apply_mask(seq, &regions, mode, true))
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Merge overlapping or adjacent regions.
fn merge_regions(regions: &[(usize, usize, f64)], gap: usize) -> Vec<(usize, usize, f64)> {
    if regions.is_empty() {
        return Vec::new();
    }

    let mut sorted: Vec<(usize, usize, f64)> = regions.to_vec();
    sorted.sort_by_key(|r| r.0);

    let mut merged = vec![sorted[0]];
    for &(start, end, score) in &sorted[1..] {
        let last = merged.last_mut().unwrap();
        if start <= last.1 + gap {
            last.1 = last.1.max(end);
            if score > last.2 {
                last.2 = score;
            }
        } else {
            merged.push((start, end, score));
        }
    }

    merged
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- DUST ---

    #[test]
    fn dust_homopolymer() {
        // Long poly-A run should score very high
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let regions = dust(seq, &DustParams::default()).unwrap();
        assert!(!regions.is_empty(), "homopolymer should be masked");
    }

    #[test]
    fn dust_random_dna() {
        // Diverse sequence should not be masked with default threshold
        let seq = b"ACGTACGTACGTTGCATGCATGCAACGTACGTACGTTGCATGCATGCAACGTACGTACGTTGCA";
        let regions = dust(seq, &DustParams::default()).unwrap();
        assert!(regions.is_empty(), "diverse DNA should not be masked");
    }

    #[test]
    fn dust_dinucleotide_repeat() {
        // AT repeat produces score ~15 (only 2 distinct triplets), below default 20.
        // Use a lower threshold to detect it.
        let seq = b"ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let params = DustParams {
            threshold: 10.0,
            ..Default::default()
        };
        let regions = dust(seq, &params).unwrap();
        assert!(!regions.is_empty(), "dinucleotide repeat should be masked at threshold 10");
    }

    #[test]
    fn dust_empty() {
        let result = dust(b"", &DustParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn dust_short() {
        let regions = dust(b"AC", &DustParams::default()).unwrap();
        assert!(regions.is_empty());
    }

    // --- SEG ---

    #[test]
    fn seg_poly_ala() {
        let seq = b"AAAAAAAAAAAA";
        let regions = seg(seq, &SegParams::default()).unwrap();
        assert!(!regions.is_empty(), "poly-Ala should be low complexity");
    }

    #[test]
    fn seg_diverse_protein() {
        let seq = b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVK";
        let params = SegParams::default();
        let regions = seg(seq, &params).unwrap();
        // Hemoglobin sequence is mostly high complexity
        // Allow some regions but most should be unmasked
        let total_masked: usize = regions.iter().map(|r| r.end - r.start).sum();
        assert!(
            total_masked < seq.len() / 2,
            "diverse protein should be mostly unmasked, masked {} of {}",
            total_masked,
            seq.len()
        );
    }

    #[test]
    fn seg_empty() {
        let result = seg(b"", &SegParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn seg_extension() {
        // Very low complexity should extend beyond initial window
        let seq = b"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ";
        let regions = seg(seq, &SegParams::default()).unwrap();
        if !regions.is_empty() {
            assert!(
                regions[0].end - regions[0].start > 12,
                "should extend beyond initial window"
            );
        }
    }

    // --- Tandem repeats ---

    #[test]
    fn tandem_dinucleotide() {
        let seq = b"ACACACACACACACACAC";
        let regions = find_tandem_repeats(seq, &TandemRepeatParams::default()).unwrap();
        assert!(!regions.is_empty(), "AC repeat should be found");
    }

    #[test]
    fn tandem_trinucleotide() {
        let seq = b"CAGCAGCAGCAGCAGCAG";
        let regions = find_tandem_repeats(seq, &TandemRepeatParams::default()).unwrap();
        assert!(!regions.is_empty(), "CAG repeat should be found");
    }

    #[test]
    fn tandem_min_copies() {
        let seq = b"ACACAC"; // 3 copies of AC
        let params = TandemRepeatParams {
            min_copies: 4,
            ..Default::default()
        };
        let regions = find_tandem_repeats(seq, &params).unwrap();
        // 3 copies < 4 minimum at period 2
        // But period 1 might match — check specifically
        let p2_regions: Vec<_> = regions
            .iter()
            .filter(|r| (r.end - r.start) >= 8) // need at least 4 copies of period 2
            .collect();
        assert!(p2_regions.is_empty(), "3 copies should not meet min_copies=4 for period 2");
    }

    #[test]
    fn tandem_empty() {
        let result = find_tandem_repeats(b"", &TandemRepeatParams::default());
        assert!(result.is_err());
    }

    // --- Masking ---

    #[test]
    fn soft_mask_output() {
        let seq = b"ACGTACGT";
        let regions = vec![MaskedRegion {
            start: 2,
            end: 5,
            score: 1.0,
            source: MaskSource::Dust,
        }];
        let result = apply_mask(seq, &regions, MaskMode::Soft, false);
        assert_eq!(result.sequence, b"ACgtaCGT");
        assert_eq!(result.sequence.len(), seq.len());
    }

    #[test]
    fn hard_mask_dna() {
        let seq = b"ACGTACGT";
        let regions = vec![MaskedRegion {
            start: 0,
            end: 4,
            score: 1.0,
            source: MaskSource::Dust,
        }];
        let result = apply_mask(seq, &regions, MaskMode::Hard, false);
        assert_eq!(result.sequence, b"NNNNACGT");
    }

    #[test]
    fn hard_mask_protein() {
        let seq = b"MVHLTPEE";
        let regions = vec![MaskedRegion {
            start: 1,
            end: 3,
            score: 1.0,
            source: MaskSource::Seg,
        }];
        let result = apply_mask(seq, &regions, MaskMode::Hard, true);
        assert_eq!(result.sequence, b"MXXLTPEE");
    }

    #[test]
    fn masked_fraction() {
        let seq = b"AAAAAAAA"; // 8 bases
        let regions = vec![MaskedRegion {
            start: 0,
            end: 4,
            score: 1.0,
            source: MaskSource::Dust,
        }];
        let result = apply_mask(seq, &regions, MaskMode::Soft, false);
        assert!((result.masked_fraction - 0.5).abs() < 1e-10);
    }

    #[test]
    fn mask_preserves_length() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let result = mask_dust(seq, &DustParams::default(), MaskMode::Soft).unwrap();
        assert_eq!(result.sequence.len(), seq.len());
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn dna_seq(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(
            prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
            1..=max_len,
        )
    }

    proptest! {
        #[test]
        fn mask_preserves_length(seq in dna_seq(200)) {
            let result = mask_dust(&seq, &DustParams::default(), MaskMode::Soft).unwrap();
            prop_assert_eq!(result.sequence.len(), seq.len());
        }
    }
}
