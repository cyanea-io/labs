//! BEDTools-style genome arithmetic on genomic intervals.
//!
//! All functions operate on slices of [`GenomicInterval`] for flexibility —
//! works with raw vecs, BED parse output, or [`IntervalSet::intervals()`].
//! Coordinates are 0-based, half-open `[start, end)`.

use std::collections::BTreeMap;

use cyanea_core::{CyaneaError, Result};

use crate::genomic::{GenomicInterval, Strand};

/// Chromosome name → size mapping.
pub type GenomeInfo = BTreeMap<String, u64>;

/// How to handle strand when comparing intervals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandMode {
    /// Ignore strand (default BEDTools behavior).
    Ignore,
    /// Only match intervals on the same strand (`-s`).
    Same,
    /// Only match intervals on the opposite strand (`-S`).
    Opposite,
}

/// Result of a closest-interval query.
#[derive(Debug, Clone, PartialEq)]
pub struct ClosestResult {
    /// The query interval.
    pub query: GenomicInterval,
    /// The closest interval in the target set, if any.
    pub closest: Option<GenomicInterval>,
    /// Distance to the closest interval (0 if overlapping, `None` if no target on chrom).
    pub distance: Option<u64>,
}

/// Extended Jaccard similarity statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct JaccardStats {
    /// Total base pairs in the intersection.
    pub intersection_bp: u64,
    /// Total base pairs in the union.
    pub union_bp: u64,
    /// Jaccard index (intersection_bp / union_bp).
    pub jaccard: f64,
    /// Number of intersecting interval pairs.
    pub n_intersections: u64,
}

/// Convenience constructor for [`GenomeInfo`].
pub fn genome_info(chroms: &[(&str, u64)]) -> GenomeInfo {
    chroms.iter().map(|(c, s)| (c.to_string(), *s)).collect()
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Group intervals by chromosome, sorting each group by start position.
fn group_by_chrom(intervals: &[GenomicInterval]) -> BTreeMap<String, Vec<&GenomicInterval>> {
    let mut groups: BTreeMap<String, Vec<&GenomicInterval>> = BTreeMap::new();
    for iv in intervals {
        groups.entry(iv.chrom.clone()).or_default().push(iv);
    }
    for group in groups.values_mut() {
        group.sort_by_key(|iv| (iv.start, iv.end));
    }
    groups
}

/// Check whether two intervals match under the given strand mode.
fn strands_match(a: &GenomicInterval, b: &GenomicInterval, mode: StrandMode) -> bool {
    match mode {
        StrandMode::Ignore => true,
        StrandMode::Same => a.strand == b.strand,
        StrandMode::Opposite => {
            matches!(
                (&a.strand, &b.strand),
                (Strand::Forward, Strand::Reverse) | (Strand::Reverse, Strand::Forward)
            )
        }
    }
}

/// Merge a slice of intervals that are already sorted by start on a single chromosome.
/// Abutting intervals (end == next.start) are merged.
fn merge_sorted(sorted: &[&GenomicInterval]) -> Vec<GenomicInterval> {
    if sorted.is_empty() {
        return Vec::new();
    }
    let mut merged = Vec::new();
    let mut cur_start = sorted[0].start;
    let mut cur_end = sorted[0].end;
    let chrom = &sorted[0].chrom;
    let strand = sorted[0].strand;

    for iv in &sorted[1..] {
        if iv.start <= cur_end {
            cur_end = cur_end.max(iv.end);
        } else {
            merged.push(GenomicInterval {
                chrom: chrom.clone(),
                start: cur_start,
                end: cur_end,
                strand,
            });
            cur_start = iv.start;
            cur_end = iv.end;
        }
    }
    merged.push(GenomicInterval {
        chrom: chrom.clone(),
        start: cur_start,
        end: cur_end,
        strand,
    });
    merged
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Merge overlapping and abutting intervals.
///
/// When `strand_mode` is [`StrandMode::Same`], intervals are grouped by
/// (chromosome, strand) before merging. Otherwise strand is ignored.
pub fn merge(intervals: &[GenomicInterval], strand_mode: StrandMode) -> Vec<GenomicInterval> {
    if intervals.is_empty() {
        return Vec::new();
    }

    if strand_mode == StrandMode::Same {
        // Group by (chrom, strand)
        let mut groups: BTreeMap<(String, Strand), Vec<&GenomicInterval>> = BTreeMap::new();
        for iv in intervals {
            groups
                .entry((iv.chrom.clone(), iv.strand))
                .or_default()
                .push(iv);
        }
        let mut result = Vec::new();
        for group in groups.values_mut() {
            group.sort_by_key(|iv| (iv.start, iv.end));
            result.extend(merge_sorted(group));
        }
        result
    } else {
        let groups = group_by_chrom(intervals);
        let mut result = Vec::new();
        for group in groups.values() {
            result.extend(merge_sorted(group));
        }
        result
    }
}

/// Union of two interval sets (concatenate + merge).
pub fn union(
    a: &[GenomicInterval],
    b: &[GenomicInterval],
    strand_mode: StrandMode,
) -> Vec<GenomicInterval> {
    let mut combined: Vec<GenomicInterval> = Vec::with_capacity(a.len() + b.len());
    combined.extend_from_slice(a);
    combined.extend_from_slice(b);
    merge(&combined, strand_mode)
}

/// Intersect two interval sets, returning the overlapping sub-regions.
///
/// For each pair of overlapping intervals (one from `a`, one from `b`),
/// emits `max(a.start, b.start)..min(a.end, b.end)`.
pub fn intersect(
    a: &[GenomicInterval],
    b: &[GenomicInterval],
    strand_mode: StrandMode,
) -> Result<Vec<GenomicInterval>> {
    let groups_a = group_by_chrom(a);
    let groups_b = group_by_chrom(b);
    let mut result = Vec::new();

    for (chrom, a_ivs) in &groups_a {
        let b_ivs = match groups_b.get(chrom) {
            Some(v) => v,
            None => continue,
        };

        for a_iv in a_ivs {
            for b_iv in b_ivs {
                if !strands_match(a_iv, b_iv, strand_mode) {
                    continue;
                }
                if a_iv.start < b_iv.end && b_iv.start < a_iv.end {
                    let start = a_iv.start.max(b_iv.start);
                    let end = a_iv.end.min(b_iv.end);
                    result.push(GenomicInterval {
                        chrom: chrom.clone(),
                        start,
                        end,
                        strand: a_iv.strand,
                    });
                }
            }
        }
    }

    Ok(result)
}

/// Return intervals from `a` that have any overlap with `b` (BEDTools `-u`).
pub fn intersect_report_a(
    a: &[GenomicInterval],
    b: &[GenomicInterval],
    strand_mode: StrandMode,
) -> Vec<GenomicInterval> {
    let groups_b = group_by_chrom(b);
    let mut result = Vec::new();

    for a_iv in a {
        let b_ivs = match groups_b.get(&a_iv.chrom) {
            Some(v) => v,
            None => continue,
        };
        let has_overlap = b_ivs.iter().any(|b_iv| {
            strands_match(a_iv, b_iv, strand_mode)
                && a_iv.start < b_iv.end
                && b_iv.start < a_iv.end
        });
        if has_overlap {
            result.push(a_iv.clone());
        }
    }

    result
}

/// Subtract intervals in `b` from intervals in `a`.
///
/// For each interval in `a`, removes any bases covered by (merged) `b`.
pub fn subtract(
    a: &[GenomicInterval],
    b: &[GenomicInterval],
    strand_mode: StrandMode,
) -> Result<Vec<GenomicInterval>> {
    if b.is_empty() {
        return Ok(a.to_vec());
    }

    let merged_b = merge(b, StrandMode::Ignore);
    let groups_b = group_by_chrom(&merged_b);
    let mut result = Vec::new();

    for a_iv in a {
        let b_ivs = match groups_b.get(&a_iv.chrom) {
            Some(v) => v,
            None => {
                result.push(a_iv.clone());
                continue;
            }
        };

        // Filter b intervals by strand
        let matching_b: Vec<&&GenomicInterval> = b_ivs
            .iter()
            .filter(|b_iv| strands_match(a_iv, b_iv, strand_mode))
            .collect();

        if matching_b.is_empty() {
            result.push(a_iv.clone());
            continue;
        }

        // Walk through matching b intervals, emitting uncovered regions of a
        let mut pos = a_iv.start;

        for b_iv in &matching_b {
            if b_iv.start >= a_iv.end {
                break;
            }
            if b_iv.end <= pos {
                continue;
            }
            // Emit gap before this b interval
            if b_iv.start > pos {
                let end = b_iv.start.min(a_iv.end);
                if end > pos {
                    result.push(GenomicInterval {
                        chrom: a_iv.chrom.clone(),
                        start: pos,
                        end,
                        strand: a_iv.strand,
                    });
                }
            }
            pos = pos.max(b_iv.end);
        }

        // Emit trailing region after last b overlap
        if pos < a_iv.end {
            result.push(GenomicInterval {
                chrom: a_iv.chrom.clone(),
                start: pos,
                end: a_iv.end,
                strand: a_iv.strand,
            });
        }
    }

    Ok(result)
}

/// Complement of intervals with respect to a genome.
///
/// Returns the gaps between (merged) intervals and chromosome boundaries `[0, size)`.
pub fn complement(
    intervals: &[GenomicInterval],
    genome: &GenomeInfo,
) -> Result<Vec<GenomicInterval>> {
    // Validate all intervals reference known chroms and are within bounds
    for iv in intervals {
        match genome.get(&iv.chrom) {
            None => {
                return Err(CyaneaError::InvalidInput(format!(
                    "chromosome '{}' not found in genome info",
                    iv.chrom
                )));
            }
            Some(&size) => {
                if iv.end > size {
                    return Err(CyaneaError::InvalidInput(format!(
                        "interval {}:{}-{} exceeds chromosome size {}",
                        iv.chrom, iv.start, iv.end, size
                    )));
                }
            }
        }
    }

    let merged = merge(intervals, StrandMode::Ignore);
    let groups = group_by_chrom(&merged);
    let mut result = Vec::new();

    for (chrom, size) in genome {
        let mut pos = 0u64;
        if let Some(ivs) = groups.get(chrom) {
            for iv in ivs {
                if iv.start > pos {
                    result.push(GenomicInterval {
                        chrom: chrom.clone(),
                        start: pos,
                        end: iv.start,
                        strand: Strand::Unknown,
                    });
                }
                pos = iv.end;
            }
        }
        if pos < *size {
            result.push(GenomicInterval {
                chrom: chrom.clone(),
                start: pos,
                end: *size,
                strand: Strand::Unknown,
            });
        }
    }

    Ok(result)
}

/// Find the closest interval in `b` for each interval in `a`.
pub fn closest(
    a: &[GenomicInterval],
    b: &[GenomicInterval],
    strand_mode: StrandMode,
) -> Vec<ClosestResult> {
    let groups_b = group_by_chrom(b);
    let mut results = Vec::with_capacity(a.len());

    for a_iv in a {
        let b_ivs = match groups_b.get(&a_iv.chrom) {
            Some(v) => v,
            None => {
                results.push(ClosestResult {
                    query: a_iv.clone(),
                    closest: None,
                    distance: None,
                });
                continue;
            }
        };

        // Filter by strand
        let matching: Vec<&GenomicInterval> = b_ivs
            .iter()
            .filter(|b_iv| strands_match(a_iv, b_iv, strand_mode))
            .copied()
            .collect();

        if matching.is_empty() {
            results.push(ClosestResult {
                query: a_iv.clone(),
                closest: None,
                distance: None,
            });
            continue;
        }

        // Binary search for the position where a_iv.start would insert
        let idx = matching.partition_point(|iv| iv.end <= a_iv.start);

        let mut best: Option<(&GenomicInterval, u64)> = None;

        // Check candidates around the insertion point
        for &candidate_idx in &[
            idx.wrapping_sub(1),
            idx,
            idx + 1,
        ] {
            if candidate_idx >= matching.len() {
                continue;
            }
            let b_iv = matching[candidate_idx];
            let dist = if a_iv.start < b_iv.end && b_iv.start < a_iv.end {
                0 // overlapping
            } else if a_iv.end <= b_iv.start {
                b_iv.start - a_iv.end
            } else {
                a_iv.start - b_iv.end
            };

            if best.is_none() || dist < best.unwrap().1 {
                best = Some((b_iv, dist));
            }
        }

        match best {
            Some((b_iv, dist)) => {
                results.push(ClosestResult {
                    query: a_iv.clone(),
                    closest: Some(b_iv.clone()),
                    distance: Some(dist),
                });
            }
            None => {
                results.push(ClosestResult {
                    query: a_iv.clone(),
                    closest: None,
                    distance: None,
                });
            }
        }
    }

    results
}

/// Generate non-overlapping tiling windows across a genome.
pub fn make_windows(genome: &GenomeInfo, window_size: u64) -> Result<Vec<GenomicInterval>> {
    if window_size == 0 {
        return Err(CyaneaError::InvalidInput(
            "window size must be > 0".into(),
        ));
    }

    let mut result = Vec::new();
    for (chrom, &size) in genome {
        let mut start = 0u64;
        while start < size {
            let end = (start + window_size).min(size);
            result.push(GenomicInterval {
                chrom: chrom.clone(),
                start,
                end,
                strand: Strand::Unknown,
            });
            start += window_size;
        }
    }

    Ok(result)
}

/// Generate sliding windows across a genome.
pub fn make_sliding_windows(
    genome: &GenomeInfo,
    window_size: u64,
    step: u64,
) -> Result<Vec<GenomicInterval>> {
    if window_size == 0 {
        return Err(CyaneaError::InvalidInput(
            "window size must be > 0".into(),
        ));
    }
    if step == 0 {
        return Err(CyaneaError::InvalidInput("step must be > 0".into()));
    }

    let mut result = Vec::new();
    for (chrom, &size) in genome {
        let mut start = 0u64;
        while start < size {
            let end = (start + window_size).min(size);
            result.push(GenomicInterval {
                chrom: chrom.clone(),
                start,
                end,
                strand: Strand::Unknown,
            });
            start += step;
        }
    }

    Ok(result)
}

/// Generate windows centered on each interval's midpoint, clipped to chromosome bounds.
pub fn windows_around(
    intervals: &[GenomicInterval],
    window_size: u64,
    genome: &GenomeInfo,
) -> Result<Vec<GenomicInterval>> {
    if window_size == 0 {
        return Err(CyaneaError::InvalidInput(
            "window size must be > 0".into(),
        ));
    }

    let half = window_size / 2;
    let mut result = Vec::with_capacity(intervals.len());

    for iv in intervals {
        let chrom_size = genome.get(&iv.chrom).ok_or_else(|| {
            CyaneaError::InvalidInput(format!(
                "chromosome '{}' not found in genome info",
                iv.chrom
            ))
        })?;

        let mid = iv.midpoint();
        let start = mid.saturating_sub(half);
        let end = (mid + half).min(*chrom_size);

        if start < end {
            result.push(GenomicInterval {
                chrom: iv.chrom.clone(),
                start,
                end,
                strand: iv.strand,
            });
        }
    }

    Ok(result)
}

/// Jaccard similarity between two interval sets.
///
/// Returns `intersection_bp / union_bp`, or 0.0 if both sets are empty.
pub fn jaccard(a: &[GenomicInterval], b: &[GenomicInterval]) -> f64 {
    jaccard_stats(a, b).jaccard
}

/// Extended Jaccard similarity statistics.
pub fn jaccard_stats(a: &[GenomicInterval], b: &[GenomicInterval]) -> JaccardStats {
    let merged_a = merge(a, StrandMode::Ignore);
    let merged_b = merge(b, StrandMode::Ignore);

    // Compute intersection bp
    let isect = intersect(&merged_a, &merged_b, StrandMode::Ignore).unwrap_or_default();
    let intersection_bp: u64 = isect.iter().map(|iv| iv.len()).sum();
    let n_intersections = isect.len() as u64;

    // Compute union bp
    let all_union = union(&merged_a, &merged_b, StrandMode::Ignore);
    let union_bp: u64 = all_union.iter().map(|iv| iv.len()).sum();

    let jaccard = if union_bp == 0 {
        0.0
    } else {
        intersection_bp as f64 / union_bp as f64
    };

    JaccardStats {
        intersection_bp,
        union_bp,
        jaccard,
        n_intersections,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> GenomicInterval {
        GenomicInterval::new(chrom, start, end).unwrap()
    }

    fn iv_strand(chrom: &str, start: u64, end: u64, strand: Strand) -> GenomicInterval {
        GenomicInterval::with_strand(chrom, start, end, strand).unwrap()
    }

    // -----------------------------------------------------------------------
    // merge
    // -----------------------------------------------------------------------

    #[test]
    fn test_merge_overlapping() {
        let intervals = vec![iv("chr1", 100, 200), iv("chr1", 150, 300)];
        let merged = merge(&intervals, StrandMode::Ignore);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 300);
    }

    #[test]
    fn test_merge_abutting() {
        let intervals = vec![iv("chr1", 100, 200), iv("chr1", 200, 300)];
        let merged = merge(&intervals, StrandMode::Ignore);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 300);
    }

    #[test]
    fn test_merge_disjoint() {
        let intervals = vec![iv("chr1", 100, 200), iv("chr1", 300, 400)];
        let merged = merge(&intervals, StrandMode::Ignore);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn test_merge_empty() {
        let merged = merge(&[], StrandMode::Ignore);
        assert!(merged.is_empty());
    }

    #[test]
    fn test_merge_multi_chrom() {
        let intervals = vec![
            iv("chr2", 100, 200),
            iv("chr1", 100, 200),
            iv("chr1", 150, 300),
        ];
        let merged = merge(&intervals, StrandMode::Ignore);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].chrom, "chr1");
        assert_eq!(merged[0].end, 300);
        assert_eq!(merged[1].chrom, "chr2");
    }

    #[test]
    fn test_merge_strand_aware() {
        let intervals = vec![
            iv_strand("chr1", 100, 200, Strand::Forward),
            iv_strand("chr1", 150, 300, Strand::Reverse),
        ];
        let merged = merge(&intervals, StrandMode::Same);
        assert_eq!(merged.len(), 2); // different strands, not merged
    }

    #[test]
    fn test_merge_strand_aware_same() {
        let intervals = vec![
            iv_strand("chr1", 100, 200, Strand::Forward),
            iv_strand("chr1", 150, 300, Strand::Forward),
        ];
        let merged = merge(&intervals, StrandMode::Same);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].end, 300);
    }

    // -----------------------------------------------------------------------
    // union
    // -----------------------------------------------------------------------

    #[test]
    fn test_union_basic() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 150, 300)];
        let u = union(&a, &b, StrandMode::Ignore);
        assert_eq!(u.len(), 1);
        assert_eq!(u[0].start, 100);
        assert_eq!(u[0].end, 300);
    }

    #[test]
    fn test_union_disjoint() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 400, 500)];
        let u = union(&a, &b, StrandMode::Ignore);
        assert_eq!(u.len(), 2);
    }

    // -----------------------------------------------------------------------
    // intersect
    // -----------------------------------------------------------------------

    #[test]
    fn test_intersect_basic() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 400)];
        let r = intersect(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 200);
        assert_eq!(r[0].end, 300);
    }

    #[test]
    fn test_intersect_no_overlap() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 300, 400)];
        let r = intersect(&a, &b, StrandMode::Ignore).unwrap();
        assert!(r.is_empty());
    }

    #[test]
    fn test_intersect_multiple_overlaps() {
        let a = vec![iv("chr1", 100, 500)];
        let b = vec![iv("chr1", 150, 200), iv("chr1", 300, 400)];
        let r = intersect(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].start, 150);
        assert_eq!(r[0].end, 200);
        assert_eq!(r[1].start, 300);
        assert_eq!(r[1].end, 400);
    }

    #[test]
    fn test_intersect_multi_chrom() {
        let a = vec![iv("chr1", 100, 200), iv("chr2", 100, 200)];
        let b = vec![iv("chr1", 150, 250)];
        let r = intersect(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].chrom, "chr1");
    }

    #[test]
    fn test_intersect_empty_input() {
        let a: Vec<GenomicInterval> = vec![];
        let b = vec![iv("chr1", 100, 200)];
        let r = intersect(&a, &b, StrandMode::Ignore).unwrap();
        assert!(r.is_empty());
    }

    #[test]
    fn test_intersect_identical() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 100, 200)];
        let r = intersect(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 100);
        assert_eq!(r[0].end, 200);
    }

    #[test]
    fn test_intersect_strand_same() {
        let a = vec![iv_strand("chr1", 100, 300, Strand::Forward)];
        let b = vec![iv_strand("chr1", 200, 400, Strand::Reverse)];
        let r = intersect(&a, &b, StrandMode::Same).unwrap();
        assert!(r.is_empty()); // different strands
    }

    #[test]
    fn test_intersect_strand_opposite() {
        let a = vec![iv_strand("chr1", 100, 300, Strand::Forward)];
        let b = vec![iv_strand("chr1", 200, 400, Strand::Reverse)];
        let r = intersect(&a, &b, StrandMode::Opposite).unwrap();
        assert_eq!(r.len(), 1);
    }

    // -----------------------------------------------------------------------
    // intersect_report_a
    // -----------------------------------------------------------------------

    #[test]
    fn test_intersect_report_a_basic() {
        let a = vec![iv("chr1", 100, 200), iv("chr1", 400, 500)];
        let b = vec![iv("chr1", 150, 250)];
        let r = intersect_report_a(&a, &b, StrandMode::Ignore);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 100);
        assert_eq!(r[0].end, 200);
    }

    // -----------------------------------------------------------------------
    // subtract
    // -----------------------------------------------------------------------

    #[test]
    fn test_subtract_no_overlap() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 300, 400)];
        let r = subtract(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 100);
        assert_eq!(r[0].end, 200);
    }

    #[test]
    fn test_subtract_complete_removal() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 50, 250)];
        let r = subtract(&a, &b, StrandMode::Ignore).unwrap();
        assert!(r.is_empty());
    }

    #[test]
    fn test_subtract_left_trim() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 50, 200)];
        let r = subtract(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 200);
        assert_eq!(r[0].end, 300);
    }

    #[test]
    fn test_subtract_right_trim() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 400)];
        let r = subtract(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 100);
        assert_eq!(r[0].end, 200);
    }

    #[test]
    fn test_subtract_split_middle() {
        let a = vec![iv("chr1", 100, 400)];
        let b = vec![iv("chr1", 200, 300)];
        let r = subtract(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].start, 100);
        assert_eq!(r[0].end, 200);
        assert_eq!(r[1].start, 300);
        assert_eq!(r[1].end, 400);
    }

    #[test]
    fn test_subtract_multiple_b() {
        let a = vec![iv("chr1", 100, 500)];
        let b = vec![iv("chr1", 150, 200), iv("chr1", 300, 350)];
        let r = subtract(&a, &b, StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 3);
        assert_eq!((r[0].start, r[0].end), (100, 150));
        assert_eq!((r[1].start, r[1].end), (200, 300));
        assert_eq!((r[2].start, r[2].end), (350, 500));
    }

    #[test]
    fn test_subtract_empty_b() {
        let a = vec![iv("chr1", 100, 200)];
        let r = subtract(&a, &[], StrandMode::Ignore).unwrap();
        assert_eq!(r.len(), 1);
    }

    #[test]
    fn test_subtract_strand_mode() {
        let a = vec![iv_strand("chr1", 100, 300, Strand::Forward)];
        let b = vec![iv_strand("chr1", 150, 250, Strand::Reverse)];
        // Same strand mode: b doesn't match, a passes through
        let r = subtract(&a, &b, StrandMode::Same).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start, 100);
        assert_eq!(r[0].end, 300);
    }

    // -----------------------------------------------------------------------
    // complement
    // -----------------------------------------------------------------------

    #[test]
    fn test_complement_basic() {
        let intervals = vec![iv("chr1", 100, 200), iv("chr1", 300, 400)];
        let genome = genome_info(&[("chr1", 500)]);
        let r = complement(&intervals, &genome).unwrap();
        assert_eq!(r.len(), 3);
        assert_eq!((r[0].start, r[0].end), (0, 100));
        assert_eq!((r[1].start, r[1].end), (200, 300));
        assert_eq!((r[2].start, r[2].end), (400, 500));
    }

    #[test]
    fn test_complement_full_coverage() {
        let intervals = vec![iv("chr1", 0, 500)];
        let genome = genome_info(&[("chr1", 500)]);
        let r = complement(&intervals, &genome).unwrap();
        assert!(r.is_empty());
    }

    #[test]
    fn test_complement_empty_input() {
        let genome = genome_info(&[("chr1", 1000)]);
        let r = complement(&[], &genome).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].start, r[0].end), (0, 1000));
    }

    #[test]
    fn test_complement_chrom_not_in_genome() {
        let intervals = vec![iv("chrX", 100, 200)];
        let genome = genome_info(&[("chr1", 1000)]);
        assert!(complement(&intervals, &genome).is_err());
    }

    #[test]
    fn test_complement_interval_exceeds_chrom() {
        let intervals = vec![iv("chr1", 900, 1100)];
        let genome = genome_info(&[("chr1", 1000)]);
        assert!(complement(&intervals, &genome).is_err());
    }

    #[test]
    fn test_complement_chrom_with_no_intervals() {
        let intervals = vec![iv("chr1", 100, 200)];
        let genome = genome_info(&[("chr1", 500), ("chr2", 300)]);
        let r = complement(&intervals, &genome).unwrap();
        // chr1: [0,100), [200,500); chr2: [0,300)
        assert_eq!(r.len(), 3);
        let chr2: Vec<_> = r.iter().filter(|i| i.chrom == "chr2").collect();
        assert_eq!(chr2.len(), 1);
        assert_eq!((chr2[0].start, chr2[0].end), (0, 300));
    }

    // -----------------------------------------------------------------------
    // closest
    // -----------------------------------------------------------------------

    #[test]
    fn test_closest_upstream() {
        let a = vec![iv("chr1", 300, 400)];
        let b = vec![iv("chr1", 100, 200)];
        let r = closest(&a, &b, StrandMode::Ignore);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].distance, Some(100));
    }

    #[test]
    fn test_closest_downstream() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 300, 400)];
        let r = closest(&a, &b, StrandMode::Ignore);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].distance, Some(100));
    }

    #[test]
    fn test_closest_overlapping() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 400)];
        let r = closest(&a, &b, StrandMode::Ignore);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].distance, Some(0));
    }

    #[test]
    fn test_closest_no_b_on_chrom() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr2", 100, 200)];
        let r = closest(&a, &b, StrandMode::Ignore);
        assert_eq!(r.len(), 1);
        assert!(r[0].closest.is_none());
        assert!(r[0].distance.is_none());
    }

    #[test]
    fn test_closest_single_b() {
        let a = vec![iv("chr1", 500, 600)];
        let b = vec![iv("chr1", 100, 200)];
        let r = closest(&a, &b, StrandMode::Ignore);
        assert_eq!(r[0].distance, Some(300));
    }

    #[test]
    fn test_closest_equidistant() {
        let a = vec![iv("chr1", 200, 300)];
        let b = vec![iv("chr1", 100, 150), iv("chr1", 350, 400)];
        let r = closest(&a, &b, StrandMode::Ignore);
        assert_eq!(r[0].distance, Some(50));
    }

    // -----------------------------------------------------------------------
    // windows
    // -----------------------------------------------------------------------

    #[test]
    fn test_make_windows_basic() {
        let genome = genome_info(&[("chr1", 100)]);
        let r = make_windows(&genome, 30).unwrap();
        assert_eq!(r.len(), 4); // [0,30), [30,60), [60,90), [90,100)
        assert_eq!((r[3].start, r[3].end), (90, 100));
    }

    #[test]
    fn test_make_windows_exact() {
        let genome = genome_info(&[("chr1", 100)]);
        let r = make_windows(&genome, 50).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!((r[0].start, r[0].end), (0, 50));
        assert_eq!((r[1].start, r[1].end), (50, 100));
    }

    #[test]
    fn test_make_windows_zero_size() {
        let genome = genome_info(&[("chr1", 100)]);
        assert!(make_windows(&genome, 0).is_err());
    }

    #[test]
    fn test_make_windows_larger_than_chrom() {
        let genome = genome_info(&[("chr1", 50)]);
        let r = make_windows(&genome, 100).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].start, r[0].end), (0, 50));
    }

    #[test]
    fn test_make_sliding_windows() {
        let genome = genome_info(&[("chr1", 100)]);
        let r = make_sliding_windows(&genome, 50, 25).unwrap();
        // [0,50), [25,75), [50,100), [75,100)
        assert_eq!(r.len(), 4);
        assert_eq!((r[1].start, r[1].end), (25, 75));
    }

    #[test]
    fn test_make_sliding_windows_zero_step() {
        let genome = genome_info(&[("chr1", 100)]);
        assert!(make_sliding_windows(&genome, 50, 0).is_err());
    }

    #[test]
    fn test_windows_around() {
        let intervals = vec![iv("chr1", 100, 200)];
        let genome = genome_info(&[("chr1", 1000)]);
        let r = windows_around(&intervals, 50, &genome).unwrap();
        assert_eq!(r.len(), 1);
        // midpoint = 150, half = 25, so [125, 175)
        assert_eq!((r[0].start, r[0].end), (125, 175));
    }

    #[test]
    fn test_windows_around_clipping() {
        let intervals = vec![iv("chr1", 0, 10)];
        let genome = genome_info(&[("chr1", 100)]);
        let r = windows_around(&intervals, 40, &genome).unwrap();
        assert_eq!(r.len(), 1);
        // midpoint = 5, half = 20, start = 5-20 saturates to 0, end = 5+20 = 25
        assert_eq!(r[0].start, 0);
        assert_eq!(r[0].end, 25);
    }

    // -----------------------------------------------------------------------
    // jaccard
    // -----------------------------------------------------------------------

    #[test]
    fn test_jaccard_identical() {
        let a = vec![iv("chr1", 100, 200)];
        let j = jaccard(&a, &a);
        assert!((j - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_jaccard_no_overlap() {
        let a = vec![iv("chr1", 100, 200)];
        let b = vec![iv("chr1", 300, 400)];
        let j = jaccard(&a, &b);
        assert!((j - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_jaccard_partial_overlap() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 400)];
        // intersection = [200,300) = 100bp, union = [100,400) = 300bp
        let j = jaccard(&a, &b);
        assert!((j - 100.0 / 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_jaccard_both_empty() {
        let a: Vec<GenomicInterval> = vec![];
        let b: Vec<GenomicInterval> = vec![];
        assert!((jaccard(&a, &b) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_jaccard_stats_bp_counts() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 400)];
        let stats = jaccard_stats(&a, &b);
        assert_eq!(stats.intersection_bp, 100);
        assert_eq!(stats.union_bp, 300);
        assert_eq!(stats.n_intersections, 1);
    }

    // -----------------------------------------------------------------------
    // property-style tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_subtract_plus_intersect_equals_a() {
        // Property holds on non-overlapping a (merged first)
        let a = merge(
            &[iv("chr1", 100, 300), iv("chr1", 250, 400)],
            StrandMode::Ignore,
        );
        let b = vec![iv("chr1", 200, 350)];

        let sub = subtract(&a, &b, StrandMode::Ignore).unwrap();
        let isect = intersect(&a, &b, StrandMode::Ignore).unwrap();

        let sub_bp: u64 = sub.iter().map(|i| i.len()).sum();
        let isect_bp: u64 = isect.iter().map(|i| i.len()).sum();
        let a_bp: u64 = a.iter().map(|i| i.len()).sum();

        assert_eq!(sub_bp + isect_bp, a_bp);
    }

    #[test]
    fn test_jaccard_symmetry() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 500)];
        let j1 = jaccard(&a, &b);
        let j2 = jaccard(&b, &a);
        assert!((j1 - j2).abs() < 1e-10);
    }

    #[test]
    fn test_jaccard_bounds() {
        let a = vec![iv("chr1", 100, 300)];
        let b = vec![iv("chr1", 200, 500)];
        let j = jaccard(&a, &b);
        assert!(j >= 0.0 && j <= 1.0);
    }

    #[test]
    fn test_complement_of_complement_equals_merge() {
        let a = vec![iv("chr1", 100, 200), iv("chr1", 300, 400)];
        let genome = genome_info(&[("chr1", 500)]);

        let comp = complement(&a, &genome).unwrap();
        let comp2 = complement(&comp, &genome).unwrap();

        let merged_a = merge(&a, StrandMode::Ignore);
        assert_eq!(comp2.len(), merged_a.len());
        for (c, m) in comp2.iter().zip(merged_a.iter()) {
            assert_eq!(c.start, m.start);
            assert_eq!(c.end, m.end);
        }
    }
}
