//! Sorted interval collection with overlap queries.
//!
//! [`IntervalSet`] stores [`GenomicInterval`]s in a `Vec` and supports
//! overlap/containment queries via linear scan. This is suitable for small to
//! moderate datasets; a future version may use an interval tree for O(log n + k)
//! queries.

use std::collections::BTreeSet;

use cyanea_core::Summarizable;

use crate::genomic::GenomicInterval;

/// A collection of genomic intervals with overlap query support.
///
/// Current query complexity is O(n). A planned upgrade will replace the
/// backing store with an interval tree for O(log n + k) overlap queries.
pub struct IntervalSet {
    intervals: Vec<GenomicInterval>,
    sorted: bool,
}

impl IntervalSet {
    /// Create an empty interval set.
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
            sorted: true,
        }
    }

    /// Create an interval set from existing intervals.
    pub fn from_intervals(intervals: Vec<GenomicInterval>) -> Self {
        Self {
            sorted: false,
            intervals,
        }
    }

    /// Add an interval. Marks the set as unsorted.
    pub fn push(&mut self, interval: GenomicInterval) {
        self.sorted = false;
        self.intervals.push(interval);
    }

    /// Number of intervals in the set.
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Sort intervals by chromosome then start position.
    pub fn sort(&mut self) {
        self.intervals
            .sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
        self.sorted = true;
    }

    /// Return all intervals that overlap the query. O(n) scan.
    pub fn overlapping(&self, query: &GenomicInterval) -> Vec<&GenomicInterval> {
        self.intervals
            .iter()
            .filter(|iv| iv.overlaps(query))
            .collect()
    }

    /// Return all intervals containing a given position. O(n) scan.
    pub fn containing(&self, chrom: &str, position: u64) -> Vec<&GenomicInterval> {
        self.intervals
            .iter()
            .filter(|iv| iv.chrom == chrom && iv.contains(position))
            .collect()
    }

    /// Merge all overlapping intervals and return a new, sorted set.
    ///
    /// Adjacent intervals on the same chromosome that overlap or abut are
    /// merged into a single interval.
    pub fn merge_overlapping(&self) -> IntervalSet {
        if self.intervals.is_empty() {
            return IntervalSet::new();
        }

        let mut sorted: Vec<GenomicInterval> = self.intervals.clone();
        sorted.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

        let mut merged = Vec::new();
        let mut current = sorted[0].clone();

        for iv in &sorted[1..] {
            if iv.chrom == current.chrom && iv.start <= current.end {
                current.end = current.end.max(iv.end);
            } else {
                merged.push(current);
                current = iv.clone();
            }
        }
        merged.push(current);

        IntervalSet {
            intervals: merged,
            sorted: true,
        }
    }

    /// Consume the set and return the inner intervals.
    pub fn into_intervals(self) -> Vec<GenomicInterval> {
        self.intervals
    }

    /// Total bases covered on a given chromosome (after merging overlaps).
    pub fn coverage(&self, chrom: &str) -> u64 {
        let chrom_intervals: Vec<GenomicInterval> = self
            .intervals
            .iter()
            .filter(|iv| iv.chrom == chrom)
            .cloned()
            .collect();

        if chrom_intervals.is_empty() {
            return 0;
        }

        let subset = IntervalSet::from_intervals(chrom_intervals);
        let merged = subset.merge_overlapping();
        merged.intervals.iter().map(|iv| iv.len()).sum()
    }
}

impl Default for IntervalSet {
    fn default() -> Self {
        Self::new()
    }
}

impl Summarizable for IntervalSet {
    fn summary(&self) -> String {
        let chroms: BTreeSet<&str> = self.intervals.iter().map(|iv| iv.chrom.as_str()).collect();
        format!(
            "IntervalSet: {} intervals across {} chromosomes",
            self.len(),
            chroms.len()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> GenomicInterval {
        GenomicInterval::new(chrom, start, end).unwrap()
    }

    #[test]
    fn test_empty_set() {
        let set = IntervalSet::new();
        assert!(set.is_empty());
        assert_eq!(set.len(), 0);
    }

    #[test]
    fn test_push_and_sort() {
        let mut set = IntervalSet::new();
        set.push(iv("chr1", 200, 300));
        set.push(iv("chr1", 100, 150));
        assert!(!set.sorted);

        set.sort();
        assert!(set.sorted);
        assert_eq!(set.intervals[0].start, 100);
        assert_eq!(set.intervals[1].start, 200);
    }

    #[test]
    fn test_overlapping() {
        let set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr1", 300, 400),
            iv("chr2", 100, 200),
        ]);

        let query = iv("chr1", 150, 250);
        let hits = set.overlapping(&query);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 100);
    }

    #[test]
    fn test_containing() {
        let set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr1", 150, 250),
            iv("chr1", 300, 400),
        ]);

        let hits = set.containing("chr1", 175);
        assert_eq!(hits.len(), 2);

        let hits = set.containing("chr1", 350);
        assert_eq!(hits.len(), 1);

        let hits = set.containing("chr2", 100);
        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn test_merge_overlapping() {
        let set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr1", 150, 300),
            iv("chr1", 500, 600),
            iv("chr2", 100, 200),
        ]);

        let merged = set.merge_overlapping();
        assert_eq!(merged.len(), 3);
        assert_eq!(merged.intervals[0].start, 100);
        assert_eq!(merged.intervals[0].end, 300);
        assert_eq!(merged.intervals[1].start, 500);
    }

    #[test]
    fn test_coverage() {
        let set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr1", 150, 300),
            iv("chr1", 500, 600),
        ]);

        // merged: [100,300) = 200bp + [500,600) = 100bp = 300bp
        assert_eq!(set.coverage("chr1"), 300);
        assert_eq!(set.coverage("chr2"), 0);
    }

    #[test]
    fn test_summary() {
        let set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr2", 100, 200),
        ]);
        assert_eq!(set.summary(), "IntervalSet: 2 intervals across 2 chromosomes");
    }

    #[test]
    fn test_merge_empty() {
        let set = IntervalSet::new();
        let merged = set.merge_overlapping();
        assert!(merged.is_empty());
    }
}
