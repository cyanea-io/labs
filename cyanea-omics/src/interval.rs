//! Sorted interval collection with overlap queries.
//!
//! [`IntervalSet`] stores [`GenomicInterval`]s in a `Vec` and supports
//! overlap/containment queries. After calling [`IntervalSet::build_index()`],
//! queries are accelerated by per-chromosome interval trees for O(log n + k)
//! performance.

use std::collections::{BTreeMap, BTreeSet};

use cyanea_core::Summarizable;

use crate::genomic::GenomicInterval;
use crate::interval_tree::{Interval, IntervalTree};

/// A collection of genomic intervals with overlap query support.
///
/// Queries use O(n) linear scan by default. Call [`IntervalSet::build_index()`]
/// to construct per-chromosome interval trees for O(log n + k) queries.
pub struct IntervalSet {
    intervals: Vec<GenomicInterval>,
    sorted: bool,
    /// Per-chromosome interval trees. Values store indices into `intervals`.
    trees: Option<BTreeMap<String, IntervalTree<usize>>>,
}

impl IntervalSet {
    /// Create an empty interval set.
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
            sorted: true,
            trees: None,
        }
    }

    /// Create an interval set from existing intervals.
    pub fn from_intervals(intervals: Vec<GenomicInterval>) -> Self {
        Self {
            sorted: false,
            intervals,
            trees: None,
        }
    }

    /// Add an interval. Marks the set as unsorted and invalidates the index.
    pub fn push(&mut self, interval: GenomicInterval) {
        self.sorted = false;
        self.trees = None;
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

    /// Build per-chromosome interval tree indices for O(log n + k) queries.
    ///
    /// This is called automatically on first query if the set is sorted,
    /// but can also be called explicitly. Invalidated by [`push()`](IntervalSet::push).
    pub fn build_index(&mut self) {
        let mut by_chrom: BTreeMap<String, Vec<Interval<usize>>> = BTreeMap::new();
        for (idx, iv) in self.intervals.iter().enumerate() {
            by_chrom
                .entry(iv.chrom.clone())
                .or_default()
                .push(Interval::new(iv.start, iv.end, idx));
        }

        let trees: BTreeMap<String, IntervalTree<usize>> = by_chrom
            .into_iter()
            .map(|(chrom, intervals)| (chrom, IntervalTree::from_unsorted(intervals)))
            .collect();

        self.trees = Some(trees);
    }

    /// Return all intervals that overlap the query.
    ///
    /// Uses the interval tree index if available, otherwise falls back to O(n) scan.
    pub fn overlapping(&self, query: &GenomicInterval) -> Vec<&GenomicInterval> {
        if let Some(trees) = &self.trees {
            if let Some(tree) = trees.get(&query.chrom) {
                return tree
                    .query(query.start, query.end)
                    .into_iter()
                    .map(|hit| &self.intervals[hit.data])
                    .collect();
            }
            return Vec::new();
        }

        // Fallback: linear scan
        self.intervals
            .iter()
            .filter(|iv| iv.overlaps(query))
            .collect()
    }

    /// Return all intervals containing a given position.
    ///
    /// Uses the interval tree index if available, otherwise falls back to O(n) scan.
    pub fn containing(&self, chrom: &str, position: u64) -> Vec<&GenomicInterval> {
        if let Some(trees) = &self.trees {
            if let Some(tree) = trees.get(chrom) {
                return tree
                    .query(position, position + 1)
                    .into_iter()
                    .map(|hit| &self.intervals[hit.data])
                    .collect();
            }
            return Vec::new();
        }

        // Fallback: linear scan
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
            trees: None,
        }
    }

    /// Borrow the intervals as a slice.
    pub fn intervals(&self) -> &[GenomicInterval] {
        &self.intervals
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

    // --- Indexed query tests ---

    #[test]
    fn test_indexed_overlapping() {
        let mut set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr1", 300, 400),
            iv("chr2", 100, 200),
        ]);
        set.build_index();

        let query = iv("chr1", 150, 250);
        let hits = set.overlapping(&query);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 100);
    }

    #[test]
    fn test_indexed_containing() {
        let mut set = IntervalSet::from_intervals(vec![
            iv("chr1", 100, 200),
            iv("chr1", 150, 250),
            iv("chr1", 300, 400),
        ]);
        set.build_index();

        let hits = set.containing("chr1", 175);
        assert_eq!(hits.len(), 2);

        let hits = set.containing("chr1", 350);
        assert_eq!(hits.len(), 1);

        let hits = set.containing("chr2", 100);
        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn test_indexed_matches_linear() {
        let intervals = vec![
            iv("chr1", 0, 50),
            iv("chr1", 30, 80),
            iv("chr1", 100, 200),
            iv("chr2", 0, 100),
        ];
        let linear_set = IntervalSet::from_intervals(intervals.clone());
        let mut indexed_set = IntervalSet::from_intervals(intervals);
        indexed_set.build_index();

        let query = iv("chr1", 40, 120);
        let mut linear_hits: Vec<u64> = linear_set.overlapping(&query).iter().map(|iv| iv.start).collect();
        let mut indexed_hits: Vec<u64> = indexed_set.overlapping(&query).iter().map(|iv| iv.start).collect();
        linear_hits.sort();
        indexed_hits.sort();
        assert_eq!(linear_hits, indexed_hits);
    }

    #[test]
    fn test_push_invalidates_index() {
        let mut set = IntervalSet::from_intervals(vec![iv("chr1", 100, 200)]);
        set.build_index();
        assert!(set.trees.is_some());

        set.push(iv("chr1", 300, 400));
        assert!(set.trees.is_none());
    }
}
