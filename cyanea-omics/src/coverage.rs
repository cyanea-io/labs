//! Run-length encoded coverage vectors for memory-efficient genome-wide depth.
//!
//! [`RleCoverage`] stores depth values as (depth, run_length) pairs, which is
//! highly efficient for typical genomic coverage profiles where long stretches
//! share the same depth.

use crate::genomic::GenomicInterval;

/// Run-length encoded coverage vector.
///
/// Stores coverage depth as a sequence of (depth, length) runs.
/// This is memory-efficient for genome-wide coverage where long runs of
/// identical depth are common (e.g., depth=0 in intergenic regions).
#[derive(Debug, Clone)]
pub struct RleCoverage {
    /// (depth, run_length) pairs in positional order.
    runs: Vec<(u32, u64)>,
    /// Total number of positions covered by this vector.
    total_length: u64,
}

impl RleCoverage {
    /// Build an RLE coverage vector from a set of intervals on a chromosome.
    ///
    /// Each interval contributes +1 to all positions it covers. The resulting
    /// coverage vector spans `[0, chrom_length)`.
    pub fn from_intervals(intervals: &[GenomicInterval], chrom_length: u64) -> Self {
        if chrom_length == 0 {
            return Self {
                runs: Vec::new(),
                total_length: 0,
            };
        }

        // Use sweep-line approach: collect +1/-1 events
        let mut events: Vec<(u64, i32)> = Vec::with_capacity(intervals.len() * 2);
        for iv in intervals {
            let start = iv.start.min(chrom_length);
            let end = iv.end.min(chrom_length);
            if start < end {
                events.push((start, 1));
                events.push((end, -1));
            }
        }
        events.sort_by_key(|&(pos, delta)| (pos, std::cmp::Reverse(delta)));

        // Sweep line
        let mut runs = Vec::new();
        let mut depth: i64 = 0;
        let mut prev_pos: u64 = 0;

        for (pos, delta) in &events {
            let pos = *pos;
            if pos > prev_pos {
                runs.push((depth as u32, pos - prev_pos));
                prev_pos = pos;
            }
            depth += *delta as i64;
        }

        // Remaining positions to end of chromosome
        if prev_pos < chrom_length {
            runs.push((depth as u32, chrom_length - prev_pos));
        }

        // Merge consecutive runs with same depth
        let merged = merge_runs(runs);

        Self {
            runs: merged,
            total_length: chrom_length,
        }
    }

    /// Build an RLE coverage vector from a dense depth array.
    pub fn from_depths(depths: &[u32]) -> Self {
        if depths.is_empty() {
            return Self {
                runs: Vec::new(),
                total_length: 0,
            };
        }

        let mut runs = Vec::new();
        let mut current_depth = depths[0];
        let mut run_len: u64 = 1;

        for &d in &depths[1..] {
            if d == current_depth {
                run_len += 1;
            } else {
                runs.push((current_depth, run_len));
                current_depth = d;
                run_len = 1;
            }
        }
        runs.push((current_depth, run_len));

        Self {
            runs,
            total_length: depths.len() as u64,
        }
    }

    /// Get the coverage depth at a specific position.
    ///
    /// Returns 0 if the position is out of range.
    pub fn get(&self, position: u64) -> u32 {
        if position >= self.total_length {
            return 0;
        }

        let mut offset: u64 = 0;
        for &(depth, length) in &self.runs {
            offset += length;
            if position < offset {
                return depth;
            }
        }

        0
    }

    /// Mean coverage across all positions.
    pub fn mean_coverage(&self) -> f64 {
        if self.total_length == 0 {
            return 0.0;
        }

        let sum: u64 = self
            .runs
            .iter()
            .map(|&(depth, length)| depth as u64 * length)
            .sum();

        sum as f64 / self.total_length as f64
    }

    /// Number of bases with depth above a threshold.
    pub fn bases_above(&self, threshold: u32) -> u64 {
        self.runs
            .iter()
            .filter(|&&(depth, _)| depth > threshold)
            .map(|&(_, length)| length)
            .sum()
    }

    /// Expand the RLE coverage to a dense vector.
    pub fn to_vec(&self) -> Vec<u32> {
        let mut result = Vec::with_capacity(self.total_length as usize);
        for &(depth, length) in &self.runs {
            result.extend(std::iter::repeat(depth).take(length as usize));
        }
        result
    }

    /// Total length of the coverage vector.
    pub fn len(&self) -> u64 {
        self.total_length
    }

    /// Whether the coverage vector is empty.
    pub fn is_empty(&self) -> bool {
        self.total_length == 0
    }

    /// Number of RLE runs.
    pub fn num_runs(&self) -> usize {
        self.runs.len()
    }
}

fn merge_runs(runs: Vec<(u32, u64)>) -> Vec<(u32, u64)> {
    if runs.is_empty() {
        return runs;
    }

    let mut merged = Vec::with_capacity(runs.len());
    let mut current = runs[0];

    for &(depth, length) in &runs[1..] {
        if depth == current.0 {
            current.1 += length;
        } else {
            merged.push(current);
            current = (depth, length);
        }
    }
    merged.push(current);

    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: u64, end: u64) -> GenomicInterval {
        GenomicInterval::new("chr1", start, end).unwrap()
    }

    #[test]
    fn empty_coverage() {
        let cov = RleCoverage::from_intervals(&[], 100);
        assert_eq!(cov.len(), 100);
        assert_eq!(cov.mean_coverage(), 0.0);
        assert_eq!(cov.bases_above(0), 0);
        assert_eq!(cov.get(50), 0);
    }

    #[test]
    fn zero_length() {
        let cov = RleCoverage::from_intervals(&[], 0);
        assert!(cov.is_empty());
        assert_eq!(cov.mean_coverage(), 0.0);
    }

    #[test]
    fn single_interval() {
        let cov = RleCoverage::from_intervals(&[iv(10, 20)], 100);
        assert_eq!(cov.get(5), 0);
        assert_eq!(cov.get(10), 1);
        assert_eq!(cov.get(15), 1);
        assert_eq!(cov.get(19), 1);
        assert_eq!(cov.get(20), 0);
        assert_eq!(cov.get(50), 0);
        assert_eq!(cov.bases_above(0), 10);
        assert!((cov.mean_coverage() - 0.1).abs() < 1e-10);
    }

    #[test]
    fn overlapping_intervals() {
        let cov = RleCoverage::from_intervals(&[iv(10, 30), iv(20, 40)], 100);
        assert_eq!(cov.get(15), 1); // only first
        assert_eq!(cov.get(25), 2); // both
        assert_eq!(cov.get(35), 1); // only second
        assert_eq!(cov.get(45), 0); // neither
    }

    #[test]
    fn from_depths_basic() {
        let depths = vec![0, 0, 1, 1, 1, 2, 2, 0];
        let cov = RleCoverage::from_depths(&depths);
        assert_eq!(cov.len(), 8);
        assert_eq!(cov.num_runs(), 4);
        assert_eq!(cov.get(0), 0);
        assert_eq!(cov.get(2), 1);
        assert_eq!(cov.get(5), 2);
        assert_eq!(cov.get(7), 0);
    }

    #[test]
    fn from_depths_empty() {
        let cov = RleCoverage::from_depths(&[]);
        assert!(cov.is_empty());
    }

    #[test]
    fn round_trip() {
        let depths = vec![0, 0, 3, 3, 3, 1, 0, 0, 0, 5];
        let cov = RleCoverage::from_depths(&depths);
        let expanded = cov.to_vec();
        assert_eq!(expanded, depths);
    }

    #[test]
    fn mean_coverage_value() {
        let cov = RleCoverage::from_depths(&[0, 1, 2, 3, 4]);
        assert!((cov.mean_coverage() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn bases_above_threshold() {
        let cov = RleCoverage::from_depths(&[0, 1, 2, 3, 4, 5]);
        assert_eq!(cov.bases_above(2), 3); // positions with depth 3, 4, 5
        assert_eq!(cov.bases_above(0), 5); // positions with depth > 0
        assert_eq!(cov.bases_above(10), 0);
    }

    #[test]
    fn get_out_of_range() {
        let cov = RleCoverage::from_depths(&[1, 2, 3]);
        assert_eq!(cov.get(3), 0);
        assert_eq!(cov.get(100), 0);
    }

    #[test]
    fn interval_coverage_to_vec() {
        let cov = RleCoverage::from_intervals(&[iv(2, 5)], 8);
        let expanded = cov.to_vec();
        assert_eq!(expanded, vec![0, 0, 1, 1, 1, 0, 0, 0]);
    }

    #[test]
    fn heavily_overlapping() {
        // 5 intervals all covering [0, 10)
        let intervals: Vec<GenomicInterval> = (0..5).map(|_| iv(0, 10)).collect();
        let cov = RleCoverage::from_intervals(&intervals, 10);
        assert_eq!(cov.get(0), 5);
        assert_eq!(cov.get(9), 5);
        assert_eq!(cov.num_runs(), 1); // one run of depth 5
    }

    #[test]
    fn uniform_depth() {
        let depths = vec![3; 1000];
        let cov = RleCoverage::from_depths(&depths);
        assert_eq!(cov.num_runs(), 1);
        assert!((cov.mean_coverage() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn intervals_beyond_chrom_length() {
        // Interval extends past chrom_length — should be clamped
        let cov = RleCoverage::from_intervals(&[iv(5, 200)], 10);
        assert_eq!(cov.len(), 10);
        assert_eq!(cov.get(5), 1);
        assert_eq!(cov.get(9), 1);
        let expanded = cov.to_vec();
        assert_eq!(expanded.len(), 10);
    }
}
