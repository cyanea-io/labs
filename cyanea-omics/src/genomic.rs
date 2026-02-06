//! Genomic primitives — strand, position, and interval types.
//!
//! These are the foundation types used across the omics crate for representing
//! genomic coordinates. All coordinates are 0-based, half-open `[start, end)`
//! unless otherwise noted.

use core::fmt;

use cyanea_core::{CyaneaError, Result};

/// Strand orientation on a reference genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    /// Returns `true` if this is the forward (+) strand.
    pub fn is_forward(&self) -> bool {
        matches!(self, Strand::Forward)
    }

    /// Returns `true` if this is the reverse (-) strand.
    pub fn is_reverse(&self) -> bool {
        matches!(self, Strand::Reverse)
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

/// A single position on a chromosome (0-based).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GenomicPosition {
    pub chrom: String,
    pub position: u64,
}

/// A half-open interval `[start, end)` on a chromosome (0-based).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GenomicInterval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

impl GenomicInterval {
    /// Create a new interval on the unknown strand.
    ///
    /// Returns an error if `start >= end`.
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Result<Self> {
        Self::with_strand(chrom, start, end, Strand::Unknown)
    }

    /// Create a new interval with an explicit strand.
    ///
    /// Returns an error if `start >= end`.
    pub fn with_strand(
        chrom: impl Into<String>,
        start: u64,
        end: u64,
        strand: Strand,
    ) -> Result<Self> {
        if start >= end {
            return Err(CyaneaError::InvalidInput(format!(
                "interval start ({start}) must be less than end ({end})"
            )));
        }
        Ok(Self {
            chrom: chrom.into(),
            start,
            end,
            strand,
        })
    }

    /// Length of the interval in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the interval has zero length (never true for valid intervals).
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// Whether a 0-based position falls within `[start, end)`.
    pub fn contains(&self, position: u64) -> bool {
        position >= self.start && position < self.end
    }

    /// Whether this interval overlaps with `other` (must be on the same chromosome).
    pub fn overlaps(&self, other: &GenomicInterval) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    /// Return the overlapping region, or `None` if the intervals don't overlap.
    pub fn intersect(&self, other: &GenomicInterval) -> Option<GenomicInterval> {
        if !self.overlaps(other) {
            return None;
        }
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        Some(GenomicInterval {
            chrom: self.chrom.clone(),
            start,
            end,
            strand: self.strand,
        })
    }

    /// Merge two overlapping intervals into their union, or `None` if they don't overlap.
    pub fn merge(&self, other: &GenomicInterval) -> Option<GenomicInterval> {
        if !self.overlaps(other) {
            return None;
        }
        let start = self.start.min(other.start);
        let end = self.end.max(other.end);
        Some(GenomicInterval {
            chrom: self.chrom.clone(),
            start,
            end,
            strand: self.strand,
        })
    }

    /// Distance between two intervals, or `None` if they are on different chromosomes.
    /// Returns 0 if the intervals overlap.
    pub fn distance(&self, other: &GenomicInterval) -> Option<u64> {
        if self.chrom != other.chrom {
            return None;
        }
        if self.overlaps(other) {
            return Some(0);
        }
        if self.end <= other.start {
            Some(other.start - self.end)
        } else {
            Some(self.start - other.end)
        }
    }

    /// The midpoint of the interval (rounded down).
    pub fn midpoint(&self) -> u64 {
        (self.start + self.end) / 2
    }
}

impl fmt::Display for GenomicInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}-{}({})",
            self.chrom, self.start, self.end, self.strand
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_display() {
        assert_eq!(Strand::Forward.to_string(), "+");
        assert_eq!(Strand::Reverse.to_string(), "-");
        assert_eq!(Strand::Unknown.to_string(), ".");
    }

    #[test]
    fn test_strand_predicates() {
        assert!(Strand::Forward.is_forward());
        assert!(!Strand::Forward.is_reverse());
        assert!(Strand::Reverse.is_reverse());
        assert!(!Strand::Unknown.is_forward());
    }

    #[test]
    fn test_interval_new() {
        let iv = GenomicInterval::new("chr1", 100, 200).unwrap();
        assert_eq!(iv.chrom, "chr1");
        assert_eq!(iv.start, 100);
        assert_eq!(iv.end, 200);
        assert_eq!(iv.strand, Strand::Unknown);
    }

    #[test]
    fn test_interval_invalid() {
        assert!(GenomicInterval::new("chr1", 200, 100).is_err());
        assert!(GenomicInterval::new("chr1", 100, 100).is_err());
    }

    #[test]
    fn test_interval_with_strand() {
        let iv = GenomicInterval::with_strand("chr1", 0, 50, Strand::Forward).unwrap();
        assert_eq!(iv.strand, Strand::Forward);
    }

    #[test]
    fn test_interval_len() {
        let iv = GenomicInterval::new("chr1", 100, 300).unwrap();
        assert_eq!(iv.len(), 200);
    }

    #[test]
    fn test_interval_contains() {
        let iv = GenomicInterval::new("chr1", 100, 200).unwrap();
        assert!(iv.contains(100));
        assert!(iv.contains(150));
        assert!(iv.contains(199));
        assert!(!iv.contains(200));
        assert!(!iv.contains(99));
    }

    #[test]
    fn test_interval_overlaps() {
        let a = GenomicInterval::new("chr1", 100, 200).unwrap();
        let b = GenomicInterval::new("chr1", 150, 250).unwrap();
        let c = GenomicInterval::new("chr1", 200, 300).unwrap();
        let d = GenomicInterval::new("chr2", 100, 200).unwrap();

        assert!(a.overlaps(&b));
        assert!(b.overlaps(&a));
        assert!(!a.overlaps(&c)); // abutting, not overlapping
        assert!(!a.overlaps(&d)); // different chrom
    }

    #[test]
    fn test_interval_intersect() {
        let a = GenomicInterval::new("chr1", 100, 200).unwrap();
        let b = GenomicInterval::new("chr1", 150, 250).unwrap();

        let isect = a.intersect(&b).unwrap();
        assert_eq!(isect.start, 150);
        assert_eq!(isect.end, 200);

        let c = GenomicInterval::new("chr1", 200, 300).unwrap();
        assert!(a.intersect(&c).is_none());
    }

    #[test]
    fn test_interval_merge() {
        let a = GenomicInterval::new("chr1", 100, 200).unwrap();
        let b = GenomicInterval::new("chr1", 150, 250).unwrap();

        let merged = a.merge(&b).unwrap();
        assert_eq!(merged.start, 100);
        assert_eq!(merged.end, 250);

        let c = GenomicInterval::new("chr1", 300, 400).unwrap();
        assert!(a.merge(&c).is_none());
    }

    #[test]
    fn test_interval_distance() {
        let a = GenomicInterval::new("chr1", 100, 200).unwrap();
        let b = GenomicInterval::new("chr1", 300, 400).unwrap();
        let c = GenomicInterval::new("chr2", 100, 200).unwrap();
        let d = GenomicInterval::new("chr1", 150, 250).unwrap();

        assert_eq!(a.distance(&b), Some(100));
        assert_eq!(b.distance(&a), Some(100));
        assert_eq!(a.distance(&c), None);
        assert_eq!(a.distance(&d), Some(0));
    }

    #[test]
    fn test_interval_midpoint() {
        let iv = GenomicInterval::new("chr1", 100, 200).unwrap();
        assert_eq!(iv.midpoint(), 150);
    }

    #[test]
    fn test_interval_display() {
        let iv = GenomicInterval::with_strand("chr1", 100, 200, Strand::Forward).unwrap();
        assert_eq!(iv.to_string(), "chr1:100-200(+)");
    }
}
