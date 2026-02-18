//! Assembly quality-control metrics.
//!
//! Compute standard assembly statistics (N50, L50, N90, GC content, auN)
//! from a set of contig sequences.

use cyanea_core::{CyaneaError, Result};

/// Assembly quality statistics.
#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct AssemblyStats {
    /// Number of contigs.
    pub n_contigs: usize,
    /// Total assembly length in bases.
    pub total_length: usize,
    /// Length of the largest contig.
    pub largest_contig: usize,
    /// Length of the smallest contig.
    pub smallest_contig: usize,
    /// GC content as a fraction (0.0–1.0), ignoring N bases.
    pub gc_content: f64,
    /// N50: the length such that contigs of this length or longer cover ≥ 50% of the assembly.
    pub n50: usize,
    /// L50: the number of contigs needed to reach N50 cumulative length.
    pub l50: usize,
    /// N90: the length such that contigs of this length or longer cover ≥ 90% of the assembly.
    pub n90: usize,
    /// L90: the number of contigs needed to reach N90 cumulative length.
    pub l90: usize,
    /// Area under the Nx curve (auN = Σ length² / total_length).
    #[allow(non_snake_case)]
    pub auN: f64,
}

/// Compute assembly statistics for a set of contigs.
///
/// # Errors
///
/// Returns an error if `contigs` is empty.
pub fn assembly_stats(contigs: &[&[u8]]) -> Result<AssemblyStats> {
    if contigs.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "at least one contig is required".into(),
        ));
    }

    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a)); // descending

    let total_length: usize = lengths.iter().sum();
    let largest_contig = lengths[0];
    let smallest_contig = *lengths.last().unwrap();
    let n_contigs = lengths.len();

    // GC content: count G+C over all non-N bases.
    let mut gc = 0usize;
    let mut non_n = 0usize;
    for contig in contigs {
        for &b in *contig {
            let upper = b.to_ascii_uppercase();
            match upper {
                b'G' | b'C' => {
                    gc += 1;
                    non_n += 1;
                }
                b'A' | b'T' => {
                    non_n += 1;
                }
                _ => {} // N and others ignored
            }
        }
    }
    let gc_content = if non_n > 0 {
        gc as f64 / non_n as f64
    } else {
        0.0
    };

    let (n50, l50) = nx_from_sorted(&lengths, total_length, 0.5);
    let (n90, l90) = nx_from_sorted(&lengths, total_length, 0.9);

    // auN = Σ length² / total_length
    #[allow(non_snake_case)]
    let auN: f64 = lengths.iter().map(|&l| (l as f64).powi(2)).sum::<f64>() / total_length as f64;

    Ok(AssemblyStats {
        n_contigs,
        total_length,
        largest_contig,
        smallest_contig,
        gc_content,
        n50,
        l50,
        n90,
        l90,
        auN,
    })
}

/// Compute Nx and Lx values for a given fraction `x` (0.0–1.0).
///
/// Nx is the contig length at which cumulative length first reaches x × total.
/// Lx is the number of contigs required to reach that threshold.
///
/// # Errors
///
/// Returns an error if `contigs` is empty or `x` is not in (0, 1].
pub fn nx_values(contigs: &[&[u8]], x: f64) -> Result<(usize, usize)> {
    if contigs.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "at least one contig is required".into(),
        ));
    }
    if x <= 0.0 || x > 1.0 {
        return Err(CyaneaError::InvalidInput(format!(
            "x must be in (0, 1], got {}",
            x
        )));
    }

    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let total: usize = lengths.iter().sum();

    Ok(nx_from_sorted(&lengths, total, x))
}

/// Internal: compute Nx/Lx from pre-sorted (descending) lengths.
fn nx_from_sorted(lengths: &[usize], total: usize, x: f64) -> (usize, usize) {
    let threshold = (total as f64 * x).ceil() as usize;
    let mut cumulative = 0usize;
    for (i, &len) in lengths.iter().enumerate() {
        cumulative += len;
        if cumulative >= threshold {
            return (len, i + 1);
        }
    }
    // Should not reach here if lengths sum to total.
    (*lengths.last().unwrap(), lengths.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembly_stats_simple() {
        let contigs: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGT", b"AC"];
        let stats = assembly_stats(&contigs).unwrap();
        assert_eq!(stats.n_contigs, 3);
        assert_eq!(stats.total_length, 14);
        assert_eq!(stats.largest_contig, 8);
        assert_eq!(stats.smallest_contig, 2);
    }

    #[test]
    fn n50_known_values() {
        // Contigs: 10, 8, 6, 4, 2 → total 30, threshold = 15
        // cumulative: 10 → 18 ≥ 15, so N50=8, wait let's re-check:
        // sorted desc: 10, 8, 6, 4, 2. cumulative after 10: 10 < 15. after 10+8=18 ≥ 15. N50=8, L50=2.
        let c1 = vec![b'A'; 10];
        let c2 = vec![b'A'; 8];
        let c3 = vec![b'A'; 6];
        let c4 = vec![b'A'; 4];
        let c5 = vec![b'A'; 2];
        let contigs: Vec<&[u8]> = vec![&c1, &c2, &c3, &c4, &c5];
        let (n50, l50) = nx_values(&contigs, 0.5).unwrap();
        assert_eq!(n50, 8);
        assert_eq!(l50, 2);
    }

    #[test]
    fn gc_content_correct() {
        let contigs: Vec<&[u8]> = vec![b"GGCC", b"AATT"];
        let stats = assembly_stats(&contigs).unwrap();
        assert!((stats.gc_content - 0.5).abs() < 1e-10);
    }

    #[test]
    fn empty_contigs_error() {
        let contigs: Vec<&[u8]> = vec![];
        assert!(assembly_stats(&contigs).is_err());
    }
}
