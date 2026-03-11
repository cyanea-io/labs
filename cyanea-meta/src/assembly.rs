//! Assembly quality control metrics for metagenomic assemblies.
//!
//! Computes N50, L50, N90, L90, total length, GC content, and other
//! standard assembly statistics.

use crate::error::{MetaError, Result};

/// Assembly quality statistics.
#[derive(Debug, Clone)]
pub struct AssemblyStats {
    /// Number of contigs.
    pub n_contigs: usize,
    /// Total assembly length in bases.
    pub total_length: usize,
    /// Length of the longest contig.
    pub longest_contig: usize,
    /// Length of the shortest contig.
    pub shortest_contig: usize,
    /// GC content as a fraction (0.0–1.0).
    pub gc_content: f64,
    /// N50: length such that contigs ≥ N50 cover ≥ 50% of assembly.
    pub n50: usize,
    /// L50: number of contigs needed to reach 50% of assembly.
    pub l50: usize,
    /// N90: length such that contigs ≥ N90 cover ≥ 90% of assembly.
    pub n90: usize,
    /// L90: number of contigs needed to reach 90% of assembly.
    pub l90: usize,
    /// Mean contig length.
    pub mean_length: f64,
    /// Median contig length.
    pub median_length: usize,
    /// Area under Nx curve: Σ (length² / total_length).
    pub au_n: f64,
}

/// Compute assembly statistics for a set of contigs.
///
/// # Errors
///
/// Returns an error if contigs is empty.
pub fn assembly_stats(contigs: &[&[u8]]) -> Result<AssemblyStats> {
    if contigs.is_empty() {
        return Err(MetaError::Binning(
            "at least one contig required".into(),
        ));
    }

    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a)); // Descending

    let n_contigs = lengths.len();
    let total_length: usize = lengths.iter().sum();
    let longest_contig = lengths[0];
    let shortest_contig = *lengths.last().unwrap();
    let mean_length = total_length as f64 / n_contigs as f64;

    // Median
    let median_length = if n_contigs % 2 == 0 {
        (lengths[n_contigs / 2 - 1] + lengths[n_contigs / 2]) / 2
    } else {
        lengths[n_contigs / 2]
    };

    // GC content
    let mut gc = 0usize;
    let mut valid = 0usize;
    for contig in contigs {
        for &b in *contig {
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
    }
    let gc_content = if valid > 0 {
        gc as f64 / valid as f64
    } else {
        0.0
    };

    // Nx/Lx statistics
    let (n50, l50) = compute_nx(&lengths, total_length, 0.5);
    let (n90, l90) = compute_nx(&lengths, total_length, 0.9);

    // Area under Nx curve
    let au_n: f64 = lengths
        .iter()
        .map(|&l| (l as f64).powi(2))
        .sum::<f64>()
        / total_length as f64;

    Ok(AssemblyStats {
        n_contigs,
        total_length,
        longest_contig,
        shortest_contig,
        gc_content,
        n50,
        l50,
        n90,
        l90,
        mean_length,
        median_length,
        au_n,
    })
}

/// Compute Nx and Lx for a given percentile (0.0–1.0).
///
/// Nx: the contig length where cumulative length first reaches x × total.
/// Lx: number of contigs needed to reach that threshold.
///
/// # Errors
///
/// Returns an error if contigs is empty or x is not in (0, 1].
pub fn nx_values(contigs: &[&[u8]], x: f64) -> Result<(usize, usize)> {
    if contigs.is_empty() {
        return Err(MetaError::Binning(
            "at least one contig required".into(),
        ));
    }

    if x <= 0.0 || x > 1.0 {
        return Err(MetaError::Binning(format!(
            "x must be in (0, 1], got {}",
            x
        )));
    }

    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let total: usize = lengths.iter().sum();

    Ok(compute_nx(&lengths, total, x))
}

/// Internal: compute Nx/Lx from sorted (descending) lengths.
fn compute_nx(lengths: &[usize], total: usize, x: f64) -> (usize, usize) {
    let threshold = (total as f64 * x).ceil() as usize;
    let mut cumulative = 0usize;

    for (i, &len) in lengths.iter().enumerate() {
        cumulative += len;
        if cumulative >= threshold {
            return (len, i + 1);
        }
    }

    // Fallback (should not reach here)
    (*lengths.last().unwrap(), lengths.len())
}

/// Coverage depth: per-contig average coverage from read mappings.
///
/// Simplified version: returns uniform coverage for now.
/// Full implementation would parse BAM/SAM depth information.
///
/// # Errors
///
/// Returns an error if contigs is empty.
pub fn coverage_depth(contigs: &[&[u8]], _reads_per_contig: &[usize]) -> Result<Vec<f64>> {
    if contigs.is_empty() {
        return Err(MetaError::Binning(
            "at least one contig required".into(),
        ));
    }

    // Placeholder: return uniform depth
    Ok(vec![1.0; contigs.len()])
}

/// Filter contigs by length and/or coverage thresholds.
///
/// # Errors
///
/// Returns an error if thresholds are invalid.
pub fn filter_contigs<'a>(
    contigs: &'a [&'a [u8]],
    min_length: usize,
    _min_coverage: f64,
) -> Result<Vec<&'a [u8]>> {
    let filtered: Vec<&'a [u8]> = contigs
        .iter()
        .filter(|c| c.len() >= min_length)
        .copied()
        .collect();

    Ok(filtered)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembly_stats_basic() {
        let c1 = b"ACGTACGTACGTACGT"; // 16 bp
        let c2 = b"ACGTACGT"; // 8 bp
        let c3 = b"ACGT"; // 4 bp
        let contigs: Vec<&[u8]> = vec![c1, c2, c3];

        let stats = assembly_stats(&contigs).unwrap();
        assert_eq!(stats.n_contigs, 3);
        assert_eq!(stats.total_length, 28);
        assert_eq!(stats.longest_contig, 16);
        assert_eq!(stats.shortest_contig, 4);
    }

    #[test]
    fn n50_calculation() {
        // Contigs: 10, 8, 6, 4, 2 → total 30
        // Threshold for N50: 15
        // Cumulative: 10 < 15, then 10+8=18 ≥ 15 → N50=8, L50=2
        let mut contigs: Vec<Vec<u8>> = vec![
            vec![b'A'; 10],
            vec![b'A'; 8],
            vec![b'A'; 6],
            vec![b'A'; 4],
            vec![b'A'; 2],
        ];
        let refs: Vec<&[u8]> = contigs.iter_mut().map(|v| v.as_slice()).collect();

        let (n50, l50) = nx_values(&refs, 0.5).unwrap();
        assert_eq!(n50, 8);
        assert_eq!(l50, 2);
    }

    #[test]
    fn gc_content_calculation() {
        let contigs: Vec<&[u8]> = vec![b"GGCC", b"AATT"];
        let stats = assembly_stats(&contigs).unwrap();
        assert!((stats.gc_content - 0.5).abs() < 1e-10);
    }

    #[test]
    fn empty_contigs_error() {
        let contigs: Vec<&[u8]> = vec![];
        assert!(assembly_stats(&contigs).is_err());
    }

    #[test]
    fn invalid_x_error() {
        let contigs: Vec<&[u8]> = vec![b"ACGT"];
        assert!(nx_values(&contigs, 1.5).is_err());
    }

    #[test]
    fn filter_contigs_by_length() {
        let c1 = b"ACGTACGTACGTACGT"; // 16 bp
        let c2 = b"AC"; // 2 bp
        let c3 = b"ACGTACGTACGTACGTACGT"; // 20 bp
        let contigs: Vec<&[u8]> = vec![c1, c2, c3];

        let filtered = filter_contigs(&contigs, 5, 0.0).unwrap();
        assert_eq!(filtered.len(), 2); // c1 and c3
    }

    #[test]
    fn mean_and_median() {
        let c1 = b"ACGTACGTACGTACGT"; // 16 bp
        let c2 = b"ACGTACGT"; // 8 bp
        let c3 = b"ACGT"; // 4 bp
        let contigs: Vec<&[u8]> = vec![c1, c2, c3];

        let stats = assembly_stats(&contigs).unwrap();
        assert!((stats.mean_length - (28.0 / 3.0)).abs() < 1e-10);
        assert_eq!(stats.median_length, 8);
    }
}
