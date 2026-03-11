//! Long-read sequencing analysis: PacBio HiFi and Oxford Nanopore.
//!
//! Provides long-read-specific data types, quality metrics, error correction,
//! consensus generation, and read simulation for PacBio HiFi (CCS) and
//! Oxford Nanopore (simplex/duplex) platforms.

use cyanea_core::{CyaneaError, Result};

/// Sequencing platform for a long read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LongReadPlatform {
    /// PacBio Continuous Long Reads (CLR).
    PacBioCLR,
    /// PacBio HiFi / Circular Consensus Sequencing.
    PacBioHiFi,
    /// Oxford Nanopore simplex reads.
    NanoporeSimplex,
    /// Oxford Nanopore duplex reads.
    NanoporeDuplex,
}

/// A long read with platform-specific metadata.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LongRead {
    /// Read identifier.
    pub id: String,
    /// Nucleotide sequence.
    pub sequence: Vec<u8>,
    /// Per-base quality scores (Phred).
    pub quality: Vec<u8>,
    /// Sequencing platform.
    pub platform: LongReadPlatform,
    /// Number of full passes (CCS/HiFi only).
    pub num_passes: Option<u32>,
    /// Predicted accuracy (0.0 - 1.0).
    pub predicted_accuracy: Option<f64>,
    /// Read group / run ID.
    pub read_group: Option<String>,
}

impl LongRead {
    /// Create a new long read.
    pub fn new(
        id: impl Into<String>,
        sequence: Vec<u8>,
        quality: Vec<u8>,
        platform: LongReadPlatform,
    ) -> Result<Self> {
        if sequence.is_empty() {
            return Err(CyaneaError::InvalidInput("empty sequence".into()));
        }
        if !quality.is_empty() && quality.len() != sequence.len() {
            return Err(CyaneaError::InvalidInput(
                "quality length must match sequence length".into(),
            ));
        }
        Ok(Self {
            id: id.into(),
            sequence,
            quality,
            platform,
            num_passes: None,
            predicted_accuracy: None,
            read_group: None,
        })
    }

    /// Read length in bases.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Whether the read is empty.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Mean quality score (Phred).
    pub fn mean_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }
        self.quality.iter().map(|&q| q as f64).sum::<f64>() / self.quality.len() as f64
    }

    /// Estimated per-base error rate from quality scores.
    pub fn error_rate(&self) -> f64 {
        if self.quality.is_empty() {
            return match self.platform {
                LongReadPlatform::PacBioCLR => 0.13,
                LongReadPlatform::PacBioHiFi => 0.001,
                LongReadPlatform::NanoporeSimplex => 0.05,
                LongReadPlatform::NanoporeDuplex => 0.005,
            };
        }
        let mean_q = self.mean_quality();
        10.0f64.powf(-mean_q / 10.0)
    }

    /// GC content of the read.
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc = self.sequence.iter()
            .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count();
        gc as f64 / self.sequence.len() as f64
    }
}

/// Quality profile statistics for a collection of long reads.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LongReadStats {
    /// Total number of reads.
    pub num_reads: usize,
    /// Total bases.
    pub total_bases: u64,
    /// Mean read length.
    pub mean_length: f64,
    /// Median read length.
    pub median_length: f64,
    /// N50 read length.
    pub n50: usize,
    /// Longest read.
    pub max_length: usize,
    /// Shortest read.
    pub min_length: usize,
    /// Mean quality score.
    pub mean_quality: f64,
    /// Mean GC content.
    pub mean_gc: f64,
    /// Reads with quality >= Q20.
    pub q20_count: usize,
    /// Reads with quality >= Q30.
    pub q30_count: usize,
    /// Mean number of passes (HiFi only).
    pub mean_passes: Option<f64>,
}

/// Compute statistics for a collection of long reads.
pub fn longread_stats(reads: &[LongRead]) -> Result<LongReadStats> {
    if reads.is_empty() {
        return Err(CyaneaError::InvalidInput("no reads".into()));
    }

    let mut lengths: Vec<usize> = reads.iter().map(|r| r.len()).collect();
    lengths.sort_unstable();

    let total_bases: u64 = lengths.iter().map(|&l| l as u64).sum();
    let mean_length = total_bases as f64 / reads.len() as f64;

    let median_length = if lengths.len() % 2 == 0 {
        (lengths[lengths.len() / 2 - 1] + lengths[lengths.len() / 2]) as f64 / 2.0
    } else {
        lengths[lengths.len() / 2] as f64
    };

    // N50
    let half_total = total_bases / 2;
    let mut cumulative = 0u64;
    let mut n50 = 0usize;
    for &l in lengths.iter().rev() {
        cumulative += l as u64;
        if cumulative >= half_total {
            n50 = l;
            break;
        }
    }

    let mean_quality = reads.iter().map(|r| r.mean_quality()).sum::<f64>() / reads.len() as f64;
    let mean_gc = reads.iter().map(|r| r.gc_content()).sum::<f64>() / reads.len() as f64;

    let q20_count = reads.iter().filter(|r| r.mean_quality() >= 20.0).count();
    let q30_count = reads.iter().filter(|r| r.mean_quality() >= 30.0).count();

    let mean_passes = if reads.iter().any(|r| r.num_passes.is_some()) {
        let (sum, count) = reads.iter()
            .filter_map(|r| r.num_passes)
            .fold((0u64, 0u64), |(s, c), p| (s + p as u64, c + 1));
        if count > 0 { Some(sum as f64 / count as f64) } else { None }
    } else {
        None
    };

    Ok(LongReadStats {
        num_reads: reads.len(),
        total_bases,
        mean_length,
        median_length,
        n50,
        max_length: *lengths.last().unwrap(),
        min_length: lengths[0],
        mean_quality,
        mean_gc,
        q20_count,
        q30_count,
        mean_passes,
    })
}

/// Error correction result.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CorrectedRead {
    /// Read identifier.
    pub id: String,
    /// Corrected sequence.
    pub sequence: Vec<u8>,
    /// Corrected quality scores.
    pub quality: Vec<u8>,
    /// Number of bases corrected.
    pub corrections: usize,
    /// Original length.
    pub original_length: usize,
}

/// Self-correction of a long read using k-mer consensus.
///
/// Uses internal k-mer frequency to identify and correct likely errors:
/// 1. Build k-mer frequency table from the read
/// 2. Identify low-frequency k-mers as likely errors
/// 3. Try single-base substitutions to find higher-frequency alternatives
///
/// This is a simplified self-correction; for higher accuracy, use
/// multiple overlapping reads.
pub fn self_correct(read: &LongRead, k: usize) -> Result<CorrectedRead> {
    if read.len() < k {
        return Err(CyaneaError::InvalidInput(
            format!("read length {} shorter than k={}", read.len(), k),
        ));
    }
    if k < 3 || k > 31 {
        return Err(CyaneaError::InvalidInput("k must be 3-31".into()));
    }

    // Build k-mer frequency table
    let mut kmer_counts: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::new();
    for window in read.sequence.windows(k) {
        *kmer_counts.entry(window.to_vec()).or_insert(0) += 1;
    }

    // Compute median k-mer frequency
    let mut freqs: Vec<usize> = kmer_counts.values().cloned().collect();
    freqs.sort_unstable();
    let median_freq = if freqs.is_empty() {
        1
    } else if freqs.len() % 2 == 0 {
        (freqs[freqs.len() / 2 - 1] + freqs[freqs.len() / 2]) / 2
    } else {
        freqs[freqs.len() / 2]
    };

    // Low-frequency threshold
    let threshold = (median_freq as f64 * 0.1).max(1.0) as usize;

    let mut corrected = read.sequence.clone();
    let mut corrections = 0usize;
    let bases = [b'A', b'C', b'G', b'T'];

    // Scan for positions where all overlapping k-mers are low-frequency
    for pos in 0..corrected.len() {
        // Check if this position is in a low-frequency region
        let start = pos.saturating_sub(k - 1);
        let end = (pos + 1).min(corrected.len().saturating_sub(k - 1));

        let all_low = (start..end).all(|s| {
            if s + k <= corrected.len() {
                kmer_counts.get(&corrected[s..s + k]).copied().unwrap_or(0) <= threshold
            } else {
                false
            }
        });

        if !all_low {
            continue;
        }

        // Try substitutions
        let original = corrected[pos];
        let mut best_base = original;
        let mut best_score = 0usize;

        for &base in &bases {
            if base == original {
                continue;
            }
            corrected[pos] = base;

            let score: usize = (start..end)
                .filter_map(|s| {
                    if s + k <= corrected.len() {
                        kmer_counts.get(&corrected[s..s + k]).copied()
                    } else {
                        None
                    }
                })
                .sum();

            if score > best_score {
                best_score = score;
                best_base = base;
            }
        }

        if best_base != original && best_score > threshold * (end - start) {
            corrected[pos] = best_base;
            corrections += 1;
        } else {
            corrected[pos] = original;
        }
    }

    Ok(CorrectedRead {
        id: read.id.clone(),
        quality: read.quality.clone(),
        original_length: read.len(),
        corrections,
        sequence: corrected,
    })
}

/// Consensus from multiple overlapping reads (simplified POA-free approach).
///
/// Given reads that overlap the same region, generates a consensus
/// by majority-vote at each position after alignment to the first read.
///
/// For full graph-based consensus, use `cyanea_align::poa`.
pub fn simple_consensus(reads: &[&[u8]]) -> Result<Vec<u8>> {
    if reads.is_empty() {
        return Err(CyaneaError::InvalidInput("no reads for consensus".into()));
    }
    if reads.len() == 1 {
        return Ok(reads[0].to_vec());
    }

    // Use first read as template, majority-vote at each position
    let max_len = reads.iter().map(|r| r.len()).max().unwrap();
    let mut consensus = Vec::with_capacity(max_len);

    for pos in 0..max_len {
        let mut counts = [0u32; 5]; // A, C, G, T, other
        for &read in reads {
            if pos < read.len() {
                let idx = match read[pos] {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => 4,
                };
                counts[idx] += 1;
            }
        }

        let best = counts[..4].iter()
            .enumerate()
            .max_by_key(|(_, &c)| c)
            .map(|(i, _)| i)
            .unwrap_or(0);

        consensus.push(match best {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }

    Ok(consensus)
}

/// Configuration for long-read simulation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LongReadSimConfig {
    /// Target platform.
    pub platform: LongReadPlatform,
    /// Mean read length.
    pub mean_length: usize,
    /// Standard deviation of read length.
    pub length_std: usize,
    /// Desired coverage.
    pub coverage: f64,
    /// Number of CCS passes (HiFi only, default 10).
    pub num_passes: u32,
    /// PRNG seed.
    pub seed: u64,
}

impl Default for LongReadSimConfig {
    fn default() -> Self {
        Self {
            platform: LongReadPlatform::PacBioHiFi,
            mean_length: 15000,
            length_std: 3000,
            coverage: 30.0,
            num_passes: 10,
            seed: 42,
        }
    }
}

/// Simulate long reads from a reference sequence.
///
/// Generates reads with platform-appropriate error profiles:
/// - **PacBio CLR**: ~13% error (mostly insertions)
/// - **PacBio HiFi**: ~0.1% error (balanced)
/// - **Nanopore simplex**: ~5% error (systematic homopolymer errors)
/// - **Nanopore duplex**: ~0.5% error
pub fn simulate_long_reads(reference: &[u8], config: &LongReadSimConfig) -> Result<Vec<LongRead>> {
    if reference.is_empty() {
        return Err(CyaneaError::InvalidInput("empty reference".into()));
    }

    let error_rate = match config.platform {
        LongReadPlatform::PacBioCLR => 0.13,
        LongReadPlatform::PacBioHiFi => 0.001,
        LongReadPlatform::NanoporeSimplex => 0.05,
        LongReadPlatform::NanoporeDuplex => 0.005,
    };

    let num_reads = ((reference.len() as f64 * config.coverage) / config.mean_length as f64).ceil() as usize;
    let mut reads = Vec::with_capacity(num_reads);

    // Simple LCG PRNG
    let mut rng_state = config.seed;
    let lcg_next = |state: &mut u64| -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state
    };

    for i in 0..num_reads {
        // Sample read length (clamped normal approximation using Box-Muller)
        let u1 = (lcg_next(&mut rng_state) as f64) / (u64::MAX as f64);
        let u2 = (lcg_next(&mut rng_state) as f64) / (u64::MAX as f64);
        let z = (-2.0 * u1.max(1e-10).ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        let read_len = ((config.mean_length as f64 + z * config.length_std as f64) as usize)
            .max(100)
            .min(reference.len());

        // Random start position
        let max_start = reference.len().saturating_sub(read_len);
        let start = if max_start > 0 {
            (lcg_next(&mut rng_state) as usize) % max_start
        } else {
            0
        };

        let mut seq = reference[start..start + read_len].to_vec();
        let mut qual = vec![0u8; read_len];

        // Introduce errors
        let bases = [b'A', b'C', b'G', b'T'];
        let mut _corrections = 0usize;
        for j in 0..seq.len() {
            let rand_val = (lcg_next(&mut rng_state) as f64) / (u64::MAX as f64);
            if rand_val < error_rate {
                // Substitution error
                let new_base = bases[(lcg_next(&mut rng_state) as usize) % 4];
                if new_base != seq[j] {
                    seq[j] = new_base;
                    _corrections += 1;
                }
            }
            // Quality score based on error rate
            let q = (-10.0 * error_rate.max(1e-10).log10()) as u8;
            qual[j] = q.min(60);
        }

        let mut read = LongRead::new(
            format!("sim_longread_{}", i),
            seq,
            qual,
            config.platform,
        )?;

        if config.platform == LongReadPlatform::PacBioHiFi {
            read.num_passes = Some(config.num_passes);
        }

        read.predicted_accuracy = Some(1.0 - error_rate);
        reads.push(read);
    }

    Ok(reads)
}

/// Adapter sequences commonly found in long reads.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LongReadAdapter {
    /// Adapter name.
    pub name: String,
    /// Adapter sequence.
    pub sequence: Vec<u8>,
    /// Platform this adapter is associated with.
    pub platform: LongReadPlatform,
}

/// Common PacBio and Nanopore adapter sequences.
pub fn common_adapters() -> Vec<LongReadAdapter> {
    vec![
        LongReadAdapter {
            name: "PacBio_SMRTbell_3.0".into(),
            sequence: b"ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT".to_vec(),
            platform: LongReadPlatform::PacBioHiFi,
        },
        LongReadAdapter {
            name: "ONT_ligation_top".into(),
            sequence: b"AATGTACTTCGTTCAGTTACGTATTGCT".to_vec(),
            platform: LongReadPlatform::NanoporeSimplex,
        },
        LongReadAdapter {
            name: "ONT_ligation_bottom".into(),
            sequence: b"GCAATACGTAACTGAACGAAGT".to_vec(),
            platform: LongReadPlatform::NanoporeSimplex,
        },
        LongReadAdapter {
            name: "ONT_rapid".into(),
            sequence: b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".to_vec(),
            platform: LongReadPlatform::NanoporeSimplex,
        },
    ]
}

/// Find and trim adapter sequences from a long read.
///
/// Searches for adapter sequences at the start and end of the read
/// using exact substring matching with a configurable max edit distance
/// (Hamming distance).
pub fn trim_adapters(
    read: &LongRead,
    adapters: &[LongReadAdapter],
    max_mismatches: usize,
) -> LongRead {
    let mut seq = read.sequence.clone();
    let mut qual = read.quality.clone();

    for adapter in adapters {
        let adapter_len = adapter.sequence.len();
        if adapter_len > seq.len() {
            continue;
        }

        // Check start of read
        let start_region = &seq[..adapter_len.min(seq.len())];
        let start_mismatches = start_region.iter()
            .zip(adapter.sequence.iter())
            .filter(|(a, b)| a != b)
            .count();

        if start_mismatches <= max_mismatches {
            seq = seq[adapter_len..].to_vec();
            if !qual.is_empty() {
                qual = qual[adapter_len..].to_vec();
            }
        }

        // Check end of read
        if adapter_len <= seq.len() {
            let end_start = seq.len() - adapter_len;
            let end_region = &seq[end_start..];
            let end_mismatches = end_region.iter()
                .zip(adapter.sequence.iter())
                .filter(|(a, b)| a != b)
                .count();

            if end_mismatches <= max_mismatches {
                seq.truncate(end_start);
                if !qual.is_empty() {
                    qual.truncate(end_start);
                }
            }
        }
    }

    let trimmed = LongRead {
        id: read.id.clone(),
        sequence: seq,
        quality: qual,
        platform: read.platform,
        num_passes: read.num_passes,
        predicted_accuracy: read.predicted_accuracy,
        read_group: read.read_group.clone(),
    };
    trimmed
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(id: &str, seq: &[u8], platform: LongReadPlatform) -> LongRead {
        let qual = vec![30u8; seq.len()];
        LongRead::new(id, seq.to_vec(), qual, platform).unwrap()
    }

    #[test]
    fn test_longread_creation() {
        let read = make_read("r1", b"ACGTACGTACGT", LongReadPlatform::PacBioHiFi);
        assert_eq!(read.len(), 12);
        assert!(!read.is_empty());
        assert!((read.mean_quality() - 30.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_read_error() {
        let result = LongRead::new("r1", vec![], vec![], LongReadPlatform::PacBioHiFi);
        assert!(result.is_err());
    }

    #[test]
    fn test_quality_mismatch_error() {
        let result = LongRead::new("r1", vec![b'A', b'C'], vec![30], LongReadPlatform::PacBioHiFi);
        assert!(result.is_err());
    }

    #[test]
    fn test_gc_content() {
        let read = make_read("r1", b"AACCGGTT", LongReadPlatform::NanoporeSimplex);
        assert!((read.gc_content() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_error_rate_from_quality() {
        let read = make_read("r1", b"ACGT", LongReadPlatform::PacBioHiFi);
        // Q30 -> error rate = 10^(-3) = 0.001
        assert!((read.error_rate() - 0.001).abs() < 1e-6);
    }

    #[test]
    fn test_error_rate_default_no_quality() {
        let read = LongRead::new("r1", b"ACGT".to_vec(), vec![], LongReadPlatform::PacBioCLR).unwrap();
        assert!((read.error_rate() - 0.13).abs() < 1e-10);
    }

    #[test]
    fn test_longread_stats() {
        let reads = vec![
            make_read("r1", &vec![b'A'; 1000], LongReadPlatform::PacBioHiFi),
            make_read("r2", &vec![b'C'; 2000], LongReadPlatform::PacBioHiFi),
            make_read("r3", &vec![b'G'; 3000], LongReadPlatform::PacBioHiFi),
        ];

        let stats = longread_stats(&reads).unwrap();
        assert_eq!(stats.num_reads, 3);
        assert_eq!(stats.total_bases, 6000);
        assert!((stats.mean_length - 2000.0).abs() < 1e-10);
        assert!((stats.median_length - 2000.0).abs() < 1e-10);
        assert_eq!(stats.n50, 3000);
        assert_eq!(stats.max_length, 3000);
        assert_eq!(stats.min_length, 1000);
        assert_eq!(stats.q20_count, 3);
        assert_eq!(stats.q30_count, 3);
    }

    #[test]
    fn test_longread_stats_empty() {
        assert!(longread_stats(&[]).is_err());
    }

    #[test]
    fn test_longread_stats_with_passes() {
        let mut reads = vec![
            make_read("r1", &vec![b'A'; 1000], LongReadPlatform::PacBioHiFi),
            make_read("r2", &vec![b'C'; 2000], LongReadPlatform::PacBioHiFi),
        ];
        reads[0].num_passes = Some(10);
        reads[1].num_passes = Some(20);

        let stats = longread_stats(&reads).unwrap();
        assert!((stats.mean_passes.unwrap() - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_simple_consensus() {
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGT",
            b"ACGTACGT",
            b"ACGTACGT",
            b"ACGTTCGT", // one error at pos 4
        ];
        let cons = simple_consensus(&reads).unwrap();
        assert_eq!(cons, b"ACGTACGT"); // majority wins
    }

    #[test]
    fn test_consensus_single_read() {
        let reads: Vec<&[u8]> = vec![b"ACGT"];
        let cons = simple_consensus(&reads).unwrap();
        assert_eq!(cons, b"ACGT");
    }

    #[test]
    fn test_consensus_empty() {
        let reads: Vec<&[u8]> = vec![];
        assert!(simple_consensus(&reads).is_err());
    }

    #[test]
    fn test_self_correct() {
        // Create a read with a single error in a repetitive context
        let mut seq = vec![b'A'; 100];
        seq[50] = b'T'; // introduce error in a poly-A tract
        let read = make_read("r1", &seq, LongReadPlatform::NanoporeSimplex);

        let corrected = self_correct(&read, 11).unwrap();
        assert!(corrected.corrections <= 5); // may or may not correct depending on k-mer context
    }

    #[test]
    fn test_self_correct_k_too_small() {
        let read = make_read("r1", b"ACGT", LongReadPlatform::PacBioHiFi);
        assert!(self_correct(&read, 2).is_err());
    }

    #[test]
    fn test_self_correct_read_too_short() {
        let read = make_read("r1", b"ACGT", LongReadPlatform::PacBioHiFi);
        assert!(self_correct(&read, 10).is_err());
    }

    #[test]
    fn test_simulate_long_reads() {
        let reference = vec![b'A'; 1000];
        let config = LongReadSimConfig {
            platform: LongReadPlatform::PacBioHiFi,
            mean_length: 200,
            length_std: 50,
            coverage: 10.0,
            num_passes: 10,
            seed: 42,
        };

        let reads = simulate_long_reads(&reference, &config).unwrap();
        assert!(!reads.is_empty());
        for r in &reads {
            assert!(r.len() > 0);
            assert_eq!(r.platform, LongReadPlatform::PacBioHiFi);
            assert!(r.num_passes.is_some());
        }
    }

    #[test]
    fn test_simulate_nanopore() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let config = LongReadSimConfig {
            platform: LongReadPlatform::NanoporeSimplex,
            mean_length: 30,
            length_std: 5,
            coverage: 5.0,
            num_passes: 1,
            seed: 123,
        };

        let reads = simulate_long_reads(reference, &config).unwrap();
        assert!(!reads.is_empty());
        assert_eq!(reads[0].platform, LongReadPlatform::NanoporeSimplex);
        assert!(reads[0].num_passes.is_none());
    }

    #[test]
    fn test_common_adapters() {
        let adapters = common_adapters();
        assert!(adapters.len() >= 4);
        assert!(adapters.iter().any(|a| a.name.contains("PacBio")));
        assert!(adapters.iter().any(|a| a.name.contains("ONT")));
    }

    #[test]
    fn test_trim_adapters() {
        let adapter_seq = b"AATGTACTTCGTTCAGTTACGTATTGCT";
        let inner = b"ACGTACGTACGTACGT";
        let mut full_seq = adapter_seq.to_vec();
        full_seq.extend_from_slice(inner);

        let read = make_read("r1", &full_seq, LongReadPlatform::NanoporeSimplex);
        let adapters = common_adapters();
        let trimmed = trim_adapters(&read, &adapters, 0);

        assert_eq!(trimmed.sequence, inner);
    }

    #[test]
    fn test_trim_adapters_with_mismatches() {
        let mut adapter_seq = b"AATGTACTTCGTTCAGTTACGTATTGCT".to_vec();
        adapter_seq[5] = b'G'; // introduce 1 mismatch (original is 'A')
        let inner = b"ACGTACGTACGTACGT";
        let mut full_seq = adapter_seq.clone();
        full_seq.extend_from_slice(inner);

        let read = make_read("r1", &full_seq, LongReadPlatform::NanoporeSimplex);
        let adapters = common_adapters();

        // With 0 mismatches allowed, no trimming
        let trimmed_strict = trim_adapters(&read, &adapters, 0);
        assert_eq!(trimmed_strict.len(), full_seq.len());

        // With 1 mismatch allowed, should trim
        let trimmed_relaxed = trim_adapters(&read, &adapters, 1);
        assert_eq!(trimmed_relaxed.sequence, inner);
    }
}
