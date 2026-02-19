//! Illumina-style short-read simulator.
//!
//! Generates synthetic paired-end or single-end reads from a reference
//! sequence, with configurable coverage, fragment sizes, and a
//! position-dependent error profile that mimics Illumina 3' quality
//! degradation.

use std::f64::consts::PI;

/// Configuration for the read simulator.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ReadSimConfig {
    /// Read length in bases (default 150).
    pub read_length: usize,
    /// Desired sequencing coverage (default 30.0).
    pub coverage: f64,
    /// Mean insert/fragment size in bases (default 300.0).
    pub fragment_mean: f64,
    /// Standard deviation of fragment size (default 50.0).
    pub fragment_std: f64,
    /// Per-base error rate at position 0 (default 0.001).
    pub error_rate: f64,
    /// Whether to produce paired-end reads (default true).
    pub paired: bool,
    /// PRNG seed for reproducibility.
    pub seed: u64,
}

impl Default for ReadSimConfig {
    fn default() -> Self {
        Self {
            read_length: 150,
            coverage: 30.0,
            fragment_mean: 300.0,
            fragment_std: 50.0,
            error_rate: 0.001,
            paired: true,
            seed: 42,
        }
    }
}

/// A single simulated read with its true genomic origin.
#[derive(Debug, Clone)]
pub struct SimulatedRead {
    /// Read name (e.g. "sim_read_0/1").
    pub name: String,
    /// Nucleotide sequence (A, C, G, T) with introduced errors.
    pub sequence: Vec<u8>,
    /// Phred+33 encoded quality scores.
    pub quality: Vec<u8>,
    /// True 0-based start position on the reference.
    pub true_position: u64,
    /// Chromosome / contig name.
    pub true_chrom: String,
    /// Whether this is read 1 of a pair (false for read 2 or single-end).
    pub is_read1: bool,
}

// ---------------------------------------------------------------------------
// Private Xorshift64 PRNG
// ---------------------------------------------------------------------------

struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Avoid the all-zero fixed point.
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        self.next_u64() as f64 / u64::MAX as f64
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Box-Muller transform: produce a standard-normal sample from two uniforms.
fn box_muller(rng: &mut Xorshift64) -> f64 {
    let u1 = rng.next_f64().max(f64::MIN_POSITIVE); // avoid ln(0)
    let u2 = rng.next_f64();
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// Reverse complement of a byte slice (A<->T, C<->G).
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Pick a random base that differs from `original`.
fn substitute_base(rng: &mut Xorshift64, original: u8) -> u8 {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    loop {
        let b = BASES[(rng.next_u64() % 4) as usize];
        if b != original {
            return b;
        }
    }
}

/// Compute Phred+33 quality byte for a given error probability.
fn phred_quality(error_prob: f64) -> u8 {
    if error_prob <= 0.0 {
        // Perfect quality — cap at Q40.
        return 73; // 40 + 33
    }
    let q = (-10.0 * error_prob.log10()).round() as i32;
    let clamped = q.clamp(2, 40) as u8; // Q2..Q40
    clamped + 33 // Phred+33 encoding: 35..=73
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Simulate Illumina-style reads from a reference sequence.
///
/// Returns an empty `Vec` when the reference is shorter than `read_length`.
///
/// # Arguments
///
/// * `reference` - The reference sequence as uppercase ACGT bytes.
/// * `chrom_name` - Name to tag on every read's `true_chrom` field.
/// * `config` - Simulation parameters (coverage, error model, etc.).
pub fn simulate_reads(
    reference: &[u8],
    chrom_name: &str,
    config: &ReadSimConfig,
) -> Vec<SimulatedRead> {
    let ref_len = reference.len();
    if ref_len < config.read_length {
        return Vec::new();
    }

    let mut rng = Xorshift64::new(config.seed);

    // Total reads needed to achieve the requested coverage.
    let n_reads =
        ((ref_len as f64 * config.coverage) / config.read_length as f64).round() as usize;

    // Pre-compute per-position quality scores.
    let qualities: Vec<u8> = (0..config.read_length)
        .map(|i| {
            let ep = config.error_rate * (1.0 + 2.0 * i as f64 / config.read_length as f64);
            phred_quality(ep)
        })
        .collect();

    let mut reads = Vec::new();

    if config.paired {
        let n_fragments = n_reads / 2;
        for frag_idx in 0..n_fragments {
            // Sample fragment size from N(mean, std), clamped.
            let frag_size = {
                let z = box_muller(&mut rng);
                let raw = config.fragment_mean + config.fragment_std * z;
                raw.round().clamp(config.read_length as f64, ref_len as f64) as usize
            };

            // Random start position.
            let max_start = if ref_len > frag_size {
                ref_len - frag_size
            } else {
                0
            };
            let start = if max_start > 0 {
                (rng.next_u64() % max_start as u64) as usize
            } else {
                0
            };

            let fragment = &reference[start..start + frag_size];

            // Read 1: first read_length bases of fragment.
            let r1_bases: Vec<u8> = fragment[..config.read_length].to_vec();

            // Read 2: last read_length bases of fragment, reverse complemented.
            let r2_fwd = &fragment[frag_size - config.read_length..];
            let r2_bases = reverse_complement(r2_fwd);

            let r1_pos = start as u64;
            let r2_pos = (start + frag_size - config.read_length) as u64;

            // Introduce errors and build reads.
            let r1 = build_read(
                &mut rng,
                &r1_bases,
                &qualities,
                config,
                frag_idx,
                true,
                r1_pos,
                chrom_name,
            );
            let r2 = build_read(
                &mut rng,
                &r2_bases,
                &qualities,
                config,
                frag_idx,
                false,
                r2_pos,
                chrom_name,
            );
            reads.push(r1);
            reads.push(r2);
        }
    } else {
        for frag_idx in 0..n_reads {
            let max_start = ref_len - config.read_length;
            let start = if max_start > 0 {
                (rng.next_u64() % max_start as u64) as usize
            } else {
                0
            };
            let bases: Vec<u8> = reference[start..start + config.read_length].to_vec();
            let r = build_read(
                &mut rng,
                &bases,
                &qualities,
                config,
                frag_idx,
                true,
                start as u64,
                chrom_name,
            );
            reads.push(r);
        }
    }

    reads
}

/// Apply position-dependent errors and package into a [`SimulatedRead`].
fn build_read(
    rng: &mut Xorshift64,
    bases: &[u8],
    qualities: &[u8],
    config: &ReadSimConfig,
    frag_idx: usize,
    is_read1: bool,
    true_position: u64,
    chrom_name: &str,
) -> SimulatedRead {
    let read_len = config.read_length;
    let mut seq = bases.to_vec();

    for i in 0..read_len {
        let error_prob =
            config.error_rate * (1.0 + 2.0 * i as f64 / read_len as f64);
        if rng.next_f64() < error_prob {
            seq[i] = substitute_base(rng, seq[i]);
        }
    }

    let suffix = if is_read1 { "/1" } else { "/2" };
    SimulatedRead {
        name: format!("sim_read_{}{}", frag_idx, suffix),
        sequence: seq,
        quality: qualities.to_vec(),
        true_position,
        true_chrom: chrom_name.to_string(),
        is_read1,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: a 1 kb reference of repeating ACGT.
    fn reference_1kb() -> Vec<u8> {
        let unit = b"ACGTACGTACGTACGT"; // 16 bp
        unit.iter().copied().cycle().take(1000).collect()
    }

    #[test]
    fn read_count_proportional_to_coverage() {
        let reference = reference_1kb();
        let config = ReadSimConfig {
            read_length: 100,
            coverage: 10.0,
            paired: false,
            seed: 123,
            ..Default::default()
        };
        let reads = simulate_reads(&reference, "chr1", &config);
        // Expected: 1000 * 10.0 / 100 = 100 reads.
        let expected = ((reference.len() as f64 * config.coverage) / config.read_length as f64)
            .round() as usize;
        assert_eq!(reads.len(), expected);
    }

    #[test]
    fn read_length_matches_config() {
        let reference = reference_1kb();
        let config = ReadSimConfig {
            read_length: 75,
            coverage: 5.0,
            paired: true,
            seed: 456,
            ..Default::default()
        };
        let reads = simulate_reads(&reference, "chr1", &config);
        assert!(!reads.is_empty());
        for r in &reads {
            assert_eq!(r.sequence.len(), 75);
            assert_eq!(r.quality.len(), 75);
        }
    }

    #[test]
    fn paired_produces_both_mates() {
        let reference = reference_1kb();
        let config = ReadSimConfig {
            read_length: 100,
            coverage: 10.0,
            paired: true,
            seed: 789,
            ..Default::default()
        };
        let reads = simulate_reads(&reference, "chrX", &config);
        assert!(!reads.is_empty());
        // Paired reads come in pairs (read1, read2, read1, read2, ...).
        assert_eq!(reads.len() % 2, 0);
        let n_read1 = reads.iter().filter(|r| r.is_read1).count();
        let n_read2 = reads.iter().filter(|r| !r.is_read1).count();
        assert_eq!(n_read1, n_read2);
        assert!(n_read1 > 0);
        // Check naming convention.
        assert!(reads[0].name.ends_with("/1"));
        assert!(reads[1].name.ends_with("/2"));
    }

    #[test]
    fn reads_are_valid_dna() {
        let reference = reference_1kb();
        let config = ReadSimConfig {
            coverage: 5.0,
            seed: 101,
            ..Default::default()
        };
        let reads = simulate_reads(&reference, "chr1", &config);
        for r in &reads {
            for &b in &r.sequence {
                assert!(
                    b == b'A' || b == b'C' || b == b'G' || b == b'T',
                    "unexpected base: {}",
                    b as char
                );
            }
        }
    }

    #[test]
    fn quality_scores_valid_phred33() {
        let reference = reference_1kb();
        let config = ReadSimConfig {
            coverage: 5.0,
            seed: 202,
            ..Default::default()
        };
        let reads = simulate_reads(&reference, "chr1", &config);
        for r in &reads {
            for &q in &r.quality {
                assert!(
                    (35..=73).contains(&q),
                    "quality byte {} out of Phred+33 range 35..=73",
                    q
                );
            }
        }
    }

    #[test]
    fn zero_error_rate_reads_match_reference() {
        let reference = reference_1kb();
        let config = ReadSimConfig {
            read_length: 100,
            coverage: 10.0,
            error_rate: 0.0,
            paired: false,
            seed: 303,
            ..Default::default()
        };
        let reads = simulate_reads(&reference, "chr1", &config);
        for r in &reads {
            let pos = r.true_position as usize;
            let expected = &reference[pos..pos + 100];
            assert_eq!(
                r.sequence, expected,
                "read at position {} does not match reference",
                pos
            );
        }
    }
}
