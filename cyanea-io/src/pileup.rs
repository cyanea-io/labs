//! Pileup generation from SAM/BAM alignment records.
//!
//! Transforms [`SamRecord`] data into per-position pileup columns showing
//! base counts, quality sums, and depth-of-coverage. Supports mpileup-format
//! output and a simple SNP caller.

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

use crate::sam::SamRecord;

// ---------------------------------------------------------------------------
// Internal: CIGAR parsing
// ---------------------------------------------------------------------------

/// Minimal CIGAR operations for pileup walking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CigarOp {
    /// M/=/X: consumes ref + query.
    Align(usize),
    /// I: consumes query only.
    Ins(usize),
    /// D/N: consumes ref only.
    Del(usize),
    /// S: consumes query only.
    SoftClip(usize),
    /// H/P: consumes neither.
    Skip(usize),
}

/// Parse a CIGAR string into a sequence of operations.
fn parse_cigar_ops(cigar: &str) -> Result<Vec<CigarOp>> {
    if cigar == "*" {
        return Ok(Vec::new());
    }

    let mut ops = Vec::new();
    let mut num_start = 0;

    for (i, c) in cigar.char_indices() {
        if c.is_ascii_digit() {
            continue;
        }
        let len: usize = cigar[num_start..i]
            .parse()
            .map_err(|_| CyaneaError::Parse(format!("invalid CIGAR length in '{cigar}'")))?;
        let op = match c {
            'M' | '=' | 'X' => CigarOp::Align(len),
            'I' => CigarOp::Ins(len),
            'D' | 'N' => CigarOp::Del(len),
            'S' => CigarOp::SoftClip(len),
            'H' | 'P' => CigarOp::Skip(len),
            _ => {
                return Err(CyaneaError::Parse(format!(
                    "unknown CIGAR operation '{c}' in '{cigar}'"
                )));
            }
        };
        ops.push(op);
        num_start = i + 1;
    }

    Ok(ops)
}

/// Map a base to an index: A=0, C=1, G=2, T=3, N=4, del=5.
fn base_index(base: u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'N' | b'n' => 4,
        _ => 5, // deletion marker or unknown
    }
}

/// FLAG bit 0x10: read is on the reverse strand.
fn is_reverse(flag: u16) -> bool {
    flag & 0x10 != 0
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Per-position pileup data.
#[derive(Debug, Clone)]
pub struct PileupColumn {
    /// Reference sequence name.
    pub rname: String,
    /// 0-based position on the reference.
    pub pos: u64,
    /// Reference base at this position (N if no reference provided).
    pub ref_base: u8,
    /// Depth of coverage.
    pub depth: u32,
    /// Observed bases (upper-case = forward strand, lower-case = reverse strand).
    pub bases: Vec<u8>,
    /// Phred quality scores corresponding to each base in `bases`.
    pub qualities: Vec<u8>,
    /// Counts per base type: [A, C, G, T, N, del].
    pub base_counts: [u32; 6],
    /// Sum of Phred quality scores per base type: [A, C, G, T, N, del].
    pub quality_sums: [u64; 6],
}

/// Pileup for a single reference sequence.
#[derive(Debug, Clone)]
pub struct Pileup {
    /// Reference sequence name.
    pub rname: String,
    /// Sorted columns (by position).
    pub columns: Vec<PileupColumn>,
}

/// Depth-of-coverage statistics for a pileup.
#[derive(Debug, Clone)]
pub struct DepthStats {
    /// Reference sequence name.
    pub rname: String,
    /// Length of the covered region (max_pos - min_pos + 1, or 0 if empty).
    pub length: u64,
    /// Number of positions with depth > 0.
    pub covered: u64,
    /// Breadth of coverage (covered / length).
    pub breadth: f64,
    /// Minimum depth across all positions.
    pub min_depth: u32,
    /// Maximum depth across all positions.
    pub max_depth: u32,
    /// Mean depth across all positions.
    pub mean_depth: f64,
    /// Median depth across all positions.
    pub median_depth: f64,
}

// ---------------------------------------------------------------------------
// Internal: column builder
// ---------------------------------------------------------------------------

/// Builder for accumulating bases at a single position.
struct ColumnBuilder {
    bases: Vec<u8>,
    qualities: Vec<u8>,
    base_counts: [u32; 6],
    quality_sums: [u64; 6],
}

impl ColumnBuilder {
    fn new() -> Self {
        Self {
            bases: Vec::new(),
            qualities: Vec::new(),
            base_counts: [0; 6],
            quality_sums: [0; 6],
        }
    }

    fn add(&mut self, base: u8, qual: u8) {
        self.bases.push(base);
        self.qualities.push(qual);
        let idx = base_index(base);
        self.base_counts[idx] += 1;
        self.quality_sums[idx] += qual as u64;
    }

    fn into_column(self, rname: &str, pos: u64, ref_base: u8) -> PileupColumn {
        PileupColumn {
            rname: rname.to_string(),
            pos,
            ref_base,
            depth: self.bases.len() as u32,
            bases: self.bases,
            qualities: self.qualities,
            base_counts: self.base_counts,
            quality_sums: self.quality_sums,
        }
    }
}

// ---------------------------------------------------------------------------
// Core: walk a single read into the column map
// ---------------------------------------------------------------------------

/// Walk a single alignment record and deposit bases/qualities into columns.
fn walk_read(
    record: &SamRecord,
    columns: &mut HashMap<u64, ColumnBuilder>,
) -> Result<()> {
    let ops = parse_cigar_ops(&record.cigar)?;
    if ops.is_empty() {
        return Ok(());
    }

    // SAM POS is 1-based; convert to 0-based.
    let mut ref_pos = record.pos.saturating_sub(1);
    let mut query_pos: usize = 0;

    let seq = record.sequence.as_bytes();
    let qual = record.quality.as_bytes();
    let reverse = is_reverse(record.flag);

    for op in &ops {
        match *op {
            CigarOp::Align(len) => {
                for _ in 0..len {
                    let base = if query_pos < seq.len() {
                        let b = seq[query_pos];
                        if reverse {
                            b.to_ascii_lowercase()
                        } else {
                            b.to_ascii_uppercase()
                        }
                    } else if reverse {
                        b'n'
                    } else {
                        b'N'
                    };

                    let q = if query_pos < qual.len() && qual != b"*" {
                        qual[query_pos].saturating_sub(33)
                    } else {
                        0
                    };

                    columns.entry(ref_pos).or_insert_with(ColumnBuilder::new).add(base, q);
                    ref_pos += 1;
                    query_pos += 1;
                }
            }
            CigarOp::Ins(len) => {
                // Insertions consume query only — no ref positions covered.
                query_pos += len;
            }
            CigarOp::Del(len) => {
                // Deletions consume ref — deposit deletion markers.
                for _ in 0..len {
                    let del_base = if reverse { b'#' } else { b'*' };
                    columns.entry(ref_pos).or_insert_with(ColumnBuilder::new).add(del_base, 0);
                    ref_pos += 1;
                }
            }
            CigarOp::SoftClip(len) => {
                // Soft clips consume query only.
                query_pos += len;
            }
            CigarOp::Skip(len) => {
                // H/P: consume neither ref nor query. Just a no-op for length.
                let _ = len;
            }
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Generate pileup from SAM records, grouped by reference.
///
/// Records are grouped by `rname`. Within each group, CIGAR strings are
/// walked to produce per-position base and quality data.
///
/// If `reference` is provided, it maps `rname` to a reference sequence
/// so that `ref_base` in each [`PileupColumn`] reflects the true reference
/// base. Otherwise, `ref_base` is set to `b'N'`.
///
/// Unmapped records (FLAG 0x4) and records with `*` CIGAR are skipped.
pub fn pileup(
    records: &[SamRecord],
    reference: Option<&HashMap<String, Vec<u8>>>,
) -> Result<Vec<Pileup>> {
    // Group records by rname.
    let mut groups: HashMap<String, Vec<&SamRecord>> = HashMap::new();
    for rec in records {
        if rec.is_unmapped() || rec.cigar == "*" {
            continue;
        }
        groups.entry(rec.rname.clone()).or_default().push(rec);
    }

    let mut pileups: Vec<Pileup> = Vec::new();
    let mut rnames: Vec<String> = groups.keys().cloned().collect();
    rnames.sort();

    for rname in rnames {
        let recs = &groups[&rname];
        let ref_seq = reference.and_then(|r| r.get(&rname));
        let p = build_pileup(&rname, recs, ref_seq.map(|v| v.as_slice()))?;
        pileups.push(p);
    }

    Ok(pileups)
}

/// Generate pileup for a single reference sequence.
///
/// Only records mapping to `rname` are included. Unmapped records and
/// records with `*` CIGAR are skipped.
pub fn pileup_region(
    records: &[SamRecord],
    rname: &str,
    reference: Option<&[u8]>,
) -> Result<Pileup> {
    let filtered: Vec<&SamRecord> = records
        .iter()
        .filter(|r| r.rname == rname && r.is_mapped() && r.cigar != "*")
        .collect();

    build_pileup(rname, &filtered, reference)
}

/// Build a Pileup from a set of records on the same reference.
fn build_pileup(
    rname: &str,
    records: &[&SamRecord],
    ref_seq: Option<&[u8]>,
) -> Result<Pileup> {
    let mut columns: HashMap<u64, ColumnBuilder> = HashMap::new();

    for rec in records {
        walk_read(rec, &mut columns)?;
    }

    // Sort by position and build PileupColumn structs.
    let mut positions: Vec<u64> = columns.keys().copied().collect();
    positions.sort_unstable();

    let cols: Vec<PileupColumn> = positions
        .into_iter()
        .map(|pos| {
            let builder = columns.remove(&pos).unwrap();
            let ref_base = ref_seq
                .and_then(|s| s.get(pos as usize).copied())
                .unwrap_or(b'N');
            builder.into_column(rname, pos, ref_base)
        })
        .collect();

    Ok(Pileup {
        rname: rname.to_string(),
        columns: cols,
    })
}

/// Format a pileup as mpileup text.
///
/// Output format (tab-separated, one line per position):
/// `RNAME  POS(1-based)  REF_BASE  DEPTH  BASES  QUALS`
///
/// Base encoding:
/// - `.` / `,` — match to reference on forward / reverse strand
/// - `ACGTN` / `acgtn` — mismatch on forward / reverse strand
/// - `*` — deletion
pub fn pileup_to_mpileup(pileup: &Pileup) -> String {
    let mut lines = Vec::with_capacity(pileup.columns.len());

    for col in &pileup.columns {
        let ref_upper = col.ref_base.to_ascii_uppercase();
        let mut bases_str = String::with_capacity(col.bases.len());
        let mut quals_str = String::with_capacity(col.qualities.len());

        for (i, &base) in col.bases.iter().enumerate() {
            // Deletion markers.
            if base == b'*' || base == b'#' {
                bases_str.push('*');
            } else {
                let base_upper = base.to_ascii_uppercase();
                if base_upper == ref_upper && ref_upper != b'N' {
                    // Match to reference.
                    if base.is_ascii_uppercase() {
                        bases_str.push('.');
                    } else {
                        bases_str.push(',');
                    }
                } else {
                    // Mismatch or unknown ref.
                    bases_str.push(base as char);
                }
            }

            let q = if i < col.qualities.len() {
                col.qualities[i]
            } else {
                0
            };
            quals_str.push((q + 33) as char);
        }

        // mpileup is 1-based.
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            col.rname,
            col.pos + 1,
            ref_upper as char,
            col.depth,
            bases_str,
            quals_str,
        ));
    }

    lines.join("\n")
}

/// Compute depth-of-coverage statistics from a pileup.
pub fn depth_stats(pileup: &Pileup) -> DepthStats {
    if pileup.columns.is_empty() {
        return DepthStats {
            rname: pileup.rname.clone(),
            length: 0,
            covered: 0,
            breadth: 0.0,
            min_depth: 0,
            max_depth: 0,
            mean_depth: 0.0,
            median_depth: 0.0,
        };
    }

    let min_pos = pileup.columns.first().unwrap().pos;
    let max_pos = pileup.columns.last().unwrap().pos;
    let length = max_pos - min_pos + 1;

    // Build depth vector including gaps (0 depth).
    let mut depths: Vec<u32> = vec![0; length as usize];
    for col in &pileup.columns {
        let idx = (col.pos - min_pos) as usize;
        depths[idx] = col.depth;
    }

    let covered = depths.iter().filter(|&&d| d > 0).count() as u64;
    let min_depth = *depths.iter().min().unwrap();
    let max_depth = *depths.iter().max().unwrap();
    let sum: u64 = depths.iter().map(|&d| d as u64).sum();
    let mean_depth = sum as f64 / length as f64;

    let mut sorted = depths.clone();
    sorted.sort_unstable();
    let median_depth = if sorted.len().is_multiple_of(2) {
        let mid = sorted.len() / 2;
        (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0
    } else {
        sorted[sorted.len() / 2] as f64
    };

    let breadth = covered as f64 / length as f64;

    DepthStats {
        rname: pileup.rname.clone(),
        length,
        covered,
        breadth,
        min_depth,
        max_depth,
        mean_depth,
        median_depth,
    }
}

/// Call simple SNPs from a pileup.
///
/// For each position with sufficient depth, finds the most common non-reference,
/// non-deletion base. Emits a [`Variant`] if the alternate allele meets the
/// frequency and count thresholds.
///
/// Quality is estimated as `-10 * log10(1 - alt_freq)`, capped at 99.
#[cfg(feature = "vcf")]
pub fn call_snps(
    pileup: &Pileup,
    min_depth: u32,
    min_alt_freq: f64,
    min_alt_count: u32,
) -> Vec<cyanea_omics::variant::Variant> {
    use cyanea_omics::variant::Variant;

    let mut variants = Vec::new();

    for col in &pileup.columns {
        if col.depth < min_depth {
            continue;
        }

        let ref_upper = col.ref_base.to_ascii_uppercase();
        let ref_idx = base_index(ref_upper);

        // Find most frequent non-ref, non-deletion base.
        let mut best_idx = usize::MAX;
        let mut best_count: u32 = 0;

        for i in 0..5 {
            // A, C, G, T, N (skip deletions at index 5)
            if i == ref_idx {
                continue;
            }
            if col.base_counts[i] > best_count {
                best_count = col.base_counts[i];
                best_idx = i;
            }
        }

        if best_idx == usize::MAX || best_count < min_alt_count {
            continue;
        }

        let alt_freq = best_count as f64 / col.depth as f64;
        if alt_freq < min_alt_freq {
            continue;
        }

        let alt_base = match best_idx {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        };

        // Quality = -10 * log10(1 - alt_freq), capped at 99.
        let qual = if alt_freq >= 1.0 {
            99.0
        } else {
            let raw = -10.0 * (1.0 - alt_freq).log10();
            raw.min(99.0)
        };

        // VCF is 1-based.
        if let Ok(mut v) = Variant::new(
            &col.rname,
            col.pos + 1,
            vec![ref_upper],
            vec![vec![alt_base]],
        ) {
            v.quality = Some(qual);
            variants.push(v);
        }
    }

    variants
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a SamRecord from individual fields for testing.
    fn make_record(
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u64,
        mapq: u8,
        cigar: &str,
        seq: &str,
        qual: &str,
    ) -> SamRecord {
        SamRecord {
            qname: qname.to_string(),
            flag,
            rname: rname.to_string(),
            pos,
            mapq,
            cigar: cigar.to_string(),
            sequence: seq.to_string(),
            quality: qual.to_string(),
        }
    }

    // -- CIGAR parsing tests --

    #[test]
    fn test_cigar_simple_match() {
        let ops = parse_cigar_ops("50M").unwrap();
        assert_eq!(ops, vec![CigarOp::Align(50)]);
    }

    #[test]
    fn test_cigar_complex() {
        let ops = parse_cigar_ops("10M2I3D5M").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::Align(10),
                CigarOp::Ins(2),
                CigarOp::Del(3),
                CigarOp::Align(5),
            ]
        );
    }

    #[test]
    fn test_cigar_star() {
        let ops = parse_cigar_ops("*").unwrap();
        assert!(ops.is_empty());
    }

    #[test]
    fn test_cigar_eq_x_mapped_to_align() {
        let ops = parse_cigar_ops("5=3X2M").unwrap();
        assert_eq!(
            ops,
            vec![CigarOp::Align(5), CigarOp::Align(3), CigarOp::Align(2)]
        );
    }

    #[test]
    fn test_cigar_invalid_char() {
        let result = parse_cigar_ops("10Z");
        assert!(result.is_err());
    }

    // -- Core pileup tests --

    #[test]
    fn test_pileup_single_read() {
        let records = vec![make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "IIII")];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups.len(), 1);
        assert_eq!(pileups[0].rname, "chr1");
        assert_eq!(pileups[0].columns.len(), 4);

        // Position 0: base A
        assert_eq!(pileups[0].columns[0].pos, 0);
        assert_eq!(pileups[0].columns[0].depth, 1);
        assert_eq!(pileups[0].columns[0].bases, vec![b'A']);
        assert_eq!(pileups[0].columns[0].base_counts[0], 1); // A count
    }

    #[test]
    fn test_pileup_two_overlapping_reads() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "IIII"),
            make_record("r2", 0, "chr1", 3, 60, "4M", "TTAA", "IIII"),
        ];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups.len(), 1);
        assert_eq!(pileups[0].columns.len(), 6); // pos 0..5

        // Position 2 (0-based): covered by r1 (G) and r2 (T)
        let col2 = &pileups[0].columns[2];
        assert_eq!(col2.pos, 2);
        assert_eq!(col2.depth, 2);
        assert_eq!(col2.bases, vec![b'G', b'T']);
    }

    #[test]
    fn test_pileup_deletion() {
        // 2M2D2M on sequence ACGT — positions 0,1 get A,C; positions 2,3 get deletions; positions 4,5 get G,T
        let records = vec![make_record("r1", 0, "chr1", 1, 60, "2M2D2M", "ACGT", "IIII")];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups[0].columns.len(), 6);

        // Positions 2 and 3 should have deletion markers.
        assert_eq!(pileups[0].columns[2].bases, vec![b'*']);
        assert_eq!(pileups[0].columns[2].base_counts[5], 1); // del count
        assert_eq!(pileups[0].columns[3].bases, vec![b'*']);
    }

    #[test]
    fn test_pileup_insertion() {
        // 3M2I3M on sequence ACGTTTAA — covers 6 ref positions (not 8)
        let records = vec![make_record(
            "r1", 0, "chr1", 1, 60, "3M2I3M", "ACGTTAAA", "IIIIIIII",
        )];
        let pileups = pileup(&records, None).unwrap();
        // 3 + 3 = 6 ref positions
        assert_eq!(pileups[0].columns.len(), 6);
    }

    #[test]
    fn test_pileup_soft_clip() {
        // 2S3M2S on sequence NNACGNN — only the 3M part covers ref positions
        let records = vec![make_record(
            "r1", 0, "chr1", 1, 60, "2S3M2S", "NNACGNN", "IIIIIII",
        )];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups[0].columns.len(), 3);
        // After 2S, query_pos=2, so bases are A, C, G
        assert_eq!(pileups[0].columns[0].bases, vec![b'A']);
        assert_eq!(pileups[0].columns[1].bases, vec![b'C']);
        assert_eq!(pileups[0].columns[2].bases, vec![b'G']);
    }

    #[test]
    fn test_pileup_unmapped_skipped() {
        let records = vec![
            make_record("r1", 4, "*", 0, 0, "*", "ACGT", "*"), // unmapped
            make_record("r2", 0, "chr1", 1, 60, "4M", "TTTT", "IIII"),
        ];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups.len(), 1);
        assert_eq!(pileups[0].columns[0].depth, 1); // only r2
    }

    #[test]
    fn test_pileup_star_cigar_skipped() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "*", "ACGT", "IIII"),
            make_record("r2", 0, "chr1", 1, 60, "4M", "TTTT", "IIII"),
        ];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups.len(), 1);
        assert_eq!(pileups[0].columns[0].depth, 1);
    }

    #[test]
    fn test_pileup_with_reference() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGT".to_vec());

        let records = vec![make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "IIII")];
        let pileups = pileup(&records, Some(&reference)).unwrap();

        assert_eq!(pileups[0].columns[0].ref_base, b'A');
        assert_eq!(pileups[0].columns[1].ref_base, b'C');
        assert_eq!(pileups[0].columns[2].ref_base, b'G');
        assert_eq!(pileups[0].columns[3].ref_base, b'T');
    }

    // -- Multi-reference tests --

    #[test]
    fn test_pileup_multi_reference() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "AAAA", "IIII"),
            make_record("r2", 0, "chr2", 1, 60, "4M", "CCCC", "IIII"),
        ];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups.len(), 2);
        assert_eq!(pileups[0].rname, "chr1");
        assert_eq!(pileups[1].rname, "chr2");
    }

    #[test]
    fn test_pileup_region_filters() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "AAAA", "IIII"),
            make_record("r2", 0, "chr2", 1, 60, "4M", "CCCC", "IIII"),
            make_record("r3", 0, "chr1", 5, 60, "4M", "GGGG", "IIII"),
        ];
        let p = pileup_region(&records, "chr1", None).unwrap();
        assert_eq!(p.rname, "chr1");
        // r1 covers 0..3, r3 covers 4..7 → 8 positions
        assert_eq!(p.columns.len(), 8);
    }

    // -- mpileup output tests --

    #[test]
    fn test_mpileup_basic_format() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGT".to_vec());

        let records = vec![make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "IIII")];
        let pileups = pileup(&records, Some(&reference)).unwrap();
        let output = pileup_to_mpileup(&pileups[0]);

        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4);
        // Each line: rname, 1-based pos, ref_base, depth, bases, quals
        assert!(lines[0].starts_with("chr1\t1\tA\t1\t"));
    }

    #[test]
    fn test_mpileup_match_dot_comma() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "AAAA", "IIII"),  // fwd
            make_record("r2", 16, "chr1", 1, 60, "4M", "AAAA", "IIII"), // rev
        ];
        let pileups = pileup(&records, Some(&reference)).unwrap();
        let output = pileup_to_mpileup(&pileups[0]);

        let lines: Vec<&str> = output.lines().collect();
        // At each position: fwd A matches ref A → '.', rev a matches ref A → ','
        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields[4], ".,");
    }

    #[test]
    fn test_mpileup_mismatch_case() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "CCCC", "IIII"),  // fwd, mismatch
            make_record("r2", 16, "chr1", 1, 60, "4M", "CCCC", "IIII"), // rev, mismatch
        ];
        let pileups = pileup(&records, Some(&reference)).unwrap();
        let output = pileup_to_mpileup(&pileups[0]);

        let lines: Vec<&str> = output.lines().collect();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        // fwd C mismatch → 'C', rev c mismatch → 'c'
        assert_eq!(fields[4], "Cc");
    }

    // -- Depth stats tests --

    #[test]
    fn test_depth_stats_basic() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "IIII"),
            make_record("r2", 0, "chr1", 3, 60, "4M", "TTAA", "IIII"),
        ];
        let pileups = pileup(&records, None).unwrap();
        let stats = depth_stats(&pileups[0]);

        assert_eq!(stats.rname, "chr1");
        assert_eq!(stats.length, 6); // pos 0..5
        assert_eq!(stats.covered, 6);
        assert_eq!(stats.min_depth, 1);
        assert_eq!(stats.max_depth, 2);
        assert!((stats.breadth - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_depth_stats_gaps() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "2M", "AC", "II"),
            make_record("r2", 0, "chr1", 5, 60, "2M", "GT", "II"),
        ];
        let pileups = pileup(&records, None).unwrap();
        let stats = depth_stats(&pileups[0]);

        // Range: 0..5, length = 6, but positions 2,3 have 0 depth.
        assert_eq!(stats.length, 6);
        assert_eq!(stats.covered, 4);
        assert_eq!(stats.min_depth, 0);
        assert_eq!(stats.max_depth, 1);
        assert!((stats.breadth - 4.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_depth_stats_empty() {
        let p = Pileup {
            rname: "chr1".to_string(),
            columns: vec![],
        };
        let stats = depth_stats(&p);
        assert_eq!(stats.length, 0);
        assert_eq!(stats.covered, 0);
        assert_eq!(stats.min_depth, 0);
        assert_eq!(stats.max_depth, 0);
        assert_eq!(stats.mean_depth, 0.0);
    }

    // -- SNP caller tests --

    #[cfg(feature = "vcf")]
    mod snp_tests {
        use super::*;

        #[test]
        fn test_call_snps_basic() {
            // 3 reads with A at pos 0, 2 reads with T at pos 0 → if ref is A, alt is T
            let mut reference = HashMap::new();
            reference.insert("chr1".to_string(), b"AAAA".to_vec());

            let records = vec![
                make_record("r1", 0, "chr1", 1, 60, "1M", "T", "I"),
                make_record("r2", 0, "chr1", 1, 60, "1M", "T", "I"),
                make_record("r3", 0, "chr1", 1, 60, "1M", "T", "I"),
                make_record("r4", 0, "chr1", 1, 60, "1M", "A", "I"),
            ];
            let pileups = pileup(&records, Some(&reference)).unwrap();
            let variants = call_snps(&pileups[0], 1, 0.2, 1);

            assert_eq!(variants.len(), 1);
            assert_eq!(variants[0].chrom, "chr1");
            assert_eq!(variants[0].position, 1); // 1-based
            assert_eq!(variants[0].ref_allele, vec![b'A']);
            assert_eq!(variants[0].alt_alleles, vec![vec![b'T']]);
        }

        #[test]
        fn test_call_snps_below_threshold() {
            let mut reference = HashMap::new();
            reference.insert("chr1".to_string(), b"AAAA".to_vec());

            // 9 A's and 1 T → alt_freq = 0.1, below 0.2 threshold
            let mut records: Vec<SamRecord> = (0..9)
                .map(|i| make_record(&format!("r{i}"), 0, "chr1", 1, 60, "1M", "A", "I"))
                .collect();
            records.push(make_record("r9", 0, "chr1", 1, 60, "1M", "T", "I"));

            let pileups = pileup(&records, Some(&reference)).unwrap();
            let variants = call_snps(&pileups[0], 1, 0.2, 1);
            assert!(variants.is_empty());
        }

        #[test]
        fn test_call_snps_min_depth() {
            let mut reference = HashMap::new();
            reference.insert("chr1".to_string(), b"AAAA".to_vec());

            let records = vec![
                make_record("r1", 0, "chr1", 1, 60, "1M", "T", "I"),
                make_record("r2", 0, "chr1", 1, 60, "1M", "T", "I"),
            ];
            let pileups = pileup(&records, Some(&reference)).unwrap();

            // min_depth = 5, only 2 reads → no calls
            let variants = call_snps(&pileups[0], 5, 0.1, 1);
            assert!(variants.is_empty());
        }

        #[test]
        fn test_call_snps_variant_quality() {
            let mut reference = HashMap::new();
            reference.insert("chr1".to_string(), b"AAAA".to_vec());

            // All reads are alt → quality should be 99 (capped)
            let records = vec![
                make_record("r1", 0, "chr1", 1, 60, "1M", "G", "I"),
                make_record("r2", 0, "chr1", 1, 60, "1M", "G", "I"),
            ];
            let pileups = pileup(&records, Some(&reference)).unwrap();
            let variants = call_snps(&pileups[0], 1, 0.1, 1);

            assert_eq!(variants.len(), 1);
            assert!(variants[0].quality.unwrap() >= 90.0);
        }
    }

    // -- Edge case tests --

    #[test]
    fn test_pileup_empty_records() {
        let records: Vec<SamRecord> = vec![];
        let pileups = pileup(&records, None).unwrap();
        assert!(pileups.is_empty());
    }

    #[test]
    fn test_pileup_quality_star() {
        // Quality "*" should be treated as 0 for all bases.
        let records = vec![make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "*")];
        let pileups = pileup(&records, None).unwrap();

        for col in &pileups[0].columns {
            for &q in &col.qualities {
                assert_eq!(q, 0);
            }
        }
    }

    #[test]
    fn test_pileup_reverse_strand_bases() {
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "ACGT", "IIII"),  // forward
            make_record("r2", 16, "chr1", 1, 60, "4M", "ACGT", "IIII"), // reverse
        ];
        let pileups = pileup(&records, None).unwrap();

        // Forward strand: uppercase; reverse strand: lowercase.
        assert_eq!(pileups[0].columns[0].bases, vec![b'A', b'a']);
        assert_eq!(pileups[0].columns[1].bases, vec![b'C', b'c']);
        assert_eq!(pileups[0].columns[2].bases, vec![b'G', b'g']);
        assert_eq!(pileups[0].columns[3].bases, vec![b'T', b't']);
    }

    #[test]
    fn test_pileup_hard_clip() {
        // 5H3M5H — hard clips consume neither ref nor query
        let records = vec![make_record("r1", 0, "chr1", 1, 60, "5H3M5H", "ACG", "III")];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups[0].columns.len(), 3);
        assert_eq!(pileups[0].columns[0].bases, vec![b'A']);
    }

    #[test]
    fn test_pileup_n_splice() {
        // N operation (splice) consumes ref only, same as D for pileup purposes
        let records = vec![make_record(
            "r1", 0, "chr1", 1, 60, "2M3N2M", "ACGT", "IIII",
        )];
        let pileups = pileup(&records, None).unwrap();
        // 2 aligned + 3 deleted + 2 aligned = 7 ref positions
        assert_eq!(pileups[0].columns.len(), 7);
        // Positions 2,3,4 should be deletion markers
        assert_eq!(pileups[0].columns[2].bases, vec![b'*']);
    }

    #[test]
    fn test_cigar_hard_soft_clip_combo() {
        // 3H2S4M — hard clip consumes nothing, soft clip consumes 2 query bases,
        // then 4M starting at query position 2
        let records = vec![make_record(
            "r1", 0, "chr1", 1, 60, "3H2S4M", "NNACGT", "IIIIII",
        )];
        let pileups = pileup(&records, None).unwrap();
        assert_eq!(pileups[0].columns.len(), 4);
        // After 2S, query_pos=2, so bases from seq are A, C, G, T
        assert_eq!(pileups[0].columns[0].bases, vec![b'A']);
        assert_eq!(pileups[0].columns[1].bases, vec![b'C']);
        assert_eq!(pileups[0].columns[2].bases, vec![b'G']);
        assert_eq!(pileups[0].columns[3].bases, vec![b'T']);
    }

    #[test]
    fn test_depth_stats_median() {
        // 3 reads at pos 1-4, 1 read at pos 5-6 → depths: [3, 3, 3, 3, 1, 1]
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "4M", "AAAA", "IIII"),
            make_record("r2", 0, "chr1", 1, 60, "4M", "AAAA", "IIII"),
            make_record("r3", 0, "chr1", 1, 60, "6M", "AAAAAA", "IIIIII"),
        ];
        let pileups = pileup(&records, None).unwrap();
        let stats = depth_stats(&pileups[0]);

        assert_eq!(stats.min_depth, 1);
        assert_eq!(stats.max_depth, 3);
        // sorted depths: [1, 1, 3, 3, 3, 3] → median = (3+3)/2 = 3.0
        assert!((stats.median_depth - 3.0).abs() < f64::EPSILON);
    }
}
