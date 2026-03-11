//! Structural variant detection from long reads.
//!
//! Detects insertions, deletions, inversions, duplications, and translocations
//! from split/supplementary alignments and read depth signals in long-read data.

use cyanea_core::Result;

/// Type of structural variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SvType {
    /// Insertion relative to reference.
    Insertion,
    /// Deletion relative to reference.
    Deletion,
    /// Inverted segment.
    Inversion,
    /// Tandem duplication.
    Duplication,
    /// Breakend / translocation.
    Breakend,
}

impl SvType {
    /// VCF-style string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            SvType::Insertion => "INS",
            SvType::Deletion => "DEL",
            SvType::Inversion => "INV",
            SvType::Duplication => "DUP",
            SvType::Breakend => "BND",
        }
    }
}

/// A detected structural variant.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StructuralVariant {
    /// Chromosome / contig name.
    pub chrom: String,
    /// Start position (0-based).
    pub start: u64,
    /// End position (0-based, exclusive).
    pub end: u64,
    /// Variant type.
    pub sv_type: SvType,
    /// Variant length (positive for insertions, negative for deletions).
    pub length: i64,
    /// Supporting read count.
    pub support: usize,
    /// Quality score (Phred-scaled).
    pub quality: f64,
    /// Genotype estimate (e.g., "0/1" or "1/1").
    pub genotype: String,
    /// Inserted sequence (for insertions).
    pub inserted_seq: Option<Vec<u8>>,
}

/// A split alignment segment from a single read.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SplitAlignment {
    /// Read name.
    pub read_name: String,
    /// Reference name for this segment.
    pub chrom: String,
    /// Reference start position.
    pub ref_start: u64,
    /// Reference end position.
    pub ref_end: u64,
    /// Query (read) start position.
    pub query_start: u64,
    /// Query (read) end position.
    pub query_end: u64,
    /// Whether this segment is on the reverse strand.
    pub is_reverse: bool,
    /// Mapping quality.
    pub mapq: u8,
}

/// Configuration for SV calling.
#[derive(Debug, Clone)]
pub struct SvCallConfig {
    /// Minimum SV length to report (default 50).
    pub min_sv_length: u64,
    /// Maximum SV length (default 100_000_000).
    pub max_sv_length: u64,
    /// Minimum number of supporting reads (default 3).
    pub min_support: usize,
    /// Maximum distance to cluster breakpoints (default 500).
    pub cluster_distance: u64,
    /// Minimum mapping quality for split segments (default 20).
    pub min_mapq: u8,
}

impl Default for SvCallConfig {
    fn default() -> Self {
        Self {
            min_sv_length: 50,
            max_sv_length: 100_000_000,
            min_support: 3,
            cluster_distance: 500,
            min_mapq: 20,
        }
    }
}

/// Detect structural variants from split alignments.
///
/// Groups split alignment segments by read name, then identifies
/// SVs from discordant segments:
///
/// - **Deletion**: Gap in reference between consecutive segments on same strand
/// - **Insertion**: Gap in query between consecutive segments
/// - **Inversion**: Consecutive segments on different strands
/// - **Duplication**: Overlapping reference positions
/// - **Breakend**: Segments on different chromosomes
pub fn call_svs(
    alignments: &[SplitAlignment],
    config: &SvCallConfig,
) -> Result<Vec<StructuralVariant>> {
    if alignments.is_empty() {
        return Ok(Vec::new());
    }

    // Group by read name
    let mut by_read: std::collections::HashMap<&str, Vec<&SplitAlignment>> =
        std::collections::HashMap::new();
    for aln in alignments {
        if aln.mapq >= config.min_mapq {
            by_read.entry(&aln.read_name).or_default().push(aln);
        }
    }

    // Collect raw SV signals
    let mut raw_svs: Vec<StructuralVariant> = Vec::new();

    for (_read_name, segments) in &by_read {
        if segments.len() < 2 {
            continue;
        }

        let mut sorted = segments.clone();
        sorted.sort_by_key(|s| s.query_start);

        for pair in sorted.windows(2) {
            let s1 = pair[0];
            let s2 = pair[1];

            // Different chromosomes → breakend
            if s1.chrom != s2.chrom {
                raw_svs.push(StructuralVariant {
                    chrom: s1.chrom.clone(),
                    start: s1.ref_end,
                    end: s1.ref_end + 1,
                    sv_type: SvType::Breakend,
                    length: 0,
                    support: 1,
                    quality: s1.mapq.min(s2.mapq) as f64,
                    genotype: "0/1".into(),
                    inserted_seq: None,
                });
                continue;
            }

            // Different strands → inversion
            if s1.is_reverse != s2.is_reverse {
                let start = s1.ref_start.min(s2.ref_start);
                let end = s1.ref_end.max(s2.ref_end);
                let length = (end - start) as i64;

                if length >= config.min_sv_length as i64 && length <= config.max_sv_length as i64 {
                    raw_svs.push(StructuralVariant {
                        chrom: s1.chrom.clone(),
                        start,
                        end,
                        sv_type: SvType::Inversion,
                        length,
                        support: 1,
                        quality: s1.mapq.min(s2.mapq) as f64,
                        genotype: "0/1".into(),
                        inserted_seq: None,
                    });
                }
                continue;
            }

            // Same strand, same chromosome
            let ref_gap = if s2.ref_start > s1.ref_end {
                (s2.ref_start - s1.ref_end) as i64
            } else {
                -((s1.ref_end - s2.ref_start) as i64)
            };

            let query_gap = if s2.query_start > s1.query_end {
                (s2.query_start - s1.query_end) as i64
            } else {
                0
            };

            let size_diff = ref_gap - query_gap;

            if size_diff.unsigned_abs() < config.min_sv_length {
                continue;
            }

            if size_diff > 0 && (size_diff as u64) <= config.max_sv_length {
                // Deletion: reference gap larger than query gap
                raw_svs.push(StructuralVariant {
                    chrom: s1.chrom.clone(),
                    start: s1.ref_end,
                    end: s2.ref_start,
                    sv_type: SvType::Deletion,
                    length: -size_diff,
                    support: 1,
                    quality: s1.mapq.min(s2.mapq) as f64,
                    genotype: "0/1".into(),
                    inserted_seq: None,
                });
            } else if size_diff < 0 && ((-size_diff) as u64) <= config.max_sv_length {
                // Insertion: query gap larger than reference gap
                raw_svs.push(StructuralVariant {
                    chrom: s1.chrom.clone(),
                    start: s1.ref_end,
                    end: s1.ref_end + 1,
                    sv_type: SvType::Insertion,
                    length: -size_diff,
                    support: 1,
                    quality: s1.mapq.min(s2.mapq) as f64,
                    genotype: "0/1".into(),
                    inserted_seq: None,
                });
            }

            // Overlapping reference → duplication
            if ref_gap < 0 && (-ref_gap as u64) >= config.min_sv_length {
                raw_svs.push(StructuralVariant {
                    chrom: s1.chrom.clone(),
                    start: s2.ref_start,
                    end: s1.ref_end,
                    sv_type: SvType::Duplication,
                    length: -ref_gap,
                    support: 1,
                    quality: s1.mapq.min(s2.mapq) as f64,
                    genotype: "0/1".into(),
                    inserted_seq: None,
                });
            }
        }
    }

    // Cluster nearby SVs of the same type
    let clustered = cluster_svs(&mut raw_svs, config);

    // Filter by support
    let filtered: Vec<StructuralVariant> = clustered.into_iter()
        .filter(|sv| sv.support >= config.min_support)
        .collect();

    Ok(filtered)
}

/// Cluster nearby SVs of the same type, merging support counts.
pub fn cluster_svs(svs: &mut [StructuralVariant], config: &SvCallConfig) -> Vec<StructuralVariant> {
    if svs.is_empty() {
        return Vec::new();
    }

    svs.sort_by(|a, b| {
        a.chrom.cmp(&b.chrom)
            .then(a.sv_type.as_str().cmp(b.sv_type.as_str()))
            .then(a.start.cmp(&b.start))
    });

    let mut clustered = Vec::new();
    let mut current = svs[0].clone();

    for sv in svs.iter().skip(1) {
        if sv.chrom == current.chrom
            && sv.sv_type == current.sv_type
            && sv.start.abs_diff(current.start) <= config.cluster_distance
        {
            current.support += sv.support;
            current.quality = current.quality.max(sv.quality);
            // Keep the breakpoint with more support
            if sv.quality > current.quality {
                current.start = sv.start;
                current.end = sv.end;
                current.length = sv.length;
            }
        } else {
            clustered.push(current);
            current = sv.clone();
        }
    }
    clustered.push(current);

    clustered
}

/// Detect SVs from CIGAR-level indels in aligned reads.
///
/// Extracts large insertions and deletions (>= min_length) directly
/// from CIGAR strings without requiring split alignments.
pub fn svs_from_cigar(
    chrom: &str,
    ref_start: u64,
    cigar_ops: &[(char, u32)],
    min_length: u64,
) -> Vec<StructuralVariant> {
    let mut svs = Vec::new();
    let mut ref_pos = ref_start;

    for &(op, len) in cigar_ops {
        match op {
            'M' | '=' | 'X' => ref_pos += len as u64,
            'D' | 'N' => {
                if len as u64 >= min_length {
                    svs.push(StructuralVariant {
                        chrom: chrom.to_string(),
                        start: ref_pos,
                        end: ref_pos + len as u64,
                        sv_type: SvType::Deletion,
                        length: -(len as i64),
                        support: 1,
                        quality: 0.0,
                        genotype: "0/1".into(),
                        inserted_seq: None,
                    });
                }
                ref_pos += len as u64;
            }
            'I' => {
                if len as u64 >= min_length {
                    svs.push(StructuralVariant {
                        chrom: chrom.to_string(),
                        start: ref_pos,
                        end: ref_pos + 1,
                        sv_type: SvType::Insertion,
                        length: len as i64,
                        support: 1,
                        quality: 0.0,
                        genotype: "0/1".into(),
                        inserted_seq: None,
                    });
                }
            }
            'S' | 'H' => {} // clips don't affect reference position
            _ => {}
        }
    }

    svs
}

/// Summary of SV calls.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SvSummary {
    pub total: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub inversions: usize,
    pub duplications: usize,
    pub breakends: usize,
    pub mean_length: f64,
    pub median_length: f64,
}

/// Summarize a set of SV calls.
pub fn sv_summary(svs: &[StructuralVariant]) -> SvSummary {
    let mut lengths: Vec<i64> = svs.iter().map(|sv| sv.length.abs()).collect();
    lengths.sort_unstable();

    let median = if lengths.is_empty() {
        0.0
    } else if lengths.len() % 2 == 0 {
        (lengths[lengths.len() / 2 - 1] + lengths[lengths.len() / 2]) as f64 / 2.0
    } else {
        lengths[lengths.len() / 2] as f64
    };

    let mean = if lengths.is_empty() {
        0.0
    } else {
        lengths.iter().sum::<i64>() as f64 / lengths.len() as f64
    };

    SvSummary {
        total: svs.len(),
        insertions: svs.iter().filter(|sv| sv.sv_type == SvType::Insertion).count(),
        deletions: svs.iter().filter(|sv| sv.sv_type == SvType::Deletion).count(),
        inversions: svs.iter().filter(|sv| sv.sv_type == SvType::Inversion).count(),
        duplications: svs.iter().filter(|sv| sv.sv_type == SvType::Duplication).count(),
        breakends: svs.iter().filter(|sv| sv.sv_type == SvType::Breakend).count(),
        mean_length: mean,
        median_length: median,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_split(
        read: &str, chrom: &str, ref_start: u64, ref_end: u64,
        query_start: u64, query_end: u64, is_reverse: bool,
    ) -> SplitAlignment {
        SplitAlignment {
            read_name: read.into(),
            chrom: chrom.into(),
            ref_start,
            ref_end,
            query_start,
            query_end,
            is_reverse,
            mapq: 60,
        }
    }

    #[test]
    fn test_sv_type_str() {
        assert_eq!(SvType::Insertion.as_str(), "INS");
        assert_eq!(SvType::Deletion.as_str(), "DEL");
        assert_eq!(SvType::Inversion.as_str(), "INV");
        assert_eq!(SvType::Duplication.as_str(), "DUP");
        assert_eq!(SvType::Breakend.as_str(), "BND");
    }

    #[test]
    fn test_deletion_from_splits() {
        // Read maps in two segments with a large gap in reference
        let alns: Vec<SplitAlignment> = (0..5).flat_map(|i| {
            vec![
                make_split(&format!("read_{}", i), "chr1", 1000, 2000, 0, 1000, false),
                make_split(&format!("read_{}", i), "chr1", 5000, 6000, 1000, 2000, false),
            ]
        }).collect();

        let config = SvCallConfig { min_support: 3, ..Default::default() };
        let svs = call_svs(&alns, &config).unwrap();

        assert!(!svs.is_empty());
        assert!(svs.iter().any(|sv| sv.sv_type == SvType::Deletion));
    }

    #[test]
    fn test_inversion_from_splits() {
        // Read maps forward then reverse on same chromosome
        let alns: Vec<SplitAlignment> = (0..4).flat_map(|i| {
            vec![
                make_split(&format!("read_{}", i), "chr1", 1000, 2000, 0, 1000, false),
                make_split(&format!("read_{}", i), "chr1", 3000, 4000, 1000, 2000, true), // reversed!
            ]
        }).collect();

        let config = SvCallConfig { min_support: 3, ..Default::default() };
        let svs = call_svs(&alns, &config).unwrap();

        assert!(svs.iter().any(|sv| sv.sv_type == SvType::Inversion));
    }

    #[test]
    fn test_breakend_from_splits() {
        // Read maps to different chromosomes
        let alns: Vec<SplitAlignment> = (0..4).flat_map(|i| {
            vec![
                make_split(&format!("read_{}", i), "chr1", 1000, 2000, 0, 1000, false),
                make_split(&format!("read_{}", i), "chr5", 3000, 4000, 1000, 2000, false),
            ]
        }).collect();

        let config = SvCallConfig { min_support: 3, ..Default::default() };
        let svs = call_svs(&alns, &config).unwrap();

        assert!(svs.iter().any(|sv| sv.sv_type == SvType::Breakend));
    }

    #[test]
    fn test_svs_from_cigar() {
        let cigar = vec![
            ('M', 1000),
            ('D', 500),  // 500bp deletion
            ('M', 500),
            ('I', 200),  // 200bp insertion
            ('M', 1000),
        ];

        let svs = svs_from_cigar("chr1", 10000, &cigar, 50);
        assert_eq!(svs.len(), 2);
        assert!(svs.iter().any(|sv| sv.sv_type == SvType::Deletion && sv.length == -500));
        assert!(svs.iter().any(|sv| sv.sv_type == SvType::Insertion && sv.length == 200));
    }

    #[test]
    fn test_svs_from_cigar_below_threshold() {
        let cigar = vec![('M', 100), ('D', 10), ('M', 100)];
        let svs = svs_from_cigar("chr1", 0, &cigar, 50);
        assert!(svs.is_empty()); // 10bp deletion below 50bp threshold
    }

    #[test]
    fn test_sv_summary() {
        let svs = vec![
            StructuralVariant {
                chrom: "chr1".into(), start: 1000, end: 2000, sv_type: SvType::Deletion,
                length: -1000, support: 5, quality: 60.0, genotype: "0/1".into(), inserted_seq: None,
            },
            StructuralVariant {
                chrom: "chr1".into(), start: 5000, end: 5001, sv_type: SvType::Insertion,
                length: 500, support: 3, quality: 40.0, genotype: "0/1".into(), inserted_seq: None,
            },
        ];

        let summary = sv_summary(&svs);
        assert_eq!(summary.total, 2);
        assert_eq!(summary.deletions, 1);
        assert_eq!(summary.insertions, 1);
        assert!((summary.mean_length - 750.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_calls() {
        let config = SvCallConfig::default();
        let svs = call_svs(&[], &config).unwrap();
        assert!(svs.is_empty());
    }

    #[test]
    fn test_low_support_filtered() {
        // Only 2 reads support this SV, but min_support is 3
        let alns: Vec<SplitAlignment> = (0..2).flat_map(|i| {
            vec![
                make_split(&format!("read_{}", i), "chr1", 1000, 2000, 0, 1000, false),
                make_split(&format!("read_{}", i), "chr1", 5000, 6000, 1000, 2000, false),
            ]
        }).collect();

        let config = SvCallConfig { min_support: 3, ..Default::default() };
        let svs = call_svs(&alns, &config).unwrap();
        assert!(svs.is_empty());
    }

    #[test]
    fn test_low_mapq_filtered() {
        let mut alns: Vec<SplitAlignment> = (0..5).flat_map(|i| {
            vec![
                make_split(&format!("read_{}", i), "chr1", 1000, 2000, 0, 1000, false),
                make_split(&format!("read_{}", i), "chr1", 5000, 6000, 1000, 2000, false),
            ]
        }).collect();
        // Set all mapq to 5 (below threshold)
        for aln in &mut alns {
            aln.mapq = 5;
        }

        let config = SvCallConfig { min_support: 3, min_mapq: 20, ..Default::default() };
        let svs = call_svs(&alns, &config).unwrap();
        assert!(svs.is_empty());
    }
}
