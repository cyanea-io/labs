//! BAM/SAM operations: sorting, merging, duplicate marking, and statistics.
//!
//! Provides samtools-equivalent operations on [`SamRecord`] collections.

use crate::bam::BamReference;
use crate::sam::SamRecord;

// ---------------------------------------------------------------------------
// FLAG constants (duplicated from sam.rs for local use)
// ---------------------------------------------------------------------------
const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_FIRST_IN_PAIR: u16 = 0x40;
const FLAG_SECOND_IN_PAIR: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

// ===========================================================================
// Sorting
// ===========================================================================

/// Sort order for BAM records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    /// Sort by reference sequence then position (like `samtools sort`).
    Coordinate,
    /// Sort by query name (like `samtools sort -n`).
    Queryname,
    /// No sorting.
    Unsorted,
}

/// Sort records by coordinate: reference sequence order then position.
///
/// `ref_order` gives the reference name ordering from the BAM header.
/// Records mapped to references not in `ref_order` sort after all known references.
/// Unmapped records (rname == "*") sort last.
pub fn coordinate_sort(records: &mut [SamRecord], ref_order: &[String]) {
    use std::collections::HashMap;
    let ref_idx: HashMap<&str, usize> = ref_order
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let fallback = ref_order.len();
    let unmapped_idx = fallback + 1;

    records.sort_by(|a, b| {
        let a_idx = if a.rname == "*" {
            unmapped_idx
        } else {
            ref_idx.get(a.rname.as_str()).copied().unwrap_or(fallback)
        };
        let b_idx = if b.rname == "*" {
            unmapped_idx
        } else {
            ref_idx.get(b.rname.as_str()).copied().unwrap_or(fallback)
        };
        a_idx.cmp(&b_idx).then(a.pos.cmp(&b.pos))
    });
}

/// Sort records by query name.
pub fn queryname_sort(records: &mut [SamRecord]) {
    records.sort_by(|a, b| a.qname.cmp(&b.qname));
}

// ===========================================================================
// Merge
// ===========================================================================

/// Merge multiple coordinate-sorted record sets into one sorted stream.
///
/// Assumes all inputs are already coordinate-sorted with the same reference order.
pub fn merge_bam_records(inputs: &[&[SamRecord]], ref_order: &[String]) -> Vec<SamRecord> {
    let mut all: Vec<SamRecord> = inputs.iter().flat_map(|s| s.iter().cloned()).collect();
    coordinate_sort(&mut all, ref_order);
    all
}

// ===========================================================================
// Duplicate marking
// ===========================================================================

/// Statistics from duplicate marking.
#[derive(Debug, Clone, Default)]
pub struct DuplicateStats {
    /// Total reads processed.
    pub total_reads: usize,
    /// Number of duplicates marked (FLAG 0x400 set).
    pub duplicates_marked: usize,
}

/// Mark duplicate reads based on mapping position and orientation.
///
/// For each group of reads with the same (rname, pos, strand), the read with
/// the highest sum of base qualities is kept; all others get FLAG 0x400 set.
/// Secondary and supplementary alignments are skipped.
pub fn mark_duplicates(records: &mut [SamRecord]) -> DuplicateStats {
    use std::collections::HashMap;

    let mut stats = DuplicateStats {
        total_reads: records.len(),
        ..Default::default()
    };

    // Group by (rname, pos, is_reverse, is_first_in_pair_or_not_paired)
    // Key: (rname, pos, strand, read_number)
    let mut groups: HashMap<(String, u64, bool, u8), Vec<usize>> = HashMap::new();

    for (i, rec) in records.iter().enumerate() {
        if rec.flag & FLAG_UNMAPPED != 0
            || rec.flag & FLAG_SECONDARY != 0
            || rec.flag & FLAG_SUPPLEMENTARY != 0
        {
            continue;
        }

        let strand = rec.flag & FLAG_REVERSE != 0;
        let read_num = if rec.flag & FLAG_PAIRED != 0 {
            if rec.flag & FLAG_FIRST_IN_PAIR != 0 {
                1
            } else {
                2
            }
        } else {
            0
        };

        let key = (rec.rname.clone(), rec.pos, strand, read_num);
        groups.entry(key).or_default().push(i);
    }

    for indices in groups.values() {
        if indices.len() <= 1 {
            continue;
        }

        // Find the record with the highest quality sum
        let best = indices
            .iter()
            .max_by_key(|&&idx| quality_sum(&records[idx]))
            .copied()
            .unwrap();

        for &idx in indices {
            if idx != best {
                records[idx].flag |= FLAG_DUPLICATE;
                stats.duplicates_marked += 1;
            }
        }
    }

    stats
}

fn quality_sum(record: &SamRecord) -> u64 {
    if record.quality == "*" {
        0
    } else {
        record
            .quality
            .bytes()
            .map(|b| (b.saturating_sub(33)) as u64)
            .sum()
    }
}

// ===========================================================================
// Fixmate
// ===========================================================================

/// Fix mate information for name-sorted paired reads.
///
/// For each pair of reads with the same query name, sets:
/// - RNEXT to mate's RNAME
/// - PNEXT to mate's POS
/// - TLEN (template length)
///
/// Records must be name-sorted. Unpaired reads are skipped.
pub fn fixmate(records: &mut [SamRecord]) {
    let mut i = 0;
    while i + 1 < records.len() {
        if records[i].qname == records[i + 1].qname
            && records[i].flag & FLAG_PAIRED != 0
            && records[i + 1].flag & FLAG_PAIRED != 0
        {
            let rname_a = records[i].rname.clone();
            let rname_b = records[i + 1].rname.clone();
            let pos_a = records[i].pos;
            let pos_b = records[i + 1].pos;

            // Set RNEXT
            records[i].rnext = if rname_b == rname_a {
                "=".to_string()
            } else {
                rname_b
            };
            records[i + 1].rnext = if rname_a == records[i + 1].rname {
                "=".to_string()
            } else {
                rname_a
            };

            // Set PNEXT
            records[i].pnext = pos_b;
            records[i + 1].pnext = pos_a;

            // Set TLEN
            if records[i].rname == records[i + 1].rname
                || (records[i].rnext == "=" || records[i + 1].rnext == "=")
            {
                let (lo, hi) = if pos_a <= pos_b {
                    (pos_a, pos_b)
                } else {
                    (pos_b, pos_a)
                };
                // Approximate: hi - lo + read_length (use seq len of downstream read)
                let tlen = (hi as i64) - (lo as i64) + (records[i].seq_len() as i64);
                records[i].tlen = tlen;
                records[i + 1].tlen = -tlen;
            }

            i += 2;
        } else {
            i += 1;
        }
    }
}

// ===========================================================================
// Idxstats
// ===========================================================================

/// Per-reference sequence alignment counts (like `samtools idxstats`).
#[derive(Debug, Clone)]
pub struct IdxStats {
    /// Reference sequence name.
    pub reference_name: String,
    /// Reference sequence length.
    pub reference_length: u32,
    /// Number of mapped reads.
    pub mapped_count: u64,
    /// Number of unmapped reads placed on this reference.
    pub unmapped_count: u64,
}

/// Compute per-reference mapped/unmapped counts.
pub fn idxstats(records: &[SamRecord], references: &[BamReference]) -> Vec<IdxStats> {
    use std::collections::HashMap;

    let mut mapped: HashMap<&str, u64> = HashMap::new();
    let mut unmapped: HashMap<&str, u64> = HashMap::new();

    for rec in records {
        if rec.rname == "*" {
            continue;
        }
        if rec.flag & FLAG_UNMAPPED != 0 {
            *unmapped.entry(rec.rname.as_str()).or_default() += 1;
        } else {
            *mapped.entry(rec.rname.as_str()).or_default() += 1;
        }
    }

    references
        .iter()
        .map(|r| IdxStats {
            reference_name: r.name.clone(),
            reference_length: r.length,
            mapped_count: mapped.get(r.name.as_str()).copied().unwrap_or(0),
            unmapped_count: unmapped.get(r.name.as_str()).copied().unwrap_or(0),
        })
        .collect()
}

// ===========================================================================
// Flagstat
// ===========================================================================

/// Comprehensive flag-based statistics (like `samtools flagstat`).
#[derive(Debug, Clone, Default)]
pub struct FlagStats {
    /// Total number of reads (QC-passed).
    pub total: usize,
    /// Primary alignments (not secondary, not supplementary).
    pub primary: usize,
    /// Secondary alignments (FLAG 0x100).
    pub secondary: usize,
    /// Supplementary alignments (FLAG 0x800).
    pub supplementary: usize,
    /// Duplicate reads (FLAG 0x400).
    pub duplicates: usize,
    /// Primary alignments that are mapped.
    pub primary_mapped: usize,
    /// Paired reads (FLAG 0x1).
    pub paired: usize,
    /// First in pair (FLAG 0x40).
    pub read1: usize,
    /// Second in pair (FLAG 0x80).
    pub read2: usize,
    /// Properly paired (FLAG 0x2 and both mapped).
    pub properly_paired: usize,
    /// Both reads in pair mapped.
    pub both_mapped: usize,
    /// Singletons: read mapped but mate unmapped.
    pub singletons: usize,
    /// Mate mapped to different chromosome.
    pub mate_different_chrom: usize,
    /// Mate mapped to different chromosome with MAPQ >= 5.
    pub mate_different_chrom_mapq5: usize,
}

/// Compute flag-based statistics.
pub fn flagstat(records: &[SamRecord]) -> FlagStats {
    let mut fs = FlagStats::default();

    for rec in records {
        fs.total += 1;

        let is_secondary = rec.flag & FLAG_SECONDARY != 0;
        let is_supplementary = rec.flag & FLAG_SUPPLEMENTARY != 0;
        let is_unmapped = rec.flag & FLAG_UNMAPPED != 0;
        let is_paired = rec.flag & FLAG_PAIRED != 0;
        let is_mate_unmapped = rec.flag & FLAG_MATE_UNMAPPED != 0;

        if is_secondary {
            fs.secondary += 1;
        } else if is_supplementary {
            fs.supplementary += 1;
        } else {
            fs.primary += 1;
            if !is_unmapped {
                fs.primary_mapped += 1;
            }
        }

        if rec.flag & FLAG_DUPLICATE != 0 {
            fs.duplicates += 1;
        }

        if is_paired {
            fs.paired += 1;

            if rec.flag & FLAG_FIRST_IN_PAIR != 0 {
                fs.read1 += 1;
            }
            if rec.flag & FLAG_SECOND_IN_PAIR != 0 {
                fs.read2 += 1;
            }
            if rec.flag & FLAG_PROPER_PAIR != 0 {
                fs.properly_paired += 1;
            }
            if !is_unmapped && !is_mate_unmapped {
                fs.both_mapped += 1;

                // Check if mate is on different chromosome
                if rec.rnext != "=" && rec.rnext != rec.rname && rec.rnext != "*" {
                    fs.mate_different_chrom += 1;
                    if rec.mapq >= 5 {
                        fs.mate_different_chrom_mapq5 += 1;
                    }
                }
            }
            if !is_unmapped && is_mate_unmapped {
                fs.singletons += 1;
            }
        }
    }

    fs
}

// ===========================================================================
// Alignment statistics
// ===========================================================================

/// Extended alignment statistics.
#[derive(Debug, Clone)]
pub struct AlignmentStats {
    /// Flag-based statistics.
    pub flagstats: FlagStats,
    /// Median insert size for proper pairs.
    pub insert_size_median: f64,
    /// Mean insert size for proper pairs.
    pub insert_size_mean: f64,
    /// Standard deviation of insert sizes for proper pairs.
    pub insert_size_std: f64,
    /// Mean GC content across mapped reads.
    pub gc_content_mean: f64,
}

/// Compute extended alignment statistics.
pub fn alignment_stats(records: &[SamRecord]) -> AlignmentStats {
    let fs = flagstat(records);

    let mut insert_sizes: Vec<f64> = Vec::new();
    let mut gc_values: Vec<f64> = Vec::new();

    for rec in records {
        // Insert sizes from proper pairs (first in pair only, positive TLEN)
        if rec.flag & FLAG_PROPER_PAIR != 0
            && rec.flag & FLAG_FIRST_IN_PAIR != 0
            && rec.tlen > 0
        {
            insert_sizes.push(rec.tlen as f64);
        }

        // GC content of mapped reads
        if rec.flag & FLAG_UNMAPPED == 0 && rec.sequence != "*" {
            let gc = rec
                .sequence
                .bytes()
                .filter(|&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                .count();
            gc_values.push(gc as f64 / rec.sequence.len() as f64);
        }
    }

    let insert_size_mean = if insert_sizes.is_empty() {
        0.0
    } else {
        insert_sizes.iter().sum::<f64>() / insert_sizes.len() as f64
    };

    let insert_size_std = if insert_sizes.len() < 2 {
        0.0
    } else {
        let mean = insert_size_mean;
        let variance =
            insert_sizes.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
                / (insert_sizes.len() - 1) as f64;
        variance.sqrt()
    };

    let insert_size_median = if insert_sizes.is_empty() {
        0.0
    } else {
        insert_sizes.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = insert_sizes.len() / 2;
        if insert_sizes.len() % 2 == 0 {
            (insert_sizes[mid - 1] + insert_sizes[mid]) / 2.0
        } else {
            insert_sizes[mid]
        }
    };

    let gc_content_mean = if gc_values.is_empty() {
        0.0
    } else {
        gc_values.iter().sum::<f64>() / gc_values.len() as f64
    };

    AlignmentStats {
        flagstats: fs,
        insert_size_median,
        insert_size_mean,
        insert_size_std,
        gc_content_mean,
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u64,
        mapq: u8,
        seq: &str,
    ) -> SamRecord {
        SamRecord {
            qname: qname.to_string(),
            flag,
            rname: rname.to_string(),
            pos,
            mapq,
            cigar: format!("{}M", seq.len()),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            sequence: seq.to_string(),
            quality: "*".to_string(),
        }
    }

    fn make_paired_record(
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u64,
        mapq: u8,
        seq: &str,
        rnext: &str,
        pnext: u64,
        tlen: i64,
    ) -> SamRecord {
        SamRecord {
            qname: qname.to_string(),
            flag,
            rname: rname.to_string(),
            pos,
            mapq,
            cigar: format!("{}M", seq.len()),
            rnext: rnext.to_string(),
            pnext,
            tlen,
            sequence: seq.to_string(),
            quality: "*".to_string(),
        }
    }

    // -----------------------------------------------------------------------
    // Sorting
    // -----------------------------------------------------------------------

    #[test]
    fn coordinate_sort_basic() {
        let ref_order = vec!["chr1".to_string(), "chr2".to_string()];
        let mut records = vec![
            make_record("r3", 0, "chr2", 100, 60, "ACGT"),
            make_record("r1", 0, "chr1", 200, 60, "ACGT"),
            make_record("r2", 0, "chr1", 100, 60, "ACGT"),
        ];

        coordinate_sort(&mut records, &ref_order);

        assert_eq!(records[0].rname, "chr1");
        assert_eq!(records[0].pos, 100);
        assert_eq!(records[1].rname, "chr1");
        assert_eq!(records[1].pos, 200);
        assert_eq!(records[2].rname, "chr2");
        assert_eq!(records[2].pos, 100);
    }

    #[test]
    fn coordinate_sort_unmapped_last() {
        let ref_order = vec!["chr1".to_string()];
        let mut records = vec![
            make_record("r1", FLAG_UNMAPPED, "*", 0, 0, "*"),
            make_record("r2", 0, "chr1", 100, 60, "ACGT"),
        ];

        coordinate_sort(&mut records, &ref_order);

        assert_eq!(records[0].rname, "chr1");
        assert_eq!(records[1].rname, "*");
    }

    #[test]
    fn queryname_sort_basic() {
        let mut records = vec![
            make_record("readC", 0, "chr1", 100, 60, "ACGT"),
            make_record("readA", 0, "chr1", 200, 60, "ACGT"),
            make_record("readB", 0, "chr1", 150, 60, "ACGT"),
        ];

        queryname_sort(&mut records);

        assert_eq!(records[0].qname, "readA");
        assert_eq!(records[1].qname, "readB");
        assert_eq!(records[2].qname, "readC");
    }

    // -----------------------------------------------------------------------
    // Merge
    // -----------------------------------------------------------------------

    #[test]
    fn merge_basic() {
        let ref_order = vec!["chr1".to_string()];
        let set1 = vec![
            make_record("r1", 0, "chr1", 100, 60, "ACGT"),
            make_record("r3", 0, "chr1", 300, 60, "ACGT"),
        ];
        let set2 = vec![make_record("r2", 0, "chr1", 200, 60, "ACGT")];

        let merged = merge_bam_records(&[&set1, &set2], &ref_order);

        assert_eq!(merged.len(), 3);
        assert_eq!(merged[0].pos, 100);
        assert_eq!(merged[1].pos, 200);
        assert_eq!(merged[2].pos, 300);
    }

    // -----------------------------------------------------------------------
    // Duplicate marking
    // -----------------------------------------------------------------------

    #[test]
    fn mark_duplicates_basic() {
        let mut records = vec![
            SamRecord {
                quality: "IIII".to_string(), // higher quality (sum = 4 * 40)
                ..make_record("r1", 0, "chr1", 100, 60, "ACGT")
            },
            SamRecord {
                quality: "####".to_string(), // lower quality (sum = 4 * 2)
                ..make_record("r2", 0, "chr1", 100, 60, "ACGT")
            },
        ];

        let stats = mark_duplicates(&mut records);

        assert_eq!(stats.total_reads, 2);
        assert_eq!(stats.duplicates_marked, 1);
        // Higher quality read should NOT be marked as duplicate
        assert_eq!(records[0].flag & FLAG_DUPLICATE, 0);
        // Lower quality read should be marked
        assert_ne!(records[1].flag & FLAG_DUPLICATE, 0);
    }

    #[test]
    fn mark_duplicates_different_positions() {
        let mut records = vec![
            make_record("r1", 0, "chr1", 100, 60, "ACGT"),
            make_record("r2", 0, "chr1", 200, 60, "ACGT"),
        ];

        let stats = mark_duplicates(&mut records);

        assert_eq!(stats.duplicates_marked, 0);
    }

    #[test]
    fn mark_duplicates_skips_unmapped() {
        let mut records = vec![
            make_record("r1", FLAG_UNMAPPED, "*", 0, 0, "ACGT"),
            make_record("r2", FLAG_UNMAPPED, "*", 0, 0, "ACGT"),
        ];

        let stats = mark_duplicates(&mut records);

        assert_eq!(stats.duplicates_marked, 0);
    }

    #[test]
    fn mark_duplicates_skips_secondary() {
        let mut records = vec![
            make_record("r1", FLAG_SECONDARY, "chr1", 100, 60, "ACGT"),
            make_record("r2", 0, "chr1", 100, 60, "ACGT"),
        ];

        let stats = mark_duplicates(&mut records);

        assert_eq!(stats.duplicates_marked, 0);
    }

    // -----------------------------------------------------------------------
    // Fixmate
    // -----------------------------------------------------------------------

    #[test]
    fn fixmate_basic() {
        let mut records = vec![
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_FIRST_IN_PAIR,
                "chr1",
                100,
                60,
                "ACGTACGT",
                "*",
                0,
                0,
            ),
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_SECOND_IN_PAIR,
                "chr1",
                200,
                60,
                "ACGTACGT",
                "*",
                0,
                0,
            ),
        ];

        fixmate(&mut records);

        assert_eq!(records[0].rnext, "=");
        assert_eq!(records[0].pnext, 200);
        assert!(records[0].tlen > 0);
        assert_eq!(records[1].rnext, "=");
        assert_eq!(records[1].pnext, 100);
        assert!(records[1].tlen < 0);
    }

    #[test]
    fn fixmate_different_chroms() {
        let mut records = vec![
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_FIRST_IN_PAIR,
                "chr1",
                100,
                60,
                "ACGT",
                "*",
                0,
                0,
            ),
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_SECOND_IN_PAIR,
                "chr2",
                200,
                60,
                "ACGT",
                "*",
                0,
                0,
            ),
        ];

        fixmate(&mut records);

        assert_eq!(records[0].rnext, "chr2");
        assert_eq!(records[0].pnext, 200);
        assert_eq!(records[1].rnext, "chr1");
        assert_eq!(records[1].pnext, 100);
    }

    // -----------------------------------------------------------------------
    // Idxstats
    // -----------------------------------------------------------------------

    #[test]
    fn idxstats_basic() {
        let records = vec![
            make_record("r1", 0, "chr1", 100, 60, "ACGT"),
            make_record("r2", 0, "chr1", 200, 60, "ACGT"),
            make_record("r3", 0, "chr2", 100, 60, "ACGT"),
            make_record("r4", FLAG_UNMAPPED, "*", 0, 0, "ACGT"),
        ];
        let refs = vec![
            BamReference {
                name: "chr1".to_string(),
                length: 1000,
            },
            BamReference {
                name: "chr2".to_string(),
                length: 2000,
            },
        ];

        let stats = idxstats(&records, &refs);

        assert_eq!(stats.len(), 2);
        assert_eq!(stats[0].reference_name, "chr1");
        assert_eq!(stats[0].mapped_count, 2);
        assert_eq!(stats[0].unmapped_count, 0);
        assert_eq!(stats[1].reference_name, "chr2");
        assert_eq!(stats[1].mapped_count, 1);
        assert_eq!(stats[1].unmapped_count, 0);
    }

    // -----------------------------------------------------------------------
    // Flagstat
    // -----------------------------------------------------------------------

    #[test]
    fn flagstat_basic() {
        let records = vec![
            make_record("r1", 0, "chr1", 100, 60, "ACGT"),
            make_record("r2", FLAG_UNMAPPED, "*", 0, 0, "ACGT"),
            make_record("r3", FLAG_SECONDARY, "chr1", 100, 30, "ACGT"),
            make_record("r4", FLAG_SUPPLEMENTARY, "chr1", 100, 30, "ACGT"),
            SamRecord {
                flag: FLAG_DUPLICATE,
                ..make_record("r5", 0, "chr1", 100, 60, "ACGT")
            },
        ];

        let fs = flagstat(&records);

        assert_eq!(fs.total, 5);
        assert_eq!(fs.primary, 3); // r1, r2, r5
        assert_eq!(fs.secondary, 1);
        assert_eq!(fs.supplementary, 1);
        assert_eq!(fs.primary_mapped, 2); // r1, r5
        assert_eq!(fs.duplicates, 1);
    }

    #[test]
    fn flagstat_paired() {
        let records = vec![
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_IN_PAIR,
                "chr1",
                100,
                60,
                "ACGT",
                "=",
                200,
                100,
            ),
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_SECOND_IN_PAIR,
                "chr1",
                200,
                60,
                "ACGT",
                "=",
                100,
                -100,
            ),
        ];

        let fs = flagstat(&records);

        assert_eq!(fs.paired, 2);
        assert_eq!(fs.read1, 1);
        assert_eq!(fs.read2, 1);
        assert_eq!(fs.properly_paired, 2);
        assert_eq!(fs.both_mapped, 2);
        assert_eq!(fs.singletons, 0);
        assert_eq!(fs.mate_different_chrom, 0);
    }

    #[test]
    fn flagstat_singleton() {
        let records = vec![make_paired_record(
            "pair1",
            FLAG_PAIRED | FLAG_MATE_UNMAPPED | FLAG_FIRST_IN_PAIR,
            "chr1",
            100,
            60,
            "ACGT",
            "*",
            0,
            0,
        )];

        let fs = flagstat(&records);

        assert_eq!(fs.singletons, 1);
        assert_eq!(fs.both_mapped, 0);
    }

    #[test]
    fn flagstat_mate_different_chrom() {
        let records = vec![make_paired_record(
            "pair1",
            FLAG_PAIRED | FLAG_FIRST_IN_PAIR,
            "chr1",
            100,
            60,
            "ACGT",
            "chr2",
            200,
            0,
        )];

        let fs = flagstat(&records);

        assert_eq!(fs.mate_different_chrom, 1);
        assert_eq!(fs.mate_different_chrom_mapq5, 1);
    }

    #[test]
    fn flagstat_mate_different_chrom_low_mapq() {
        let records = vec![make_paired_record(
            "pair1",
            FLAG_PAIRED | FLAG_FIRST_IN_PAIR,
            "chr1",
            100,
            3,
            "ACGT",
            "chr2",
            200,
            0,
        )];

        let fs = flagstat(&records);

        assert_eq!(fs.mate_different_chrom, 1);
        assert_eq!(fs.mate_different_chrom_mapq5, 0);
    }

    // -----------------------------------------------------------------------
    // Alignment stats
    // -----------------------------------------------------------------------

    #[test]
    fn alignment_stats_insert_size() {
        let records = vec![
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_IN_PAIR,
                "chr1",
                100,
                60,
                "ACGT",
                "=",
                200,
                150,
            ),
            make_paired_record(
                "pair1",
                FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_SECOND_IN_PAIR,
                "chr1",
                200,
                60,
                "ACGT",
                "=",
                100,
                -150,
            ),
            make_paired_record(
                "pair2",
                FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_IN_PAIR,
                "chr1",
                300,
                60,
                "ACGT",
                "=",
                500,
                250,
            ),
            make_paired_record(
                "pair2",
                FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_SECOND_IN_PAIR,
                "chr1",
                500,
                60,
                "ACGT",
                "=",
                300,
                -250,
            ),
        ];

        let stats = alignment_stats(&records);

        // Insert sizes from first-in-pair with positive TLEN: 150, 250
        assert!((stats.insert_size_mean - 200.0).abs() < 0.01);
        assert!((stats.insert_size_median - 200.0).abs() < 0.01);
        assert!(stats.insert_size_std > 0.0);
    }

    #[test]
    fn alignment_stats_gc_content() {
        let records = vec![
            make_record("r1", 0, "chr1", 100, 60, "GGCC"), // 100% GC
            make_record("r2", 0, "chr1", 200, 60, "AATT"), // 0% GC
        ];

        let stats = alignment_stats(&records);

        assert!((stats.gc_content_mean - 0.5).abs() < 0.01);
    }

    #[test]
    fn alignment_stats_empty() {
        let stats = alignment_stats(&[]);

        assert_eq!(stats.flagstats.total, 0);
        assert_eq!(stats.insert_size_mean, 0.0);
        assert_eq!(stats.gc_content_mean, 0.0);
    }

    // -----------------------------------------------------------------------
    // Merge multi-chrom
    // -----------------------------------------------------------------------

    #[test]
    fn merge_multi_chrom() {
        let ref_order = vec!["chr1".to_string(), "chr2".to_string()];
        let set1 = vec![make_record("r1", 0, "chr2", 100, 60, "ACGT")];
        let set2 = vec![make_record("r2", 0, "chr1", 100, 60, "ACGT")];

        let merged = merge_bam_records(&[&set1, &set2], &ref_order);

        assert_eq!(merged[0].rname, "chr1");
        assert_eq!(merged[1].rname, "chr2");
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn mark_duplicates_three_at_same_pos() {
        let mut records = vec![
            SamRecord {
                quality: "AAAA".to_string(), // sum = 4 * 32
                ..make_record("r1", 0, "chr1", 100, 60, "ACGT")
            },
            SamRecord {
                quality: "IIII".to_string(), // sum = 4 * 40 (highest)
                ..make_record("r2", 0, "chr1", 100, 60, "ACGT")
            },
            SamRecord {
                quality: "!!!!".to_string(), // sum = 4 * 0 (lowest)
                ..make_record("r3", 0, "chr1", 100, 60, "ACGT")
            },
        ];

        let stats = mark_duplicates(&mut records);

        assert_eq!(stats.duplicates_marked, 2);
        // r2 (highest quality) should be kept
        assert_eq!(records[1].flag & FLAG_DUPLICATE, 0);
        assert_ne!(records[0].flag & FLAG_DUPLICATE, 0);
        assert_ne!(records[2].flag & FLAG_DUPLICATE, 0);
    }

    #[test]
    fn coordinate_sort_same_chrom_same_pos() {
        let ref_order = vec!["chr1".to_string()];
        let mut records = vec![
            make_record("r1", 0, "chr1", 100, 60, "ACGT"),
            make_record("r2", 0, "chr1", 100, 60, "ACGT"),
        ];

        coordinate_sort(&mut records, &ref_order);
        // Both have same position, order is stable-ish
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn fixmate_unpaired_skipped() {
        let mut records = vec![
            make_record("r1", 0, "chr1", 100, 60, "ACGT"),
            make_record("r2", 0, "chr1", 200, 60, "ACGT"),
        ];

        fixmate(&mut records);

        // Unpaired reads should not be modified
        assert_eq!(records[0].rnext, "*");
        assert_eq!(records[1].rnext, "*");
    }
}
