//! SAM (Sequence Alignment/Map) parser.
//!
//! Parses SAM text format files into [`SamRecord`] records. SAM is a
//! tab-delimited text format with 11 mandatory fields per alignment line.
//! Header lines (starting with `@`) are skipped.
//!
//! This is a lightweight, pure-Rust parser for the text SAM format.
//! BAM (binary) format is not supported.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

/// Flag bit: read is unmapped.
const FLAG_UNMAPPED: u16 = 0x4;

/// A single SAM alignment record.
///
/// Contains the 11 mandatory fields from the SAM specification.
/// Optional fields (columns 12+) are currently ignored.
#[derive(Debug, Clone)]
pub struct SamRecord {
    /// Query template name (QNAME).
    pub qname: String,
    /// Bitwise flag (FLAG).
    pub flag: u16,
    /// Reference sequence name (RNAME). `"*"` if unmapped.
    pub rname: String,
    /// 1-based leftmost mapping position (POS). 0 if unmapped.
    pub pos: u64,
    /// Mapping quality (MAPQ). 255 if unavailable.
    pub mapq: u8,
    /// CIGAR string. `"*"` if unavailable.
    pub cigar: String,
    /// Query sequence (SEQ). `"*"` if not stored.
    pub sequence: String,
    /// ASCII of base quality plus 33 (QUAL). `"*"` if not stored.
    pub quality: String,
}

impl SamRecord {
    /// Returns `true` if the read is mapped (FLAG bit 0x4 is NOT set).
    pub fn is_mapped(&self) -> bool {
        self.flag & FLAG_UNMAPPED == 0
    }

    /// Returns `true` if the read is unmapped (FLAG bit 0x4 is set).
    pub fn is_unmapped(&self) -> bool {
        self.flag & FLAG_UNMAPPED != 0
    }

    /// Returns the length of the query sequence, or 0 if SEQ is `"*"`.
    pub fn seq_len(&self) -> usize {
        if self.sequence == "*" {
            0
        } else {
            self.sequence.len()
        }
    }
}

/// Summary statistics for a collection of SAM records.
#[derive(Debug, Clone)]
pub struct SamStats {
    /// Total number of alignment records.
    pub total_reads: usize,
    /// Number of mapped reads (FLAG bit 0x4 NOT set).
    pub mapped: usize,
    /// Number of unmapped reads (FLAG bit 0x4 set).
    pub unmapped: usize,
    /// Mean mapping quality of mapped reads.
    pub avg_mapq: f64,
    /// Mean sequence length across all reads.
    pub avg_length: f64,
    /// Distribution of MAPQ values as `(mapq_value, count)` pairs,
    /// sorted by MAPQ value. Only includes values with count > 0.
    pub mapq_distribution: Vec<(u8, usize)>,
}

/// Parse a SAM file from a file path and return all alignment records.
///
/// Header lines (starting with `@`) and empty lines are skipped.
/// Each data line is parsed into a [`SamRecord`].
pub fn parse_sam(path: impl AsRef<Path>) -> Result<Vec<SamRecord>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    parse_sam_reader(reader, path)
}

/// Parse SAM records from a buffered reader.
///
/// Header lines (starting with `@`) and empty lines are skipped.
fn parse_sam_reader(reader: impl BufRead, path: &Path) -> Result<Vec<SamRecord>> {
    let mut data_lines: Vec<(usize, String)> = Vec::new();
    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let trimmed = line.trim().to_string();
        if trimmed.is_empty() || trimmed.starts_with('@') {
            continue;
        }
        data_lines.push((line_num + 1, trimmed));
    }

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data_lines
            .par_iter()
            .map(|(line_num, line)| parse_sam_record(line, *line_num, path))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    data_lines
        .iter()
        .map(|(line_num, line)| parse_sam_record(line, *line_num, path))
        .collect()
}

/// Parse a single SAM alignment line into a [`SamRecord`].
///
/// The SAM format requires exactly 11 mandatory tab-separated fields:
/// QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
///
/// Optional fields (columns 12+) are silently ignored.
fn parse_sam_record(line: &str, line_num: usize, path: &Path) -> Result<SamRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected at least 11 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let qname = fields[0].to_string();

    let flag: u16 = fields[1].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid FLAG '{}'",
            path.display(),
            line_num,
            fields[1]
        ))
    })?;

    let rname = fields[2].to_string();

    let pos: u64 = fields[3].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid POS '{}'",
            path.display(),
            line_num,
            fields[3]
        ))
    })?;

    let mapq: u8 = fields[4].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid MAPQ '{}'",
            path.display(),
            line_num,
            fields[4]
        ))
    })?;

    let cigar = fields[5].to_string();

    // fields[6] = RNEXT, fields[7] = PNEXT, fields[8] = TLEN — skipped for now

    let sequence = fields[9].to_string();
    let quality = fields[10].to_string();

    Ok(SamRecord {
        qname,
        flag,
        rname,
        pos,
        mapq,
        cigar,
        sequence,
        quality,
    })
}

/// Compute summary statistics from a slice of SAM records.
///
/// - `mapped` = reads where FLAG bit 0x4 is NOT set
/// - `unmapped` = reads where FLAG bit 0x4 is set
/// - `avg_mapq` = mean MAPQ of mapped reads only (0.0 if no mapped reads)
/// - `avg_length` = mean sequence length across all reads
/// - `mapq_distribution` = histogram of MAPQ values (mapped reads only)
pub fn sam_stats(records: &[SamRecord]) -> SamStats {
    let total_reads = records.len();
    let mut mapped: usize = 0;
    let mut unmapped: usize = 0;
    let mut mapq_sum: u64 = 0;
    let mut length_sum: u64 = 0;
    let mut mapq_counts = [0usize; 256];

    for record in records {
        length_sum += record.seq_len() as u64;

        if record.is_mapped() {
            mapped += 1;
            mapq_sum += record.mapq as u64;
            mapq_counts[record.mapq as usize] += 1;
        } else {
            unmapped += 1;
        }
    }

    let avg_mapq = if mapped > 0 {
        mapq_sum as f64 / mapped as f64
    } else {
        0.0
    };

    let avg_length = if total_reads > 0 {
        length_sum as f64 / total_reads as f64
    } else {
        0.0
    };

    let mapq_distribution: Vec<(u8, usize)> = mapq_counts
        .iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
        .map(|(val, &count)| (val as u8, count))
        .collect();

    SamStats {
        total_reads,
        mapped,
        unmapped,
        avg_mapq,
        avg_length,
        mapq_distribution,
    }
}

/// Parse a SAM file and return summary statistics without storing all records.
///
/// This is more memory-efficient than calling [`parse_sam`] followed by
/// [`sam_stats`] for large files, as it computes statistics in a streaming fashion.
pub fn sam_stats_from_path(path: impl AsRef<Path>) -> Result<SamStats> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut total_reads: usize = 0;
    let mut mapped: usize = 0;
    let mut unmapped: usize = 0;
    let mut mapq_sum: u64 = 0;
    let mut length_sum: u64 = 0;
    let mut mapq_counts = [0usize; 256];

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(CyaneaError::Io)?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('@') {
            continue;
        }

        let record = parse_sam_record(line, line_num + 1, path)?;
        total_reads += 1;
        length_sum += record.seq_len() as u64;

        if record.is_mapped() {
            mapped += 1;
            mapq_sum += record.mapq as u64;
            mapq_counts[record.mapq as usize] += 1;
        } else {
            unmapped += 1;
        }
    }

    let avg_mapq = if mapped > 0 {
        mapq_sum as f64 / mapped as f64
    } else {
        0.0
    };

    let avg_length = if total_reads > 0 {
        length_sum as f64 / total_reads as f64
    } else {
        0.0
    };

    let mapq_distribution: Vec<(u8, usize)> = mapq_counts
        .iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
        .map(|(val, &count)| (val as u8, count))
        .collect();

    Ok(SamStats {
        total_reads,
        mapped,
        unmapped,
        avg_mapq,
        avg_length,
        mapq_distribution,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_sam(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".sam").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    /// Minimal SAM with header and three alignment records.
    /// read1: 50 bp mapped, read2: 50 bp mapped, read3: 50 bp unmapped.
    const BASIC_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:248956422
@SQ\tSN:chr2\tLN:242193529
read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*
read2\t16\tchr1\t200\t30\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*
read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*
";

    #[test]
    fn test_parse_basic_sam() {
        let file = write_sam(BASIC_SAM);
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records.len(), 3);

        // First record: mapped, forward strand
        assert_eq!(records[0].qname, "read1");
        assert_eq!(records[0].flag, 0);
        assert_eq!(records[0].rname, "chr1");
        assert_eq!(records[0].pos, 100);
        assert_eq!(records[0].mapq, 60);
        assert_eq!(records[0].cigar, "50M");
        assert_eq!(records[0].seq_len(), 50);
        assert!(records[0].is_mapped());
        assert!(!records[0].is_unmapped());

        // Second record: mapped, reverse strand (FLAG 16)
        assert_eq!(records[1].qname, "read2");
        assert_eq!(records[1].flag, 16);
        assert_eq!(records[1].rname, "chr1");
        assert_eq!(records[1].pos, 200);
        assert_eq!(records[1].mapq, 30);
        assert!(records[1].is_mapped());

        // Third record: unmapped (FLAG 4)
        assert_eq!(records[2].qname, "read3");
        assert_eq!(records[2].flag, 4);
        assert_eq!(records[2].rname, "*");
        assert_eq!(records[2].pos, 0);
        assert_eq!(records[2].mapq, 0);
        assert_eq!(records[2].cigar, "*");
        assert!(records[2].is_unmapped());
        assert!(!records[2].is_mapped());
    }

    #[test]
    fn test_sam_stats_basic() {
        let file = write_sam(BASIC_SAM);
        let records = parse_sam(file.path()).unwrap();
        let stats = sam_stats(&records);

        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.mapped, 2);
        assert_eq!(stats.unmapped, 1);
        // avg_mapq = (60 + 30) / 2 = 45.0
        assert!((stats.avg_mapq - 45.0).abs() < f64::EPSILON);
        // All three reads have 50 bp sequences
        assert!((stats.avg_length - 50.0).abs() < f64::EPSILON);

        // MAPQ distribution: 30 -> 1, 60 -> 1
        assert_eq!(stats.mapq_distribution.len(), 2);
        assert_eq!(stats.mapq_distribution[0], (30, 1));
        assert_eq!(stats.mapq_distribution[1], (60, 1));
    }

    #[test]
    fn test_sam_stats_from_path() {
        let file = write_sam(BASIC_SAM);
        let stats = sam_stats_from_path(file.path()).unwrap();

        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.mapped, 2);
        assert_eq!(stats.unmapped, 1);
        assert!((stats.avg_mapq - 45.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_sam_header_only() {
        let file = write_sam(
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:248956422\n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_sam_empty_file() {
        let file = write_sam("");
        let records = parse_sam(file.path()).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_sam_stats_empty() {
        let stats = sam_stats(&[]);
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.mapped, 0);
        assert_eq!(stats.unmapped, 0);
        assert_eq!(stats.avg_mapq, 0.0);
        assert_eq!(stats.avg_length, 0.0);
        assert!(stats.mapq_distribution.is_empty());
    }

    #[test]
    fn test_sam_all_unmapped() {
        let file = write_sam(
            "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n\
             read2\t4\t*\t0\t0\t*\t*\t0\t0\tGGCC\t*\n",
        );
        let records = parse_sam(file.path()).unwrap();
        let stats = sam_stats(&records);

        assert_eq!(stats.total_reads, 2);
        assert_eq!(stats.mapped, 0);
        assert_eq!(stats.unmapped, 2);
        assert_eq!(stats.avg_mapq, 0.0); // no mapped reads, default to 0
        assert!((stats.avg_length - 4.0).abs() < f64::EPSILON);
        assert!(stats.mapq_distribution.is_empty());
    }

    #[test]
    fn test_sam_all_mapped() {
        let file = write_sam(
            "read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
             read2\t0\tchr1\t200\t60\t10M\t*\t0\t0\tGCTAGCTAGC\tIIIIIIIIII\n\
             read3\t16\tchr2\t300\t40\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n",
        );
        let records = parse_sam(file.path()).unwrap();
        let stats = sam_stats(&records);

        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.mapped, 3);
        assert_eq!(stats.unmapped, 0);
        // avg_mapq = (60 + 60 + 40) / 3 = 53.333...
        assert!((stats.avg_mapq - 160.0 / 3.0).abs() < 0.001);
        assert!((stats.avg_length - 10.0).abs() < f64::EPSILON);

        // MAPQ distribution: 40 -> 1, 60 -> 2
        assert_eq!(stats.mapq_distribution.len(), 2);
        assert_eq!(stats.mapq_distribution[0], (40, 1));
        assert_eq!(stats.mapq_distribution[1], (60, 2));
    }

    #[test]
    fn test_sam_too_few_columns() {
        let file = write_sam("read1\t0\tchr1\t100\t60\n");
        let result = parse_sam(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("11 tab-separated columns"));
    }

    #[test]
    fn test_sam_invalid_flag() {
        let file = write_sam(
            "read1\tNOTANUMBER\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\t*\n",
        );
        let result = parse_sam(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("invalid FLAG"));
    }

    #[test]
    fn test_sam_invalid_pos() {
        let file = write_sam(
            "read1\t0\tchr1\tXYZ\t60\t50M\t*\t0\t0\tACGT\t*\n",
        );
        let result = parse_sam(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("invalid POS"));
    }

    #[test]
    fn test_sam_invalid_mapq() {
        let file = write_sam(
            "read1\t0\tchr1\t100\t999\t50M\t*\t0\t0\tACGT\t*\n",
        );
        let result = parse_sam(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("invalid MAPQ"));
    }

    #[test]
    fn test_sam_file_not_found() {
        let result = parse_sam("/nonexistent/file.sam");
        assert!(result.is_err());
    }

    #[test]
    fn test_sam_with_optional_fields() {
        // SAM records can have optional fields after the 11 mandatory ones.
        // These should be silently ignored.
        let file = write_sam(
            "read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\tNM:i:2\tMD:Z:48A1\n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].qname, "read1");
        assert_eq!(records[0].mapq, 60);
    }

    #[test]
    fn test_sam_paired_end_flags() {
        // Paired-end reads have FLAG bits: 0x1 (paired), 0x40 (first), 0x80 (second)
        let file = write_sam(
            "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\n\
             read1\t147\tchr1\t200\t60\t50M\t=\t100\t-150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records.len(), 2);

        // Both should be mapped (FLAG 99 = 0x63, FLAG 147 = 0x93, neither has 0x4)
        assert!(records[0].is_mapped());
        assert!(records[1].is_mapped());
        assert_eq!(records[0].flag, 99);
        assert_eq!(records[1].flag, 147);
    }

    #[test]
    fn test_sam_seq_len_star() {
        // When SEQ is "*", seq_len should return 0
        let file = write_sam(
            "read1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq_len(), 0);
        assert_eq!(records[0].sequence, "*");
    }

    #[test]
    fn test_sam_mapq_255() {
        // MAPQ 255 means unavailable
        let file = write_sam(
            "read1\t0\tchr1\t100\t255\t50M\t*\t0\t0\tACGT\t*\n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records[0].mapq, 255);
    }

    #[test]
    fn test_sam_skip_blank_lines() {
        let file = write_sam(
            "\n\
             @HD\tVN:1.6\n\
             \n\
             read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n\
             \n\
             read2\t0\tchr1\t200\t40\t10M\t*\t0\t0\tGCTAGCTAGC\t*\n\
             \n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn test_sam_quality_string() {
        let file = write_sam(
            "read1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        );
        let records = parse_sam(file.path()).unwrap();
        assert_eq!(records[0].quality, "IIII");
    }

    #[test]
    fn test_sam_stats_mapq_distribution_sorted() {
        let file = write_sam(
            "r1\t0\tchr1\t100\t10\t5M\t*\t0\t0\tACGTA\t*\n\
             r2\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\t*\n\
             r3\t0\tchr1\t300\t10\t5M\t*\t0\t0\tACGTA\t*\n\
             r4\t0\tchr1\t400\t30\t5M\t*\t0\t0\tACGTA\t*\n\
             r5\t0\tchr1\t500\t60\t5M\t*\t0\t0\tACGTA\t*\n",
        );
        let records = parse_sam(file.path()).unwrap();
        let stats = sam_stats(&records);

        assert_eq!(stats.total_reads, 5);
        assert_eq!(stats.mapped, 5);
        // MAPQ distribution sorted by value: 10->2, 30->1, 60->2
        assert_eq!(stats.mapq_distribution.len(), 3);
        assert_eq!(stats.mapq_distribution[0], (10, 2));
        assert_eq!(stats.mapq_distribution[1], (30, 1));
        assert_eq!(stats.mapq_distribution[2], (60, 2));

        // avg_mapq = (10+60+10+30+60)/5 = 34.0
        assert!((stats.avg_mapq - 34.0).abs() < f64::EPSILON);
    }
}
