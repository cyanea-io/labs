//! BLAST tabular output parser (`-outfmt 6` and `-outfmt 7`).
//!
//! Parses the standard 12-column BLAST tabular format produced by
//! `blastn`, `blastp`, `blastx`, etc. with `-outfmt 6` (no comments)
//! or `-outfmt 7` (with `#` comment lines, which are skipped).

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

/// A single BLAST hit from tabular output.
///
/// Corresponds to the standard 12 columns of `-outfmt 6`:
/// `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`
#[derive(Debug, Clone)]
pub struct BlastRecord {
    /// Query sequence ID.
    pub query_id: String,
    /// Subject (database) sequence ID.
    pub subject_id: String,
    /// Percentage of identical matches.
    pub pct_identity: f64,
    /// Alignment length.
    pub alignment_length: u64,
    /// Number of mismatches.
    pub mismatches: u64,
    /// Number of gap openings.
    pub gap_opens: u64,
    /// Start of alignment in query (1-based).
    pub query_start: u64,
    /// End of alignment in query (1-based).
    pub query_end: u64,
    /// Start of alignment in subject (1-based).
    pub subject_start: u64,
    /// End of alignment in subject (1-based).
    pub subject_end: u64,
    /// Expect value.
    pub evalue: f64,
    /// Bit score.
    pub bit_score: f64,
}

/// Summary statistics for a BLAST tabular output file.
#[derive(Debug, Clone)]
pub struct BlastStats {
    /// Total number of hits.
    pub hit_count: u64,
    /// Number of unique query sequences.
    pub unique_queries: u64,
    /// Number of unique subject sequences.
    pub unique_subjects: u64,
    /// Average percent identity across all hits.
    pub avg_identity: f64,
    /// Average E-value across all hits.
    pub avg_evalue: f64,
}

/// Parse a BLAST tabular output file (`-outfmt 6` or `-outfmt 7`).
///
/// Lines starting with `#` (outfmt 7 comment/header lines) are skipped.
/// Returns an error if the file is empty (no hits) or contains malformed lines.
pub fn parse_blast(path: impl AsRef<Path>) -> Result<Vec<BlastRecord>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let record = parse_blast_line(line, line_num + 1, path)?;
        records.push(record);
    }

    if records.is_empty() {
        return Err(CyaneaError::Parse(format!(
            "{}: no BLAST hits found",
            path.display()
        )));
    }

    Ok(records)
}

/// Compute summary statistics from a BLAST tabular output file.
pub fn blast_stats(path: impl AsRef<Path>) -> Result<BlastStats> {
    let records = parse_blast(path)?;
    let mut queries = std::collections::HashSet::new();
    let mut subjects = std::collections::HashSet::new();
    let mut sum_identity = 0.0;
    let mut sum_evalue = 0.0;

    for rec in &records {
        queries.insert(&rec.query_id);
        subjects.insert(&rec.subject_id);
        sum_identity += rec.pct_identity;
        sum_evalue += rec.evalue;
    }

    let n = records.len() as f64;
    Ok(BlastStats {
        hit_count: records.len() as u64,
        unique_queries: queries.len() as u64,
        unique_subjects: subjects.len() as u64,
        avg_identity: sum_identity / n,
        avg_evalue: sum_evalue / n,
    })
}

/// Parse a single BLAST tabular line.
fn parse_blast_line(line: &str, line_num: usize, path: &Path) -> Result<BlastRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected 12 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let parse_f64 = |idx: usize, name: &str| -> Result<f64> {
        fields[idx].parse::<f64>().map_err(|_| {
            CyaneaError::Parse(format!(
                "{}: line {}: invalid {} '{}'",
                path.display(),
                line_num,
                name,
                fields[idx]
            ))
        })
    };

    let parse_u64 = |idx: usize, name: &str| -> Result<u64> {
        fields[idx].parse::<u64>().map_err(|_| {
            CyaneaError::Parse(format!(
                "{}: line {}: invalid {} '{}'",
                path.display(),
                line_num,
                name,
                fields[idx]
            ))
        })
    };

    Ok(BlastRecord {
        query_id: fields[0].to_string(),
        subject_id: fields[1].to_string(),
        pct_identity: parse_f64(2, "pct_identity")?,
        alignment_length: parse_u64(3, "alignment_length")?,
        mismatches: parse_u64(4, "mismatches")?,
        gap_opens: parse_u64(5, "gap_opens")?,
        query_start: parse_u64(6, "query_start")?,
        query_end: parse_u64(7, "query_end")?,
        subject_start: parse_u64(8, "subject_start")?,
        subject_end: parse_u64(9, "subject_end")?,
        evalue: parse_f64(10, "evalue")?,
        bit_score: parse_f64(11, "bit_score")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_blast(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".blast").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn blast_parse_simple() {
        let file = write_blast(
            "query1\tsubj1\t99.5\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200.0\n\
             query1\tsubj2\t85.0\t80\t12\t0\t1\t80\t50\t129\t1e-20\t150.0\n\
             query2\tsubj1\t92.3\t60\t5\t1\t10\t69\t200\t259\t1e-30\t170.5\n",
        );

        let records = parse_blast(file.path()).unwrap();
        assert_eq!(records.len(), 3);

        assert_eq!(records[0].query_id, "query1");
        assert_eq!(records[0].subject_id, "subj1");
        assert!((records[0].pct_identity - 99.5).abs() < f64::EPSILON);
        assert_eq!(records[0].alignment_length, 100);
        assert_eq!(records[0].mismatches, 0);
        assert_eq!(records[0].gap_opens, 0);
        assert_eq!(records[0].query_start, 1);
        assert_eq!(records[0].query_end, 100);
        assert_eq!(records[0].subject_start, 1);
        assert_eq!(records[0].subject_end, 100);
        assert!((records[0].evalue - 1e-50).abs() < 1e-60);
        assert!((records[0].bit_score - 200.0).abs() < f64::EPSILON);

        assert_eq!(records[2].query_id, "query2");
        assert_eq!(records[2].gap_opens, 1);
    }

    #[test]
    fn blast_outfmt7_comments() {
        let file = write_blast(
            "# BLASTN 2.14.0+\n\
             # Query: query1\n\
             # Database: nr\n\
             # Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n\
             # 2 hits found\n\
             query1\tsubj1\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-25\t100.0\n\
             query1\tsubj2\t95.0\t50\t2\t1\t1\t50\t100\t149\t1e-15\t80.0\n",
        );

        let records = parse_blast(file.path()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].query_id, "query1");
        assert_eq!(records[1].subject_id, "subj2");
    }

    #[test]
    fn blast_stats_computed() {
        let file = write_blast(
            "q1\ts1\t90.0\t100\t10\t0\t1\t100\t1\t100\t1e-40\t180.0\n\
             q1\ts2\t80.0\t80\t16\t0\t1\t80\t1\t80\t1e-20\t120.0\n\
             q2\ts1\t95.0\t50\t2\t1\t1\t50\t1\t50\t1e-30\t160.0\n",
        );

        let stats = blast_stats(file.path()).unwrap();
        assert_eq!(stats.hit_count, 3);
        assert_eq!(stats.unique_queries, 2);
        assert_eq!(stats.unique_subjects, 2);
        assert!((stats.avg_identity - 88.333_333_333_333_33).abs() < 1e-6);
    }

    #[test]
    fn blast_empty_file_error() {
        let file = write_blast("# BLASTN output\n# 0 hits found\n");
        let result = parse_blast(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn blast_malformed_line() {
        let file = write_blast("query1\tsubj1\t99.0\n");
        let result = parse_blast(file.path());
        assert!(result.is_err());
    }
}
