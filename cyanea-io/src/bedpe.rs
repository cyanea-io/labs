//! BEDPE (paired-end BED) parser.
//!
//! Parses BEDPE files into paired [`GenomicInterval`] records.
//! BEDPE format uses 0-based, half-open coordinates `[start, end)`.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::genomic::{GenomicInterval, Strand};

/// A BEDPE record with two genomic intervals and optional name/score.
#[derive(Debug, Clone)]
pub struct BedpeRecord {
    /// First interval (chrom1:start1-end1).
    pub interval1: GenomicInterval,
    /// Second interval (chrom2:start2-end2).
    pub interval2: GenomicInterval,
    /// Column 7: feature name (optional).
    pub name: Option<String>,
    /// Column 8: score 0–1000 (optional).
    pub score: Option<u32>,
}

/// Parse a BEDPE file and return all records.
///
/// Supports BEDPE6 (6 required columns) through BEDPE10 (all optional columns).
/// Lines starting with `#`, `track`, or `browser` are treated as headers
/// and skipped. Empty lines are also skipped.
pub fn parse_bedpe(path: impl AsRef<Path>) -> Result<Vec<BedpeRecord>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut data_lines: Vec<(usize, String)> = Vec::new();
    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let trimmed = line.trim().to_string();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }
        data_lines.push((line_num + 1, trimmed));
    }

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data_lines
            .par_iter()
            .map(|(line_num, line)| parse_bedpe_line(line, *line_num, path))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    data_lines
        .iter()
        .map(|(line_num, line)| parse_bedpe_line(line, *line_num, path))
        .collect()
}

/// Summary statistics for a BEDPE file.
#[derive(Debug, Clone)]
pub struct BedpeStats {
    pub record_count: u64,
    pub total_span: u64,
    pub inter_chromosomal: u64,
    pub intra_chromosomal: u64,
    pub chromosomes: Vec<String>,
}

/// Parse a BEDPE file and return summary statistics.
pub fn bedpe_stats(path: impl AsRef<Path>) -> Result<BedpeStats> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut record_count: u64 = 0;
    let mut total_span: u64 = 0;
    let mut inter_chromosomal: u64 = 0;
    let mut intra_chromosomal: u64 = 0;
    let mut chroms = Vec::new();
    let mut seen_chroms = std::collections::HashSet::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(CyaneaError::Io)?;
        let line = line.trim();

        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let record = parse_bedpe_line(line, line_num + 1, path)?;
        record_count += 1;
        total_span += record.interval1.len() + record.interval2.len();

        if record.interval1.chrom == record.interval2.chrom {
            intra_chromosomal += 1;
        } else {
            inter_chromosomal += 1;
        }

        if seen_chroms.insert(record.interval1.chrom.clone()) {
            chroms.push(record.interval1.chrom.clone());
        }
        if seen_chroms.insert(record.interval2.chrom.clone()) {
            chroms.push(record.interval2.chrom.clone());
        }
    }

    Ok(BedpeStats {
        record_count,
        total_span,
        inter_chromosomal,
        intra_chromosomal,
        chromosomes: chroms,
    })
}

/// Parse a single BEDPE line into a BedpeRecord.
fn parse_bedpe_line(line: &str, line_num: usize, path: &Path) -> Result<BedpeRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected at least 6 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let chrom1 = fields[0].to_string();
    let start1: u64 = fields[1].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid start1 coordinate '{}'",
            path.display(),
            line_num,
            fields[1]
        ))
    })?;
    let end1: u64 = fields[2].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid end1 coordinate '{}'",
            path.display(),
            line_num,
            fields[2]
        ))
    })?;
    if start1 > end1 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: start1 ({}) > end1 ({})",
            path.display(),
            line_num,
            start1,
            end1
        )));
    }

    let chrom2 = fields[3].to_string();
    let start2: u64 = fields[4].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid start2 coordinate '{}'",
            path.display(),
            line_num,
            fields[4]
        ))
    })?;
    let end2: u64 = fields[5].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid end2 coordinate '{}'",
            path.display(),
            line_num,
            fields[5]
        ))
    })?;
    if start2 > end2 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: start2 ({}) > end2 ({})",
            path.display(),
            line_num,
            start2,
            end2
        )));
    }

    let name = fields.get(6).and_then(|s| {
        if *s == "." { None } else { Some(s.to_string()) }
    });

    let score = fields.get(7).and_then(|s| {
        if *s == "." { None } else { s.parse().ok() }
    });

    let strand1 = fields.get(8).map_or(Strand::Unknown, |s| match *s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    });

    let strand2 = fields.get(9).map_or(Strand::Unknown, |s| match *s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    });

    let interval1 = GenomicInterval::with_strand(chrom1, start1, end1, strand1)?;
    let interval2 = GenomicInterval::with_strand(chrom2, start2, end2, strand2)?;

    Ok(BedpeRecord {
        interval1,
        interval2,
        name,
        score,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_bedpe(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".bedpe").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_parse_bedpe6() {
        let file = write_bedpe(
            "chr1\t100\t200\tchr1\t300\t400\n\
             chr2\t500\t600\tchr3\t700\t800\n",
        );
        let records = parse_bedpe(file.path()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].interval1.chrom, "chr1");
        assert_eq!(records[0].interval1.start, 100);
        assert_eq!(records[0].interval1.end, 200);
        assert_eq!(records[0].interval2.chrom, "chr1");
        assert_eq!(records[0].interval2.start, 300);
        assert_eq!(records[0].interval2.end, 400);
        assert!(records[0].name.is_none());
        assert!(records[0].score.is_none());
    }

    #[test]
    fn test_parse_bedpe10() {
        let file = write_bedpe(
            "chr1\t100\t200\tchr1\t300\t400\tpair1\t500\t+\t-\n",
        );
        let records = parse_bedpe(file.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, Some("pair1".to_string()));
        assert_eq!(records[0].score, Some(500));
        assert_eq!(records[0].interval1.strand, Strand::Forward);
        assert_eq!(records[0].interval2.strand, Strand::Reverse);
    }

    #[test]
    fn test_bedpe_skip_headers() {
        let file = write_bedpe(
            "# comment\n\
             track name=test\n\
             browser position chr1:1-1000\n\
             chr1\t100\t200\tchr1\t300\t400\n",
        );
        let records = parse_bedpe(file.path()).unwrap();
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn test_bedpe_too_few_columns() {
        let file = write_bedpe("chr1\t100\t200\tchr1\t300\n");
        let result = parse_bedpe(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_bedpe_invalid_coordinates() {
        let file = write_bedpe("chr1\t200\t100\tchr1\t300\t400\n");
        assert!(parse_bedpe(file.path()).is_err());

        let file2 = write_bedpe("chr1\t100\t200\tchr1\t400\t300\n");
        assert!(parse_bedpe(file2.path()).is_err());
    }

    #[test]
    fn test_bedpe_stats() {
        let file = write_bedpe(
            "chr1\t0\t100\tchr1\t200\t300\n\
             chr1\t0\t50\tchr2\t0\t50\n\
             chr3\t0\t25\tchr3\t50\t75\n",
        );
        let stats = bedpe_stats(file.path()).unwrap();
        assert_eq!(stats.record_count, 3);
        assert_eq!(stats.total_span, 350); // (100+100) + (50+50) + (25+25)
        assert_eq!(stats.intra_chromosomal, 2);
        assert_eq!(stats.inter_chromosomal, 1);
        assert_eq!(stats.chromosomes.len(), 3);
    }

    #[test]
    fn test_bedpe_empty_file() {
        let file = write_bedpe("# header only\n");
        let records = parse_bedpe(file.path()).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_bedpe_file_not_found() {
        let result = parse_bedpe("/nonexistent/file.bedpe");
        assert!(result.is_err());
    }

    #[test]
    fn test_bedpe_dot_name_score() {
        let file = write_bedpe(
            "chr1\t100\t200\tchr1\t300\t400\t.\t.\t+\t-\n",
        );
        let records = parse_bedpe(file.path()).unwrap();
        assert!(records[0].name.is_none());
        assert!(records[0].score.is_none());
    }

    #[test]
    fn test_bedpe_stats_empty() {
        let file = write_bedpe("");
        let stats = bedpe_stats(file.path()).unwrap();
        assert_eq!(stats.record_count, 0);
        assert_eq!(stats.total_span, 0);
        assert_eq!(stats.inter_chromosomal, 0);
        assert_eq!(stats.intra_chromosomal, 0);
    }
}
