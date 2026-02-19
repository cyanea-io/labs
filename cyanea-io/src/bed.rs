//! BED (Browser Extensible Data) parser.
//!
//! Parses BED3–BED6 files into [`GenomicInterval`] records.
//! BED format uses 0-based, half-open coordinates `[start, end)`,
//! which matches [`GenomicInterval`] directly.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::genomic::{GenomicInterval, Strand};

/// A BED record with optional name and score fields (BED3–BED6).
#[derive(Debug, Clone)]
pub struct BedRecord {
    pub interval: GenomicInterval,
    /// Column 4: feature name (BED4+).
    pub name: Option<String>,
    /// Column 5: score 0–1000 (BED5+).
    pub score: Option<u32>,
}

/// Parse a BED file and return all records as [`BedRecord`]s.
///
/// Supports BED3 through BED6. Lines starting with `#`, `track`, or `browser`
/// are treated as headers and skipped. Empty lines are also skipped.
pub fn parse_bed(path: impl AsRef<Path>) -> Result<Vec<BedRecord>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    // Read all data lines (skip headers)
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

    // Parse data lines (optionally in parallel)
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data_lines
            .par_iter()
            .map(|(line_num, line)| parse_bed_line(line, *line_num, path))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    data_lines
        .iter()
        .map(|(line_num, line)| parse_bed_line(line, *line_num, path))
        .collect()
}

/// Parse a BED file and return only the genomic intervals (discards name/score).
pub fn parse_bed_intervals(path: impl AsRef<Path>) -> Result<Vec<GenomicInterval>> {
    parse_bed(path).map(|records| records.into_iter().map(|r| r.interval).collect())
}

/// Summary statistics for a BED file.
#[derive(Debug, Clone)]
pub struct BedStats {
    pub record_count: u64,
    pub total_bases: u64,
    pub chromosomes: Vec<String>,
}

/// Parse a BED file and return summary statistics.
pub fn bed_stats(path: impl AsRef<Path>) -> Result<BedStats> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut record_count: u64 = 0;
    let mut total_bases: u64 = 0;
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

        let record = parse_bed_line(line, line_num + 1, path)?;
        record_count += 1;
        total_bases += record.interval.len();

        if seen_chroms.insert(record.interval.chrom.clone()) {
            chroms.push(record.interval.chrom);
        }
    }

    Ok(BedStats {
        record_count,
        total_bases,
        chromosomes: chroms,
    })
}

/// Parse BED text from a string and return all records.
///
/// Behaves like [`parse_bed`] but reads from an in-memory string instead of a file.
pub fn parse_bed_str(text: &str) -> Result<Vec<BedRecord>> {
    let dummy = Path::new("<string>");
    text.lines()
        .enumerate()
        .filter(|(_, line)| {
            let t = line.trim();
            !t.is_empty()
                && !t.starts_with('#')
                && !t.starts_with("track")
                && !t.starts_with("browser")
        })
        .map(|(i, line)| parse_bed_line(line.trim(), i + 1, dummy))
        .collect()
}

/// Parse a single BED line into a BedRecord.
fn parse_bed_line(line: &str, line_num: usize, path: &Path) -> Result<BedRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected at least 3 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let chrom = fields[0].to_string();

    let start: u64 = fields[1].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid start coordinate '{}'",
            path.display(),
            line_num,
            fields[1]
        ))
    })?;

    let end: u64 = fields[2].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid end coordinate '{}'",
            path.display(),
            line_num,
            fields[2]
        ))
    })?;

    if start > end {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: start ({}) > end ({})",
            path.display(),
            line_num,
            start,
            end
        )));
    }

    let name = fields.get(3).and_then(|s| {
        if *s == "." { None } else { Some(s.to_string()) }
    });

    let score = fields.get(4).and_then(|s| {
        if *s == "." { None } else { s.parse().ok() }
    });

    let strand = fields.get(5).map_or(Strand::Unknown, |s| match *s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    });

    let interval = GenomicInterval::with_strand(chrom, start, end, strand)?;

    Ok(BedRecord {
        interval,
        name,
        score,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_bed(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".bed").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_parse_bed3() {
        let file = write_bed("chr1\t100\t200\nchr1\t300\t400\nchr2\t500\t600\n");
        let records = parse_bed(file.path()).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].interval.chrom, "chr1");
        assert_eq!(records[0].interval.start, 100);
        assert_eq!(records[0].interval.end, 200);
        assert_eq!(records[0].interval.strand, Strand::Unknown);
        assert!(records[0].name.is_none());
        assert!(records[0].score.is_none());
    }

    #[test]
    fn test_parse_bed6() {
        let file = write_bed("chr1\t100\t200\tgene1\t500\t+\nchr1\t300\t400\tgene2\t0\t-\n");
        let records = parse_bed(file.path()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[0].score, Some(500));
        assert_eq!(records[0].interval.strand, Strand::Forward);
        assert_eq!(records[1].interval.strand, Strand::Reverse);
    }

    #[test]
    fn test_parse_bed_intervals() {
        let file = write_bed("chr1\t0\t100\nchr1\t200\t300\n");
        let intervals = parse_bed_intervals(file.path()).unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].len(), 100);
        assert_eq!(intervals[1].start, 200);
    }

    #[test]
    fn test_bed_stats() {
        let file = write_bed("chr1\t0\t100\nchr1\t200\t350\nchr2\t0\t50\n");
        let stats = bed_stats(file.path()).unwrap();
        assert_eq!(stats.record_count, 3);
        assert_eq!(stats.total_bases, 300); // 100 + 150 + 50
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn test_bed_skip_headers() {
        let file = write_bed(
            "# comment\n\
             track name=test\n\
             browser position chr1:1-1000\n\
             chr1\t100\t200\n",
        );
        let records = parse_bed(file.path()).unwrap();
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn test_bed_empty_file() {
        let file = write_bed("# header only\n");
        let records = parse_bed(file.path()).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_bed_invalid_start_end() {
        let file = write_bed("chr1\t200\t100\n");
        let result = parse_bed(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_bed_too_few_columns() {
        let file = write_bed("chr1\t100\n");
        let result = parse_bed(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_bed_str() {
        let text = "chr1\t100\t200\tgene1\t500\t+\nchr2\t300\t400\n";
        let records = super::parse_bed_str(text).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].interval.chrom, "chr1");
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[1].interval.chrom, "chr2");
    }

    #[test]
    fn test_bed_file_not_found() {
        let result = parse_bed("/nonexistent/file.bed");
        assert!(result.is_err());
    }
}
