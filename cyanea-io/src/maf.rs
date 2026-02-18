//! MAF (Multiple Alignment Format) parser.
//!
//! Parses both UCSC Multiple Alignment Format and pairwise MAF output
//! from LAST/minimap2. Both use the same `a`/`s`/`q` block structure.
//!
//! An MAF file consists of alignment blocks beginning with `a` (alignment)
//! lines, followed by `s` (sequence) lines with aligned text.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

/// A single MAF alignment block.
///
/// Each block starts with an `a` line (optionally containing a score)
/// and contains one or more `s` lines describing aligned sequences.
#[derive(Debug, Clone)]
pub struct MafBlock {
    /// Alignment score from the `a` line, if present.
    pub score: Option<f64>,
    /// Sequences participating in this alignment block.
    pub sequences: Vec<MafSequence>,
}

/// An aligned sequence within a MAF block.
#[derive(Debug, Clone)]
pub struct MafSequence {
    /// Source sequence name (e.g. `hg38.chr1`).
    pub src: String,
    /// 0-based start position in the source sequence.
    pub start: u64,
    /// Number of non-gap bases in this alignment.
    pub size: u64,
    /// Strand: `+` or `-`.
    pub strand: char,
    /// Total length of the source sequence.
    pub src_size: u64,
    /// Alignment text including gap characters (`-`).
    pub text: String,
}

/// Summary statistics for an MAF file.
#[derive(Debug, Clone)]
pub struct MafStats {
    /// Number of alignment blocks.
    pub block_count: u64,
    /// Total aligned bases (sum of all sequence sizes).
    pub total_aligned_bases: u64,
    /// Unique species/source names across all blocks.
    pub species: Vec<String>,
}

/// Parse an MAF file and return all alignment blocks.
///
/// Handles both multi-species (UCSC) and pairwise (LAST/minimap2) MAF files.
/// Lines starting with `#` are treated as comments and skipped.
pub fn parse_maf(path: impl AsRef<Path>) -> Result<Vec<MafBlock>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    let mut blocks = Vec::new();
    let mut current_block: Option<MafBlock> = None;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let line = line.trim();

        if line.is_empty() {
            // Empty line ends current block
            if let Some(block) = current_block.take() {
                if !block.sequences.is_empty() {
                    blocks.push(block);
                }
            }
            continue;
        }

        if line.starts_with('#') {
            continue;
        }

        if line.starts_with("a ") || line == "a" {
            // Finish previous block if any
            if let Some(block) = current_block.take() {
                if !block.sequences.is_empty() {
                    blocks.push(block);
                }
            }
            let score = parse_score(line);
            current_block = Some(MafBlock {
                score,
                sequences: Vec::new(),
            });
        } else if line.starts_with("s ") {
            let seq = parse_s_line(line, line_num + 1, path)?;
            if let Some(ref mut block) = current_block {
                block.sequences.push(seq);
            }
        }
        // Skip q, i, e, and other line types
    }

    // Don't forget the last block
    if let Some(block) = current_block.take() {
        if !block.sequences.is_empty() {
            blocks.push(block);
        }
    }

    Ok(blocks)
}

/// Compute summary statistics from an MAF file.
pub fn maf_stats(path: impl AsRef<Path>) -> Result<MafStats> {
    let blocks = parse_maf(path)?;
    let mut total_aligned_bases: u64 = 0;
    let mut species_set = std::collections::HashSet::new();
    let mut species_order = Vec::new();

    for block in &blocks {
        for seq in &block.sequences {
            total_aligned_bases += seq.size;
            if species_set.insert(seq.src.clone()) {
                species_order.push(seq.src.clone());
            }
        }
    }

    Ok(MafStats {
        block_count: blocks.len() as u64,
        total_aligned_bases,
        species: species_order,
    })
}

/// Extract score from an `a` line.
///
/// Format: `a score=1234.5 ...`
fn parse_score(line: &str) -> Option<f64> {
    for part in line.split_whitespace() {
        if let Some(val) = part.strip_prefix("score=") {
            return val.parse().ok();
        }
    }
    None
}

/// Parse an `s` (sequence) line.
///
/// Format: `s src start size strand srcSize text`
fn parse_s_line(line: &str, line_num: usize, path: &Path) -> Result<MafSequence> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 7 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected 7 fields in 's' line, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let src = fields[1].to_string();

    let start: u64 = fields[2].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid start '{}'",
            path.display(),
            line_num,
            fields[2]
        ))
    })?;

    let size: u64 = fields[3].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid size '{}'",
            path.display(),
            line_num,
            fields[3]
        ))
    })?;

    let strand = fields[4].chars().next().unwrap_or('+');
    if strand != '+' && strand != '-' {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: invalid strand '{}'",
            path.display(),
            line_num,
            fields[4]
        )));
    }

    let src_size: u64 = fields[5].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid srcSize '{}'",
            path.display(),
            line_num,
            fields[5]
        ))
    })?;

    let text = fields[6].to_string();

    Ok(MafSequence {
        src,
        start,
        size,
        strand,
        src_size,
        text,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_maf(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".maf").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn maf_parse_two_blocks() {
        let file = write_maf(
            "##maf version=1 scoring=tba.v8\n\
             \n\
             a score=23262.0\n\
             s hg38.chr7       27578828 38 + 159345973 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
             s panTro4.chr6    28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
             s ponAbe2.chr6    28594507 38 + 174210431 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
             \n\
             a score=5000.0\n\
             s hg38.chr7       27699739 6 + 159345973 TAAAGA\n\
             s panTro4.chr6    28862317 6 + 161576975 TAAAGA\n\
             \n",
        );

        let blocks = parse_maf(file.path()).unwrap();
        assert_eq!(blocks.len(), 2);

        assert!((blocks[0].score.unwrap() - 23262.0).abs() < f64::EPSILON);
        assert_eq!(blocks[0].sequences.len(), 3);
        assert_eq!(blocks[0].sequences[0].src, "hg38.chr7");
        assert_eq!(blocks[0].sequences[0].start, 27578828);
        assert_eq!(blocks[0].sequences[0].size, 38);
        assert_eq!(blocks[0].sequences[0].strand, '+');
        assert_eq!(blocks[0].sequences[0].src_size, 159345973);

        assert_eq!(blocks[1].sequences.len(), 2);
        assert!((blocks[1].score.unwrap() - 5000.0).abs() < f64::EPSILON);
    }

    #[test]
    fn maf_pairwise() {
        let file = write_maf(
            "a score=100.0\n\
             s ref.chr1   100 50 + 1000 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n\
             s query.ctg1   0 50 + 500  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n\
             \n",
        );

        let blocks = parse_maf(file.path()).unwrap();
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].sequences.len(), 2);
        assert_eq!(blocks[0].sequences[0].src, "ref.chr1");
        assert_eq!(blocks[0].sequences[1].src, "query.ctg1");
    }

    #[test]
    fn maf_minus_strand() {
        let file = write_maf(
            "a score=200.0\n\
             s hg38.chr1   1000 50 + 248956422 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n\
             s mm10.chr5   5000 50 - 151834684 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n\
             \n",
        );

        let blocks = parse_maf(file.path()).unwrap();
        assert_eq!(blocks[0].sequences[0].strand, '+');
        assert_eq!(blocks[0].sequences[1].strand, '-');
        assert_eq!(blocks[0].sequences[1].src_size, 151834684);
    }

    #[test]
    fn maf_stats_computed() {
        let file = write_maf(
            "a score=100.0\n\
             s species1.chr1  0 30 + 1000 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\
             s species2.chr2  0 30 + 2000 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\
             \n\
             a score=200.0\n\
             s species1.chr1  50 20 + 1000 AAAAAAAAAAAAAAAAAAAA\n\
             s species3.chr3  10 20 + 3000 AAAAAAAAAAAAAAAAAAAA\n\
             \n",
        );

        let stats = maf_stats(file.path()).unwrap();
        assert_eq!(stats.block_count, 2);
        assert_eq!(stats.total_aligned_bases, 100); // 30+30+20+20
        assert_eq!(stats.species.len(), 3);
        assert!(stats.species.contains(&"species1.chr1".to_string()));
        assert!(stats.species.contains(&"species2.chr2".to_string()));
        assert!(stats.species.contains(&"species3.chr3".to_string()));
    }

    #[test]
    fn maf_empty_block_skipped() {
        let file = write_maf(
            "a score=0\n\
             \n\
             a score=100.0\n\
             s src1 0 10 + 100 AAAAAAAAAA\n\
             \n",
        );

        let blocks = parse_maf(file.path()).unwrap();
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].sequences.len(), 1);
    }
}
