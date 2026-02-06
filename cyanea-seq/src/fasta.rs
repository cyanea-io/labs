//! FASTA/FASTQ parsing and statistics.

use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use needletail::parse_fastx_file;

/// Summary statistics for a FASTA/FASTQ file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FastaStats {
    pub sequence_count: u64,
    pub total_bases: u64,
    pub gc_content: f64,
    pub avg_length: f64,
}

/// Parse a FASTA/FASTQ file and compute summary statistics.
pub fn parse_fasta_stats(path: impl AsRef<Path>) -> Result<FastaStats> {
    let path = path.as_ref();
    let mut reader = parse_fastx_file(path).map_err(|e| CyaneaError::Parse(e.to_string()))?;

    let mut sequence_count: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut gc_count: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| CyaneaError::Parse(e.to_string()))?;
        let seq = record.seq();

        sequence_count += 1;
        total_bases += seq.len() as u64;

        for &base in seq.iter() {
            match base {
                b'G' | b'g' | b'C' | b'c' => gc_count += 1,
                _ => {}
            }
        }
    }

    let gc_content = if total_bases > 0 {
        (gc_count as f64 / total_bases as f64) * 100.0
    } else {
        0.0
    };

    let avg_length = if sequence_count > 0 {
        total_bases as f64 / sequence_count as f64
    } else {
        0.0
    };

    Ok(FastaStats {
        sequence_count,
        total_bases,
        gc_content,
        avg_length,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_fasta_parsing() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ATCGATCG").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "GCGCGCGC").unwrap();
        file.flush().unwrap();

        let stats = parse_fasta_stats(file.path()).unwrap();
        assert_eq!(stats.sequence_count, 2);
        assert_eq!(stats.total_bases, 16);
        assert!((stats.gc_content - 75.0).abs() < 0.01);
        assert!((stats.avg_length - 8.0).abs() < 0.01);
    }

    #[test]
    fn test_empty_fasta() {
        let file = NamedTempFile::new().unwrap();

        // needletail returns an error for empty files
        let result = parse_fasta_stats(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_fasta_file_not_found() {
        let result = parse_fasta_stats("/nonexistent/file.fasta");
        assert!(result.is_err());
    }
}
