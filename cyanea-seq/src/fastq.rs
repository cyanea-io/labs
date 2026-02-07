//! FASTQ record type and parsing.
//!
//! [`FastqRecord`] wraps a [`DnaSequence`] with quality scores and metadata.
//! Parsing uses needletail for streaming FASTQ reading.

use std::path::Path;

use cyanea_core::{Annotated, CyaneaError, Result, Sequence, Summarizable};

use crate::quality::{PhredEncoding, QualityScores};
use crate::types::DnaSequence;

/// A single FASTQ record: name, sequence, and quality scores.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FastqRecord {
    name: String,
    description: Option<String>,
    sequence: DnaSequence,
    quality: QualityScores,
}

impl FastqRecord {
    /// Create a new FASTQ record.
    ///
    /// Returns an error if the sequence and quality lengths don't match.
    pub fn new(
        name: String,
        description: Option<String>,
        sequence: DnaSequence,
        quality: QualityScores,
    ) -> Result<Self> {
        if sequence.len() != quality.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence length ({}) does not match quality length ({})",
                sequence.len(),
                quality.len()
            )));
        }
        Ok(Self {
            name,
            description,
            sequence,
            quality,
        })
    }

    /// The underlying DNA sequence.
    pub fn sequence(&self) -> &DnaSequence {
        &self.sequence
    }

    /// The quality scores.
    pub fn quality(&self) -> &QualityScores {
        &self.quality
    }
}

impl Sequence for FastqRecord {
    fn as_bytes(&self) -> &[u8] {
        self.sequence.as_bytes()
    }
}

impl Annotated for FastqRecord {
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
}

impl Summarizable for FastqRecord {
    fn summary(&self) -> String {
        format!(
            "FASTQ {} ({} bp, mean Q{:.1})",
            self.name,
            self.sequence.len(),
            self.quality.mean()
        )
    }
}

/// Summary statistics for a FASTQ file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FastqStats {
    pub sequence_count: u64,
    pub total_bases: u64,
    pub gc_content: f64,
    pub avg_length: f64,
    pub mean_quality: f64,
    pub q20_fraction: f64,
    pub q30_fraction: f64,
}

/// Parse a FASTQ file into a vector of [`FastqRecord`]s.
///
/// Uses Phred+33 encoding for quality scores. All sequences are treated as DNA.
pub fn parse_fastq_file(path: impl AsRef<Path>) -> Result<Vec<FastqRecord>> {
    let path = path.as_ref();
    let mut reader =
        needletail::parse_fastx_file(path).map_err(|e| CyaneaError::Parse(e.to_string()))?;

    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|e| CyaneaError::Parse(e.to_string()))?;

        let raw_id = std::str::from_utf8(record.id())
            .map_err(|e| CyaneaError::Parse(e.to_string()))?;

        // Split "name description" on first whitespace
        let (name, description) = match raw_id.split_once(char::is_whitespace) {
            Some((n, d)) => (n.to_string(), Some(d.to_string())),
            None => (raw_id.to_string(), None),
        };

        let sequence = DnaSequence::new(record.seq())?;

        let qual_bytes = record
            .qual()
            .ok_or_else(|| CyaneaError::Parse("missing quality scores".into()))?;
        let quality = QualityScores::from_ascii(qual_bytes, PhredEncoding::Phred33)?;

        records.push(FastqRecord::new(name, description, sequence, quality)?);
    }

    Ok(records)
}

/// Parse a FASTQ file and compute summary statistics without storing records.
pub fn parse_fastq_stats(path: impl AsRef<Path>) -> Result<FastqStats> {
    let path = path.as_ref();
    let mut reader =
        needletail::parse_fastx_file(path).map_err(|e| CyaneaError::Parse(e.to_string()))?;

    let mut sequence_count: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut gc_count: u64 = 0;
    let mut quality_sum: u64 = 0;
    let mut q20_count: u64 = 0;
    let mut q30_count: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| CyaneaError::Parse(e.to_string()))?;
        let seq = record.seq();

        sequence_count += 1;
        total_bases += seq.len() as u64;

        for &base in seq.iter() {
            match base.to_ascii_uppercase() {
                b'G' | b'C' => gc_count += 1,
                _ => {}
            }
        }

        if let Some(qual) = record.qual() {
            for &q in qual {
                let phred = q.saturating_sub(33);
                quality_sum += phred as u64;
                if phred >= 20 {
                    q20_count += 1;
                }
                if phred >= 30 {
                    q30_count += 1;
                }
            }
        }
    }

    let gc_content = if total_bases > 0 {
        gc_count as f64 / total_bases as f64
    } else {
        0.0
    };

    let avg_length = if sequence_count > 0 {
        total_bases as f64 / sequence_count as f64
    } else {
        0.0
    };

    let mean_quality = if total_bases > 0 {
        quality_sum as f64 / total_bases as f64
    } else {
        0.0
    };

    let q20_fraction = if total_bases > 0 {
        q20_count as f64 / total_bases as f64
    } else {
        0.0
    };

    let q30_fraction = if total_bases > 0 {
        q30_count as f64 / total_bases as f64
    } else {
        0.0
    };

    Ok(FastqStats {
        sequence_count,
        total_bases,
        gc_content,
        avg_length,
        mean_quality,
        q20_fraction,
        q30_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn length_mismatch_error() {
        let seq = DnaSequence::new(b"ACGT").unwrap();
        let qual = QualityScores::from_raw(vec![30, 30, 30]); // length 3 vs 4
        let result = FastqRecord::new("test".into(), None, seq, qual);
        assert!(result.is_err());
    }

    #[test]
    fn annotated_trait() {
        let seq = DnaSequence::new(b"ACGT").unwrap();
        let qual = QualityScores::from_raw(vec![30, 30, 30, 30]);
        let record =
            FastqRecord::new("read1".into(), Some("sample A".into()), seq, qual).unwrap();
        assert_eq!(record.name(), "read1");
        assert_eq!(record.description(), Some("sample A"));
    }

    #[test]
    fn parse_temp_fastq() {
        let mut file = NamedTempFile::new().unwrap();
        // Write a minimal FASTQ record
        writeln!(file, "@read1 test read").unwrap();
        writeln!(file, "ACGTACGT").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "IIIIIIII").unwrap(); // Q40 in Phred+33
        file.flush().unwrap();

        let records = parse_fastq_file(file.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name(), "read1");
        assert_eq!(records[0].description(), Some("test read"));
        assert_eq!(records[0].sequence().as_bytes(), b"ACGTACGT");
        assert_eq!(records[0].quality().as_slice(), &[40; 8]);
    }
}
