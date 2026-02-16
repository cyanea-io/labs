//! Paired-end FASTQ support.
//!
//! Types and functions for working with paired-end Illumina reads (R1/R2).
//! Supports separate files, interleaved files, and conversions between them.
//!
//! # Example
//!
//! ```no_run
//! use cyanea_seq::paired::{parse_paired_fastq_files, MateValidation};
//!
//! let pairs = parse_paired_fastq_files("reads_R1.fq", "reads_R2.fq", MateValidation::Relaxed)
//!     .unwrap();
//! for pair in &pairs {
//!     println!("{}: R1={} bp, R2={} bp",
//!         pair.pair_name(),
//!         pair.r1().sequence().len(),
//!         pair.r2().sequence().len());
//! }
//! ```

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use cyanea_core::{Annotated, CyaneaError, Result, Sequence};

use crate::fastq::{parse_fastq_stats, FastqRecord, FastqStats};
use crate::quality::{PhredEncoding, QualityScores};
use crate::types::DnaSequence;

// ---------------------------------------------------------------------------
// Core types
// ---------------------------------------------------------------------------

/// How to validate that two reads form a proper mate pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MateValidation {
    /// Names must match after stripping suffixes, and reads must have `/1` and `/2` suffixes.
    Strict,
    /// Names must match after stripping any read-number suffix (`/1`, `/2`, `_1`, `_2`).
    Relaxed,
    /// No validation — trust file order.
    None,
}

/// A paired-end FASTQ record containing R1 (forward) and R2 (reverse) reads.
#[derive(Debug, Clone)]
pub struct PairedFastqRecord {
    r1: FastqRecord,
    r2: FastqRecord,
}

impl PairedFastqRecord {
    /// Create a new paired record with mate validation.
    pub fn new(r1: FastqRecord, r2: FastqRecord, validation: MateValidation) -> Result<Self> {
        match validation {
            MateValidation::Strict => validate_mate_pair_strict(&r1, &r2)?,
            MateValidation::Relaxed => validate_mate_pair(&r1, &r2)?,
            MateValidation::None => {}
        }
        Ok(Self { r1, r2 })
    }

    /// Create a new paired record without any validation.
    pub fn new_unchecked(r1: FastqRecord, r2: FastqRecord) -> Self {
        Self { r1, r2 }
    }

    /// The forward (R1) read.
    pub fn r1(&self) -> &FastqRecord {
        &self.r1
    }

    /// The reverse (R2) read.
    pub fn r2(&self) -> &FastqRecord {
        &self.r2
    }

    /// Consume the pair and return both reads.
    pub fn into_reads(self) -> (FastqRecord, FastqRecord) {
        (self.r1, self.r2)
    }

    /// The shared pair name (R1 name with read-number suffix stripped).
    pub fn pair_name(&self) -> &str {
        strip_read_suffix(self.r1.name())
    }
}

/// Summary statistics for a paired-end FASTQ dataset.
#[derive(Debug, Clone)]
pub struct PairedFastqStats {
    /// Number of read pairs.
    pub pair_count: u64,
    /// Statistics for R1 reads.
    pub r1_stats: FastqStats,
    /// Statistics for R2 reads.
    pub r2_stats: FastqStats,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Strip read-number suffixes from a read name.
///
/// Recognizes `/1`, `/2`, `_1`, `_2` at the end of the name.
/// Illumina-style read numbers (e.g. `1:N:0:ATCACG`) are stored in the
/// description field by the parser, so the name is already clean.
pub fn strip_read_suffix(name: &str) -> &str {
    if name.len() >= 2 {
        let suffix = &name[name.len() - 2..];
        if suffix == "/1" || suffix == "/2" || suffix == "_1" || suffix == "_2" {
            return &name[..name.len() - 2];
        }
    }
    name
}

/// Validate that two reads form a mate pair (relaxed: name prefixes must match).
pub fn validate_mate_pair(r1: &FastqRecord, r2: &FastqRecord) -> Result<()> {
    let name1 = strip_read_suffix(r1.name());
    let name2 = strip_read_suffix(r2.name());
    if name1 != name2 {
        return Err(CyaneaError::InvalidInput(format!(
            "mate pair name mismatch: '{}' vs '{}'",
            r1.name(),
            r2.name()
        )));
    }
    Ok(())
}

/// Validate that two reads form a mate pair (strict: requires `/1` and `/2` suffixes).
pub fn validate_mate_pair_strict(r1: &FastqRecord, r2: &FastqRecord) -> Result<()> {
    if !r1.name().ends_with("/1") {
        return Err(CyaneaError::InvalidInput(format!(
            "R1 read name '{}' does not end with /1",
            r1.name()
        )));
    }
    if !r2.name().ends_with("/2") {
        return Err(CyaneaError::InvalidInput(format!(
            "R2 read name '{}' does not end with /2",
            r2.name()
        )));
    }
    validate_mate_pair(r1, r2)
}

/// Convert raw needletail record fields into a [`FastqRecord`].
fn record_from_parts(id: &[u8], seq: &[u8], qual: Option<&[u8]>) -> Result<FastqRecord> {
    let raw_id =
        std::str::from_utf8(id).map_err(|e| CyaneaError::Parse(e.to_string()))?;
    let (name, description) = match raw_id.split_once(char::is_whitespace) {
        Some((n, d)) => (n.to_string(), Some(d.to_string())),
        None => (raw_id.to_string(), None),
    };
    let sequence = DnaSequence::new(seq)?;
    let qual_bytes =
        qual.ok_or_else(|| CyaneaError::Parse("missing quality scores".into()))?;
    let quality = QualityScores::from_ascii(qual_bytes, PhredEncoding::Phred33)?;
    FastqRecord::new(name, description, sequence, quality)
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

/// Parse paired FASTQ files (separate R1 and R2 files) into paired records.
///
/// Opens both files with needletail and iterates in lockstep. Returns an error
/// if one file has more records than the other or if validation fails.
pub fn parse_paired_fastq_files(
    r1_path: impl AsRef<Path>,
    r2_path: impl AsRef<Path>,
    validation: MateValidation,
) -> Result<Vec<PairedFastqRecord>> {
    let r1_empty = std::fs::metadata(r1_path.as_ref())
        .map(|m| m.len() == 0)
        .unwrap_or(false);
    let r2_empty = std::fs::metadata(r2_path.as_ref())
        .map(|m| m.len() == 0)
        .unwrap_or(false);
    if r1_empty && r2_empty {
        return Ok(Vec::new());
    }
    if r1_empty != r2_empty {
        return Err(CyaneaError::InvalidInput(
            "one file is empty but the other is not".into(),
        ));
    }

    let mut r1_reader = needletail::parse_fastx_file(r1_path.as_ref())
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;
    let mut r2_reader = needletail::parse_fastx_file(r2_path.as_ref())
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;

    let mut pairs = Vec::new();
    loop {
        let r1_next = r1_reader.next();
        let r2_next = r2_reader.next();

        match (r1_next, r2_next) {
            (Some(r1), Some(r2)) => {
                let r1 = r1.map_err(|e| CyaneaError::Parse(e.to_string()))?;
                let r2 = r2.map_err(|e| CyaneaError::Parse(e.to_string()))?;
                let r1_rec = record_from_parts(r1.id(), &r1.seq(), r1.qual())?;
                let r2_rec = record_from_parts(r2.id(), &r2.seq(), r2.qual())?;
                pairs.push(PairedFastqRecord::new(r1_rec, r2_rec, validation)?);
            }
            (None, None) => break,
            (Some(_), None) => {
                return Err(CyaneaError::InvalidInput(
                    "R1 file has more records than R2 file".into(),
                ));
            }
            (None, Some(_)) => {
                return Err(CyaneaError::InvalidInput(
                    "R2 file has more records than R1 file".into(),
                ));
            }
        }
    }

    Ok(pairs)
}

/// Parse an interleaved FASTQ file (alternating R1/R2 records) into paired records.
///
/// Reads two records at a time from a single file. Returns an error on odd
/// record count or if validation fails.
pub fn parse_interleaved_fastq(
    path: impl AsRef<Path>,
    validation: MateValidation,
) -> Result<Vec<PairedFastqRecord>> {
    if std::fs::metadata(path.as_ref())
        .map(|m| m.len() == 0)
        .unwrap_or(false)
    {
        return Ok(Vec::new());
    }

    let mut reader = needletail::parse_fastx_file(path.as_ref())
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;

    let mut pairs = Vec::new();
    loop {
        let r1 = match reader.next() {
            Some(r1) => r1.map_err(|e| CyaneaError::Parse(e.to_string()))?,
            None => break,
        };
        let r1_rec = record_from_parts(r1.id(), &r1.seq(), r1.qual())?;

        let r2 = reader
            .next()
            .ok_or_else(|| {
                CyaneaError::InvalidInput(
                    "odd number of records in interleaved FASTQ".into(),
                )
            })?
            .map_err(|e| CyaneaError::Parse(e.to_string()))?;
        let r2_rec = record_from_parts(r2.id(), &r2.seq(), r2.qual())?;

        pairs.push(PairedFastqRecord::new(r1_rec, r2_rec, validation)?);
    }

    Ok(pairs)
}

/// Compute summary statistics for paired FASTQ files without storing records.
pub fn parse_paired_fastq_stats(
    r1_path: impl AsRef<Path>,
    r2_path: impl AsRef<Path>,
) -> Result<PairedFastqStats> {
    let r1_stats = parse_fastq_stats(r1_path)?;
    let r2_stats = parse_fastq_stats(r2_path)?;
    let pair_count = r1_stats.sequence_count.min(r2_stats.sequence_count);
    Ok(PairedFastqStats {
        pair_count,
        r1_stats,
        r2_stats,
    })
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// Write a single FASTQ record to a writer.
fn write_fastq_record_to(
    writer: &mut impl Write,
    record: &FastqRecord,
    encoding: PhredEncoding,
) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(record.name().as_bytes())?;
    if let Some(desc) = record.description() {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    writer.write_all(b"\n")?;
    writer.write_all(record.sequence().as_bytes())?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&record.quality().to_ascii(encoding))?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// Write paired records to separate R1 and R2 FASTQ files.
pub fn write_paired_fastq(
    pairs: &[PairedFastqRecord],
    r1_path: impl AsRef<Path>,
    r2_path: impl AsRef<Path>,
    encoding: PhredEncoding,
) -> Result<()> {
    let mut r1_writer = BufWriter::new(File::create(r1_path.as_ref())?);
    let mut r2_writer = BufWriter::new(File::create(r2_path.as_ref())?);

    for pair in pairs {
        write_fastq_record_to(&mut r1_writer, &pair.r1, encoding)?;
        write_fastq_record_to(&mut r2_writer, &pair.r2, encoding)?;
    }

    r1_writer.flush()?;
    r2_writer.flush()?;
    Ok(())
}

/// Write paired records to a single interleaved FASTQ file.
pub fn write_interleaved_fastq(
    pairs: &[PairedFastqRecord],
    path: impl AsRef<Path>,
    encoding: PhredEncoding,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path.as_ref())?);

    for pair in pairs {
        write_fastq_record_to(&mut writer, &pair.r1, encoding)?;
        write_fastq_record_to(&mut writer, &pair.r2, encoding)?;
    }

    writer.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Interleave / Deinterleave (streaming)
// ---------------------------------------------------------------------------

/// Interleave two separate FASTQ files into one output file.
///
/// Streaming — does not load all records into memory. Returns the number of
/// pairs written.
pub fn interleave_fastq_files(
    r1_path: impl AsRef<Path>,
    r2_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    validation: MateValidation,
) -> Result<u64> {
    let mut r1_reader = needletail::parse_fastx_file(r1_path.as_ref())
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;
    let mut r2_reader = needletail::parse_fastx_file(r2_path.as_ref())
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;
    let mut writer = BufWriter::new(File::create(output_path.as_ref())?);

    let encoding = PhredEncoding::Phred33;
    let mut count = 0u64;

    loop {
        let r1_next = r1_reader.next();
        let r2_next = r2_reader.next();

        match (r1_next, r2_next) {
            (Some(r1), Some(r2)) => {
                let r1 = r1.map_err(|e| CyaneaError::Parse(e.to_string()))?;
                let r2 = r2.map_err(|e| CyaneaError::Parse(e.to_string()))?;
                let r1_rec = record_from_parts(r1.id(), &r1.seq(), r1.qual())?;
                let r2_rec = record_from_parts(r2.id(), &r2.seq(), r2.qual())?;

                match validation {
                    MateValidation::Strict => validate_mate_pair_strict(&r1_rec, &r2_rec)?,
                    MateValidation::Relaxed => validate_mate_pair(&r1_rec, &r2_rec)?,
                    MateValidation::None => {}
                }

                write_fastq_record_to(&mut writer, &r1_rec, encoding)?;
                write_fastq_record_to(&mut writer, &r2_rec, encoding)?;
                count += 1;
            }
            (None, None) => break,
            (Some(_), None) => {
                return Err(CyaneaError::InvalidInput(
                    "R1 file has more records than R2 file".into(),
                ));
            }
            (None, Some(_)) => {
                return Err(CyaneaError::InvalidInput(
                    "R2 file has more records than R1 file".into(),
                ));
            }
        }
    }

    writer.flush()?;
    Ok(count)
}

/// Deinterleave an interleaved FASTQ file into separate R1 and R2 files.
///
/// Streaming — does not load all records into memory. Returns the number of
/// pairs written.
pub fn deinterleave_fastq_file(
    input_path: impl AsRef<Path>,
    r1_path: impl AsRef<Path>,
    r2_path: impl AsRef<Path>,
    validation: MateValidation,
) -> Result<u64> {
    let mut reader = needletail::parse_fastx_file(input_path.as_ref())
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;
    let mut r1_writer = BufWriter::new(File::create(r1_path.as_ref())?);
    let mut r2_writer = BufWriter::new(File::create(r2_path.as_ref())?);

    let encoding = PhredEncoding::Phred33;
    let mut count = 0u64;

    loop {
        let r1 = match reader.next() {
            Some(r1) => r1.map_err(|e| CyaneaError::Parse(e.to_string()))?,
            None => break,
        };
        let r1_rec = record_from_parts(r1.id(), &r1.seq(), r1.qual())?;

        let r2 = reader
            .next()
            .ok_or_else(|| {
                CyaneaError::InvalidInput(
                    "odd number of records in interleaved FASTQ".into(),
                )
            })?
            .map_err(|e| CyaneaError::Parse(e.to_string()))?;
        let r2_rec = record_from_parts(r2.id(), &r2.seq(), r2.qual())?;

        match validation {
            MateValidation::Strict => validate_mate_pair_strict(&r1_rec, &r2_rec)?,
            MateValidation::Relaxed => validate_mate_pair(&r1_rec, &r2_rec)?,
            MateValidation::None => {}
        }

        write_fastq_record_to(&mut r1_writer, &r1_rec, encoding)?;
        write_fastq_record_to(&mut r2_writer, &r2_rec, encoding)?;
        count += 1;
    }

    r1_writer.flush()?;
    r2_writer.flush()?;
    Ok(count)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_record(name: &str, seq: &[u8], quals: &[u8]) -> FastqRecord {
        let sequence = DnaSequence::new(seq).unwrap();
        let quality = QualityScores::from_raw(quals.to_vec());
        FastqRecord::new(name.to_string(), None, sequence, quality).unwrap()
    }

    fn make_record_with_desc(
        name: &str,
        desc: &str,
        seq: &[u8],
        quals: &[u8],
    ) -> FastqRecord {
        let sequence = DnaSequence::new(seq).unwrap();
        let quality = QualityScores::from_raw(quals.to_vec());
        FastqRecord::new(name.to_string(), Some(desc.to_string()), sequence, quality).unwrap()
    }

    fn write_temp_fastq(records: &[(&str, &str, &str)]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        for (name, seq, qual) in records {
            writeln!(file, "@{}", name).unwrap();
            writeln!(file, "{}", seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", qual).unwrap();
        }
        file.flush().unwrap();
        file
    }

    // --- strip_read_suffix ---

    #[test]
    fn strip_suffix_slash_1() {
        assert_eq!(strip_read_suffix("read1/1"), "read1");
    }

    #[test]
    fn strip_suffix_slash_2() {
        assert_eq!(strip_read_suffix("read1/2"), "read1");
    }

    #[test]
    fn strip_suffix_underscore_1() {
        assert_eq!(strip_read_suffix("read1_1"), "read1");
    }

    #[test]
    fn strip_suffix_underscore_2() {
        assert_eq!(strip_read_suffix("read1_2"), "read1");
    }

    #[test]
    fn strip_suffix_none() {
        assert_eq!(strip_read_suffix("read1"), "read1");
    }

    #[test]
    fn strip_suffix_illumina_name() {
        // Illumina names don't have suffixes in the name field
        assert_eq!(
            strip_read_suffix("INSTRUMENT:1:2:3:4:5:6"),
            "INSTRUMENT:1:2:3:4:5:6"
        );
    }

    #[test]
    fn strip_suffix_short_name() {
        assert_eq!(strip_read_suffix("r"), "r");
    }

    // --- Mate validation ---

    #[test]
    fn validate_relaxed_matching_slash() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1/2", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair(&r1, &r2).is_ok());
    }

    #[test]
    fn validate_relaxed_matching_underscore() {
        let r1 = make_record("read1_1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1_2", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair(&r1, &r2).is_ok());
    }

    #[test]
    fn validate_relaxed_identical_names() {
        let r1 = make_record("read1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair(&r1, &r2).is_ok());
    }

    #[test]
    fn validate_relaxed_mismatch() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read2/2", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair(&r1, &r2).is_err());
    }

    #[test]
    fn validate_strict_valid() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1/2", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair_strict(&r1, &r2).is_ok());
    }

    #[test]
    fn validate_strict_missing_suffix_r1() {
        let r1 = make_record("read1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1/2", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair_strict(&r1, &r2).is_err());
    }

    #[test]
    fn validate_strict_missing_suffix_r2() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair_strict(&r1, &r2).is_err());
    }

    #[test]
    fn validate_strict_wrong_suffix_order() {
        let r1 = make_record("read1/2", b"ACGT", &[30; 4]);
        let r2 = make_record("read1/1", b"ACGT", &[30; 4]);
        assert!(validate_mate_pair_strict(&r1, &r2).is_err());
    }

    // --- PairedFastqRecord ---

    #[test]
    fn pair_creation_valid() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1/2", b"TGCA", &[25; 4]);
        let pair = PairedFastqRecord::new(r1, r2, MateValidation::Relaxed).unwrap();
        assert_eq!(pair.r1().sequence().as_bytes(), b"ACGT");
        assert_eq!(pair.r2().sequence().as_bytes(), b"TGCA");
        assert_eq!(pair.pair_name(), "read1");
    }

    #[test]
    fn pair_creation_invalid() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read2/2", b"TGCA", &[25; 4]);
        assert!(PairedFastqRecord::new(r1, r2, MateValidation::Relaxed).is_err());
    }

    #[test]
    fn pair_creation_unchecked() {
        let r1 = make_record("read1", b"ACGT", &[30; 4]);
        let r2 = make_record("read2", b"TGCA", &[25; 4]);
        let pair = PairedFastqRecord::new_unchecked(r1, r2);
        assert_eq!(pair.r1().name(), "read1");
        assert_eq!(pair.r2().name(), "read2");
    }

    #[test]
    fn pair_into_reads() {
        let r1 = make_record("read1/1", b"ACGT", &[30; 4]);
        let r2 = make_record("read1/2", b"TGCA", &[25; 4]);
        let pair = PairedFastqRecord::new_unchecked(r1, r2);
        let (r1, r2) = pair.into_reads();
        assert_eq!(r1.sequence().as_bytes(), b"ACGT");
        assert_eq!(r2.sequence().as_bytes(), b"TGCA");
    }

    #[test]
    fn pair_no_validation() {
        let r1 = make_record("foo", b"ACGT", &[30; 4]);
        let r2 = make_record("bar", b"TGCA", &[25; 4]);
        assert!(PairedFastqRecord::new(r1, r2, MateValidation::None).is_ok());
    }

    // --- Parse paired files ---

    #[test]
    fn parse_paired_files_matching() {
        let r1_file = write_temp_fastq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TGCATGCA", "IIIIIIII"),
        ]);
        let r2_file = write_temp_fastq(&[
            ("read1/2", "GGGGCCCC", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "IIIIIIII"),
        ]);
        let pairs = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::Relaxed,
        )
        .unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].r1().sequence().as_bytes(), b"ACGTACGT");
        assert_eq!(pairs[0].r2().sequence().as_bytes(), b"GGGGCCCC");
        assert_eq!(pairs[1].pair_name(), "read2");
    }

    #[test]
    fn parse_paired_files_unequal_r1_longer() {
        let r1_file = write_temp_fastq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TGCATGCA", "IIIIIIII"),
        ]);
        let r2_file = write_temp_fastq(&[("read1/2", "GGGGCCCC", "IIIIIIII")]);
        let result = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::Relaxed,
        );
        assert!(result.is_err());
    }

    #[test]
    fn parse_paired_files_unequal_r2_longer() {
        let r1_file = write_temp_fastq(&[("read1/1", "ACGTACGT", "IIIIIIII")]);
        let r2_file = write_temp_fastq(&[
            ("read1/2", "GGGGCCCC", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "IIIIIIII"),
        ]);
        let result = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::Relaxed,
        );
        assert!(result.is_err());
    }

    #[test]
    fn parse_paired_files_no_validation() {
        let r1_file = write_temp_fastq(&[("foo", "ACGTACGT", "IIIIIIII")]);
        let r2_file = write_temp_fastq(&[("bar", "GGGGCCCC", "IIIIIIII")]);
        let pairs = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::None,
        )
        .unwrap();
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn parse_paired_files_validation_failure() {
        let r1_file = write_temp_fastq(&[("read1/1", "ACGTACGT", "IIIIIIII")]);
        let r2_file = write_temp_fastq(&[("read2/2", "GGGGCCCC", "IIIIIIII")]);
        let result = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::Relaxed,
        );
        assert!(result.is_err());
    }

    #[test]
    fn parse_paired_files_strict_validation() {
        let r1_file = write_temp_fastq(&[("read1/1", "ACGT", "IIII")]);
        let r2_file = write_temp_fastq(&[("read1/2", "TGCA", "IIII")]);
        let pairs = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::Strict,
        )
        .unwrap();
        assert_eq!(pairs.len(), 1);
    }

    // --- Parse interleaved ---

    #[test]
    fn parse_interleaved_valid() {
        let file = write_temp_fastq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read1/2", "TGCATGCA", "IIIIIIII"),
            ("read2/1", "GGGGCCCC", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "IIIIIIII"),
        ]);
        let pairs = parse_interleaved_fastq(file.path(), MateValidation::Relaxed).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].r1().sequence().as_bytes(), b"ACGTACGT");
        assert_eq!(pairs[0].r2().sequence().as_bytes(), b"TGCATGCA");
    }

    #[test]
    fn parse_interleaved_odd_count() {
        let file = write_temp_fastq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read1/2", "TGCATGCA", "IIIIIIII"),
            ("read2/1", "GGGGCCCC", "IIIIIIII"),
        ]);
        let result = parse_interleaved_fastq(file.path(), MateValidation::None);
        assert!(result.is_err());
    }

    // --- Write + re-parse round-trips ---

    #[test]
    fn write_parse_roundtrip_separate() {
        let r1 = make_record("read1/1", b"ACGTACGT", &[30; 8]);
        let r2 = make_record("read1/2", b"TGCATGCA", &[25; 8]);
        let pairs = vec![PairedFastqRecord::new_unchecked(r1, r2)];

        let r1_file = NamedTempFile::new().unwrap();
        let r2_file = NamedTempFile::new().unwrap();
        write_paired_fastq(&pairs, r1_file.path(), r2_file.path(), PhredEncoding::Phred33)
            .unwrap();

        let parsed = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::None,
        )
        .unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].r1().sequence().as_bytes(), b"ACGTACGT");
        assert_eq!(parsed[0].r2().sequence().as_bytes(), b"TGCATGCA");
        assert_eq!(parsed[0].r1().quality().as_slice(), &[30; 8]);
        assert_eq!(parsed[0].r2().quality().as_slice(), &[25; 8]);
    }

    #[test]
    fn write_parse_roundtrip_interleaved() {
        let r1 = make_record("read1/1", b"ACGTACGT", &[30; 8]);
        let r2 = make_record("read1/2", b"TGCATGCA", &[25; 8]);
        let pairs = vec![PairedFastqRecord::new_unchecked(r1, r2)];

        let file = NamedTempFile::new().unwrap();
        write_interleaved_fastq(&pairs, file.path(), PhredEncoding::Phred33).unwrap();

        let parsed = parse_interleaved_fastq(file.path(), MateValidation::None).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].r1().sequence().as_bytes(), b"ACGTACGT");
        assert_eq!(parsed[0].r2().sequence().as_bytes(), b"TGCATGCA");
    }

    #[test]
    fn write_parse_roundtrip_with_description() {
        let r1 = make_record_with_desc("read1", "1:N:0:ATCACG", b"ACGTACGT", &[30; 8]);
        let r2 = make_record_with_desc("read1", "2:N:0:ATCACG", b"TGCATGCA", &[25; 8]);
        let pairs = vec![PairedFastqRecord::new_unchecked(r1, r2)];

        let file = NamedTempFile::new().unwrap();
        write_interleaved_fastq(&pairs, file.path(), PhredEncoding::Phred33).unwrap();

        let parsed = parse_interleaved_fastq(file.path(), MateValidation::None).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].r1().name(), "read1");
        assert_eq!(parsed[0].r1().description(), Some("1:N:0:ATCACG"));
        assert_eq!(parsed[0].r2().description(), Some("2:N:0:ATCACG"));
    }

    #[test]
    fn write_parse_roundtrip_multiple_pairs() {
        let pairs = vec![
            PairedFastqRecord::new_unchecked(
                make_record("r1/1", b"AAAA", &[30; 4]),
                make_record("r1/2", b"CCCC", &[25; 4]),
            ),
            PairedFastqRecord::new_unchecked(
                make_record("r2/1", b"GGGG", &[35; 4]),
                make_record("r2/2", b"TTTT", &[20; 4]),
            ),
        ];

        let r1_file = NamedTempFile::new().unwrap();
        let r2_file = NamedTempFile::new().unwrap();
        write_paired_fastq(&pairs, r1_file.path(), r2_file.path(), PhredEncoding::Phred33)
            .unwrap();

        let parsed = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::None,
        )
        .unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].r1().sequence().as_bytes(), b"AAAA");
        assert_eq!(parsed[1].r2().sequence().as_bytes(), b"TTTT");
    }

    // --- Interleave / deinterleave ---

    #[test]
    fn interleave_deinterleave_roundtrip() {
        let r1_file = write_temp_fastq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TGCATGCA", "IIIIIIII"),
        ]);
        let r2_file = write_temp_fastq(&[
            ("read1/2", "GGGGCCCC", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "IIIIIIII"),
        ]);

        let interleaved = NamedTempFile::new().unwrap();
        let count = interleave_fastq_files(
            r1_file.path(),
            r2_file.path(),
            interleaved.path(),
            MateValidation::Relaxed,
        )
        .unwrap();
        assert_eq!(count, 2);

        let out_r1 = NamedTempFile::new().unwrap();
        let out_r2 = NamedTempFile::new().unwrap();
        let count = deinterleave_fastq_file(
            interleaved.path(),
            out_r1.path(),
            out_r2.path(),
            MateValidation::Relaxed,
        )
        .unwrap();
        assert_eq!(count, 2);

        let pairs = parse_paired_fastq_files(
            out_r1.path(),
            out_r2.path(),
            MateValidation::None,
        )
        .unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].r1().sequence().as_bytes(), b"ACGTACGT");
        assert_eq!(pairs[0].r2().sequence().as_bytes(), b"GGGGCCCC");
    }

    #[test]
    fn interleave_unequal_files() {
        let r1_file = write_temp_fastq(&[
            ("read1/1", "ACGT", "IIII"),
            ("read2/1", "TGCA", "IIII"),
        ]);
        let r2_file = write_temp_fastq(&[("read1/2", "GGGG", "IIII")]);
        let output = NamedTempFile::new().unwrap();
        let result = interleave_fastq_files(
            r1_file.path(),
            r2_file.path(),
            output.path(),
            MateValidation::None,
        );
        assert!(result.is_err());
    }

    // --- Empty / single pair edge cases ---

    #[test]
    fn parse_paired_empty_files() {
        let r1_file = NamedTempFile::new().unwrap();
        let r2_file = NamedTempFile::new().unwrap();
        let pairs = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::None,
        )
        .unwrap();
        assert!(pairs.is_empty());
    }

    #[test]
    fn parse_interleaved_empty_file() {
        let file = NamedTempFile::new().unwrap();
        let pairs = parse_interleaved_fastq(file.path(), MateValidation::None).unwrap();
        assert!(pairs.is_empty());
    }

    #[test]
    fn parse_paired_single_pair() {
        let r1_file = write_temp_fastq(&[("read1/1", "ACGT", "IIII")]);
        let r2_file = write_temp_fastq(&[("read1/2", "TGCA", "IIII")]);
        let pairs = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::Strict,
        )
        .unwrap();
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn write_empty_pairs() {
        let pairs: Vec<PairedFastqRecord> = vec![];
        let r1_file = NamedTempFile::new().unwrap();
        let r2_file = NamedTempFile::new().unwrap();
        write_paired_fastq(&pairs, r1_file.path(), r2_file.path(), PhredEncoding::Phred33)
            .unwrap();
        let parsed = parse_paired_fastq_files(
            r1_file.path(),
            r2_file.path(),
            MateValidation::None,
        )
        .unwrap();
        assert!(parsed.is_empty());
    }

    // --- Stats ---

    #[test]
    fn paired_stats() {
        let r1_file = write_temp_fastq(&[
            ("read1", "ACGTACGT", "IIIIIIII"),
            ("read2", "TGCA", "IIII"),
        ]);
        let r2_file = write_temp_fastq(&[
            ("read1", "GGGGCCCC", "IIIIIIII"),
            ("read2", "AAAA", "IIII"),
        ]);
        let stats = parse_paired_fastq_stats(r1_file.path(), r2_file.path()).unwrap();
        assert_eq!(stats.pair_count, 2);
        assert_eq!(stats.r1_stats.sequence_count, 2);
        assert_eq!(stats.r2_stats.sequence_count, 2);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use crate::types::DnaSequence;
    use proptest::prelude::*;
    use tempfile::NamedTempFile;

    fn dna_and_quality(max_len: usize) -> impl Strategy<Value = (Vec<u8>, Vec<u8>)> {
        (1..=max_len).prop_flat_map(|len| {
            let seq = proptest::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
                len,
            );
            let qual = proptest::collection::vec(0..=41u8, len);
            (seq, qual)
        })
    }

    proptest! {
        #[test]
        fn roundtrip_write_parse(
            (seq1, qual1) in dna_and_quality(100),
            (seq2, qual2) in dna_and_quality(100),
        ) {
            let r1 = {
                let s = DnaSequence::new(&seq1).unwrap();
                let q = QualityScores::from_raw(qual1.clone());
                FastqRecord::new("read/1".into(), None, s, q).unwrap()
            };
            let r2 = {
                let s = DnaSequence::new(&seq2).unwrap();
                let q = QualityScores::from_raw(qual2.clone());
                FastqRecord::new("read/2".into(), None, s, q).unwrap()
            };
            let pairs = vec![PairedFastqRecord::new_unchecked(r1, r2)];

            let file = NamedTempFile::new().unwrap();
            write_interleaved_fastq(&pairs, file.path(), PhredEncoding::Phred33).unwrap();

            let parsed = parse_interleaved_fastq(file.path(), MateValidation::None).unwrap();
            prop_assert_eq!(parsed.len(), 1);
            prop_assert_eq!(parsed[0].r1().sequence().as_bytes(), seq1.as_slice());
            prop_assert_eq!(parsed[0].r2().sequence().as_bytes(), seq2.as_slice());
            prop_assert_eq!(parsed[0].r1().quality().as_slice(), qual1.as_slice());
            prop_assert_eq!(parsed[0].r2().quality().as_slice(), qual2.as_slice());
        }

        #[test]
        fn interleave_deinterleave_identity(
            (seq1, qual1) in dna_and_quality(50),
            (seq2, qual2) in dna_and_quality(50),
        ) {
            use std::io::Write as _;

            let r1_file = NamedTempFile::new().unwrap();
            let r2_file = NamedTempFile::new().unwrap();

            // Write R1 FASTQ
            {
                let ascii1: Vec<u8> = qual1.iter().map(|&q| q + 33).collect();
                let seq_str = std::str::from_utf8(&seq1).unwrap();
                let qual_str = std::str::from_utf8(&ascii1).unwrap();
                write!(r1_file.as_file(), "@read/1\n{}\n+\n{}\n", seq_str, qual_str).unwrap();
            }
            // Write R2 FASTQ
            {
                let ascii2: Vec<u8> = qual2.iter().map(|&q| q + 33).collect();
                let seq_str = std::str::from_utf8(&seq2).unwrap();
                let qual_str = std::str::from_utf8(&ascii2).unwrap();
                write!(r2_file.as_file(), "@read/2\n{}\n+\n{}\n", seq_str, qual_str).unwrap();
            }

            let interleaved = NamedTempFile::new().unwrap();
            let count = interleave_fastq_files(
                r1_file.path(),
                r2_file.path(),
                interleaved.path(),
                MateValidation::None,
            )
            .unwrap();
            prop_assert_eq!(count, 1);

            let out_r1 = NamedTempFile::new().unwrap();
            let out_r2 = NamedTempFile::new().unwrap();
            let count2 = deinterleave_fastq_file(
                interleaved.path(),
                out_r1.path(),
                out_r2.path(),
                MateValidation::None,
            )
            .unwrap();
            prop_assert_eq!(count2, 1);

            let pairs = parse_paired_fastq_files(
                out_r1.path(),
                out_r2.path(),
                MateValidation::None,
            )
            .unwrap();
            prop_assert_eq!(pairs.len(), 1);
            prop_assert_eq!(pairs[0].r1().sequence().as_bytes(), seq1.as_slice());
            prop_assert_eq!(pairs[0].r2().sequence().as_bytes(), seq2.as_slice());
        }
    }
}
