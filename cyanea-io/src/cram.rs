//! CRAM (Compressed Reference-oriented Alignment Map) parser.
//!
//! Parses CRAM files into [`SamRecord`] records, reusing the same types
//! as the SAM/BAM modules. CRAM achieves better compression than BAM
//! through reference-based encoding.
//!
//! Requires the `cram` feature flag. The `cram` feature implies `sam`,
//! so [`SamRecord`] and [`SamStats`] are always available.

use std::path::{Path, PathBuf};

use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;

use cyanea_core::{CyaneaError, Result};

use crate::sam::{sam_stats, SamRecord, SamStats};

/// Configuration for CRAM reading.
#[derive(Debug, Clone, Default)]
pub struct CramConfig {
    /// Path to reference FASTA. `None` uses embedded references.
    pub reference_path: Option<PathBuf>,
}

/// Parse a CRAM file and return all alignment records.
///
/// Uses the provided [`CramConfig`] to locate reference sequences.
/// A reference FASTA is required for CRAM files containing mapped reads.
pub fn parse_cram(path: impl AsRef<Path>, config: &CramConfig) -> Result<Vec<SamRecord>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;

    let repository = build_repository(config)?;

    let mut reader = noodles_cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_reader(std::io::BufReader::new(file));

    let header = reader.read_header().map_err(|e| {
        CyaneaError::Parse(format!("{}: failed to read CRAM header: {e}", path.display()))
    })?;

    let mut records = Vec::new();
    for result in reader.records(&header) {
        let cram_record = result.map_err(|e| {
            CyaneaError::Parse(format!(
                "{}: error reading CRAM record: {e}",
                path.display()
            ))
        })?;
        let sam = convert_record(&cram_record, &header)?;
        records.push(sam);
    }

    Ok(records)
}

/// Parse a CRAM file with default configuration (no external reference).
pub fn parse_cram_default(path: impl AsRef<Path>) -> Result<Vec<SamRecord>> {
    parse_cram(path, &CramConfig::default())
}

/// Compute summary statistics from a CRAM file.
pub fn cram_stats(path: impl AsRef<Path>, config: &CramConfig) -> Result<SamStats> {
    let records = parse_cram(path, config)?;
    Ok(sam_stats(&records))
}

/// Compute summary statistics from a CRAM file with default configuration.
pub fn cram_stats_default(path: impl AsRef<Path>) -> Result<SamStats> {
    cram_stats(path, &CramConfig::default())
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Build a FASTA reference repository from the config.
fn build_repository(config: &CramConfig) -> Result<noodles_fasta::Repository> {
    match &config.reference_path {
        Some(ref_path) => {
            let file = std::fs::File::open(ref_path).map_err(|e| {
                CyaneaError::Io(std::io::Error::new(
                    e.kind(),
                    format!("reference FASTA {}: {}", ref_path.display(), e),
                ))
            })?;
            let mut fasta_reader =
                noodles_fasta::io::Reader::new(std::io::BufReader::new(file));
            let records: Vec<noodles_fasta::Record> = fasta_reader
                .records()
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(|e| {
                    CyaneaError::Parse(format!(
                        "reference FASTA {}: {e}",
                        ref_path.display()
                    ))
                })?;
            Ok(noodles_fasta::Repository::new(records))
        }
        None => Ok(noodles_fasta::Repository::default()),
    }
}

fn cigar_kind_char(kind: CigarKind) -> char {
    match kind {
        CigarKind::Match => 'M',
        CigarKind::Insertion => 'I',
        CigarKind::Deletion => 'D',
        CigarKind::Skip => 'N',
        CigarKind::SoftClip => 'S',
        CigarKind::HardClip => 'H',
        CigarKind::Pad => 'P',
        CigarKind::SequenceMatch => '=',
        CigarKind::SequenceMismatch => 'X',
    }
}

/// Convert a noodles CRAM record into our SamRecord.
fn convert_record(
    record: &noodles_cram::Record,
    header: &noodles_sam::Header,
) -> Result<SamRecord> {
    // Read name
    let qname = record
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
        .unwrap_or_else(|| "*".to_string());

    // Flags (returns Flags directly, not Result)
    let flags = record.flags();
    let flag = u16::from(flags);

    // Reference name
    let rname = record
        .reference_sequence(header.reference_sequences())
        .and_then(|r| r.ok())
        .map(|(name, _)| String::from_utf8_lossy(name).to_string())
        .unwrap_or_else(|| "*".to_string());

    // Position (1-based, returns Option<Position> directly)
    let pos = record
        .alignment_start()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(0);

    // Mapping quality (returns Option<MappingQuality> directly)
    let mapq = record
        .mapping_quality()
        .map(|m| u8::from(m))
        .unwrap_or(255);

    // CIGAR string — reconstruct from features
    let cigar = {
        use noodles_sam::alignment::record::Cigar as CigarTrait;
        let cigar_obj = <noodles_cram::Record as noodles_sam::alignment::Record>::cigar(record);
        if cigar_obj.is_empty() {
            "*".to_string()
        } else {
            let mut s = String::new();
            for op_result in cigar_obj.iter() {
                match op_result {
                    Ok(op) => {
                        s.push_str(&format!("{}{}", op.len(), cigar_kind_char(op.kind())));
                    }
                    Err(_) => {
                        return Err(CyaneaError::Parse("failed to read CIGAR op".into()));
                    }
                }
            }
            if s.is_empty() { "*".to_string() } else { s }
        }
    };

    // Sequence — access via record_buf Sequence type
    let sequence = {
        let seq = record.sequence();
        if seq.is_empty() {
            "*".to_string()
        } else {
            // Sequence is a noodles_sam::alignment::record_buf::Sequence
            // which implements AsRef<[u8]>
            let bytes: &[u8] = seq.as_ref();
            String::from_utf8_lossy(bytes).to_string()
        }
    };

    // Quality scores
    let quality = {
        let quals = record.quality_scores();
        let bytes: &[u8] = quals.as_ref();
        if bytes.is_empty() || bytes.iter().all(|&q| q == 0xFF) {
            "*".to_string()
        } else {
            let phred: Vec<u8> = bytes.iter().map(|&q| q + 33).collect();
            String::from_utf8(phred).unwrap_or_else(|_| "*".to_string())
        }
    };

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    use noodles_core::Position;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record_buf::{
        QualityScores as QualBuf, Sequence as SeqBuf,
    };
    use noodles_sam::header::record::value::map::ReferenceSequence;
    use noodles_sam::header::record::value::Map;
    use std::num::NonZeroUsize;

    const REF_LEN: usize = 1000;

    fn build_header() -> noodles_sam::Header {
        let len1 = NonZeroUsize::try_from(REF_LEN).unwrap();
        let len2 = NonZeroUsize::try_from(REF_LEN).unwrap();
        noodles_sam::Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len1))
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(len2))
            .build()
    }

    fn build_repository() -> noodles_fasta::Repository {
        use noodles_fasta::record::{Definition, Sequence};

        let ref_seqs = vec![
            noodles_fasta::Record::new(
                Definition::new("chr1", None),
                Sequence::from(vec![b'A'; REF_LEN]),
            ),
            noodles_fasta::Record::new(
                Definition::new("chr2", None),
                Sequence::from(vec![b'A'; REF_LEN]),
            ),
        ];
        noodles_fasta::Repository::new(ref_seqs)
    }

    /// Write a temp FASTA reference file matching our test header.
    fn write_ref_fasta() -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".fasta").unwrap();
        write!(file, ">chr1\n{}\n", "A".repeat(REF_LEN)).unwrap();
        write!(file, ">chr2\n{}\n", "A".repeat(REF_LEN)).unwrap();
        file.flush().unwrap();
        file
    }

    /// Returns (cram_file, ref_fasta_file).
    fn write_test_cram(recs: &[TestRec]) -> (NamedTempFile, NamedTempFile) {
        let header = build_header();
        let repository = build_repository();

        let mut writer = noodles_cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_writer(Vec::new());

        writer.write_file_definition().unwrap();
        writer.write_file_header(&header).unwrap();

        for rec in recs {
            let cram_rec = rec.to_cram_record();
            writer.write_record(&header, cram_rec).unwrap();
        }
        writer.try_finish(&header).unwrap();

        let bytes = writer.get_ref().clone();
        let mut cram_file = NamedTempFile::with_suffix(".cram").unwrap();
        cram_file.write_all(&bytes).unwrap();
        cram_file.flush().unwrap();

        let ref_file = write_ref_fasta();
        (cram_file, ref_file)
    }

    fn config_with_ref(ref_file: &NamedTempFile) -> CramConfig {
        CramConfig {
            reference_path: Some(ref_file.path().to_path_buf()),
        }
    }

    struct TestRec {
        qname: &'static str,
        flag: u16,
        ref_id: Option<usize>,
        pos: Option<usize>, // 1-based
        mapq: Option<u8>,
        seq: &'static str,
        qual: Vec<u8>, // raw phred (not +33)
    }

    impl TestRec {
        fn mapped(qname: &'static str, ref_id: usize, pos: usize, mapq: u8, seq: &'static str) -> Self {
            Self {
                qname,
                flag: 0,
                ref_id: Some(ref_id),
                pos: Some(pos),
                mapq: Some(mapq),
                seq,
                qual: vec![30; seq.len()],
            }
        }

        fn unmapped(qname: &'static str, seq: &'static str) -> Self {
            Self {
                qname,
                flag: 4,
                ref_id: None,
                pos: None,
                mapq: None,
                seq,
                qual: vec![],
            }
        }

        fn to_cram_record(&self) -> noodles_cram::Record {
            let mut builder = noodles_cram::Record::builder();

            builder = builder.set_bam_flags(Flags::from(self.flag));

            if let Some(ref_id) = self.ref_id {
                builder = builder.set_reference_sequence_id(ref_id);
            }
            if let Some(pos) = self.pos {
                builder = builder.set_alignment_start(Position::try_from(pos).unwrap());
            }
            if let Some(mapq) = self.mapq {
                builder = builder.set_mapping_quality(MappingQuality::new(mapq).unwrap());
            }

            builder = builder.set_read_length(self.seq.len());

            if !self.seq.is_empty() {
                builder = builder.set_bases(SeqBuf::from(self.seq.as_bytes().to_vec()));
            }
            if !self.qual.is_empty() {
                builder = builder.set_quality_scores(QualBuf::from(self.qual.clone()));
            }

            builder = builder.set_name(self.qname);

            builder.build()
        }
    }

    #[test]
    fn roundtrip_with_ref() {
        let recs = vec![TestRec::mapped("read1", 0, 100, 60, "ACGTACGTAC")];
        let (cram_file, ref_file) = write_test_cram(&recs);
        let config = config_with_ref(&ref_file);
        let result = parse_cram(cram_file.path(), &config).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].qname, "read1");
        assert!(result[0].is_mapped());
        assert_eq!(result[0].pos, 100);
        assert_eq!(result[0].mapq, 60);
    }

    #[test]
    fn cram_stats_basic() {
        let recs = vec![
            TestRec::mapped("read1", 0, 100, 60, "ACGT"),
            TestRec::mapped("read2", 0, 200, 30, "TGCA"),
            TestRec::unmapped("read3", "NNNN"),
        ];
        let (cram_file, ref_file) = write_test_cram(&recs);
        let config = config_with_ref(&ref_file);
        let stats = cram_stats(cram_file.path(), &config).unwrap();
        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.mapped, 2);
        assert_eq!(stats.unmapped, 1);
        assert!((stats.avg_mapq - 45.0).abs() < f64::EPSILON);
    }

    #[test]
    fn unmapped_reads() {
        let recs = vec![
            TestRec::unmapped("u1", "AAAA"),
            TestRec::unmapped("u2", "CCCC"),
        ];
        let (cram_file, _ref_file) = write_test_cram(&recs);
        let result = parse_cram_default(cram_file.path()).unwrap();
        assert_eq!(result.len(), 2);
        assert!(result[0].is_unmapped());
        assert!(result[1].is_unmapped());
    }

    #[test]
    fn multiple_references() {
        let recs = vec![
            TestRec::mapped("r1", 0, 100, 60, "ACGT"),
            TestRec::mapped("r2", 1, 200, 40, "TGCA"),
        ];
        let (cram_file, ref_file) = write_test_cram(&recs);
        let config = config_with_ref(&ref_file);
        let result = parse_cram(cram_file.path(), &config).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].rname, "chr1");
        assert_eq!(result[1].rname, "chr2");
    }

    #[test]
    fn empty_cram() {
        let (cram_file, ref_file) = write_test_cram(&[]);
        let config = config_with_ref(&ref_file);
        let result = parse_cram(cram_file.path(), &config).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn quality_scores() {
        let recs = vec![TestRec::mapped("r1", 0, 100, 60, "ACGT")];
        let (cram_file, ref_file) = write_test_cram(&recs);
        let config = config_with_ref(&ref_file);
        let result = parse_cram(cram_file.path(), &config).unwrap();
        // Quality should be phred+33 encoded
        assert_ne!(result[0].quality, "*");
    }

    #[test]
    fn nonexistent_file_error() {
        let result = parse_cram_default("/nonexistent/file.cram");
        assert!(result.is_err());
    }
}
