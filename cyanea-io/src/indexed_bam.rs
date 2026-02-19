//! Indexed BAM reader for random-access region queries.
//!
//! Uses noodles for BAI/CSI index parsing and BGZF seeking.
//! Records are converted to the same [`SamRecord`] type used by the sequential BAM parser.

use std::path::Path;

use cyanea_core::{CyaneaError, Result};

use crate::bam::BamReference;
use crate::sam::SamRecord;

use noodles_bam as noodles_bam_crate;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;

/// An indexed BAM reader for random-access region queries.
pub struct IndexedBamReader {
    reader: noodles_bam_crate::io::IndexedReader<noodles_bgzf::Reader<std::fs::File>>,
    header: noodles_sam::Header,
}

impl IndexedBamReader {
    /// Open a BAM file with its index.
    ///
    /// Looks for the index at `<bam_path>.bai` automatically.
    pub fn open(bam_path: impl AsRef<Path>) -> Result<Self> {
        let bam_path = bam_path.as_ref();
        let reader = noodles_bam_crate::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .map_err(|e| {
                CyaneaError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("{}: {}", bam_path.display(), e),
                ))
            })?;

        Self::from_reader(reader)
    }

    /// Open a BAM file with an explicitly specified index path.
    pub fn open_with_index(
        bam_path: impl AsRef<Path>,
        bai_path: impl AsRef<Path>,
    ) -> Result<Self> {
        let bam_path = bam_path.as_ref();
        let bai_path = bai_path.as_ref();

        let index = noodles_bam_crate::bai::read(bai_path).map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("{}: {}", bai_path.display(), e),
            ))
        })?;

        let file = std::fs::File::open(bam_path).map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {}", bam_path.display(), e),
            ))
        })?;

        let reader = noodles_bam_crate::io::IndexedReader::new(file, index);

        Self::from_reader(reader)
    }

    fn from_reader(
        mut reader: noodles_bam_crate::io::IndexedReader<
            noodles_bgzf::Reader<std::fs::File>,
        >,
    ) -> Result<Self> {
        let header = reader.read_header().map_err(|e| {
            CyaneaError::Parse(format!("failed to read BAM header: {e}"))
        })?;

        Ok(Self { reader, header })
    }

    /// Fetch all records overlapping a genomic region.
    ///
    /// Coordinates are 0-based, half-open `[start, end)`.
    pub fn fetch(&mut self, chrom: &str, start: u64, end: u64) -> Result<Vec<SamRecord>> {
        // Find reference sequence ID
        self.header
            .reference_sequences()
            .get_index_of(&chrom.as_bytes()[..])
            .ok_or_else(|| {
                CyaneaError::InvalidInput(format!(
                    "chromosome '{}' not found in BAM header",
                    chrom
                ))
            })?;

        let start_pos = noodles_core::Position::try_from((start + 1) as usize).map_err(|_| {
            CyaneaError::InvalidInput(format!("invalid start position: {}", start))
        })?;
        let end_pos = noodles_core::Position::try_from(end as usize).map_err(|_| {
            CyaneaError::InvalidInput(format!("invalid end position: {}", end))
        })?;

        let region = noodles_core::Region::new(chrom, start_pos..=end_pos);

        let query = self
            .reader
            .query(&self.header, &region)
            .map_err(|e| CyaneaError::Parse(format!("BAM query failed: {e}")))?;

        let mut records = Vec::new();
        for result in query {
            let noodles_rec = result.map_err(|e| {
                CyaneaError::Parse(format!("error reading BAM record: {e}"))
            })?;

            let sam = convert_noodles_record(&noodles_rec, &self.header)?;
            records.push(sam);
        }

        Ok(records)
    }

    /// Get reference sequence information from the BAM header.
    pub fn references(&self) -> Vec<BamReference> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, map)| BamReference {
                name: String::from_utf8_lossy(name).to_string(),
                length: usize::from(map.length()) as u32,
            })
            .collect()
    }
}

/// Create a BAI index from a coordinate-sorted BAM file.
///
/// Reads the BAM file, tracks virtual positions and reference sequence placements,
/// then writes a standard BAI index to `bai_path`.
pub fn create_bai_index(bam_path: &Path, bai_path: &Path) -> Result<()> {
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles_csi::binning_index::index::reference_sequence::index::LinearIndex;
    use noodles_csi::binning_index::Indexer;
    use noodles_sam::alignment::Record;

    let file = std::fs::File::open(bam_path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", bam_path.display(), e),
        ))
    })?;

    let mut reader = noodles_bam_crate::io::Reader::new(file);
    let header = reader.read_header().map_err(|e| {
        CyaneaError::Parse(format!("failed to read BAM header: {e}"))
    })?;

    let mut indexer = Indexer::<LinearIndex>::new(14, 5);
    let mut record = noodles_bam_crate::Record::default();

    loop {
        let start_vpos = reader.get_ref().virtual_position();
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => {
                return Err(CyaneaError::Parse(format!(
                    "failed to read BAM record during indexing: {e}"
                )));
            }
        }
        let end_vpos = reader.get_ref().virtual_position();

        let chunk = Chunk::new(start_vpos, end_vpos);

        let ref_seq_id = record.reference_sequence_id().and_then(|r| r.ok());
        let start = record.alignment_start().and_then(|r| r.ok());
        let end = Record::alignment_end(&record).and_then(|r| r.ok());
        let is_mapped = !record.flags().is_unmapped();

        if let (Some(ref_id), Some(s)) = (ref_seq_id, start) {
            let e = end.unwrap_or(s);
            indexer
                .add_record(Some((ref_id, s, e, is_mapped)), chunk)
                .map_err(|e| CyaneaError::Parse(format!("indexer error: {e}")))?;
        } else {
            indexer
                .add_record(None, chunk)
                .map_err(|e| CyaneaError::Parse(format!("indexer error: {e}")))?;
        }
    }

    let ref_count = header.reference_sequences().len();
    let index = indexer.build(ref_count);

    noodles_bam_crate::bai::write(bai_path, &index).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("{}: {}", bai_path.display(), e),
        ))
    })?;

    Ok(())
}

/// Convenience function: open a BAM file and fetch records in a region.
pub fn fetch_bam(
    path: impl AsRef<Path>,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<SamRecord>> {
    let mut reader = IndexedBamReader::open(path)?;
    reader.fetch(chrom, start, end)
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

/// Convert a noodles BAM record to our SamRecord type.
fn convert_noodles_record(
    record: &noodles_bam_crate::Record,
    header: &noodles_sam::Header,
) -> Result<SamRecord> {
    use noodles_sam::alignment::Record;

    let qname = record
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
        .unwrap_or_else(|| "*".to_string());

    let flags = record.flags();
    let flag = u16::from(flags);

    let rname = record
        .reference_sequence(header)
        .and_then(|r| r.ok())
        .map(|(name, _)| String::from_utf8_lossy(name).to_string())
        .unwrap_or_else(|| "*".to_string());

    let pos = record
        .alignment_start()
        .and_then(|r| r.ok())
        .map(|p| usize::from(p) as u64)
        .unwrap_or(0);

    let mapq = record
        .mapping_quality()
        .map(|m| u8::from(m))
        .unwrap_or(255);

    let cigar = {
        let cigar_obj = record.cigar();
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
            if s.is_empty() {
                "*".to_string()
            } else {
                s
            }
        }
    };

    let sequence = {
        let seq = record.sequence();
        if seq.is_empty() {
            "*".to_string()
        } else {
            // BAM sequences are packed 2 bases per byte; iterate to decode
            let decoded: Vec<u8> = seq.iter().collect();
            String::from_utf8_lossy(&decoded).to_string()
        }
    };

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

    let rnext = record
        .mate_reference_sequence(header)
        .and_then(|r| r.ok())
        .map(|(name, _)| {
            let mate_name = String::from_utf8_lossy(name).to_string();
            if mate_name == rname {
                "=".to_string()
            } else {
                mate_name
            }
        })
        .unwrap_or_else(|| "*".to_string());

    let pnext = record
        .mate_alignment_start()
        .and_then(|r| r.ok())
        .map(|p| usize::from(p) as u64)
        .unwrap_or(0);

    let tlen = record.template_length() as i64;

    Ok(SamRecord {
        qname,
        flag,
        rname,
        pos,
        mapq,
        cigar,
        rnext,
        pnext,
        tlen,
        sequence,
        quality,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::io::Write as SamIoWrite;
    use std::io::Write;

    /// Write a BAM file using noodles and create a BAI index.
    fn write_indexed_bam(
        records: &[noodles_sam::alignment::RecordBuf],
    ) -> (tempfile::NamedTempFile, tempfile::NamedTempFile) {
        use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles_csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles_csi::binning_index::Indexer;
        use noodles_sam::alignment::Record;
        use noodles_sam::header::record::value::map::ReferenceSequence;
        use noodles_sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        let len1 = NonZeroUsize::try_from(248956422usize).unwrap();
        let len2 = NonZeroUsize::try_from(242193529usize).unwrap();
        let header = noodles_sam::Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len1))
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(len2))
            .build();

        // Write BAM to buffer
        let mut bam_buf = Vec::new();
        {
            let mut writer = noodles_bam_crate::io::Writer::new(&mut bam_buf);
            writer.write_header(&header).unwrap();

            for rec in records {
                writer
                    .write_alignment_record(&header, rec)
                    .unwrap();
            }
            writer.try_finish().unwrap();
        }

        // Write BAM to temp file
        let mut bam_file = tempfile::NamedTempFile::with_suffix(".bam").unwrap();
        bam_file.write_all(&bam_buf).unwrap();
        bam_file.flush().unwrap();

        // Build BAI index by reading back the BAM and tracking virtual positions
        let index = {
            let mut reader =
                noodles_bam_crate::io::Reader::new(std::io::Cursor::new(&bam_buf));
            reader.read_header().unwrap();

            let mut indexer = Indexer::<LinearIndex>::new(14, 5);
            let mut record = noodles_bam_crate::Record::default();

            loop {
                let start_vpos = reader.get_ref().virtual_position();
                match reader.read_record(&mut record) {
                    Ok(0) => break,
                    Ok(_) => {}
                    Err(e) => panic!("failed to read BAM record: {e}"),
                }
                let end_vpos = reader.get_ref().virtual_position();

                let chunk = Chunk::new(start_vpos, end_vpos);

                let ref_seq_id = record.reference_sequence_id().and_then(|r| r.ok());
                let start = record.alignment_start().and_then(|r| r.ok());
                let end = Record::alignment_end(&record).and_then(|r| r.ok());
                let is_mapped = !record.flags().is_unmapped();

                if let (Some(ref_id), Some(s)) = (ref_seq_id, start) {
                    let e = end.unwrap_or(s);
                    indexer
                        .add_record(Some((ref_id, s, e, is_mapped)), chunk)
                        .unwrap();
                } else {
                    indexer.add_record(None, chunk).unwrap();
                }
            }

            let ref_count = header.reference_sequences().len();
            indexer.build(ref_count)
        };

        let bai_file = tempfile::NamedTempFile::with_suffix(".bai").unwrap();
        noodles_bam_crate::bai::write(bai_file.path(), &index).unwrap();

        (bam_file, bai_file)
    }

    fn make_record(
        ref_id: usize,
        pos: usize,
        mapq: u8,
        name: &str,
        seq: &str,
    ) -> noodles_sam::alignment::RecordBuf {
        use noodles_sam::alignment::record::cigar::op::Kind;
        use noodles_sam::alignment::record::cigar::Op;
        use noodles_sam::alignment::record::Flags;
        use noodles_sam::alignment::record::MappingQuality;
        use noodles_sam::alignment::record_buf::{
            Cigar as CigarBuf, QualityScores as QualBuf, Sequence as SeqBuf,
        };

        let cigar = CigarBuf::from(vec![Op::new(Kind::Match, seq.len())]);

        noodles_sam::alignment::RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::from(0u16))
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(noodles_core::Position::try_from(pos).unwrap())
            .set_mapping_quality(MappingQuality::new(mapq).unwrap())
            .set_cigar(cigar)
            .set_sequence(SeqBuf::from(seq.as_bytes().to_vec()))
            .set_quality_scores(QualBuf::from(vec![30u8; seq.len()]))
            .build()
    }

    #[test]
    fn indexed_bam_fetch_region() {
        let records = vec![
            make_record(0, 100, 60, "read1", "ACGTACGTAC"),
            make_record(0, 200, 40, "read2", "TGCATGCATG"),
            make_record(0, 500, 50, "read3", "GGCCGGCCGG"),
        ];

        let (bam_file, bai_file) = write_indexed_bam(&records);
        let mut reader =
            IndexedBamReader::open_with_index(bam_file.path(), bai_file.path()).unwrap();

        // Fetch region covering read1 and read2 but not read3
        let results = reader.fetch("chr1", 50, 250).unwrap();
        assert!(results.len() >= 2);

        let names: Vec<&str> = results.iter().map(|r| r.qname.as_str()).collect();
        assert!(names.contains(&"read1"));
        assert!(names.contains(&"read2"));
    }

    #[test]
    fn indexed_bam_empty_region() {
        let records = vec![make_record(0, 100, 60, "read1", "ACGTACGTAC")];

        let (bam_file, bai_file) = write_indexed_bam(&records);
        let mut reader =
            IndexedBamReader::open_with_index(bam_file.path(), bai_file.path()).unwrap();

        let results = reader.fetch("chr1", 10000, 20000).unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn indexed_bam_references() {
        let (bam_file, bai_file) = write_indexed_bam(&[]);
        let reader =
            IndexedBamReader::open_with_index(bam_file.path(), bai_file.path()).unwrap();
        let refs = reader.references();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs[0].name, "chr1");
        assert_eq!(refs[1].name, "chr2");
    }

    #[test]
    fn create_bai_roundtrip() {
        let records = vec![
            make_record(0, 100, 60, "read1", "ACGTACGTAC"),
            make_record(0, 200, 40, "read2", "TGCATGCATG"),
            make_record(0, 500, 50, "read3", "GGCCGGCCGG"),
        ];

        // Write BAM using the test helper (which creates its own BAI)
        let (bam_file, _original_bai) = write_indexed_bam(&records);

        // Now create a new BAI using the public function
        let new_bai = tempfile::NamedTempFile::with_suffix(".bai").unwrap();
        create_bai_index(bam_file.path(), new_bai.path()).unwrap();

        // Verify we can fetch records using the newly created index
        let mut reader =
            IndexedBamReader::open_with_index(bam_file.path(), new_bai.path()).unwrap();
        let results = reader.fetch("chr1", 50, 250).unwrap();
        assert!(results.len() >= 2);
        let names: Vec<&str> = results.iter().map(|r| r.qname.as_str()).collect();
        assert!(names.contains(&"read1"));
        assert!(names.contains(&"read2"));
    }

    #[test]
    fn create_bai_empty() {
        // Empty BAM file → valid empty index
        let (bam_file, _) = write_indexed_bam(&[]);

        let new_bai = tempfile::NamedTempFile::with_suffix(".bai").unwrap();
        create_bai_index(bam_file.path(), new_bai.path()).unwrap();

        // Verify the index works (no records to fetch)
        let mut reader =
            IndexedBamReader::open_with_index(bam_file.path(), new_bai.path()).unwrap();
        let results = reader.fetch("chr1", 0, 1000).unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn indexed_bam_invalid_chrom() {
        let (bam_file, bai_file) = write_indexed_bam(&[]);
        let mut reader =
            IndexedBamReader::open_with_index(bam_file.path(), bai_file.path()).unwrap();
        let result = reader.fetch("chrINVALID", 0, 100);
        assert!(result.is_err());
    }
}
