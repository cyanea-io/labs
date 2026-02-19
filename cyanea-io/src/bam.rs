//! BAM (Binary Alignment/Map) parser with BGZF decompression.
//!
//! BAM is the binary equivalent of SAM, compressed with BGZF (Blocked GNU Zip Format).
//! BGZF is a series of concatenated gzip blocks, each at most 64 KiB uncompressed,
//! enabling random access when paired with a BAI index.
//!
//! This module provides sequential reading of BAM files, reusing the [`SamRecord`]
//! type from the SAM module.

use std::io::{self, Read};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

#[cfg(feature = "sam")]
use crate::sam::{sam_stats, SamRecord, SamStats};

use crate::bgzf;

/// BAM magic bytes at the start of decompressed data.
const BAM_MAGIC: [u8; 4] = [b'B', b'A', b'M', 0x01];

// ---------------------------------------------------------------------------
// BAM reader
// ---------------------------------------------------------------------------

/// A reference sequence from the BAM header.
#[derive(Debug, Clone)]
pub struct BamReference {
    /// Reference sequence name.
    pub name: String,
    /// Reference sequence length.
    pub length: u32,
}

/// Parse a BAM file and return all alignment records as [`SamRecord`]s.
///
/// Decompresses BGZF blocks, parses the BAM header and reference dictionary,
/// then reads all alignment records.
///
/// # Errors
///
/// Returns an error if the file cannot be opened, is not valid BAM/BGZF,
/// or contains malformed alignment records.
#[cfg(feature = "sam")]
pub fn parse_bam(path: impl AsRef<Path>) -> Result<Vec<SamRecord>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path).map_err(|e| {
        CyaneaError::Io(io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let mut reader = io::BufReader::new(file);
    parse_bam_reader(&mut reader)
}

/// Parse BAM records from a reader.
#[cfg(feature = "sam")]
fn parse_bam_reader(reader: &mut impl Read) -> Result<Vec<SamRecord>> {
    // Decompress all BGZF blocks into a contiguous buffer
    let data = bgzf::decompress_all(reader)?;

    if data.len() < 4 {
        return Err(CyaneaError::Parse("BAM data too short".into()));
    }

    let mut pos = 0;

    // Validate BAM magic
    if data[pos..pos + 4] != BAM_MAGIC {
        return Err(CyaneaError::Parse("not a valid BAM file (bad magic)".into()));
    }
    pos += 4;

    // Header text length + header text
    let header_len = read_u32_le(&data, &mut pos)? as usize;
    if pos + header_len > data.len() {
        return Err(CyaneaError::Parse("BAM header truncated".into()));
    }
    pos += header_len; // Skip header text

    // Reference sequences
    let n_ref = read_u32_le(&data, &mut pos)? as usize;
    let mut references = Vec::with_capacity(n_ref);
    for _ in 0..n_ref {
        let name_len = read_u32_le(&data, &mut pos)? as usize;
        if pos + name_len > data.len() {
            return Err(CyaneaError::Parse("BAM reference name truncated".into()));
        }
        // name_len includes NUL terminator
        let name = std::str::from_utf8(&data[pos..pos + name_len - 1])
            .map_err(|_| CyaneaError::Parse("invalid UTF-8 in reference name".into()))?
            .to_string();
        pos += name_len;
        let length = read_u32_le(&data, &mut pos)?;
        references.push(BamReference { name, length });
    }

    // Read alignment records
    let mut records = Vec::new();
    while pos < data.len() {
        if pos + 4 > data.len() {
            break;
        }
        let record = parse_bam_record(&data, &mut pos, &references)?;
        records.push(record);
    }

    Ok(records)
}

/// Parse a single BAM alignment record.
#[cfg(feature = "sam")]
fn parse_bam_record(
    data: &[u8],
    pos: &mut usize,
    references: &[BamReference],
) -> Result<SamRecord> {
    let block_size = read_u32_le(data, pos)? as usize;
    if *pos + block_size > data.len() {
        return Err(CyaneaError::Parse("BAM record truncated".into()));
    }

    let record_start = *pos;

    let ref_id = read_i32_le(data, pos)?;
    let position = read_i32_le(data, pos)?; // 0-based
    let bin_mq_nl = read_u32_le(data, pos)?;
    let _bin = (bin_mq_nl >> 16) as u16;
    let mapq = ((bin_mq_nl >> 8) & 0xFF) as u8;
    let name_len = (bin_mq_nl & 0xFF) as usize;

    let flag_nc = read_u32_le(data, pos)?;
    let flag = (flag_nc >> 16) as u16;
    let n_cigar_op = (flag_nc & 0xFFFF) as usize;

    let seq_len = read_u32_le(data, pos)? as usize;
    let next_ref_id = read_i32_le(data, pos)?;
    let next_pos = read_i32_le(data, pos)?;
    let tlen = read_i32_le(data, pos)?;

    // Read name (NUL-terminated)
    if *pos + name_len > data.len() {
        return Err(CyaneaError::Parse("BAM record name truncated".into()));
    }
    let qname = std::str::from_utf8(&data[*pos..*pos + name_len - 1])
        .map_err(|_| CyaneaError::Parse("invalid UTF-8 in read name".into()))?
        .to_string();
    *pos += name_len;

    // Read CIGAR
    if *pos + n_cigar_op * 4 > data.len() {
        return Err(CyaneaError::Parse("BAM CIGAR truncated".into()));
    }
    let cigar = decode_cigar(data, pos, n_cigar_op);

    // Read sequence (4-bit encoded, 2 bases per byte)
    let seq_bytes = (seq_len + 1) / 2;
    if *pos + seq_bytes > data.len() {
        return Err(CyaneaError::Parse("BAM sequence truncated".into()));
    }
    let sequence = decode_sequence(data, *pos, seq_len);
    *pos += seq_bytes;

    // Read quality (one byte per base)
    if *pos + seq_len > data.len() {
        return Err(CyaneaError::Parse("BAM quality truncated".into()));
    }
    let quality = decode_quality(data, *pos, seq_len);
    *pos += seq_len;

    // Skip auxiliary data (remaining bytes in block)
    let consumed = *pos - record_start;
    if consumed < block_size {
        *pos = record_start + block_size;
    }

    // Convert reference ID to name
    let rname = if ref_id < 0 {
        "*".to_string()
    } else {
        references
            .get(ref_id as usize)
            .map(|r| r.name.clone())
            .unwrap_or_else(|| "*".to_string())
    };

    // Convert mate reference ID to RNEXT
    let rnext = if next_ref_id < 0 {
        "*".to_string()
    } else if next_ref_id == ref_id {
        "=".to_string()
    } else {
        references
            .get(next_ref_id as usize)
            .map(|r| r.name.clone())
            .unwrap_or_else(|| "*".to_string())
    };

    // BAM position is 0-based; SAM position is 1-based
    let sam_pos = if position < 0 {
        0u64
    } else {
        (position + 1) as u64
    };

    // Mate position: 0-based to 1-based, 0 if unavailable (-1 in BAM)
    let pnext = if next_pos < 0 { 0u64 } else { (next_pos + 1) as u64 };

    Ok(SamRecord {
        qname,
        flag,
        rname,
        pos: sam_pos,
        mapq,
        cigar,
        rnext,
        pnext,
        tlen: tlen as i64,
        sequence,
        quality,
    })
}

// ---------------------------------------------------------------------------
// BAM field decoders
// ---------------------------------------------------------------------------

/// Decode CIGAR operations from BAM binary format.
fn decode_cigar(data: &[u8], pos: &mut usize, n_ops: usize) -> String {
    if n_ops == 0 {
        return "*".to_string();
    }

    let mut cigar = String::new();
    for _ in 0..n_ops {
        let val = u32::from_le_bytes([
            data[*pos],
            data[*pos + 1],
            data[*pos + 2],
            data[*pos + 3],
        ]);
        *pos += 4;

        let op_len = val >> 4;
        let op_code = val & 0xF;
        let op_char = match op_code {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            _ => '?',
        };
        cigar.push_str(&format!("{}{}", op_len, op_char));
    }
    cigar
}

/// Decode BAM 4-bit encoded sequence to ASCII.
const SEQ_DECODE: [u8; 16] = [
    b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'N',
    b'N',
];

fn decode_sequence(data: &[u8], offset: usize, seq_len: usize) -> String {
    if seq_len == 0 {
        return "*".to_string();
    }

    let mut seq = Vec::with_capacity(seq_len);
    for i in 0..seq_len {
        let byte = data[offset + i / 2];
        let base = if i % 2 == 0 {
            SEQ_DECODE[(byte >> 4) as usize]
        } else {
            SEQ_DECODE[(byte & 0x0F) as usize]
        };
        seq.push(base);
    }
    String::from_utf8(seq).unwrap_or_else(|_| "*".to_string())
}

/// Decode BAM quality scores to Phred+33 ASCII.
fn decode_quality(data: &[u8], offset: usize, seq_len: usize) -> String {
    if seq_len == 0 {
        return "*".to_string();
    }

    // Check if all quality values are 0xFF (unavailable)
    let all_missing = data[offset..offset + seq_len].iter().all(|&q| q == 0xFF);
    if all_missing {
        return "*".to_string();
    }

    let qual: Vec<u8> = data[offset..offset + seq_len]
        .iter()
        .map(|&q| q + 33) // Convert to Phred+33 ASCII
        .collect();
    String::from_utf8(qual).unwrap_or_else(|_| "*".to_string())
}

// ---------------------------------------------------------------------------
// Binary reading helpers
// ---------------------------------------------------------------------------

fn read_u32_le(data: &[u8], pos: &mut usize) -> Result<u32> {
    if *pos + 4 > data.len() {
        return Err(CyaneaError::Parse("unexpected end of BAM data".into()));
    }
    let val = u32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
    *pos += 4;
    Ok(val)
}

fn read_i32_le(data: &[u8], pos: &mut usize) -> Result<i32> {
    if *pos + 4 > data.len() {
        return Err(CyaneaError::Parse("unexpected end of BAM data".into()));
    }
    let val = i32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
    *pos += 4;
    Ok(val)
}

/// Compute BAM statistics from a BAM file path.
///
/// Parses the BAM file and computes summary statistics using the same
/// [`SamStats`] type as the SAM module.
#[cfg(feature = "sam")]
pub fn bam_stats(path: impl AsRef<Path>) -> Result<SamStats> {
    let records = parse_bam(path)?;
    Ok(sam_stats(&records))
}

// ---------------------------------------------------------------------------
// Write a minimal valid BAM for testing
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test_helpers {
    pub(crate) use crate::bgzf::test_helpers::{bgzf_compress, bgzf_eof_block};

    /// Build a minimal valid BAM binary (uncompressed content).
    ///
    /// Creates a BAM with one reference (chr1, length 248956422) and
    /// the given alignment records as raw BAM binary.
    pub(super) fn build_bam_content(records: &[BamTestRecord]) -> Vec<u8> {
        let mut data = Vec::new();

        // BAM magic
        data.extend_from_slice(b"BAM\x01");

        // Header text
        let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
        data.extend_from_slice(&(header.len() as u32).to_le_bytes());
        data.extend_from_slice(header);

        // Reference sequences
        data.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        let ref_name = b"chr1\0";
        data.extend_from_slice(&(ref_name.len() as u32).to_le_bytes());
        data.extend_from_slice(ref_name);
        data.extend_from_slice(&248956422u32.to_le_bytes());

        // Alignment records
        for rec in records {
            let record_bytes = encode_bam_record(rec);
            data.extend_from_slice(&(record_bytes.len() as u32).to_le_bytes());
            data.extend_from_slice(&record_bytes);
        }

        data
    }

    pub(super) struct BamTestRecord {
        pub qname: &'static str,
        pub flag: u16,
        pub ref_id: i32,
        pub pos: i32, // 0-based
        pub mapq: u8,
        pub cigar: Vec<(u32, u8)>, // (length, op_code)
        pub seq: &'static [u8],
        pub qual: &'static [u8], // raw Phred values (not +33)
        pub next_ref_id: i32,
        pub next_pos: i32, // 0-based
        pub tlen: i32,
    }

    fn encode_bam_record(rec: &BamTestRecord) -> Vec<u8> {
        let mut data = Vec::new();

        // refID
        data.extend_from_slice(&rec.ref_id.to_le_bytes());
        // pos (0-based)
        data.extend_from_slice(&rec.pos.to_le_bytes());

        // bin_mq_nl
        let name_len = rec.qname.len() + 1; // includes NUL
        let bin = reg2bin(rec.pos, rec.pos + rec.seq.len() as i32) as u32;
        let bin_mq_nl = (bin << 16) | ((rec.mapq as u32) << 8) | (name_len as u32);
        data.extend_from_slice(&bin_mq_nl.to_le_bytes());

        // flag_nc
        let flag_nc = ((rec.flag as u32) << 16) | (rec.cigar.len() as u32);
        data.extend_from_slice(&flag_nc.to_le_bytes());

        // seq length
        data.extend_from_slice(&(rec.seq.len() as u32).to_le_bytes());

        // next_refID, next_pos, tlen
        data.extend_from_slice(&rec.next_ref_id.to_le_bytes());
        data.extend_from_slice(&rec.next_pos.to_le_bytes());
        data.extend_from_slice(&rec.tlen.to_le_bytes());

        // read name (NUL-terminated)
        data.extend_from_slice(rec.qname.as_bytes());
        data.push(0);

        // CIGAR
        for &(len, op) in &rec.cigar {
            let val = (len << 4) | (op as u32);
            data.extend_from_slice(&val.to_le_bytes());
        }

        // Sequence (4-bit encoded)
        let seq_bytes = (rec.seq.len() + 1) / 2;
        for i in 0..seq_bytes {
            let high = encode_base(rec.seq.get(i * 2).copied().unwrap_or(0));
            let low = if i * 2 + 1 < rec.seq.len() {
                encode_base(rec.seq[i * 2 + 1])
            } else {
                0
            };
            data.push((high << 4) | low);
        }

        // Quality
        if rec.qual.is_empty() {
            // All 0xFF means unavailable
            data.extend(std::iter::repeat(0xFF).take(rec.seq.len()));
        } else {
            data.extend_from_slice(rec.qual);
        }

        data
    }

    fn encode_base(b: u8) -> u8 {
        match b.to_ascii_uppercase() {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            b'N' => 15,
            _ => 0,
        }
    }

    /// Compute BAM bin for a 0-based [beg, end) region.
    fn reg2bin(beg: i32, end: i32) -> u16 {
        let beg = beg as u32;
        let end = (end - 1) as u32;
        if beg >> 14 == end >> 14 {
            return (((1 << 15) - 1) / 7 + (beg >> 14)) as u16;
        }
        if beg >> 17 == end >> 17 {
            return (((1 << 12) - 1) / 7 + (beg >> 17)) as u16;
        }
        if beg >> 20 == end >> 20 {
            return (((1 << 9) - 1) / 7 + (beg >> 20)) as u16;
        }
        if beg >> 23 == end >> 23 {
            return (((1 << 6) - 1) / 7 + (beg >> 23)) as u16;
        }
        if beg >> 26 == end >> 26 {
            return (((1 << 3) - 1) / 7 + (beg >> 26)) as u16;
        }
        0
    }

    /// Build a complete BAM file (BGZF-compressed) and write to a temp file.
    pub(super) fn write_test_bam(
        records: &[BamTestRecord],
    ) -> tempfile::NamedTempFile {
        let content = build_bam_content(records);
        let compressed = bgzf_compress(&content);
        let eof = bgzf_eof_block();

        let mut file = tempfile::NamedTempFile::with_suffix(".bam").unwrap();
        use std::io::Write;
        file.write_all(&compressed).unwrap();
        file.write_all(&eof).unwrap();
        file.flush().unwrap();
        file
    }
}

#[cfg(test)]
#[cfg(feature = "sam")]
mod tests {
    use super::test_helpers::*;
    use super::*;

    #[test]
    fn parse_minimal_bam() {
        let records = vec![BamTestRecord {
            qname: "read1",
            flag: 0,
            ref_id: 0,
            pos: 99, // 0-based → SAM pos 100
            mapq: 60,
            cigar: vec![(10, 0)], // 10M
            seq: b"ACGTACGTAC",
            qual: &[],
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        }];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].qname, "read1");
        assert_eq!(result[0].flag, 0);
        assert_eq!(result[0].rname, "chr1");
        assert_eq!(result[0].pos, 100); // 1-based
        assert_eq!(result[0].mapq, 60);
        assert_eq!(result[0].cigar, "10M");
        assert_eq!(result[0].sequence, "ACGTACGTAC");
        assert!(result[0].is_mapped());
        assert_eq!(result[0].rnext, "*");
        assert_eq!(result[0].pnext, 0);
        assert_eq!(result[0].tlen, 0);
    }

    #[test]
    fn parse_bam_multiple_records() {
        let records = vec![
            BamTestRecord {
                qname: "read1",
                flag: 0,
                ref_id: 0,
                pos: 99,
                mapq: 60,
                cigar: vec![(4, 0)], // 4M
                seq: b"ACGT",
                qual: &[30, 30, 30, 30],
                next_ref_id: -1,
                next_pos: -1,
                tlen: 0,
            },
            BamTestRecord {
                qname: "read2",
                flag: 16, // reverse strand
                ref_id: 0,
                pos: 199,
                mapq: 40,
                cigar: vec![(4, 0)],
                seq: b"TGCA",
                qual: &[],
                next_ref_id: -1,
                next_pos: -1,
                tlen: 0,
            },
        ];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result.len(), 2);

        assert_eq!(result[0].qname, "read1");
        assert_eq!(result[0].pos, 100);
        assert_eq!(result[0].mapq, 60);
        assert_eq!(result[0].quality, "????"); // 30+33=63='?'

        assert_eq!(result[1].qname, "read2");
        assert_eq!(result[1].flag, 16);
        assert_eq!(result[1].pos, 200);
        assert_eq!(result[1].mapq, 40);
    }

    #[test]
    fn parse_bam_unmapped_read() {
        let records = vec![BamTestRecord {
            qname: "unmapped1",
            flag: 4, // unmapped
            ref_id: -1,
            pos: -1,
            mapq: 0,
            cigar: vec![],
            seq: b"AAAA",
            qual: &[],
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        }];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].qname, "unmapped1");
        assert!(result[0].is_unmapped());
        assert_eq!(result[0].rname, "*");
        assert_eq!(result[0].cigar, "*");
    }

    #[test]
    fn bam_stats_basic() {
        let records = vec![
            BamTestRecord {
                qname: "read1",
                flag: 0,
                ref_id: 0,
                pos: 99,
                mapq: 60,
                cigar: vec![(4, 0)],
                seq: b"ACGT",
                qual: &[],
                next_ref_id: -1,
                next_pos: -1,
                tlen: 0,
            },
            BamTestRecord {
                qname: "read2",
                flag: 0,
                ref_id: 0,
                pos: 199,
                mapq: 30,
                cigar: vec![(4, 0)],
                seq: b"TGCA",
                qual: &[],
                next_ref_id: -1,
                next_pos: -1,
                tlen: 0,
            },
            BamTestRecord {
                qname: "read3",
                flag: 4,
                ref_id: -1,
                pos: -1,
                mapq: 0,
                cigar: vec![],
                seq: b"NNNN",
                qual: &[],
                next_ref_id: -1,
                next_pos: -1,
                tlen: 0,
            },
        ];
        let file = write_test_bam(&records);
        let stats = bam_stats(file.path()).unwrap();
        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.mapped, 2);
        assert_eq!(stats.unmapped, 1);
        assert!((stats.avg_mapq - 45.0).abs() < f64::EPSILON);
    }

    #[test]
    fn bam_cigar_operations() {
        // Test various CIGAR operations
        let records = vec![BamTestRecord {
            qname: "read_complex",
            flag: 0,
            ref_id: 0,
            pos: 99,
            mapq: 60,
            cigar: vec![(3, 0), (1, 1), (2, 0), (1, 2), (1, 0)], // 3M1I2M1D1M
            seq: b"ACGTTGA",                                       // 3+1+2+1 = 7 bases
            qual: &[],
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        }];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result[0].cigar, "3M1I2M1D1M");
    }

    #[test]
    fn bam_empty_file() {
        // BAM with no records
        let file = write_test_bam(&[]);
        let result = parse_bam(file.path()).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn bam_sequence_encoding() {
        // Verify all standard bases decode correctly
        let records = vec![BamTestRecord {
            qname: "bases",
            flag: 0,
            ref_id: 0,
            pos: 0,
            mapq: 60,
            cigar: vec![(5, 0)],
            seq: b"ACGTN",
            qual: &[],
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        }];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result[0].sequence, "ACGTN");
    }

    #[test]
    fn bam_mate_fields_same_ref() {
        // Paired reads on same reference: next_ref_id == ref_id → rnext = "="
        let records = vec![
            BamTestRecord {
                qname: "read1",
                flag: 99, // paired, proper, mate_reverse, first_in_pair
                ref_id: 0,
                pos: 99,
                mapq: 60,
                cigar: vec![(10, 0)],
                seq: b"ACGTACGTAC",
                qual: &[],
                next_ref_id: 0, // same as ref_id
                next_pos: 199,  // 0-based
                tlen: 150,
            },
            BamTestRecord {
                qname: "read1",
                flag: 147, // paired, proper, reverse, second_in_pair
                ref_id: 0,
                pos: 199,
                mapq: 60,
                cigar: vec![(10, 0)],
                seq: b"ACGTACGTAC",
                qual: &[],
                next_ref_id: 0,
                next_pos: 99,
                tlen: -150,
            },
        ];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result.len(), 2);

        assert_eq!(result[0].rnext, "=");
        assert_eq!(result[0].pnext, 200); // 0-based 199 → 1-based 200
        assert_eq!(result[0].tlen, 150);

        assert_eq!(result[1].rnext, "=");
        assert_eq!(result[1].pnext, 100);
        assert_eq!(result[1].tlen, -150);
    }

    #[test]
    fn bam_mate_fields_unmapped_mate() {
        let records = vec![BamTestRecord {
            qname: "read1",
            flag: 73, // paired + mate_unmapped + first_in_pair
            ref_id: 0,
            pos: 99,
            mapq: 60,
            cigar: vec![(4, 0)],
            seq: b"ACGT",
            qual: &[],
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        }];
        let file = write_test_bam(&records);
        let result = parse_bam(file.path()).unwrap();
        assert_eq!(result[0].rnext, "*");
        assert_eq!(result[0].pnext, 0);
        assert_eq!(result[0].tlen, 0);
    }

    #[test]
    fn bgzf_block_decompression() {
        let original = b"Hello, BGZF world!";
        let compressed = bgzf_compress(original);
        let mut reader = io::Cursor::new(compressed);
        let decompressed = crate::bgzf::read_bgzf_block(&mut reader).unwrap().unwrap();
        assert_eq!(decompressed, original);
    }

    #[test]
    fn bgzf_eof_detection() {
        let eof = bgzf_eof_block();
        let mut reader = io::Cursor::new(eof);
        let block = crate::bgzf::read_bgzf_block(&mut reader).unwrap().unwrap();
        assert!(block.is_empty());
    }

    #[test]
    fn bam_file_not_found() {
        let result = parse_bam("/nonexistent/file.bam");
        assert!(result.is_err());
    }
}
