//! BCF (Binary Call Format) reader.
//!
//! BCF2 is the binary encoding of VCF, using BGZF compression (same as BAM).
//! Records contain integer-encoded contig and allele information, with the
//! VCF text header stored in the first BGZF block for field definitions.
//!
//! This module reuses [`Variant`] and [`VcfStats`] from the VCF module.

use std::io::{self, Read};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::variant::{Variant, VariantFilter};
use flate2::read::DeflateDecoder;

use crate::vcf::VcfStats;

/// BCF2 magic: "BCF\2\1"
const BCF_MAGIC: [u8; 5] = [b'B', b'C', b'F', 0x02, 0x01];

/// BGZF magic bytes: standard gzip header with FEXTRA flag.
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Parse a BCF2 file and return variant records.
///
/// Decompresses BGZF blocks, parses the BCF header for contig and filter
/// definitions, then reads all variant records.
pub fn parse_bcf(path: impl AsRef<Path>) -> Result<Vec<Variant>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path).map_err(|e| {
        CyaneaError::Io(io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let mut reader = io::BufReader::new(file);

    // Decompress all BGZF blocks
    let mut data = Vec::new();
    loop {
        match read_bgzf_block(&mut reader)? {
            None => break,
            Some(block) => {
                if block.is_empty() {
                    break;
                }
                data.extend_from_slice(&block);
            }
        }
    }

    if data.len() < 9 {
        return Err(CyaneaError::Parse(format!(
            "{}: BCF data too short",
            path.display()
        )));
    }

    let mut pos = 0;

    // Validate BCF magic
    if data[pos..pos + 5] != BCF_MAGIC {
        return Err(CyaneaError::Parse(format!(
            "{}: not a valid BCF file (bad magic)",
            path.display()
        )));
    }
    pos += 5;

    // Header text length (u32 LE)
    let header_len = read_u32_le(&data, &mut pos)? as usize;
    if pos + header_len > data.len() {
        return Err(CyaneaError::Parse(format!(
            "{}: BCF header truncated",
            path.display()
        )));
    }

    // Parse header text for contig names and filter definitions
    let header_text = std::str::from_utf8(&data[pos..pos + header_len])
        .map_err(|_| CyaneaError::Parse(format!("{}: invalid UTF-8 in BCF header", path.display())))?;
    pos += header_len;

    let contigs = parse_header_contigs(header_text);
    let filters = parse_header_filters(header_text);

    // Read variant records
    let mut variants = Vec::new();
    while pos < data.len() {
        if pos + 8 > data.len() {
            break;
        }
        let variant = parse_bcf_record(&data, &mut pos, &contigs, &filters)?;
        variants.push(variant);
    }

    Ok(variants)
}

/// Compute variant statistics from a BCF file.
///
/// Returns the same [`VcfStats`] type used by the VCF module.
pub fn bcf_stats(path: impl AsRef<Path>) -> Result<VcfStats> {
    let variants = parse_bcf(path)?;

    let mut snv_count: u64 = 0;
    let mut indel_count: u64 = 0;
    let mut pass_count: u64 = 0;
    let mut chroms = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for v in &variants {
        if v.is_snv() {
            snv_count += 1;
        }
        if v.is_indel() {
            indel_count += 1;
        }
        if v.filter == VariantFilter::Pass {
            pass_count += 1;
        }
        if seen.insert(v.chrom.clone()) {
            chroms.push(v.chrom.clone());
        }
    }

    Ok(VcfStats {
        variant_count: variants.len() as u64,
        snv_count,
        indel_count,
        pass_count,
        chromosomes: chroms,
    })
}

// ---------------------------------------------------------------------------
// BGZF block reader (same pattern as bam.rs)
// ---------------------------------------------------------------------------

fn read_bgzf_block(reader: &mut impl Read) -> Result<Option<Vec<u8>>> {
    let mut header = [0u8; 18];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(CyaneaError::Io(e)),
    }

    if header[0] != BGZF_MAGIC[0]
        || header[1] != BGZF_MAGIC[1]
        || header[2] != BGZF_MAGIC[2]
        || (header[3] & 0x04) == 0
    {
        return Err(CyaneaError::Parse("not a valid BGZF block".into()));
    }

    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;

    let mut extra = vec![0u8; xlen];
    extra[..6].copy_from_slice(&header[12..18]);
    if xlen > 6 {
        reader
            .read_exact(&mut extra[6..])
            .map_err(CyaneaError::Io)?;
    }

    let bsize = find_bsize(&extra)?;
    let block_size = bsize as usize + 1;
    let header_size = 12 + xlen;
    if block_size < header_size + 8 {
        return Err(CyaneaError::Parse("BGZF block too small".into()));
    }
    let data_size = block_size - header_size - 8;

    let mut cdata = vec![0u8; data_size];
    reader.read_exact(&mut cdata).map_err(CyaneaError::Io)?;

    let mut trailer = [0u8; 8];
    reader.read_exact(&mut trailer).map_err(CyaneaError::Io)?;
    let isize = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;

    if isize == 0 {
        return Ok(Some(Vec::new()));
    }

    let mut decompressed = vec![0u8; isize];
    let mut decoder = DeflateDecoder::new(&cdata[..]);
    decoder
        .read_exact(&mut decompressed)
        .map_err(|e| CyaneaError::Parse(format!("BGZF decompression failed: {e}")))?;

    Ok(Some(decompressed))
}

fn find_bsize(extra: &[u8]) -> Result<u16> {
    let mut pos = 0;
    while pos + 4 <= extra.len() {
        let si1 = extra[pos];
        let si2 = extra[pos + 1];
        let slen = u16::from_le_bytes([extra[pos + 2], extra[pos + 3]]) as usize;

        if si1 == b'B' && si2 == b'C' && slen == 2 && pos + 6 <= extra.len() {
            return Ok(u16::from_le_bytes([extra[pos + 4], extra[pos + 5]]));
        }
        pos += 4 + slen;
    }
    Err(CyaneaError::Parse(
        "BGZF extra field missing BSIZE subfield".into(),
    ))
}

// ---------------------------------------------------------------------------
// Binary reading helpers
// ---------------------------------------------------------------------------

fn read_u32_le(data: &[u8], pos: &mut usize) -> Result<u32> {
    if *pos + 4 > data.len() {
        return Err(CyaneaError::Parse("unexpected end of BCF data".into()));
    }
    let val = u32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
    *pos += 4;
    Ok(val)
}

fn read_i32_le(data: &[u8], pos: &mut usize) -> Result<i32> {
    if *pos + 4 > data.len() {
        return Err(CyaneaError::Parse("unexpected end of BCF data".into()));
    }
    let val = i32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
    *pos += 4;
    Ok(val)
}

fn read_f32_le(data: &[u8], pos: &mut usize) -> Result<f32> {
    if *pos + 4 > data.len() {
        return Err(CyaneaError::Parse("unexpected end of BCF data".into()));
    }
    let val = f32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
    *pos += 4;
    Ok(val)
}

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

fn parse_header_contigs(header: &str) -> Vec<String> {
    let mut contigs = Vec::new();
    for line in header.lines() {
        if line.starts_with("##contig=<") {
            if let Some(rest) = line.strip_prefix("##contig=<") {
                let rest = rest.trim_end_matches('>');
                for field in rest.split(',') {
                    if let Some(id) = field.strip_prefix("ID=") {
                        contigs.push(id.to_string());
                        break;
                    }
                }
            }
        }
    }
    contigs
}

fn parse_header_filters(header: &str) -> Vec<String> {
    let mut filters = Vec::new();
    for line in header.lines() {
        if line.starts_with("##FILTER=<") {
            if let Some(rest) = line.strip_prefix("##FILTER=<") {
                let rest = rest.trim_end_matches('>');
                for field in rest.split(',') {
                    if let Some(id) = field.strip_prefix("ID=") {
                        filters.push(id.to_string());
                        break;
                    }
                }
            }
        }
    }
    filters
}

// ---------------------------------------------------------------------------
// BCF record parsing
// ---------------------------------------------------------------------------

/// BCF2 typed value encoding.
///
/// Type byte: low 4 bits = type, high 4 bits = count (if < 15).
/// Types: 0=missing, 1=int8, 2=int16, 3=int32, 5=float, 7=char/string
fn read_typed_string(data: &[u8], pos: &mut usize) -> Result<String> {
    if *pos >= data.len() {
        return Ok(String::new());
    }

    let type_byte = data[*pos];
    *pos += 1;

    let type_id = type_byte & 0x0F;
    let mut count = ((type_byte >> 4) & 0x0F) as usize;

    // If count == 15, next byte(s) encode the actual count
    if count == 15 {
        if *pos >= data.len() {
            return Ok(String::new());
        }
        let count_type_byte = data[*pos];
        *pos += 1;
        let count_type = count_type_byte & 0x0F;
        count = match count_type {
            1 => {
                // int8
                if *pos >= data.len() {
                    return Ok(String::new());
                }
                let v = data[*pos] as usize;
                *pos += 1;
                v
            }
            2 => {
                // int16
                if *pos + 2 > data.len() {
                    return Ok(String::new());
                }
                let v = u16::from_le_bytes([data[*pos], data[*pos + 1]]) as usize;
                *pos += 2;
                v
            }
            3 => {
                // int32
                if *pos + 4 > data.len() {
                    return Ok(String::new());
                }
                let v = u32::from_le_bytes([
                    data[*pos],
                    data[*pos + 1],
                    data[*pos + 2],
                    data[*pos + 3],
                ]) as usize;
                *pos += 4;
                v
            }
            _ => 0,
        };
    }

    if type_id == 7 {
        // char/string
        if *pos + count > data.len() {
            count = data.len() - *pos;
        }
        let s = std::str::from_utf8(&data[*pos..*pos + count])
            .unwrap_or("")
            .trim_end_matches('\0')
            .to_string();
        *pos += count;
        Ok(s)
    } else {
        Ok(String::new())
    }
}

fn read_typed_int_vec(data: &[u8], pos: &mut usize) -> Result<Vec<i32>> {
    if *pos >= data.len() {
        return Ok(Vec::new());
    }

    let type_byte = data[*pos];
    *pos += 1;

    let type_id = type_byte & 0x0F;
    let mut count = ((type_byte >> 4) & 0x0F) as usize;

    if count == 15 {
        if *pos >= data.len() {
            return Ok(Vec::new());
        }
        let count_type_byte = data[*pos];
        *pos += 1;
        let count_type = count_type_byte & 0x0F;
        count = match count_type {
            1 => {
                if *pos >= data.len() { return Ok(Vec::new()); }
                let v = data[*pos] as usize;
                *pos += 1;
                v
            }
            2 => {
                if *pos + 2 > data.len() { return Ok(Vec::new()); }
                let v = u16::from_le_bytes([data[*pos], data[*pos + 1]]) as usize;
                *pos += 2;
                v
            }
            3 => {
                if *pos + 4 > data.len() { return Ok(Vec::new()); }
                let v = u32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]) as usize;
                *pos += 4;
                v
            }
            _ => 0,
        };
    }

    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        match type_id {
            1 => {
                if *pos >= data.len() { break; }
                values.push(data[*pos] as i8 as i32);
                *pos += 1;
            }
            2 => {
                if *pos + 2 > data.len() { break; }
                let v = i16::from_le_bytes([data[*pos], data[*pos + 1]]);
                values.push(v as i32);
                *pos += 2;
            }
            3 => {
                if *pos + 4 > data.len() { break; }
                let v = i32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
                values.push(v);
                *pos += 4;
            }
            _ => break,
        }
    }

    Ok(values)
}

fn parse_bcf_record(
    data: &[u8],
    pos: &mut usize,
    contigs: &[String],
    filters: &[String],
) -> Result<Variant> {
    // l_shared (u32) + l_indiv (u32)
    let l_shared = read_u32_le(data, pos)? as usize;
    let l_indiv = read_u32_le(data, pos)? as usize;

    let record_start = *pos;

    // CHROM (i32), POS (i32), rlen (i32), QUAL (f32)
    let chrom_idx = read_i32_le(data, pos)?;
    let position = read_i32_le(data, pos)?; // 0-based in BCF
    let _rlen = read_i32_le(data, pos)?;
    let qual = read_f32_le(data, pos)?;

    // n_info (u16) + n_allele (u16)
    if *pos + 4 > data.len() {
        return Err(CyaneaError::Parse("BCF record truncated".into()));
    }
    let n_info_allele = read_u32_le(data, pos)?;
    let _n_info = (n_info_allele & 0xFFFF) as u16;
    let n_allele = (n_info_allele >> 16) as u16;

    // n_fmt_sample (u32)
    let _n_fmt_sample = read_u32_le(data, pos)?;

    // ID (typed string)
    let id_str = read_typed_string(data, pos)?;
    let id = if id_str.is_empty() || id_str == "." {
        None
    } else {
        Some(id_str)
    };

    // Alleles: n_allele typed strings
    let mut alleles = Vec::with_capacity(n_allele as usize);
    for _ in 0..n_allele {
        let allele = read_typed_string(data, pos)?;
        alleles.push(allele);
    }

    // FILTER: typed int vector of filter indices
    let filter_indices = read_typed_int_vec(data, pos)?;

    // Skip remaining shared data (INFO fields) and individual data
    let consumed = *pos - record_start;
    let total = l_shared + l_indiv;
    if consumed < total {
        *pos = record_start + total;
    }

    // Build Variant
    let chrom = if chrom_idx >= 0 && (chrom_idx as usize) < contigs.len() {
        contigs[chrom_idx as usize].clone()
    } else {
        format!("chr{}", chrom_idx)
    };

    let ref_allele = alleles
        .first()
        .map(|a| a.as_bytes().to_vec())
        .unwrap_or_default();

    let alt_alleles: Vec<Vec<u8>> = alleles
        .iter()
        .skip(1)
        .filter(|a| !a.is_empty() && *a != ".")
        .map(|a| a.as_bytes().to_vec())
        .collect();

    let quality = if qual.is_nan() || qual == f32::from_bits(0x7F80_0001) {
        None
    } else {
        Some(qual as f64)
    };

    let filter = if filter_indices.is_empty() {
        VariantFilter::Missing
    } else if filter_indices == [0] {
        VariantFilter::Pass
    } else {
        let filter_names: Vec<String> = filter_indices
            .iter()
            .map(|&idx| {
                if idx >= 0 && (idx as usize) < filters.len() {
                    filters[idx as usize].clone()
                } else {
                    format!("FILTER{}", idx)
                }
            })
            .collect();
        if filter_names.iter().any(|f| f == "PASS") {
            VariantFilter::Pass
        } else {
            VariantFilter::Fail(filter_names)
        }
    };

    // BCF position is 0-based; convert to 1-based for consistency with VCF
    let vcf_position = if position < 0 { 0u64 } else { (position + 1) as u64 };

    Ok(Variant {
        chrom,
        position: vcf_position,
        id,
        ref_allele,
        alt_alleles,
        quality,
        filter,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test_helpers {
    use flate2::write::DeflateEncoder;
    use flate2::Compression;
    use std::io::Write;

    pub fn bgzf_compress(data: &[u8]) -> Vec<u8> {
        let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        let cdata = encoder.finish().unwrap();

        let cdata_len = cdata.len();
        let bsize = (18 + cdata_len + 8 - 1) as u16;

        let mut block = Vec::new();
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]);
        block.extend_from_slice(&[0u8; 6]);
        block.extend_from_slice(&6u16.to_le_bytes());
        block.push(b'B');
        block.push(b'C');
        block.extend_from_slice(&2u16.to_le_bytes());
        block.extend_from_slice(&bsize.to_le_bytes());
        block.extend_from_slice(&cdata);

        let crc = crc32(data);
        block.extend_from_slice(&crc.to_le_bytes());
        block.extend_from_slice(&(data.len() as u32).to_le_bytes());

        block
    }

    pub fn bgzf_eof_block() -> Vec<u8> {
        bgzf_compress(&[])
    }

    fn crc32(data: &[u8]) -> u32 {
        let mut crc = 0xFFFFFFFFu32;
        for &byte in data {
            crc ^= byte as u32;
            for _ in 0..8 {
                if crc & 1 != 0 {
                    crc = (crc >> 1) ^ 0xEDB88320;
                } else {
                    crc >>= 1;
                }
            }
        }
        !crc
    }

    /// Encode a BCF typed string.
    pub fn encode_typed_string(s: &str) -> Vec<u8> {
        let bytes = s.as_bytes();
        let mut out = Vec::new();
        let len = bytes.len();
        if len < 15 {
            out.push(((len as u8) << 4) | 7); // type=7 (char), count in high nibble
        } else {
            out.push((15 << 4) | 7); // count=15 means overflow
            if len < 128 {
                out.push((1 << 4) | 1); // count type: 1 value of int8
                out.push(len as u8);
            } else {
                out.push((1 << 4) | 2); // count type: 1 value of int16
                out.extend_from_slice(&(len as u16).to_le_bytes());
            }
        }
        out.extend_from_slice(bytes);
        out
    }

    /// Encode a BCF typed int8 vector (for filter indices).
    pub fn encode_typed_int8_vec(values: &[i8]) -> Vec<u8> {
        let mut out = Vec::new();
        let len = values.len();
        if len < 15 {
            out.push(((len as u8) << 4) | 1); // type=1 (int8)
        } else {
            out.push((15 << 4) | 1);
            out.push((1 << 4) | 1);
            out.push(len as u8);
        }
        for &v in values {
            out.push(v as u8);
        }
        out
    }

    /// Build a minimal BCF record (shared fields only, no indiv data).
    pub fn build_bcf_record(
        chrom_idx: i32,
        pos: i32, // 0-based
        qual: f32,
        id: &str,
        alleles: &[&str],
        filter_indices: &[i8],
    ) -> Vec<u8> {
        let mut shared = Vec::new();

        // CHROM, POS, RLEN, QUAL
        shared.extend_from_slice(&chrom_idx.to_le_bytes());
        shared.extend_from_slice(&pos.to_le_bytes());
        let rlen = alleles.first().map(|a| a.len() as i32).unwrap_or(0);
        shared.extend_from_slice(&rlen.to_le_bytes());
        shared.extend_from_slice(&qual.to_le_bytes());

        // n_info (low 16) + n_allele (high 16)
        let n_allele = alleles.len() as u32;
        let n_info_allele = (n_allele << 16) | 0;
        shared.extend_from_slice(&n_info_allele.to_le_bytes());

        // n_fmt_sample (0)
        shared.extend_from_slice(&0u32.to_le_bytes());

        // ID
        if id.is_empty() || id == "." {
            shared.push((0 << 4) | 7); // empty string
        } else {
            shared.extend_from_slice(&encode_typed_string(id));
        }

        // Alleles
        for allele in alleles {
            shared.extend_from_slice(&encode_typed_string(allele));
        }

        // FILTER
        shared.extend_from_slice(&encode_typed_int8_vec(filter_indices));

        // Wrap with l_shared + l_indiv
        let l_shared = shared.len() as u32;
        let l_indiv = 0u32;
        let mut record = Vec::new();
        record.extend_from_slice(&l_shared.to_le_bytes());
        record.extend_from_slice(&l_indiv.to_le_bytes());
        record.extend_from_slice(&shared);
        record
    }

    /// Build BCF binary content (magic + header + records).
    pub fn build_bcf_content(
        contigs: &[&str],
        filter_defs: &[&str],
        records_data: &[Vec<u8>],
    ) -> Vec<u8> {
        let mut data = Vec::new();

        // BCF magic
        data.extend_from_slice(b"BCF\x02\x01");

        // Build VCF header text
        let mut header = String::new();
        header.push_str("##fileformat=VCFv4.3\n");
        for (i, contig) in contigs.iter().enumerate() {
            header.push_str(&format!("##contig=<ID={},length={}>\n", contig, 1000 * (i + 1)));
        }
        for filter_name in filter_defs {
            header.push_str(&format!("##FILTER=<ID={},Description=\"{}\">\n", filter_name, filter_name));
        }
        header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        // Add NUL terminator
        let header_bytes = header.as_bytes();
        let header_len = header_bytes.len() + 1; // +1 for NUL
        data.extend_from_slice(&(header_len as u32).to_le_bytes());
        data.extend_from_slice(header_bytes);
        data.push(0); // NUL terminator

        // Records
        for rec in records_data {
            data.extend_from_slice(rec);
        }

        data
    }

    pub fn write_test_bcf(
        contigs: &[&str],
        filter_defs: &[&str],
        records_data: &[Vec<u8>],
    ) -> tempfile::NamedTempFile {
        let content = build_bcf_content(contigs, filter_defs, records_data);
        let compressed = bgzf_compress(&content);
        let eof = bgzf_eof_block();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&compressed).unwrap();
        file.write_all(&eof).unwrap();
        file.flush().unwrap();
        file
    }
}

#[cfg(test)]
mod tests {
    use super::test_helpers::*;
    use super::*;

    #[test]
    fn bcf_parse_simple() {
        let records = vec![
            build_bcf_record(0, 99, 30.0, "rs123", &["A", "G"], &[0]),     // PASS
            build_bcf_record(0, 199, f32::NAN, ".", &["AC", "A"], &[]),     // Missing filter
            build_bcf_record(1, 299, 50.5, ".", &["T", "TA", "TG"], &[1]), // LowQual
        ];

        let file = write_test_bcf(
            &["chr1", "chr2"],
            &["PASS", "LowQual"],
            &records,
        );

        let variants = parse_bcf(file.path()).unwrap();
        assert_eq!(variants.len(), 3);

        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].position, 100); // 0-based 99 → 1-based 100
        assert_eq!(variants[0].id, Some("rs123".to_string()));
        assert_eq!(variants[0].ref_allele, b"A");
        assert_eq!(variants[0].alt_alleles, vec![b"G".to_vec()]);
        assert!((variants[0].quality.unwrap() - 30.0).abs() < 0.1);
        assert_eq!(variants[0].filter, VariantFilter::Pass);

        assert_eq!(variants[1].position, 200);
        assert_eq!(variants[1].id, None);
        assert_eq!(variants[1].quality, None);

        assert_eq!(variants[2].chrom, "chr2");
        assert_eq!(variants[2].alt_alleles.len(), 2);
    }

    #[test]
    fn bcf_matches_vcf() {
        // Verify BCF and VCF produce compatible Variant records
        let records = vec![
            build_bcf_record(0, 99, 30.0, ".", &["A", "G"], &[0]),
        ];
        let file = write_test_bcf(&["chr1"], &["PASS"], &records);
        let variants = parse_bcf(file.path()).unwrap();

        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].position, 100);
        assert_eq!(variants[0].ref_allele, b"A");
        assert_eq!(variants[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(variants[0].filter, VariantFilter::Pass);
    }

    #[test]
    fn bcf_stats_computed() {
        let records = vec![
            build_bcf_record(0, 99, 30.0, ".", &["A", "G"], &[0]),
            build_bcf_record(0, 199, 40.0, ".", &["AC", "A"], &[0]),
            build_bcf_record(1, 299, 50.0, ".", &["T", "C"], &[1]),
        ];
        let file = write_test_bcf(&["chr1", "chr2"], &["PASS", "LowQual"], &records);

        let stats = bcf_stats(file.path()).unwrap();
        assert_eq!(stats.variant_count, 3);
        assert_eq!(stats.snv_count, 2);
        assert_eq!(stats.indel_count, 1);
        assert_eq!(stats.pass_count, 2);
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn bcf_invalid_magic() {
        let mut bad_data = vec![0u8; 100];
        bad_data[0..5].copy_from_slice(b"XXXXX");
        let compressed = bgzf_compress(&bad_data);
        let eof = bgzf_eof_block();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        use std::io::Write;
        file.write_all(&compressed).unwrap();
        file.write_all(&eof).unwrap();
        file.flush().unwrap();

        let result = parse_bcf(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("not a valid BCF"));
    }

    #[test]
    fn bcf_multiallelic() {
        let records = vec![
            build_bcf_record(0, 99, 30.0, ".", &["A", "G", "T", "C"], &[0]),
        ];
        let file = write_test_bcf(&["chr1"], &["PASS"], &records);
        let variants = parse_bcf(file.path()).unwrap();

        assert_eq!(variants[0].alt_alleles.len(), 3);
        assert_eq!(variants[0].alt_alleles[0], b"G");
        assert_eq!(variants[0].alt_alleles[1], b"T");
        assert_eq!(variants[0].alt_alleles[2], b"C");
    }
}
