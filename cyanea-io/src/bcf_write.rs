//! BCF2 binary writer.
//!
//! Writes [`Variant`] records in BCF2 format with BGZF compression.
//! Reuses BCF typed value encoding patterns from the reader module.

use std::io::Write;
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::variant::{Variant, VariantFilter};

use crate::vcf_header::VcfHeader;

/// BCF2 magic: "BCF\2\1"
const BCF_MAGIC: [u8; 5] = [b'B', b'C', b'F', 0x02, 0x01];

/// Write variants in BCF2 format to a file.
///
/// The header defines contig and filter ordering for integer encoding.
pub fn write_bcf(header: &VcfHeader, variants: &[Variant], path: impl AsRef<Path>) -> Result<()> {
    let data = write_bcf_bytes(header, variants)?;
    let path = path.as_ref();
    let mut file = std::fs::File::create(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    file.write_all(&data).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    Ok(())
}

/// Write variants in BCF2 format to a byte vector (BGZF-compressed).
pub fn write_bcf_bytes(header: &VcfHeader, variants: &[Variant]) -> Result<Vec<u8>> {
    let mut raw = Vec::new();

    // BCF magic
    raw.extend_from_slice(&BCF_MAGIC);

    // Build contig and filter lookup maps
    let contig_map: std::collections::HashMap<&str, i32> = header
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.id.as_str(), i as i32))
        .collect();

    let filter_map: std::collections::HashMap<&str, i32> = {
        let mut m = std::collections::HashMap::new();
        // PASS is always filter index 0 in BCF
        m.insert("PASS", 0i32);
        for (i, f) in header.filter_fields.iter().enumerate() {
            m.insert(f.id.as_str(), (i + 1) as i32);
        }
        m
    };

    // Build header text — ensure PASS is the first FILTER definition (BCF spec: index 0)
    let header_with_pass = {
        let mut h = header.clone();
        // Insert PASS at the front if not already present
        if !h.filter_fields.iter().any(|f| f.id == "PASS") {
            h.filter_fields.insert(
                0,
                crate::vcf_header::FilterDef {
                    id: "PASS".to_string(),
                    description: "All filters passed".to_string(),
                },
            );
        }
        h
    };
    let header_text = header_with_pass.to_vcf_string();
    let header_len = (header_text.len() + 1) as u32; // +1 for NUL
    raw.extend_from_slice(&header_len.to_le_bytes());
    raw.extend_from_slice(header_text.as_bytes());
    raw.push(0); // NUL terminator

    // Write variant records
    for v in variants {
        let record = encode_bcf_record(v, &contig_map, &filter_map)?;
        raw.extend_from_slice(&record);
    }

    // BGZF-compress the raw data
    bgzf_compress_all(&raw)
}

/// Encode a single variant as a BCF2 record.
fn encode_bcf_record(
    v: &Variant,
    contig_map: &std::collections::HashMap<&str, i32>,
    filter_map: &std::collections::HashMap<&str, i32>,
) -> Result<Vec<u8>> {
    let mut shared = Vec::new();

    // CHROM index
    let chrom_idx = contig_map.get(v.chrom.as_str()).copied().unwrap_or(-1);
    shared.extend_from_slice(&chrom_idx.to_le_bytes());

    // POS (0-based in BCF, our Variant uses 1-based)
    let pos = if v.position > 0 {
        (v.position - 1) as i32
    } else {
        0i32
    };
    shared.extend_from_slice(&pos.to_le_bytes());

    // RLEN (reference allele length)
    let rlen = v.ref_allele.len() as i32;
    shared.extend_from_slice(&rlen.to_le_bytes());

    // QUAL
    let qual = v
        .quality
        .map(|q| q as f32)
        .unwrap_or(f32::from_bits(0x7F80_0001));
    shared.extend_from_slice(&qual.to_le_bytes());

    // n_info (low 16 bits) + n_allele (high 16 bits)
    let n_allele = (1 + v.alt_alleles.len()) as u32;
    let n_info_allele = (n_allele << 16) | 0;
    shared.extend_from_slice(&n_info_allele.to_le_bytes());

    // n_fmt_sample (0 — no genotype data)
    shared.extend_from_slice(&0u32.to_le_bytes());

    // ID
    if let Some(ref id) = v.id {
        encode_typed_string(&mut shared, id);
    } else {
        shared.push((0 << 4) | 7); // empty string
    }

    // Alleles: ref first, then alts
    let ref_str = std::str::from_utf8(&v.ref_allele).unwrap_or("N");
    encode_typed_string(&mut shared, ref_str);
    for alt in &v.alt_alleles {
        let alt_str = std::str::from_utf8(alt).unwrap_or(".");
        encode_typed_string(&mut shared, alt_str);
    }

    // FILTER
    match &v.filter {
        VariantFilter::Pass => {
            encode_typed_int8_vec(&mut shared, &[0]); // PASS = index 0
        }
        VariantFilter::Missing => {
            shared.push((0 << 4) | 1); // empty int8 vec
        }
        VariantFilter::Fail(reasons) => {
            let indices: Vec<i8> = reasons
                .iter()
                .map(|r| filter_map.get(r.as_str()).copied().unwrap_or(0) as i8)
                .collect();
            encode_typed_int8_vec(&mut shared, &indices);
        }
    }

    // Wrap with l_shared + l_indiv
    let l_shared = shared.len() as u32;
    let l_indiv = 0u32;
    let mut record = Vec::new();
    record.extend_from_slice(&l_shared.to_le_bytes());
    record.extend_from_slice(&l_indiv.to_le_bytes());
    record.extend_from_slice(&shared);

    Ok(record)
}

/// Encode a BCF typed string.
fn encode_typed_string(buf: &mut Vec<u8>, s: &str) {
    let bytes = s.as_bytes();
    let len = bytes.len();
    if len < 15 {
        buf.push(((len as u8) << 4) | 7);
    } else if len < 128 {
        buf.push((15 << 4) | 7);
        buf.push((1 << 4) | 1);
        buf.push(len as u8);
    } else {
        buf.push((15 << 4) | 7);
        buf.push((1 << 4) | 2);
        buf.extend_from_slice(&(len as u16).to_le_bytes());
    }
    buf.extend_from_slice(bytes);
}

/// Encode a BCF typed int8 vector.
fn encode_typed_int8_vec(buf: &mut Vec<u8>, values: &[i8]) {
    let len = values.len();
    if len < 15 {
        buf.push(((len as u8) << 4) | 1);
    } else {
        buf.push((15 << 4) | 1);
        buf.push((1 << 4) | 1);
        buf.push(len as u8);
    }
    for &v in values {
        buf.push(v as u8);
    }
}

// ---------------------------------------------------------------------------
// BGZF compression
// ---------------------------------------------------------------------------

/// BGZF-compress data into blocks + EOF marker.
fn bgzf_compress_all(data: &[u8]) -> Result<Vec<u8>> {
    use flate2::write::DeflateEncoder;
    use flate2::Compression;

    let mut output = Vec::new();
    const MAX_BLOCK: usize = 60 * 1024;

    let mut offset = 0;
    while offset < data.len() {
        let end = (offset + MAX_BLOCK).min(data.len());
        let chunk = &data[offset..end];

        let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(chunk).map_err(|e| {
            CyaneaError::Io(std::io::Error::new(e.kind(), format!("BGZF compress: {e}")))
        })?;
        let cdata = encoder.finish().map_err(|e| {
            CyaneaError::Io(std::io::Error::new(e.kind(), format!("BGZF finish: {e}")))
        })?;

        write_bgzf_block(&mut output, &cdata, chunk);
        offset = end;
    }

    // EOF block
    let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(&[]).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(e.kind(), format!("BGZF EOF: {e}")))
    })?;
    let eof_cdata = encoder.finish().map_err(|e| {
        CyaneaError::Io(std::io::Error::new(e.kind(), format!("BGZF EOF finish: {e}")))
    })?;
    write_bgzf_block(&mut output, &eof_cdata, &[]);

    Ok(output)
}

/// Write a single BGZF block with proper CRC32.
fn write_bgzf_block(output: &mut Vec<u8>, cdata: &[u8], original: &[u8]) {
    let bsize = (18 + cdata.len() + 8 - 1) as u16;

    // Gzip header with FEXTRA
    output.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]);
    output.extend_from_slice(&[0u8; 6]); // MTIME, XFL, OS
    output.extend_from_slice(&6u16.to_le_bytes()); // XLEN = 6

    // BGZF extra subfield
    output.push(b'B');
    output.push(b'C');
    output.extend_from_slice(&2u16.to_le_bytes()); // SLEN = 2
    output.extend_from_slice(&bsize.to_le_bytes()); // BSIZE

    // Compressed data
    output.extend_from_slice(cdata);

    // CRC32 of original uncompressed data + ISIZE
    let crc = crc32(original);
    output.extend_from_slice(&crc.to_le_bytes());
    output.extend_from_slice(&(original.len() as u32).to_le_bytes());
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::parse_bcf;

    fn snv(chrom: &str, pos: u64, ref_a: &str, alt: &str, qual: Option<f64>) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            position: pos,
            id: None,
            ref_allele: ref_a.as_bytes().to_vec(),
            alt_alleles: vec![alt.as_bytes().to_vec()],
            quality: qual,
            filter: VariantFilter::Pass,
        }
    }

    #[test]
    fn bcf_write_roundtrip_single() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));

        let variants = vec![snv("chr1", 100, "A", "G", Some(30.0))];

        let bytes = write_bcf_bytes(&header, &variants).unwrap();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let parsed = parse_bcf(file.path()).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].chrom, "chr1");
        assert_eq!(parsed[0].position, 100);
        assert_eq!(parsed[0].ref_allele, b"A");
        assert_eq!(parsed[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(parsed[0].filter, VariantFilter::Pass);
    }

    #[test]
    fn bcf_write_roundtrip_multiple() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));
        header.add_contig("chr2", Some(2000));
        header.add_filter("LowQual", "Low quality");

        let variants = vec![
            snv("chr1", 100, "A", "G", Some(30.0)),
            Variant {
                chrom: "chr1".to_string(),
                position: 200,
                id: Some("rs456".to_string()),
                ref_allele: b"C".to_vec(),
                alt_alleles: vec![b"T".to_vec()],
                quality: Some(50.0),
                filter: VariantFilter::Pass,
            },
            Variant {
                chrom: "chr2".to_string(),
                position: 300,
                id: None,
                ref_allele: b"G".to_vec(),
                alt_alleles: vec![b"A".to_vec()],
                quality: None,
                filter: VariantFilter::Missing,
            },
        ];

        let bytes = write_bcf_bytes(&header, &variants).unwrap();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let parsed = parse_bcf(file.path()).unwrap();
        assert_eq!(parsed.len(), 3);

        assert_eq!(parsed[0].chrom, "chr1");
        assert_eq!(parsed[0].position, 100);

        assert_eq!(parsed[1].chrom, "chr1");
        assert_eq!(parsed[1].position, 200);
        assert_eq!(parsed[1].id, Some("rs456".to_string()));

        assert_eq!(parsed[2].chrom, "chr2");
        assert_eq!(parsed[2].position, 300);
        assert_eq!(parsed[2].quality, None);
    }

    #[test]
    fn bcf_write_multiallelic() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));

        let v = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: None,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec(), b"T".to_vec()],
            quality: Some(40.0),
            filter: VariantFilter::Pass,
        };

        let bytes = write_bcf_bytes(&header, &[v]).unwrap();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let parsed = parse_bcf(file.path()).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].alt_alleles.len(), 2);
        assert_eq!(parsed[0].alt_alleles[0], b"G");
        assert_eq!(parsed[0].alt_alleles[1], b"T");
    }

    #[test]
    fn bcf_write_filter_fail() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));
        header.add_filter("LowQual", "Low quality");

        let v = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: None,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            quality: Some(5.0),
            filter: VariantFilter::Fail(vec!["LowQual".to_string()]),
        };

        let bytes = write_bcf_bytes(&header, &[v]).unwrap();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let parsed = parse_bcf(file.path()).unwrap();
        assert_eq!(parsed.len(), 1);
        match &parsed[0].filter {
            VariantFilter::Fail(reasons) => {
                assert!(reasons.contains(&"LowQual".to_string()));
            }
            other => panic!("expected Fail filter, got {:?}", other),
        }
    }

    #[test]
    fn bcf_write_empty() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));

        let bytes = write_bcf_bytes(&header, &[]).unwrap();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let parsed = parse_bcf(file.path()).unwrap();
        assert!(parsed.is_empty());
    }

    #[test]
    fn bcf_write_with_id() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));

        let v = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: Some("rs12345".to_string()),
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            quality: Some(30.0),
            filter: VariantFilter::Pass,
        };

        let bytes = write_bcf_bytes(&header, &[v]).unwrap();

        let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let parsed = parse_bcf(file.path()).unwrap();
        assert_eq!(parsed[0].id, Some("rs12345".to_string()));
    }

    #[test]
    fn bcf_write_to_file_path() {
        let mut header = VcfHeader::new();
        header.add_contig("chr1", Some(1000));

        let v = snv("chr1", 100, "A", "G", Some(30.0));

        let tmp = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
        write_bcf(&header, &[v], tmp.path()).unwrap();

        let parsed = parse_bcf(tmp.path()).unwrap();
        assert_eq!(parsed.len(), 1);
    }
}
