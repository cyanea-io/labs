//! Shared BGZF (Blocked GNU Zip Format) utilities.
//!
//! BGZF is a series of concatenated gzip blocks, each at most 64 KiB uncompressed,
//! enabling random access when paired with an index (BAI, TBI, CSI).
//!
//! This module provides common BGZF reading functions used by both the BAM and
//! BCF parsers.

use std::io::{self, Read};

use cyanea_core::{CyaneaError, Result};
use flate2::read::DeflateDecoder;

/// BGZF magic bytes: standard gzip header with FEXTRA flag.
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// A BGZF virtual file offset: block offset (upper 48 bits) + within-block offset (lower 16 bits).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualOffset(pub u64);

impl VirtualOffset {
    /// Create a virtual offset from block offset and within-block offset.
    pub fn new(block_offset: u64, within_block_offset: u16) -> Self {
        Self((block_offset << 16) | within_block_offset as u64)
    }

    /// The compressed file offset of the BGZF block.
    pub fn block_offset(&self) -> u64 {
        self.0 >> 16
    }

    /// The offset within the uncompressed BGZF block.
    pub fn within_block_offset(&self) -> u16 {
        (self.0 & 0xFFFF) as u16
    }
}

/// Read and decompress the next BGZF block from a reader.
///
/// Returns `Ok(None)` at EOF, `Ok(Some(data))` for a valid block.
/// An empty `Vec` signals the BGZF EOF marker block.
pub fn read_bgzf_block(reader: &mut impl Read) -> Result<Option<Vec<u8>>> {
    // Read gzip header (at least 18 bytes for basic header + BGZF extra)
    let mut header = [0u8; 18];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(CyaneaError::Io(e)),
    }

    // Validate gzip magic + FEXTRA flag
    if header[0] != BGZF_MAGIC[0]
        || header[1] != BGZF_MAGIC[1]
        || header[2] != BGZF_MAGIC[2]
        || (header[3] & 0x04) == 0
    {
        return Err(CyaneaError::Parse("not a valid BGZF block".into()));
    }

    // XLEN is at offset 10-11 (little-endian)
    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;

    // We already read 18 bytes total (12 gzip header + 6 of extra field)
    // The extra field starts at offset 12. We've read 6 bytes of it (offsets 12-17).
    let mut extra = vec![0u8; xlen];
    extra[..6].copy_from_slice(&header[12..18]);
    if xlen > 6 {
        reader
            .read_exact(&mut extra[6..])
            .map_err(CyaneaError::Io)?;
    }

    // Find BSIZE in extra field subfields
    let bsize = find_bsize(&extra)?;

    // Total block size = BSIZE + 1 (BSIZE is 0-based)
    let block_size = bsize as usize + 1;

    // Bytes remaining after the gzip header (12 bytes) + extra field (xlen bytes):
    let header_size = 12 + xlen;
    if block_size < header_size + 8 {
        return Err(CyaneaError::Parse("BGZF block too small".into()));
    }
    let data_size = block_size - header_size - 8; // subtract CRC32(4) + ISIZE(4)

    // Read compressed data
    let mut cdata = vec![0u8; data_size];
    reader.read_exact(&mut cdata).map_err(CyaneaError::Io)?;

    // Read CRC32 and ISIZE
    let mut trailer = [0u8; 8];
    reader.read_exact(&mut trailer).map_err(CyaneaError::Io)?;
    let isize = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;

    // Empty block (EOF marker)
    if isize == 0 {
        return Ok(Some(Vec::new()));
    }

    // Decompress using raw DEFLATE
    let mut decompressed = vec![0u8; isize];
    let mut decoder = DeflateDecoder::new(&cdata[..]);
    decoder
        .read_exact(&mut decompressed)
        .map_err(|e| CyaneaError::Parse(format!("BGZF decompression failed: {e}")))?;

    Ok(Some(decompressed))
}

/// Find the BSIZE field in BGZF extra data.
///
/// BGZF stores BSIZE in a subfield with SI1=66 ('B'), SI2=67 ('C'), SLEN=2.
pub fn find_bsize(extra: &[u8]) -> Result<u16> {
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

/// Decompress all BGZF blocks from a reader into a contiguous buffer.
pub fn decompress_all(reader: &mut impl Read) -> Result<Vec<u8>> {
    let mut data = Vec::new();
    loop {
        match read_bgzf_block(reader)? {
            None => break,
            Some(block) => {
                if block.is_empty() {
                    break; // EOF marker
                }
                data.extend_from_slice(&block);
            }
        }
    }
    Ok(data)
}

// ---------------------------------------------------------------------------
// Test helpers (shared between BAM and BCF tests)
// ---------------------------------------------------------------------------

#[cfg(test)]
pub(crate) mod test_helpers {
    use flate2::write::DeflateEncoder;
    use flate2::Compression;
    use std::io::Write;

    /// Create a BGZF-compressed block from uncompressed data.
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

    /// Create a BGZF EOF marker block.
    pub fn bgzf_eof_block() -> Vec<u8> {
        bgzf_compress(&[])
    }

    /// Simple CRC32 computation.
    pub fn crc32(data: &[u8]) -> u32 {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_helpers::*;

    #[test]
    fn virtual_offset_roundtrip() {
        let vo = VirtualOffset::new(12345, 678);
        assert_eq!(vo.block_offset(), 12345);
        assert_eq!(vo.within_block_offset(), 678);
    }

    #[test]
    fn bgzf_block_roundtrip() {
        let original = b"Hello, BGZF world!";
        let compressed = bgzf_compress(original);
        let mut reader = std::io::Cursor::new(compressed);
        let decompressed = read_bgzf_block(&mut reader).unwrap().unwrap();
        assert_eq!(decompressed, original);
    }

    #[test]
    fn bgzf_eof_detection() {
        let eof = bgzf_eof_block();
        let mut reader = std::io::Cursor::new(eof);
        let block = read_bgzf_block(&mut reader).unwrap().unwrap();
        assert!(block.is_empty());
    }

    #[test]
    fn decompress_all_multiple_blocks() {
        let block1 = bgzf_compress(b"block1");
        let block2 = bgzf_compress(b"block2");
        let eof = bgzf_eof_block();

        let mut data = Vec::new();
        data.extend_from_slice(&block1);
        data.extend_from_slice(&block2);
        data.extend_from_slice(&eof);

        let mut reader = std::io::Cursor::new(data);
        let result = decompress_all(&mut reader).unwrap();
        assert_eq!(result, b"block1block2");
    }
}
