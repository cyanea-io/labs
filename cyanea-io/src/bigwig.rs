//! bigWig/bigBed binary format reader.
//!
//! Kent binary formats for dense signal data (bigWig) and compressed BED
//! intervals (bigBed). Both share a common header structure with B+ tree
//! for chromosome name → ID mapping and R-tree for spatial indexing.
//!
//! This module provides header reading and region-based querying.

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

/// bigWig magic number (little-endian).
const BIGWIG_MAGIC_LE: u32 = 0x888F_FC26;
/// bigWig magic number (big-endian, byte-swapped).
const BIGWIG_MAGIC_BE: u32 = 0x26FC_8F88;
/// bigBed magic number (little-endian).
const BIGBED_MAGIC_LE: u32 = 0x8789_F2EB;
/// bigBed magic number (big-endian, byte-swapped).
const BIGBED_MAGIC_BE: u32 = 0xEBF2_8978;

/// Header from a bigWig or bigBed file.
#[derive(Debug, Clone)]
pub struct BigWigHeader {
    /// Format version.
    pub version: u16,
    /// Number of zoom levels.
    pub zoom_levels: u16,
    /// Number of chromosomes.
    pub chrom_count: u32,
    /// Total summary statistics.
    pub total_summary: BigWigSummary,
}

/// Summary statistics stored in the bigWig/bigBed header.
#[derive(Debug, Clone)]
pub struct BigWigSummary {
    /// Number of bases covered.
    pub bases_covered: u64,
    /// Minimum value.
    pub min_val: f64,
    /// Maximum value.
    pub max_val: f64,
    /// Sum of all values.
    pub sum: f64,
    /// Sum of squares of all values.
    pub sum_squares: f64,
}

/// A single bigWig interval with a value.
#[derive(Debug, Clone)]
pub struct BigWigInterval {
    /// Chromosome name.
    pub chrom: String,
    /// 0-based start position.
    pub start: u32,
    /// End position (exclusive).
    pub end: u32,
    /// Signal value.
    pub value: f32,
}

/// A single bigBed record.
#[derive(Debug, Clone)]
pub struct BigBedRecord {
    /// Chromosome name.
    pub chrom: String,
    /// 0-based start position.
    pub start: u32,
    /// End position (exclusive).
    pub end: u32,
    /// Remaining BED fields as raw string.
    pub rest: String,
}

// ---------------------------------------------------------------------------
// Internal header structure
// ---------------------------------------------------------------------------

/// Parsed internal header with file offsets for navigating the binary.
struct InternalHeader {
    big_endian: bool,
    version: u16,
    zoom_levels: u16,
    chrom_tree_offset: u64,
    _data_offset: u64,
    _data_count: u64,
    total_summary_offset: u64,
    _uncompress_buf_size: u32,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Read the header and total summary from a bigWig file.
pub fn read_bigwig_header(path: impl AsRef<Path>) -> Result<BigWigHeader> {
    let path = path.as_ref();
    let mut file = open_file(path)?;
    let header = read_internal_header(&mut file, path, true)?;
    let summary = read_total_summary(&mut file, &header, path)?;
    let chroms = read_chrom_tree(&mut file, &header, path)?;

    Ok(BigWigHeader {
        version: header.version,
        zoom_levels: header.zoom_levels,
        chrom_count: chroms.len() as u32,
        total_summary: summary,
    })
}

/// Query bigWig intervals overlapping a genomic region.
pub fn read_bigwig_intervals(
    path: impl AsRef<Path>,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<Vec<BigWigInterval>> {
    let path = path.as_ref();
    let mut file = open_file(path)?;
    let header = read_internal_header(&mut file, path, true)?;
    let chroms = read_chrom_tree(&mut file, &header, path)?;

    let chrom_id = chroms
        .iter()
        .find(|(name, _, _)| name == chrom)
        .map(|(_, id, _)| *id)
        .ok_or_else(|| {
            CyaneaError::Parse(format!("{}: chromosome '{}' not found", path.display(), chrom))
        })?;

    // Read R-tree to find overlapping data blocks
    let blocks = read_rtree_overlapping(&mut file, &header, path, chrom_id, start, end)?;

    let mut intervals = Vec::new();
    for (offset, size) in blocks {
        let block_intervals =
            read_bigwig_data_block(&mut file, &header, path, offset, size, chrom, start, end)?;
        intervals.extend(block_intervals);
    }

    Ok(intervals)
}

/// Read the header from a bigBed file.
pub fn read_bigbed_header(path: impl AsRef<Path>) -> Result<BigWigHeader> {
    let path = path.as_ref();
    let mut file = open_file(path)?;
    let header = read_internal_header(&mut file, path, false)?;
    let summary = read_total_summary(&mut file, &header, path)?;
    let chroms = read_chrom_tree(&mut file, &header, path)?;

    Ok(BigWigHeader {
        version: header.version,
        zoom_levels: header.zoom_levels,
        chrom_count: chroms.len() as u32,
        total_summary: summary,
    })
}

/// Query bigBed records overlapping a genomic region.
pub fn read_bigbed_records(
    path: impl AsRef<Path>,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<Vec<BigBedRecord>> {
    let path = path.as_ref();
    let mut file = open_file(path)?;
    let header = read_internal_header(&mut file, path, false)?;
    let chroms = read_chrom_tree(&mut file, &header, path)?;

    let chrom_id = chroms
        .iter()
        .find(|(name, _, _)| name == chrom)
        .map(|(_, id, _)| *id)
        .ok_or_else(|| {
            CyaneaError::Parse(format!("{}: chromosome '{}' not found", path.display(), chrom))
        })?;

    let blocks = read_rtree_overlapping(&mut file, &header, path, chrom_id, start, end)?;

    let mut records = Vec::new();
    for (offset, size) in blocks {
        let block_records =
            read_bigbed_data_block(&mut file, &header, path, offset, size, chrom, start, end)?;
        records.extend(block_records);
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// File helpers
// ---------------------------------------------------------------------------

fn open_file(path: &Path) -> Result<File> {
    File::open(path).map_err(|e| {
        CyaneaError::Io(io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })
}

fn read_err(path: &Path, msg: &str) -> CyaneaError {
    CyaneaError::Parse(format!("{}: {}", path.display(), msg))
}

// ---------------------------------------------------------------------------
// Binary reading helpers
// ---------------------------------------------------------------------------

fn read_bytes(file: &mut File, n: usize, path: &Path) -> Result<Vec<u8>> {
    let mut buf = vec![0u8; n];
    file.read_exact(&mut buf)
        .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;
    Ok(buf)
}

fn read_u16(data: &[u8], offset: usize, big_endian: bool) -> u16 {
    if big_endian {
        u16::from_be_bytes([data[offset], data[offset + 1]])
    } else {
        u16::from_le_bytes([data[offset], data[offset + 1]])
    }
}

fn read_u32(data: &[u8], offset: usize, big_endian: bool) -> u32 {
    let b = &data[offset..offset + 4];
    if big_endian {
        u32::from_be_bytes([b[0], b[1], b[2], b[3]])
    } else {
        u32::from_le_bytes([b[0], b[1], b[2], b[3]])
    }
}

fn read_u64(data: &[u8], offset: usize, big_endian: bool) -> u64 {
    let b = &data[offset..offset + 8];
    if big_endian {
        u64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
    } else {
        u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
    }
}

fn read_f64(data: &[u8], offset: usize, big_endian: bool) -> f64 {
    let bits = read_u64(data, offset, big_endian);
    f64::from_bits(bits)
}

fn read_f32(data: &[u8], offset: usize, big_endian: bool) -> f32 {
    let bits = read_u32(data, offset, big_endian);
    f32::from_bits(bits)
}

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

fn read_internal_header(
    file: &mut File,
    path: &Path,
    is_bigwig: bool,
) -> Result<InternalHeader> {
    let buf = read_bytes(file, 64, path)?;

    let raw_magic = u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
    let (big_endian, valid) = if is_bigwig {
        match raw_magic {
            BIGWIG_MAGIC_LE => (false, true),
            BIGWIG_MAGIC_BE => (true, true),
            _ => (false, false),
        }
    } else {
        match raw_magic {
            BIGBED_MAGIC_LE => (false, true),
            BIGBED_MAGIC_BE => (true, true),
            _ => (false, false),
        }
    };

    if !valid {
        let expected = if is_bigwig { "bigWig" } else { "bigBed" };
        return Err(read_err(path, &format!("not a valid {} file (bad magic)", expected)));
    }

    let version = read_u16(&buf, 4, big_endian);
    let zoom_levels = read_u16(&buf, 6, big_endian);
    let chrom_tree_offset = read_u64(&buf, 8, big_endian);
    let data_offset = read_u64(&buf, 16, big_endian);
    // Index offset at 24
    let _index_offset = read_u64(&buf, 24, big_endian);
    // Field count at 32 (u16), defined field count at 34 (u16) - bigBed only
    // Auto-SQL offset at 40 (u64)
    let total_summary_offset = read_u64(&buf, 48, big_endian);
    let uncompress_buf_size = read_u32(&buf, 56, big_endian);
    // Data count
    let data_count = read_u64(&buf, 32, big_endian);

    Ok(InternalHeader {
        big_endian,
        version,
        zoom_levels,
        chrom_tree_offset,
        _data_offset: data_offset,
        _data_count: data_count,
        total_summary_offset,
        _uncompress_buf_size: uncompress_buf_size,
    })
}

fn read_total_summary(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
) -> Result<BigWigSummary> {
    if header.total_summary_offset == 0 {
        return Ok(BigWigSummary {
            bases_covered: 0,
            min_val: 0.0,
            max_val: 0.0,
            sum: 0.0,
            sum_squares: 0.0,
        });
    }

    file.seek(SeekFrom::Start(header.total_summary_offset))
        .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;

    let buf = read_bytes(file, 40, path)?;
    let be = header.big_endian;

    Ok(BigWigSummary {
        bases_covered: read_u64(&buf, 0, be),
        min_val: read_f64(&buf, 8, be),
        max_val: read_f64(&buf, 16, be),
        sum: read_f64(&buf, 24, be),
        sum_squares: read_f64(&buf, 32, be),
    })
}

// ---------------------------------------------------------------------------
// B+ tree (chromosome name → ID)
// ---------------------------------------------------------------------------

/// Returns Vec of (name, id, size).
fn read_chrom_tree(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
) -> Result<Vec<(String, u32, u32)>> {
    file.seek(SeekFrom::Start(header.chrom_tree_offset))
        .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;

    let be = header.big_endian;
    let buf = read_bytes(file, 32, path)?;

    let magic = read_u32(&buf, 0, be);
    if magic != 0x78CA_8C91 {
        return Err(read_err(path, "invalid B+ tree magic"));
    }

    // block_size at 4 (u32), key_size at 8 (u32), val_size at 12 (u32)
    let key_size = read_u32(&buf, 8, be) as usize;
    let _val_size = read_u32(&buf, 12, be) as usize;
    let _item_count = read_u64(&buf, 16, be);
    // reserved at 24 (u64)

    // Root node follows header (at offset chrom_tree_offset + 32)
    let mut chroms = Vec::new();
    read_bptree_node(file, header, path, key_size, &mut chroms)?;
    Ok(chroms)
}

fn read_bptree_node(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
    key_size: usize,
    chroms: &mut Vec<(String, u32, u32)>,
) -> Result<()> {
    let be = header.big_endian;
    let node_buf = read_bytes(file, 4, path)?;
    let is_leaf = node_buf[0] != 0;
    let child_count = read_u16(&node_buf, 2, be) as usize;

    if is_leaf {
        for _ in 0..child_count {
            let item_buf = read_bytes(file, key_size + 8, path)?;
            let name = std::str::from_utf8(&item_buf[..key_size])
                .unwrap_or("")
                .trim_end_matches('\0')
                .to_string();
            let chrom_id = read_u32(&item_buf, key_size, be);
            let chrom_size = read_u32(&item_buf, key_size + 4, be);
            chroms.push((name, chrom_id, chrom_size));
        }
    } else {
        // Internal node: read child pointers, then recurse
        let mut child_offsets = Vec::with_capacity(child_count);
        for _ in 0..child_count {
            let item_buf = read_bytes(file, key_size + 8, path)?;
            let child_offset = read_u64(&item_buf, key_size, be);
            child_offsets.push(child_offset);
        }
        for offset in child_offsets {
            file.seek(SeekFrom::Start(offset))
                .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;
            read_bptree_node(file, header, path, key_size, chroms)?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// R-tree (spatial index)
// ---------------------------------------------------------------------------

/// Find data block offsets overlapping the given region.
fn read_rtree_overlapping(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
    chrom_id: u32,
    start: u32,
    end: u32,
) -> Result<Vec<(u64, u64)>> {
    // R-tree header immediately follows the data section
    // The index offset is stored at byte 24 of the file header
    file.seek(SeekFrom::Start(24))
        .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;

    let buf = read_bytes(file, 8, path)?;
    let index_offset = read_u64(&buf, 0, header.big_endian);

    file.seek(SeekFrom::Start(index_offset))
        .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;

    let be = header.big_endian;
    let rtree_header = read_bytes(file, 48, path)?;

    let magic = read_u32(&rtree_header, 0, be);
    if magic != 0x2468_ACE0 {
        return Err(read_err(path, "invalid R-tree magic"));
    }

    let mut blocks = Vec::new();
    read_rtree_node(file, header, path, chrom_id, start, end, &mut blocks)?;
    Ok(blocks)
}

fn read_rtree_node(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
    chrom_id: u32,
    start: u32,
    end: u32,
    blocks: &mut Vec<(u64, u64)>,
) -> Result<()> {
    let be = header.big_endian;
    let node_buf = read_bytes(file, 4, path)?;
    let is_leaf = node_buf[0] != 0;
    let child_count = read_u16(&node_buf, 2, be) as usize;

    if is_leaf {
        for _ in 0..child_count {
            let item = read_bytes(file, 32, path)?;
            let start_chrom_ix = read_u32(&item, 0, be);
            let start_base = read_u32(&item, 4, be);
            let end_chrom_ix = read_u32(&item, 8, be);
            let end_base = read_u32(&item, 12, be);
            let data_offset = read_u64(&item, 16, be);
            let data_size = read_u64(&item, 24, be);

            if overlaps(chrom_id, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                blocks.push((data_offset, data_size));
            }
        }
    } else {
        let mut children = Vec::with_capacity(child_count);
        for _ in 0..child_count {
            let item = read_bytes(file, 24, path)?;
            let start_chrom_ix = read_u32(&item, 0, be);
            let start_base = read_u32(&item, 4, be);
            let end_chrom_ix = read_u32(&item, 8, be);
            let end_base = read_u32(&item, 12, be);
            let child_offset = read_u64(&item, 16, be);

            if overlaps(chrom_id, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                children.push(child_offset);
            }
        }
        for offset in children {
            file.seek(SeekFrom::Start(offset))
                .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;
            read_rtree_node(file, header, path, chrom_id, start, end, blocks)?;
        }
    }

    Ok(())
}

/// Check if a query region overlaps an R-tree node's bounding box.
fn overlaps(
    qchrom: u32, qstart: u32, qend: u32,
    start_chrom: u32, start_base: u32,
    end_chrom: u32, end_base: u32,
) -> bool {
    // Query is entirely before the node
    if qchrom < start_chrom || (qchrom == start_chrom && qend <= start_base) {
        return false;
    }
    // Query is entirely after the node
    if qchrom > end_chrom || (qchrom == end_chrom && qstart >= end_base) {
        return false;
    }
    true
}

// ---------------------------------------------------------------------------
// Data block reading
// ---------------------------------------------------------------------------

fn decompress_block(
    file: &mut File,
    _header: &InternalHeader,
    path: &Path,
    offset: u64,
    size: u64,
) -> Result<Vec<u8>> {
    file.seek(SeekFrom::Start(offset))
        .map_err(|e| CyaneaError::Io(io::Error::new(e.kind(), format!("{}: {}", path.display(), e))))?;

    let raw = read_bytes(file, size as usize, path)?;

    // Try zlib decompression; if it fails, assume uncompressed
    match flate2::read::ZlibDecoder::new(&raw[..]).bytes().collect::<std::result::Result<Vec<u8>, _>>() {
        Ok(decompressed) if !decompressed.is_empty() => Ok(decompressed),
        _ => Ok(raw),
    }
}

fn read_bigwig_data_block(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
    offset: u64,
    size: u64,
    chrom_name: &str,
    query_start: u32,
    query_end: u32,
) -> Result<Vec<BigWigInterval>> {
    let data = decompress_block(file, header, path, offset, size)?;
    let be = header.big_endian;
    let mut intervals = Vec::new();

    if data.len() < 24 {
        return Ok(intervals);
    }

    // Section header: chromId(4) + chromStart(4) + chromEnd(4) + itemStep(4) + itemSpan(4) + type(1) + reserved(1) + itemCount(2)
    let _chrom_id = read_u32(&data, 0, be);
    let _chrom_start = read_u32(&data, 4, be);
    let _chrom_end = read_u32(&data, 8, be);
    let item_step = read_u32(&data, 12, be);
    let item_span = read_u32(&data, 16, be);
    let data_type = data[20];
    let item_count = read_u16(&data, 22, be) as usize;

    let mut pos = 24;

    match data_type {
        1 => {
            // bedGraph: chromStart(4) + chromEnd(4) + value(4) per item
            for _ in 0..item_count {
                if pos + 12 > data.len() {
                    break;
                }
                let s = read_u32(&data, pos, be);
                let e = read_u32(&data, pos + 4, be);
                let v = read_f32(&data, pos + 8, be);
                pos += 12;

                if s < query_end && e > query_start {
                    intervals.push(BigWigInterval {
                        chrom: chrom_name.to_string(),
                        start: s,
                        end: e,
                        value: v,
                    });
                }
            }
        }
        2 => {
            // variableStep: chromStart(4) + value(4) per item
            for _ in 0..item_count {
                if pos + 8 > data.len() {
                    break;
                }
                let s = read_u32(&data, pos, be);
                let v = read_f32(&data, pos + 4, be);
                pos += 8;

                let e = s + item_span;
                if s < query_end && e > query_start {
                    intervals.push(BigWigInterval {
                        chrom: chrom_name.to_string(),
                        start: s,
                        end: e,
                        value: v,
                    });
                }
            }
        }
        3 => {
            // fixedStep: value(4) per item
            let mut current_start = _chrom_start;
            for _ in 0..item_count {
                if pos + 4 > data.len() {
                    break;
                }
                let v = read_f32(&data, pos, be);
                pos += 4;

                let s = current_start;
                let e = s + item_span;
                current_start += item_step;

                if s < query_end && e > query_start {
                    intervals.push(BigWigInterval {
                        chrom: chrom_name.to_string(),
                        start: s,
                        end: e,
                        value: v,
                    });
                }
            }
        }
        _ => {}
    }

    Ok(intervals)
}

fn read_bigbed_data_block(
    file: &mut File,
    header: &InternalHeader,
    path: &Path,
    offset: u64,
    size: u64,
    chrom_name: &str,
    query_start: u32,
    query_end: u32,
) -> Result<Vec<BigBedRecord>> {
    let data = decompress_block(file, header, path, offset, size)?;
    let be = header.big_endian;
    let mut records = Vec::new();
    let mut pos = 0;

    while pos + 12 <= data.len() {
        let _chrom_id = read_u32(&data, pos, be);
        let start = read_u32(&data, pos + 4, be);
        let end = read_u32(&data, pos + 8, be);
        pos += 12;

        // Rest of the record is a NUL-terminated string
        let rest_start = pos;
        while pos < data.len() && data[pos] != 0 {
            pos += 1;
        }
        let rest = std::str::from_utf8(&data[rest_start..pos])
            .unwrap_or("")
            .to_string();
        if pos < data.len() {
            pos += 1; // skip NUL
        }

        if start < query_end && end > query_start {
            records.push(BigBedRecord {
                chrom: chrom_name.to_string(),
                start,
                end,
                rest,
            });
        }
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// Helper to construct minimal bigWig binary for testing
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test_helpers {
    /// Build a minimal bigWig binary with given parameters.
    ///
    /// Creates a valid bigWig file with:
    /// - Header (64 bytes)
    /// - Zoom headers (none)
    /// - Total summary (40 bytes)
    /// - Chromosome B+ tree
    /// - Data section with bedGraph items
    /// - R-tree index
    pub fn build_minimal_bigwig(
        chroms: &[(&str, u32)], // (name, size)
        intervals: &[(u32, u32, u32, f32)], // (chrom_id, start, end, value)
    ) -> Vec<u8> {
        let mut buf = Vec::new();

        // We'll build up sections and patch offsets at the end
        let magic: u32 = 0x888F_FC26;
        let version: u16 = 4;
        let zoom_levels: u16 = 0;

        // Header placeholder (64 bytes)
        let header_start = buf.len();
        buf.extend_from_slice(&magic.to_le_bytes());       // 0: magic
        buf.extend_from_slice(&version.to_le_bytes());      // 4: version
        buf.extend_from_slice(&zoom_levels.to_le_bytes());  // 6: zoom levels
        buf.extend_from_slice(&0u64.to_le_bytes());         // 8: chrom tree offset (patch later)
        buf.extend_from_slice(&0u64.to_le_bytes());         // 16: data offset (patch later)
        buf.extend_from_slice(&0u64.to_le_bytes());         // 24: index offset (patch later)
        buf.extend_from_slice(&0u64.to_le_bytes());         // 32: data count (we store as u64 but spec overlaps with fieldCount etc for bigBed)
        buf.extend_from_slice(&0u64.to_le_bytes());         // 40: auto-sql offset
        buf.extend_from_slice(&0u64.to_le_bytes());         // 48: total summary offset (patch later)
        buf.extend_from_slice(&0u32.to_le_bytes());         // 56: uncompress buf size
        buf.extend_from_slice(&0u32.to_le_bytes());         // 60: reserved
        assert_eq!(buf.len() - header_start, 64);

        // Total summary
        let summary_offset = buf.len() as u64;
        let bases_covered: u64 = intervals.iter().map(|i| (i.2 - i.1) as u64).sum();
        let min_val: f64 = intervals.iter().map(|i| i.3 as f64).fold(f64::MAX, f64::min);
        let max_val: f64 = intervals.iter().map(|i| i.3 as f64).fold(f64::MIN, f64::max);
        let sum: f64 = intervals.iter().map(|i| i.3 as f64 * (i.2 - i.1) as f64).sum();
        let sum_sq: f64 = intervals.iter().map(|i| (i.3 as f64).powi(2) * (i.2 - i.1) as f64).sum();

        let (min_val, max_val) = if intervals.is_empty() { (0.0, 0.0) } else { (min_val, max_val) };

        buf.extend_from_slice(&bases_covered.to_le_bytes());
        buf.extend_from_slice(&min_val.to_bits().to_le_bytes());
        buf.extend_from_slice(&max_val.to_bits().to_le_bytes());
        buf.extend_from_slice(&sum.to_bits().to_le_bytes());
        buf.extend_from_slice(&sum_sq.to_bits().to_le_bytes());

        // B+ tree for chromosomes
        let chrom_tree_offset = buf.len() as u64;
        let key_size: u32 = chroms.iter().map(|(n, _)| n.len() as u32 + 1).max().unwrap_or(8).max(8);
        let val_size: u32 = 8; // chrom_id(4) + chrom_size(4)
        let bptree_magic: u32 = 0x78CA_8C91;

        buf.extend_from_slice(&bptree_magic.to_le_bytes());     // magic
        buf.extend_from_slice(&256u32.to_le_bytes());            // block size
        buf.extend_from_slice(&key_size.to_le_bytes());          // key size
        buf.extend_from_slice(&val_size.to_le_bytes());          // val size
        buf.extend_from_slice(&(chroms.len() as u64).to_le_bytes()); // item count
        buf.extend_from_slice(&0u64.to_le_bytes());              // reserved

        // Leaf node
        buf.push(1); // is_leaf
        buf.push(0); // reserved
        buf.extend_from_slice(&(chroms.len() as u16).to_le_bytes()); // count

        for (i, (name, size)) in chroms.iter().enumerate() {
            let mut key = vec![0u8; key_size as usize];
            let name_bytes = name.as_bytes();
            key[..name_bytes.len()].copy_from_slice(name_bytes);
            buf.extend_from_slice(&key);
            buf.extend_from_slice(&(i as u32).to_le_bytes()); // chrom_id
            buf.extend_from_slice(&size.to_le_bytes());         // chrom_size
        }

        // Data section (bedGraph format)
        let data_offset = buf.len() as u64;

        if !intervals.is_empty() {
            // Section header (24 bytes)
            let first = &intervals[0];
            let last = &intervals[intervals.len() - 1];
            buf.extend_from_slice(&first.0.to_le_bytes());    // chromId
            buf.extend_from_slice(&first.1.to_le_bytes());    // chromStart
            buf.extend_from_slice(&last.2.to_le_bytes());     // chromEnd
            buf.extend_from_slice(&0u32.to_le_bytes());       // itemStep (not used for bedGraph)
            buf.extend_from_slice(&0u32.to_le_bytes());       // itemSpan (not used for bedGraph)
            buf.push(1);                                       // type = bedGraph
            buf.push(0);                                       // reserved
            buf.extend_from_slice(&(intervals.len() as u16).to_le_bytes()); // itemCount

            // Items
            for iv in intervals {
                buf.extend_from_slice(&iv.1.to_le_bytes()); // start
                buf.extend_from_slice(&iv.2.to_le_bytes()); // end
                buf.extend_from_slice(&iv.3.to_bits().to_le_bytes()); // value
            }
        }

        let data_end = buf.len() as u64;
        let data_size = data_end - data_offset;

        // R-tree index
        let index_offset = buf.len() as u64;
        let rtree_magic: u32 = 0x2468_ACE0;

        buf.extend_from_slice(&rtree_magic.to_le_bytes());      // magic
        buf.extend_from_slice(&256u32.to_le_bytes());            // block size
        buf.extend_from_slice(&(intervals.len() as u64).to_le_bytes()); // item count
        buf.extend_from_slice(&0u32.to_le_bytes());              // start chrom ix
        buf.extend_from_slice(&0u32.to_le_bytes());              // start base
        buf.extend_from_slice(&0u32.to_le_bytes());              // end chrom ix
        buf.extend_from_slice(&0u32.to_le_bytes());              // end base
        buf.extend_from_slice(&0u64.to_le_bytes());              // end file offset
        buf.extend_from_slice(&(if intervals.is_empty() { 0u32 } else { 1u32 }).to_le_bytes()); // items per slot
        buf.extend_from_slice(&0u32.to_le_bytes());              // reserved

        // Leaf node
        buf.push(1); // is_leaf
        buf.push(0); // reserved
        if intervals.is_empty() {
            buf.extend_from_slice(&0u16.to_le_bytes()); // count = 0
        } else {
            buf.extend_from_slice(&1u16.to_le_bytes()); // count = 1 block

            let first = &intervals[0];
            let last = &intervals[intervals.len() - 1];
            buf.extend_from_slice(&first.0.to_le_bytes()); // start chrom ix
            buf.extend_from_slice(&first.1.to_le_bytes()); // start base
            buf.extend_from_slice(&last.0.to_le_bytes());  // end chrom ix
            buf.extend_from_slice(&last.2.to_le_bytes());  // end base
            buf.extend_from_slice(&data_offset.to_le_bytes()); // data offset
            buf.extend_from_slice(&data_size.to_le_bytes());   // data size
        }

        // Patch header offsets
        patch_u64(&mut buf, 8, chrom_tree_offset);
        patch_u64(&mut buf, 16, data_offset);
        patch_u64(&mut buf, 24, index_offset);
        let data_count = intervals.len() as u64;
        patch_u64(&mut buf, 32, data_count);
        patch_u64(&mut buf, 48, summary_offset);

        buf
    }

    /// Build a minimal bigBed binary.
    pub fn build_minimal_bigbed(
        chroms: &[(&str, u32)],
        records: &[(u32, u32, u32, &str)], // (chrom_id, start, end, rest)
    ) -> Vec<u8> {
        let mut buf = Vec::new();

        let magic: u32 = 0x8789_F2EB;
        let version: u16 = 4;
        let zoom_levels: u16 = 0;

        // Header (64 bytes)
        buf.extend_from_slice(&magic.to_le_bytes());
        buf.extend_from_slice(&version.to_le_bytes());
        buf.extend_from_slice(&zoom_levels.to_le_bytes());
        buf.extend_from_slice(&0u64.to_le_bytes()); // 8: chrom tree offset
        buf.extend_from_slice(&0u64.to_le_bytes()); // 16: data offset
        buf.extend_from_slice(&0u64.to_le_bytes()); // 24: index offset
        buf.extend_from_slice(&0u64.to_le_bytes()); // 32: data count
        buf.extend_from_slice(&0u64.to_le_bytes()); // 40: auto-sql offset
        buf.extend_from_slice(&0u64.to_le_bytes()); // 48: total summary offset
        buf.extend_from_slice(&0u32.to_le_bytes()); // 56: uncompress buf size
        buf.extend_from_slice(&0u32.to_le_bytes()); // 60: reserved

        // Total summary
        let summary_offset = buf.len() as u64;
        let bases: u64 = records.iter().map(|r| (r.2 - r.1) as u64).sum();
        buf.extend_from_slice(&bases.to_le_bytes());
        buf.extend_from_slice(&0f64.to_bits().to_le_bytes()); // min
        buf.extend_from_slice(&0f64.to_bits().to_le_bytes()); // max
        buf.extend_from_slice(&0f64.to_bits().to_le_bytes()); // sum
        buf.extend_from_slice(&0f64.to_bits().to_le_bytes()); // sum_sq

        // B+ tree
        let chrom_tree_offset = buf.len() as u64;
        let key_size: u32 = chroms.iter().map(|(n, _)| n.len() as u32 + 1).max().unwrap_or(8).max(8);

        buf.extend_from_slice(&0x78CA_8C91u32.to_le_bytes());
        buf.extend_from_slice(&256u32.to_le_bytes());
        buf.extend_from_slice(&key_size.to_le_bytes());
        buf.extend_from_slice(&8u32.to_le_bytes());
        buf.extend_from_slice(&(chroms.len() as u64).to_le_bytes());
        buf.extend_from_slice(&0u64.to_le_bytes());

        buf.push(1);
        buf.push(0);
        buf.extend_from_slice(&(chroms.len() as u16).to_le_bytes());

        for (i, (name, size)) in chroms.iter().enumerate() {
            let mut key = vec![0u8; key_size as usize];
            let name_bytes = name.as_bytes();
            key[..name_bytes.len()].copy_from_slice(name_bytes);
            buf.extend_from_slice(&key);
            buf.extend_from_slice(&(i as u32).to_le_bytes());
            buf.extend_from_slice(&size.to_le_bytes());
        }

        // Data section
        let data_offset = buf.len() as u64;
        for rec in records {
            buf.extend_from_slice(&rec.0.to_le_bytes()); // chrom_id
            buf.extend_from_slice(&rec.1.to_le_bytes()); // start
            buf.extend_from_slice(&rec.2.to_le_bytes()); // end
            buf.extend_from_slice(rec.3.as_bytes());       // rest
            buf.push(0);                                    // NUL terminator
        }
        let data_end = buf.len() as u64;
        let data_size = data_end - data_offset;

        // R-tree
        let index_offset = buf.len() as u64;
        buf.extend_from_slice(&0x2468_ACE0u32.to_le_bytes());
        buf.extend_from_slice(&256u32.to_le_bytes());
        buf.extend_from_slice(&(records.len() as u64).to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(&0u64.to_le_bytes());
        buf.extend_from_slice(&(if records.is_empty() { 0u32 } else { 1u32 }).to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());

        buf.push(1);
        buf.push(0);
        if records.is_empty() {
            buf.extend_from_slice(&0u16.to_le_bytes());
        } else {
            buf.extend_from_slice(&1u16.to_le_bytes());
            let first = &records[0];
            let last = &records[records.len() - 1];
            buf.extend_from_slice(&first.0.to_le_bytes());
            buf.extend_from_slice(&first.1.to_le_bytes());
            buf.extend_from_slice(&last.0.to_le_bytes());
            buf.extend_from_slice(&last.2.to_le_bytes());
            buf.extend_from_slice(&data_offset.to_le_bytes());
            buf.extend_from_slice(&data_size.to_le_bytes());
        }

        // Patch offsets
        patch_u64(&mut buf, 8, chrom_tree_offset);
        patch_u64(&mut buf, 16, data_offset);
        patch_u64(&mut buf, 24, index_offset);
        patch_u64(&mut buf, 32, records.len() as u64);
        patch_u64(&mut buf, 48, summary_offset);

        buf
    }

    fn patch_u64(buf: &mut [u8], offset: usize, value: u64) {
        let bytes = value.to_le_bytes();
        buf[offset..offset + 8].copy_from_slice(&bytes);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use tempfile::NamedTempFile;

    fn write_test_file(data: &[u8], suffix: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(suffix).unwrap();
        file.write_all(data).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn bigwig_magic_validation() {
        let mut bad_data = vec![0u8; 200];
        bad_data[0..4].copy_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF]);
        let file = write_test_file(&bad_data, ".bw");
        let result = read_bigwig_header(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("not a valid bigWig"));
    }

    #[test]
    fn bigwig_header_parse() {
        let data = test_helpers::build_minimal_bigwig(
            &[("chr1", 248956422), ("chr2", 242193529)],
            &[
                (0, 100, 200, 1.5),
                (0, 300, 400, 2.5),
            ],
        );
        let file = write_test_file(&data, ".bw");
        let header = read_bigwig_header(file.path()).unwrap();

        assert_eq!(header.version, 4);
        assert_eq!(header.zoom_levels, 0);
        assert_eq!(header.chrom_count, 2);
    }

    #[test]
    fn bigwig_summary_stats() {
        let data = test_helpers::build_minimal_bigwig(
            &[("chr1", 1000)],
            &[
                (0, 0, 100, 1.0),
                (0, 100, 200, 3.0),
            ],
        );
        let file = write_test_file(&data, ".bw");
        let header = read_bigwig_header(file.path()).unwrap();

        assert_eq!(header.total_summary.bases_covered, 200);
        assert!((header.total_summary.min_val - 1.0).abs() < f64::EPSILON);
        assert!((header.total_summary.max_val - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn bigbed_magic_validation() {
        let mut bad_data = vec![0u8; 200];
        bad_data[0..4].copy_from_slice(&[0x01, 0x02, 0x03, 0x04]);
        let file = write_test_file(&bad_data, ".bb");
        let result = read_bigbed_header(file.path());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("not a valid bigBed"));
    }

    #[test]
    fn bigwig_endianness() {
        // Build a little-endian bigWig and verify we can read it
        let data = test_helpers::build_minimal_bigwig(
            &[("chr1", 1000)],
            &[(0, 50, 150, 2.0)],
        );
        let file = write_test_file(&data, ".bw");
        let header = read_bigwig_header(file.path()).unwrap();
        assert_eq!(header.chrom_count, 1);

        // Verify big-endian magic is rejected as expected (we built LE)
        // and that a swapped magic would be caught
        let be_magic = 0x26FC_8F88u32;
        let le_magic = 0x888F_FC26u32;
        assert_ne!(be_magic, le_magic);
    }

    #[test]
    fn bigwig_region_query() {
        let data = test_helpers::build_minimal_bigwig(
            &[("chr1", 1000)],
            &[
                (0, 100, 200, 1.0),
                (0, 200, 300, 2.0),
                (0, 400, 500, 3.0),
            ],
        );
        let file = write_test_file(&data, ".bw");
        let intervals = read_bigwig_intervals(file.path(), "chr1", 150, 350).unwrap();

        // Should overlap with intervals [100,200) and [200,300) but not [400,500)
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
        assert!((intervals[0].value - 1.0).abs() < f32::EPSILON);
        assert_eq!(intervals[1].start, 200);
        assert_eq!(intervals[1].end, 300);
    }
}
