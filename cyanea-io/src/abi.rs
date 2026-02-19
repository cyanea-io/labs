//! ABI chromatogram (`.ab1`) binary file parser.
//!
//! Parses Applied Biosystems ABIF (ABI) Sanger sequencing trace files.
//! Extracts called bases, quality scores, raw trace data for all four
//! channels (A, C, G, T), peak positions, and the sample name.
//!
//! The ABIF format stores data in a tagged directory structure where each
//! entry is identified by a 4-character tag name and a tag number.

use cyanea_core::{CyaneaError, Result};

/// A parsed ABI chromatogram record.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AbiRecord {
    /// Called base sequence (ASCII: A, C, G, T, N).
    pub sequence: Vec<u8>,
    /// Phred-style quality values (one per called base).
    pub quality: Vec<u8>,
    /// Raw fluorescence traces for the four channels.
    pub traces: AbiTraces,
    /// Sample name from the instrument.
    pub sample_name: String,
    /// Peak positions (trace data point index for each called base).
    pub peak_positions: Vec<u16>,
}

/// Raw fluorescence trace data for the four nucleotide channels.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AbiTraces {
    /// Adenine channel intensities.
    pub a: Vec<i16>,
    /// Cytosine channel intensities.
    pub c: Vec<i16>,
    /// Guanine channel intensities.
    pub g: Vec<i16>,
    /// Thymine channel intensities.
    pub t: Vec<i16>,
}

/// ABIF magic bytes.
const ABIF_MAGIC: &[u8; 4] = b"ABIF";

/// Size of each directory entry in bytes.
const DIR_ENTRY_SIZE: usize = 28;

/// A parsed directory entry from the ABIF file.
#[derive(Debug)]
struct DirEntry {
    /// 4-character tag name.
    tag_name: [u8; 4],
    /// Tag number (distinguishes multiple entries with the same tag name).
    tag_number: i32,
    /// Element type code.
    _element_type: i16,
    /// Size of each element in bytes.
    _element_size: i16,
    /// Number of elements.
    num_elements: i32,
    /// Total data size in bytes.
    data_size: i32,
    /// File offset to the data (or inline data if data_size <= 4).
    data_offset: i32,
}

/// Read a big-endian i16 from a byte slice at the given offset.
fn read_i16_be(data: &[u8], offset: usize) -> Result<i16> {
    if offset + 2 > data.len() {
        return Err(CyaneaError::Parse(format!(
            "ABI: unexpected EOF reading i16 at offset {}",
            offset
        )));
    }
    Ok(i16::from_be_bytes([data[offset], data[offset + 1]]))
}

/// Read a big-endian i32 from a byte slice at the given offset.
fn read_i32_be(data: &[u8], offset: usize) -> Result<i32> {
    if offset + 4 > data.len() {
        return Err(CyaneaError::Parse(format!(
            "ABI: unexpected EOF reading i32 at offset {}",
            offset
        )));
    }
    Ok(i32::from_be_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ]))
}

/// Parse a single directory entry from the raw bytes.
fn parse_dir_entry(data: &[u8], offset: usize) -> Result<DirEntry> {
    if offset + DIR_ENTRY_SIZE > data.len() {
        return Err(CyaneaError::Parse(format!(
            "ABI: directory entry at offset {} exceeds file size",
            offset
        )));
    }

    let mut tag_name = [0u8; 4];
    tag_name.copy_from_slice(&data[offset..offset + 4]);

    Ok(DirEntry {
        tag_name,
        tag_number: read_i32_be(data, offset + 4)?,
        _element_type: read_i16_be(data, offset + 8)?,
        _element_size: read_i16_be(data, offset + 10)?,
        num_elements: read_i32_be(data, offset + 12)?,
        data_size: read_i32_be(data, offset + 16)?,
        data_offset: read_i32_be(data, offset + 20)?,
    })
}

/// Get the data region for a directory entry.
///
/// If `data_size <= 4`, the data is stored inline in the offset field
/// (bytes 20..24 of the directory entry). Otherwise, `data_offset` is
/// a file offset.
fn get_entry_data<'a>(
    file_data: &'a [u8],
    entry: &DirEntry,
    entry_file_offset: usize,
) -> Result<&'a [u8]> {
    let size = entry.data_size as usize;
    if size <= 4 {
        // Data is stored inline at the offset field position (bytes 20..24)
        let inline_start = entry_file_offset + 20;
        if inline_start + size > file_data.len() {
            return Err(CyaneaError::Parse(
                "ABI: inline data exceeds file size".to_string(),
            ));
        }
        Ok(&file_data[inline_start..inline_start + size])
    } else {
        let offset = entry.data_offset as usize;
        if offset + size > file_data.len() {
            return Err(CyaneaError::Parse(format!(
                "ABI: data at offset {} with size {} exceeds file size {}",
                offset,
                size,
                file_data.len()
            )));
        }
        Ok(&file_data[offset..offset + size])
    }
}

/// Read an array of big-endian i16 values from raw bytes.
fn read_i16_array(data: &[u8], count: usize) -> Result<Vec<i16>> {
    if data.len() < count * 2 {
        return Err(CyaneaError::Parse(format!(
            "ABI: expected {} bytes for {} i16 values, got {}",
            count * 2,
            count,
            data.len()
        )));
    }
    let mut values = Vec::with_capacity(count);
    for i in 0..count {
        let offset = i * 2;
        values.push(i16::from_be_bytes([data[offset], data[offset + 1]]));
    }
    Ok(values)
}

/// Read an array of big-endian u16 values from raw bytes.
fn read_u16_array(data: &[u8], count: usize) -> Result<Vec<u16>> {
    if data.len() < count * 2 {
        return Err(CyaneaError::Parse(format!(
            "ABI: expected {} bytes for {} u16 values, got {}",
            count * 2,
            count,
            data.len()
        )));
    }
    let mut values = Vec::with_capacity(count);
    for i in 0..count {
        let offset = i * 2;
        values.push(u16::from_be_bytes([data[offset], data[offset + 1]]));
    }
    Ok(values)
}

/// Parse an ABI chromatogram file from raw bytes.
///
/// Extracts called bases (PBAS), quality scores (PCON), trace data
/// (DATA channels 9-12), peak positions (PLOC), and sample name (SMPL).
/// The channel-to-base mapping is determined from the filter wheel order
/// tag (FWO_), defaulting to `ACGT` if not present.
///
/// # Errors
///
/// Returns an error if the data does not start with the `ABIF` magic bytes,
/// the version is unsupported, or required tags are missing.
pub fn parse_abi_bytes(data: &[u8]) -> Result<AbiRecord> {
    // Check minimum size and magic bytes
    if data.len() < 128 {
        return Err(CyaneaError::Parse(
            "ABI: file too small for ABIF header".to_string(),
        ));
    }

    if &data[0..4] != ABIF_MAGIC {
        return Err(CyaneaError::Parse(format!(
            "ABI: invalid magic bytes (expected 'ABIF', got '{}')",
            String::from_utf8_lossy(&data[0..4])
        )));
    }

    // Version check (big-endian i16 at offset 4)
    let version = read_i16_be(data, 4)?;
    if version < 100 {
        return Err(CyaneaError::Parse(format!(
            "ABI: unsupported version {} (expected >= 100)",
            version
        )));
    }

    // The header directory entry starts at offset 6 and is 28 bytes.
    // It tells us where the actual directory entries are.
    // num_elements (at offset 6+12=18) = number of directory entries
    // data_offset (at offset 6+20=26) = file offset of directory
    let num_entries = read_i32_be(data, 18)? as usize;
    let dir_offset = read_i32_be(data, 26)? as usize;

    if dir_offset + num_entries * DIR_ENTRY_SIZE > data.len() {
        return Err(CyaneaError::Parse(
            "ABI: directory extends beyond file".to_string(),
        ));
    }

    // Parse all directory entries
    let mut entries = Vec::with_capacity(num_entries);
    for i in 0..num_entries {
        let offset = dir_offset + i * DIR_ENTRY_SIZE;
        entries.push((offset, parse_dir_entry(data, offset)?));
    }

    // Find entries by tag name and number
    let find_entry = |tag: &[u8; 4], number: i32| -> Option<usize> {
        entries
            .iter()
            .position(|(_, e)| &e.tag_name == tag && e.tag_number == number)
    };

    // Determine channel order from FWO_ tag (filter wheel order)
    let channel_order: [u8; 4] = if let Some(idx) = find_entry(b"FWO_", 1) {
        let (entry_offset, ref entry) = entries[idx];
        let fwo_data = get_entry_data(data, entry, entry_offset)?;
        if fwo_data.len() >= 4 {
            [fwo_data[0], fwo_data[1], fwo_data[2], fwo_data[3]]
        } else {
            *b"ACGT"
        }
    } else {
        *b"ACGT"
    };

    // Extract called bases (PBAS, tag number 1 or 2)
    let pbas_idx = find_entry(b"PBAS", 2)
        .or_else(|| find_entry(b"PBAS", 1))
        .ok_or_else(|| CyaneaError::Parse("ABI: missing PBAS (called bases) tag".to_string()))?;
    let (pbas_offset, ref pbas_entry) = entries[pbas_idx];
    let pbas_data = get_entry_data(data, pbas_entry, pbas_offset)?;
    let sequence = pbas_data[..pbas_entry.num_elements as usize].to_vec();

    // Extract quality values (PCON, tag number 1 or 2)
    let pcon_idx = find_entry(b"PCON", 2)
        .or_else(|| find_entry(b"PCON", 1))
        .ok_or_else(|| {
            CyaneaError::Parse("ABI: missing PCON (quality values) tag".to_string())
        })?;
    let (pcon_offset, ref pcon_entry) = entries[pcon_idx];
    let pcon_data = get_entry_data(data, pcon_entry, pcon_offset)?;
    let quality = pcon_data[..pcon_entry.num_elements as usize].to_vec();

    // Extract trace data (DATA tags 9, 10, 11, 12)
    let mut trace_channels: [Vec<i16>; 4] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    for (i, data_num) in [9i32, 10, 11, 12].iter().enumerate() {
        let idx = find_entry(b"DATA", *data_num).ok_or_else(|| {
            CyaneaError::Parse(format!("ABI: missing DATA.{} (trace channel) tag", data_num))
        })?;
        let (entry_offset, ref entry) = entries[idx];
        let channel_data = get_entry_data(data, entry, entry_offset)?;
        trace_channels[i] = read_i16_array(channel_data, entry.num_elements as usize)?;
    }

    // Map channels to bases according to filter wheel order
    let mut trace_a = Vec::new();
    let mut trace_c = Vec::new();
    let mut trace_g = Vec::new();
    let mut trace_t = Vec::new();

    for (i, &base) in channel_order.iter().enumerate() {
        match base {
            b'A' | b'a' => trace_a = trace_channels[i].clone(),
            b'C' | b'c' => trace_c = trace_channels[i].clone(),
            b'G' | b'g' => trace_g = trace_channels[i].clone(),
            b'T' | b't' => trace_t = trace_channels[i].clone(),
            _ => {}
        }
    }

    // Extract peak positions (PLOC, tag number 1 or 2)
    let ploc_idx = find_entry(b"PLOC", 2)
        .or_else(|| find_entry(b"PLOC", 1))
        .ok_or_else(|| {
            CyaneaError::Parse("ABI: missing PLOC (peak positions) tag".to_string())
        })?;
    let (ploc_offset, ref ploc_entry) = entries[ploc_idx];
    let ploc_data = get_entry_data(data, ploc_entry, ploc_offset)?;
    let peak_positions = read_u16_array(ploc_data, ploc_entry.num_elements as usize)?;

    // Extract sample name (SMPL, tag number 1)
    let sample_name = if let Some(idx) = find_entry(b"SMPL", 1) {
        let (entry_offset, ref entry) = entries[idx];
        let smpl_data = get_entry_data(data, entry, entry_offset)?;
        String::from_utf8_lossy(&smpl_data[..entry.num_elements as usize]).to_string()
    } else {
        String::new()
    };

    Ok(AbiRecord {
        sequence,
        quality,
        traces: AbiTraces {
            a: trace_a,
            c: trace_c,
            g: trace_g,
            t: trace_t,
        },
        sample_name,
        peak_positions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal synthetic ABIF file for testing.
    ///
    /// Creates a valid ABIF binary with directory entries for PBAS, PCON,
    /// DATA.9-12, PLOC, SMPL, and FWO_.
    fn build_test_abi(sequence: &[u8], quality: &[u8], sample_name: &str) -> Vec<u8> {
        let seq_len = sequence.len();
        let trace_len: usize = 50; // 50 data points per trace channel
        let trace_bytes = trace_len * 2; // i16 = 2 bytes each
        let peak_bytes = seq_len * 2; // u16 = 2 bytes each

        // We'll build the file as: header (128 bytes) + data section + directory
        // Data section layout (starting at offset 128):
        //   - sequence bytes (PBAS)
        //   - quality bytes (PCON)
        //   - trace channel 0 (DATA.9) — trace_bytes
        //   - trace channel 1 (DATA.10) — trace_bytes
        //   - trace channel 2 (DATA.11) — trace_bytes
        //   - trace channel 3 (DATA.12) — trace_bytes
        //   - peak positions (PLOC) — peak_bytes
        //   - sample name (SMPL)
        let data_start: usize = 128;
        let pbas_offset = data_start;
        let pcon_offset = pbas_offset + seq_len;
        let data9_offset = pcon_offset + seq_len;
        let data10_offset = data9_offset + trace_bytes;
        let data11_offset = data10_offset + trace_bytes;
        let data12_offset = data11_offset + trace_bytes;
        let ploc_offset = data12_offset + trace_bytes;
        let smpl_offset = ploc_offset + peak_bytes;
        let dir_offset = smpl_offset + sample_name.len();

        // Number of directory entries: PBAS, PCON, DATA.9-12, PLOC, SMPL, FWO_ = 9
        let num_entries: i32 = 9;

        // Total file size
        let total_size = dir_offset + (num_entries as usize) * DIR_ENTRY_SIZE;
        let mut buf = vec![0u8; total_size];

        // --- Header ---
        // Magic
        buf[0..4].copy_from_slice(b"ABIF");
        // Version = 101
        buf[4..6].copy_from_slice(&101i16.to_be_bytes());

        // Header directory entry (at offset 6, 28 bytes)
        // tag_name = "tdir"
        buf[6..10].copy_from_slice(b"tdir");
        // tag_number = 1
        buf[10..14].copy_from_slice(&1i32.to_be_bytes());
        // element_type = 1023
        buf[14..16].copy_from_slice(&1023i16.to_be_bytes());
        // element_size = 28
        buf[16..18].copy_from_slice(&28i16.to_be_bytes());
        // num_elements = number of directory entries
        buf[18..22].copy_from_slice(&num_entries.to_be_bytes());
        // data_size = num_entries * 28
        let dir_data_size = num_entries * DIR_ENTRY_SIZE as i32;
        buf[22..26].copy_from_slice(&dir_data_size.to_be_bytes());
        // data_offset = directory offset
        buf[26..30].copy_from_slice(&(dir_offset as i32).to_be_bytes());

        // --- Data section ---
        // PBAS data (called bases)
        buf[pbas_offset..pbas_offset + seq_len].copy_from_slice(sequence);

        // PCON data (quality values)
        buf[pcon_offset..pcon_offset + seq_len].copy_from_slice(quality);

        // DATA.9-12 (trace channels — fill with simple ramp patterns)
        for ch in 0..4u16 {
            let ch_offset = data9_offset + (ch as usize) * trace_bytes;
            for i in 0..trace_len {
                let value = ((ch as i16 + 1) * 100 + i as i16) as i16;
                let byte_offset = ch_offset + i * 2;
                buf[byte_offset..byte_offset + 2].copy_from_slice(&value.to_be_bytes());
            }
        }

        // PLOC data (peak positions — evenly spaced)
        for i in 0..seq_len {
            let pos = ((i * trace_len) / seq_len.max(1)) as u16;
            let byte_offset = ploc_offset + i * 2;
            buf[byte_offset..byte_offset + 2].copy_from_slice(&pos.to_be_bytes());
        }

        // SMPL data (sample name)
        buf[smpl_offset..smpl_offset + sample_name.len()]
            .copy_from_slice(sample_name.as_bytes());

        // --- Directory entries ---
        fn write_entry(
            buf: &mut [u8],
            dir_offset: usize,
            idx: usize,
            tag: &[u8; 4],
            number: i32,
            elem_type: i16,
            elem_size: i16,
            num_elems: i32,
            d_size: i32,
            d_offset: i32,
        ) {
            let base = dir_offset + idx * 28;
            buf[base..base + 4].copy_from_slice(tag);
            buf[base + 4..base + 8].copy_from_slice(&number.to_be_bytes());
            buf[base + 8..base + 10].copy_from_slice(&elem_type.to_be_bytes());
            buf[base + 10..base + 12].copy_from_slice(&elem_size.to_be_bytes());
            buf[base + 12..base + 16].copy_from_slice(&num_elems.to_be_bytes());
            buf[base + 16..base + 20].copy_from_slice(&d_size.to_be_bytes());
            buf[base + 20..base + 24].copy_from_slice(&d_offset.to_be_bytes());
        }

        // Entry 0: FWO_ (filter wheel order) — 4 bytes inline = "ACGT"
        // For inline data (size <= 4), the data_offset field contains the data itself.
        // We'll write "ACGT" as bytes into the offset field.
        {
            let base = dir_offset;
            buf[base..base + 4].copy_from_slice(b"FWO_");
            buf[base + 4..base + 8].copy_from_slice(&1i32.to_be_bytes());
            buf[base + 8..base + 10].copy_from_slice(&2i16.to_be_bytes()); // char type
            buf[base + 10..base + 12].copy_from_slice(&1i16.to_be_bytes()); // 1 byte each
            buf[base + 12..base + 16].copy_from_slice(&4i32.to_be_bytes()); // 4 elements
            buf[base + 16..base + 20].copy_from_slice(&4i32.to_be_bytes()); // 4 bytes total
            // Inline data: "ACGT"
            buf[base + 20] = b'A';
            buf[base + 21] = b'C';
            buf[base + 22] = b'G';
            buf[base + 23] = b'T';
        }

        // Entry 1: PBAS (called bases)
        write_entry(
            &mut buf, dir_offset,
            1,
            b"PBAS",
            2,
            2, // char
            1, // 1 byte per element
            seq_len as i32,
            seq_len as i32,
            pbas_offset as i32,
        );

        // Entry 2: PCON (quality values)
        write_entry(
            &mut buf, dir_offset,
            2,
            b"PCON",
            2,
            1, // byte
            1,
            seq_len as i32,
            seq_len as i32,
            pcon_offset as i32,
        );

        // Entries 3-6: DATA.9, DATA.10, DATA.11, DATA.12
        for ch in 0..4u32 {
            let ch_offset = data9_offset + (ch as usize) * trace_bytes;
            write_entry(
                &mut buf, dir_offset,
                3 + ch as usize,
                b"DATA",
                9 + ch as i32,
                4, // short (i16)
                2, // 2 bytes per element
                trace_len as i32,
                trace_bytes as i32,
                ch_offset as i32,
            );
        }

        // Entry 7: PLOC (peak locations)
        write_entry(
            &mut buf, dir_offset,
            7,
            b"PLOC",
            2,
            5, // unsigned short (u16)
            2,
            seq_len as i32,
            peak_bytes as i32,
            ploc_offset as i32,
        );

        // Entry 8: SMPL (sample name)
        write_entry(
            &mut buf, dir_offset,
            8,
            b"SMPL",
            1,
            18, // pString
            1,
            sample_name.len() as i32,
            sample_name.len() as i32,
            smpl_offset as i32,
        );

        buf
    }

    #[test]
    fn abi_valid_sequence() {
        let seq = b"ACGTACGT";
        let qual = vec![40u8; 8];
        let data = build_test_abi(seq, &qual, "TestSample");

        let record = parse_abi_bytes(&data).unwrap();
        assert_eq!(record.sequence, b"ACGTACGT");
    }

    #[test]
    fn abi_quality_values() {
        let seq = b"ACGTNN";
        let qual = vec![30, 35, 40, 45, 10, 5];
        let data = build_test_abi(seq, &qual, "QualTest");

        let record = parse_abi_bytes(&data).unwrap();
        assert_eq!(record.quality.len(), seq.len());
        assert_eq!(record.quality, qual);
    }

    #[test]
    fn abi_trace_consistency() {
        let seq = b"ATGC";
        let qual = vec![40u8; 4];
        let data = build_test_abi(seq, &qual, "TraceTest");

        let record = parse_abi_bytes(&data).unwrap();

        // All four channels should have the same length (50 data points)
        assert_eq!(record.traces.a.len(), 50);
        assert_eq!(record.traces.c.len(), 50);
        assert_eq!(record.traces.g.len(), 50);
        assert_eq!(record.traces.t.len(), 50);

        // Verify channel data is non-trivial (our ramp pattern)
        assert!(record.traces.a.iter().any(|&v| v > 0));
        assert!(record.traces.c.iter().any(|&v| v > 0));
        assert!(record.traces.g.iter().any(|&v| v > 0));
        assert!(record.traces.t.iter().any(|&v| v > 0));
    }

    #[test]
    fn abi_sample_name() {
        let seq = b"ACGT";
        let qual = vec![40u8; 4];
        let data = build_test_abi(seq, &qual, "MySample_2024");

        let record = parse_abi_bytes(&data).unwrap();
        assert_eq!(record.sample_name, "MySample_2024");
    }

    #[test]
    fn abi_invalid_magic() {
        let mut data = vec![0u8; 256];
        data[0..4].copy_from_slice(b"NOPE");

        let result = parse_abi_bytes(&data);
        assert!(result.is_err());

        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("magic"));
    }
}
