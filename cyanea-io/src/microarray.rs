//! Microarray file format parsing — Affymetrix CEL, GenePix GPR, Illumina IDAT.
//!
//! Provides parsers for legacy microarray data formats, extracting probe-level
//! intensities and metadata.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// CEL (Affymetrix)
// ---------------------------------------------------------------------------

/// Parsed Affymetrix CEL file data.
#[derive(Debug, Clone)]
pub struct CelFile {
    /// CEL file version (3 = text, 4 = binary v4, CC = Command Console).
    pub version: CelVersion,
    /// Chip type / array type (e.g., "HG-U133_Plus_2").
    pub chip_type: String,
    /// Algorithm used for intensity extraction.
    pub algorithm: String,
    /// Number of columns in the array.
    pub cols: u32,
    /// Number of rows in the array.
    pub rows: u32,
    /// Number of cells (probes).
    pub n_cells: usize,
    /// Intensity values per cell.
    pub intensities: Vec<f32>,
    /// Standard deviation per cell (if available).
    pub std_devs: Option<Vec<f32>>,
    /// Number of pixels per cell (if available).
    pub n_pixels: Option<Vec<u16>>,
    /// Outlier cell indices.
    pub outliers: Vec<usize>,
    /// Masked cell indices.
    pub masked: Vec<usize>,
    /// Header key-value pairs.
    pub header: Vec<(String, String)>,
}

/// CEL file version.
#[derive(Debug, Clone, PartialEq)]
pub enum CelVersion {
    /// Version 3 (text-based).
    V3,
    /// Version 4 (binary).
    V4,
    /// Command Console (CC/Calvin).
    CommandConsole,
}

/// Parse an Affymetrix CEL v3 (text format) file.
///
/// Extracts header metadata, intensity values, standard deviations,
/// and outlier/masked cell lists.
pub fn parse_cel_v3(data: &str) -> Result<CelFile> {
    let mut chip_type = String::new();
    let mut algorithm = String::new();
    let mut cols = 0u32;
    let mut rows = 0u32;
    let mut n_cells = 0usize;
    let mut intensities = Vec::new();
    let mut std_devs = Vec::new();
    let mut n_pixels = Vec::new();
    let mut outliers = Vec::new();
    let mut masked = Vec::new();
    let mut header = Vec::new();

    #[derive(PartialEq)]
    enum Section {
        None,
        Header,
        Intensity,
        Outliers,
        Masked,
    }
    let mut section = Section::None;

    for line in data.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        match line {
            "[HEADER]" => {
                section = Section::Header;
                continue;
            }
            "[INTENSITY]" => {
                section = Section::Intensity;
                continue;
            }
            "[OUTLIERS]" => {
                section = Section::Outliers;
                continue;
            }
            "[MASKS]" => {
                section = Section::Masked;
                continue;
            }
            _ => {}
        }

        if line.starts_with('[') {
            section = Section::None;
            continue;
        }

        match section {
            Section::Header => {
                if let Some((key, val)) = line.split_once('=') {
                    let key = key.trim();
                    let val = val.trim();
                    match key {
                        "DatHeader" | "ChipType" => chip_type = val.to_string(),
                        "Algorithm" => algorithm = val.to_string(),
                        "Cols" => cols = val.parse().unwrap_or(0),
                        "Rows" => rows = val.parse().unwrap_or(0),
                        "NumberCells" | "NumCells" => n_cells = val.parse().unwrap_or(0),
                        _ => {}
                    }
                    header.push((key.to_string(), val.to_string()));
                }
            }
            Section::Intensity => {
                // Skip header line "CellHeader=X Y MEAN STDV NPIXELS"
                if line.starts_with("CellHeader") || line.starts_with("NumberCells") {
                    if line.starts_with("NumberCells") {
                        if let Some((_, v)) = line.split_once('=') {
                            n_cells = v.trim().parse().unwrap_or(0);
                        }
                    }
                    continue;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    // X Y MEAN [STDV] [NPIXELS]
                    if let Ok(mean) = parts[2].parse::<f32>() {
                        intensities.push(mean);
                    }
                    if parts.len() >= 4 {
                        if let Ok(sd) = parts[3].parse::<f32>() {
                            std_devs.push(sd);
                        }
                    }
                    if parts.len() >= 5 {
                        if let Ok(np) = parts[4].parse::<u16>() {
                            n_pixels.push(np);
                        }
                    }
                }
            }
            Section::Outliers => {
                if line.starts_with("CellHeader") || line.starts_with("NumberCells") {
                    continue;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let (Ok(x), Ok(y)) = (parts[0].parse::<u32>(), parts[1].parse::<u32>()) {
                        outliers.push((y * cols + x) as usize);
                    }
                }
            }
            Section::Masked => {
                if line.starts_with("CellHeader") || line.starts_with("NumberCells") {
                    continue;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let (Ok(x), Ok(y)) = (parts[0].parse::<u32>(), parts[1].parse::<u32>()) {
                        masked.push((y * cols + x) as usize);
                    }
                }
            }
            Section::None => {}
        }
    }

    if intensities.is_empty() {
        return Err(CyaneaError::Parse("no intensity data found in CEL file".into()));
    }

    if n_cells == 0 {
        n_cells = intensities.len();
    }

    Ok(CelFile {
        version: CelVersion::V3,
        chip_type,
        algorithm,
        cols,
        rows,
        n_cells,
        intensities,
        std_devs: if std_devs.is_empty() { None } else { Some(std_devs) },
        n_pixels: if n_pixels.is_empty() { None } else { Some(n_pixels) },
        outliers,
        masked,
        header,
    })
}

// ---------------------------------------------------------------------------
// GPR (GenePix)
// ---------------------------------------------------------------------------

/// Parsed GenePix GPR (GenePix Results) file data.
#[derive(Debug, Clone)]
pub struct GprFile {
    /// ATF version.
    pub atf_version: String,
    /// Header key-value pairs.
    pub header: Vec<(String, String)>,
    /// Column names from the data section.
    pub columns: Vec<String>,
    /// Spot data rows.
    pub spots: Vec<GprSpot>,
}

/// A single spot in a GenePix GPR file.
#[derive(Debug, Clone)]
pub struct GprSpot {
    /// Block number.
    pub block: u32,
    /// Row within block.
    pub row: u32,
    /// Column within block.
    pub column: u32,
    /// Probe/gene name.
    pub name: String,
    /// Probe/gene ID.
    pub id: String,
    /// Foreground median (wavelength 1 / Cy5).
    pub f635_median: f64,
    /// Background median (wavelength 1 / Cy5).
    pub b635_median: f64,
    /// Foreground median (wavelength 2 / Cy3).
    pub f532_median: f64,
    /// Background median (wavelength 2 / Cy3).
    pub b532_median: f64,
    /// Ratio of medians (635/532).
    pub ratio_of_medians: Option<f64>,
    /// Log ratio (base 2).
    pub log_ratio: Option<f64>,
    /// Flags (negative = bad).
    pub flags: i32,
}

impl GprSpot {
    /// Background-subtracted intensity for channel 635 (Cy5).
    pub fn corrected_635(&self) -> f64 {
        (self.f635_median - self.b635_median).max(0.0)
    }

    /// Background-subtracted intensity for channel 532 (Cy3).
    pub fn corrected_532(&self) -> f64 {
        (self.f532_median - self.b532_median).max(0.0)
    }

    /// Whether the spot is flagged as bad.
    pub fn is_flagged(&self) -> bool {
        self.flags < 0
    }
}

/// Parse a GenePix GPR file (ATF format).
///
/// GPR files use the Axon Text File (ATF) format with a header section
/// followed by tab-delimited data.
pub fn parse_gpr(data: &str) -> Result<GprFile> {
    let mut lines = data.lines();

    // First line: ATF version
    let first_line = lines
        .next()
        .ok_or_else(|| CyaneaError::Parse("empty GPR file".into()))?;
    let atf_version = if first_line.starts_with("ATF") {
        first_line.trim().to_string()
    } else {
        return Err(CyaneaError::Parse(
            "GPR file must start with ATF header".into(),
        ));
    };

    // Second line: number of header lines and data columns
    let counts_line = lines
        .next()
        .ok_or_else(|| CyaneaError::Parse("missing ATF counts line".into()))?;
    let counts: Vec<usize> = counts_line
        .split('\t')
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    let n_header_lines = counts.first().copied().unwrap_or(0);

    // Read header key-value pairs
    let mut header = Vec::new();
    for _ in 0..n_header_lines {
        if let Some(line) = lines.next() {
            if let Some((key, val)) = line.split_once('=') {
                header.push((key.trim().trim_matches('"').to_string(), val.trim().trim_matches('"').to_string()));
            }
        }
    }

    // Column header line
    let col_line = lines
        .next()
        .ok_or_else(|| CyaneaError::Parse("missing column header".into()))?;
    let columns: Vec<String> = col_line
        .split('\t')
        .map(|s| s.trim().trim_matches('"').to_string())
        .collect();

    // Find column indices
    let find_col = |name: &str| -> Option<usize> {
        columns.iter().position(|c| c == name)
    };

    let idx_block = find_col("Block");
    let idx_row = find_col("Row");
    let idx_col = find_col("Column");
    let idx_name = find_col("Name");
    let idx_id = find_col("ID");
    let idx_f635 = find_col("F635 Median");
    let idx_b635 = find_col("B635 Median");
    let idx_f532 = find_col("F532 Median");
    let idx_b532 = find_col("B532 Median");
    let idx_ratio = find_col("Ratio of Medians (635/532)");
    let idx_log_ratio = find_col("Log Ratio (635/532)");
    let idx_flags = find_col("Flags");

    // Parse data rows
    let mut spots = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let get_str = |idx: Option<usize>| -> String {
            idx.and_then(|i| fields.get(i))
                .map(|s| s.trim().trim_matches('"').to_string())
                .unwrap_or_default()
        };
        let get_u32 = |idx: Option<usize>| -> u32 {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0)
        };
        let get_f64 = |idx: Option<usize>| -> f64 {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0.0)
        };
        let get_opt_f64 = |idx: Option<usize>| -> Option<f64> {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| s.trim().parse().ok())
        };
        let get_i32 = |idx: Option<usize>| -> i32 {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0)
        };

        spots.push(GprSpot {
            block: get_u32(idx_block),
            row: get_u32(idx_row),
            column: get_u32(idx_col),
            name: get_str(idx_name),
            id: get_str(idx_id),
            f635_median: get_f64(idx_f635),
            b635_median: get_f64(idx_b635),
            f532_median: get_f64(idx_f532),
            b532_median: get_f64(idx_b532),
            ratio_of_medians: get_opt_f64(idx_ratio),
            log_ratio: get_opt_f64(idx_log_ratio),
            flags: get_i32(idx_flags),
        });
    }

    if spots.is_empty() {
        return Err(CyaneaError::Parse("no spots found in GPR file".into()));
    }

    Ok(GprFile {
        atf_version,
        header,
        columns,
        spots,
    })
}

// ---------------------------------------------------------------------------
// IDAT (Illumina)
// ---------------------------------------------------------------------------

/// Parsed Illumina IDAT file data.
#[derive(Debug, Clone)]
pub struct IdatFile {
    /// IDAT file version.
    pub version: u64,
    /// Number of entries (probes).
    pub n_entries: usize,
    /// Illumina IDs (probe addresses).
    pub illumina_ids: Vec<u32>,
    /// Mean intensities per probe.
    pub mean_intensities: Vec<u16>,
    /// Number of beads per probe.
    pub n_beads: Vec<u8>,
    /// Standard deviations (if available).
    pub std_devs: Option<Vec<u16>>,
    /// Barcode (if available).
    pub barcode: Option<String>,
    /// Chip type (if available).
    pub chip_type: Option<String>,
    /// Manifest (if available).
    pub manifest: Option<String>,
    /// Run info entries.
    pub run_info: Vec<(String, String)>,
}

/// Parse an Illumina IDAT file from binary data.
///
/// IDAT is a binary format containing bead-level intensity data from
/// Illumina microarrays (expression and methylation arrays).
///
/// The format starts with magic bytes "IDAT", followed by version and
/// field count, then tagged data fields.
pub fn parse_idat(data: &[u8]) -> Result<IdatFile> {
    if data.len() < 12 {
        return Err(CyaneaError::Parse("IDAT file too short".into()));
    }

    // Magic bytes: "IDAT"
    if &data[0..4] != b"IDAT" {
        return Err(CyaneaError::Parse("not an IDAT file (missing magic bytes)".into()));
    }

    // Version (u64 LE at offset 4)
    let version = u64::from_le_bytes(data[4..12].try_into().map_err(|_| {
        CyaneaError::Parse("failed to read IDAT version".into())
    })?);

    if version != 3 && version != 1 {
        return Err(CyaneaError::Parse(format!(
            "unsupported IDAT version: {}",
            version
        )));
    }

    // Number of fields (u32 LE at offset 12)
    if data.len() < 16 {
        return Err(CyaneaError::Parse("IDAT truncated before field count".into()));
    }
    let n_fields = u32::from_le_bytes(data[12..16].try_into().unwrap()) as usize;

    // Parse field directory: each entry is (field_code: u16, byte_offset: u64)
    let dir_start = 16;
    let dir_entry_size = 10; // 2 + 8
    let dir_end = dir_start + n_fields * dir_entry_size;
    if data.len() < dir_end {
        return Err(CyaneaError::Parse("IDAT truncated in field directory".into()));
    }

    let mut field_offsets: Vec<(u16, u64)> = Vec::with_capacity(n_fields);
    for i in 0..n_fields {
        let base = dir_start + i * dir_entry_size;
        let code = u16::from_le_bytes(data[base..base + 2].try_into().unwrap());
        let offset = u64::from_le_bytes(data[base + 2..base + 10].try_into().unwrap());
        field_offsets.push((code, offset));
    }

    // Field codes:
    // 1000 = nSNPsRead (n_entries)
    // 102  = IlluminaID (probe addresses)
    // 104  = Mean
    // 107  = NBeads
    // 103  = SD
    // 402  = Barcode
    // 403  = ChipType
    // 401  = Manifest

    let mut n_entries = 0usize;
    let mut illumina_ids = Vec::new();
    let mut mean_intensities = Vec::new();
    let mut n_beads = Vec::new();
    let mut std_devs_raw = Vec::new();
    let mut barcode = None;
    let mut chip_type = None;
    let mut manifest = None;
    let mut run_info = Vec::new();

    for &(code, offset) in &field_offsets {
        let off = offset as usize;
        if off >= data.len() {
            continue;
        }

        match code {
            1000 => {
                // nSNPsRead: u32
                if off + 4 <= data.len() {
                    n_entries = u32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as usize;
                }
            }
            102 => {
                // IlluminaID: array of u32
                if n_entries > 0 && off + n_entries * 4 <= data.len() {
                    for i in 0..n_entries {
                        let base = off + i * 4;
                        let id = u32::from_le_bytes(data[base..base + 4].try_into().unwrap());
                        illumina_ids.push(id);
                    }
                }
            }
            104 => {
                // Mean: array of u16
                if n_entries > 0 && off + n_entries * 2 <= data.len() {
                    for i in 0..n_entries {
                        let base = off + i * 2;
                        let val = u16::from_le_bytes(data[base..base + 2].try_into().unwrap());
                        mean_intensities.push(val);
                    }
                }
            }
            107 => {
                // NBeads: array of u8
                if n_entries > 0 && off + n_entries <= data.len() {
                    n_beads = data[off..off + n_entries].to_vec();
                }
            }
            103 => {
                // SD: array of u16
                if n_entries > 0 && off + n_entries * 2 <= data.len() {
                    for i in 0..n_entries {
                        let base = off + i * 2;
                        let val = u16::from_le_bytes(data[base..base + 2].try_into().unwrap());
                        std_devs_raw.push(val);
                    }
                }
            }
            402 => {
                // Barcode: null-terminated string
                barcode = read_idat_string(data, off);
            }
            403 => {
                // ChipType
                chip_type = read_idat_string(data, off);
            }
            401 => {
                // Manifest
                manifest = read_idat_string(data, off);
            }
            300 => {
                // RunInfo: array of (key, value) strings
                if off + 4 <= data.len() {
                    let n_run = u32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as usize;
                    let mut pos = off + 4;
                    for _ in 0..n_run {
                        if let Some((k, new_pos)) = read_idat_string_at(data, pos) {
                            pos = new_pos;
                            if let Some((v, new_pos2)) = read_idat_string_at(data, pos) {
                                pos = new_pos2;
                                run_info.push((k, v));
                            }
                        } else {
                            break;
                        }
                    }
                }
            }
            _ => {}
        }
    }

    if n_entries == 0 && !illumina_ids.is_empty() {
        n_entries = illumina_ids.len();
    }

    if mean_intensities.is_empty() {
        return Err(CyaneaError::Parse("no intensity data found in IDAT file".into()));
    }

    Ok(IdatFile {
        version,
        n_entries,
        illumina_ids,
        mean_intensities,
        n_beads,
        std_devs: if std_devs_raw.is_empty() {
            None
        } else {
            Some(std_devs_raw)
        },
        barcode,
        chip_type,
        manifest,
        run_info,
    })
}

/// Read a length-prefixed string from IDAT data at the given offset.
fn read_idat_string(data: &[u8], off: usize) -> Option<String> {
    read_idat_string_at(data, off).map(|(s, _)| s)
}

/// Read a length-prefixed string and return the string and the new position.
fn read_idat_string_at(data: &[u8], off: usize) -> Option<(String, usize)> {
    if off + 1 > data.len() {
        return None;
    }
    // IDAT strings: first byte is length (or u16 LE for longer strings)
    let len = data[off] as usize;
    let start = off + 1;
    if start + len > data.len() {
        return None;
    }
    let s = String::from_utf8_lossy(&data[start..start + len]).to_string();
    Some((s, start + len))
}

// ---------------------------------------------------------------------------
// Write helpers
// ---------------------------------------------------------------------------

/// Write a CEL v3 text format from a CelFile.
pub fn write_cel_v3(cel: &CelFile) -> String {
    let mut lines = Vec::new();

    lines.push("[CEL]".to_string());
    lines.push("Version=3".to_string());
    lines.push(String::new());

    lines.push("[HEADER]".to_string());
    lines.push(format!("Cols={}", cel.cols));
    lines.push(format!("Rows={}", cel.rows));
    lines.push(format!("NumberCells={}", cel.n_cells));
    if !cel.chip_type.is_empty() {
        lines.push(format!("ChipType={}", cel.chip_type));
    }
    if !cel.algorithm.is_empty() {
        lines.push(format!("Algorithm={}", cel.algorithm));
    }
    lines.push(String::new());

    lines.push("[INTENSITY]".to_string());
    lines.push(format!("NumberCells={}", cel.intensities.len()));
    lines.push("CellHeader=X\tY\tMEAN\tSTDV\tNPIXELS".to_string());
    for (i, &intensity) in cel.intensities.iter().enumerate() {
        let x = i as u32 % cel.cols;
        let y = i as u32 / cel.cols;
        let stdv = cel
            .std_devs
            .as_ref()
            .and_then(|v| v.get(i))
            .copied()
            .unwrap_or(0.0);
        let npix = cel
            .n_pixels
            .as_ref()
            .and_then(|v| v.get(i))
            .copied()
            .unwrap_or(1);
        lines.push(format!("{}\t{}\t{:.1}\t{:.1}\t{}", x, y, intensity, stdv, npix));
    }
    lines.push(String::new());

    lines.push("[OUTLIERS]".to_string());
    lines.push(format!("NumberCells={}", cel.outliers.len()));
    lines.push("CellHeader=X\tY".to_string());
    for &idx in &cel.outliers {
        let x = idx as u32 % cel.cols;
        let y = idx as u32 / cel.cols;
        lines.push(format!("{}\t{}", x, y));
    }
    lines.push(String::new());

    lines.push("[MASKS]".to_string());
    lines.push(format!("NumberCells={}", cel.masked.len()));
    lines.push("CellHeader=X\tY".to_string());
    for &idx in &cel.masked {
        let x = idx as u32 % cel.cols;
        let y = idx as u32 / cel.cols;
        lines.push(format!("{}\t{}", x, y));
    }

    lines.join("\n")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn demo_cel_v3() -> String {
        "\
[CEL]
Version=3

[HEADER]
Cols=5
Rows=4
NumberCells=20
ChipType=HG-U133_Plus_2
Algorithm=Plier

[INTENSITY]
NumberCells=20
CellHeader=X\tY\tMEAN\tSTDV\tNPIXELS
0\t0\t150.3\t12.1\t9
1\t0\t200.5\t15.3\t9
2\t0\t80.1\t8.2\t9
3\t0\t312.7\t22.0\t9
4\t0\t95.4\t9.8\t9
0\t1\t180.2\t14.5\t9
1\t1\t220.8\t18.1\t9
2\t1\t45.3\t5.1\t9
3\t1\t500.1\t35.2\t9
4\t1\t110.9\t11.0\t9
0\t2\t160.0\t13.0\t9
1\t2\t190.5\t16.2\t9
2\t2\t70.8\t7.5\t9
3\t2\t280.3\t20.1\t9
4\t2\t125.6\t12.8\t9
0\t3\t145.1\t11.9\t9
1\t3\t210.3\t17.5\t9
2\t3\t55.2\t6.3\t9
3\t3\t350.8\t25.0\t9
4\t3\t100.0\t10.5\t9

[OUTLIERS]
NumberCells=2
CellHeader=X\tY
2\t1
4\t3

[MASKS]
NumberCells=1
CellHeader=X\tY
0\t0
"
        .to_string()
    }

    #[test]
    fn test_parse_cel_v3() {
        let cel = parse_cel_v3(&demo_cel_v3()).unwrap();
        assert_eq!(cel.version, CelVersion::V3);
        assert_eq!(cel.chip_type, "HG-U133_Plus_2");
        assert_eq!(cel.algorithm, "Plier");
        assert_eq!(cel.cols, 5);
        assert_eq!(cel.rows, 4);
        assert_eq!(cel.n_cells, 20);
        assert_eq!(cel.intensities.len(), 20);
        assert!(cel.std_devs.is_some());
        assert_eq!(cel.std_devs.as_ref().unwrap().len(), 20);
        assert!(cel.n_pixels.is_some());
    }

    #[test]
    fn test_cel_outliers_masks() {
        let cel = parse_cel_v3(&demo_cel_v3()).unwrap();
        assert_eq!(cel.outliers.len(), 2);
        // (2,1) = 1*5+2 = 7
        assert!(cel.outliers.contains(&7));
        // (4,3) = 3*5+4 = 19
        assert!(cel.outliers.contains(&19));
        assert_eq!(cel.masked.len(), 1);
        // (0,0) = 0
        assert!(cel.masked.contains(&0));
    }

    #[test]
    fn test_cel_intensity_values() {
        let cel = parse_cel_v3(&demo_cel_v3()).unwrap();
        assert!((cel.intensities[0] - 150.3).abs() < 0.1);
        assert!((cel.intensities[8] - 500.1).abs() < 0.1); // (3,1) cell
    }

    #[test]
    fn test_cel_roundtrip() {
        let cel = parse_cel_v3(&demo_cel_v3()).unwrap();
        let output = write_cel_v3(&cel);
        let cel2 = parse_cel_v3(&output).unwrap();
        assert_eq!(cel2.cols, cel.cols);
        assert_eq!(cel2.rows, cel.rows);
        assert_eq!(cel2.intensities.len(), cel.intensities.len());
        assert_eq!(cel2.outliers.len(), cel.outliers.len());
    }

    fn demo_gpr() -> String {
        "ATF\t1.0\n\
         2\t13\n\
         \"Type\"=\"GenePix Results\"\n\
         \"DateTime\"=\"2023-01-15\"\n\
         \"Block\"\t\"Row\"\t\"Column\"\t\"Name\"\t\"ID\"\t\"F635 Median\"\t\"B635 Median\"\t\"F532 Median\"\t\"B532 Median\"\t\"Ratio of Medians (635/532)\"\t\"Log Ratio (635/532)\"\t\"Flags\"\t\"Diameter\"\n\
         1\t1\t1\t\"GeneA\"\t\"PROBE001\"\t5000\t200\t3000\t150\t1.68\t0.75\t0\t100\n\
         1\t1\t2\t\"GeneB\"\t\"PROBE002\"\t8000\t250\t2000\t100\t4.08\t2.03\t0\t100\n\
         1\t2\t1\t\"Empty\"\t\"EMPTY01\"\t300\t280\t290\t270\t1.0\t0.0\t-50\t100\n\
         1\t2\t2\t\"GeneC\"\t\"PROBE003\"\t12000\t300\t6000\t200\t2.02\t1.01\t100\t100\n"
        .to_string()
    }

    #[test]
    fn test_parse_gpr() {
        let gpr = parse_gpr(&demo_gpr()).unwrap();
        assert_eq!(gpr.atf_version, "ATF\t1.0");
        assert_eq!(gpr.header.len(), 2);
        assert_eq!(gpr.spots.len(), 4);
    }

    #[test]
    fn test_gpr_spot_values() {
        let gpr = parse_gpr(&demo_gpr()).unwrap();
        let s = &gpr.spots[0];
        assert_eq!(s.name, "GeneA");
        assert_eq!(s.id, "PROBE001");
        assert!((s.f635_median - 5000.0).abs() < 0.1);
        assert!((s.b635_median - 200.0).abs() < 0.1);
        assert!((s.corrected_635() - 4800.0).abs() < 0.1);
    }

    #[test]
    fn test_gpr_flagged() {
        let gpr = parse_gpr(&demo_gpr()).unwrap();
        assert!(!gpr.spots[0].is_flagged());
        assert!(gpr.spots[2].is_flagged()); // Empty spot, flags=-50
        assert!(!gpr.spots[3].is_flagged()); // flags=100 (positive = good)
    }

    #[test]
    fn test_gpr_background_correction() {
        let gpr = parse_gpr(&demo_gpr()).unwrap();
        let s = &gpr.spots[1]; // GeneB
        assert!((s.corrected_635() - 7750.0).abs() < 0.1);
        assert!((s.corrected_532() - 1900.0).abs() < 0.1);
    }

    #[test]
    fn test_parse_idat() {
        // Build a minimal IDAT binary
        let mut data = Vec::new();
        // Magic
        data.extend_from_slice(b"IDAT");
        // Version (u64 LE) = 3
        data.extend_from_slice(&3u64.to_le_bytes());
        // Number of fields (u32 LE) = 3 (nSNPsRead, Mean, IlluminaID)
        data.extend_from_slice(&3u32.to_le_bytes());

        // Field directory:
        // Each entry: code (u16 LE) + offset (u64 LE)
        let dir_end = 16 + 3 * 10; // 46

        // Field 0: nSNPsRead (1000) at offset dir_end
        data.extend_from_slice(&1000u16.to_le_bytes());
        data.extend_from_slice(&(dir_end as u64).to_le_bytes());

        // Field 1: IlluminaID (102) at offset dir_end + 4
        let id_offset = dir_end + 4;
        data.extend_from_slice(&102u16.to_le_bytes());
        data.extend_from_slice(&(id_offset as u64).to_le_bytes());

        // Field 2: Mean (104) at offset dir_end + 4 + 3*4 = dir_end + 16
        let mean_offset = id_offset + 3 * 4;
        data.extend_from_slice(&104u16.to_le_bytes());
        data.extend_from_slice(&(mean_offset as u64).to_le_bytes());

        // Data:
        // nSNPsRead = 3
        data.extend_from_slice(&3u32.to_le_bytes());
        // IlluminaIDs: 100, 200, 300
        data.extend_from_slice(&100u32.to_le_bytes());
        data.extend_from_slice(&200u32.to_le_bytes());
        data.extend_from_slice(&300u32.to_le_bytes());
        // Mean intensities: 1000, 2000, 3000
        data.extend_from_slice(&1000u16.to_le_bytes());
        data.extend_from_slice(&2000u16.to_le_bytes());
        data.extend_from_slice(&3000u16.to_le_bytes());

        let idat = parse_idat(&data).unwrap();
        assert_eq!(idat.version, 3);
        assert_eq!(idat.n_entries, 3);
        assert_eq!(idat.illumina_ids, vec![100, 200, 300]);
        assert_eq!(idat.mean_intensities, vec![1000, 2000, 3000]);
    }

    #[test]
    fn test_idat_bad_magic() {
        let data = b"NOTX\x03\x00\x00\x00\x00\x00\x00\x00";
        assert!(parse_idat(data).is_err());
    }

    #[test]
    fn test_idat_too_short() {
        let data = b"IDAT";
        assert!(parse_idat(data).is_err());
    }
}
