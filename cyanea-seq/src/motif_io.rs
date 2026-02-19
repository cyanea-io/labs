//! Motif format I/O — MEME, TRANSFAC, and JASPAR parsers/writers plus PWM comparison.
//!
//! Supports reading and writing position weight matrices in three standard formats
//! used by transcription factor binding motif databases, plus a Pearson correlation-based
//! similarity measure for comparing motifs.

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

/// Alphabet for a position weight matrix.
#[derive(Debug, Clone, PartialEq)]
pub enum MotifAlphabet {
    Dna,
    Rna,
    Protein,
}

/// A sequence motif represented as a position weight matrix.
#[derive(Debug, Clone)]
pub struct Motif {
    /// Name or identifier of the motif.
    pub name: String,
    /// Alphabet type (DNA, RNA, or Protein).
    pub alphabet: MotifAlphabet,
    /// Position weight matrix: rows = positions, cols = A/C/G/T (or U).
    pub matrix: Vec<[f64; 4]>,
    /// Arbitrary metadata key-value pairs from the source format.
    pub metadata: HashMap<String, String>,
}

// ---------------------------------------------------------------------------
// MEME format
// ---------------------------------------------------------------------------

/// Parse motifs from MEME minimal format.
///
/// Expects a header starting with `MEME version`, followed by one or more
/// `MOTIF` blocks each containing a `letter-probability matrix:` line and
/// rows of four space-separated probability values.
///
/// # Errors
///
/// Returns [`CyaneaError::Parse`] on malformed input (missing header, bad
/// numbers, row count mismatch).
pub fn parse_meme(input: &str) -> Result<Vec<Motif>> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    let lines: Vec<&str> = input.lines().collect();

    // Validate header.
    let has_header = lines.iter().any(|l| l.trim().starts_with("MEME version"));
    if !has_header {
        return Err(CyaneaError::Parse(
            "MEME format must start with 'MEME version' header".into(),
        ));
    }

    let mut motifs = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        if line.starts_with("MOTIF") {
            let name = line
                .strip_prefix("MOTIF")
                .unwrap()
                .trim()
                .split_whitespace()
                .next()
                .unwrap_or("unnamed")
                .to_string();

            // Advance to find the letter-probability matrix line.
            i += 1;
            let mut width: Option<usize> = None;
            while i < lines.len() {
                let mline = lines[i].trim();
                if mline.starts_with("letter-probability matrix:") {
                    // Parse "w=" param from the line, e.g. "alength= 4 w= 8"
                    if let Some(w_idx) = mline.find("w=") {
                        let after_w = &mline[w_idx + 2..];
                        let w_str = after_w.trim().split_whitespace().next().unwrap_or("0");
                        if let Ok(w) = w_str.parse::<usize>() {
                            width = Some(w);
                        }
                    }
                    i += 1;
                    break;
                }
                i += 1;
            }

            // Read matrix rows.
            let mut matrix = Vec::new();
            while i < lines.len() {
                let rline = lines[i].trim();
                if rline.is_empty() || rline.starts_with("MOTIF") || rline.starts_with("URL") {
                    break;
                }
                let vals: std::result::Result<Vec<f64>, _> = rline
                    .split_whitespace()
                    .map(|s| s.parse::<f64>())
                    .collect();
                match vals {
                    Ok(v) if v.len() == 4 => {
                        matrix.push([v[0], v[1], v[2], v[3]]);
                        i += 1;
                    }
                    _ => break,
                }
            }

            // Validate width if specified.
            if let Some(w) = width {
                if matrix.len() != w {
                    return Err(CyaneaError::Parse(format!(
                        "MEME motif '{}': expected {} rows but found {}",
                        name,
                        w,
                        matrix.len()
                    )));
                }
            }

            if matrix.is_empty() {
                return Err(CyaneaError::Parse(format!(
                    "MEME motif '{}': no matrix rows found",
                    name
                )));
            }

            motifs.push(Motif {
                name,
                alphabet: MotifAlphabet::Dna,
                matrix,
                metadata: HashMap::new(),
            });
        } else {
            i += 1;
        }
    }

    Ok(motifs)
}

/// Write motifs in MEME minimal format.
pub fn write_meme(motifs: &[Motif]) -> String {
    let mut out = String::new();
    out.push_str("MEME version 5\n\n");
    out.push_str("ALPHABET= ACGT\n\n");
    out.push_str("strands: + -\n\n");

    for motif in motifs {
        out.push_str(&format!("MOTIF {}\n", motif.name));
        out.push_str(&format!(
            "letter-probability matrix: alength= 4 w= {}\n",
            motif.matrix.len()
        ));
        for row in &motif.matrix {
            out.push_str(&format!(
                " {:.6}  {:.6}  {:.6}  {:.6}\n",
                row[0], row[1], row[2], row[3]
            ));
        }
        out.push('\n');
    }

    out
}

// ---------------------------------------------------------------------------
// TRANSFAC format
// ---------------------------------------------------------------------------

/// Parse motifs from TRANSFAC format.
///
/// Each record starts with `ID  name`, contains a `P0` header line followed
/// by numbered rows of four values, and ends with `//`.
///
/// # Errors
///
/// Returns [`CyaneaError::Parse`] on malformed rows or missing structure.
pub fn parse_transfac(input: &str) -> Result<Vec<Motif>> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    let mut motifs = Vec::new();
    let mut name = String::new();
    let mut metadata = HashMap::new();
    let mut matrix: Vec<[f64; 4]> = Vec::new();
    let mut in_matrix = false;

    for line in input.lines() {
        let trimmed = line.trim();

        if trimmed.starts_with("ID") {
            // Start of a new record.
            name = trimmed
                .strip_prefix("ID")
                .unwrap()
                .trim()
                .to_string();
            metadata.clear();
            matrix.clear();
            in_matrix = false;
        } else if trimmed.starts_with("BF") {
            let val = trimmed.strip_prefix("BF").unwrap().trim().to_string();
            if !val.is_empty() {
                metadata.insert("BF".to_string(), val);
            }
        } else if trimmed.starts_with("P0") || trimmed.starts_with("PO") {
            in_matrix = true;
        } else if trimmed == "//" {
            // End of record.
            if !matrix.is_empty() {
                motifs.push(Motif {
                    name: name.clone(),
                    alphabet: MotifAlphabet::Dna,
                    matrix: matrix.clone(),
                    metadata: metadata.clone(),
                });
            }
            name.clear();
            metadata.clear();
            matrix.clear();
            in_matrix = false;
        } else if trimmed == "XX" {
            // Section separator — stop reading matrix rows.
            in_matrix = false;
        } else if in_matrix && !trimmed.is_empty() {
            // Numbered matrix row, e.g. "01   0.25   0.25   0.25   0.25"
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 5 {
                // First element is the row number, next 4 are values.
                let vals: std::result::Result<Vec<f64>, _> =
                    parts[1..5].iter().map(|s| s.parse::<f64>()).collect();
                match vals {
                    Ok(v) => {
                        // Normalize to probabilities if they look like counts (sum > 2).
                        let sum: f64 = v.iter().sum();
                        if sum > 2.0 {
                            let s = if sum > 0.0 { sum } else { 1.0 };
                            matrix.push([v[0] / s, v[1] / s, v[2] / s, v[3] / s]);
                        } else {
                            matrix.push([v[0], v[1], v[2], v[3]]);
                        }
                    }
                    Err(e) => {
                        return Err(CyaneaError::Parse(format!(
                            "TRANSFAC: failed to parse matrix row '{}': {}",
                            trimmed, e
                        )));
                    }
                }
            }
        }
    }

    // Handle file without trailing //.
    if !matrix.is_empty() {
        motifs.push(Motif {
            name,
            alphabet: MotifAlphabet::Dna,
            matrix,
            metadata,
        });
    }

    Ok(motifs)
}

/// Write motifs in TRANSFAC format.
pub fn write_transfac(motifs: &[Motif]) -> String {
    let mut out = String::new();

    for motif in motifs {
        out.push_str(&format!("ID  {}\n", motif.name));
        if let Some(bf) = motif.metadata.get("BF") {
            out.push_str(&format!("BF  {}\n", bf));
        }
        out.push_str("P0      A      C      G      T\n");
        for (i, row) in motif.matrix.iter().enumerate() {
            out.push_str(&format!(
                "{:02}   {:.4}   {:.4}   {:.4}   {:.4}\n",
                i + 1,
                row[0],
                row[1],
                row[2],
                row[3]
            ));
        }
        out.push_str("XX\n");
        out.push_str("//\n");
    }

    out
}

// ---------------------------------------------------------------------------
// JASPAR format
// ---------------------------------------------------------------------------

/// Parse motifs from JASPAR format.
///
/// Each motif starts with `>ID name` followed by four rows labelled A, C, G, T
/// containing bracket-delimited count values.
///
/// # Errors
///
/// Returns [`CyaneaError::Parse`] on malformed rows (missing brackets, wrong
/// row count, or unparseable numbers).
pub fn parse_jaspar(input: &str) -> Result<Vec<Motif>> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    let lines: Vec<&str> = input.lines().collect();
    let mut motifs = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        if line.starts_with('>') {
            let header = line[1..].trim();
            let mut parts = header.splitn(2, char::is_whitespace);
            let id = parts.next().unwrap_or("unnamed").to_string();
            let description = parts.next().unwrap_or("").trim().to_string();

            // Read exactly 4 rows (A, C, G, T).
            let mut raw_rows: Vec<Vec<f64>> = Vec::new();
            let expected_labels = ['A', 'C', 'G', 'T'];
            i += 1;

            for expected_label in &expected_labels {
                if i >= lines.len() {
                    return Err(CyaneaError::Parse(format!(
                        "JASPAR motif '{}': expected row '{}' but reached end of input",
                        id, expected_label
                    )));
                }
                let rline = lines[i].trim();

                // Extract the part between [ and ].
                let bracket_start = rline.find('[').ok_or_else(|| {
                    CyaneaError::Parse(format!(
                        "JASPAR motif '{}': missing '[' in row '{}'",
                        id, expected_label
                    ))
                })?;
                let bracket_end = rline.find(']').ok_or_else(|| {
                    CyaneaError::Parse(format!(
                        "JASPAR motif '{}': missing ']' in row '{}'",
                        id, expected_label
                    ))
                })?;

                let inner = &rline[bracket_start + 1..bracket_end];
                let vals: std::result::Result<Vec<f64>, _> = inner
                    .split_whitespace()
                    .map(|s| s.parse::<f64>())
                    .collect();
                match vals {
                    Ok(v) => raw_rows.push(v),
                    Err(e) => {
                        return Err(CyaneaError::Parse(format!(
                            "JASPAR motif '{}': failed to parse row '{}': {}",
                            id, expected_label, e
                        )));
                    }
                }
                i += 1;
            }

            // Validate all rows have the same length.
            let width = raw_rows[0].len();
            for (ri, row) in raw_rows.iter().enumerate() {
                if row.len() != width {
                    return Err(CyaneaError::Parse(format!(
                        "JASPAR motif '{}': row {} has {} values but expected {}",
                        id,
                        expected_labels[ri],
                        row.len(),
                        width
                    )));
                }
            }

            // Transpose: raw_rows[base][position] -> matrix[position][base].
            let mut matrix = Vec::with_capacity(width);
            for pos in 0..width {
                let col = [
                    raw_rows[0][pos],
                    raw_rows[1][pos],
                    raw_rows[2][pos],
                    raw_rows[3][pos],
                ];
                // Normalize to probabilities.
                let sum: f64 = col.iter().sum();
                if sum > 0.0 {
                    matrix.push([col[0] / sum, col[1] / sum, col[2] / sum, col[3] / sum]);
                } else {
                    matrix.push([0.25, 0.25, 0.25, 0.25]);
                }
            }

            let mut metadata = HashMap::new();
            if !description.is_empty() {
                metadata.insert("description".to_string(), description);
            }

            motifs.push(Motif {
                name: id,
                alphabet: MotifAlphabet::Dna,
                matrix,
                metadata,
            });
        } else {
            i += 1;
        }
    }

    Ok(motifs)
}

/// Write motifs in JASPAR format.
///
/// Outputs counts derived by scaling probabilities to sum to 100 per position.
pub fn write_jaspar(motifs: &[Motif]) -> String {
    let mut out = String::new();

    for motif in motifs {
        // Header line.
        let desc = motif.metadata.get("description").map(|s| s.as_str()).unwrap_or("");
        if desc.is_empty() {
            out.push_str(&format!(">{}\n", motif.name));
        } else {
            out.push_str(&format!(">{} {}\n", motif.name, desc));
        }

        let width = motif.matrix.len();
        let labels = ['A', 'C', 'G', 'T'];

        for (base_idx, label) in labels.iter().enumerate() {
            out.push_str(&format!("{} [", label));
            for pos in 0..width {
                let count = (motif.matrix[pos][base_idx] * 100.0).round() as i64;
                if pos > 0 {
                    out.push(' ');
                }
                out.push_str(&format!("{:3}", count));
            }
            out.push_str(" ]\n");
        }
    }

    out
}

// ---------------------------------------------------------------------------
// PWM comparison — Pearson correlation
// ---------------------------------------------------------------------------

/// Compute the similarity between two motifs using Pearson correlation.
///
/// For each possible ungapped alignment offset of two motifs, the function
/// computes the per-column Pearson correlation between the overlapping
/// portions and averages them. The maximum average correlation across all
/// offsets is returned.
///
/// Returns `0.0` if the motifs do not overlap (overlap < 1 column).
pub fn motif_similarity(a: &Motif, b: &Motif) -> f64 {
    let la = a.matrix.len() as isize;
    let lb = b.matrix.len() as isize;

    if la == 0 || lb == 0 {
        return 0.0;
    }

    let mut best = f64::NEG_INFINITY;

    // Slide b across a at every possible offset.
    // offset = start position of b relative to a.
    // Range: -(lb - 1) .. la - 1  (inclusive on both sides).
    for offset in (-(lb - 1))..la {
        // Determine overlap region.
        let col_start_a = offset.max(0) as usize;
        let col_start_b = (-offset).max(0) as usize;
        let col_end_a = (la).min(offset + lb) as usize;
        let overlap = col_end_a as isize - col_start_a as isize;

        if overlap < 1 {
            continue;
        }

        let n = overlap as usize;
        let mut corr_sum = 0.0;

        for k in 0..n {
            let col_a = &a.matrix[col_start_a + k];
            let col_b = &b.matrix[col_start_b + k];
            corr_sum += pearson4(col_a, col_b);
        }

        let avg = corr_sum / n as f64;
        if avg > best {
            best = avg;
        }
    }

    if best == f64::NEG_INFINITY {
        0.0
    } else {
        best
    }
}

/// Pearson correlation between two 4-element arrays.
fn pearson4(x: &[f64; 4], y: &[f64; 4]) -> f64 {
    let n = 4.0;
    let sum_x: f64 = x.iter().sum();
    let sum_y: f64 = y.iter().sum();
    let mean_x = sum_x / n;
    let mean_y = sum_y / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;

    for i in 0..4 {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom < 1e-15 {
        // If both are constant (e.g. uniform), they are "perfectly similar".
        if var_x < 1e-15 && var_y < 1e-15 {
            return 1.0;
        }
        return 0.0;
    }

    cov / denom
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_motif(name: &str, matrix: Vec<[f64; 4]>) -> Motif {
        Motif {
            name: name.to_string(),
            alphabet: MotifAlphabet::Dna,
            matrix,
            metadata: HashMap::new(),
        }
    }

    // 1. MEME round-trip
    #[test]
    fn meme_round_trip() {
        let motif = sample_motif(
            "test1",
            vec![
                [0.9, 0.03, 0.04, 0.03],
                [0.1, 0.7, 0.1, 0.1],
                [0.05, 0.05, 0.85, 0.05],
            ],
        );
        let written = write_meme(&[motif]);
        let parsed = parse_meme(&written).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "test1");
        assert_eq!(parsed[0].matrix.len(), 3);
        for (orig, read) in [
            [0.9, 0.03, 0.04, 0.03],
            [0.1, 0.7, 0.1, 0.1],
            [0.05, 0.05, 0.85, 0.05],
        ]
        .iter()
        .zip(parsed[0].matrix.iter())
        {
            for j in 0..4 {
                assert!((orig[j] - read[j]).abs() < 1e-4);
            }
        }
    }

    // 2. TRANSFAC round-trip
    #[test]
    fn transfac_round_trip() {
        let motif = sample_motif(
            "TF_MOTIF",
            vec![
                [0.5, 0.2, 0.2, 0.1],
                [0.1, 0.6, 0.2, 0.1],
            ],
        );
        let written = write_transfac(&[motif]);
        let parsed = parse_transfac(&written).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "TF_MOTIF");
        assert_eq!(parsed[0].matrix.len(), 2);
        for j in 0..4 {
            assert!((parsed[0].matrix[0][j] - [0.5, 0.2, 0.2, 0.1][j]).abs() < 1e-3);
        }
    }

    // 3. JASPAR round-trip
    #[test]
    fn jaspar_round_trip() {
        let motif = sample_motif(
            "MA0001.1",
            vec![
                [0.25, 0.25, 0.25, 0.25],
                [0.8, 0.1, 0.05, 0.05],
                [0.1, 0.1, 0.7, 0.1],
            ],
        );
        let written = write_jaspar(&[motif]);
        let parsed = parse_jaspar(&written).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "MA0001.1");
        assert_eq!(parsed[0].matrix.len(), 3);
        // JASPAR round-trips through integer counts (scaled x100), so check
        // that values are within 0.02 tolerance.
        for (orig, read) in [
            [0.25, 0.25, 0.25, 0.25],
            [0.8, 0.1, 0.05, 0.05],
            [0.1, 0.1, 0.7, 0.1],
        ]
        .iter()
        .zip(parsed[0].matrix.iter())
        {
            for j in 0..4 {
                assert!(
                    (orig[j] - read[j]).abs() < 0.02,
                    "expected ~{} got {} at col {}",
                    orig[j],
                    read[j],
                    j
                );
            }
        }
    }

    // 4. motif_similarity of identical motifs = 1.0
    #[test]
    fn similarity_identical() {
        let m = sample_motif(
            "id",
            vec![
                [0.9, 0.03, 0.04, 0.03],
                [0.05, 0.8, 0.1, 0.05],
                [0.1, 0.1, 0.7, 0.1],
            ],
        );
        let sim = motif_similarity(&m, &m);
        assert!(
            (sim - 1.0).abs() < 1e-10,
            "identical motifs should have similarity 1.0, got {}",
            sim
        );
    }

    // 5. motif_similarity of different motifs < 1.0
    #[test]
    fn similarity_different() {
        let a = sample_motif(
            "a",
            vec![
                [0.9, 0.03, 0.04, 0.03],
                [0.03, 0.9, 0.04, 0.03],
                [0.03, 0.03, 0.9, 0.04],
            ],
        );
        let b = sample_motif(
            "b",
            vec![
                [0.1, 0.2, 0.3, 0.4],
                [0.4, 0.3, 0.2, 0.1],
                [0.1, 0.4, 0.2, 0.3],
            ],
        );
        let sim = motif_similarity(&a, &b);
        assert!(
            sim < 1.0,
            "different motifs should have similarity < 1.0, got {}",
            sim
        );
    }

    // 6. Multi-motif MEME parsing
    #[test]
    fn multi_motif_meme() {
        let input = "\
MEME version 5

ALPHABET= ACGT

strands: + -

MOTIF first
letter-probability matrix: alength= 4 w= 2
 0.25  0.25  0.25  0.25
 0.90  0.03  0.04  0.03

MOTIF second
letter-probability matrix: alength= 4 w= 3
 0.10  0.10  0.70  0.10
 0.25  0.25  0.25  0.25
 0.05  0.80  0.10  0.05
";
        let motifs = parse_meme(input).unwrap();
        assert_eq!(motifs.len(), 2);
        assert_eq!(motifs[0].name, "first");
        assert_eq!(motifs[0].matrix.len(), 2);
        assert_eq!(motifs[1].name, "second");
        assert_eq!(motifs[1].matrix.len(), 3);
    }

    // 7. Multi-motif TRANSFAC parsing
    #[test]
    fn multi_motif_transfac() {
        let input = "\
ID  motif_A
P0      A      C      G      T
01   0.2500   0.2500   0.2500   0.2500
XX
//
ID  motif_B
P0      A      C      G      T
01   0.9000   0.0300   0.0400   0.0300
02   0.1000   0.7000   0.1000   0.1000
XX
//
";
        let motifs = parse_transfac(input).unwrap();
        assert_eq!(motifs.len(), 2);
        assert_eq!(motifs[0].name, "motif_A");
        assert_eq!(motifs[0].matrix.len(), 1);
        assert_eq!(motifs[1].name, "motif_B");
        assert_eq!(motifs[1].matrix.len(), 2);
    }

    // 8. Empty input returns empty vec
    #[test]
    fn empty_input_returns_empty() {
        assert!(parse_meme("").unwrap().is_empty());
        assert!(parse_transfac("").unwrap().is_empty());
        assert!(parse_jaspar("").unwrap().is_empty());
        assert!(parse_meme("   \n  \n").unwrap().is_empty());
    }

    // 9. Malformed MEME returns error
    #[test]
    fn malformed_meme_error() {
        // Missing MEME version header.
        let input = "MOTIF test\nletter-probability matrix: alength= 4 w= 2\n0.25 0.25 0.25 0.25\n";
        assert!(parse_meme(input).is_err());
    }

    // 10. Malformed JASPAR returns error
    #[test]
    fn malformed_jaspar_error() {
        // Missing brackets.
        let input = ">MA0001.1 test\nA  1  2  3  4\nC  5  6  7  8\nG  9 10 11 12\nT 13 14 15 16\n";
        assert!(parse_jaspar(input).is_err());
    }
}
