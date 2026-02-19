//! PIR/NBRF sequence format parser and writer.
//!
//! Parses the PIR (Protein Information Resource) / NBRF format used
//! by sequence databases. Each record has a header line starting with
//! `>` followed by a two-character entry type code and a semicolon-separated
//! name, a description line, and sequence data terminated by `*`.

use cyanea_core::{CyaneaError, Result};

/// A parsed PIR/NBRF record.
#[derive(Debug, Clone)]
pub struct PirRecord {
    /// Entry type code (e.g. `P1` for protein, `F1` for fragment, `DL` for DNA).
    pub entry_type: String,
    /// Sequence name/identifier.
    pub name: String,
    /// Description from the second header line.
    pub description: String,
    /// Sequence data (without the terminating `*`).
    pub sequence: String,
}

/// Parse one or more PIR records from a string.
///
/// Each record begins with `>type;name`, followed by a description line,
/// then sequence data terminated by `*`. Returns an error if the sequence
/// terminator `*` is missing.
///
/// # Examples
///
/// ```
/// # use cyanea_io::pir::parse_pir;
/// let input = ">P1;CRAB_ANAPL\nAlpha crystallin\nMASHAKE*\n";
/// let records = parse_pir(input).unwrap();
/// assert_eq!(records.len(), 1);
/// assert_eq!(records[0].sequence, "MASHAKE");
/// ```
pub fn parse_pir(input: &str) -> Result<Vec<PirRecord>> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    let mut records = Vec::new();
    let mut lines = input.lines().peekable();

    while let Some(line) = lines.next() {
        let line = line.trim();

        // Skip empty lines between records
        if line.is_empty() {
            continue;
        }

        // Header line: >type;name
        if !line.starts_with('>') {
            continue;
        }

        let header = &line[1..]; // strip '>'
        let (entry_type, name) = match header.split_once(';') {
            Some((et, n)) => (et.trim().to_string(), n.trim().to_string()),
            None => {
                return Err(CyaneaError::Parse(format!(
                    "invalid PIR header (missing ';'): '{}'",
                    line
                )));
            }
        };

        // Description line
        let description = match lines.next() {
            Some(desc) => desc.trim().to_string(),
            None => {
                return Err(CyaneaError::Parse(format!(
                    "unexpected end of input after PIR header for '{}'",
                    name
                )));
            }
        };

        // Sequence lines until we encounter '*'
        let mut sequence = String::new();
        let mut found_terminator = false;

        while let Some(seq_line) = lines.next() {
            let seq_line = seq_line.trim();
            if seq_line.is_empty() {
                continue;
            }

            // Check if this line contains the terminator
            if let Some(pos) = seq_line.find('*') {
                // Take everything before the *
                let fragment: String = seq_line[..pos]
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                sequence.push_str(&fragment);
                found_terminator = true;
                break;
            } else {
                let fragment: String = seq_line
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                sequence.push_str(&fragment);
            }
        }

        if !found_terminator {
            return Err(CyaneaError::Parse(format!(
                "missing sequence terminator '*' for PIR record '{}'",
                name
            )));
        }

        records.push(PirRecord {
            entry_type,
            name,
            description,
            sequence,
        });
    }

    Ok(records)
}

/// Write PIR records as a string.
///
/// Each record is formatted with `>type;name`, description, sequence
/// (60 characters per line), and a `*` terminator.
///
/// # Examples
///
/// ```
/// # use cyanea_io::pir::{PirRecord, write_pir};
/// let rec = PirRecord {
///     entry_type: "P1".to_string(),
///     name: "TEST".to_string(),
///     description: "Test protein".to_string(),
///     sequence: "MASHAKE".to_string(),
/// };
/// let output = write_pir(&[rec]);
/// assert!(output.contains(">P1;TEST"));
/// assert!(output.contains("MASHAKE*"));
/// ```
pub fn write_pir(records: &[PirRecord]) -> String {
    const LINE_WIDTH: usize = 60;

    let mut out = String::new();

    for rec in records {
        // Header
        out.push_str(&format!(">{};{}\n", rec.entry_type, rec.name));
        // Description
        out.push_str(&rec.description);
        out.push('\n');

        // Sequence in lines of LINE_WIDTH, with * at the end
        let seq = &rec.sequence;
        let mut pos = 0;
        while pos < seq.len() {
            let end = std::cmp::min(pos + LINE_WIDTH, seq.len());
            out.push_str(&seq[pos..end]);
            if end == seq.len() {
                out.push('*');
            }
            out.push('\n');
            pos = end;
        }

        if seq.is_empty() {
            out.push_str("*\n");
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pir_round_trip() {
        let rec = PirRecord {
            entry_type: "P1".to_string(),
            name: "CRAB_ANAPL".to_string(),
            description: "Alpha A crystallin - Atlantic crab".to_string(),
            sequence: "MFNIFFHDFKAPLNASSQVIFEKGESK".to_string(),
        };

        let written = write_pir(&[rec]);
        let parsed = parse_pir(&written).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].entry_type, "P1");
        assert_eq!(parsed[0].name, "CRAB_ANAPL");
        assert_eq!(parsed[0].description, "Alpha A crystallin - Atlantic crab");
        assert_eq!(parsed[0].sequence, "MFNIFFHDFKAPLNASSQVIFEKGESK");
    }

    #[test]
    fn pir_multi_record() {
        let input = "\
>P1;CRAB_ANAPL
Alpha A crystallin - Atlantic crab
MFNIFFHDFK*
>P1;CRAB_BOVIN
Alpha crystallin - Bovine
MDIAIHHPWI*
>DL;DNA_TEST
Test DNA sequence
ACGTACGT*
";
        let records = parse_pir(input).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].name, "CRAB_ANAPL");
        assert_eq!(records[0].sequence, "MFNIFFHDFK");
        assert_eq!(records[1].name, "CRAB_BOVIN");
        assert_eq!(records[1].sequence, "MDIAIHHPWI");
        assert_eq!(records[2].entry_type, "DL");
        assert_eq!(records[2].name, "DNA_TEST");
        assert_eq!(records[2].sequence, "ACGTACGT");
    }

    #[test]
    fn pir_various_entry_types() {
        let input = "\
>F1;FRAG1
Protein fragment
ACDEF*
>DL;DNA1
Linear DNA
ACGT*
>N1;RNA1
RNA sequence
ACGU*
";
        let records = parse_pir(input).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].entry_type, "F1");
        assert_eq!(records[1].entry_type, "DL");
        assert_eq!(records[2].entry_type, "N1");
    }

    #[test]
    fn pir_missing_terminator_error() {
        let input = ">P1;TEST\nTest protein\nMASHAKE\n";
        let result = parse_pir(input);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("terminator"));
    }

    #[test]
    fn pir_empty_input() {
        let records = parse_pir("").unwrap();
        assert!(records.is_empty());

        let records = parse_pir("  \n\n  ").unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn pir_multiline_sequence() {
        let input = "\
>P1;LONGSEQ
A long protein
MFNIFFHDFKAPLNASSQVIFEKGESKFVNQDVKEMEEQTQAFNKHSSGFGRNYQD
ADVDVAYGFSSGGPTIAAGK*
";
        let records = parse_pir(input).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence.len(), 76);
        assert!(records[0].sequence.starts_with("MFNIFFHDFK"));
        assert!(records[0].sequence.ends_with("AAGK"));
    }
}
