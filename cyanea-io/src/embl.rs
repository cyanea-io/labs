//! EMBL flat file format parser and writer.
//!
//! Parses EMBL/ENA format files (`.embl`, `.dat`) containing one or more
//! records separated by `//`. Each record uses two-letter line-type
//! prefixes: `ID`, `AC`, `DE`, `FT`, `SQ`, and sequence data lines.

use cyanea_core::Result;

/// A parsed EMBL record.
#[derive(Debug, Clone)]
pub struct EmblRecord {
    /// Entry name from the ID line.
    pub id: String,
    /// Primary accession number from the AC line.
    pub accession: String,
    /// Description from DE lines (concatenated).
    pub description: String,
    /// Raw sequence with only alphabetic characters.
    pub sequence: String,
    /// Feature table entries as `(key, location)` pairs from FT lines.
    pub features: Vec<(String, String)>,
}

/// Parse one or more EMBL records from a string.
///
/// Records are terminated by `//`. Returns an empty `Vec` for empty input.
///
/// # Examples
///
/// ```
/// # use cyanea_io::embl::parse_embl;
/// let input = "ID   TEST; SV 1; linear; DNA; STD; HUM; 10 BP.\n\
///              AC   X00001;\n\
///              DE   Test sequence.\n\
///              SQ   Sequence 10 BP;\n\
///              \x20    acgtacgtac                                                          10\n\
///              //\n";
/// let records = parse_embl(input).unwrap();
/// assert_eq!(records.len(), 1);
/// assert_eq!(records[0].sequence, "acgtacgtac");
/// ```
pub fn parse_embl(input: &str) -> Result<Vec<EmblRecord>> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    let mut records = Vec::new();
    let mut builder = RecordBuilder::new();
    let mut in_sequence = false;

    for line in input.lines() {
        // Record terminator
        if line.starts_with("//") {
            if builder.has_data() {
                records.push(builder.build());
            }
            builder = RecordBuilder::new();
            in_sequence = false;
            continue;
        }

        if in_sequence {
            // Sequence data lines: spaces + bases + optional trailing number
            // Extract only alphabetic characters
            for ch in line.chars() {
                if ch.is_ascii_alphabetic() {
                    builder.sequence.push(ch);
                }
            }
            continue;
        }

        // ID line: "ID   name; ..."
        if line.starts_with("ID   ") || line.starts_with("ID\t") {
            let rest = &line[5..];
            let rest = rest.trim_start();
            // Name is the first token (before ';' or whitespace)
            let name = rest
                .split(|c: char| c == ';' || c.is_whitespace())
                .next()
                .unwrap_or("")
                .trim();
            builder.id = name.to_string();
            continue;
        }

        // AC line: "AC   accession;"
        if line.starts_with("AC   ") || line.starts_with("AC\t") {
            let rest = &line[5..];
            let rest = rest.trim();
            let acc = rest.trim_end_matches(';').trim();
            if builder.accession.is_empty() {
                builder.accession = acc.to_string();
            }
            continue;
        }

        // DE line: "DE   description"
        if line.starts_with("DE   ") || line.starts_with("DE\t") {
            let rest = &line[5..];
            let rest = rest.trim();
            if !builder.description.is_empty() {
                builder.description.push(' ');
            }
            builder.description.push_str(rest);
            continue;
        }

        // FT line: "FT   key             location" or "FT                   /qualifier"
        if line.starts_with("FT   ") || line.starts_with("FT\t") {
            let rest = &line[5..];
            // A new feature starts with a non-space character after FT prefix
            let trimmed = rest.trim_start();
            if !trimmed.is_empty() && !trimmed.starts_with('/') {
                // Try to split key and location
                let parts: Vec<&str> = trimmed.splitn(2, char::is_whitespace).collect();
                if parts.len() >= 2 {
                    let key = parts[0].trim().to_string();
                    let location = parts[1].trim().to_string();
                    builder.features.push((key, location));
                } else if parts.len() == 1 {
                    builder.features.push((parts[0].trim().to_string(), String::new()));
                }
            }
            // Qualifier lines (starting with /) are skipped for simplicity
            continue;
        }

        // SQ line starts the sequence block
        if line.starts_with("SQ   ") || line.starts_with("SQ\t") || line == "SQ" {
            in_sequence = true;
            continue;
        }

        // Other lines (OS, OC, OX, RN, etc.) are skipped
    }

    // Handle records without trailing //
    if builder.has_data() {
        records.push(builder.build());
    }

    Ok(records)
}

/// Write EMBL records as a string.
///
/// Each record is terminated by `//`.
///
/// # Examples
///
/// ```
/// # use cyanea_io::embl::{EmblRecord, write_embl};
/// let rec = EmblRecord {
///     id: "TEST".to_string(),
///     accession: "X00001".to_string(),
///     description: "Test sequence.".to_string(),
///     sequence: "acgtacgtac".to_string(),
///     features: vec![],
/// };
/// let output = write_embl(&[rec]);
/// assert!(output.contains("ID   TEST"));
/// ```
pub fn write_embl(records: &[EmblRecord]) -> String {
    let mut out = String::new();

    for rec in records {
        // ID line
        out.push_str(&format!("ID   {}; SV 1; linear; DNA; STD; UNC; {} BP.\n",
            rec.id, rec.sequence.len()));

        // AC line
        out.push_str(&format!("AC   {};\n", rec.accession));

        // DE line
        out.push_str(&format!("DE   {}\n", rec.description));

        // FT lines
        for (key, location) in &rec.features {
            out.push_str(&format!("FT   {:<16}{}\n", key, location));
        }

        // SQ header
        out.push_str(&format!("SQ   Sequence {} BP;\n", rec.sequence.len()));

        // Sequence data: 60 chars per line in groups of 10, with trailing count
        let seq_bytes = rec.sequence.as_bytes();
        let mut pos = 0;
        while pos < seq_bytes.len() {
            out.push_str("     ");
            let line_end = std::cmp::min(pos + 60, seq_bytes.len());
            let mut col = 0;
            for i in pos..line_end {
                out.push(seq_bytes[i] as char);
                col += 1;
                if col % 10 == 0 && i + 1 < line_end {
                    out.push(' ');
                }
            }
            // Pad to align the count
            let chars_written = (line_end - pos) + ((line_end - pos).saturating_sub(1)) / 10;
            let target_width = 60 + 5; // 60 bases + 5 spaces between groups
            for _ in chars_written..target_width {
                out.push(' ');
            }
            out.push_str(&format!(" {}\n", line_end));
            pos = line_end;
        }

        out.push_str("//\n");
    }

    out
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

struct RecordBuilder {
    id: String,
    accession: String,
    description: String,
    sequence: String,
    features: Vec<(String, String)>,
}

impl RecordBuilder {
    fn new() -> Self {
        Self {
            id: String::new(),
            accession: String::new(),
            description: String::new(),
            sequence: String::new(),
            features: Vec::new(),
        }
    }

    fn has_data(&self) -> bool {
        !self.id.is_empty() || !self.sequence.is_empty()
    }

    fn build(self) -> EmblRecord {
        EmblRecord {
            id: self.id,
            accession: self.accession,
            description: self.description,
            sequence: self.sequence,
            features: self.features,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embl_round_trip() {
        let rec = EmblRecord {
            id: "HSBGLOBIN".to_string(),
            accession: "V00497".to_string(),
            description: "Human beta-globin gene.".to_string(),
            sequence: "acgtacgtacgtacgt".to_string(),
            features: vec![
                ("gene".to_string(), "1..16".to_string()),
            ],
        };

        let written = write_embl(&[rec]);
        let parsed = parse_embl(&written).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].id, "HSBGLOBIN");
        assert_eq!(parsed[0].accession, "V00497");
        assert_eq!(parsed[0].description, "Human beta-globin gene.");
        assert_eq!(parsed[0].sequence, "acgtacgtacgtacgt");
        assert_eq!(parsed[0].features.len(), 1);
        assert_eq!(parsed[0].features[0].0, "gene");
        assert_eq!(parsed[0].features[0].1, "1..16");
    }

    #[test]
    fn embl_multi_record() {
        let input = "\
ID   REC1; SV 1; linear; DNA; STD; HUM; 10 BP.
AC   X00001;
DE   First record.
SQ   Sequence 10 BP;
     acgtacgtac                                                          10
//
ID   REC2; SV 1; linear; DNA; STD; HUM; 8 BP.
AC   X00002;
DE   Second record.
SQ   Sequence 8 BP;
     tgcatgca                                                            8
//
";
        let records = parse_embl(input).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "REC1");
        assert_eq!(records[0].accession, "X00001");
        assert_eq!(records[0].sequence, "acgtacgtac");
        assert_eq!(records[1].id, "REC2");
        assert_eq!(records[1].accession, "X00002");
        assert_eq!(records[1].sequence, "tgcatgca");
    }

    #[test]
    fn embl_feature_extraction() {
        let input = "\
ID   FEAT1; SV 1; linear; DNA; STD; HUM; 20 BP.
AC   Y00001;
DE   Feature test.
FT   gene            1..20
FT                   /gene=\"TP53\"
FT   CDS             join(1..10,15..20)
FT                   /protein_id=\"AAA001.1\"
SQ   Sequence 20 BP;
     acgtacgtac acgtacgtac                                               20
//
";
        let records = parse_embl(input).unwrap();
        assert_eq!(records[0].features.len(), 2);
        assert_eq!(records[0].features[0].0, "gene");
        assert_eq!(records[0].features[0].1, "1..20");
        assert_eq!(records[0].features[1].0, "CDS");
        assert_eq!(records[0].features[1].1, "join(1..10,15..20)");
    }

    #[test]
    fn embl_empty_input() {
        let records = parse_embl("").unwrap();
        assert!(records.is_empty());

        let records = parse_embl("   \n\n  ").unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn embl_sequence_whitespace_extraction() {
        // Sequence data lines have spaces and trailing numbers that must be stripped
        let input = "\
ID   SEQTEST; SV 1; linear; DNA; STD; UNC; 30 BP.
AC   Z99999;
DE   Whitespace test.
SQ   Sequence 30 BP;
     aaaaaaaaaa cccccccccc gggggggggg                                    30
//
";
        let records = parse_embl(input).unwrap();
        assert_eq!(records[0].sequence, "aaaaaaaaaaccccccccccgggggggggg");
        assert_eq!(records[0].sequence.len(), 30);
    }
}
