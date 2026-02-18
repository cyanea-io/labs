//! GenBank flat file parser.
//!
//! Parses GenBank format files (`.gb`, `.gbk`) containing one or more records
//! separated by `//`. Each record includes LOCUS, DEFINITION, ACCESSION,
//! FEATURES table, and ORIGIN sequence sections.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

/// A parsed GenBank record.
#[derive(Debug, Clone)]
pub struct GenbankRecord {
    /// Locus name from the LOCUS line.
    pub locus: String,
    /// Definition/description line.
    pub definition: String,
    /// Primary accession number.
    pub accession: String,
    /// Version (accession.version).
    pub version: Option<String>,
    /// Organism name from the SOURCE/ORGANISM section.
    pub organism: Option<String>,
    /// Raw DNA/RNA/protein sequence.
    pub sequence: Vec<u8>,
    /// Feature table entries.
    pub features: Vec<GenbankFeature>,
}

/// A feature from the GenBank FEATURES table.
#[derive(Debug, Clone)]
pub struct GenbankFeature {
    /// Feature type (e.g. "gene", "CDS", "mRNA").
    pub feature_type: String,
    /// Raw location string (e.g. "join(100..200,300..400)", "complement(1..100)").
    pub location: String,
    /// Qualifier key-value pairs from `/key="value"` lines.
    pub qualifiers: Vec<(String, String)>,
}

/// Summary statistics for a GenBank file.
#[derive(Debug, Clone)]
pub struct GenbankStats {
    /// Number of records in the file.
    pub record_count: u64,
    /// Total number of sequence bases/residues.
    pub total_bases: u64,
    /// Count of each feature type.
    pub feature_counts: HashMap<String, u64>,
}

/// Parse a GenBank flat file and return all records.
///
/// Multi-record files (separated by `//`) are fully supported.
pub fn parse_genbank(path: impl AsRef<Path>) -> Result<Vec<GenbankRecord>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut builder = RecordBuilder::new();
    let mut in_features = false;
    let mut in_origin = false;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;

        // Record terminator
        if line.starts_with("//") {
            if builder.has_data() {
                records.push(builder.build());
            }
            builder = RecordBuilder::new();
            in_features = false;
            in_origin = false;
            continue;
        }

        if in_origin {
            // ORIGIN sequence lines: "   1 atgcatgcat gcatgcatgc ..."
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            for part in trimmed.split_whitespace() {
                // Skip the line number at the start
                if part.chars().all(|c| c.is_ascii_digit()) {
                    continue;
                }
                for ch in part.chars() {
                    if ch.is_ascii_alphabetic() {
                        builder.sequence.push(ch.to_ascii_uppercase() as u8);
                    }
                }
            }
            continue;
        }

        if line.starts_with("ORIGIN") {
            in_origin = true;
            in_features = false;
            continue;
        }

        if in_features {
            parse_feature_line(&line, &mut builder);
            continue;
        }

        if line.starts_with("FEATURES") {
            in_features = true;
            continue;
        }

        // Header sections
        if line.starts_with("LOCUS") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                builder.locus = parts[1].to_string();
            }
        } else if line.starts_with("DEFINITION") {
            builder.definition = line[12..].trim().to_string();
        } else if line.starts_with("            ") && !builder.definition.is_empty() && builder.accession.is_empty() {
            // Continuation of DEFINITION
            builder.definition.push(' ');
            builder.definition.push_str(line.trim());
        } else if line.starts_with("ACCESSION") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                builder.accession = parts[1].to_string();
            }
        } else if line.starts_with("VERSION") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                builder.version = Some(parts[1].to_string());
            }
        } else if line.trim_start().starts_with("ORGANISM") {
            let trimmed = line.trim_start();
            if trimmed.len() > 8 {
                builder.organism = Some(trimmed[8..].trim().to_string());
            }
        }
    }

    // Handle file without trailing //
    if builder.has_data() {
        records.push(builder.build());
    }

    if records.is_empty() {
        return Err(CyaneaError::Parse(format!(
            "{}: no GenBank records found",
            path.display()
        )));
    }

    Ok(records)
}

/// Compute summary statistics from a GenBank file.
pub fn genbank_stats(path: impl AsRef<Path>) -> Result<GenbankStats> {
    let records = parse_genbank(path)?;
    let mut total_bases: u64 = 0;
    let mut feature_counts: HashMap<String, u64> = HashMap::new();

    for rec in &records {
        total_bases += rec.sequence.len() as u64;
        for feat in &rec.features {
            *feature_counts.entry(feat.feature_type.clone()).or_insert(0) += 1;
        }
    }

    Ok(GenbankStats {
        record_count: records.len() as u64,
        total_bases,
        feature_counts,
    })
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

struct RecordBuilder {
    locus: String,
    definition: String,
    accession: String,
    version: Option<String>,
    organism: Option<String>,
    sequence: Vec<u8>,
    features: Vec<GenbankFeature>,
    current_feature: Option<FeatureBuilder>,
}

struct FeatureBuilder {
    feature_type: String,
    location: String,
    qualifiers: Vec<(String, String)>,
    current_qualifier_key: Option<String>,
    current_qualifier_value: String,
}

impl RecordBuilder {
    fn new() -> Self {
        Self {
            locus: String::new(),
            definition: String::new(),
            accession: String::new(),
            version: None,
            organism: None,
            sequence: Vec::new(),
            features: Vec::new(),
            current_feature: None,
        }
    }

    fn has_data(&self) -> bool {
        !self.locus.is_empty() || !self.sequence.is_empty()
    }

    fn finish_feature(&mut self) {
        if let Some(mut fb) = self.current_feature.take() {
            fb.finish_qualifier();
            self.features.push(GenbankFeature {
                feature_type: fb.feature_type,
                location: fb.location,
                qualifiers: fb.qualifiers,
            });
        }
    }

    fn build(mut self) -> GenbankRecord {
        self.finish_feature();
        GenbankRecord {
            locus: self.locus,
            definition: self.definition,
            accession: self.accession,
            version: self.version,
            organism: self.organism,
            sequence: self.sequence,
            features: self.features,
        }
    }
}

impl FeatureBuilder {
    fn finish_qualifier(&mut self) {
        if let Some(key) = self.current_qualifier_key.take() {
            let value = self.current_qualifier_value.trim().to_string();
            // Strip surrounding quotes
            let value = value
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .unwrap_or(&value)
                .to_string();
            self.qualifiers.push((key, value));
            self.current_qualifier_value.clear();
        }
    }
}

/// Parse a line within the FEATURES section.
///
/// Feature key lines start at column 5 (after 5 spaces),
/// qualifier/continuation lines start at column 21 (after 21 spaces).
fn parse_feature_line(line: &str, builder: &mut RecordBuilder) {
    if line.len() < 6 {
        return;
    }

    // Check if this is a new feature key (starts at column 5, non-space at column 5)
    // Feature keys: "     gene            complement(1..100)"
    let prefix = &line[..std::cmp::min(21, line.len())];
    let trimmed = prefix.trim_start();

    if !line.starts_with("                     ") && !trimmed.is_empty() {
        // New feature line
        builder.finish_feature();
        let parts: Vec<&str> = line.trim().splitn(2, char::is_whitespace).collect();
        if parts.len() >= 2 {
            builder.current_feature = Some(FeatureBuilder {
                feature_type: parts[0].to_string(),
                location: parts[1].trim().to_string(),
                qualifiers: Vec::new(),
                current_qualifier_key: None,
                current_qualifier_value: String::new(),
            });
        }
    } else if line.starts_with("                     ") {
        // Qualifier or continuation line (21 spaces)
        let content = line[21..].trim_end();
        if let Some(ref mut feat) = builder.current_feature {
            if let Some(stripped) = content.strip_prefix('/') {
                // New qualifier: /key="value" or /key=value or /key
                feat.finish_qualifier();
                if let Some((key, value)) = stripped.split_once('=') {
                    feat.current_qualifier_key = Some(key.to_string());
                    feat.current_qualifier_value = value.to_string();
                    // Check if the value is complete (ends with quote)
                    if value.starts_with('"') && value.ends_with('"') && value.len() > 1 {
                        feat.finish_qualifier();
                    }
                } else {
                    // Flag qualifier with no value
                    feat.qualifiers.push((stripped.to_string(), String::new()));
                }
            } else {
                // Continuation of previous qualifier value
                if feat.current_qualifier_key.is_some() {
                    feat.current_qualifier_value.push(' ');
                    feat.current_qualifier_value.push_str(content);
                    // Check if value is now complete
                    if feat.current_qualifier_value.starts_with('"')
                        && feat.current_qualifier_value.ends_with('"')
                    {
                        feat.finish_qualifier();
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_genbank(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".gb").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn genbank_parse_simple() {
        let file = write_genbank(
            "LOCUS       AB000001             500 bp    DNA     linear   PRI 01-JAN-2020\n\
             DEFINITION  Homo sapiens test gene.\n\
             ACCESSION   AB000001\n\
             VERSION     AB000001.1\n\
             SOURCE      Homo sapiens (human)\n\
               ORGANISM  Homo sapiens\n\
             FEATURES             Location/Qualifiers\n\
             ORIGIN\n\
             \x20       1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat\n\
             \x20      51 atgcatgcat gcatgcatgc\n\
             //\n",
        );

        let records = parse_genbank(file.path()).unwrap();
        assert_eq!(records.len(), 1);

        let rec = &records[0];
        assert_eq!(rec.locus, "AB000001");
        assert_eq!(rec.definition, "Homo sapiens test gene.");
        assert_eq!(rec.accession, "AB000001");
        assert_eq!(rec.version, Some("AB000001.1".to_string()));
        assert_eq!(rec.organism, Some("Homo sapiens".to_string()));
        assert_eq!(rec.sequence.len(), 70);
        assert_eq!(&rec.sequence[..4], b"ATGC");
    }

    #[test]
    fn genbank_multi_record() {
        let file = write_genbank(
            "LOCUS       REC1                  30 bp    DNA     linear   PRI\n\
             DEFINITION  Record one.\n\
             ACCESSION   REC1\n\
             FEATURES             Location/Qualifiers\n\
             ORIGIN\n\
             \x20       1 atgcatgcat gcatgcatgc atgcatgcat\n\
             //\n\
             LOCUS       REC2                  20 bp    DNA     linear   PRI\n\
             DEFINITION  Record two.\n\
             ACCESSION   REC2\n\
             FEATURES             Location/Qualifiers\n\
             ORIGIN\n\
             \x20       1 gggggggggg aaaaaaaaaa\n\
             //\n",
        );

        let records = parse_genbank(file.path()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].locus, "REC1");
        assert_eq!(records[0].sequence.len(), 30);
        assert_eq!(records[1].locus, "REC2");
        assert_eq!(records[1].sequence.len(), 20);
    }

    #[test]
    fn genbank_features() {
        let file = write_genbank(
            "LOCUS       TEST                 100 bp    DNA     linear   PRI\n\
             DEFINITION  Test record.\n\
             ACCESSION   TEST001\n\
             FEATURES             Location/Qualifiers\n\
             \x20    gene            1..100\n\
             \x20                    /gene=\"TP53\"\n\
             \x20    CDS             join(1..50,60..100)\n\
             \x20                    /gene=\"TP53\"\n\
             \x20                    /protein_id=\"AAA001.1\"\n\
             \x20    mRNA            1..100\n\
             \x20                    /gene=\"TP53\"\n\
             ORIGIN\n\
             \x20       1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat\n\
             \x20      51 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat\n\
             //\n",
        );

        let records = parse_genbank(file.path()).unwrap();
        assert_eq!(records[0].features.len(), 3);
        assert_eq!(records[0].features[0].feature_type, "gene");
        assert_eq!(records[0].features[0].location, "1..100");
        assert_eq!(records[0].features[1].feature_type, "CDS");
        assert_eq!(records[0].features[1].location, "join(1..50,60..100)");
        assert_eq!(records[0].features[2].feature_type, "mRNA");
    }

    #[test]
    fn genbank_qualifiers() {
        let file = write_genbank(
            "LOCUS       TEST                  10 bp    DNA     linear   PRI\n\
             DEFINITION  Test.\n\
             ACCESSION   TEST001\n\
             FEATURES             Location/Qualifiers\n\
             \x20    gene            1..10\n\
             \x20                    /gene=\"TP53\"\n\
             \x20                    /note=\"tumor protein p53\"\n\
             ORIGIN\n\
             \x20       1 atgcatgcat\n\
             //\n",
        );

        let records = parse_genbank(file.path()).unwrap();
        let feat = &records[0].features[0];
        assert_eq!(feat.qualifiers.len(), 2);
        assert_eq!(feat.qualifiers[0], ("gene".to_string(), "TP53".to_string()));
        assert_eq!(
            feat.qualifiers[1],
            ("note".to_string(), "tumor protein p53".to_string())
        );
    }

    #[test]
    fn genbank_sequence() {
        let file = write_genbank(
            "LOCUS       SEQ1                  20 bp    DNA     linear\n\
             DEFINITION  Sequence test.\n\
             ACCESSION   SEQ1\n\
             FEATURES             Location/Qualifiers\n\
             ORIGIN\n\
             \x20       1 acgtacgtac gtacgtacgt\n\
             //\n",
        );

        let records = parse_genbank(file.path()).unwrap();
        assert_eq!(records[0].sequence, b"ACGTACGTACGTACGTACGT");
        assert_eq!(records[0].sequence.len(), 20);
    }

    #[test]
    fn genbank_stats_computed() {
        let file = write_genbank(
            "LOCUS       R1                    10 bp    DNA     linear\n\
             DEFINITION  Record 1.\n\
             ACCESSION   R1\n\
             FEATURES             Location/Qualifiers\n\
             \x20    gene            1..10\n\
             \x20                    /gene=\"A\"\n\
             \x20    CDS             1..10\n\
             \x20                    /gene=\"A\"\n\
             ORIGIN\n\
             \x20       1 atgcatgcat\n\
             //\n\
             LOCUS       R2                    20 bp    DNA     linear\n\
             DEFINITION  Record 2.\n\
             ACCESSION   R2\n\
             FEATURES             Location/Qualifiers\n\
             \x20    gene            1..20\n\
             \x20                    /gene=\"B\"\n\
             ORIGIN\n\
             \x20       1 atgcatgcat gcatgcatgc\n\
             //\n",
        );

        let stats = genbank_stats(file.path()).unwrap();
        assert_eq!(stats.record_count, 2);
        assert_eq!(stats.total_bases, 30);
        assert_eq!(*stats.feature_counts.get("gene").unwrap(), 2);
        assert_eq!(*stats.feature_counts.get("CDS").unwrap(), 1);
    }
}
