//! Structured VCF header construction and parsing.
//!
//! Provides [`VcfHeader`] for building, serializing, and parsing VCF 4.3 headers
//! with contig, INFO, FORMAT, and FILTER field definitions.

use cyanea_core::{CyaneaError, Result};

/// A contig (reference sequence) definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigLine {
    /// Contig identifier (e.g. "chr1").
    pub id: String,
    /// Optional contig length.
    pub length: Option<u64>,
}

/// An INFO or FORMAT field definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldDef {
    /// Field identifier.
    pub id: String,
    /// Number of values (e.g. "1", "A", "R", "G", ".").
    pub number: String,
    /// Value type (e.g. "Integer", "Float", "String", "Flag", "Character").
    pub field_type: String,
    /// Human-readable description.
    pub description: String,
}

/// A FILTER field definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterDef {
    /// Filter identifier.
    pub id: String,
    /// Human-readable description.
    pub description: String,
}

/// A structured VCF header.
///
/// Build headers programmatically with [`VcfHeader::new`] and the `add_*` methods,
/// or parse from text with [`VcfHeader::parse`].
#[derive(Debug, Clone, Default)]
pub struct VcfHeader {
    /// Contig (reference sequence) definitions.
    pub contigs: Vec<ContigLine>,
    /// INFO field definitions.
    pub info_fields: Vec<FieldDef>,
    /// FORMAT field definitions.
    pub format_fields: Vec<FieldDef>,
    /// FILTER definitions.
    pub filter_fields: Vec<FilterDef>,
    /// Sample names (columns after FORMAT).
    pub samples: Vec<String>,
    /// Extra header lines not covered above (stored verbatim).
    pub extra_lines: Vec<String>,
}

impl VcfHeader {
    /// Create an empty VCF header.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a contig definition.
    pub fn add_contig(&mut self, id: &str, length: Option<u64>) {
        self.contigs.push(ContigLine {
            id: id.to_string(),
            length,
        });
    }

    /// Add an INFO field definition.
    pub fn add_info(&mut self, id: &str, number: &str, field_type: &str, desc: &str) {
        self.info_fields.push(FieldDef {
            id: id.to_string(),
            number: number.to_string(),
            field_type: field_type.to_string(),
            description: desc.to_string(),
        });
    }

    /// Add a FORMAT field definition.
    pub fn add_format(&mut self, id: &str, number: &str, field_type: &str, desc: &str) {
        self.format_fields.push(FieldDef {
            id: id.to_string(),
            number: number.to_string(),
            field_type: field_type.to_string(),
            description: desc.to_string(),
        });
    }

    /// Add a FILTER definition.
    pub fn add_filter(&mut self, id: &str, desc: &str) {
        self.filter_fields.push(FilterDef {
            id: id.to_string(),
            description: desc.to_string(),
        });
    }

    /// Add a sample name.
    pub fn add_sample(&mut self, name: &str) {
        self.samples.push(name.to_string());
    }

    /// Serialize the header to VCF text (including the `#CHROM` line).
    pub fn to_vcf_string(&self) -> String {
        let mut out = String::new();
        out.push_str("##fileformat=VCFv4.3\n");

        for c in &self.contigs {
            if let Some(len) = c.length {
                out.push_str(&format!("##contig=<ID={},length={}>\n", c.id, len));
            } else {
                out.push_str(&format!("##contig=<ID={}>\n", c.id));
            }
        }

        for f in &self.info_fields {
            out.push_str(&format!(
                "##INFO=<ID={},Number={},Type={},Description=\"{}\">\n",
                f.id, f.number, f.field_type, f.description
            ));
        }

        for f in &self.format_fields {
            out.push_str(&format!(
                "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">\n",
                f.id, f.number, f.field_type, f.description
            ));
        }

        for f in &self.filter_fields {
            out.push_str(&format!(
                "##FILTER=<ID={},Description=\"{}\">\n",
                f.id, f.description
            ));
        }

        for line in &self.extra_lines {
            out.push_str(line);
            if !line.ends_with('\n') {
                out.push('\n');
            }
        }

        out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            out.push_str("\tFORMAT");
            for s in &self.samples {
                out.push('\t');
                out.push_str(s);
            }
        }
        out.push('\n');

        out
    }

    /// Parse a VCF header from text.
    ///
    /// Expects lines starting with `##` (meta-information) and one `#CHROM` header line.
    /// Only header lines are parsed; data lines are ignored.
    pub fn parse(text: &str) -> Result<Self> {
        let mut header = Self::new();

        for line in text.lines() {
            let line = line.trim();
            if line.starts_with("##contig=<") {
                let inner = extract_angle_bracket_content(line, "##contig=")?;
                let fields = parse_meta_fields(&inner);
                let id = fields
                    .get("ID")
                    .ok_or_else(|| CyaneaError::Parse("contig missing ID".into()))?
                    .to_string();
                let length = fields.get("length").and_then(|v| v.parse().ok());
                header.contigs.push(ContigLine { id, length });
            } else if line.starts_with("##INFO=<") {
                let def = parse_field_def(line, "##INFO=")?;
                header.info_fields.push(def);
            } else if line.starts_with("##FORMAT=<") {
                let def = parse_field_def(line, "##FORMAT=")?;
                header.format_fields.push(def);
            } else if line.starts_with("##FILTER=<") {
                let inner = extract_angle_bracket_content(line, "##FILTER=")?;
                let fields = parse_meta_fields(&inner);
                let id = fields
                    .get("ID")
                    .ok_or_else(|| CyaneaError::Parse("FILTER missing ID".into()))?
                    .to_string();
                let description = fields
                    .get("Description")
                    .map(|s| s.trim_matches('"').to_string())
                    .unwrap_or_default();
                header.filter_fields.push(FilterDef { id, description });
            } else if line.starts_with("#CHROM") {
                // Parse sample names from the #CHROM line
                let cols: Vec<&str> = line.split('\t').collect();
                // Standard columns: CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT sample1 ...]
                if cols.len() > 9 {
                    for sample in &cols[9..] {
                        header.samples.push(sample.to_string());
                    }
                }
            } else if line.starts_with("##fileformat=") {
                // Skip fileformat line (we always write VCFv4.3)
            } else if line.starts_with("##") {
                header.extra_lines.push(line.to_string());
            }
            // Skip data lines (don't start with #)
        }

        Ok(header)
    }
}

/// Extract content between `<` and `>` after a prefix.
fn extract_angle_bracket_content(line: &str, prefix: &str) -> Result<String> {
    let rest = line
        .strip_prefix(prefix)
        .ok_or_else(|| CyaneaError::Parse(format!("expected prefix '{}'", prefix)))?;
    let inner = rest
        .strip_prefix('<')
        .and_then(|s| s.strip_suffix('>'))
        .ok_or_else(|| CyaneaError::Parse("missing angle brackets".into()))?;
    Ok(inner.to_string())
}

/// Parse comma-separated key=value fields from a meta-information line.
///
/// Handles quoted values containing commas (e.g. `Description="some,text"`).
fn parse_meta_fields(inner: &str) -> std::collections::HashMap<String, String> {
    let mut fields = std::collections::HashMap::new();
    let bytes = inner.as_bytes();
    let mut pos = 0;
    let len = bytes.len();

    while pos < len {
        // Find '='
        let eq_pos = match bytes[pos..].iter().position(|&b| b == b'=') {
            Some(p) => pos + p,
            None => break,
        };
        let key = inner[pos..eq_pos].to_string();
        pos = eq_pos + 1;

        if pos < len && bytes[pos] == b'"' {
            // Quoted value — find closing quote
            pos += 1; // skip opening quote
            let close = match bytes[pos..].iter().position(|&b| b == b'"') {
                Some(p) => pos + p,
                None => len,
            };
            let value = inner[pos..close].to_string();
            fields.insert(key, value);
            pos = close + 1; // skip closing quote
            if pos < len && bytes[pos] == b',' {
                pos += 1; // skip comma
            }
        } else {
            // Unquoted value — find next comma
            let comma = bytes[pos..]
                .iter()
                .position(|&b| b == b',')
                .map(|p| pos + p)
                .unwrap_or(len);
            let value = inner[pos..comma].to_string();
            fields.insert(key, value);
            pos = comma + 1;
        }
    }

    fields
}

/// Parse an INFO or FORMAT field definition.
fn parse_field_def(line: &str, prefix: &str) -> Result<FieldDef> {
    let inner = extract_angle_bracket_content(line, prefix)?;
    let fields = parse_meta_fields(&inner);

    let id = fields
        .get("ID")
        .ok_or_else(|| CyaneaError::Parse(format!("{} missing ID", prefix.trim_end_matches('='))))?
        .to_string();
    let number = fields.get("Number").cloned().unwrap_or_else(|| ".".to_string());
    let field_type = fields.get("Type").cloned().unwrap_or_else(|| "String".to_string());
    let description = fields
        .get("Description")
        .map(|s| s.trim_matches('"').to_string())
        .unwrap_or_default();

    Ok(FieldDef {
        id,
        number,
        field_type,
        description,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_empty_header() {
        let h = VcfHeader::new();
        let s = h.to_vcf_string();
        assert!(s.starts_with("##fileformat=VCFv4.3\n"));
        assert!(s.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    }

    #[test]
    fn build_header_with_contigs() {
        let mut h = VcfHeader::new();
        h.add_contig("chr1", Some(248956422));
        h.add_contig("chrM", None);
        let s = h.to_vcf_string();
        assert!(s.contains("##contig=<ID=chr1,length=248956422>"));
        assert!(s.contains("##contig=<ID=chrM>"));
    }

    #[test]
    fn build_header_with_info() {
        let mut h = VcfHeader::new();
        h.add_info("DP", "1", "Integer", "Total Depth");
        h.add_info("AF", "A", "Float", "Allele Frequency");
        let s = h.to_vcf_string();
        assert!(s.contains("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"));
        assert!(s.contains("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
    }

    #[test]
    fn build_header_with_format() {
        let mut h = VcfHeader::new();
        h.add_format("GT", "1", "String", "Genotype");
        let s = h.to_vcf_string();
        assert!(s.contains("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    }

    #[test]
    fn build_header_with_filters() {
        let mut h = VcfHeader::new();
        h.add_filter("LowQual", "Low quality");
        h.add_filter("StrandBias", "Strand bias detected");
        let s = h.to_vcf_string();
        assert!(s.contains("##FILTER=<ID=LowQual,Description=\"Low quality\">"));
        assert!(s.contains("##FILTER=<ID=StrandBias,Description=\"Strand bias detected\">"));
    }

    #[test]
    fn build_header_with_samples() {
        let mut h = VcfHeader::new();
        h.add_format("GT", "1", "String", "Genotype");
        h.add_sample("SAMPLE1");
        h.add_sample("SAMPLE2");
        let s = h.to_vcf_string();
        assert!(s.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\n"));
    }

    #[test]
    fn build_header_no_samples_no_format() {
        let h = VcfHeader::new();
        let s = h.to_vcf_string();
        assert!(!s.contains("FORMAT"));
    }

    #[test]
    fn parse_empty_header() {
        let text = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let h = VcfHeader::parse(text).unwrap();
        assert!(h.contigs.is_empty());
        assert!(h.info_fields.is_empty());
        assert!(h.samples.is_empty());
    }

    #[test]
    fn parse_header_contigs() {
        let text = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##contig=<ID=chrM>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";
        let h = VcfHeader::parse(text).unwrap();
        assert_eq!(h.contigs.len(), 2);
        assert_eq!(h.contigs[0].id, "chr1");
        assert_eq!(h.contigs[0].length, Some(248956422));
        assert_eq!(h.contigs[1].id, "chrM");
        assert_eq!(h.contigs[1].length, None);
    }

    #[test]
    fn parse_header_info_format_filter() {
        let text = "\
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FILTER=<ID=LowQual,Description=\"Low quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
";
        let h = VcfHeader::parse(text).unwrap();
        assert_eq!(h.info_fields.len(), 1);
        assert_eq!(h.info_fields[0].id, "DP");
        assert_eq!(h.info_fields[0].number, "1");
        assert_eq!(h.info_fields[0].field_type, "Integer");
        assert_eq!(h.info_fields[0].description, "Total Depth");

        assert_eq!(h.format_fields.len(), 1);
        assert_eq!(h.format_fields[0].id, "GT");

        assert_eq!(h.filter_fields.len(), 1);
        assert_eq!(h.filter_fields[0].id, "LowQual");

        assert_eq!(h.samples, vec!["SAMPLE1"]);
    }

    #[test]
    fn parse_roundtrip() {
        let mut h = VcfHeader::new();
        h.add_contig("chr1", Some(1000));
        h.add_contig("chr2", Some(2000));
        h.add_info("DP", "1", "Integer", "Total Depth");
        h.add_format("GT", "1", "String", "Genotype");
        h.add_filter("LowQual", "Low quality");
        h.add_sample("S1");

        let text = h.to_vcf_string();
        let parsed = VcfHeader::parse(&text).unwrap();

        assert_eq!(parsed.contigs.len(), 2);
        assert_eq!(parsed.contigs[0].id, "chr1");
        assert_eq!(parsed.contigs[0].length, Some(1000));
        assert_eq!(parsed.info_fields.len(), 1);
        assert_eq!(parsed.info_fields[0].id, "DP");
        assert_eq!(parsed.format_fields.len(), 1);
        assert_eq!(parsed.filter_fields.len(), 1);
        assert_eq!(parsed.samples, vec!["S1"]);
    }

    #[test]
    fn parse_extra_lines_preserved() {
        let text = "\
##fileformat=VCFv4.3
##source=Cyanea
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";
        let h = VcfHeader::parse(text).unwrap();
        assert_eq!(h.extra_lines.len(), 2);
        assert!(h.extra_lines[0].contains("source=Cyanea"));
        assert!(h.extra_lines[1].contains("reference=GRCh38"));
    }

    #[test]
    fn parse_multiple_samples() {
        let text = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891\tNA12892
";
        let h = VcfHeader::parse(text).unwrap();
        assert_eq!(h.samples, vec!["NA12878", "NA12891", "NA12892"]);
    }

    #[test]
    fn parse_description_with_comma() {
        let text = "\
##fileformat=VCFv4.3
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";
        let h = VcfHeader::parse(text).unwrap();
        assert_eq!(h.info_fields[0].description, "Allele Frequency, for each ALT allele");
    }

    #[test]
    fn contig_line_without_length() {
        let mut h = VcfHeader::new();
        h.add_contig("chrUn", None);
        let s = h.to_vcf_string();
        assert!(s.contains("##contig=<ID=chrUn>\n"));
        assert!(!s.contains("length"));
    }
}
