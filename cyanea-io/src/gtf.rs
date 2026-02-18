//! GTF (Gene Transfer Format / GFF2) parser.
//!
//! GTF uses the same 9-column tab format as GFF3 but with different attribute
//! syntax. GFF3 uses `key=value;` pairs, while GTF uses space-separated
//! `key "value";` pairs with quoted values.
//!
//! Parses GTF files into the same hierarchical [`Gene`] → [`Transcript`] →
//! [`Exon`] structure as the GFF3 parser. Coordinates are converted from
//! GTF's 1-based closed intervals to 0-based half-open `[start, end)`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::annotation::{Exon, Gene, GeneType, Transcript};
use cyanea_omics::genomic::Strand;

/// Parse a GTF file and return assembled gene models.
///
/// Features are grouped by `gene_id` and `transcript_id` attributes into
/// the hierarchy: `gene` → `transcript` → `exon`/`CDS`.
///
/// Coordinates are converted from GTF's 1-based closed intervals to
/// 0-based half-open `[start, end)`.
pub fn parse_gtf(path: impl AsRef<Path>) -> Result<Vec<Gene>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    let mut builder = GtfBuilder::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let record = parse_gtf_line(line, line_num + 1, path)?;
        builder.add_record(record);
    }

    Ok(builder.build())
}

/// Summary statistics for a GTF file.
#[derive(Debug, Clone)]
pub struct GtfStats {
    /// Number of genes.
    pub gene_count: u64,
    /// Number of transcripts.
    pub transcript_count: u64,
    /// Number of exons.
    pub exon_count: u64,
    /// Number of CDS features.
    pub cds_count: u64,
}

/// Parse a GTF file and return summary statistics.
pub fn gtf_stats(path: impl AsRef<Path>) -> Result<GtfStats> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut gene_ids = std::collections::HashSet::new();
    let mut transcript_ids = std::collections::HashSet::new();
    let mut exon_count: u64 = 0;
    let mut cds_count: u64 = 0;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature_type = fields[2];
        let attrs = parse_gtf_attributes(fields[8]);

        match feature_type {
            "gene" => {
                if let Some(gid) = attrs.get("gene_id") {
                    gene_ids.insert(gid.clone());
                }
            }
            "transcript" => {
                if let Some(tid) = attrs.get("transcript_id") {
                    transcript_ids.insert(tid.clone());
                }
            }
            "exon" => {
                exon_count += 1;
                // Also capture gene/transcript IDs from exon lines
                // in case there are no explicit gene/transcript features
                if let Some(gid) = attrs.get("gene_id") {
                    gene_ids.insert(gid.clone());
                }
                if let Some(tid) = attrs.get("transcript_id") {
                    transcript_ids.insert(tid.clone());
                }
            }
            "CDS" => {
                cds_count += 1;
                if let Some(gid) = attrs.get("gene_id") {
                    gene_ids.insert(gid.clone());
                }
                if let Some(tid) = attrs.get("transcript_id") {
                    transcript_ids.insert(tid.clone());
                }
            }
            _ => {}
        }
    }

    Ok(GtfStats {
        gene_count: gene_ids.len() as u64,
        transcript_count: transcript_ids.len() as u64,
        exon_count,
        cds_count,
    })
}

// ---------------------------------------------------------------------------
// Internal types and helpers
// ---------------------------------------------------------------------------

struct GtfRecord {
    seqid: String,
    feature_type: String,
    start: u64, // already 0-based
    end: u64,   // already half-open
    strand: Strand,
    attributes: HashMap<String, String>,
}

struct GtfBuilder {
    genes: HashMap<String, GeneBuilder>,
    transcripts: HashMap<String, (String, TranscriptBuilder)>,
    exons: Vec<(String, ExonRecord)>,
    cds_records: Vec<(String, CdsRecord)>,
    gene_order: Vec<String>,
}

struct GeneBuilder {
    gene_id: String,
    gene_name: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: Strand,
    gene_type: GeneType,
}

struct TranscriptBuilder {
    transcript_id: String,
    start: u64,
    end: u64,
    exons: Vec<Exon>,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
}

struct ExonRecord {
    start: u64,
    end: u64,
    number: Option<u32>,
}

struct CdsRecord {
    start: u64,
    end: u64,
}

impl GtfBuilder {
    fn new() -> Self {
        Self {
            genes: HashMap::new(),
            transcripts: HashMap::new(),
            exons: Vec::new(),
            cds_records: Vec::new(),
            gene_order: Vec::new(),
        }
    }

    fn add_record(&mut self, record: GtfRecord) {
        let ft = record.feature_type.as_str();
        match ft {
            "gene" => self.add_gene(record),
            "transcript" | "mRNA" => self.add_transcript(record),
            "exon" => self.add_exon(record),
            "CDS" => self.add_cds(record),
            _ => {}
        }
    }

    fn add_gene(&mut self, record: GtfRecord) {
        let gene_id = record
            .attributes
            .get("gene_id")
            .cloned()
            .unwrap_or_default();
        if gene_id.is_empty() {
            return;
        }

        let gene_name = record
            .attributes
            .get("gene_name")
            .cloned()
            .unwrap_or_else(|| gene_id.clone());

        let gene_type = record
            .attributes
            .get("gene_biotype")
            .or_else(|| record.attributes.get("gene_type"))
            .map(|bt| parse_gene_type(bt))
            .unwrap_or(GeneType::Other("unknown".to_string()));

        if !self.genes.contains_key(&gene_id) {
            self.gene_order.push(gene_id.clone());
        }
        self.genes.insert(
            gene_id.clone(),
            GeneBuilder {
                gene_id,
                gene_name,
                chrom: record.seqid,
                start: record.start,
                end: record.end,
                strand: record.strand,
                gene_type,
            },
        );
    }

    fn add_transcript(&mut self, record: GtfRecord) {
        let transcript_id = record
            .attributes
            .get("transcript_id")
            .cloned()
            .unwrap_or_default();
        let gene_id = record
            .attributes
            .get("gene_id")
            .cloned()
            .unwrap_or_default();
        if transcript_id.is_empty() || gene_id.is_empty() {
            return;
        }

        // Auto-create gene if not yet seen (GTF may lack explicit gene lines)
        if !self.genes.contains_key(&gene_id) {
            let gene_name = record
                .attributes
                .get("gene_name")
                .cloned()
                .unwrap_or_else(|| gene_id.clone());
            let gene_type = record
                .attributes
                .get("gene_biotype")
                .or_else(|| record.attributes.get("gene_type"))
                .map(|bt| parse_gene_type(bt))
                .unwrap_or(GeneType::Other("unknown".to_string()));

            self.gene_order.push(gene_id.clone());
            self.genes.insert(
                gene_id.clone(),
                GeneBuilder {
                    gene_id: gene_id.clone(),
                    gene_name,
                    chrom: record.seqid.clone(),
                    start: record.start,
                    end: record.end,
                    strand: record.strand,
                    gene_type,
                },
            );
        }

        self.transcripts.insert(
            transcript_id.clone(),
            (
                gene_id,
                TranscriptBuilder {
                    transcript_id,
                    start: record.start,
                    end: record.end,
                    exons: Vec::new(),
                    cds_start: None,
                    cds_end: None,
                },
            ),
        );
    }

    fn add_exon(&mut self, record: GtfRecord) {
        let transcript_id = record
            .attributes
            .get("transcript_id")
            .cloned()
            .unwrap_or_default();
        let gene_id = record
            .attributes
            .get("gene_id")
            .cloned()
            .unwrap_or_default();
        if transcript_id.is_empty() {
            return;
        }

        // Auto-create gene and transcript if not yet seen
        if !gene_id.is_empty() && !self.genes.contains_key(&gene_id) {
            let gene_name = record
                .attributes
                .get("gene_name")
                .cloned()
                .unwrap_or_else(|| gene_id.clone());
            self.gene_order.push(gene_id.clone());
            self.genes.insert(
                gene_id.clone(),
                GeneBuilder {
                    gene_id: gene_id.clone(),
                    gene_name,
                    chrom: record.seqid.clone(),
                    start: record.start,
                    end: record.end,
                    strand: record.strand,
                    gene_type: GeneType::Other("unknown".to_string()),
                },
            );
        }
        if !self.transcripts.contains_key(&transcript_id) && !gene_id.is_empty() {
            self.transcripts.insert(
                transcript_id.clone(),
                (
                    gene_id,
                    TranscriptBuilder {
                        transcript_id: transcript_id.clone(),
                        start: record.start,
                        end: record.end,
                        exons: Vec::new(),
                        cds_start: None,
                        cds_end: None,
                    },
                ),
            );
        }

        let number = record
            .attributes
            .get("exon_number")
            .and_then(|s| s.parse().ok());

        self.exons.push((
            transcript_id,
            ExonRecord {
                start: record.start,
                end: record.end,
                number,
            },
        ));
    }

    fn add_cds(&mut self, record: GtfRecord) {
        let transcript_id = record
            .attributes
            .get("transcript_id")
            .cloned()
            .unwrap_or_default();
        if transcript_id.is_empty() {
            return;
        }

        self.cds_records.push((
            transcript_id,
            CdsRecord {
                start: record.start,
                end: record.end,
            },
        ));
    }

    fn build(mut self) -> Vec<Gene> {
        // Attach exons to transcripts
        for (parent_id, exon) in self.exons {
            if let Some((_, tx)) = self.transcripts.get_mut(&parent_id) {
                let exon_number = exon.number.unwrap_or((tx.exons.len() + 1) as u32);
                tx.exons.push(Exon {
                    exon_number,
                    start: exon.start,
                    end: exon.end,
                });
            }
        }

        // Attach CDS boundaries to transcripts
        for (parent_id, cds) in &self.cds_records {
            if let Some((_, tx)) = self.transcripts.get_mut(parent_id) {
                let cds_start = tx.cds_start.map_or(cds.start, |s| s.min(cds.start));
                let cds_end = tx.cds_end.map_or(cds.end, |e| e.max(cds.end));
                tx.cds_start = Some(cds_start);
                tx.cds_end = Some(cds_end);
            }
        }

        // Build transcripts and attach to genes
        let mut gene_transcripts: HashMap<String, Vec<Transcript>> = HashMap::new();
        for (_, (parent_gene_id, mut tx)) in self.transcripts {
            tx.exons.sort_by_key(|e| e.start);
            for (i, exon) in tx.exons.iter_mut().enumerate() {
                exon.exon_number = (i + 1) as u32;
            }

            let transcript = Transcript {
                transcript_id: tx.transcript_id,
                start: tx.start,
                end: tx.end,
                exons: tx.exons,
                cds_start: tx.cds_start,
                cds_end: tx.cds_end,
            };

            gene_transcripts
                .entry(parent_gene_id)
                .or_default()
                .push(transcript);
        }

        // Assemble genes in insertion order
        let mut genes = Vec::new();
        for gene_id in &self.gene_order {
            if let Some(gb) = self.genes.remove(gene_id) {
                let transcripts = gene_transcripts.remove(gene_id).unwrap_or_default();
                genes.push(Gene {
                    gene_id: gb.gene_id,
                    gene_name: gb.gene_name,
                    chrom: gb.chrom,
                    start: gb.start,
                    end: gb.end,
                    strand: gb.strand,
                    gene_type: gb.gene_type,
                    transcripts,
                });
            }
        }

        genes
    }
}

/// Parse a single GTF data line.
fn parse_gtf_line(line: &str, line_num: usize, path: &Path) -> Result<GtfRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 9 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected 9 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let seqid = fields[0].to_string();
    let feature_type = fields[2].to_string();

    // GTF is 1-based, closed. Convert to 0-based, half-open.
    let start: u64 = fields[3]
        .parse::<u64>()
        .map_err(|_| {
            CyaneaError::Parse(format!(
                "{}: line {}: invalid start '{}'",
                path.display(),
                line_num,
                fields[3]
            ))
        })?
        .saturating_sub(1); // 1-based → 0-based

    let end: u64 = fields[4].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid end '{}'",
            path.display(),
            line_num,
            fields[4]
        ))
    })?; // 1-based closed → 0-based half-open (no adjustment)

    let strand = match fields[6] {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    };

    let attributes = parse_gtf_attributes(fields[8]);

    Ok(GtfRecord {
        seqid,
        feature_type,
        start,
        end,
        strand,
        attributes,
    })
}

/// Parse GTF attribute column.
///
/// GTF format: `gene_id "ENSG00000141510"; gene_name "TP53"; gene_biotype "protein_coding";`
/// Key-value pairs separated by `;`, values quoted with `"`.
fn parse_gtf_attributes(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        // Split on first whitespace: key "value"
        if let Some(space_pos) = pair.find(|c: char| c.is_whitespace()) {
            let key = pair[..space_pos].trim();
            let value = pair[space_pos..].trim();
            // Strip surrounding quotes
            let value = value
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .unwrap_or(value);
            attrs.insert(key.to_string(), value.to_string());
        }
    }
    attrs
}

/// Map biotype strings to GeneType enum.
fn parse_gene_type(biotype: &str) -> GeneType {
    match biotype {
        "protein_coding" => GeneType::ProteinCoding,
        "lncRNA" | "lnc_RNA" => GeneType::LncRNA,
        "miRNA" => GeneType::MiRNA,
        "rRNA" => GeneType::RRNA,
        "tRNA" => GeneType::TRNA,
        "pseudogene" | "processed_pseudogene" | "unprocessed_pseudogene"
        | "transcribed_pseudogene" => GeneType::Pseudogene,
        other => GeneType::Other(other.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_gtf(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".gtf").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn gtf_parse_simple() {
        let file = write_gtf(
            "#!genome-build GRCh38\n\
             chr1\tENSEMBL\tgene\t1000\t5000\t.\t+\t.\tgene_id \"ENSG001\"; gene_name \"TP53\"; gene_biotype \"protein_coding\";\n\
             chr1\tENSEMBL\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"ENSG001\"; transcript_id \"ENST001\";\n\
             chr1\tENSEMBL\texon\t1000\t1200\t.\t+\t.\tgene_id \"ENSG001\"; transcript_id \"ENST001\"; exon_number \"1\";\n\
             chr2\tENSEMBL\tgene\t2000\t3000\t.\t-\t.\tgene_id \"ENSG002\"; gene_name \"BRCA1\"; gene_biotype \"protein_coding\";\n\
             chr2\tENSEMBL\ttranscript\t2000\t3000\t.\t-\t.\tgene_id \"ENSG002\"; transcript_id \"ENST002\";\n\
             chr2\tENSEMBL\texon\t2000\t2500\t.\t-\t.\tgene_id \"ENSG002\"; transcript_id \"ENST002\"; exon_number \"1\";\n\
             chr3\tENSEMBL\tgene\t5000\t8000\t.\t+\t.\tgene_id \"ENSG003\"; gene_name \"MYC\"; gene_biotype \"protein_coding\";\n",
        );

        let genes = parse_gtf(file.path()).unwrap();
        assert_eq!(genes.len(), 3);
        assert_eq!(genes[0].gene_name, "TP53");
        assert_eq!(genes[1].gene_name, "BRCA1");
        assert_eq!(genes[2].gene_name, "MYC");
    }

    #[test]
    fn gtf_attributes_quoted() {
        let attrs = parse_gtf_attributes(
            "gene_id \"ENSG00000141510\"; gene_name \"TP53\"; gene_biotype \"protein_coding\";",
        );
        assert_eq!(attrs.get("gene_id").unwrap(), "ENSG00000141510");
        assert_eq!(attrs.get("gene_name").unwrap(), "TP53");
        assert_eq!(attrs.get("gene_biotype").unwrap(), "protein_coding");
    }

    #[test]
    fn gtf_hierarchy() {
        let file = write_gtf(
            "chr1\tENSEMBL\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"GeneA\"; gene_biotype \"protein_coding\";\n\
             chr1\tENSEMBL\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tENSEMBL\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number \"1\";\n\
             chr1\tENSEMBL\texon\t2000\t2500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number \"2\";\n\
             chr1\tENSEMBL\texon\t4000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number \"3\";\n\
             chr1\tENSEMBL\tCDS\t1050\t1200\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tENSEMBL\tCDS\t2000\t2400\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(file.path()).unwrap();
        assert_eq!(genes.len(), 1);

        let gene = &genes[0];
        assert_eq!(gene.gene_id, "G1");
        assert_eq!(gene.gene_name, "GeneA");
        assert_eq!(gene.gene_type, GeneType::ProteinCoding);
        assert_eq!(gene.n_transcripts(), 1);

        let tx = &gene.transcripts[0];
        assert_eq!(tx.transcript_id, "T1");
        assert_eq!(tx.n_exons(), 3);
        assert_eq!(tx.cds_start, Some(1049));
        assert_eq!(tx.cds_end, Some(2400));
    }

    #[test]
    fn gtf_coordinates() {
        let file = write_gtf(
            "chr1\tENSEMBL\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"A\";\n\
             chr1\tENSEMBL\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tENSEMBL\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(file.path()).unwrap();
        // 1-based 1000 → 0-based 999
        assert_eq!(genes[0].start, 999);
        // 1-based closed 5000 → 0-based half-open 5000
        assert_eq!(genes[0].end, 5000);

        let exon = &genes[0].transcripts[0].exons[0];
        assert_eq!(exon.start, 999);
        assert_eq!(exon.end, 1200);
    }

    #[test]
    fn gtf_stats_computed() {
        let file = write_gtf(
            "chr1\tENSEMBL\tgene\t1\t1000\t.\t+\t.\tgene_id \"G1\"; gene_name \"A\";\n\
             chr1\tENSEMBL\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tENSEMBL\texon\t1\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tENSEMBL\texon\t600\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tENSEMBL\tCDS\t100\t500\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr2\tENSEMBL\tgene\t1\t500\t.\t-\t.\tgene_id \"G2\"; gene_name \"B\";\n\
             chr2\tENSEMBL\ttranscript\t1\t500\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";\n\
             chr2\tENSEMBL\texon\t1\t500\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";\n",
        );

        let stats = gtf_stats(file.path()).unwrap();
        assert_eq!(stats.gene_count, 2);
        assert_eq!(stats.transcript_count, 2);
        assert_eq!(stats.exon_count, 3);
        assert_eq!(stats.cds_count, 1);
    }

    #[test]
    fn gtf_comment_lines_skipped() {
        let file = write_gtf(
            "#!genome-build GRCh38.p14\n\
             #!genome-version GRCh38\n\
             #!genome-date 2013-12\n\
             chr1\tENSEMBL\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_name \"A\";\n",
        );

        let genes = parse_gtf(file.path()).unwrap();
        assert_eq!(genes.len(), 1);
    }
}
