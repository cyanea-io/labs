//! GFF3 (General Feature Format) parser.
//!
//! Parses GFF3 files into a hierarchical [`Gene`] → [`Transcript`] → [`Exon`]
//! structure. GFF3 uses 1-based, closed coordinates; this parser converts them
//! to 0-based, half-open `[start, end)` to match the internal coordinate system.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::annotation::{Exon, Gene, GeneType, Transcript};
use cyanea_omics::genomic::Strand;

/// Parse a GFF3 file and return assembled gene models.
///
/// Features are grouped by `Parent` attribute into the hierarchy:
/// `gene` → `mRNA`/`transcript` → `exon`/`CDS`.
///
/// Coordinates are converted from GFF3's 1-based closed intervals to
/// 0-based half-open `[start, end)`.
pub fn parse_gff3(path: impl AsRef<Path>) -> Result<Vec<Gene>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    let mut builder = Gff3Builder::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        // Stop at FASTA section (must check before generic comment skip).
        if line == "##FASTA" {
            break;
        }

        // Skip comments and pragmas.
        if line.starts_with('#') {
            continue;
        }

        let record = parse_gff3_line(line, line_num + 1, path)?;
        builder.add_record(record);
    }

    Ok(builder.build())
}

/// Summary statistics for a GFF3 file.
#[derive(Debug, Clone)]
pub struct GffStats {
    pub gene_count: u64,
    pub transcript_count: u64,
    pub exon_count: u64,
    pub protein_coding_count: u64,
    pub chromosomes: Vec<String>,
}

/// Parse a GFF3 file and return summary statistics.
pub fn gff3_stats(path: impl AsRef<Path>) -> Result<GffStats> {
    let genes = parse_gff3(path)?;
    let mut chroms = Vec::new();
    let mut seen_chroms = std::collections::HashSet::new();
    let mut transcript_count: u64 = 0;
    let mut exon_count: u64 = 0;
    let mut protein_coding_count: u64 = 0;

    for gene in &genes {
        if seen_chroms.insert(gene.chrom.clone()) {
            chroms.push(gene.chrom.clone());
        }
        if gene.is_protein_coding() {
            protein_coding_count += 1;
        }
        transcript_count += gene.n_transcripts() as u64;
        for tx in &gene.transcripts {
            exon_count += tx.n_exons() as u64;
        }
    }

    Ok(GffStats {
        gene_count: genes.len() as u64,
        transcript_count,
        exon_count,
        protein_coding_count,
        chromosomes: chroms,
    })
}

// ---------------------------------------------------------------------------
// Internal types and helpers
// ---------------------------------------------------------------------------

/// A raw parsed GFF3 line before hierarchical assembly.
struct Gff3Record {
    seqid: String,
    feature_type: String,
    start: u64, // already 0-based
    end: u64,   // already half-open
    strand: Strand,
    attributes: HashMap<String, String>,
}

/// Assembles GFF3 records into hierarchical Gene models.
struct Gff3Builder {
    /// Genes keyed by their GFF3 ID attribute.
    genes: HashMap<String, GeneBuilder>,
    /// Transcripts keyed by their GFF3 ID, mapped to parent gene ID.
    transcripts: HashMap<String, (String, TranscriptBuilder)>,
    /// Exons mapped to parent transcript ID.
    exons: Vec<(String, ExonRecord)>,
    /// CDS records mapped to parent transcript ID.
    cds_records: Vec<(String, CdsRecord)>,
    /// Insertion order for genes.
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

impl Gff3Builder {
    fn new() -> Self {
        Self {
            genes: HashMap::new(),
            transcripts: HashMap::new(),
            exons: Vec::new(),
            cds_records: Vec::new(),
            gene_order: Vec::new(),
        }
    }

    fn add_record(&mut self, record: Gff3Record) {
        let ft = record.feature_type.as_str();
        match ft {
            "gene" | "pseudogene" => self.add_gene(record),
            "mRNA" | "transcript" | "lnc_RNA" | "miRNA" | "rRNA" | "tRNA" | "ncRNA"
            | "snRNA" | "snoRNA" => self.add_transcript(record),
            "exon" => self.add_exon(record),
            "CDS" => self.add_cds(record),
            _ => {} // skip other feature types
        }
    }

    fn add_gene(&mut self, record: Gff3Record) {
        let id = record
            .attributes
            .get("ID")
            .cloned()
            .unwrap_or_default();
        if id.is_empty() {
            return;
        }

        let gene_name = record
            .attributes
            .get("Name")
            .cloned()
            .unwrap_or_else(|| id.clone());

        let gene_type = record
            .attributes
            .get("biotype")
            .or_else(|| record.attributes.get("gene_biotype"))
            .map(|bt| parse_gene_type(bt))
            .unwrap_or_else(|| {
                if record.feature_type == "pseudogene" {
                    GeneType::Pseudogene
                } else {
                    GeneType::Other("unknown".to_string())
                }
            });

        self.gene_order.push(id.clone());
        self.genes.insert(
            id.clone(),
            GeneBuilder {
                gene_id: id,
                gene_name,
                chrom: record.seqid,
                start: record.start,
                end: record.end,
                strand: record.strand,
                gene_type,
            },
        );
    }

    fn add_transcript(&mut self, record: Gff3Record) {
        let id = record
            .attributes
            .get("ID")
            .cloned()
            .unwrap_or_default();
        let parent = record
            .attributes
            .get("Parent")
            .cloned()
            .unwrap_or_default();
        if id.is_empty() || parent.is_empty() {
            return;
        }

        self.transcripts.insert(
            id.clone(),
            (
                parent,
                TranscriptBuilder {
                    transcript_id: id,
                    start: record.start,
                    end: record.end,
                    exons: Vec::new(),
                    cds_start: None,
                    cds_end: None,
                },
            ),
        );
    }

    fn add_exon(&mut self, record: Gff3Record) {
        let parent = record
            .attributes
            .get("Parent")
            .cloned()
            .unwrap_or_default();
        if parent.is_empty() {
            return;
        }

        let number = record
            .attributes
            .get("exon_number")
            .or_else(|| record.attributes.get("rank"))
            .and_then(|s| s.parse().ok());

        self.exons.push((
            parent,
            ExonRecord {
                start: record.start,
                end: record.end,
                number,
            },
        ));
    }

    fn add_cds(&mut self, record: Gff3Record) {
        let parent = record
            .attributes
            .get("Parent")
            .cloned()
            .unwrap_or_default();
        if parent.is_empty() {
            return;
        }

        self.cds_records.push((
            parent,
            CdsRecord {
                start: record.start,
                end: record.end,
            },
        ));
    }

    fn build(mut self) -> Vec<Gene> {
        // Attach exons to transcripts.
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

        // Attach CDS boundaries to transcripts.
        for (parent_id, cds) in &self.cds_records {
            if let Some((_, tx)) = self.transcripts.get_mut(parent_id) {
                let cds_start = tx.cds_start.map_or(cds.start, |s| s.min(cds.start));
                let cds_end = tx.cds_end.map_or(cds.end, |e| e.max(cds.end));
                tx.cds_start = Some(cds_start);
                tx.cds_end = Some(cds_end);
            }
        }

        // Build transcripts and attach to genes.
        let mut gene_transcripts: HashMap<String, Vec<Transcript>> = HashMap::new();
        for (_, (parent_gene_id, mut tx)) in self.transcripts {
            // Sort exons by start position.
            tx.exons.sort_by_key(|e| e.start);
            // Re-number if needed.
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

        // Assemble genes in insertion order.
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

/// Parse a single GFF3 data line.
fn parse_gff3_line(line: &str, line_num: usize, path: &Path) -> Result<Gff3Record> {
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

    // GFF3 is 1-based, closed. Convert to 0-based, half-open.
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
    })?; // closed → half-open (no adjustment needed since 1-based closed [s,e] = 0-based half-open [s-1, e))

    let strand = match fields[6] {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    };

    let attributes = parse_gff3_attributes(fields[8]);

    Ok(Gff3Record {
        seqid,
        feature_type,
        start,
        end,
        strand,
        attributes,
    })
}

/// Parse GFF3 attribute column (key=value pairs separated by `;`).
fn parse_gff3_attributes(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        if let Some((key, value)) = pair.split_once('=') {
            // URL-decode percent-encoded characters.
            let value = url_decode(value);
            attrs.insert(key.to_string(), value);
        }
    }
    attrs
}

/// Minimal URL decoding for common GFF3 percent-encoded characters.
fn url_decode(s: &str) -> String {
    s.replace("%3B", ";")
        .replace("%3D", "=")
        .replace("%26", "&")
        .replace("%2C", ",")
        .replace("%25", "%")
        .replace("%09", "\t")
        .replace("%0A", "\n")
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

    fn write_gff3(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".gff3").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_parse_simple_gene() {
        let file = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=TP53;biotype=protein_coding\n\
             chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=tx1;Parent=gene1\n\
             chr1\t.\texon\t1000\t1200\t.\t+\t.\tParent=tx1;exon_number=1\n\
             chr1\t.\texon\t2000\t2500\t.\t+\t.\tParent=tx1;exon_number=2\n\
             chr1\t.\texon\t4000\t5000\t.\t+\t.\tParent=tx1;exon_number=3\n\
             chr1\t.\tCDS\t1050\t1200\t.\t+\t0\tParent=tx1\n\
             chr1\t.\tCDS\t2000\t2400\t.\t+\t0\tParent=tx1\n",
        );

        let genes = parse_gff3(file.path()).unwrap();
        assert_eq!(genes.len(), 1);

        let gene = &genes[0];
        assert_eq!(gene.gene_id, "gene1");
        assert_eq!(gene.gene_name, "TP53");
        assert_eq!(gene.chrom, "chr1");
        // 1-based 1000 → 0-based 999
        assert_eq!(gene.start, 999);
        // 1-based closed 5000 → 0-based half-open 5000
        assert_eq!(gene.end, 5000);
        assert_eq!(gene.strand, Strand::Forward);
        assert_eq!(gene.gene_type, GeneType::ProteinCoding);
        assert!(gene.is_protein_coding());

        assert_eq!(gene.n_transcripts(), 1);
        let tx = &gene.transcripts[0];
        assert_eq!(tx.transcript_id, "tx1");
        assert_eq!(tx.n_exons(), 3);

        // Exon coordinates converted to 0-based half-open.
        assert_eq!(tx.exons[0].start, 999);
        assert_eq!(tx.exons[0].end, 1200);
        assert_eq!(tx.exons[1].start, 1999);
        assert_eq!(tx.exons[1].end, 2500);

        // CDS boundaries.
        assert_eq!(tx.cds_start, Some(1049));
        assert_eq!(tx.cds_end, Some(2400));
    }

    #[test]
    fn test_multiple_transcripts() {
        let file = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\t100\t1000\t.\t-\t.\tID=g1;Name=GeneA;biotype=protein_coding\n\
             chr1\t.\tmRNA\t100\t1000\t.\t-\t.\tID=t1;Parent=g1\n\
             chr1\t.\texon\t100\t300\t.\t-\t.\tParent=t1\n\
             chr1\t.\tmRNA\t100\t800\t.\t-\t.\tID=t2;Parent=g1\n\
             chr1\t.\texon\t100\t200\t.\t-\t.\tParent=t2\n",
        );

        let genes = parse_gff3(file.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].n_transcripts(), 2);
        assert_eq!(genes[0].strand, Strand::Reverse);
    }

    #[test]
    fn test_gff3_stats() {
        let file = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;Name=A;biotype=protein_coding\n\
             chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1\n\
             chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t1\n\
             chr1\t.\texon\t600\t1000\t.\t+\t.\tParent=t1\n\
             chr2\t.\tgene\t1\t500\t.\t-\t.\tID=g2;Name=B;biotype=lncRNA\n\
             chr2\t.\tmRNA\t1\t500\t.\t-\t.\tID=t2;Parent=g2\n\
             chr2\t.\texon\t1\t500\t.\t-\t.\tParent=t2\n",
        );

        let stats = gff3_stats(file.path()).unwrap();
        assert_eq!(stats.gene_count, 2);
        assert_eq!(stats.transcript_count, 2);
        assert_eq!(stats.exon_count, 3);
        assert_eq!(stats.protein_coding_count, 1);
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn test_gff3_stops_at_fasta() {
        let file = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Name=A;biotype=protein_coding\n\
             ##FASTA\n\
             >chr1\n\
             ATCGATCG\n",
        );

        let genes = parse_gff3(file.path()).unwrap();
        assert_eq!(genes.len(), 1);
    }

    #[test]
    fn test_gff3_empty_file() {
        let file = write_gff3("##gff-version 3\n# no features\n");
        let genes = parse_gff3(file.path()).unwrap();
        assert!(genes.is_empty());
    }

    #[test]
    fn test_gff3_url_decode() {
        assert_eq!(url_decode("hello%3Bworld"), "hello;world");
        assert_eq!(url_decode("key%3Dvalue"), "key=value");
    }

    #[test]
    fn test_gff3_wrong_columns() {
        let file = write_gff3("chr1\tgene\n");
        let result = parse_gff3(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_gff3_file_not_found() {
        let result = parse_gff3("/nonexistent/file.gff3");
        assert!(result.is_err());
    }

    #[test]
    fn test_pseudogene_type() {
        let file = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tpseudogene\t1\t100\t.\t+\t.\tID=pg1;Name=PseudoA\n",
        );
        let genes = parse_gff3(file.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].gene_type, GeneType::Pseudogene);
    }
}
