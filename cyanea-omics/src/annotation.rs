//! Gene, transcript, and exon annotation types.
//!
//! Hierarchical gene → transcript → exon model for representing genome
//! annotations from sources like GENCODE, Ensembl, or RefSeq.

use cyanea_core::{Annotated, Summarizable};

use crate::genomic::{GenomicInterval, Strand};

/// Classification of a gene's biotype.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GeneType {
    ProteinCoding,
    LncRNA,
    MiRNA,
    RRNA,
    TRNA,
    Pseudogene,
    Other(String),
}

impl core::fmt::Display for GeneType {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            GeneType::ProteinCoding => write!(f, "protein_coding"),
            GeneType::LncRNA => write!(f, "lncRNA"),
            GeneType::MiRNA => write!(f, "miRNA"),
            GeneType::RRNA => write!(f, "rRNA"),
            GeneType::TRNA => write!(f, "tRNA"),
            GeneType::Pseudogene => write!(f, "pseudogene"),
            GeneType::Other(s) => write!(f, "{s}"),
        }
    }
}

/// An exon within a transcript (0-based coordinates).
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Exon {
    pub exon_number: u32,
    /// 0-based start (inclusive).
    pub start: u64,
    /// 0-based end (exclusive).
    pub end: u64,
}

impl Exon {
    /// Length of the exon in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the exon has zero length.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

/// A transcript with optional CDS boundaries.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Transcript {
    pub transcript_id: String,
    /// 0-based start (inclusive).
    pub start: u64,
    /// 0-based end (exclusive).
    pub end: u64,
    pub exons: Vec<Exon>,
    /// CDS start (0-based inclusive), if protein-coding.
    pub cds_start: Option<u64>,
    /// CDS end (0-based exclusive), if protein-coding.
    pub cds_end: Option<u64>,
}

impl Transcript {
    /// Length of the transcript span in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the transcript has zero length.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// Number of exons.
    pub fn n_exons(&self) -> usize {
        self.exons.len()
    }

    /// Total exonic length (sum of individual exon lengths).
    pub fn exonic_length(&self) -> u64 {
        self.exons.iter().map(|e| e.len()).sum()
    }

    /// Convert this transcript to a [`GenomicInterval`] on the given chromosome and strand.
    pub fn to_genomic_interval(&self, chrom: &str, strand: Strand) -> GenomicInterval {
        GenomicInterval {
            chrom: chrom.into(),
            start: self.start,
            end: self.end,
            strand,
        }
    }
}

/// A gene with its transcripts.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Gene {
    pub gene_id: String,
    pub gene_name: String,
    pub chrom: String,
    /// 0-based start (inclusive).
    pub start: u64,
    /// 0-based end (exclusive).
    pub end: u64,
    pub strand: Strand,
    pub gene_type: GeneType,
    pub transcripts: Vec<Transcript>,
}

impl Gene {
    /// Length of the gene span in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the gene has zero length.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// Number of transcripts.
    pub fn n_transcripts(&self) -> usize {
        self.transcripts.len()
    }

    /// Convert to a [`GenomicInterval`].
    pub fn to_genomic_interval(&self) -> GenomicInterval {
        GenomicInterval {
            chrom: self.chrom.clone(),
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }

    /// Whether this gene is protein-coding.
    pub fn is_protein_coding(&self) -> bool {
        self.gene_type == GeneType::ProteinCoding
    }
}

impl Annotated for Gene {
    fn name(&self) -> &str {
        &self.gene_name
    }
}

impl Summarizable for Gene {
    fn summary(&self) -> String {
        format!(
            "Gene: {} ({}:{}-{}, {}, {}, {} transcripts)",
            self.gene_name,
            self.chrom,
            self.start,
            self.end,
            self.strand,
            self.gene_type,
            self.n_transcripts()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_gene() -> Gene {
        Gene {
            gene_id: "ENSG00000141510".into(),
            gene_name: "TP53".into(),
            chrom: "chr17".into(),
            start: 7668421,
            end: 7687490,
            strand: Strand::Reverse,
            gene_type: GeneType::ProteinCoding,
            transcripts: vec![
                Transcript {
                    transcript_id: "ENST00000269305".into(),
                    start: 7668421,
                    end: 7687490,
                    exons: vec![
                        Exon { exon_number: 1, start: 7668421, end: 7668586 },
                        Exon { exon_number: 2, start: 7670609, end: 7670715 },
                        Exon { exon_number: 3, start: 7673534, end: 7673608 },
                    ],
                    cds_start: Some(7668421),
                    cds_end: Some(7687490),
                },
                Transcript {
                    transcript_id: "ENST00000413465".into(),
                    start: 7669608,
                    end: 7687490,
                    exons: vec![
                        Exon { exon_number: 1, start: 7669608, end: 7669690 },
                    ],
                    cds_start: None,
                    cds_end: None,
                },
            ],
        }
    }

    #[test]
    fn test_exon_len() {
        let exon = Exon { exon_number: 1, start: 100, end: 300 };
        assert_eq!(exon.len(), 200);
    }

    #[test]
    fn test_transcript_exonic_length() {
        let gene = sample_gene();
        let tx = &gene.transcripts[0];
        // (7668586-7668421) + (7670715-7670609) + (7673608-7673534)
        // = 165 + 106 + 74 = 345
        assert_eq!(tx.exonic_length(), 345);
    }

    #[test]
    fn test_transcript_n_exons() {
        let gene = sample_gene();
        assert_eq!(gene.transcripts[0].n_exons(), 3);
        assert_eq!(gene.transcripts[1].n_exons(), 1);
    }

    #[test]
    fn test_transcript_to_interval() {
        let gene = sample_gene();
        let tx = &gene.transcripts[0];
        let iv = tx.to_genomic_interval("chr17", Strand::Reverse);
        assert_eq!(iv.chrom, "chr17");
        assert_eq!(iv.start, 7668421);
        assert_eq!(iv.end, 7687490);
        assert_eq!(iv.strand, Strand::Reverse);
    }

    #[test]
    fn test_gene_len() {
        let gene = sample_gene();
        assert_eq!(gene.len(), 7687490 - 7668421);
    }

    #[test]
    fn test_gene_n_transcripts() {
        let gene = sample_gene();
        assert_eq!(gene.n_transcripts(), 2);
    }

    #[test]
    fn test_gene_to_interval() {
        let gene = sample_gene();
        let iv = gene.to_genomic_interval();
        assert_eq!(iv.chrom, "chr17");
        assert_eq!(iv.strand, Strand::Reverse);
    }

    #[test]
    fn test_gene_is_protein_coding() {
        let gene = sample_gene();
        assert!(gene.is_protein_coding());
    }

    #[test]
    fn test_annotated() {
        let gene = sample_gene();
        assert_eq!(gene.name(), "TP53");
    }

    #[test]
    fn test_summary() {
        let gene = sample_gene();
        assert_eq!(
            gene.summary(),
            "Gene: TP53 (chr17:7668421-7687490, -, protein_coding, 2 transcripts)"
        );
    }

    #[test]
    fn test_gene_type_display() {
        assert_eq!(GeneType::ProteinCoding.to_string(), "protein_coding");
        assert_eq!(GeneType::LncRNA.to_string(), "lncRNA");
        assert_eq!(GeneType::Other("snRNA".into()).to_string(), "snRNA");
    }
}
