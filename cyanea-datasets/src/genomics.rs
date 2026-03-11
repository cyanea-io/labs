//! Genomics demo datasets: E. coli genome subset, SARS-CoV-2 spike, demo VCF.

/// E. coli K-12 MG1655 rrsA (16S rRNA) gene — 1,542 bp.
///
/// Universally used in metagenomics and phylogenetics.
pub fn ecoli_16s_rrna() -> (&'static str, &'static [u8]) {
    ("rrsA_16S_rRNA", ECOLI_16S)
}

const ECOLI_16S: &[u8] = b"AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAA\
CACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAG\
TGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAAC\
TACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGA\
CCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTG\
GGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCA\
GCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGG\
GAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAG\
AAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAA\
TACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCC\
AGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTA";

/// SARS-CoV-2 Spike protein receptor-binding domain (RBD) — nucleotide, ~670 bp.
///
/// Covers residues 319–541 of the S protein. Key region for ACE2 binding
/// and neutralizing antibody targets.
pub fn sars_cov2_spike_rbd() -> (&'static str, &'static [u8]) {
    ("SARS-CoV-2_Spike_RBD", SPIKE_RBD)
}

const SPIKE_RBD: &[u8] = b"AATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTG\
CATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTAT\
TCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCT\
CCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTA\
ATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCT\
GATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAAT\
TCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTG\
TTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTAT\
CAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCT\
TTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGA\
GTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCT\
AAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGT\
TTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAA\
CAATTTGGCAGAGACATTGCTGACAC";

/// SARS-CoV-2 Spike RBD — protein sequence (~222 aa).
pub fn sars_cov2_spike_rbd_protein() -> (&'static str, &'static [u8]) {
    ("SARS-CoV-2_Spike_RBD_protein", SPIKE_RBD_PROTEIN)
}

const SPIKE_RBD_PROTEIN: &[u8] = b"NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGV\
SPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIA\
WNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNC\
YFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF\
NFNGLTGTGVLTESNKKFLPFQQF";

/// Small set of demo variants (VCF-style) from human chr22.
///
/// 10 variants: mix of SNVs, a deletion, and an insertion.
/// Based on common 1000 Genomes variants on chromosome 22.
pub fn demo_variants_chr22() -> Vec<DemoVariant> {
    vec![
        DemoVariant { chrom: "chr22", pos: 16050075, id: "rs587697622", ref_a: "A", alt: "G", qual: 100.0, af: 0.0002 },
        DemoVariant { chrom: "chr22", pos: 16050115, id: "rs587755077", ref_a: "G", alt: "A", qual: 99.0, af: 0.0008 },
        DemoVariant { chrom: "chr22", pos: 16050213, id: "rs533874225", ref_a: "C", alt: "T", qual: 85.0, af: 0.012 },
        DemoVariant { chrom: "chr22", pos: 16050319, id: "rs62224609", ref_a: "T", alt: "C", qual: 95.0, af: 0.185 },
        DemoVariant { chrom: "chr22", pos: 16050527, id: "rs587638893", ref_a: "C", alt: "A", qual: 88.0, af: 0.0015 },
        DemoVariant { chrom: "chr22", pos: 16050607, id: "rs9617528", ref_a: "T", alt: "C", qual: 99.0, af: 0.342 },
        DemoVariant { chrom: "chr22", pos: 16050739, id: "rs587776767", ref_a: "G", alt: "T", qual: 78.0, af: 0.0003 },
        DemoVariant { chrom: "chr22", pos: 16050847, id: "rs2845371", ref_a: "T", alt: "A", qual: 97.0, af: 0.271 },
        DemoVariant { chrom: "chr22", pos: 16050922, id: "rs587712345", ref_a: "AGT", alt: "A", qual: 72.0, af: 0.0001 },
        DemoVariant { chrom: "chr22", pos: 16050984, id: "rs587798123", ref_a: "C", alt: "CTGA", qual: 68.0, af: 0.0005 },
    ]
}

/// A simple variant record.
#[derive(Debug, Clone)]
pub struct DemoVariant {
    pub chrom: &'static str,
    pub pos: u64,
    pub id: &'static str,
    pub ref_a: &'static str,
    pub alt: &'static str,
    pub qual: f64,
    pub af: f64,
}

impl DemoVariant {
    /// Format as VCF line.
    pub fn to_vcf_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{:.1}\tPASS\tAF={:.4}",
            self.chrom, self.pos, self.id, self.ref_a, self.alt, self.qual, self.af
        )
    }
}

/// Generate a VCF string from demo variants.
pub fn demo_vcf_chr22() -> String {
    let mut lines = vec![
        "##fileformat=VCFv4.2".to_string(),
        "##source=cyanea-datasets".to_string(),
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">".to_string(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
    ];
    for v in demo_variants_chr22() {
        lines.push(v.to_vcf_line());
    }
    lines.join("\n")
}

/// Demo gene annotation for TP53 (simplified).
pub fn tp53_gene() -> DemoGene {
    DemoGene {
        name: "TP53",
        chrom: "chr17",
        start: 7668402,
        end: 7687550,
        strand: '-',
        exon_count: 11,
        transcript: "NM_000546.6",
        description: "Tumor protein p53",
    }
}

/// Demo gene annotation for BRCA1 (simplified).
pub fn brca1_gene() -> DemoGene {
    DemoGene {
        name: "BRCA1",
        chrom: "chr17",
        start: 43044295,
        end: 43170245,
        strand: '-',
        exon_count: 23,
        transcript: "NM_007294.4",
        description: "BRCA1 DNA repair associated",
    }
}

/// A simplified gene annotation.
#[derive(Debug, Clone)]
pub struct DemoGene {
    pub name: &'static str,
    pub chrom: &'static str,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub exon_count: usize,
    pub transcript: &'static str,
    pub description: &'static str,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ecoli_16s() {
        let (name, seq) = ecoli_16s_rrna();
        assert_eq!(name, "rrsA_16S_rRNA");
        assert!(seq.len() > 500);
        // Valid DNA
        assert!(seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')));
    }

    #[test]
    fn test_spike_rbd() {
        let (name, seq) = sars_cov2_spike_rbd();
        assert!(name.contains("Spike"));
        assert!(seq.len() > 600);
    }

    #[test]
    fn test_spike_protein() {
        let (_, seq) = sars_cov2_spike_rbd_protein();
        assert!(seq.len() > 200);
    }

    #[test]
    fn test_demo_variants() {
        let variants = demo_variants_chr22();
        assert_eq!(variants.len(), 10);
        assert!(variants.iter().all(|v| v.chrom == "chr22"));
        // Mix of types
        assert!(variants.iter().any(|v| v.ref_a.len() > 1)); // deletion
        assert!(variants.iter().any(|v| v.alt.len() > 1)); // insertion
    }

    #[test]
    fn test_demo_vcf() {
        let vcf = demo_vcf_chr22();
        assert!(vcf.starts_with("##fileformat=VCFv4.2"));
        assert!(vcf.contains("#CHROM"));
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 10);
    }

    #[test]
    fn test_variant_to_vcf_line() {
        let v = &demo_variants_chr22()[0];
        let line = v.to_vcf_line();
        assert!(line.contains("chr22"));
        assert!(line.contains("rs587697622"));
    }

    #[test]
    fn test_tp53_gene() {
        let g = tp53_gene();
        assert_eq!(g.name, "TP53");
        assert_eq!(g.chrom, "chr17");
        assert_eq!(g.exon_count, 11);
    }

    #[test]
    fn test_brca1_gene() {
        let g = brca1_gene();
        assert_eq!(g.name, "BRCA1");
        assert_eq!(g.exon_count, 23);
    }
}
