//! Alignment demo datasets: pre-aligned sequences, reference pairs.

/// Three SARS-CoV-2 spike RBD sequences for multiple alignment demo.
///
/// Wuhan-Hu-1 (reference), Alpha (B.1.1.7), and Delta (B.1.617.2) variants.
/// Short subsequence (~60 aa) around key mutation sites.
pub fn spike_alignment_seqs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("Wuhan-Hu-1", b"YQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKSTN"),
        ("Alpha_B117",  b"YQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKSTN"),
        ("Delta_B1617", b"YQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKSTN"),
        ("Omicron_BA1", b"YQAGSTPCNGVKGFNCYFPLQSYGFHPTNGVGYQPYRVVVLSFELLHAPATVCGPKSTN"),
    ]
}

/// Hemoglobin alpha sequences from 5 species for cross-species alignment.
///
/// ~141 amino acids each — classic phylogenetics dataset.
pub fn hemoglobin_alpha() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("Human",      b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"),
        ("Chimpanzee", b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"),
        ("Dog",        b"MVLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKAHGKKVADALTLAVAHVDDMPQALSALSDLHAHKLRVDPVNFKLLSHCLLSTLAVHLPNDFTPAVHASLDKFLASVSTVLTSKYR"),
        ("Chicken",    b"MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSHGSAQIKAHGKKVANALIEAVNHIDDIAGALSKLSDLHAQKLRVDPVNFKFLGHCFLVVVAIHHPSVLTPEVHASLDKFLCAVGNVLTAKYR"),
        ("Zebrafish",  b"MSLSDKDKAAVKGLWAKISPKADDIGAEALGRMLTVYPQTKTYFAHWADLSPGSAPVKKHGITIMNQIDDCVGHMDDLFGFLTKLSELHATKLRVDPANFKILAHNLIVVIAAYFPAEFTPEIHLSVDKFLQQLALALAEKYR"),
    ]
}

/// A pre-computed pairwise alignment (Needleman-Wunsch) for demo purposes.
///
/// Human vs. chicken hemoglobin alpha, with gaps.
pub fn demo_pairwise_alignment() -> DemoAlignment {
    DemoAlignment {
        seq_a_name: "Human_HBA",
        seq_b_name: "Chicken_HBA",
        aligned_a: b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL-SFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        aligned_b: b"MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSHGSAQIKAHGKKVANALIEAVNHIDDIAGALSKLSDLHAQKLRVDPVNFKFLGHCFLVVVAI-HHPSVLTPEVHASLDKFLCAVGNVLTAKYR",
        score: 425,
        identity: 0.62,
        gaps: 2,
    }
}

/// A pre-computed alignment result.
#[derive(Debug, Clone)]
pub struct DemoAlignment {
    pub seq_a_name: &'static str,
    pub seq_b_name: &'static str,
    pub aligned_a: &'static [u8],
    pub aligned_b: &'static [u8],
    pub score: i64,
    pub identity: f64,
    pub gaps: usize,
}

/// Demo CIGAR strings for SAM/BAM examples.
pub fn demo_cigar_strings() -> Vec<(&'static str, &'static str)> {
    vec![
        ("perfect_match", "100M"),
        ("with_insertion", "50M2I48M"),
        ("with_deletion", "30M5D70M"),
        ("soft_clipped", "5S90M5S"),
        ("complex_realign", "20M1I10M2D15M1I52M"),
        ("spliced_rna", "75M5000N25M"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spike_alignment_seqs() {
        let seqs = spike_alignment_seqs();
        assert_eq!(seqs.len(), 4);
        // All same length (pre-aligned region)
        let len = seqs[0].1.len();
        assert!(seqs.iter().all(|s| s.1.len() == len || (s.1.len() as i64 - len as i64).abs() <= 1));
    }

    #[test]
    fn test_hemoglobin_alpha() {
        let seqs = hemoglobin_alpha();
        assert_eq!(seqs.len(), 5);
        assert_eq!(seqs[0].0, "Human");
        // All ~141 aa
        assert!(seqs.iter().all(|s| s.1.len() >= 139 && s.1.len() <= 145));
    }

    #[test]
    fn test_demo_pairwise_alignment() {
        let aln = demo_pairwise_alignment();
        assert_eq!(aln.aligned_a.len(), aln.aligned_b.len());
        assert!(aln.identity > 0.5);
    }

    #[test]
    fn test_demo_cigar_strings() {
        let cigars = demo_cigar_strings();
        assert_eq!(cigars.len(), 6);
        assert!(cigars[0].1.contains('M'));
    }
}
