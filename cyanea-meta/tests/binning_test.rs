/// Integration tests for binning module.

use cyanea_meta::binning::{tetranucleotide_frequency, bin_contigs, Contig, filter_bins};

#[test]
fn tetranucleotide_frequency_basic() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let tnf = tetranucleotide_frequency(seq).unwrap();

    assert_eq!(tnf.len(), 256);
    // All frequencies should be non-negative
    for &f in &tnf {
        assert!(f >= 0.0);
    }
    // Sum should be approximately 1.0
    let sum: f64 = tnf.iter().sum();
    assert!((sum - 1.0).abs() < 0.01);
}

#[test]
fn tnf_too_short() {
    let seq = b"ACG";
    assert!(tetranucleotide_frequency(seq).is_err());
}

#[test]
fn tnf_minimum_length() {
    let seq = b"ACGT";
    let tnf = tetranucleotide_frequency(seq).unwrap();
    assert_eq!(tnf.len(), 256);
}

#[test]
fn tnf_repeated_sequence() {
    let seq = b"AAAAAAAAAAAAAAAA";
    let tnf = tetranucleotide_frequency(seq).unwrap();

    // Only AAAA should have non-zero frequency
    let mut non_zero = 0;
    for (i, &f) in tnf.iter().enumerate() {
        if f > 0.0 {
            non_zero = i;
            assert!(f > 0.999); // Should be ~1.0
        }
    }
    // AAAA is encoded as 0
    assert_eq!(non_zero, 0);
}

#[test]
fn contig_creation_basic() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let contig = Contig::new("ctg1", seq, 10.5).unwrap();

    assert_eq!(contig.id, "ctg1");
    assert_eq!(contig.length, 32);
    assert!((contig.coverage - 10.5).abs() < 1e-10);
    assert_eq!(contig.tnf.len(), 256);
}

#[test]
fn contig_gc_content() {
    let seq = b"GGCCAATT".to_vec();
    let contig = Contig::new("ctg1", seq, 1.0).unwrap();

    // 50% GC
    assert!((contig.gc_content - 0.5).abs() < 1e-10);
}

#[test]
fn contig_empty_sequence_error() {
    let seq = vec![];
    assert!(Contig::new("ctg1", seq, 1.0).is_err());
}

#[test]
fn bin_contigs_basic() {
    let seq1 = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq2 = b"TGCATGCATGCATGCATGCATGCATGCATGCA".to_vec();
    let seq3 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();

    let contigs = vec![
        Contig::new("ctg1", seq1, 10.0).unwrap(),
        Contig::new("ctg2", seq2, 5.0).unwrap(),
        Contig::new("ctg3", seq3, 15.0).unwrap(),
    ];

    let bins = bin_contigs(&contigs, 3).unwrap();
    assert_eq!(bins.len(), 3);
    assert!(bins.iter().all(|b| !b.contig_ids.is_empty()));
}

#[test]
fn bin_contigs_fewer_bins_than_contigs() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let contigs = vec![
        Contig::new("ctg1", seq.clone(), 10.0).unwrap(),
        Contig::new("ctg2", seq.clone(), 5.0).unwrap(),
        Contig::new("ctg3", seq.clone(), 15.0).unwrap(),
        Contig::new("ctg4", seq.clone(), 8.0).unwrap(),
    ];

    let bins = bin_contigs(&contigs, 2).unwrap();
    assert_eq!(bins.len(), 2);
    // All contigs should be assigned
    let total_contigs: usize = bins.iter().map(|b| b.contig_ids.len()).sum();
    assert_eq!(total_contigs, 4);
}

#[test]
fn bin_contigs_zero_clusters_error() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let contigs = vec![Contig::new("ctg1", seq, 10.0).unwrap()];

    assert!(bin_contigs(&contigs, 0).is_err());
}

#[test]
fn bin_contigs_empty_error() {
    assert!(bin_contigs(&[], 2).is_err());
}

#[test]
fn filter_bins_basic() {
    use cyanea_meta::binning::Bin;

    let bins = vec![
        Bin::new("bin1", vec![], 0.9, 0.02),
        Bin::new("bin2", vec![], 0.5, 0.1),
        Bin::new("bin3", vec![], 0.95, 0.08),
    ];

    let filtered = filter_bins(&bins, 0.8, 0.05).unwrap();
    // Only bin1 should pass: completeness ≥ 0.8 and contamination ≤ 0.05
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].id, "bin1");
}

#[test]
fn filter_bins_all_pass() {
    use cyanea_meta::binning::Bin;

    let bins = vec![
        Bin::new("bin1", vec![], 0.9, 0.02),
        Bin::new("bin2", vec![], 0.85, 0.03),
    ];

    let filtered = filter_bins(&bins, 0.8, 0.05).unwrap();
    assert_eq!(filtered.len(), 2);
}

#[test]
fn filter_bins_none_pass() {
    use cyanea_meta::binning::Bin;

    let bins = vec![
        Bin::new("bin1", vec![], 0.5, 0.1),
        Bin::new("bin2", vec![], 0.4, 0.15),
    ];

    let filtered = filter_bins(&bins, 0.8, 0.05).unwrap();
    assert_eq!(filtered.len(), 0);
}

#[test]
fn filter_bins_invalid_completeness_error() {
    use cyanea_meta::binning::Bin;

    let bins = vec![Bin::new("bin1", vec![], 0.9, 0.02)];
    assert!(filter_bins(&bins, 1.5, 0.05).is_err());
}

#[test]
fn filter_bins_invalid_contamination_error() {
    use cyanea_meta::binning::Bin;

    let bins = vec![Bin::new("bin1", vec![], 0.9, 0.02)];
    assert!(filter_bins(&bins, 0.8, -0.1).is_err());
}

#[test]
fn bin_quality_score() {
    use cyanea_meta::binning::Bin;

    let bin = Bin::new("bin1", vec![], 0.8, 0.05);
    let score = bin.quality_score();
    // Score = completeness - 5 * contamination = 0.8 - 0.25 = 0.55
    assert!((score - 0.55).abs() < 1e-10);
}
