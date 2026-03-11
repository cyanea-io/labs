//! Integration tests for motif discovery and scanning.

use cyanea_epi::motifs::{
    Motif, DiscoveryParams, scan_sequence, discover_motifs, parse_meme, write_meme,
    compare_motifs, motif_enrichment,
};

#[test]
fn test_motif_scanning_forward_strand() {
    // Create a motif (ACGT consensus)
    let pwm = vec![
        [0.9, 0.05, 0.03, 0.02],
        [0.05, 0.9, 0.03, 0.02],
        [0.03, 0.05, 0.9, 0.02],
        [0.02, 0.03, 0.05, 0.9],
    ];

    let motif = Motif::new("test_motif", pwm);

    // Scan a sequence
    let seq = b"NNNNACGTNNNNACGTNNNN";
    let matches = scan_sequence(seq, &motif, -5.0);

    assert!(!matches.is_empty());
    // Both forward and reverse strand matches are possible
    assert!(matches.iter().any(|m| m.position > 0));
}

#[test]
fn test_motif_scanning_reverse_strand() {
    let pwm = vec![
        [0.9, 0.05, 0.03, 0.02],
        [0.05, 0.9, 0.03, 0.02],
        [0.03, 0.05, 0.9, 0.02],
        [0.02, 0.03, 0.05, 0.9],
    ];

    let motif = Motif::new("test_motif", pwm);

    // ACGT reverse complement is ACGT
    let seq = b"NNNNACGTNNNN";
    let matches = scan_sequence(seq, &motif, -5.0);

    // Should have both strands represented
    assert!(!matches.is_empty());
}

#[test]
fn test_motif_discovery_from_peaks() {
    let peak_seqs = vec![
        b"ACGTACGTACGTACGTACGT".as_slice(),
        b"ACGTACGTACGTACGTACGT".as_slice(),
        b"ACGTACGTNNNNNNNNNNNN".as_slice(),
        b"NNNNACGTACGTACGTACGT".as_slice(),
    ];

    let params = DiscoveryParams {
        motif_width: 4,
        n_motifs: 3,
        background_freq: [0.25; 4],
    };

    let motifs = discover_motifs(&peak_seqs, &params).unwrap();

    assert!(!motifs.is_empty());
    assert!(motifs.len() <= params.n_motifs);
}

#[test]
fn test_meme_format_roundtrip() {
    let pwm = vec![
        [0.8, 0.1, 0.05, 0.05],
        [0.1, 0.8, 0.05, 0.05],
        [0.05, 0.05, 0.8, 0.1],
        [0.05, 0.05, 0.1, 0.8],
    ];

    let motif = Motif::new("test_motif", pwm);

    // Write to MEME format
    let meme_str = write_meme(&[motif.clone()]);
    assert!(meme_str.contains("MOTIF"));
    assert!(meme_str.contains("test_motif"));

    // Parse back
    let parsed = parse_meme(&meme_str).unwrap();
    assert_eq!(parsed.len(), 1);
    assert_eq!(parsed[0].name, "test_motif");
    assert_eq!(parsed[0].width(), motif.width());
}

#[test]
fn test_meme_parse_multiple_motifs() {
    let meme_content = r#"MEME version 4.0

MOTIF motif1 ACGT
letter-probability matrix: alength= 4 w= 4
 0.900000 0.050000 0.030000 0.020000
 0.050000 0.900000 0.030000 0.020000
 0.030000 0.050000 0.900000 0.020000
 0.020000 0.030000 0.050000 0.900000
//

MOTIF motif2 TGCA
letter-probability matrix: alength= 4 w= 4
 0.020000 0.030000 0.050000 0.900000
 0.050000 0.900000 0.030000 0.020000
 0.900000 0.050000 0.030000 0.020000
 0.030000 0.050000 0.900000 0.020000
//
"#;

    let motifs = parse_meme(meme_content).unwrap();
    assert_eq!(motifs.len(), 2);
    assert_eq!(motifs[0].name, "motif1");
    assert_eq!(motifs[1].name, "motif2");
}

#[test]
fn test_motif_comparison() {
    let pwm1 = vec![
        [0.9, 0.05, 0.03, 0.02],
        [0.05, 0.9, 0.03, 0.02],
        [0.03, 0.05, 0.9, 0.02],
        [0.02, 0.03, 0.05, 0.9],
    ];

    let pwm2 = vec![
        [0.9, 0.05, 0.03, 0.02],
        [0.05, 0.9, 0.03, 0.02],
        [0.03, 0.05, 0.9, 0.02],
        [0.02, 0.03, 0.05, 0.9],
    ];

    let pwm3 = vec![
        [0.25, 0.25, 0.25, 0.25],
        [0.25, 0.25, 0.25, 0.25],
        [0.25, 0.25, 0.25, 0.25],
        [0.25, 0.25, 0.25, 0.25],
    ];

    let motif1 = Motif::new("m1", pwm1);
    let motif2 = Motif::new("m2", pwm2);
    let motif3 = Motif::new("m3", pwm3);

    let sim_same = compare_motifs(&motif1, &motif2);
    let sim_diff = compare_motifs(&motif1, &motif3);

    // Identical motifs should have higher similarity
    assert!(sim_same > sim_diff);
}

#[test]
fn test_motif_enrichment_in_peaks() {
    let pwm = vec![
        [0.9, 0.05, 0.03, 0.02],
        [0.05, 0.9, 0.03, 0.02],
        [0.03, 0.05, 0.9, 0.02],
        [0.02, 0.03, 0.05, 0.9],
    ];

    let motif = Motif::new("test", pwm);

    // Peak sequences enriched for ACGT
    let peak_seqs = vec![
        b"ACGTACGTACGTACGT".as_slice(),
        b"AAACGTACGTACGTAA".as_slice(),
        b"ACGTNNNNACGTACGT".as_slice(),
    ];

    // Background sequences depleted for ACGT
    let bg_seqs = vec![
        b"GGGGGGGGGGGGGGGG".as_slice(),
        b"TTTTTTTTTTTTTTTT".as_slice(),
        b"GGGGTTTTGGGGTTTT".as_slice(),
    ];

    let (fold_enrichment, p_value) = motif_enrichment(&motif, &peak_seqs, &bg_seqs, -5.0).unwrap();

    assert!(fold_enrichment > 0.0);
    assert!(p_value >= 0.0 && p_value <= 1.0);
}

#[test]
fn test_motif_width_handling() {
    let pwm = vec![
        [0.9, 0.05, 0.03, 0.02],
        [0.05, 0.9, 0.03, 0.02],
    ];

    let motif = Motif::new("short", pwm);
    assert_eq!(motif.width(), 2);

    // Scanning shorter sequence
    let seq = b"AC";
    let matches = scan_sequence(seq, &motif, -5.0);
    assert!(!matches.is_empty());

    // Scanning longer sequence
    let seq_long = b"ACGTACGTACGT";
    let matches_long = scan_sequence(seq_long, &motif, -5.0);
    assert!(!matches_long.is_empty());
}

#[test]
fn test_empty_discovery() {
    let params = DiscoveryParams::default();
    let motifs = discover_motifs(&[], &params).unwrap();
    assert!(motifs.is_empty());
}
