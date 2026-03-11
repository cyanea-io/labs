use cyanea_proteomics::protein::*;
use cyanea_proteomics::quantification::*;
use cyanea_proteomics::search::Psm;
use cyanea_proteomics::spectrum::*;

fn make_psm(seq: &str, spec_id: &str, score: f64) -> Psm {
    Psm {
        spectrum_id: spec_id.to_string(),
        peptide_sequence: seq.to_string(),
        xcorr: 1.0,
        hyperscore: score,
        matched_b: 5,
        matched_y: 5,
        total_b: 6,
        total_y: 6,
        delta_mass: 0.0,
        charge: 2,
        is_decoy: false,
    }
}

#[test]
fn test_spectral_counting_multiple_groups() {
    let groups = vec![
        ProteinGroup {
            accessions: vec!["P1".into()],
            unique_peptides: vec!["AAAK".into(), "BBBK".into()],
            shared_peptides: vec![],
            psm_count: 0,
            best_score: 0.0,
            coverage: 0.0,
            is_decoy: false,
        },
        ProteinGroup {
            accessions: vec!["P2".into()],
            unique_peptides: vec!["CCCK".into()],
            shared_peptides: vec![],
            psm_count: 0,
            best_score: 0.0,
            coverage: 0.0,
            is_decoy: false,
        },
    ];

    let psms = vec![
        make_psm("AAAK", "s1", 20.0),
        make_psm("AAAK", "s2", 18.0),
        make_psm("BBBK", "s3", 15.0),
        make_psm("CCCK", "s4", 10.0),
        make_psm("CCCK", "s5", 8.0),
        make_psm("CCCK", "s6", 6.0),
    ];

    let quants = spectral_counting(&groups, &psms);
    assert_eq!(quants.len(), 2);
    assert_eq!(quants[0].raw_value, 3.0); // P1: AAAK(2) + BBBK(1)
    assert_eq!(quants[1].raw_value, 3.0); // P2: CCCK(3)
}

#[test]
fn test_tmt_quantification_pipeline() {
    // Create spectra with TMT-like reporter ions
    let reporters = TmtPlex::Tmt6.reporter_mzs();
    let channel_values = vec![100.0, 200.0, 300.0, 150.0, 250.0, 175.0];

    let peaks: Vec<Peak> = reporters.iter()
        .zip(channel_values.iter())
        .map(|(&mz, &intensity)| Peak { mz, intensity })
        .chain(std::iter::once(Peak { mz: 300.0, intensity: 5000.0 }))
        .collect();

    let spec = MassSpectrum::new("tmt_scan", MsLevel::Ms2, 100.0, peaks).unwrap();
    let quants = quantify_tmt(&[spec], TmtPlex::Tmt6, 0.02);

    assert_eq!(quants.len(), 1);
    assert_eq!(quants[0].channel_intensities.len(), 6);
    assert!((quants[0].channel_intensities[0] - 100.0).abs() < 1e-5);
    assert!((quants[0].channel_intensities[2] - 300.0).abs() < 1e-5);
}

#[test]
fn test_normalization_pipeline() {
    let mut quants: Vec<ProteinQuant> = (1..=10)
        .map(|i| ProteinQuant {
            accessions: vec![format!("P{}", i)],
            raw_value: i as f64 * 100.0,
            normalized_value: 0.0,
            peptide_count: 2,
            psm_count: 3,
        })
        .collect();

    // Median normalize
    median_normalize(&mut quants);
    // Median of 100,200,...,1000 = 550
    assert!((quants[0].normalized_value - 100.0 / 550.0).abs() < 1e-10);
    assert!((quants[4].normalized_value - 500.0 / 550.0).abs() < 1e-10);

    // Log2 transform
    log2_transform(&mut quants);
    assert!((quants[0].normalized_value - (100.0f64).log2()).abs() < 1e-10);
}

#[test]
fn test_protein_inference_full() {
    let psms = vec![
        make_psm("AAAK", "s1", 50.0),
        make_psm("BBBK", "s2", 40.0),
        make_psm("CCCK", "s3", 30.0),
        make_psm("DDDK", "s4", 20.0),
    ];

    let proteins = vec![
        ProteinEntry::new("Prot_A", b"MAAAKBBBK"),     // has AAAK, BBBK
        ProteinEntry::new("Prot_B", b"MBBBKCCCKDDDK"), // has BBBK, CCCK, DDDK
    ];

    let groups = infer_proteins(&psms, &proteins).unwrap();
    assert!(!groups.is_empty());

    // Both proteins should be needed (parsimony)
    let total_accessions: usize = groups.iter().map(|g| g.accessions.len()).sum();
    assert_eq!(total_accessions, 2);

    // Quantify by spectral counting
    let quants = spectral_counting(&groups, &psms);
    assert_eq!(quants.len(), groups.len());
    for q in &quants {
        assert!(q.raw_value > 0.0);
    }
}
