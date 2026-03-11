use cyanea_proteomics::fdr::*;
use cyanea_proteomics::search::Psm;
use cyanea_proteomics::mztab::{write_mztab_psms, write_mztab_proteins};
use cyanea_proteomics::protein::ProteinGroup;

fn make_psm(id: &str, seq: &str, score: f64) -> Psm {
    Psm {
        spectrum_id: id.to_string(),
        peptide_sequence: seq.to_string(),
        xcorr: score * 0.1,
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
fn test_fdr_with_realistic_distribution() {
    // Simulate a realistic score distribution:
    // - 100 targets with scores 10-110
    // - 20 decoys with scores 5-25
    let targets: Vec<Psm> = (0..100)
        .map(|i| make_psm(&format!("t{}", i), &format!("PEP{}", i), 10.0 + i as f64))
        .collect();

    let decoys: Vec<Psm> = (0..20)
        .map(|i| make_psm(&format!("d{}", i), &format!("DEC{}", i), 5.0 + i as f64))
        .collect();

    let config = FdrConfig { threshold: 0.01, ..Default::default() };
    let results = estimate_fdr(&targets, &decoys, &config).unwrap();

    // High-scoring targets should pass
    let passing: Vec<&FdrResult> = results.iter().filter(|r| r.passes).collect();
    assert!(!passing.is_empty());

    // All passing should be targets
    assert!(passing.iter().all(|r| !r.is_decoy));

    // Q-values should be non-decreasing (walking down by score)
    for i in 1..results.len() {
        assert!(results[i].q_value >= results[i - 1].q_value - 1e-10);
    }
}

#[test]
fn test_fdr_stringency_levels() {
    let targets: Vec<Psm> = (0..50)
        .map(|i| make_psm(&format!("t{}", i), &format!("PEP{}", i), 10.0 + i as f64 * 2.0))
        .collect();

    let decoys: Vec<Psm> = (0..10)
        .map(|i| make_psm(&format!("d{}", i), &format!("DEC{}", i), 5.0 + i as f64 * 3.0))
        .collect();

    let passing_1 = filter_fdr(&targets, &decoys, 0.01).unwrap();
    let passing_5 = filter_fdr(&targets, &decoys, 0.05).unwrap();
    let passing_10 = filter_fdr(&targets, &decoys, 0.10).unwrap();

    // More stringent FDR should yield fewer PSMs
    assert!(passing_1.len() <= passing_5.len());
    assert!(passing_5.len() <= passing_10.len());
}

#[test]
fn test_fdr_summary_complete() {
    let targets: Vec<Psm> = (0..30)
        .map(|i| make_psm(&format!("t{}", i), &format!("PEP{}", i), 20.0 + i as f64))
        .collect();

    let decoys: Vec<Psm> = (0..5)
        .map(|i| make_psm(&format!("d{}", i), &format!("DEC{}", i), 15.0 + i as f64))
        .collect();

    let summary = fdr_summary(&targets, &decoys).unwrap();
    assert_eq!(summary.total_targets, 30);
    assert_eq!(summary.total_decoys, 5);
    assert!(summary.passing_1pct <= summary.passing_5pct);
    assert!(summary.passing_5pct <= 30);
}

#[test]
fn test_mztab_output_from_fdr() {
    let targets = vec![
        make_psm("s1", "PEPTIDEK", 50.0),
        make_psm("s2", "SEQUENCER", 40.0),
    ];
    let decoys = vec![
        make_psm("d1", "REVKDITP", 20.0),
    ];

    let config = FdrConfig { threshold: 0.05, ..Default::default() };
    let fdr_results = estimate_fdr(&targets, &decoys, &config).unwrap();

    let mztab = write_mztab_psms(&fdr_results);
    assert!(mztab.contains("PSH\tsequence"));
    assert!(mztab.contains("PEPTIDEK"));
    assert!(mztab.contains("SEQUENCER"));

    // Test protein mzTab output
    let groups = vec![
        ProteinGroup {
            accessions: vec!["P12345".into()],
            unique_peptides: vec!["PEPTIDEK".into(), "SEQUENCER".into()],
            shared_peptides: vec![],
            psm_count: 2,
            best_score: 50.0,
            coverage: 0.45,
            is_decoy: false,
        },
    ];

    let prot_mztab = write_mztab_proteins(&groups, None);
    assert!(prot_mztab.contains("PRT\tP12345"));
}

#[test]
fn test_all_decoys_below_targets() {
    // Edge case: all decoys score below all targets
    let targets: Vec<Psm> = (0..10)
        .map(|i| make_psm(&format!("t{}", i), &format!("PEP{}", i), 50.0 + i as f64))
        .collect();

    let decoys: Vec<Psm> = (0..5)
        .map(|i| make_psm(&format!("d{}", i), &format!("DEC{}", i), 1.0 + i as f64))
        .collect();

    let config = FdrConfig { threshold: 0.01, ..Default::default() };
    let results = estimate_fdr(&targets, &decoys, &config).unwrap();

    // All targets should pass since all decoys are below
    let passing: Vec<&FdrResult> = results.iter().filter(|r| r.passes).collect();
    assert_eq!(passing.len(), 10);
}
