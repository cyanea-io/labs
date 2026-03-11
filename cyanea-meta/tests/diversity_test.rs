/// Integration tests for diversity module.

use cyanea_meta::diversity::{
    alpha_diversity, beta_diversity_matrix, rarefaction_curve, rarefy,
};

#[test]
fn alpha_diversity_basic() {
    let counts = vec![100, 50, 25, 10];
    let ad = alpha_diversity(&counts).unwrap();

    assert_eq!(ad.observed_species, 4);
    assert!(ad.shannon > 0.0);
    assert!(ad.simpson >= 0.0 && ad.simpson <= 1.0);
    assert!(ad.chao1 >= 4.0); // At least the observed count
    assert!(ad.ace >= 4.0);
    assert!(ad.fisher_alpha > 0.0);
}

#[test]
fn shannon_uniform_distribution() {
    // Uniform distribution with 5 species of equal abundance
    let counts = vec![10, 10, 10, 10, 10];
    let ad = alpha_diversity(&counts).unwrap();

    let expected_shannon = (5.0f64).ln();
    assert!((ad.shannon - expected_shannon).abs() < 1e-10);
}

#[test]
fn shannon_single_species() {
    let counts = vec![100];
    let ad = alpha_diversity(&counts).unwrap();
    assert!((ad.shannon - 0.0).abs() < 1e-10);
}

#[test]
fn simpson_single_species_is_one() {
    let counts = vec![100];
    let ad = alpha_diversity(&counts).unwrap();
    assert!((ad.simpson - 1.0).abs() < 1e-10);
}

#[test]
fn simpson_two_equal_species() {
    let counts = vec![50, 50];
    let ad = alpha_diversity(&counts).unwrap();
    // D = 2 * 49 * 49 / (100 * 99) = 4802 / 9900 ≈ 0.485
    assert!(ad.simpson > 0.4 && ad.simpson < 0.6);
}

#[test]
fn chao1_with_singletons() {
    let counts = vec![10, 5, 3, 2, 1, 1];
    let ad = alpha_diversity(&counts).unwrap();
    // f1 = 2 (two singletons), f2 = 1 (one doubleton)
    // Chao1 = 6 + 2²/(2*1) = 6 + 2 = 8
    assert!(ad.chao1 >= 8.0);
}

#[test]
fn beta_diversity_identical_samples() {
    let s1 = vec![10, 20, 30];
    let s2 = vec![10, 20, 30];

    let matrix = beta_diversity_matrix(&[&s1, &s2]).unwrap();
    assert_eq!(matrix.distances.len(), 2);
    assert!(matrix.distances[0][1].abs() < 1e-10);
}

#[test]
fn beta_diversity_disjoint_samples() {
    let s1 = vec![10, 0, 0];
    let s2 = vec![0, 20, 30];

    let matrix = beta_diversity_matrix(&[&s1, &s2]).unwrap();
    assert!((matrix.distances[0][1] - 1.0).abs() < 1e-10);
}

#[test]
fn beta_diversity_symmetric() {
    let s1 = vec![10, 20, 0, 5];
    let s2 = vec![0, 15, 30, 0];
    let s3 = vec![5, 10, 10, 5];

    let matrix = beta_diversity_matrix(&[&s1, &s2, &s3]).unwrap();

    // Matrix should be symmetric
    for i in 0..3 {
        for j in 0..3 {
            assert!((matrix.distances[i][j] - matrix.distances[j][i]).abs() < 1e-10);
        }
    }
}

#[test]
fn rarefaction_curve_monotonic() {
    let counts = vec![100, 50, 25, 10];
    let total: usize = counts.iter().map(|&c| c as usize).sum();
    let steps: Vec<usize> = (1..=5).map(|i| i * total / 5).collect();

    let curve = rarefaction_curve(&counts, &steps).unwrap();

    // Rarefaction curve should be non-decreasing
    for w in curve.windows(2) {
        assert!(
            w[1].1 >= w[0].1 - 1e-10,
            "rarefaction not monotonic: {:?}",
            curve
        );
    }
}

#[test]
fn rarefaction_curve_at_full_depth() {
    let counts = vec![100, 50, 25];
    let total: usize = counts.iter().map(|&c| c as usize).sum();

    let curve = rarefaction_curve(&counts, &[total]).unwrap();
    assert_eq!(curve.len(), 1);
    assert_eq!(curve[0].0, total);
    // At full depth, expected richness should equal observed richness
    assert!((curve[0].1 - 3.0).abs() < 1e-10);
}

#[test]
fn rarefy_to_target_depth() {
    let counts = vec![100, 50, 25];
    let target = 100;

    let rarefied = rarefy(&counts, target).unwrap();
    let total: usize = rarefied.iter().map(|&c| c as usize).sum();
    assert_eq!(total, target);
}

#[test]
fn rarefy_preserves_proportions() {
    let counts = vec![100, 50];
    let target = 75;

    let rarefied = rarefy(&counts, target).unwrap();
    // Approximate proportions should be similar: 100:50 ≈ 2:1
    let ratio = rarefied[0] as f64 / rarefied[1].max(1) as f64;
    assert!(ratio > 1.5 && ratio < 2.5); // Some tolerance due to rounding
}

#[test]
fn rarefy_depth_exceeds_total_error() {
    let counts = vec![100, 50];
    assert!(rarefy(&counts, 200).is_err());
}

#[test]
fn empty_alpha_diversity_error() {
    assert!(alpha_diversity(&[]).is_err());
}

#[test]
fn zero_counts_alpha_diversity_error() {
    assert!(alpha_diversity(&[0, 0, 0]).is_err());
}

#[test]
fn insufficient_samples_beta_diversity_error() {
    let s1 = vec![10, 20];
    assert!(beta_diversity_matrix(&[&s1]).is_err());
}

#[test]
fn empty_steps_rarefaction_error() {
    let counts = vec![100, 50];
    assert!(rarefaction_curve(&counts, &[]).is_err());
}

#[test]
fn step_exceeds_total_rarefaction_error() {
    let counts = vec![100, 50]; // total = 150
    assert!(rarefaction_curve(&counts, &[200]).is_err());
}
