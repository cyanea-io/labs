/// Integration tests for composition module.

use cyanea_meta::composition::{clr_transform, ilr_transform, differential_abundance, ancom};

#[test]
fn clr_transform_basic() {
    let abundances = vec![vec![10.0, 20.0, 30.0]];
    let result = clr_transform(&abundances).unwrap();

    assert_eq!(result.len(), 1);
    // Sum of CLR values should be approximately 0 (property of CLR)
    let sum: f64 = result[0].values.values().sum();
    assert!(sum.abs() < 1e-10);
}

#[test]
fn clr_transform_multiple_samples() {
    let abundances = vec![
        vec![100.0, 50.0, 25.0],
        vec![10.0, 20.0, 30.0],
    ];
    let result = clr_transform(&abundances).unwrap();

    assert_eq!(result.len(), 2);
    for transform in &result {
        let sum: f64 = transform.values.values().sum();
        assert!(sum.abs() < 1e-10, "sum should be ~0: {}", sum);
    }
}

#[test]
fn clr_rejects_zero_abundances() {
    let abundances = vec![vec![0.0, 20.0, 30.0]];
    assert!(clr_transform(&abundances).is_err());
}

#[test]
fn clr_rejects_negative_abundances() {
    let abundances = vec![vec![-10.0, 20.0, 30.0]];
    assert!(clr_transform(&abundances).is_err());
}

#[test]
fn ilr_transform_valid() {
    let abundances = vec![vec![10.0, 20.0, 30.0]];
    let result = ilr_transform(&abundances).unwrap();

    assert_eq!(result.len(), 1);
    assert_eq!(result[0].len(), 3);
}

#[test]
fn ilr_transform_multiple_samples() {
    let abundances = vec![
        vec![100.0, 50.0, 25.0],
        vec![10.0, 20.0, 30.0],
    ];
    let result = ilr_transform(&abundances).unwrap();

    assert_eq!(result.len(), 2);
    assert_eq!(result[0].len(), 3);
    assert_eq!(result[1].len(), 3);
}

#[test]
fn differential_abundance_basic() {
    let group1 = vec![
        vec![100.0, 10.0, 5.0],
        vec![95.0, 15.0, 10.0],
    ];
    let group2 = vec![
        vec![10.0, 100.0, 50.0],
        vec![15.0, 95.0, 45.0],
    ];

    let results = differential_abundance(&group1, &group2).unwrap();

    assert!(!results.is_empty());
    assert_eq!(results.len(), 3); // 3 taxa

    // All p-values should be in [0, 1]
    for r in &results {
        assert!(r.p_value >= 0.0 && r.p_value <= 1.0);
        assert!(r.t_statistic.is_finite());
    }
}

#[test]
fn differential_abundance_sorted_by_pvalue() {
    let group1 = vec![vec![100.0, 10.0]];
    let group2 = vec![vec![10.0, 100.0]];

    let results = differential_abundance(&group1, &group2).unwrap();

    // Results should be sorted by p-value (ascending)
    for w in results.windows(2) {
        assert!(w[0].p_value <= w[1].p_value);
    }
}

#[test]
fn differential_abundance_empty_group_error() {
    let group1 = vec![vec![100.0, 10.0]];
    let group2: Vec<Vec<f64>> = vec![];

    assert!(differential_abundance(&group1, &group2).is_err());
}

#[test]
fn differential_abundance_mismatched_dimensions_error() {
    let group1 = vec![vec![100.0, 10.0]];
    let group2 = vec![vec![10.0, 100.0, 50.0]];

    assert!(differential_abundance(&group1, &group2).is_err());
}

#[test]
fn ancom_basic() {
    let group1 = vec![
        vec![100.0, 10.0],
        vec![95.0, 15.0],
    ];
    let group2 = vec![
        vec![10.0, 100.0],
        vec![15.0, 95.0],
    ];

    let results = ancom(&group1, &group2).unwrap();

    assert_eq!(results.len(), 2); // 2 taxa
    for r in &results {
        assert!(r.w_statistic <= 1); // W should be 0 or 1 for 2 taxa
    }
}

#[test]
fn ancom_sorted_by_w_statistic() {
    let group1 = vec![
        vec![100.0, 10.0, 5.0],
        vec![95.0, 15.0, 10.0],
    ];
    let group2 = vec![
        vec![10.0, 100.0, 50.0],
        vec![15.0, 95.0, 45.0],
    ];

    let results = ancom(&group1, &group2).unwrap();

    // Results should be sorted by W-statistic (descending)
    for w in results.windows(2) {
        assert!(w[0].w_statistic >= w[1].w_statistic);
    }
}

#[test]
fn ancom_empty_group_error() {
    let group1 = vec![vec![100.0, 10.0]];
    let group2: Vec<Vec<f64>> = vec![];

    assert!(ancom(&group1, &group2).is_err());
}

#[test]
fn ancom_mismatched_taxa_error() {
    let group1 = vec![vec![100.0, 10.0, 50.0]];
    let group2 = vec![vec![10.0, 100.0]];

    assert!(ancom(&group1, &group2).is_err());
}
