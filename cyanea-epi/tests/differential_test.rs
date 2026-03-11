//! Integration tests for differential binding analysis.

use cyanea_epi::differential::{
    differential_peaks, count_reads_in_peaks, ma_plot_data,
};

#[test]
fn test_differential_analysis_clear_signal() {
    // Create count matrix with clear differential signal
    // Samples 0-1: group 1 (high counts in peaks 0-1, low in peak 2)
    // Samples 2-3: group 2 (low counts in peaks 0-1, high in peak 2)
    let count_matrix = vec![
        vec![100, 110, 10, 5],   // Peak 1: enriched in group 1
        vec![50, 60, 10, 15],    // Peak 2: enriched in group 1
        vec![5, 10, 100, 95],    // Peak 3: enriched in group 2
    ];

    let conditions = vec![0, 0, 1, 1];
    let region_ids = vec![
        "peak1".to_string(),
        "peak2".to_string(),
        "peak3".to_string(),
    ];

    let results = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0).unwrap();

    assert_eq!(results.len(), 3);

    // Results should contain valid log2_fc values
    for result in &results {
        assert!(!result.log2_fc.is_nan());
        assert!(result.log2_fc.is_finite());
    }
}

#[test]
fn test_differential_no_difference() {
    // Same counts across conditions
    let count_matrix = vec![
        vec![50, 50, 50, 50],
        vec![100, 100, 100, 100],
    ];

    let conditions = vec![0, 0, 1, 1];
    let region_ids = vec!["peak1".to_string(), "peak2".to_string()];

    let results = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0).unwrap();

    assert_eq!(results.len(), 2);

    // log2_fc should be ~0 (fold change = 1)
    for result in &results {
        assert!(result.log2_fc.abs() < 0.5);
        // p-value should be high (not significant)
        assert!(result.p_value > 0.05);
    }
}

#[test]
fn test_count_reads_in_peaks_basic() {
    let reads = vec![
        ("chr1".to_string(), 100, 150),
        ("chr1".to_string(), 110, 160),
        ("chr1".to_string(), 200, 250),
        ("chr1".to_string(), 600, 650),
    ];

    let peaks = vec![(100, 200), (200, 300), (500, 700)];
    let counts = count_reads_in_peaks(&reads, &peaks);

    assert_eq!(counts.len(), 3);
    assert_eq!(counts[0], 2); // Two reads in peak 1
    assert_eq!(counts[1], 1); // One read in peak 2
    assert_eq!(counts[2], 1); // One read in peak 3
}

#[test]
fn test_count_reads_partial_overlap() {
    let reads = vec![
        ("chr1".to_string(), 50, 150),  // Spans into peak
        ("chr1".to_string(), 150, 250), // Spans out of peak
    ];

    let peaks = vec![(100, 200)];
    let counts = count_reads_in_peaks(&reads, &peaks);

    assert_eq!(counts[0], 2); // Both reads overlap
}

#[test]
fn test_ma_plot_data_generation() {
    let results = vec![
        cyanea_epi::differential::DiffResult {
            region_id: "peak1".to_string(),
            log2_fc: 2.0,
            p_value: 0.001,
            q_value: 0.01,
            mean_count_1: 100.0,
            mean_count_2: 25.0,
        },
        cyanea_epi::differential::DiffResult {
            region_id: "peak2".to_string(),
            log2_fc: -1.0,
            p_value: 0.01,
            q_value: 0.05,
            mean_count_1: 10.0,
            mean_count_2: 20.0,
        },
    ];

    let data = ma_plot_data(&results);

    assert_eq!(data.len(), 2);
    assert!(data[0].0 > 0.0); // M (log-ratio)
    assert!(data[0].1 > 0.0); // A (mean)
}

#[test]
fn test_size_factor_normalization() {
    // Count matrix with different sequencing depths
    let count_matrix = vec![
        vec![10, 100],  // Peak 1: low counts in sample 1, high in sample 2
        vec![20, 200],  // Peak 2: proportionally same as peak 1
        vec![100, 50],  // Peak 3: opposite pattern
    ];

    let conditions = vec![0, 1];
    let region_ids = vec!["peak1".to_string(), "peak2".to_string(), "peak3".to_string()];

    let results = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0).unwrap();

    // Peaks 1 and 2 have same proportions (should give similar results)
    let fc1 = results[0].log2_fc;
    let fc2 = results[1].log2_fc;

    // They should be similar after normalization
    assert!((fc1 - fc2).abs() < 1.0);
}

#[test]
fn test_bh_correction() {
    // Create differential results with multiple tests
    let count_matrix = vec![
        vec![100, 110, 10, 5],
        vec![50, 60, 10, 15],
        vec![5, 10, 100, 95],
        vec![100, 100, 100, 100],
        vec![50, 50, 50, 50],
    ];

    let conditions = vec![0, 0, 1, 1];
    let region_ids = vec![
        "peak1".to_string(),
        "peak2".to_string(),
        "peak3".to_string(),
        "peak4".to_string(),
        "peak5".to_string(),
    ];

    let results = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0).unwrap();

    // Q-values should be monotonically non-decreasing when sorted by p-value
    let mut sorted_results = results.clone();
    sorted_results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap());

    for i in 1..sorted_results.len() {
        assert!(
            sorted_results[i].q_value >= sorted_results[i - 1].q_value,
            "Q-values not monotonic"
        );
    }
}

#[test]
fn test_invalid_conditions_error() {
    let count_matrix = vec![vec![10, 20, 30]];
    let conditions = vec![0, 2]; // Invalid condition (must be 0 or 1)
    let region_ids = vec!["peak1".to_string()];

    let result = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0);
    assert!(result.is_err());
}

#[test]
fn test_mismatched_condition_count_error() {
    let count_matrix = vec![vec![10, 20, 30]];
    let conditions = vec![0, 1]; // Only 2 conditions but 3 samples
    let region_ids = vec!["peak1".to_string()];

    let result = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0);
    assert!(result.is_err());
}

#[test]
fn test_empty_group_error() {
    let count_matrix = vec![vec![10, 20, 30]];
    let conditions = vec![0, 0, 0]; // All group 0, no group 1
    let region_ids = vec!["peak1".to_string()];

    let result = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0);
    assert!(result.is_err());
}

#[test]
fn test_pseudocount_handling() {
    // Very low count data
    let count_matrix = vec![
        vec![1, 2, 0, 1],
        vec![0, 1, 1, 0],
    ];

    let conditions = vec![0, 0, 1, 1];
    let region_ids = vec!["peak1".to_string(), "peak2".to_string()];

    // Should not crash with low counts or zeros
    let results = differential_peaks(&count_matrix, &conditions, &region_ids, 1.0).unwrap();
    assert_eq!(results.len(), 2);

    for result in &results {
        assert!(!result.log2_fc.is_nan());
        assert!(!result.p_value.is_nan());
    }
}
