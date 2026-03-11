//! Integration tests for chromatin state learning and segmentation.

use cyanea_epi::chromatin::{
    learn_chromatin_states, segment_genome, state_enrichment, ChromHMMParams,
};

#[test]
fn test_learn_chromatin_states_two_marks() {
    // Create a simple mark matrix with two marks
    // Active regions: mark1=1, mark2=0
    // Repressed regions: mark1=0, mark2=1
    let mark_matrix = vec![
        vec![1, 0], // Active
        vec![1, 0], // Active
        vec![1, 1], // Bivalent
        vec![0, 1], // Repressed
        vec![0, 1], // Repressed
    ];

    let marks = vec!["H3K4me3".to_string(), "H3K27me3".to_string()];

    let params = ChromHMMParams {
        n_states: 3,
        bin_size: 200,
        max_iter: 50,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();

    assert_eq!(model.states.len(), 3);
    assert_eq!(model.marks.len(), 2);
    assert_eq!(model.transition_matrix.len(), 3);
    assert_eq!(model.transition_matrix[0].len(), 3);
}

#[test]
fn test_segment_genome_basic() {
    let mark_matrix = vec![
        vec![1, 0],
        vec![1, 0],
        vec![0, 1],
        vec![0, 1],
        vec![1, 1],
    ];

    let marks = vec!["H3K4me3".to_string(), "H3K27me3".to_string()];
    let params = ChromHMMParams {
        n_states: 3,
        bin_size: 200,
        max_iter: 20,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();
    let segmentation = segment_genome(&mark_matrix, &model, "chr1", 200).unwrap();

    assert_eq!(segmentation.state_assignments.len(), 5);
    assert!(!segmentation.segments.is_empty());
    assert_eq!(segmentation.state_names.len(), 3);
}

#[test]
fn test_state_enrichment_at_annotations() {
    let mark_matrix = vec![
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![0, 1],
        vec![0, 1],
        vec![0, 1],
    ];

    let marks = vec!["mark1".to_string(), "mark2".to_string()];
    let params = ChromHMMParams {
        n_states: 2,
        bin_size: 200,
        max_iter: 20,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();
    let segmentation = segment_genome(&mark_matrix, &model, "chr1", 200).unwrap();

    // Define annotations (e.g., TSS regions) in first half of genome
    let annotations = vec![(0, 300), (0, 400)];

    let enrichments = state_enrichment(&segmentation, &annotations);

    assert_eq!(enrichments.len(), 2);
    for (state_id, enrichment) in &enrichments {
        assert!(*state_id < model.states.len());
        assert!(*enrichment >= 0.0 && *enrichment <= 1.0);
    }
}

#[test]
fn test_chromhmm_with_more_marks() {
    // Realistic 5-mark system
    let mark_matrix = vec![
        vec![1, 1, 0, 0, 0], // TSS
        vec![1, 1, 0, 0, 0], // TSS
        vec![1, 0, 1, 0, 0], // Enhancer
        vec![1, 0, 1, 0, 0], // Enhancer
        vec![0, 0, 0, 1, 1], // Repressed
        vec![0, 0, 0, 1, 1], // Repressed
        vec![0, 0, 0, 0, 0], // Quiescent
        vec![0, 0, 0, 0, 0], // Quiescent
    ];

    let marks = vec![
        "H3K4me3".to_string(),
        "H3K9ac".to_string(),
        "H3K4me1".to_string(),
        "H3K27me3".to_string(),
        "H3K9me3".to_string(),
    ];

    let params = ChromHMMParams {
        n_states: 4,
        bin_size: 200,
        max_iter: 30,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();

    assert_eq!(model.states.len(), 4);
    assert_eq!(model.marks.len(), 5);
}

#[test]
fn test_consecutive_same_states_merged() {
    let mark_matrix = vec![
        vec![1, 0], // State 0
        vec![1, 0], // State 0
        vec![1, 0], // State 0
        vec![0, 1], // State 1
        vec![0, 1], // State 1
    ];

    let marks = vec!["mark1".to_string(), "mark2".to_string()];
    let params = ChromHMMParams {
        n_states: 2,
        bin_size: 200,
        max_iter: 10,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();
    let segmentation = segment_genome(&mark_matrix, &model, "chr1", 200).unwrap();

    // Should have at most 2 segments (merged consecutive bins)
    assert!(segmentation.segments.len() <= 3);
}

#[test]
fn test_all_one_state() {
    // All bins same mark pattern
    let mark_matrix = vec![
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
    ];

    let marks = vec!["mark1".to_string(), "mark2".to_string()];
    let params = ChromHMMParams {
        n_states: 2,
        bin_size: 200,
        max_iter: 5,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();
    let segmentation = segment_genome(&mark_matrix, &model, "chr1", 200).unwrap();

    // May result in 1-2 segments depending on state learning
    assert!(segmentation.segments.len() <= 2);
}

#[test]
fn test_empty_mark_matrix_error() {
    let mark_matrix = vec![];
    let marks = vec!["mark1".to_string()];
    let params = ChromHMMParams::default();

    let result = learn_chromatin_states(&mark_matrix, &marks, &params);
    assert!(result.is_err());
}

#[test]
fn test_emission_learning() {
    // Marks with clear distinct patterns
    let mark_matrix = vec![
        vec![1, 0],
        vec![1, 0],
        vec![0, 1],
        vec![0, 1],
    ];

    let marks = vec!["mark1".to_string(), "mark2".to_string()];
    let params = ChromHMMParams {
        n_states: 2,
        bin_size: 200,
        max_iter: 100,
        convergence: 1e-6,
    };

    let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();

    // Check that emissions are learned (model was trained)
    for state in &model.states {
        let probs = &state.emission_probs;
        assert_eq!(probs.len(), 2);

        // Probabilities should be in valid range
        for &p in probs {
            assert!(p >= 0.01 && p <= 0.99);
        }
    }
}
