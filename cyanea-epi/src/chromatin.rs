//! Chromatin state analysis and HMM-based segmentation.
//!
//! This module provides simplified ChromHMM-like learning and genomic state segmentation.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::Result;
use crate::EpiError;

/// A chromatin state in the model.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChromatinState {
    /// State ID (0-based).
    pub id: usize,
    /// State name.
    pub name: String,
    /// Emission probabilities for marks (mark index -> probability).
    pub emission_probs: Vec<f64>,
    /// Color for visualization (hex string).
    pub color: String,
}

/// Learned ChromHMM model.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChromHMMModel {
    /// Chromatin states.
    pub states: Vec<ChromatinState>,
    /// Transition matrix: [from_state][to_state].
    pub transition_matrix: Vec<Vec<f64>>,
    /// Histone marks used in the model.
    pub marks: Vec<String>,
}

/// Parameters for ChromHMM learning.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChromHMMParams {
    /// Number of states (default: 15).
    pub n_states: usize,
    /// Bin size in bp (default: 200).
    pub bin_size: u64,
    /// Maximum EM iterations (default: 200).
    pub max_iter: usize,
    /// Convergence threshold (default: 1e-6).
    pub convergence: f64,
}

impl Default for ChromHMMParams {
    fn default() -> Self {
        Self {
            n_states: 15,
            bin_size: 200,
            max_iter: 200,
            convergence: 1e-6,
        }
    }
}

/// Genomic segmentation into chromatin states.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChromatinSegmentation {
    /// State assignments per bin.
    pub state_assignments: Vec<usize>,
    /// State names for reference.
    pub state_names: Vec<String>,
    /// Genome locations (start, end, chrom).
    pub segments: Vec<(u64, u64, String)>,
}

/// Learn chromatin states from mark presence matrix using EM algorithm.
///
/// Input: binary presence/absence matrix of marks at genomic bins.
/// Each row is a bin, each column is a mark.
/// Learns transition probabilities and emission probabilities via EM.
pub fn learn_chromatin_states(
    mark_matrix: &[Vec<u32>],
    marks: &[String],
    params: &ChromHMMParams,
) -> Result<ChromHMMModel> {
    if mark_matrix.is_empty() {
        return Err(EpiError::InsufficientData(
            "Mark matrix cannot be empty".to_string(),
        ));
    }

    let n_bins = mark_matrix.len();
    let n_marks = marks.len();
    let n_states = params.n_states;

    // Initialize transition matrix uniformly
    let transition = vec![vec![1.0 / (n_states as f64); n_states]; n_states];

    // Initialize emission probabilities randomly
    let mut emission = vec![vec![0.5; n_marks]; n_states];

    // EM iterations
    for _iter in 0..params.max_iter {
        // E-step: compute forward-backward probabilities
        let forward = forward_algorithm(mark_matrix, &emission, &transition, n_states)?;
        let backward = backward_algorithm(mark_matrix, &emission, &transition, n_states)?;

        // M-step: update parameters
        let old_emission = emission.clone();

        // Update emission probabilities
        for s in 0..n_states {
            for m in 0..n_marks {
                let mut numerator = 0.0;
                let mut denominator = 0.0;

                for b in 0..n_bins {
                    let gamma = if !forward.is_empty() && !backward.is_empty() && !forward[b].is_empty() {
                        let norm = forward[b].iter().sum::<f64>().max(1e-10);
                        forward[b][s] / norm
                    } else {
                        0.0
                    };

                    numerator += gamma * (mark_matrix[b][m] as f64);
                    denominator += gamma;
                }

                emission[s][m] = (numerator + 0.01) / (denominator + 0.02);
                emission[s][m] = emission[s][m].clamp(0.01, 0.99);
            }
        }

        // Check convergence
        let diff: f64 = emission
            .iter()
            .zip(old_emission.iter())
            .flat_map(|(e1, e2)| e1.iter().zip(e2.iter()))
            .map(|(a, b)| (a - b).abs())
            .sum();

        if diff < params.convergence {
            break;
        }
    }

    // Create states
    let colors = generate_colors(n_states);
    let states = (0..n_states)
        .map(|i| ChromatinState {
            id: i,
            name: format!("S{}", i),
            emission_probs: emission[i].clone(),
            color: colors[i].clone(),
        })
        .collect();

    Ok(ChromHMMModel {
        states,
        transition_matrix: transition,
        marks: marks.to_vec(),
    })
}

/// Generate distinct colors for state visualization.
fn generate_colors(n_states: usize) -> Vec<String> {
    let hues = [
        0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 15, 45, 75, 105, 135, 165, 195,
        225, 255, 285, 315, 345,
    ];

    (0..n_states)
        .map(|i| {
            let hue = hues[i % hues.len()];
            format!("hsl({},70%,50%)", hue)
        })
        .collect()
}

/// Forward algorithm for HMM (simplified; no explicit initial state).
fn forward_algorithm(
    mark_matrix: &[Vec<u32>],
    emission: &[Vec<f64>],
    transition: &[Vec<f64>],
    n_states: usize,
) -> Result<Vec<Vec<f64>>> {
    let n_bins = mark_matrix.len();
    let n_marks = emission[0].len();

    let mut forward = vec![vec![0.0; n_states]; n_bins];

    // Initialize first bin
    for s in 0..n_states {
        let mut prob = 1.0 / (n_states as f64);
        for m in 0..n_marks {
            let mark_val = mark_matrix[0][m] as f64;
            prob *= if mark_val > 0.0 {
                emission[s][m]
            } else {
                1.0 - emission[s][m]
            };
        }
        forward[0][s] = prob;
    }

    // Iterate through bins
    for b in 1..n_bins {
        for s in 0..n_states {
            let mut sum = 0.0;
            for prev_s in 0..n_states {
                sum += forward[b - 1][prev_s] * transition[prev_s][s];
            }

            let mut obs_prob = 1.0;
            for m in 0..n_marks {
                let mark_val = mark_matrix[b][m] as f64;
                obs_prob *= if mark_val > 0.0 {
                    emission[s][m]
                } else {
                    1.0 - emission[s][m]
                };
            }

            forward[b][s] = sum * obs_prob;
        }
    }

    Ok(forward)
}

/// Backward algorithm for HMM.
fn backward_algorithm(
    mark_matrix: &[Vec<u32>],
    emission: &[Vec<f64>],
    transition: &[Vec<f64>],
    n_states: usize,
) -> Result<Vec<Vec<f64>>> {
    let n_bins = mark_matrix.len();
    let n_marks = emission[0].len();

    let mut backward = vec![vec![1.0; n_states]; n_bins];

    // Backward from last bin
    for b in (0..n_bins - 1).rev() {
        for s in 0..n_states {
            let mut sum = 0.0;
            for next_s in 0..n_states {
                let mut obs_prob = 1.0;
                for m in 0..n_marks {
                    let mark_val = mark_matrix[b + 1][m] as f64;
                    obs_prob *= if mark_val > 0.0 {
                        emission[next_s][m]
                    } else {
                        1.0 - emission[next_s][m]
                    };
                }

                sum += transition[s][next_s] * obs_prob * backward[b + 1][next_s];
            }
            backward[b][s] = sum;
        }
    }

    Ok(backward)
}

/// Segment genome using Viterbi decoding.
pub fn segment_genome(
    mark_matrix: &[Vec<u32>],
    model: &ChromHMMModel,
    chrom: &str,
    bin_size: u64,
) -> Result<ChromatinSegmentation> {
    let n_bins = mark_matrix.len();
    let n_states = model.states.len();

    // Viterbi algorithm: find most likely state sequence
    let mut viterbi = vec![vec![0.0; n_states]; n_bins];
    let mut backpointer = vec![vec![0usize; n_states]; n_bins];

    // Initialize first bin
    for s in 0..n_states {
        let mut prob = (1.0 / (n_states as f64)).ln();
        for (m, &mark_val) in mark_matrix[0].iter().enumerate() {
            let emit_prob = model.states[s].emission_probs[m];
            prob += if mark_val > 0 {
                emit_prob.max(1e-10).ln()
            } else {
                (1.0 - emit_prob).max(1e-10).ln()
            };
        }
        viterbi[0][s] = prob;
    }

    // Fill in remaining bins
    for b in 1..n_bins {
        for s in 0..n_states {
            let mut best_prob = f64::NEG_INFINITY;
            let mut best_prev = 0;

            for prev_s in 0..n_states {
                let trans_prob = model.transition_matrix[prev_s][s].max(1e-10).ln();
                let path_prob = viterbi[b - 1][prev_s] + trans_prob;

                if path_prob > best_prob {
                    best_prob = path_prob;
                    best_prev = prev_s;
                }
            }

            // Add emission probability
            let mut emit_prob = 0.0;
            for (m, &mark_val) in mark_matrix[b].iter().enumerate() {
                let prob = model.states[s].emission_probs[m];
                emit_prob += if mark_val > 0 {
                    prob.max(1e-10).ln()
                } else {
                    (1.0 - prob).max(1e-10).ln()
                };
            }

            viterbi[b][s] = best_prob + emit_prob;
            backpointer[b][s] = best_prev;
        }
    }

    // Traceback
    let mut path = vec![0; n_bins];
    path[n_bins - 1] = viterbi[n_bins - 1]
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    for b in (1..n_bins).rev() {
        path[b - 1] = backpointer[b][path[b]];
    }

    // Build segments (merge adjacent same states)
    let mut segments = Vec::new();
    let state_names: Vec<_> = model.states.iter().map(|s| s.name.clone()).collect();

    let mut current_state = path[0];
    let mut start_bin = 0;

    for b in 1..n_bins {
        if path[b] != current_state {
            let start_bp = start_bin as u64 * bin_size;
            let end_bp = b as u64 * bin_size;
            segments.push((start_bp, end_bp, chrom.to_string()));
            current_state = path[b];
            start_bin = b;
        }
    }

    // Final segment
    let start_bp = start_bin as u64 * bin_size;
    let end_bp = n_bins as u64 * bin_size;
    segments.push((start_bp, end_bp, chrom.to_string()));

    Ok(ChromatinSegmentation {
        state_assignments: path,
        state_names,
        segments,
    })
}

/// Compute enrichment of states at genomic annotations.
pub fn state_enrichment(
    segmentation: &ChromatinSegmentation,
    annotations: &[(u64, u64)],
) -> Vec<(usize, f64)> {
    let mut state_overlap_counts = vec![0u64; segmentation.state_names.len()];
    let mut total_overlap = 0u64;

    for (ann_start, ann_end) in annotations {
        for (seg_start, seg_end, _) in &segmentation.segments {
            let overlap_start = ann_start.max(seg_start);
            let overlap_end = ann_end.min(seg_end);

            if overlap_start < overlap_end {
                let overlap = overlap_end - overlap_start;
                total_overlap += overlap;

                // Find state for this segment
                let bin_idx = (seg_start / 200).min(segmentation.state_assignments.len() as u64 - 1) as usize;
                if bin_idx < segmentation.state_assignments.len() {
                    let state = segmentation.state_assignments[bin_idx];
                    state_overlap_counts[state] += overlap;
                }
            }
        }
    }

    let mut enrichments = Vec::new();
    for (state_id, count) in state_overlap_counts.into_iter().enumerate() {
        let enrichment = if total_overlap == 0 {
            0.0
        } else {
            (count as f64) / (total_overlap as f64)
        };
        enrichments.push((state_id, enrichment));
    }

    enrichments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chromatin_state_creation() {
        let state = ChromatinState {
            id: 0,
            name: "Active".to_string(),
            emission_probs: vec![0.9, 0.1],
            color: "#FF0000".to_string(),
        };

        assert_eq!(state.name, "Active");
        assert_eq!(state.id, 0);
    }

    #[test]
    fn test_learn_chromatin_states() {
        let mark_matrix = vec![
            vec![1, 0],
            vec![1, 1],
            vec![0, 1],
        ];

        let marks = vec!["H3K4me3".to_string(), "H3K27me3".to_string()];
        let params = ChromHMMParams {
            n_states: 3,
            bin_size: 200,
            max_iter: 10,
            convergence: 1e-6,
        };

        let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();
        assert_eq!(model.states.len(), 3);
        assert_eq!(model.marks.len(), 2);
    }

    #[test]
    fn test_segment_genome() {
        let mark_matrix = vec![
            vec![1, 0],
            vec![1, 1],
            vec![0, 1],
            vec![0, 0],
        ];

        let marks = vec!["H3K4me3".to_string(), "H3K27me3".to_string()];
        let params = ChromHMMParams::default();

        let model = learn_chromatin_states(&mark_matrix, &marks, &params).unwrap();
        let segmentation = segment_genome(&mark_matrix, &model, "chr1", 200).unwrap();

        assert_eq!(segmentation.state_assignments.len(), 4);
        assert!(!segmentation.segments.is_empty());
    }

    #[test]
    fn test_state_enrichment() {
        let segmentation = ChromatinSegmentation {
            state_assignments: vec![0, 0, 1, 1],
            state_names: vec!["Active".to_string(), "Repressed".to_string()],
            segments: vec![
                (0, 200, "chr1".to_string()),
                (200, 400, "chr1".to_string()),
            ],
        };

        let annotations = vec![(0, 300)];
        let enrichments = state_enrichment(&segmentation, &annotations);

        assert!(!enrichments.is_empty());
    }
}
