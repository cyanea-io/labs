//! Substitution models for maximum likelihood phylogenetics.
//!
//! Provides the JC69 (Jukes-Cantor 1969) nucleotide substitution model,
//! which assumes equal base frequencies and a single substitution rate.

/// Map a nucleotide byte to an index (A=0, C=1, G=2, T=3).
///
/// Accepts both upper and lower case. Returns `None` for non-standard bases.
pub fn nucleotide_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// JC69 transition probability matrix for a given branch length `t`.
///
/// Under JC69, the probability of observing nucleotide `j` at time `t`
/// given nucleotide `i` at time 0 is:
///
/// - P(same) = 1/4 + 3/4 * e^{-4t/3}
/// - P(diff) = 1/4 - 1/4 * e^{-4t/3}
///
/// The matrix is 4x4, indexed by `nucleotide_index` (A=0, C=1, G=2, T=3).
pub fn jc69_probability(t: f64) -> [[f64; 4]; 4] {
    let e = (-4.0 * t / 3.0).exp();
    let p_same = 0.25 + 0.75 * e;
    let p_diff = 0.25 - 0.25 * e;

    [
        [p_same, p_diff, p_diff, p_diff],
        [p_diff, p_same, p_diff, p_diff],
        [p_diff, p_diff, p_same, p_diff],
        [p_diff, p_diff, p_diff, p_same],
    ]
}

/// Number of nucleotide states in DNA models.
pub const NUM_STATES: usize = 4;

/// Equilibrium frequency for JC69 (uniform: 1/4 for each base).
pub const JC69_FREQ: f64 = 0.25;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nucleotide_index_standard() {
        assert_eq!(nucleotide_index(b'A'), Some(0));
        assert_eq!(nucleotide_index(b'C'), Some(1));
        assert_eq!(nucleotide_index(b'G'), Some(2));
        assert_eq!(nucleotide_index(b'T'), Some(3));
    }

    #[test]
    fn nucleotide_index_lowercase() {
        assert_eq!(nucleotide_index(b'a'), Some(0));
        assert_eq!(nucleotide_index(b'c'), Some(1));
        assert_eq!(nucleotide_index(b'g'), Some(2));
        assert_eq!(nucleotide_index(b't'), Some(3));
    }

    #[test]
    fn nucleotide_index_uracil() {
        assert_eq!(nucleotide_index(b'U'), Some(3));
        assert_eq!(nucleotide_index(b'u'), Some(3));
    }

    #[test]
    fn nucleotide_index_invalid() {
        assert_eq!(nucleotide_index(b'N'), None);
        assert_eq!(nucleotide_index(b'-'), None);
        assert_eq!(nucleotide_index(b'X'), None);
    }

    #[test]
    fn jc69_t_zero_is_identity() {
        let p = jc69_probability(0.0);
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    assert!(
                        (p[i][j] - 1.0).abs() < 1e-12,
                        "P[{}][{}] = {}, expected 1.0",
                        i,
                        j,
                        p[i][j]
                    );
                } else {
                    assert!(
                        p[i][j].abs() < 1e-12,
                        "P[{}][{}] = {}, expected 0.0",
                        i,
                        j,
                        p[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn jc69_large_t_approaches_uniform() {
        let p = jc69_probability(1000.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p[i][j] - 0.25).abs() < 1e-10,
                    "P[{}][{}] = {}, expected ~0.25",
                    i,
                    j,
                    p[i][j]
                );
            }
        }
    }

    #[test]
    fn jc69_rows_sum_to_one() {
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 100.0] {
            let p = jc69_probability(t);
            for i in 0..4 {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "row {} sum = {} at t = {}",
                    i,
                    row_sum,
                    t
                );
            }
        }
    }

    #[test]
    fn jc69_symmetric() {
        let p = jc69_probability(0.3);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p[i][j] - p[j][i]).abs() < 1e-12,
                    "P[{}][{}] != P[{}][{}]",
                    i,
                    j,
                    j,
                    i
                );
            }
        }
    }

    #[test]
    fn jc69_probabilities_non_negative() {
        for &t in &[0.0, 0.001, 0.01, 0.1, 1.0, 10.0] {
            let p = jc69_probability(t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        p[i][j] >= 0.0,
                        "P[{}][{}] = {} < 0 at t = {}",
                        i,
                        j,
                        p[i][j],
                        t
                    );
                }
            }
        }
    }

    #[test]
    fn jc69_known_value() {
        // At t = 0.3: e^{-4*0.3/3} = e^{-0.4}
        let e = (-0.4_f64).exp();
        let p_same = 0.25 + 0.75 * e;
        let p_diff = 0.25 - 0.25 * e;
        let p = jc69_probability(0.3);
        assert!((p[0][0] - p_same).abs() < 1e-12);
        assert!((p[0][1] - p_diff).abs() < 1e-12);
    }
}
