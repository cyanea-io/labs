//! Ranking methods for numeric data.
//!
//! Provides [`rank`] which assigns ranks to data values using various
//! tie-breaking strategies via [`RankMethod`].

/// Strategy for handling tied values when ranking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RankMethod {
    /// Tied values receive the average of their would-be ranks (default for
    /// most statistical methods).
    Average,
    /// Tied values receive the minimum rank in the group.
    Min,
    /// Tied values receive the maximum rank in the group.
    Max,
    /// Tied values receive sequential ranks (first occurrence gets lower rank).
    Ordinal,
    /// Like `Min`, but ranks always increase by 1 (no gaps).
    Dense,
}

/// Assign ranks to `data` using the given [`RankMethod`].
///
/// Returns a `Vec<f64>` of the same length as `data`, where each element is
/// the rank of the corresponding input value.
///
/// Empty input produces empty output.
pub fn rank(data: &[f64], method: RankMethod) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }

    // Build (value, original_index) and sort by value.
    let mut indexed: Vec<(f64, usize)> = data.iter().copied().enumerate().map(|(i, v)| (v, i)).collect();
    indexed.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut ranks = vec![0.0; n];

    if method == RankMethod::Ordinal {
        for (rank_minus_1, &(_, orig_idx)) in indexed.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_1 + 1) as f64;
        }
    } else {
        let mut dense_rank = 0.0;
        let mut i = 0;
        while i < n {
            // Find the end of the tie group.
            let mut j = i + 1;
            while j < n && indexed[j].0.total_cmp(&indexed[i].0).is_eq() {
                j += 1;
            }
            dense_rank += 1.0;

            // Ranks in the group are (i+1)..=(j) (1-based).
            let group_len = j - i;
            let rank_val = match method {
                RankMethod::Average => {
                    let sum: f64 = (i + 1..=j).map(|r| r as f64).sum();
                    sum / group_len as f64
                }
                RankMethod::Min => (i + 1) as f64,
                RankMethod::Max => j as f64,
                RankMethod::Dense => dense_rank,
                RankMethod::Ordinal => unreachable!(),
            };

            for k in i..j {
                ranks[indexed[k].1] = rank_val;
            }

            i = j;
        }
    }

    ranks
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_average_no_ties() {
        let data = [3.0, 1.0, 2.0];
        assert_eq!(rank(&data, RankMethod::Average), vec![3.0, 1.0, 2.0]);
    }

    #[test]
    fn rank_average_with_ties() {
        let data = [3.0, 1.0, 2.0, 2.0];
        // sorted: 1(1), 2(2), 2(3), 3(4) → ties at 2 get (2+3)/2 = 2.5
        assert_eq!(rank(&data, RankMethod::Average), vec![4.0, 1.0, 2.5, 2.5]);
    }

    #[test]
    fn rank_min_with_ties() {
        let data = [3.0, 1.0, 2.0, 2.0];
        assert_eq!(rank(&data, RankMethod::Min), vec![4.0, 1.0, 2.0, 2.0]);
    }

    #[test]
    fn rank_max_with_ties() {
        let data = [3.0, 1.0, 2.0, 2.0];
        assert_eq!(rank(&data, RankMethod::Max), vec![4.0, 1.0, 3.0, 3.0]);
    }

    #[test]
    fn rank_dense_with_ties() {
        let data = [3.0, 1.0, 2.0, 2.0];
        assert_eq!(rank(&data, RankMethod::Dense), vec![3.0, 1.0, 2.0, 2.0]);
    }

    #[test]
    fn rank_ordinal() {
        let data = [3.0, 1.0, 2.0, 2.0];
        // Ties broken by original position in sorted order.
        let r = rank(&data, RankMethod::Ordinal);
        assert_eq!(r[1], 1.0); // 1.0 is smallest
        assert_eq!(r[0], 4.0); // 3.0 is largest
        // The two 2.0s get ranks 2 and 3 (order depends on stable sort).
        assert!(r[2] == 2.0 || r[2] == 3.0);
        assert!(r[3] == 2.0 || r[3] == 3.0);
        assert!((r[2] - r[3]).abs() > 0.5); // they must differ
    }

    #[test]
    fn rank_all_equal() {
        let data = [5.0, 5.0, 5.0];
        assert_eq!(rank(&data, RankMethod::Average), vec![2.0, 2.0, 2.0]);
        assert_eq!(rank(&data, RankMethod::Dense), vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn rank_empty() {
        let empty: Vec<f64> = vec![];
        assert_eq!(rank(&empty, RankMethod::Average), Vec::<f64>::new());
    }
}
