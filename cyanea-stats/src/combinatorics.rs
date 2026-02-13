//! Combinatorics utilities.
//!
//! Provides exact and log-space variants of factorials, binomial coefficients,
//! permutations, and multinomial coefficients, plus a combinations iterator.

use crate::distribution::ln_gamma;
use cyanea_core::{CyaneaError, Result};

/// Exact factorial. Returns `None` if `n > 20` (overflow for u64).
pub fn factorial(n: u64) -> Option<u64> {
    if n > 20 {
        return None;
    }
    let mut result = 1u64;
    for i in 2..=n {
        result = result.checked_mul(i)?;
    }
    Some(result)
}

/// Log-factorial via `ln(Γ(n + 1))`.
pub fn ln_factorial(n: u64) -> f64 {
    ln_gamma(n as f64 + 1.0)
}

/// Exact binomial coefficient C(n, k). Returns `None` on overflow.
pub fn binomial(n: u64, k: u64) -> Option<u64> {
    if k > n {
        return Some(0);
    }
    // Use the smaller of k and n-k for efficiency
    let k = k.min(n - k);
    let mut result = 1u64;
    for i in 0..k {
        result = result.checked_mul(n - i)?;
        result /= i + 1;
    }
    Some(result)
}

/// Log-space binomial coefficient ln(C(n, k)).
///
/// # Errors
///
/// Returns an error if `k > n`.
pub fn ln_binomial(n: u64, k: u64) -> Result<f64> {
    if k > n {
        return Err(CyaneaError::InvalidInput(
            "ln_binomial: k must be <= n".into(),
        ));
    }
    Ok(ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0))
}

/// Exact permutations P(n, k) = n! / (n-k)!. Returns `None` on overflow.
pub fn permutations(n: u64, k: u64) -> Option<u64> {
    if k > n {
        return Some(0);
    }
    let mut result = 1u64;
    for i in 0..k {
        result = result.checked_mul(n - i)?;
    }
    Some(result)
}

/// Log-space permutations ln(P(n, k)).
///
/// # Errors
///
/// Returns an error if `k > n`.
pub fn ln_permutations(n: u64, k: u64) -> Result<f64> {
    if k > n {
        return Err(CyaneaError::InvalidInput(
            "ln_permutations: k must be <= n".into(),
        ));
    }
    Ok(ln_gamma(n as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0))
}

/// Exact multinomial coefficient n! / (c₁! · c₂! · ... · cₖ!).
///
/// Returns `None` on overflow. The counts must sum to `n`.
pub fn multinomial(n: u64, counts: &[u64]) -> Option<u64> {
    let sum: u64 = counts.iter().sum();
    if sum != n {
        return None;
    }
    let mut result = 1u64;
    let mut remaining = n;
    for &c in counts {
        result = result.checked_mul(binomial(remaining, c)?)?;
        remaining -= c;
    }
    Some(result)
}

/// Log-space multinomial coefficient.
///
/// # Errors
///
/// Returns an error if counts don't sum to `n`.
pub fn ln_multinomial(n: u64, counts: &[u64]) -> Result<f64> {
    let sum: u64 = counts.iter().sum();
    if sum != n {
        return Err(CyaneaError::InvalidInput(
            "ln_multinomial: counts must sum to n".into(),
        ));
    }
    let mut result = ln_gamma(n as f64 + 1.0);
    for &c in counts {
        result -= ln_gamma(c as f64 + 1.0);
    }
    Ok(result)
}

/// Iterator over all k-element combinations of indices `[0, n)`.
///
/// Yields combinations in lexicographic order. Each combination is a
/// `Vec<usize>` of length `k` with strictly increasing indices.
///
/// # Example
///
/// ```
/// use cyanea_stats::combinatorics::combinations;
///
/// let combos: Vec<Vec<usize>> = combinations(4, 2).collect();
/// assert_eq!(combos.len(), 6); // C(4, 2) = 6
/// assert_eq!(combos[0], vec![0, 1]);
/// assert_eq!(combos[5], vec![2, 3]);
/// ```
pub fn combinations(n: usize, k: usize) -> Combinations {
    let first = if k == 0 || k > n {
        None
    } else {
        Some((0..k).collect())
    };
    Combinations { n, k, current: first }
}

/// Iterator over k-element combinations of `[0, n)`.
#[derive(Debug, Clone)]
pub struct Combinations {
    n: usize,
    k: usize,
    current: Option<Vec<usize>>,
}

impl Iterator for Combinations {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current.clone()?;

        // Advance to next combination
        let mut next = result.clone();
        let mut i = self.k;
        while i > 0 {
            i -= 1;
            next[i] += 1;
            if next[i] <= self.n - self.k + i {
                // Fill remaining positions
                for j in (i + 1)..self.k {
                    next[j] = next[j - 1] + 1;
                }
                self.current = Some(next);
                return Some(result);
            }
        }

        self.current = None;
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorial_small() {
        assert_eq!(factorial(0), Some(1));
        assert_eq!(factorial(1), Some(1));
        assert_eq!(factorial(5), Some(120));
        assert_eq!(factorial(10), Some(3_628_800));
        assert_eq!(factorial(20), Some(2_432_902_008_176_640_000));
    }

    #[test]
    fn factorial_overflow() {
        assert_eq!(factorial(21), None);
    }

    #[test]
    fn ln_factorial_matches_exact() {
        for n in 0..=20 {
            let exact = factorial(n).unwrap() as f64;
            let ln_val = ln_factorial(n);
            assert!((ln_val - exact.ln()).abs() < 1e-8 || (n == 0 && ln_val.abs() < 1e-10),
                "ln_factorial({}) = {} but expected {}", n, ln_val, exact.ln());
        }
    }

    #[test]
    fn binomial_known_values() {
        assert_eq!(binomial(5, 0), Some(1));
        assert_eq!(binomial(5, 1), Some(5));
        assert_eq!(binomial(5, 2), Some(10));
        assert_eq!(binomial(5, 3), Some(10));
        assert_eq!(binomial(5, 5), Some(1));
        assert_eq!(binomial(10, 3), Some(120));
        assert_eq!(binomial(20, 10), Some(184_756));
    }

    #[test]
    fn binomial_k_greater_than_n() {
        assert_eq!(binomial(3, 5), Some(0));
    }

    #[test]
    fn ln_binomial_accuracy() {
        let ln_val = ln_binomial(10, 3).unwrap();
        let expected = (120.0_f64).ln();
        assert!((ln_val - expected).abs() < 1e-8);
    }

    #[test]
    fn ln_binomial_invalid() {
        assert!(ln_binomial(3, 5).is_err());
    }

    #[test]
    fn permutations_known() {
        assert_eq!(permutations(5, 3), Some(60)); // 5*4*3
        assert_eq!(permutations(5, 0), Some(1));
        assert_eq!(permutations(5, 5), Some(120));
    }

    #[test]
    fn permutations_k_greater_than_n() {
        assert_eq!(permutations(3, 5), Some(0));
    }

    #[test]
    fn ln_permutations_accuracy() {
        let ln_val = ln_permutations(5, 3).unwrap();
        let expected = (60.0_f64).ln();
        assert!((ln_val - expected).abs() < 1e-8);
    }

    #[test]
    fn multinomial_known() {
        // 4! / (2! * 1! * 1!) = 12
        assert_eq!(multinomial(4, &[2, 1, 1]), Some(12));
        // 6! / (3! * 2! * 1!) = 60
        assert_eq!(multinomial(6, &[3, 2, 1]), Some(60));
    }

    #[test]
    fn multinomial_bad_sum() {
        assert_eq!(multinomial(5, &[2, 1]), None);
    }

    #[test]
    fn ln_multinomial_accuracy() {
        let ln_val = ln_multinomial(4, &[2, 1, 1]).unwrap();
        let expected = (12.0_f64).ln();
        assert!((ln_val - expected).abs() < 1e-8);
    }

    #[test]
    fn ln_multinomial_invalid() {
        assert!(ln_multinomial(5, &[2, 1]).is_err());
    }

    #[test]
    fn combinations_count() {
        let combos: Vec<Vec<usize>> = combinations(5, 2).collect();
        assert_eq!(combos.len(), 10); // C(5,2) = 10
    }

    #[test]
    fn combinations_values() {
        let combos: Vec<Vec<usize>> = combinations(4, 2).collect();
        assert_eq!(combos[0], vec![0, 1]);
        assert_eq!(combos[1], vec![0, 2]);
        assert_eq!(combos[2], vec![0, 3]);
        assert_eq!(combos[3], vec![1, 2]);
        assert_eq!(combos[4], vec![1, 3]);
        assert_eq!(combos[5], vec![2, 3]);
    }

    #[test]
    fn combinations_k_zero() {
        let combos: Vec<Vec<usize>> = combinations(5, 0).collect();
        assert!(combos.is_empty());
    }

    #[test]
    fn combinations_k_equals_n() {
        let combos: Vec<Vec<usize>> = combinations(3, 3).collect();
        assert_eq!(combos.len(), 1);
        assert_eq!(combos[0], vec![0, 1, 2]);
    }

    #[test]
    fn combinations_k_greater_than_n() {
        let combos: Vec<Vec<usize>> = combinations(2, 5).collect();
        assert!(combos.is_empty());
    }
}
