//! Log-space probability types for numerically stable computation.
//!
//! [`LogProb`] represents probabilities as natural logarithms and [`PhredProb`]
//! as Phred quality scores. Both prevent underflow in chains of small
//! probabilities and provide safe conversions between representations.

use crate::{CyaneaError, Result};

/// A probability stored as its natural logarithm: `ln(p)`.
///
/// All values are ≤ 0 (since 0 < p ≤ 1), with 0.0 representing certainty
/// (p = 1) and negative infinity representing impossibility (p = 0).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct LogProb(pub f64);

impl LogProb {
    /// Create a [`LogProb`] from a raw probability in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns an error if `p` is not in `(0, 1]`.
    pub fn from_prob(p: f64) -> Result<Self> {
        if p <= 0.0 || p > 1.0 {
            return Err(CyaneaError::InvalidInput(
                "LogProb::from_prob: p must be in (0, 1]".into(),
            ));
        }
        Ok(Self(p.ln()))
    }

    /// Convert back to a raw probability.
    pub fn to_prob(self) -> f64 {
        self.0.exp()
    }

    /// Log-sum-exp: compute `ln(exp(self) + exp(other))` without overflow.
    ///
    /// This is the log-space equivalent of addition in probability space.
    pub fn ln_add(self, other: Self) -> Self {
        if self.0 == f64::NEG_INFINITY {
            return other;
        }
        if other.0 == f64::NEG_INFINITY {
            return self;
        }
        let (max, min) = if self.0 >= other.0 {
            (self.0, other.0)
        } else {
            (other.0, self.0)
        };
        Self(max + (min - max).exp().ln_1p())
    }

    /// Multiply two probabilities in log-space (addition of log values).
    pub fn ln_mul(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }

    /// Certain event: `ln(1) = 0`.
    pub const fn certain() -> Self {
        Self(0.0)
    }

    /// Impossible event: `ln(0) = -∞`.
    pub const fn impossible() -> Self {
        Self(f64::NEG_INFINITY)
    }
}

/// A probability stored as a Phred quality score: `Q = -10 · log₁₀(p)`.
///
/// Higher values indicate lower error probability. Phred 30 ≈ 0.001.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct PhredProb(pub f64);

/// Conversion constant: `10 / ln(10)`.
const PHRED_SCALE: f64 = 10.0 / core::f64::consts::LN_10;

impl PhredProb {
    /// Create a [`PhredProb`] from a Phred score (must be ≥ 0).
    ///
    /// # Errors
    ///
    /// Returns an error if `q` is negative.
    pub fn from_phred(q: f64) -> Result<Self> {
        if q < 0.0 {
            return Err(CyaneaError::InvalidInput(
                "PhredProb::from_phred: q must be non-negative".into(),
            ));
        }
        Ok(Self(q))
    }

    /// Create a [`PhredProb`] from a raw probability in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns an error if `p` is not in `(0, 1]`.
    pub fn from_prob(p: f64) -> Result<Self> {
        if p <= 0.0 || p > 1.0 {
            return Err(CyaneaError::InvalidInput(
                "PhredProb::from_prob: p must be in (0, 1]".into(),
            ));
        }
        Ok(Self(-10.0 * p.log10()))
    }

    /// The Phred quality score.
    pub fn to_phred(self) -> f64 {
        self.0
    }

    /// Convert to a raw probability: `p = 10^(-Q/10)`.
    pub fn to_prob(self) -> f64 {
        10.0_f64.powf(-self.0 / 10.0)
    }
}

impl From<PhredProb> for LogProb {
    /// Convert Phred → LogProb: `ln(p) = -Q · ln(10) / 10`.
    fn from(phred: PhredProb) -> Self {
        Self(-phred.0 / PHRED_SCALE)
    }
}

impl From<LogProb> for PhredProb {
    /// Convert LogProb → Phred: `Q = -10 · log₁₀(p) = -ln(p) · 10 / ln(10)`.
    fn from(lp: LogProb) -> Self {
        Self(-lp.0 * PHRED_SCALE)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn logprob_from_prob_one() {
        let lp = LogProb::from_prob(1.0).unwrap();
        assert!((lp.0 - 0.0).abs() < TOL);
    }

    #[test]
    fn logprob_from_prob_half() {
        let lp = LogProb::from_prob(0.5).unwrap();
        assert!((lp.0 - 0.5_f64.ln()).abs() < TOL);
    }

    #[test]
    fn logprob_roundtrip() {
        let p = 0.001;
        let lp = LogProb::from_prob(p).unwrap();
        assert!((lp.to_prob() - p).abs() < TOL);
    }

    #[test]
    fn logprob_invalid() {
        assert!(LogProb::from_prob(0.0).is_err());
        assert!(LogProb::from_prob(-0.5).is_err());
        assert!(LogProb::from_prob(1.5).is_err());
    }

    #[test]
    fn logprob_certain_impossible() {
        assert_eq!(LogProb::certain().0, 0.0);
        assert_eq!(LogProb::certain().to_prob(), 1.0);
        assert_eq!(LogProb::impossible().0, f64::NEG_INFINITY);
        assert_eq!(LogProb::impossible().to_prob(), 0.0);
    }

    #[test]
    fn logprob_ln_mul() {
        let a = LogProb::from_prob(0.5).unwrap();
        let b = LogProb::from_prob(0.5).unwrap();
        let product = a.ln_mul(b);
        assert!((product.to_prob() - 0.25).abs() < TOL);
    }

    #[test]
    fn logprob_ln_add() {
        let a = LogProb::from_prob(0.3).unwrap();
        let b = LogProb::from_prob(0.2).unwrap();
        let sum = a.ln_add(b);
        assert!((sum.to_prob() - 0.5).abs() < TOL);
    }

    #[test]
    fn logprob_ln_add_identity() {
        let a = LogProb::from_prob(0.7).unwrap();
        let sum = a.ln_add(LogProb::impossible());
        assert!((sum.to_prob() - 0.7).abs() < TOL);

        let sum2 = LogProb::impossible().ln_add(a);
        assert!((sum2.to_prob() - 0.7).abs() < TOL);
    }

    #[test]
    fn phredprob_from_phred() {
        let q = PhredProb::from_phred(30.0).unwrap();
        assert!((q.to_prob() - 0.001).abs() < 1e-10);
    }

    #[test]
    fn phredprob_from_prob() {
        let q = PhredProb::from_prob(0.001).unwrap();
        assert!((q.to_phred() - 30.0).abs() < 1e-8);
    }

    #[test]
    fn phredprob_roundtrip() {
        let q = PhredProb::from_phred(20.0).unwrap();
        let p = q.to_prob();
        let q2 = PhredProb::from_prob(p).unwrap();
        assert!((q2.to_phred() - 20.0).abs() < 1e-8);
    }

    #[test]
    fn phredprob_invalid() {
        assert!(PhredProb::from_phred(-1.0).is_err());
        assert!(PhredProb::from_prob(0.0).is_err());
        assert!(PhredProb::from_prob(1.5).is_err());
    }

    #[test]
    fn convert_phred_to_logprob() {
        let phred = PhredProb::from_phred(30.0).unwrap();
        let lp: LogProb = phred.into();
        assert!((lp.to_prob() - 0.001).abs() < 1e-10);
    }

    #[test]
    fn convert_logprob_to_phred() {
        let lp = LogProb::from_prob(0.001).unwrap();
        let phred: PhredProb = lp.into();
        assert!((phred.to_phred() - 30.0).abs() < 1e-6);
    }

    #[test]
    fn phred_logprob_roundtrip_conversion() {
        let original = PhredProb::from_phred(25.0).unwrap();
        let lp: LogProb = original.into();
        let back: PhredProb = lp.into();
        assert!((back.to_phred() - 25.0).abs() < 1e-10);
    }
}
