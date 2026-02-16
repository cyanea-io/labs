//! Statistical methods for the Cyanea bioinformatics ecosystem.
//!
//! `cyanea-stats` provides pure-Rust implementations of common statistical
//! methods used in bioinformatics workflows:
//!
//! - **Descriptive statistics** — [`descriptive::describe`], [`descriptive::mean`],
//!   [`descriptive::variance`], quantiles, IQR, MAD
//! - **Ranking** — [`rank::rank`] with multiple tie-breaking strategies
//! - **Correlation** — [`correlation::pearson`], [`correlation::spearman`],
//!   [`correlation::CorrelationMatrix`]
//! - **Hypothesis testing** — [`testing::t_test_one_sample`],
//!   [`testing::t_test_two_sample`], [`testing::mann_whitney_u`]
//! - **Multiple testing correction** — [`correction::bonferroni`],
//!   [`correction::benjamini_hochberg`]
//! - **Distributions** — [`distribution::Normal`], [`distribution::Poisson`],
//!   plus numerical helpers [`distribution::erf`], [`distribution::ln_gamma`],
//!   [`distribution::betai`]
//!
//! # Example
//!
//! ```
//! use cyanea_stats::descriptive::describe;
//!
//! let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
//! let stats = describe(&data).unwrap();
//! assert_eq!(stats.count, 8);
//! assert!((stats.mean - 5.0).abs() < 1e-10);
//! ```

pub mod bayesian;
pub mod combinatorics;
pub mod correction;
pub mod correlation;
pub mod descriptive;
pub mod distribution;
pub mod effect_size;
pub mod popgen;
pub mod rank;
pub mod reduction;
pub mod testing;

pub use bayesian::{Beta, Dirichlet, Gamma, NormalConjugate};
pub use combinatorics::{
    binomial, combinations, factorial, ln_binomial, ln_factorial, ln_multinomial, ln_permutations,
    multinomial, permutations, Combinations,
};
pub use correction::CorrectionMethod;
pub use correlation::CorrelationMatrix;
pub use descriptive::DescriptiveStats;
pub use distribution::{Binomial, ChiSquared, Distribution, FDistribution, Normal, Poisson};
pub use popgen::{
    AlleleFrequencies, DiversityStats, FstMethod, FstResult, HweResult, LdResult, TajimaD,
};
pub use rank::RankMethod;
pub use testing::TestResult;
