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
pub mod diffexpr;
pub mod distribution;
pub mod diversity;
pub mod enrichment;
pub mod effect_size;
pub mod multivariate;
pub mod normalization;
pub mod null_model;
pub mod ordination;
pub mod popgen;
pub mod rank;
pub mod reduction;
pub mod survival;
pub mod testing;

pub use bayesian::{Beta, Dirichlet, Gamma, NormalConjugate};
pub use combinatorics::{
    binomial, combinations, factorial, ln_binomial, ln_factorial, ln_multinomial, ln_permutations,
    multinomial, permutations, Combinations,
};
pub use correction::CorrectionMethod;
pub use correlation::CorrelationMatrix;
pub use descriptive::DescriptiveStats;
pub use diffexpr::{DeGeneResult, DeMethod, DeResults, VolcanoPoint};
pub use enrichment::{
    GeneSet, GoAnnotation, GoEnrichmentConfig, GoEnrichmentResult, GoNamespace, GoTerm, GseaResult,
    OraResult,
};
pub use distribution::{
    Binomial, ChiSquared, Distribution, FDistribution, NegativeBinomial, Normal, Poisson,
};
pub use popgen::{
    AlleleFrequencies, DiversityStats, FstMethod, FstResult, HweResult, LdResult, TajimaD,
};
pub use rank::RankMethod;
pub use survival::{CoxPhResult, KmResult, KmStep, LogRankResult};
pub use testing::TestResult;
pub use diversity::{
    alpha_diversity, alpha_rarefaction, bray_curtis, bray_curtis_matrix, chao1, hill_numbers,
    jaccard, jaccard_matrix, rarefaction_curve, shannon_index, simpson_index, weighted_jaccard,
    AlphaDiversity,
};
pub use multivariate::{
    amova, anosim, bioenv, mantel_test, permanova, AmovaResult, AnosimResult, BioenvResult,
    MantelResult, PermanovaResult,
};
pub use ordination::{
    cca, nmds, pcoa, procrustes, rda, ConstrainedOrdinationResult, NmdsConfig, NmdsResult,
    PcoaResult, ProcrustesResult,
};
