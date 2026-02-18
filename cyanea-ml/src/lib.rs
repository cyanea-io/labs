//! ML primitives for bioinformatics in the Cyanea ecosystem.
//!
//! This crate provides building blocks for machine-learning workflows on
//! biological data — distance metrics, sequence encoding, k-mer feature
//! extraction, normalization, clustering, and cluster evaluation.
//!
//! # Quick start
//!
//! ```
//! use cyanea_ml::kmer::KmerCounter;
//! use cyanea_ml::encoding::Alphabet;
//!
//! let counter = KmerCounter::new(3).unwrap();
//! let counts = counter.count_sequence(b"ACGTACGT");
//! assert!(counts.total() > 0);
//!
//! let freq = counts.to_frequency_vector(Alphabet::Dna);
//! assert_eq!(freq.len(), 64); // 4^3
//! ```

pub mod cluster;
pub mod cross_validation;
pub mod distance;
pub mod embedding;
pub mod encoding;
pub mod evaluate;
pub mod forest;
pub mod gbdt;
pub mod hmm;
pub mod inference;
pub mod kmer;
pub mod metrics;
pub mod normalize;
pub mod reduction;
pub mod tree;
pub mod umap;

pub use cluster::{
    dbscan, hierarchical, kmeans, DbscanConfig, DbscanResult, HierarchicalConfig,
    HierarchicalResult, KMeansConfig, KMeansResult, Linkage, MergeStep,
};
pub use distance::{
    compute_distance, cosine_distance, cosine_similarity, euclidean, hamming, manhattan,
    pairwise_distances, DistanceMatrix, DistanceMetric,
};
pub use encoding::{label_encode, one_hot_encode, Alphabet};
pub use evaluate::{silhouette_samples, silhouette_score};
pub use metrics::{
    accuracy, f1_macro, f1_score, f1_weighted, matthews_corrcoef, pr_auc, pr_curve, roc_auc,
    roc_curve, ConfusionMatrix, PrCurve, PrPoint, RocCurve, RocPoint,
};
pub use cross_validation::{
    cross_validate_kfold, cross_validate_loo, cross_validate_stratified, CvResult, FoldResult,
};
pub use kmer::{KmerCounter, KmerCounts};
pub use normalize::{
    l2_normalize, l2_normalize_columns, min_max, min_max_columns, z_score, z_score_columns,
};
pub use reduction::{pca, tsne, PcaConfig, PcaResult, TsneConfig, TsneResult};
pub use umap::{umap, UmapConfig, UmapInit, UmapResult};
pub use hmm::HmmModel;
pub use inference::{KnnConfig, KnnModel, LinearRegression};
pub use tree::DecisionTree;
pub use forest::{RandomForest, RandomForestConfig};
pub use gbdt::{GbdtConfig, GradientBoostedTrees};
