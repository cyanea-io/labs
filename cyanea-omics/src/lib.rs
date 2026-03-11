//! Omics data structures for the Cyanea bioinformatics ecosystem.
//!
//! This crate provides core types for working with omics data:
//!
//! - **Genomic coordinates** — [`Strand`], [`GenomicPosition`], [`GenomicInterval`]
//! - **Interval collections** — [`IntervalSet`] with overlap queries
//! - **Expression matrices** — Dense [`ExpressionMatrix`] (features × samples)
//! - **Sparse matrices** — [`SparseMatrix`] in COO format
//! - **Variants** — VCF-style [`Variant`] representation
//! - **Gene annotations** — [`Gene`], [`Transcript`], [`Exon`] hierarchy
//!
//! # Quick start
//!
//! ```
//! use cyanea_omics::ExpressionMatrix;
//! use cyanea_core::Summarizable;
//!
//! let matrix = ExpressionMatrix::new(
//!     vec![vec![1.0, 2.0], vec![3.0, 4.0]],
//!     vec!["gene1".into(), "gene2".into()],
//!     vec!["sample_a".into(), "sample_b".into()],
//! ).unwrap();
//!
//! assert_eq!(matrix.shape(), (2, 2));
//! assert_eq!(matrix.get(0, 1), Some(2.0));
//! assert_eq!(matrix.summary(), "ExpressionMatrix: 2 features \u{00d7} 2 samples");
//! ```

pub mod acmg;
pub mod clinical;
pub mod genomic;
pub mod haplotype;
pub mod interval;
pub mod interval_tree;
pub mod coverage;
pub mod expr;
pub mod network;
pub mod otu;
pub mod sparse;
pub mod variant;
pub mod annotation;
pub mod variant_annotation;
pub mod single_cell;
pub mod spatial;
pub mod spatial_platforms;
pub mod spatial_segmentation;
pub mod spatial_domains;
pub mod spatial_cellchat;
pub mod spatial_deconvolution;
pub mod genome_arithmetic;
pub mod cnv;
pub mod liftover;
pub mod methylation;
pub mod pharmacogenomics;
pub mod microarray;
#[cfg(feature = "h5ad")]
pub mod h5ad;
#[cfg(feature = "zarr")]
pub mod zarr;
#[cfg(feature = "single-cell")]
pub mod sc_preprocess;
#[cfg(feature = "single-cell")]
pub mod sc_cluster;
#[cfg(feature = "single-cell")]
pub mod sc_trajectory;
#[cfg(feature = "single-cell")]
pub mod sc_markers;
#[cfg(feature = "single-cell")]
pub mod sc_integrate;
#[cfg(feature = "single-cell")]
pub mod sc_mtx;
#[cfg(feature = "single-cell")]
pub mod sc_velocity;
#[cfg(feature = "single-cell")]
pub mod sc_batch;

pub use cnv::{
    BafSegment, CbsConfig, CnvSegment, SvBreakpoint, SvType,
    baf_segmentation, circular_binary_segmentation, detect_sv_breakpoints, merge_cnv_segments,
};
pub use genomic::{GenomicInterval, GenomicPosition, Strand};
pub use interval::IntervalSet;
pub use interval_tree::{Interval, IntervalTree};
pub use coverage::RleCoverage;
pub use expr::ExpressionMatrix;
pub use sparse::SparseMatrix;
pub use variant::{Variant, VariantFilter, VariantType, Zygosity};
pub use annotation::{Exon, Gene, GeneType, Transcript};
pub use variant_annotation::{
    AnnotationConfig, Consequence, SpliceScore, VariantEffect,
    annotate_variant, annotate_variants, score_splice_disruption,
};
pub use single_cell::ColumnData;
pub use genome_arithmetic::{
    ClosestResult, GenomeInfo, JaccardStats, StrandMode,
    closest, complement, genome_info, intersect, intersect_report_a,
    jaccard, jaccard_stats, make_sliding_windows, make_windows,
    merge, subtract, union, windows_around,
};
pub use liftover::{ChainFile, LiftoverResult, liftover, liftover_batch, parse_chain};
pub use methylation::{
    CpgIsland, CpgSite, DmRegion, DmrConfig,
    bisulfite_convert, call_methylation, find_cpg_islands, find_dmrs,
};
pub use spatial::{
    CooccurrenceResult, GearysC, LrInteraction, SpatialAutocorrelation, SpatialGraph, SpatialPoint,
    cooccurrence, delaunay_neighbors, gearys_c, knn_spatial_neighbors, ligand_receptor_score,
    morans_i,
};
pub use spatial_platforms::{
    MerfishData, SlideseqData, VisiumData, VisiumScaleFactors,
    merfish_to_spatial_points, slideseq_to_spatial_points, visium_to_spatial_points,
};
pub use spatial_segmentation::{
    ExpansionParams, SegmentationResult, SegmentedCell,
    expansion_segmentation, voronoi_segmentation, watershed_grid,
};
pub use spatial_domains::{
    DomainParams, DomainResult, SpatialDomain, SpatiallyVariableGene,
    detect_domains, find_spatially_variable_genes, hmrf_smooth,
};
pub use spatial_cellchat::{
    CommParams, CommunicationResult, LrPair, PathwayCommunication,
    aggregate_pathways, analyze_communication, demo_lr_database,
};
pub use spatial_deconvolution::{
    CellTypeSignature, DeconvolutionResult, EnrichmentScore, SpotDeconvolution,
    nnls_deconvolve, score_enrichment,
};
pub use otu::OtuTable;
pub use network::{CentralityScores, Community, Graph};
pub use haplotype::{
    haplotype_blocks, haplotype_diversity, phase_em, Haplotype, HaplotypeBlock, PhasedGenotypes,
};
// Re-export ACMG/ClinVar types
pub use acmg::{
    auto_evidence, match_clinvar, parse_clinvar_tsv, AcmgClass, AcmgClassification,
    AcmgCriterion, AcmgEvidence, ClinVarAnnotation, EvidenceStrength,
};

// Re-export clinical genomics types
pub use clinical::{
    bethesda_markers, call_msi, compute_tmb, hla_compatibility, parse_hla_typing,
    HlaAllele, HlaTypingResult, MsiLocus, MsiResult, MsiStatus, TmbCategory, TmbResult,
};

// Re-export microarray analysis types
pub use microarray::{
    compute_beta, beta_to_m_value, diff_methylation, limma_diff_expr, m_value_to_beta,
    median_polish, quantile_normalize, rma_normalize, swan_normalize,
    DiffExprResult, DiffMethResult, InfiniumType, MethylationProbe,
};

// Re-export pharmacogenomics types
pub use pharmacogenomics::{
    activity_to_phenotype, call_star_alleles, demo_cyp2d6_database, lookup_drug_interactions,
    AlleleFunction, DrugGeneInteraction, MetabolizerPhenotype, PgxDatabase, StarAllele,
    StarAlleleCall,
};

#[cfg(feature = "h5ad")]
pub use h5ad::{read_h5ad, write_h5ad};
#[cfg(feature = "zarr")]
pub use zarr::{read_zarr, write_zarr};
