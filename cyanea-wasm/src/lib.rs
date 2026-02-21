//! WASM bindings and browser runtime for the Cyanea bioinformatics ecosystem.
//!
//! This crate provides in-memory, JSON-based wrappers around the Cyanea labs
//! crates, designed for environments where file I/O is unavailable (browsers,
//! sandboxed workers). Every public function accepts simple types (`&str`,
//! `f64`, `usize`) and returns a JSON `String`:
//!
//! - Success: `{"ok": <value>}`
//! - Failure: `{"error": "<message>"}`
//!
//! No `wasm-bindgen` dependency is included yet — `#[wasm_bindgen]` annotations
//! are a thin layer added when building for `wasm32`.
//!
//! # Modules
//!
//! - [`seq`] — In-memory FASTA parsing, GC content
//! - [`align`] — DNA and protein sequence alignment
//! - [`stats`] — Descriptive statistics, correlation, hypothesis testing
//! - [`ml`] — K-mer counting, distance metrics
//! - [`chem`] — SMILES parsing, molecular properties, fingerprints, similarity
//! - [`struct_bio`] — PDB parsing, secondary structure, RMSD
//! - [`phylo`] — Newick trees, evolutionary distances, UPGMA/NJ, RF distance
//!
//! # Example
//!
//! ```
//! let json = cyanea_wasm::align_dna("ACGT", "ACGT", "global");
//! let v: serde_json::Value = serde_json::from_str(&json).unwrap();
//! assert!(v["ok"]["score"].as_i64().unwrap() > 0);
//! ```

pub mod error;
pub mod seq;
pub mod align;
pub mod stats;
pub mod ml;
pub mod core_utils;
pub mod chem;
pub mod struct_bio;
pub mod phylo;
pub mod io;
pub mod omics;

/// Crate version (set from Cargo.toml at compile time).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ── Re-exports ───────────────────────────────────────────────────────────

// seq
pub use seq::{
    parse_fasta, parse_fasta_bytes, gc_content, gc_content_json,
    reverse_complement, transcribe, translate, validate, parse_fastq,
    parse_paired_fastq, parse_interleaved_fastq, trim_fastq, trim_paired_fastq,
    minhash_sketch, minhash_compare, JsMinHashSketch, JsMinHashComparison,
    rna_fold_nussinov, rna_fold_zuker, protein_props, simulate_reads,
    codon_usage, assembly_stats_json,
    JsRnaStructure, JsProteinProperties, JsSimulatedRead, JsCodonUsage, JsAssemblyStats,
};

// align
pub use align::{align_dna, align_dna_custom, align_protein, align_batch};
pub use align::{
    parse_cigar, validate_cigar, cigar_stats, cigar_to_alignment, alignment_to_cigar,
    generate_md_tag, merge_cigar, reverse_cigar, collapse_cigar, hard_clip_to_soft, split_cigar,
    progressive_msa, poa_consensus, align_banded,
    JsMsaResult, JsPoaConsensus,
};

// stats
pub use stats::{
    describe, pearson, spearman, t_test, t_test_two_sample,
    mann_whitney_u, bonferroni, benjamini_hochberg,
    kaplan_meier, log_rank_test, cox_ph,
    wright_fisher, permutation_test, bootstrap_ci,
    shannon_index, simpson_index, bray_curtis,
    fst_hudson, tajimas_d,
    JsDescriptiveStats, JsTestResult,
    JsKmStep, JsKmResult, JsLogRankResult, JsCoxPhResult,
    JsWrightFisherResult, JsTajimaD, JsFstResult,
};

// ml
pub use ml::{
    kmer_count, euclidean_distance, manhattan_distance, hamming_distance, cosine_similarity,
    umap, pca, tsne, kmeans,
    random_forest_classify, gbdt_regression, gbdt_classify,
    hmm_viterbi, hmm_likelihood,
    confusion_matrix, roc_curve, pr_curve,
    cross_validate_rf, feature_importance_variance,
    JsKmerCounts, JsUmapResult, JsPcaResult, JsTsneResult, JsKmeansResult,
    JsRandomForestResult, JsGbdtRegressionResult, JsGbdtClassifyResult,
    JsHmmViterbiResult, JsConfusionMatrix, JsRocCurve, JsRocPoint, JsPrCurve, JsPrPoint,
    JsCvResult, JsFeatureSelection,
};

// core
pub use core_utils::sha256;
#[cfg(feature = "compress")]
pub use core_utils::{zstd_compress, zstd_decompress};

// chem
pub use chem::{
    smiles_properties, canonical, smiles_fingerprint, tanimoto, smiles_substructure,
    parse_sdf, maccs_fingerprint, tanimoto_maccs,
    JsMolecularProperties, JsFingerprint, JsSubstructureResult,
    JsSdfMolecule, JsMaccsFingerprint,
};

// struct_bio
pub use struct_bio::{
    pdb_info, pdb_secondary_structure, rmsd,
    contact_map, ramachandran_analysis, parse_mmcif, kabsch_align,
    JsStructureInfo, JsChainInfo, JsSecondaryStructure, JsSSAssignment,
    JsContactMap, JsRamachandranEntry, JsMmcifInfo, JsKabschResult,
};

// phylo
pub use phylo::{
    newick_info, evolutionary_distance, build_upgma, build_nj, rf_distance,
    parse_nexus, write_nexus, simulate_evolution,
    simulate_coalescent, simulate_coalescent_growth,
    JsTreeInfo, JsRFDistance,
    JsNexusFile, JsNamedTree, JsSimulatedAlignment, JsCoalescentTree,
};

// io
pub use io::{
    pileup_from_sam, depth_stats_from_sam, pileup_to_mpileup_text,
    parse_vcf_text, parse_bed_text, parse_gff3_text, parse_blast_xml,
    parse_bedgraph, parse_gfa, ncbi_fetch_url,
    JsPileupColumn, JsPileup, JsDepthStats,
    JsVcfVariant, JsBedRecord, JsGff3Gene,
    JsBlastXmlHit, JsBlastXmlResult, JsBedGraphRecord,
    JsGfaGraph, JsGfaSegment,
};

// omics
pub use omics::{
    merge_intervals, intersect_intervals, subtract_intervals, complement_intervals,
    closest_intervals, jaccard_intervals, make_windows,
    liftover_interval, annotate_variant, cbs_segment,
    bisulfite_convert, find_cpg_islands, morans_i, gearys_c,
    JsGenomicInterval, JsVariantEffect, JsCnvSegment, JsCpgIsland,
    JsSpatialAutocorrelation, JsLiftoverResult, JsJaccard, JsClosestResult, JsGearysC,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_set() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn reexports_accessible() {
        // Verify that key functions are accessible from the crate root.
        let _ = parse_fasta(">s\nACGT\n");
        let _ = align_dna("A", "A", "global");
        let _ = describe("[1,2,3]");
        let _ = kmer_count("ACGT", 2);
        let _ = reverse_complement("ACGT");
        let _ = transcribe("ACGT");
        let _ = translate("ATGAAA");
        let _ = validate("ACGT", "dna");
        let _ = sha256("hello");
        let _ = spearman("[1,2,3]", "[1,2,3]");
        let _ = bonferroni("[0.01,0.04]");
        let _ = smiles_properties("CCO");
        let _ = canonical("CCO");
        let _ = smiles_fingerprint("CCO", 2, 2048);
        let _ = tanimoto("CCO", "CCO");
        let _ = smiles_substructure("CCO", "CC");
        let _ = pdb_info("HEADER\nEND\n");
        let _ = rmsd("[[0,0,0]]", "[[0,0,0]]");
        let _ = newick_info("(A,B);");
        let _ = evolutionary_distance("ACGT", "ACGT", "p");
        let _ = build_upgma(r#"["A","B"]"#, "[[0,1],[1,0]]");
        let _ = build_nj(r#"["A","B"]"#, "[[0,1],[1,0]]");
        let _ = rf_distance("(A,B);", "(A,B);");
    }
}
