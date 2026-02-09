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

/// Crate version (set from Cargo.toml at compile time).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ── Re-exports ───────────────────────────────────────────────────────────

// seq
pub use seq::{
    parse_fasta, parse_fasta_bytes, gc_content, gc_content_json,
    reverse_complement, transcribe, translate, validate, parse_fastq,
};

// align
pub use align::{align_dna, align_dna_custom, align_protein, align_batch};

// stats
pub use stats::{
    describe, pearson, spearman, t_test, t_test_two_sample,
    mann_whitney_u, bonferroni, benjamini_hochberg,
    JsDescriptiveStats, JsTestResult,
};

// ml
pub use ml::{
    kmer_count, euclidean_distance, manhattan_distance, hamming_distance, cosine_similarity,
    JsKmerCounts,
};

// core
pub use core_utils::{sha256, zstd_compress, zstd_decompress};

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
    }
}
