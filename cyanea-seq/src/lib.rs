//! Sequence I/O and manipulation for the Cyanea bioinformatics ecosystem.
//!
//! Provides strongly-typed, validated biological sequence types with full IUPAC
//! alphabet support, plus FASTA/FASTQ parsing:
//!
//! - **Alphabets** — [`DnaAlphabet`], [`RnaAlphabet`], [`ProteinAlphabet`]
//! - **Sequences** — [`DnaSequence`], [`RnaSequence`], [`ProteinSequence`]
//! - **Codon translation** — Standard genetic code (NCBI Table 1)
//! - **K-mer iteration** — Zero-allocation [`KmerIter`]
//! - **Quality scores** — [`QualityScores`] with Phred+33/64 support
//! - **FASTA parsing** — [`FastaStats`] via [`parse_fasta_stats`]
//! - **FASTQ parsing** — [`FastqRecord`], [`FastqStats`] via [`parse_fastq_file`]
//!
//! # Example
//!
//! ```
//! use cyanea_seq::{DnaSequence, RnaSequence, ProteinSequence};
//! use cyanea_core::Sequence;
//!
//! // Create a DNA sequence (lowercased input is normalized)
//! let dna = DnaSequence::new(b"atgaaagcttaa").unwrap();
//! assert_eq!(dna.as_bytes(), b"ATGAAAGCTTAA");
//!
//! // Reverse complement
//! let rc = dna.reverse_complement();
//! assert_eq!(rc.as_bytes(), b"TTAAGCTTTCAT");
//!
//! // Transcribe DNA → RNA
//! let rna = dna.transcribe();
//! assert_eq!(rna.as_bytes(), b"AUGAAAGCUUAA");
//!
//! // Translate RNA → Protein (stops at UAA)
//! let protein = rna.translate().unwrap();
//! assert_eq!(protein.as_bytes(), b"MKA");
//! ```

pub mod alphabet;
pub mod assembly;
pub mod bwt;
pub mod codon;
pub mod debruijn;
pub mod fasta;
#[cfg(feature = "std")]
pub mod fasta_index;
pub mod fastq;
#[cfg(feature = "std")]
pub mod paired;
pub mod fm_index;
pub mod fmd_index;
pub mod kmer;
pub mod masking;
pub mod minhash;
pub mod motif;
pub mod orf;
pub mod pattern;
pub mod protein_properties;
pub mod pssm;
pub mod restriction;
pub mod rna_structure;
pub mod quality;
pub mod seq;
pub mod suffix;
pub mod taxonomy;
pub mod trim;
pub mod twobit;
pub mod types;

// Re-export alphabet types
pub use alphabet::{Alphabet, DnaAlphabet, ProteinAlphabet, RnaAlphabet};

// Re-export the generic sequence type
pub use seq::ValidatedSeq;

// Re-export concrete type aliases and their methods
pub use types::{DnaSequence, ProteinSequence, RnaSequence};

// Re-export codon translation and analysis
pub use codon::{
    classify_substitution, codon_adaptation_index, count_syn_nonsyn_sites, translate_codon,
    translate_sequence, CodonUsage, GeneticCode, GeneticCodeId, SubstitutionClass,
};

// Re-export k-mer iterator
pub use kmer::KmerIter;

// Re-export quality scores
pub use quality::{PhredEncoding, QualityScores};

// Re-export FASTA types
pub use fasta::{parse_fasta_stats, FastaStats};

// Re-export indexed FASTA types
#[cfg(feature = "std")]
pub use fasta_index::{FastaIndex, FastaIndexEntry, IndexedFastaReader};

// Re-export FASTQ types
pub use fastq::{parse_fastq_file, parse_fastq_stats, FastqRecord, FastqStats};

// Re-export compact encoding and indexing types
pub use twobit::TwoBitSequence;
pub use suffix::SuffixArray;
pub use fm_index::FmIndex;
pub use fmd_index::{FmdIndex, BiInterval};

// Re-export MinHash sketching types
pub use minhash::{MinHash, FracMinHash};

// Re-export pattern matching algorithms
pub use pattern::{bndm, bom, horspool, kmp, myers_bitparallel, shift_and, ukkonen};

// Re-export PSSM types and helpers
pub use pssm::{dna_mapping, protein_mapping, Pssm, PssmDna, PssmProtein};

// Re-export ORF finder
pub use orf::{find_orfs, find_orfs_both_strands, find_orfs_with_codons, OrfResult, Strand};

// Re-export BWT
pub use bwt::Bwt;

// Re-export quality trimming and filtering
pub use trim::{TrimPipeline, TrimRange, TrimReport};
#[cfg(feature = "std")]
pub use trim::{OrphanPolicy, PairedTrimReport, PairedTrimResult};

// Re-export paired-end FASTQ types
#[cfg(feature = "std")]
pub use paired::{
    deinterleave_fastq_file, interleave_fastq_files, parse_interleaved_fastq,
    parse_paired_fastq_files, parse_paired_fastq_stats, strip_read_suffix, validate_mate_pair,
    validate_mate_pair_strict, write_interleaved_fastq, write_paired_fastq, MateValidation,
    PairedFastqRecord, PairedFastqStats,
};

// Re-export protein sequence properties
pub use protein_properties::{
    amino_acid_composition, chou_fasman, extinction_coefficient, gor, gravy,
    hydrophobicity_profile, isoelectric_point, predict_disorder, AminoAcidComposition,
    DisorderPrediction, ExtinctionCoefficient, HydrophobicityScale, SecondaryStructure,
    SecondaryStructurePrediction,
};

// Re-export RNA secondary structure prediction
pub use rna_structure::{
    base_pair_distance, mccaskill, mountain_distance, nussinov, zuker_mfe, MfeResult,
    NussinovResult, PartitionResult, RnaSecondaryStructure,
};

// Re-export masking types
pub use masking::{
    apply_mask, dust, find_tandem_repeats, mask_dust, mask_seg, seg, DustParams, MaskMode,
    MaskResult, MaskSource, MaskedRegion, SegParams, TandemRepeatParams,
};

// Re-export de Bruijn graph types
pub use debruijn::{DeBruijnGraph, Unitig};

// Re-export assembly QC
pub use assembly::{assembly_stats, nx_values, AssemblyStats};

// Re-export taxonomy types
pub use taxonomy::{KmerClassifier, TaxonRank, TaxonomyNode, TaxonomyTree};

// Re-export restriction enzyme types
pub use restriction::{
    common_enzymes, digest, find_cut_sites, fragment_sizes, CutSite, Fragment, Overhang,
    RestrictionEnzyme,
};

// Re-export motif discovery types
pub use motif::{discover_motifs, DiscoveredMotif, Pwm};
