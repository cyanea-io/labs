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
pub mod codon;
pub mod fasta;
#[cfg(feature = "std")]
pub mod fasta_index;
pub mod fastq;
pub mod fm_index;
pub mod fmd_index;
pub mod kmer;
pub mod minhash;
pub mod orf;
pub mod pattern;
pub mod pssm;
pub mod quality;
pub mod seq;
pub mod suffix;
pub mod twobit;
pub mod types;

// Re-export alphabet types
pub use alphabet::{Alphabet, DnaAlphabet, ProteinAlphabet, RnaAlphabet};

// Re-export the generic sequence type
pub use seq::ValidatedSeq;

// Re-export concrete type aliases and their methods
pub use types::{DnaSequence, ProteinSequence, RnaSequence};

// Re-export codon translation
pub use codon::{translate_codon, translate_sequence};

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
