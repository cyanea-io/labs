//! Concrete sequence type aliases and biologically meaningful operations.
//!
//! - [`DnaSequence`] — reverse complement, transcription, translation, GC content
//! - [`RnaSequence`] — reverse complement, reverse transcription, translation
//! - [`ProteinSequence`] — molecular weight

use cyanea_core::Result;

use crate::alphabet::{DnaAlphabet, ProteinAlphabet, RnaAlphabet};
use crate::codon;
use crate::kmer::KmerIter;
use crate::seq::ValidatedSeq;

/// A validated DNA sequence (IUPAC alphabet).
pub type DnaSequence = ValidatedSeq<DnaAlphabet>;

/// A validated RNA sequence (IUPAC alphabet).
pub type RnaSequence = ValidatedSeq<RnaAlphabet>;

/// A validated protein/amino acid sequence.
pub type ProteinSequence = ValidatedSeq<ProteinAlphabet>;

// ---------------------------------------------------------------------------
// DNA complement table (full IUPAC)
// ---------------------------------------------------------------------------

fn dna_complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y', // A|G → T|C
        b'Y' => b'R',
        b'S' => b'S', // G|C → C|G
        b'W' => b'W', // A|T → T|A
        b'K' => b'M', // G|T → C|A
        b'M' => b'K',
        b'B' => b'V', // C|G|T → G|C|A
        b'V' => b'B',
        b'D' => b'H', // A|G|T → T|C|A
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

fn rna_complement(b: u8) -> u8 {
    match b {
        b'A' => b'U',
        b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

// ---------------------------------------------------------------------------
// DNA methods
// ---------------------------------------------------------------------------

impl DnaSequence {
    /// Return the reverse complement.
    pub fn reverse_complement(&self) -> DnaSequence {
        let rc: Vec<u8> = self.iter().rev().map(|&b| dna_complement(b)).collect();
        DnaSequence::from_validated(rc)
    }

    /// Transcribe DNA to RNA (T → U).
    pub fn transcribe(&self) -> RnaSequence {
        let rna: Vec<u8> = self
            .iter()
            .map(|&b| if b == b'T' { b'U' } else { b })
            .collect();
        RnaSequence::from_validated(rna)
    }

    /// Translate DNA to protein (transcribes first, then translates).
    pub fn translate(&self) -> Result<ProteinSequence> {
        self.transcribe().translate()
    }

    /// GC content as a fraction in [0.0, 1.0].
    ///
    /// Only counts unambiguous G and C bases. Returns 0.0 for empty sequences.
    pub fn gc_content(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let gc = self.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc as f64 / self.len() as f64
    }

    /// Iterate over k-mers of length `k`.
    pub fn kmers(&self, k: usize) -> Result<KmerIter<'_>> {
        KmerIter::new(self, k)
    }
}

// ---------------------------------------------------------------------------
// RNA methods
// ---------------------------------------------------------------------------

impl RnaSequence {
    /// Return the reverse complement.
    pub fn reverse_complement(&self) -> RnaSequence {
        let rc: Vec<u8> = self.iter().rev().map(|&b| rna_complement(b)).collect();
        RnaSequence::from_validated(rc)
    }

    /// Reverse-transcribe RNA to DNA (U → T).
    pub fn reverse_transcribe(&self) -> DnaSequence {
        let dna: Vec<u8> = self
            .iter()
            .map(|&b| if b == b'U' { b'T' } else { b })
            .collect();
        DnaSequence::from_validated(dna)
    }

    /// Translate RNA to protein using the standard genetic code.
    pub fn translate(&self) -> Result<ProteinSequence> {
        let protein_bytes = codon::translate_sequence(self.as_ref());
        ProteinSequence::new(&protein_bytes)
    }

    /// Translate in all three forward reading frames.
    pub fn translate_frames(&self) -> [Result<ProteinSequence>; 3] {
        let bytes: &[u8] = self.as_ref();
        [
            ProteinSequence::new(&codon::translate_sequence(bytes)),
            ProteinSequence::new(&codon::translate_sequence(if bytes.len() > 1 {
                &bytes[1..]
            } else {
                &[]
            })),
            ProteinSequence::new(&codon::translate_sequence(if bytes.len() > 2 {
                &bytes[2..]
            } else {
                &[]
            })),
        ]
    }

    /// Iterate over k-mers of length `k`.
    pub fn kmers(&self, k: usize) -> Result<KmerIter<'_>> {
        KmerIter::new(self, k)
    }
}

// ---------------------------------------------------------------------------
// Protein methods
// ---------------------------------------------------------------------------

/// Average molecular weights (Da) for each amino acid.
fn amino_acid_weight(aa: u8) -> f64 {
    match aa {
        b'A' => 89.09, b'R' => 174.20, b'N' => 132.12, b'D' => 133.10,
        b'C' => 121.16, b'E' => 147.13, b'Q' => 146.15, b'G' => 75.03,
        b'H' => 155.16, b'I' => 131.17, b'L' => 131.17, b'K' => 146.19,
        b'M' => 149.21, b'F' => 165.19, b'P' => 115.13, b'S' => 105.09,
        b'T' => 119.12, b'W' => 204.23, b'Y' => 181.19, b'V' => 117.15,
        // Ambiguous / non-standard — use average of all 20
        _ => 128.16,
    }
}

impl ProteinSequence {
    /// Estimated molecular weight in Daltons.
    ///
    /// Sum of residue weights minus (n-1) water molecules lost in peptide bonds.
    pub fn molecular_weight(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.iter().map(|&aa| amino_acid_weight(aa)).sum();
        // Subtract water lost per peptide bond
        let water = 18.015;
        sum - (self.len() as f64 - 1.0) * water
    }

    /// Iterate over k-mers of length `k`.
    pub fn kmers(&self, k: usize) -> Result<KmerIter<'_>> {
        KmerIter::new(self, k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Reverse complement ---

    #[test]
    fn revcomp_palindromic() {
        let seq = DnaSequence::new(b"ACGT").unwrap();
        assert_eq!(seq.reverse_complement().as_ref(), b"ACGT");
    }

    #[test]
    fn revcomp_asymmetric() {
        let seq = DnaSequence::new(b"AACG").unwrap();
        assert_eq!(seq.reverse_complement().as_ref(), b"CGTT");
    }

    #[test]
    fn revcomp_iupac_ambiguity() {
        let seq = DnaSequence::new(b"RYSWKMBDHVN").unwrap();
        let rc = seq.reverse_complement();
        assert_eq!(rc.as_ref(), b"NBDHVKMWSRY");
    }

    // --- Transcription ---

    #[test]
    fn dna_to_rna() {
        let dna = DnaSequence::new(b"ATCG").unwrap();
        let rna = dna.transcribe();
        assert_eq!(rna.as_ref(), b"AUCG");
    }

    #[test]
    fn rna_to_dna() {
        let rna = RnaSequence::new(b"AUCG").unwrap();
        let dna = rna.reverse_transcribe();
        assert_eq!(dna.as_ref(), b"ATCG");
    }

    #[test]
    fn transcription_roundtrip() {
        let dna = DnaSequence::new(b"ATCGATCG").unwrap();
        let roundtrip = dna.transcribe().reverse_transcribe();
        assert_eq!(dna, roundtrip);
    }

    // --- Translation ---

    #[test]
    fn translate_aug_uuu_uaa() {
        // AUG=M, UUU=F, UAA=stop → "MF"
        let rna = RnaSequence::new(b"AUGUUUUAA").unwrap();
        let protein = rna.translate().unwrap();
        assert_eq!(protein.as_ref(), b"MF");
    }

    #[test]
    fn translate_stops_at_first_stop() {
        let rna = RnaSequence::new(b"AUGUUUUAAGCU").unwrap();
        let protein = rna.translate().unwrap();
        assert_eq!(protein.as_ref(), b"MF");
    }

    #[test]
    fn translate_incomplete_codon_ignored() {
        let rna = RnaSequence::new(b"AUGUUUAU").unwrap();
        let protein = rna.translate().unwrap();
        assert_eq!(protein.as_ref(), b"MF");
    }

    // --- GC content ---

    #[test]
    fn gc_content_basic() {
        let seq = DnaSequence::new(b"ATGC").unwrap();
        assert!((seq.gc_content() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn gc_content_empty() {
        let seq = DnaSequence::new(b"").unwrap();
        assert_eq!(seq.gc_content(), 0.0);
    }
}
