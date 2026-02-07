//! Alphabet definitions for biological sequence validation.
//!
//! Each alphabet is a zero-sized marker type that implements [`Alphabet`],
//! defining the set of valid bytes (uppercase) for a sequence type.

/// Trait for biological sequence alphabets.
///
/// Implementors define a fixed set of valid uppercase bytes. Sequence
/// constructors uppercase input first, then validate against the alphabet.
pub trait Alphabet: Clone + 'static {
    /// Human-readable name (e.g. "DNA").
    const NAME: &'static str;

    /// The set of valid uppercase bytes.
    const VALID_BYTES: &'static [u8];

    /// Check whether a byte (assumed already uppercased) is valid.
    fn is_valid(b: u8) -> bool {
        Self::VALID_BYTES.contains(&b)
    }
}

/// IUPAC DNA alphabet: `ACGTNRYSWKMBDHV`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DnaAlphabet;

impl Alphabet for DnaAlphabet {
    const NAME: &'static str = "DNA";
    const VALID_BYTES: &'static [u8] = b"ACGTNRYSWKMBDHV";
}

/// IUPAC RNA alphabet: `ACGUNRYSWKMBDHV`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RnaAlphabet;

impl Alphabet for RnaAlphabet {
    const NAME: &'static str = "RNA";
    const VALID_BYTES: &'static [u8] = b"ACGUNRYSWKMBDHV";
}

/// Protein alphabet: 20 standard amino acids plus `XBZJUO*`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ProteinAlphabet;

impl Alphabet for ProteinAlphabet {
    const NAME: &'static str = "Protein";
    const VALID_BYTES: &'static [u8] = b"ACDEFGHIKLMNPQRSTVWYXBZJUO*";
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_accepts_all_iupac_bases() {
        for &b in b"ACGTNRYSWKMBDHV" {
            assert!(DnaAlphabet::is_valid(b), "DNA should accept {}", b as char);
        }
    }

    #[test]
    fn dna_rejects_u() {
        assert!(!DnaAlphabet::is_valid(b'U'));
    }

    #[test]
    fn rna_accepts_all_iupac_bases() {
        for &b in b"ACGUNRYSWKMBDHV" {
            assert!(RnaAlphabet::is_valid(b), "RNA should accept {}", b as char);
        }
    }

    #[test]
    fn rna_rejects_t() {
        assert!(!RnaAlphabet::is_valid(b'T'));
    }

    #[test]
    fn protein_accepts_full_set() {
        for &b in b"ACDEFGHIKLMNPQRSTVWYXBZJUO*" {
            assert!(
                ProteinAlphabet::is_valid(b),
                "Protein should accept {}",
                b as char
            );
        }
    }

    #[test]
    fn protein_rejects_invalid() {
        assert!(!ProteinAlphabet::is_valid(b'1'));
        assert!(!ProteinAlphabet::is_valid(b' '));
    }
}
