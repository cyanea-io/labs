//! Sequence encoding for ML feature extraction.

/// Biological sequence alphabets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Alphabet {
    /// DNA: A, C, G, T
    Dna,
    /// RNA: A, C, G, U
    Rna,
    /// Protein: 20 standard amino acids
    Protein,
}

const DNA_SYMBOLS: &[u8] = b"ACGT";
const RNA_SYMBOLS: &[u8] = b"ACGU";
const PROTEIN_SYMBOLS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

impl Alphabet {
    /// The ordered byte symbols for this alphabet.
    pub fn symbols(&self) -> &[u8] {
        match self {
            Alphabet::Dna => DNA_SYMBOLS,
            Alphabet::Rna => RNA_SYMBOLS,
            Alphabet::Protein => PROTEIN_SYMBOLS,
        }
    }

    /// The number of symbols in this alphabet.
    pub fn size(&self) -> usize {
        self.symbols().len()
    }
}

/// One-hot encode a sequence using the given alphabet.
///
/// Each position produces `alphabet.size()` values: 1.0 at the matching index,
/// 0.0 elsewhere. Unknown characters produce an all-zero row.
///
/// The output is a flat vector of length `seq.len() * alphabet.size()`.
pub fn one_hot_encode(seq: &[u8], alphabet: Alphabet) -> Vec<f64> {
    let syms = alphabet.symbols();
    let k = syms.len();
    let mut out = vec![0.0; seq.len() * k];
    for (i, &base) in seq.iter().enumerate() {
        let upper = base.to_ascii_uppercase();
        if let Some(pos) = syms.iter().position(|&s| s == upper) {
            out[i * k + pos] = 1.0;
        }
    }
    out
}

/// Label-encode a sequence using the given alphabet.
///
/// Each position maps to its index in the alphabet (0-based) as `f64`.
/// Unknown characters are encoded as -1.0.
///
/// The output has one value per position.
pub fn label_encode(seq: &[u8], alphabet: Alphabet) -> Vec<f64> {
    let syms = alphabet.symbols();
    seq.iter()
        .map(|&base| {
            let upper = base.to_ascii_uppercase();
            match syms.iter().position(|&s| s == upper) {
                Some(pos) => pos as f64,
                None => -1.0,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_symbols() {
        assert_eq!(Alphabet::Dna.symbols(), b"ACGT");
        assert_eq!(Alphabet::Dna.size(), 4);
    }

    #[test]
    fn rna_symbols() {
        assert_eq!(Alphabet::Rna.symbols(), b"ACGU");
        assert_eq!(Alphabet::Rna.size(), 4);
    }

    #[test]
    fn protein_symbols() {
        assert_eq!(Alphabet::Protein.size(), 20);
    }

    #[test]
    fn one_hot_dna() {
        let enc = one_hot_encode(b"ACG", Alphabet::Dna);
        // A -> [1,0,0,0], C -> [0,1,0,0], G -> [0,0,1,0]
        assert_eq!(enc.len(), 12);
        assert_eq!(enc[0], 1.0); // A
        assert_eq!(enc[5], 1.0); // C
        assert_eq!(enc[10], 1.0); // G
    }

    #[test]
    fn one_hot_unknown_char() {
        let enc = one_hot_encode(b"N", Alphabet::Dna);
        assert_eq!(enc, vec![0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn one_hot_lowercase() {
        let enc = one_hot_encode(b"a", Alphabet::Dna);
        assert_eq!(enc[0], 1.0);
    }

    #[test]
    fn one_hot_empty() {
        let enc = one_hot_encode(b"", Alphabet::Dna);
        assert!(enc.is_empty());
    }

    #[test]
    fn label_encode_dna() {
        let enc = label_encode(b"ACGT", Alphabet::Dna);
        assert_eq!(enc, vec![0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn label_encode_unknown() {
        let enc = label_encode(b"N", Alphabet::Dna);
        assert_eq!(enc, vec![-1.0]);
    }

    #[test]
    fn label_encode_protein() {
        let enc = label_encode(b"AW", Alphabet::Protein);
        assert_eq!(enc[0], 0.0); // A is first
        // W is at index 18 in ACDEFGHIKLMNPQRSTVWY
        assert_eq!(enc[1], 18.0);
    }
}
