//! Codon translation using the NCBI standard genetic code (Table 1).
//!
//! Operates on raw `&[u8]` triplets. Handles both DNA (T) and RNA (U) codons.

/// Translate a single codon (3-byte slice) to an amino acid.
///
/// Accepts both DNA (T) and RNA (U) codons. Returns `None` for stop codons
/// (UAA/TAA, UAG/TAG, UGA/TGA) and unrecognized/ambiguous codons.
pub fn translate_codon(codon: &[u8]) -> Option<u8> {
    if codon.len() != 3 {
        return None;
    }
    // Normalize: uppercase and T → U
    let normalize = |b: u8| match b.to_ascii_uppercase() {
        b'T' => b'U',
        other => other,
    };
    let c = [normalize(codon[0]), normalize(codon[1]), normalize(codon[2])];
    match &c {
        // Phe
        b"UUU" | b"UUC" => Some(b'F'),
        // Leu
        b"UUA" | b"UUG" | b"CUU" | b"CUC" | b"CUA" | b"CUG" => Some(b'L'),
        // Ile
        b"AUU" | b"AUC" | b"AUA" => Some(b'I'),
        // Met (start)
        b"AUG" => Some(b'M'),
        // Val
        b"GUU" | b"GUC" | b"GUA" | b"GUG" => Some(b'V'),
        // Ser (UCx)
        b"UCU" | b"UCC" | b"UCA" | b"UCG" => Some(b'S'),
        // Pro
        b"CCU" | b"CCC" | b"CCA" | b"CCG" => Some(b'P'),
        // Thr
        b"ACU" | b"ACC" | b"ACA" | b"ACG" => Some(b'T'),
        // Ala
        b"GCU" | b"GCC" | b"GCA" | b"GCG" => Some(b'A'),
        // Tyr
        b"UAU" | b"UAC" => Some(b'Y'),
        // Stop
        b"UAA" | b"UAG" | b"UGA" => None,
        // His
        b"CAU" | b"CAC" => Some(b'H'),
        // Gln
        b"CAA" | b"CAG" => Some(b'Q'),
        // Asn
        b"AAU" | b"AAC" => Some(b'N'),
        // Lys
        b"AAA" | b"AAG" => Some(b'K'),
        // Asp
        b"GAU" | b"GAC" => Some(b'D'),
        // Glu
        b"GAA" | b"GAG" => Some(b'E'),
        // Cys
        b"UGU" | b"UGC" => Some(b'C'),
        // Trp
        b"UGG" => Some(b'W'),
        // Arg (CGx)
        b"CGU" | b"CGC" | b"CGA" | b"CGG" => Some(b'R'),
        // Ser (AGC/AGU)
        b"AGU" | b"AGC" => Some(b'S'),
        // Arg (AGR)
        b"AGA" | b"AGG" => Some(b'R'),
        // Gly
        b"GGU" | b"GGC" | b"GGA" | b"GGG" => Some(b'G'),
        // Ambiguous or invalid
        _ => None,
    }
}

/// Translate a nucleotide sequence into a protein sequence.
///
/// Processes codons (3-byte chunks) from left to right. Stops at the first
/// stop codon or end of sequence. Incomplete trailing codons are ignored.
pub fn translate_sequence(seq: &[u8]) -> Vec<u8> {
    let mut protein = Vec::with_capacity(seq.len() / 3);
    for codon in seq.chunks_exact(3) {
        match translate_codon(codon) {
            Some(aa) => protein.push(aa),
            None => break, // stop codon or ambiguous
        }
    }
    protein
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_64_standard_codons() {
        // NCBI standard genetic code (Table 1)
        let expected: &[(&[u8], Option<u8>)] = &[
            (b"UUU", Some(b'F')), (b"UUC", Some(b'F')),
            (b"UUA", Some(b'L')), (b"UUG", Some(b'L')),
            (b"CUU", Some(b'L')), (b"CUC", Some(b'L')),
            (b"CUA", Some(b'L')), (b"CUG", Some(b'L')),
            (b"AUU", Some(b'I')), (b"AUC", Some(b'I')),
            (b"AUA", Some(b'I')), (b"AUG", Some(b'M')),
            (b"GUU", Some(b'V')), (b"GUC", Some(b'V')),
            (b"GUA", Some(b'V')), (b"GUG", Some(b'V')),
            (b"UCU", Some(b'S')), (b"UCC", Some(b'S')),
            (b"UCA", Some(b'S')), (b"UCG", Some(b'S')),
            (b"CCU", Some(b'P')), (b"CCC", Some(b'P')),
            (b"CCA", Some(b'P')), (b"CCG", Some(b'P')),
            (b"ACU", Some(b'T')), (b"ACC", Some(b'T')),
            (b"ACA", Some(b'T')), (b"ACG", Some(b'T')),
            (b"GCU", Some(b'A')), (b"GCC", Some(b'A')),
            (b"GCA", Some(b'A')), (b"GCG", Some(b'A')),
            (b"UAU", Some(b'Y')), (b"UAC", Some(b'Y')),
            (b"UAA", None),       (b"UAG", None),
            (b"CAU", Some(b'H')), (b"CAC", Some(b'H')),
            (b"CAA", Some(b'Q')), (b"CAG", Some(b'Q')),
            (b"AAU", Some(b'N')), (b"AAC", Some(b'N')),
            (b"AAA", Some(b'K')), (b"AAG", Some(b'K')),
            (b"GAU", Some(b'D')), (b"GAC", Some(b'D')),
            (b"GAA", Some(b'E')), (b"GAG", Some(b'E')),
            (b"UGU", Some(b'C')), (b"UGC", Some(b'C')),
            (b"UGA", None),       (b"UGG", Some(b'W')),
            (b"CGU", Some(b'R')), (b"CGC", Some(b'R')),
            (b"CGA", Some(b'R')), (b"CGG", Some(b'R')),
            (b"AGU", Some(b'S')), (b"AGC", Some(b'S')),
            (b"AGA", Some(b'R')), (b"AGG", Some(b'R')),
            (b"GGU", Some(b'G')), (b"GGC", Some(b'G')),
            (b"GGA", Some(b'G')), (b"GGG", Some(b'G')),
        ];
        for &(codon, expected_aa) in expected {
            let result = translate_codon(codon);
            assert_eq!(
                result, expected_aa,
                "codon {} expected {:?}, got {:?}",
                std::str::from_utf8(codon).unwrap(),
                expected_aa.map(|b| b as char),
                result.map(|b| b as char),
            );
        }
    }

    #[test]
    fn translate_dna_codons() {
        // T-based codons should work too
        assert_eq!(translate_codon(b"ATG"), Some(b'M'));
        assert_eq!(translate_codon(b"TAA"), None);
    }

    #[test]
    fn translate_sequence_basic() {
        // AUG UUU UAA → MF (stop)
        let protein = translate_sequence(b"AUGUUUUAA");
        assert_eq!(protein, b"MF");
    }

    #[test]
    fn translate_sequence_incomplete_codon_ignored() {
        // AUG UUU AU (incomplete) → MF
        let protein = translate_sequence(b"AUGUUUAU");
        assert_eq!(protein, b"MF");
    }
}
