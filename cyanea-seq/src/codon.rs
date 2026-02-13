//! Codon translation, genetic codes, codon usage analysis, and CAI.
//!
//! Supports 7 NCBI genetic code tables (1, 2, 3, 4, 5, 6, 11) with
//! full translation, start/stop codon queries, codon usage frequency
//! tables, RSCU, CAI, and synonymous/non-synonymous site counting.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Base encoding: A=0, C=1, G=2, T/U=3
// ---------------------------------------------------------------------------

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Convert a 3-base codon to an index in [0, 64).
fn codon_index(codon: &[u8]) -> Option<usize> {
    if codon.len() != 3 {
        return None;
    }
    let b1 = base_index(codon[0])?;
    let b2 = base_index(codon[1])?;
    let b3 = base_index(codon[2])?;
    Some(b1 * 16 + b2 * 4 + b3)
}

/// Convert an index in [0, 64) back to a codon (as DNA: A/C/G/T).
fn index_to_codon(idx: usize) -> [u8; 3] {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    [
        BASES[idx >> 4],
        BASES[(idx >> 2) & 3],
        BASES[idx & 3],
    ]
}

// ---------------------------------------------------------------------------
// Genetic code table identifier
// ---------------------------------------------------------------------------

/// NCBI genetic code table identifier.
///
/// Covers the 7 most commonly used tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GeneticCodeId {
    Standard = 1,
    VertebrateMitochondrial = 2,
    YeastMitochondrial = 3,
    MycoplasmaSpiroplasma = 4,
    InvertebrateMitochondrial = 5,
    CiliateNuclear = 6,
    BacterialPlastid = 11,
}

// ---------------------------------------------------------------------------
// Genetic code tables (const arrays)
// ---------------------------------------------------------------------------

// Codon order: AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT,
//              ATA, ATC, ATG, ATT, CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT,
//              CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT, GAA, GAC, GAG, GAT,
//              GCA, GCC, GCG, GCT, GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT,
//              TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT, TGA, TGC, TGG, TGT,
//              TTA, TTC, TTG, TTT

/// Standard genetic code (NCBI Table 1).
const TABLE1_AA: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S',
    b'I', b'I', b'M', b'I', b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P',
    b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L', b'E', b'D', b'E', b'D',
    b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C',
    b'L', b'F', b'L', b'F',
];

const TABLE1_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    // ATG
    s[14] = true;
    s
};

/// Vertebrate mitochondrial (NCBI Table 2).
/// UGA=Trp, AGA=Stop, AGG=Stop, AUA=Met
const TABLE2_AA: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'*', b'S', b'*', b'S',
    b'M', b'I', b'M', b'I', b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P',
    b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L', b'E', b'D', b'E', b'D',
    b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'W', b'C', b'W', b'C',
    b'L', b'F', b'L', b'F',
];

const TABLE2_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    s[14] = true; // ATG
    s[12] = true; // ATA
    s[13] = true; // ATC
    s[15] = true; // ATT
    s[46] = true; // GTG
    s
};

/// Yeast mitochondrial (NCBI Table 3).
/// CUG=Thr (not Leu), UGA=Trp, CUA=Thr (not Leu)
const TABLE3_AA: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S',
    b'M', b'I', b'M', b'I', b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P',
    b'R', b'R', b'R', b'R', b'T', b'L', b'T', b'L', b'E', b'D', b'E', b'D',
    b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'W', b'C', b'W', b'C',
    b'L', b'F', b'L', b'F',
];

const TABLE3_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    s[14] = true; // ATG
    s[12] = true; // ATA
    s
};

/// Mycoplasma/Spiroplasma (NCBI Table 4).
/// UGA=Trp
const TABLE4_AA: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S',
    b'I', b'I', b'M', b'I', b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P',
    b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L', b'E', b'D', b'E', b'D',
    b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'W', b'C', b'W', b'C',
    b'L', b'F', b'L', b'F',
];

const TABLE4_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    s[14] = true; // ATG
    s[62] = true; // TTG
    s[46] = true; // GTG
    s
};

/// Invertebrate mitochondrial (NCBI Table 5).
/// AGA=Ser, AGG=Ser, UGA=Trp, AUA=Met
const TABLE5_AA: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'S', b'S', b'S', b'S',
    b'M', b'I', b'M', b'I', b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P',
    b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L', b'E', b'D', b'E', b'D',
    b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'W', b'C', b'W', b'C',
    b'L', b'F', b'L', b'F',
];

const TABLE5_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    s[14] = true; // ATG
    s[15] = true; // ATT
    s[46] = true; // GTG
    s
};

/// Ciliate nuclear (NCBI Table 6).
/// UAA=Gln, UAG=Gln (not Stop)
const TABLE6_AA: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S',
    b'I', b'I', b'M', b'I', b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P',
    b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L', b'E', b'D', b'E', b'D',
    b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'Q', b'Y', b'Q', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C',
    b'L', b'F', b'L', b'F',
];

const TABLE6_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    s[14] = true; // ATG
    s
};

/// Bacterial/plant plastid (NCBI Table 11).
/// Same AAs as standard, different starts.
const TABLE11_AA: [u8; 64] = TABLE1_AA;

const TABLE11_STARTS: [bool; 64] = {
    let mut s = [false; 64];
    s[14] = true; // ATG
    s[62] = true; // TTG
    s[46] = true; // GTG
    s[13] = true; // ATC
    s[15] = true; // ATT
    s[28] = true; // CTA (CTG in some refs — use the common bacterial ones)
    s
};

// ---------------------------------------------------------------------------
// GeneticCode
// ---------------------------------------------------------------------------

/// A genetic code translation table.
///
/// Wraps a 64-element amino acid lookup array and a 64-element start codon mask.
/// Use [`GeneticCode::from_id`] or [`GeneticCode::standard`] to create.
#[derive(Debug, Clone)]
pub struct GeneticCode {
    id: GeneticCodeId,
    name: &'static str,
    table: [u8; 64],
    starts: [bool; 64],
}

impl GeneticCode {
    /// Create a genetic code table from an NCBI table identifier.
    pub fn from_id(id: GeneticCodeId) -> Self {
        match id {
            GeneticCodeId::Standard => Self {
                id,
                name: "Standard",
                table: TABLE1_AA,
                starts: TABLE1_STARTS,
            },
            GeneticCodeId::VertebrateMitochondrial => Self {
                id,
                name: "Vertebrate Mitochondrial",
                table: TABLE2_AA,
                starts: TABLE2_STARTS,
            },
            GeneticCodeId::YeastMitochondrial => Self {
                id,
                name: "Yeast Mitochondrial",
                table: TABLE3_AA,
                starts: TABLE3_STARTS,
            },
            GeneticCodeId::MycoplasmaSpiroplasma => Self {
                id,
                name: "Mycoplasma/Spiroplasma",
                table: TABLE4_AA,
                starts: TABLE4_STARTS,
            },
            GeneticCodeId::InvertebrateMitochondrial => Self {
                id,
                name: "Invertebrate Mitochondrial",
                table: TABLE5_AA,
                starts: TABLE5_STARTS,
            },
            GeneticCodeId::CiliateNuclear => Self {
                id,
                name: "Ciliate Nuclear",
                table: TABLE6_AA,
                starts: TABLE6_STARTS,
            },
            GeneticCodeId::BacterialPlastid => Self {
                id,
                name: "Bacterial/Plant Plastid",
                table: TABLE11_AA,
                starts: TABLE11_STARTS,
            },
        }
    }

    /// Create the standard genetic code (NCBI Table 1).
    pub fn standard() -> Self {
        Self::from_id(GeneticCodeId::Standard)
    }

    /// Table identifier.
    pub fn id(&self) -> GeneticCodeId {
        self.id
    }

    /// Human-readable name.
    pub fn name(&self) -> &str {
        self.name
    }

    /// Translate a single codon (3-byte slice) to an amino acid.
    ///
    /// Returns `None` for stop codons (encoded as `b'*'` internally) and
    /// for unrecognizable/ambiguous codons.
    pub fn translate_codon(&self, codon: &[u8]) -> Option<u8> {
        let idx = codon_index(codon)?;
        let aa = self.table[idx];
        if aa == b'*' {
            None
        } else {
            Some(aa)
        }
    }

    /// Translate a nucleotide sequence, stopping at the first stop codon.
    ///
    /// Incomplete trailing codons are ignored.
    pub fn translate_sequence(&self, seq: &[u8]) -> Vec<u8> {
        let mut protein = Vec::with_capacity(seq.len() / 3);
        for codon in seq.chunks_exact(3) {
            match self.translate_codon(codon) {
                Some(aa) => protein.push(aa),
                None => break,
            }
        }
        protein
    }

    /// Translate a nucleotide sequence, including stop codons as `'*'`.
    ///
    /// Does not stop at stop codons — translates the entire sequence.
    /// Incomplete trailing codons are ignored.
    pub fn translate_sequence_full(&self, seq: &[u8]) -> Vec<u8> {
        let mut protein = Vec::with_capacity(seq.len() / 3);
        for codon in seq.chunks_exact(3) {
            if let Some(idx) = codon_index(codon) {
                protein.push(self.table[idx]);
            }
        }
        protein
    }

    /// Check whether a codon is a start codon in this table.
    pub fn is_start(&self, codon: &[u8]) -> bool {
        codon_index(codon).map_or(false, |idx| self.starts[idx])
    }

    /// Check whether a codon is a stop codon in this table.
    pub fn is_stop(&self, codon: &[u8]) -> bool {
        codon_index(codon).map_or(false, |idx| self.table[idx] == b'*')
    }

    /// Return all stop codons for this table (as DNA).
    pub fn stop_codons(&self) -> Vec<[u8; 3]> {
        (0..64)
            .filter(|&i| self.table[i] == b'*')
            .map(index_to_codon)
            .collect()
    }

    /// Return all start codons for this table (as DNA).
    pub fn start_codons(&self) -> Vec<[u8; 3]> {
        (0..64)
            .filter(|&i| self.starts[i])
            .map(index_to_codon)
            .collect()
    }

    /// Get the amino acid for a codon index (0..64), including `b'*'` for stops.
    fn aa_at(&self, idx: usize) -> u8 {
        self.table[idx]
    }
}

// ---------------------------------------------------------------------------
// Backward-compatible free functions (delegate to standard code)
// ---------------------------------------------------------------------------

/// Translate a single codon using the standard genetic code (NCBI Table 1).
///
/// Accepts both DNA (T) and RNA (U) codons. Returns `None` for stop codons
/// and unrecognized/ambiguous codons.
pub fn translate_codon(codon: &[u8]) -> Option<u8> {
    GeneticCode::standard().translate_codon(codon)
}

/// Translate a nucleotide sequence using the standard genetic code.
///
/// Stops at the first stop codon. Incomplete trailing codons are ignored.
pub fn translate_sequence(seq: &[u8]) -> Vec<u8> {
    GeneticCode::standard().translate_sequence(seq)
}

// ---------------------------------------------------------------------------
// Codon usage
// ---------------------------------------------------------------------------

/// Codon usage frequency table.
///
/// Counts occurrences of each of the 64 codons in a coding sequence.
#[derive(Debug, Clone)]
pub struct CodonUsage {
    counts: [u64; 64],
    total: u64,
}

impl CodonUsage {
    /// Build a codon usage table from a coding sequence.
    ///
    /// The sequence should be in-frame coding DNA/RNA. Non-standard bases
    /// and incomplete trailing codons are skipped.
    pub fn from_sequence(seq: &[u8]) -> Self {
        let mut counts = [0u64; 64];
        let mut total = 0u64;
        for codon in seq.chunks_exact(3) {
            if let Some(idx) = codon_index(codon) {
                counts[idx] += 1;
                total += 1;
            }
        }
        CodonUsage { counts, total }
    }

    /// Merge another usage table into this one.
    pub fn merge(&mut self, other: &CodonUsage) {
        for i in 0..64 {
            self.counts[i] += other.counts[i];
        }
        self.total += other.total;
    }

    /// Raw count for a codon.
    pub fn count(&self, codon: &[u8]) -> u64 {
        codon_index(codon).map_or(0, |idx| self.counts[idx])
    }

    /// Relative frequency of a codon (count / total codons).
    pub fn frequency(&self, codon: &[u8]) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.count(codon) as f64 / self.total as f64
    }

    /// Relative Synonymous Codon Usage for a codon.
    ///
    /// RSCU = observed / expected, where expected = (total for AA) / (number of synonymous codons).
    /// An RSCU of 1.0 means no bias.
    pub fn rscu(&self, codon: &[u8], code: &GeneticCode) -> f64 {
        let idx = match codon_index(codon) {
            Some(i) => i,
            None => return 0.0,
        };
        let aa = code.aa_at(idx);
        if aa == b'*' {
            return 0.0;
        }

        // Find all synonymous codons
        let synonymous: Vec<usize> = (0..64)
            .filter(|&i| code.aa_at(i) == aa)
            .collect();
        let n_syn = synonymous.len() as f64;
        let total_syn: u64 = synonymous.iter().map(|&i| self.counts[i]).sum();

        if total_syn == 0 {
            return 0.0;
        }

        let observed = self.counts[idx] as f64;
        let expected = total_syn as f64 / n_syn;
        observed / expected
    }

    /// Compute the relative adaptiveness (w_i) for all 64 codons.
    ///
    /// For each amino acid, w_i = RSCU_i / max(RSCU for that AA).
    /// Stop codons get w = 0.0.
    pub fn relative_adaptiveness(&self, code: &GeneticCode) -> [f64; 64] {
        let mut w = [0.0f64; 64];

        // Group codons by amino acid
        for aa in b"ACDEFGHIKLMNPQRSTVWY".iter() {
            let synonymous: Vec<usize> = (0..64)
                .filter(|&i| code.aa_at(i) == *aa)
                .collect();
            if synonymous.is_empty() {
                continue;
            }

            let rscu_vals: Vec<f64> = synonymous
                .iter()
                .map(|&i| {
                    let total_syn: u64 = synonymous.iter().map(|&j| self.counts[j]).sum();
                    if total_syn == 0 {
                        return 0.0;
                    }
                    let n_syn = synonymous.len() as f64;
                    self.counts[i] as f64 / (total_syn as f64 / n_syn)
                })
                .collect();

            let max_rscu = rscu_vals
                .iter()
                .cloned()
                .fold(0.0f64, f64::max);

            if max_rscu > 0.0 {
                for (k, &idx) in synonymous.iter().enumerate() {
                    w[idx] = rscu_vals[k] / max_rscu;
                }
            }
        }

        w
    }
}

// ---------------------------------------------------------------------------
// CAI
// ---------------------------------------------------------------------------

/// Compute the Codon Adaptation Index of a query sequence against a reference usage table.
///
/// CAI is the geometric mean of relative adaptiveness values for each codon.
/// Returns a value in (0, 1]. Higher values indicate better adaptation.
///
/// # Errors
///
/// Returns an error if the query sequence has no valid sense codons.
pub fn codon_adaptation_index(
    query: &[u8],
    reference: &CodonUsage,
    code: &GeneticCode,
) -> Result<f64> {
    let w = reference.relative_adaptiveness(code);
    let mut log_sum = 0.0f64;
    let mut n = 0u64;

    for codon in query.chunks_exact(3) {
        if let Some(idx) = codon_index(codon) {
            let aa = code.aa_at(idx);
            if aa == b'*' {
                continue; // skip stop codons
            }
            let wi = w[idx];
            if wi > 0.0 {
                log_sum += wi.ln();
                n += 1;
            }
        }
    }

    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "no valid sense codons in query".into(),
        ));
    }

    Ok((log_sum / n as f64).exp())
}

// ---------------------------------------------------------------------------
// Synonymous / non-synonymous classification
// ---------------------------------------------------------------------------

/// Classification of a single-nucleotide substitution between two codons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstitutionClass {
    /// Both codons encode the same amino acid.
    Synonymous,
    /// The codons encode different amino acids (neither is a stop).
    NonSynonymous,
    /// At least one codon is a stop codon.
    StopInvolved,
}

/// Classify the substitution between two codons that differ at exactly one position.
///
/// If the codons are identical or differ at more than one position, the
/// classification is based on the amino acid change (or lack thereof).
pub fn classify_substitution(
    codon1: &[u8],
    codon2: &[u8],
    code: &GeneticCode,
) -> SubstitutionClass {
    let idx1 = match codon_index(codon1) {
        Some(i) => i,
        None => return SubstitutionClass::StopInvolved,
    };
    let idx2 = match codon_index(codon2) {
        Some(i) => i,
        None => return SubstitutionClass::StopInvolved,
    };

    let aa1 = code.aa_at(idx1);
    let aa2 = code.aa_at(idx2);

    if aa1 == b'*' || aa2 == b'*' {
        SubstitutionClass::StopInvolved
    } else if aa1 == aa2 {
        SubstitutionClass::Synonymous
    } else {
        SubstitutionClass::NonSynonymous
    }
}

/// Count the number of synonymous and non-synonymous sites in a coding sequence.
///
/// Uses the Nei-Gojobori method: for each codon, consider all possible
/// single-nucleotide changes and classify as synonymous or non-synonymous.
/// Each codon contributes 3 sites total.
pub fn count_syn_nonsyn_sites(seq: &[u8], code: &GeneticCode) -> (f64, f64) {
    let mut syn = 0.0f64;
    let mut nonsyn = 0.0f64;
    let bases: [u8; 4] = [b'A', b'C', b'G', b'T'];

    for codon in seq.chunks_exact(3) {
        let idx = match codon_index(codon) {
            Some(i) => i,
            None => continue,
        };
        let aa = code.aa_at(idx);
        if aa == b'*' {
            continue;
        }

        // For each position in the codon
        for pos in 0..3 {
            let mut s = 0.0f64;
            let mut n = 0.0f64;
            let orig = codon[pos].to_ascii_uppercase();
            let orig = if orig == b'U' { b'T' } else { orig };

            for &base in &bases {
                if base == orig {
                    continue;
                }
                let mut mutant = [codon[0], codon[1], codon[2]];
                mutant[pos] = base;
                let midx = match codon_index(&mutant) {
                    Some(i) => i,
                    None => continue,
                };
                let maa = code.aa_at(midx);
                if maa == b'*' {
                    // Stop mutation — count as non-synonymous
                    n += 1.0;
                } else if maa == aa {
                    s += 1.0;
                } else {
                    n += 1.0;
                }
            }

            let total = s + n;
            if total > 0.0 {
                syn += s / total;
                nonsyn += n / total;
            }
        }
    }

    (syn, nonsyn)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_all_64_codons() {
        let code = GeneticCode::standard();
        // AAA=K, AAC=N, AAG=K, AAT=N
        assert_eq!(code.translate_codon(b"AAA"), Some(b'K'));
        assert_eq!(code.translate_codon(b"AAC"), Some(b'N'));
        assert_eq!(code.translate_codon(b"AAG"), Some(b'K'));
        assert_eq!(code.translate_codon(b"AAT"), Some(b'N'));
        // ATG=M (start)
        assert_eq!(code.translate_codon(b"ATG"), Some(b'M'));
        // TAA=stop, TAG=stop, TGA=stop
        assert_eq!(code.translate_codon(b"TAA"), None);
        assert_eq!(code.translate_codon(b"TAG"), None);
        assert_eq!(code.translate_codon(b"TGA"), None);
        // TGG=W
        assert_eq!(code.translate_codon(b"TGG"), Some(b'W'));
        // TTT=F, TTC=F
        assert_eq!(code.translate_codon(b"TTT"), Some(b'F'));
        assert_eq!(code.translate_codon(b"TTC"), Some(b'F'));
        // GGG=G
        assert_eq!(code.translate_codon(b"GGG"), Some(b'G'));
    }

    #[test]
    fn backward_compat_translate_codon() {
        assert_eq!(translate_codon(b"ATG"), Some(b'M'));
        assert_eq!(translate_codon(b"TAA"), None);
        assert_eq!(translate_codon(b"UGG"), Some(b'W'));
    }

    #[test]
    fn backward_compat_translate_sequence() {
        let protein = translate_sequence(b"AUGUUUUAA");
        assert_eq!(protein, b"MF");
    }

    #[test]
    fn translate_dna_codons() {
        assert_eq!(translate_codon(b"ATG"), Some(b'M'));
        assert_eq!(translate_codon(b"TAA"), None);
    }

    #[test]
    fn translate_sequence_basic() {
        let protein = translate_sequence(b"AUGUUUUAA");
        assert_eq!(protein, b"MF");
    }

    #[test]
    fn translate_sequence_incomplete_codon_ignored() {
        let protein = translate_sequence(b"AUGUUUAU");
        assert_eq!(protein, b"MF");
    }

    #[test]
    fn translate_sequence_full_includes_stops() {
        let code = GeneticCode::standard();
        let protein = code.translate_sequence_full(b"ATGTAAGCG");
        assert_eq!(protein, &[b'M', b'*', b'A']);
    }

    #[test]
    fn rna_codons_work() {
        let code = GeneticCode::standard();
        assert_eq!(code.translate_codon(b"AUG"), Some(b'M'));
        assert_eq!(code.translate_codon(b"UAA"), None);
    }

    #[test]
    fn vertebrate_mito_differences() {
        let code = GeneticCode::from_id(GeneticCodeId::VertebrateMitochondrial);
        // UGA=Trp (not stop)
        assert_eq!(code.translate_codon(b"TGA"), Some(b'W'));
        // AGA=Stop (not Arg)
        assert_eq!(code.translate_codon(b"AGA"), None);
        // AGG=Stop (not Arg)
        assert_eq!(code.translate_codon(b"AGG"), None);
        // AUA=Met (not Ile)
        assert_eq!(code.translate_codon(b"ATA"), Some(b'M'));
    }

    #[test]
    fn yeast_mito_differences() {
        let code = GeneticCode::from_id(GeneticCodeId::YeastMitochondrial);
        // CUG=Thr (not Leu)
        assert_eq!(code.translate_codon(b"CTG"), Some(b'T'));
        // CUA=Thr (not Leu)
        assert_eq!(code.translate_codon(b"CTA"), Some(b'T'));
        // UGA=Trp
        assert_eq!(code.translate_codon(b"TGA"), Some(b'W'));
        // AUA=Met
        assert_eq!(code.translate_codon(b"ATA"), Some(b'M'));
    }

    #[test]
    fn mycoplasma_differences() {
        let code = GeneticCode::from_id(GeneticCodeId::MycoplasmaSpiroplasma);
        // UGA=Trp (not stop)
        assert_eq!(code.translate_codon(b"TGA"), Some(b'W'));
        // Standard stops still apply
        assert_eq!(code.translate_codon(b"TAA"), None);
        assert_eq!(code.translate_codon(b"TAG"), None);
    }

    #[test]
    fn invertebrate_mito_differences() {
        let code = GeneticCode::from_id(GeneticCodeId::InvertebrateMitochondrial);
        // AGA=Ser (not Arg)
        assert_eq!(code.translate_codon(b"AGA"), Some(b'S'));
        // AGG=Ser (not Arg)
        assert_eq!(code.translate_codon(b"AGG"), Some(b'S'));
        // UGA=Trp
        assert_eq!(code.translate_codon(b"TGA"), Some(b'W'));
    }

    #[test]
    fn ciliate_differences() {
        let code = GeneticCode::from_id(GeneticCodeId::CiliateNuclear);
        // UAA=Gln (not stop)
        assert_eq!(code.translate_codon(b"TAA"), Some(b'Q'));
        // UAG=Gln (not stop)
        assert_eq!(code.translate_codon(b"TAG"), Some(b'Q'));
        // UGA still stop
        assert_eq!(code.translate_codon(b"TGA"), None);
    }

    #[test]
    fn bacterial_plastid_same_aas_different_starts() {
        let code = GeneticCode::from_id(GeneticCodeId::BacterialPlastid);
        // Same AAs as standard
        assert_eq!(code.translate_codon(b"ATG"), Some(b'M'));
        assert_eq!(code.translate_codon(b"TGA"), None);
        // GTG and TTG are starts
        assert!(code.is_start(b"GTG"));
        assert!(code.is_start(b"TTG"));
        assert!(code.is_start(b"ATG"));
    }

    #[test]
    fn start_stop_codons() {
        let code = GeneticCode::standard();
        assert!(code.is_start(b"ATG"));
        assert!(!code.is_start(b"GTG"));
        assert!(code.is_stop(b"TAA"));
        assert!(code.is_stop(b"TAG"));
        assert!(code.is_stop(b"TGA"));
        assert!(!code.is_stop(b"ATG"));

        let stops = code.stop_codons();
        assert_eq!(stops.len(), 3);
        let starts = code.start_codons();
        assert_eq!(starts.len(), 1);
    }

    #[test]
    fn codon_usage_basic() {
        // ATGATG = two ATG codons
        let usage = CodonUsage::from_sequence(b"ATGATG");
        assert_eq!(usage.count(b"ATG"), 2);
        assert_eq!(usage.total, 2);
        assert!((usage.frequency(b"ATG") - 1.0).abs() < 1e-10);
    }

    #[test]
    fn codon_usage_merge() {
        let mut u1 = CodonUsage::from_sequence(b"ATGATG");
        let u2 = CodonUsage::from_sequence(b"ATGAAA");
        u1.merge(&u2);
        assert_eq!(u1.count(b"ATG"), 3);
        assert_eq!(u1.count(b"AAA"), 1);
        assert_eq!(u1.total, 4);
    }

    #[test]
    fn rscu_single_codon_family() {
        let code = GeneticCode::standard();
        // ATG is the only Met codon → RSCU = 1.0
        let usage = CodonUsage::from_sequence(b"ATGATGATG");
        let rscu = usage.rscu(b"ATG", &code);
        assert!((rscu - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rscu_two_fold_degenerate() {
        let code = GeneticCode::standard();
        // Phe: TTC and TTT. If all are TTT, RSCU(TTT)=2.0, RSCU(TTC)=0.0
        let usage = CodonUsage::from_sequence(b"TTTTTTTTT");
        let rscu_ttt = usage.rscu(b"TTT", &code);
        assert!((rscu_ttt - 2.0).abs() < 1e-10);
        let rscu_ttc = usage.rscu(b"TTC", &code);
        assert!(rscu_ttc.abs() < 1e-10);
    }

    #[test]
    fn cai_high_adaptation() {
        let code = GeneticCode::standard();
        // Reference and query use same codons → CAI ≈ 1.0
        let ref_seq = b"ATGTTTAAA";
        let reference = CodonUsage::from_sequence(ref_seq);
        let cai = codon_adaptation_index(ref_seq, &reference, &code).unwrap();
        assert!((cai - 1.0).abs() < 0.01, "CAI={}", cai);
    }

    #[test]
    fn cai_error_on_empty() {
        let code = GeneticCode::standard();
        let reference = CodonUsage::from_sequence(b"ATGATG");
        let result = codon_adaptation_index(b"", &reference, &code);
        assert!(result.is_err());
    }

    #[test]
    fn classify_synonymous() {
        let code = GeneticCode::standard();
        // TTT→TTC: both Phe
        assert_eq!(
            classify_substitution(b"TTT", b"TTC", &code),
            SubstitutionClass::Synonymous
        );
    }

    #[test]
    fn classify_nonsynonymous() {
        let code = GeneticCode::standard();
        // ATG(Met)→GTG(Val)
        assert_eq!(
            classify_substitution(b"ATG", b"GTG", &code),
            SubstitutionClass::NonSynonymous
        );
    }

    #[test]
    fn classify_stop_involved() {
        let code = GeneticCode::standard();
        assert_eq!(
            classify_substitution(b"ATG", b"TAG", &code),
            SubstitutionClass::StopInvolved
        );
    }

    #[test]
    fn count_sites_basic() {
        let code = GeneticCode::standard();
        // ATG has 0 synonymous sites (only Met codon), 3 non-synonymous
        let (s, n) = count_syn_nonsyn_sites(b"ATG", &code);
        assert!(s < 0.5, "ATG syn sites should be ~0, got {}", s);
        assert!(n > 2.5, "ATG nonsyn sites should be ~3, got {}", n);
    }

    #[test]
    fn count_sites_four_fold() {
        let code = GeneticCode::standard();
        // GCN (Ala) — 4-fold degenerate at 3rd position
        // Third position: 3 of 3 possible changes are synonymous → 1 syn site
        let (s, _n) = count_syn_nonsyn_sites(b"GCA", &code);
        assert!(s > 0.5, "Ala codon should have >0.5 syn sites, got {}", s);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn coding_seq(max_codons: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(
            prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
            3..=(max_codons * 3),
        )
        .prop_map(|v| {
            let len = v.len() - (v.len() % 3);
            v[..len].to_vec()
        })
    }

    proptest! {
        #[test]
        fn cai_in_unit_interval(seq in coding_seq(20)) {
            let code = GeneticCode::standard();
            let reference = CodonUsage::from_sequence(&seq);
            if let Ok(cai) = codon_adaptation_index(&seq, &reference, &code) {
                prop_assert!(cai > 0.0 && cai <= 1.0 + 1e-10,
                    "CAI should be in (0,1], got {}", cai);
            }
        }
    }
}
