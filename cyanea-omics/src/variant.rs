//! Genomic variant representation (VCF-style).
//!
//! Types for representing single-nucleotide variants, insertions, deletions,
//! and complex variants with quality and filter information.

use cyanea_core::{Annotated, CyaneaError, Result, Scored};

use crate::genomic::{GenomicInterval, Strand};

/// The class of a genomic variant.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum VariantType {
    /// Single nucleotide variant (ref and alt are both 1 base).
    Snv,
    /// Insertion (alt is longer than ref).
    Insertion,
    /// Deletion (ref is longer than alt).
    Deletion,
    /// Multi-nucleotide variant (ref and alt are equal length > 1).
    Mnv,
    /// Complex variant (none of the above).
    Complex,
}

/// Filter status for a variant call.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum VariantFilter {
    /// Passed all filters.
    Pass,
    /// Failed one or more filters.
    Fail(Vec<String>),
    /// Filter status not available.
    Missing,
}

/// Zygosity of a variant call.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Zygosity {
    Homozygous,
    Heterozygous,
    Hemizygous,
    Unknown,
}

/// A genomic variant (VCF-style representation).
///
/// Position is 1-based following VCF convention.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Variant {
    pub chrom: String,
    /// 1-based position (VCF convention).
    pub position: u64,
    /// Optional identifier (e.g. rs12345).
    pub id: Option<String>,
    /// Reference allele.
    pub ref_allele: Vec<u8>,
    /// Alternate alleles.
    pub alt_alleles: Vec<Vec<u8>>,
    /// Phred-scaled quality score.
    pub quality: Option<f64>,
    /// Filter status.
    pub filter: VariantFilter,
}

impl Variant {
    /// Create a new variant with minimal fields.
    ///
    /// Validates that reference and alternate alleles are non-empty.
    pub fn new(
        chrom: impl Into<String>,
        position: u64,
        ref_allele: Vec<u8>,
        alt_alleles: Vec<Vec<u8>>,
    ) -> Result<Self> {
        if ref_allele.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "reference allele must not be empty".into(),
            ));
        }
        if alt_alleles.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "at least one alternate allele is required".into(),
            ));
        }
        for (i, alt) in alt_alleles.iter().enumerate() {
            if alt.is_empty() {
                return Err(CyaneaError::InvalidInput(format!(
                    "alternate allele {i} must not be empty"
                )));
            }
        }
        Ok(Self {
            chrom: chrom.into(),
            position,
            id: None,
            ref_allele,
            alt_alleles,
            quality: None,
            filter: VariantFilter::Missing,
        })
    }

    /// Infer the variant type from the first alternate allele.
    pub fn variant_type(&self) -> VariantType {
        let ref_len = self.ref_allele.len();
        let alt_len = self.alt_alleles[0].len();

        if ref_len == 1 && alt_len == 1 {
            VariantType::Snv
        } else if ref_len == alt_len {
            VariantType::Mnv
        } else if ref_len < alt_len {
            VariantType::Insertion
        } else if ref_len > alt_len {
            VariantType::Deletion
        } else {
            VariantType::Complex
        }
    }

    /// Whether this is a single-nucleotide variant.
    pub fn is_snv(&self) -> bool {
        self.variant_type() == VariantType::Snv
    }

    /// Whether this is an insertion or deletion.
    pub fn is_indel(&self) -> bool {
        matches!(
            self.variant_type(),
            VariantType::Insertion | VariantType::Deletion
        )
    }

    /// Whether this SNV is a transition (A↔G or C↔T).
    ///
    /// Returns `false` for non-SNV variants.
    pub fn is_transition(&self) -> bool {
        if !self.is_snv() {
            return false;
        }
        let r = self.ref_allele[0].to_ascii_uppercase();
        let a = self.alt_alleles[0][0].to_ascii_uppercase();
        matches!(
            (r, a),
            (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
        )
    }

    /// Whether this SNV is a transversion (complement of transition for SNVs).
    ///
    /// Returns `false` for non-SNV variants.
    pub fn is_transversion(&self) -> bool {
        self.is_snv() && !self.is_transition()
    }

    /// Convert the 1-based VCF position to a 0-based [`GenomicInterval`].
    ///
    /// The interval spans the reference allele.
    pub fn to_genomic_interval(&self) -> GenomicInterval {
        let start = self.position - 1; // VCF is 1-based
        let end = start + self.ref_allele.len() as u64;
        // This is safe because ref_allele is validated non-empty, so start < end.
        GenomicInterval {
            chrom: self.chrom.clone(),
            start,
            end,
            strand: Strand::Unknown,
        }
    }
}

impl Annotated for Variant {
    fn name(&self) -> &str {
        // Can't return a computed string from &str, so use id if available.
        // For the fallback we'd need to store a formatted name; return "" instead.
        match &self.id {
            Some(id) => id.as_str(),
            None => "",
        }
    }
}

impl Scored for Variant {
    fn score(&self) -> f64 {
        self.quality.unwrap_or(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snv_a_to_g() -> Variant {
        Variant::new("chr1", 100, vec![b'A'], vec![vec![b'G']]).unwrap()
    }

    #[test]
    fn test_snv_construction() {
        let v = snv_a_to_g();
        assert_eq!(v.chrom, "chr1");
        assert_eq!(v.position, 100);
        assert!(v.is_snv());
    }

    #[test]
    fn test_empty_ref_allele() {
        assert!(Variant::new("chr1", 100, vec![], vec![vec![b'G']]).is_err());
    }

    #[test]
    fn test_empty_alt_alleles() {
        assert!(Variant::new("chr1", 100, vec![b'A'], vec![]).is_err());
    }

    #[test]
    fn test_empty_alt_allele_entry() {
        assert!(Variant::new("chr1", 100, vec![b'A'], vec![vec![]]).is_err());
    }

    #[test]
    fn test_variant_type_snv() {
        let v = snv_a_to_g();
        assert_eq!(v.variant_type(), VariantType::Snv);
    }

    #[test]
    fn test_variant_type_insertion() {
        let v = Variant::new("chr1", 100, vec![b'A'], vec![vec![b'A', b'T', b'C']]).unwrap();
        assert_eq!(v.variant_type(), VariantType::Insertion);
        assert!(v.is_indel());
    }

    #[test]
    fn test_variant_type_deletion() {
        let v = Variant::new("chr1", 100, vec![b'A', b'T', b'C'], vec![vec![b'A']]).unwrap();
        assert_eq!(v.variant_type(), VariantType::Deletion);
        assert!(v.is_indel());
    }

    #[test]
    fn test_variant_type_mnv() {
        let v = Variant::new("chr1", 100, vec![b'A', b'T'], vec![vec![b'G', b'C']]).unwrap();
        assert_eq!(v.variant_type(), VariantType::Mnv);
    }

    #[test]
    fn test_transition() {
        // A→G is a transition
        assert!(snv_a_to_g().is_transition());
        assert!(!snv_a_to_g().is_transversion());

        // C→T is a transition
        let v = Variant::new("chr1", 100, vec![b'C'], vec![vec![b'T']]).unwrap();
        assert!(v.is_transition());
    }

    #[test]
    fn test_transversion() {
        // A→C is a transversion
        let v = Variant::new("chr1", 100, vec![b'A'], vec![vec![b'C']]).unwrap();
        assert!(v.is_transversion());
        assert!(!v.is_transition());
    }

    #[test]
    fn test_transition_non_snv() {
        let v = Variant::new("chr1", 100, vec![b'A', b'T'], vec![vec![b'G', b'C']]).unwrap();
        assert!(!v.is_transition());
        assert!(!v.is_transversion());
    }

    #[test]
    fn test_to_genomic_interval() {
        let v = snv_a_to_g();
        let iv = v.to_genomic_interval();
        assert_eq!(iv.chrom, "chr1");
        assert_eq!(iv.start, 99); // 1-based 100 → 0-based 99
        assert_eq!(iv.end, 100);
    }

    #[test]
    fn test_to_genomic_interval_deletion() {
        let v = Variant::new("chr1", 100, vec![b'A', b'T', b'C'], vec![vec![b'A']]).unwrap();
        let iv = v.to_genomic_interval();
        assert_eq!(iv.start, 99);
        assert_eq!(iv.end, 102); // 3-base ref
    }

    #[test]
    fn test_annotated() {
        let mut v = snv_a_to_g();
        assert_eq!(v.name(), "");

        v.id = Some("rs12345".into());
        assert_eq!(v.name(), "rs12345");
    }

    #[test]
    fn test_scored() {
        let mut v = snv_a_to_g();
        assert_eq!(v.score(), 0.0);

        v.quality = Some(30.0);
        assert_eq!(v.score(), 30.0);
    }
}
