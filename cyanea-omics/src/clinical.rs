//! Clinical genomics: HLA typing, tumor mutational burden, microsatellite instability.
//!
//! Provides simplified implementations of key clinical genomics metrics used
//! in oncology and transplant medicine.

use crate::variant::{Variant, VariantType};
use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

// ──────────────────────────── HLA Typing ────────────────────────────

/// An HLA allele (e.g., "A*02:01:01").
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HlaAllele {
    /// Gene (e.g., "A", "B", "C", "DRB1").
    pub gene: String,
    /// Full allele designation (e.g., "02:01:01").
    pub allele: String,
}

impl HlaAllele {
    pub fn new(gene: &str, allele: &str) -> Self {
        Self {
            gene: gene.to_string(),
            allele: allele.to_string(),
        }
    }

    /// Two-digit resolution (e.g., "A*02").
    pub fn two_digit(&self) -> String {
        let first_field = self.allele.split(':').next().unwrap_or(&self.allele);
        format!("{}*{}", self.gene, first_field)
    }

    /// Four-digit resolution (e.g., "A*02:01").
    pub fn four_digit(&self) -> String {
        let fields: Vec<&str> = self.allele.split(':').collect();
        if fields.len() >= 2 {
            format!("{}*{}:{}", self.gene, fields[0], fields[1])
        } else {
            format!("{}*{}", self.gene, self.allele)
        }
    }

    /// Full notation (e.g., "HLA-A*02:01:01").
    pub fn full_name(&self) -> String {
        format!("HLA-{}*{}", self.gene, self.allele)
    }
}

/// HLA typing result for an individual.
#[derive(Debug, Clone)]
pub struct HlaTypingResult {
    /// Gene → (allele1, allele2) for each typed locus.
    pub genotypes: HashMap<String, (HlaAllele, HlaAllele)>,
}

impl HlaTypingResult {
    /// Get all typed genes.
    pub fn genes(&self) -> Vec<&str> {
        self.genotypes.keys().map(|s| s.as_str()).collect()
    }

    /// Get the diplotype for a gene.
    pub fn diplotype(&self, gene: &str) -> Option<(&HlaAllele, &HlaAllele)> {
        self.genotypes.get(gene).map(|(a, b)| (a, b))
    }

    /// Check if two individuals share any alleles at a locus.
    pub fn shared_alleles(&self, other: &HlaTypingResult, gene: &str) -> usize {
        let (a1, a2) = match self.genotypes.get(gene) {
            Some(pair) => pair,
            None => return 0,
        };
        let (b1, b2) = match other.genotypes.get(gene) {
            Some(pair) => pair,
            None => return 0,
        };

        let mut count = 0;
        if a1.allele == b1.allele || a1.allele == b2.allele {
            count += 1;
        }
        if a2.allele == b1.allele || a2.allele == b2.allele {
            count += 1;
        }
        count
    }
}

/// Compute HLA compatibility between donor and recipient.
///
/// Returns the number of matched alleles across the specified loci (typically
/// HLA-A, -B, -C, -DRB1 for transplant matching). Maximum = 2 × num_loci.
pub fn hla_compatibility(
    donor: &HlaTypingResult,
    recipient: &HlaTypingResult,
    loci: &[&str],
) -> usize {
    loci.iter()
        .map(|gene| donor.shared_alleles(recipient, gene))
        .sum()
}

/// Parse HLA typing from a simple text format.
///
/// Each line: `GENE\tALLELE1\tALLELE2`
/// e.g., `A\t02:01\t03:01`
pub fn parse_hla_typing(content: &str) -> Result<HlaTypingResult> {
    let mut genotypes = HashMap::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let gene = fields[0];
        let allele1 = HlaAllele::new(gene, fields[1]);
        let allele2 = HlaAllele::new(gene, fields[2]);
        genotypes.insert(gene.to_string(), (allele1, allele2));
    }

    Ok(HlaTypingResult { genotypes })
}

// ──────────────────────────── TMB ────────────────────────────

/// Tumor mutational burden (TMB) result.
#[derive(Debug, Clone)]
pub struct TmbResult {
    /// Total somatic mutations counted.
    pub mutation_count: usize,
    /// Exome size in megabases.
    pub exome_size_mb: f64,
    /// TMB in mutations per megabase.
    pub tmb: f64,
    /// TMB category.
    pub category: TmbCategory,
    /// Breakdown by variant type.
    pub by_type: HashMap<String, usize>,
}

/// TMB category following common clinical thresholds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TmbCategory {
    /// < 6 mut/Mb.
    Low,
    /// 6–19 mut/Mb.
    Intermediate,
    /// ≥ 20 mut/Mb (FDA threshold for pembrolizumab).
    High,
}

impl TmbCategory {
    pub fn as_str(&self) -> &'static str {
        match self {
            TmbCategory::Low => "TMB-Low",
            TmbCategory::Intermediate => "TMB-Intermediate",
            TmbCategory::High => "TMB-High",
        }
    }
}

/// Compute tumor mutational burden from somatic variants.
///
/// # Arguments
///
/// * `variants` - Somatic variant calls (should be filtered for quality)
/// * `exome_size_mb` - Size of the assessed coding region in megabases (typically ~30–50 Mb)
/// * `count_indels` - Whether to include indels (some panels count only SNVs)
pub fn compute_tmb(
    variants: &[Variant],
    exome_size_mb: f64,
    count_indels: bool,
) -> Result<TmbResult> {
    if exome_size_mb <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "exome size must be positive".into(),
        ));
    }

    let mut count = 0usize;
    let mut by_type: HashMap<String, usize> = HashMap::new();

    for v in variants {
        let vtype = v.variant_type();
        let include = match vtype {
            VariantType::Snv | VariantType::Mnv => true,
            VariantType::Insertion | VariantType::Deletion => count_indels,
            VariantType::Complex => count_indels,
        };

        if include {
            count += 1;
            let type_name = format!("{:?}", vtype);
            *by_type.entry(type_name).or_insert(0) += 1;
        }
    }

    let tmb = count as f64 / exome_size_mb;
    let category = if tmb >= 20.0 {
        TmbCategory::High
    } else if tmb >= 6.0 {
        TmbCategory::Intermediate
    } else {
        TmbCategory::Low
    };

    Ok(TmbResult {
        mutation_count: count,
        exome_size_mb,
        tmb,
        category,
        by_type,
    })
}

// ──────────────────────────── MSI ────────────────────────────

/// A microsatellite locus for MSI testing.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MsiLocus {
    /// Locus name (e.g., "BAT25", "BAT26", "NR21", "NR24", "MONO27").
    pub name: String,
    /// Chromosome.
    pub chrom: String,
    /// Start position.
    pub start: u64,
    /// End position.
    pub end: u64,
    /// Repeat unit (e.g., "A", "CA").
    pub repeat_unit: String,
    /// Reference repeat count.
    pub reference_count: usize,
}

/// MSI status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MsiStatus {
    /// Microsatellite stable (0 unstable loci).
    MSS,
    /// MSI-Low (1 unstable locus out of 5 Bethesda markers).
    MSILow,
    /// MSI-High (≥2 unstable loci out of 5 Bethesda markers, or ≥30% unstable).
    MSIHigh,
}

impl MsiStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            MsiStatus::MSS => "MSS",
            MsiStatus::MSILow => "MSI-L",
            MsiStatus::MSIHigh => "MSI-H",
        }
    }
}

/// MSI analysis result.
#[derive(Debug, Clone)]
pub struct MsiResult {
    /// Overall MSI status.
    pub status: MsiStatus,
    /// Total loci tested.
    pub total_loci: usize,
    /// Number of unstable loci.
    pub unstable_loci: usize,
    /// Fraction of unstable loci.
    pub instability_fraction: f64,
    /// Per-locus instability calls.
    pub locus_calls: Vec<(String, bool)>,
}

/// Standard Bethesda markers for MSI testing.
pub fn bethesda_markers() -> Vec<MsiLocus> {
    vec![
        MsiLocus {
            name: "BAT25".into(),
            chrom: "chr4".into(),
            start: 55598212,
            end: 55598236,
            repeat_unit: "A".into(),
            reference_count: 25,
        },
        MsiLocus {
            name: "BAT26".into(),
            chrom: "chr2".into(),
            start: 47641560,
            end: 47641586,
            repeat_unit: "A".into(),
            reference_count: 26,
        },
        MsiLocus {
            name: "NR21".into(),
            chrom: "chr14".into(),
            start: 23652347,
            end: 23652367,
            repeat_unit: "A".into(),
            reference_count: 21,
        },
        MsiLocus {
            name: "NR24".into(),
            chrom: "chr2".into(),
            start: 95849362,
            end: 95849385,
            repeat_unit: "A".into(),
            reference_count: 24,
        },
        MsiLocus {
            name: "MONO27".into(),
            chrom: "chr4".into(),
            start: 2866492,
            end: 2866519,
            repeat_unit: "A".into(),
            reference_count: 27,
        },
    ]
}

/// Determine MSI status from observed repeat counts at marker loci.
///
/// A locus is called unstable if its observed repeat count differs from
/// the reference by more than `shift_threshold` repeats.
///
/// # Arguments
///
/// * `observed_counts` - Locus name → observed repeat count in tumor
/// * `markers` - Microsatellite loci to test
/// * `shift_threshold` - Minimum repeat count shift to call instability (typically 2–3)
pub fn call_msi(
    observed_counts: &HashMap<String, usize>,
    markers: &[MsiLocus],
    shift_threshold: usize,
) -> MsiResult {
    let mut unstable = 0usize;
    let mut locus_calls = Vec::new();
    let total = markers.len();

    for locus in markers {
        let is_unstable = if let Some(&observed) = observed_counts.get(&locus.name) {
            observed.abs_diff(locus.reference_count) > shift_threshold
        } else {
            false // no data → assume stable
        };

        locus_calls.push((locus.name.clone(), is_unstable));
        if is_unstable {
            unstable += 1;
        }
    }

    let fraction = if total > 0 {
        unstable as f64 / total as f64
    } else {
        0.0
    };

    let status = if unstable >= 2 || fraction >= 0.3 {
        MsiStatus::MSIHigh
    } else if unstable == 1 {
        MsiStatus::MSILow
    } else {
        MsiStatus::MSS
    };

    MsiResult {
        status,
        total_loci: total,
        unstable_loci: unstable,
        instability_fraction: fraction,
        locus_calls,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn make_variant(chrom: &str, pos: u64, ref_a: &[u8], alt_a: &[u8]) -> Variant {
        Variant::new(chrom, pos, ref_a.to_vec(), vec![alt_a.to_vec()]).unwrap()
    }

    // ── HLA ──

    #[test]
    fn test_hla_allele_resolution() {
        let a = HlaAllele::new("A", "02:01:01:01");
        assert_eq!(a.two_digit(), "A*02");
        assert_eq!(a.four_digit(), "A*02:01");
        assert_eq!(a.full_name(), "HLA-A*02:01:01:01");
    }

    #[test]
    fn test_hla_compatibility_perfect() {
        let typing = HlaTypingResult {
            genotypes: vec![
                ("A".into(), (HlaAllele::new("A", "02:01"), HlaAllele::new("A", "03:01"))),
                ("B".into(), (HlaAllele::new("B", "07:02"), HlaAllele::new("B", "44:02"))),
            ]
            .into_iter()
            .collect(),
        };
        let match_count = hla_compatibility(&typing, &typing, &["A", "B"]);
        assert_eq!(match_count, 4); // 2 per locus × 2 loci
    }

    #[test]
    fn test_hla_compatibility_mismatch() {
        let donor = HlaTypingResult {
            genotypes: vec![(
                "A".into(),
                (HlaAllele::new("A", "02:01"), HlaAllele::new("A", "03:01")),
            )]
            .into_iter()
            .collect(),
        };
        let recipient = HlaTypingResult {
            genotypes: vec![(
                "A".into(),
                (HlaAllele::new("A", "24:02"), HlaAllele::new("A", "68:01")),
            )]
            .into_iter()
            .collect(),
        };
        let match_count = hla_compatibility(&donor, &recipient, &["A"]);
        assert_eq!(match_count, 0);
    }

    #[test]
    fn test_hla_compatibility_partial() {
        let donor = HlaTypingResult {
            genotypes: vec![(
                "A".into(),
                (HlaAllele::new("A", "02:01"), HlaAllele::new("A", "03:01")),
            )]
            .into_iter()
            .collect(),
        };
        let recipient = HlaTypingResult {
            genotypes: vec![(
                "A".into(),
                (HlaAllele::new("A", "02:01"), HlaAllele::new("A", "68:01")),
            )]
            .into_iter()
            .collect(),
        };
        let match_count = hla_compatibility(&donor, &recipient, &["A"]);
        assert_eq!(match_count, 1);
    }

    #[test]
    fn test_parse_hla_typing() {
        let content = "A\t02:01\t03:01\nB\t07:02\t44:02\nC\t07:01\t05:01\n";
        let result = parse_hla_typing(content).unwrap();
        assert_eq!(result.genes().len(), 3);
        let (a1, a2) = result.diplotype("A").unwrap();
        assert_eq!(a1.allele, "02:01");
        assert_eq!(a2.allele, "03:01");
    }

    // ── TMB ──

    #[test]
    fn test_tmb_low() {
        let variants: Vec<Variant> = (0..50)
            .map(|i| make_variant("chr1", i * 1000, b"A", b"G"))
            .collect();
        let result = compute_tmb(&variants, 30.0, false).unwrap();
        assert!((result.tmb - 50.0 / 30.0).abs() < 0.01);
        assert_eq!(result.category, TmbCategory::Low);
    }

    #[test]
    fn test_tmb_high() {
        let variants: Vec<Variant> = (0..1000)
            .map(|i| make_variant("chr1", i * 100, b"A", b"G"))
            .collect();
        let result = compute_tmb(&variants, 30.0, false).unwrap();
        assert!(result.tmb > 20.0);
        assert_eq!(result.category, TmbCategory::High);
    }

    #[test]
    fn test_tmb_with_indels() {
        let mut variants: Vec<Variant> = (0..10)
            .map(|i| make_variant("chr1", i * 1000, b"A", b"G"))
            .collect();
        // Add indels
        variants.push(make_variant("chr1", 50000, b"AC", b"A"));
        variants.push(make_variant("chr1", 60000, b"A", b"AC"));

        let without_indels = compute_tmb(&variants, 1.0, false).unwrap();
        let with_indels = compute_tmb(&variants, 1.0, true).unwrap();
        assert_eq!(without_indels.mutation_count, 10);
        assert_eq!(with_indels.mutation_count, 12);
    }

    #[test]
    fn test_tmb_invalid_size() {
        assert!(compute_tmb(&[], 0.0, false).is_err());
    }

    #[test]
    fn test_tmb_empty() {
        let result = compute_tmb(&[], 30.0, false).unwrap();
        assert_eq!(result.tmb, 0.0);
        assert_eq!(result.category, TmbCategory::Low);
    }

    // ── MSI ──

    #[test]
    fn test_msi_stable() {
        let markers = bethesda_markers();
        let observed: HashMap<String, usize> = markers
            .iter()
            .map(|m| (m.name.clone(), m.reference_count))
            .collect();
        let result = call_msi(&observed, &markers, 2);
        assert_eq!(result.status, MsiStatus::MSS);
        assert_eq!(result.unstable_loci, 0);
    }

    #[test]
    fn test_msi_high() {
        let markers = bethesda_markers();
        let mut observed: HashMap<String, usize> = markers
            .iter()
            .map(|m| (m.name.clone(), m.reference_count))
            .collect();
        // Make BAT25 and BAT26 unstable
        *observed.get_mut("BAT25").unwrap() = 20; // shifted by 5
        *observed.get_mut("BAT26").unwrap() = 20; // shifted by 6
        let result = call_msi(&observed, &markers, 2);
        assert_eq!(result.status, MsiStatus::MSIHigh);
        assert_eq!(result.unstable_loci, 2);
    }

    #[test]
    fn test_msi_low() {
        let markers = bethesda_markers();
        let mut observed: HashMap<String, usize> = markers
            .iter()
            .map(|m| (m.name.clone(), m.reference_count))
            .collect();
        *observed.get_mut("BAT25").unwrap() = 20; // shifted by 5
        let result = call_msi(&observed, &markers, 2);
        assert_eq!(result.status, MsiStatus::MSILow);
        assert_eq!(result.unstable_loci, 1);
    }

    #[test]
    fn test_bethesda_markers() {
        let markers = bethesda_markers();
        assert_eq!(markers.len(), 5);
        assert_eq!(markers[0].name, "BAT25");
    }

    #[test]
    fn test_msi_no_data() {
        let markers = bethesda_markers();
        let observed = HashMap::new();
        let result = call_msi(&observed, &markers, 2);
        assert_eq!(result.status, MsiStatus::MSS);
    }
}
