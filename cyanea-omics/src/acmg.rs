//! ACMG/AMP variant classification.
//!
//! Implements the American College of Medical Genetics (ACMG) and Association
//! for Molecular Pathology (AMP) 2015 guidelines for interpretation of
//! sequence variants. Classifies variants as Pathogenic, Likely Pathogenic,
//! VUS, Likely Benign, or Benign based on weighted evidence criteria.

use crate::variant::{Variant, VariantType};
use cyanea_core::Result;

/// ACMG variant classification (5-tier).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AcmgClass {
    Benign,
    LikelyBenign,
    Vus,
    LikelyPathogenic,
    Pathogenic,
}

impl AcmgClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            AcmgClass::Benign => "Benign",
            AcmgClass::LikelyBenign => "Likely Benign",
            AcmgClass::Vus => "VUS",
            AcmgClass::LikelyPathogenic => "Likely Pathogenic",
            AcmgClass::Pathogenic => "Pathogenic",
        }
    }
}

/// Strength of evidence for ACMG criteria.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EvidenceStrength {
    /// Stand-alone (PVS1).
    VeryStrong,
    /// Strong (PS1–PS4).
    Strong,
    /// Moderate (PM1–PM6).
    Moderate,
    /// Supporting (PP1–PP5, BP1–BP7).
    Supporting,
}

/// An ACMG evidence criterion applied to a variant.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AcmgCriterion {
    /// Criterion code (e.g., "PVS1", "PS1", "PM2", "PP3", "BA1", "BS1", "BP4").
    pub code: String,
    /// Whether this is pathogenic (true) or benign (false) evidence.
    pub is_pathogenic: bool,
    /// Evidence strength.
    pub strength: EvidenceStrength,
    /// Brief description of why this criterion applies.
    pub description: String,
}

impl AcmgCriterion {
    pub fn new(code: &str, is_pathogenic: bool, strength: EvidenceStrength, desc: &str) -> Self {
        Self {
            code: code.to_string(),
            is_pathogenic,
            strength,
            description: desc.to_string(),
        }
    }
}

/// Result of ACMG classification for a single variant.
#[derive(Debug, Clone)]
pub struct AcmgClassification {
    /// Final classification.
    pub classification: AcmgClass,
    /// Evidence criteria that were applied.
    pub criteria: Vec<AcmgCriterion>,
    /// Number of pathogenic very-strong criteria.
    pub pvs_count: usize,
    /// Number of pathogenic strong criteria.
    pub ps_count: usize,
    /// Number of pathogenic moderate criteria.
    pub pm_count: usize,
    /// Number of pathogenic supporting criteria.
    pub pp_count: usize,
    /// Number of benign stand-alone criteria.
    pub ba_count: usize,
    /// Number of benign strong criteria.
    pub bs_count: usize,
    /// Number of benign supporting criteria.
    pub bp_count: usize,
}

/// Evidence collected for a variant prior to classification.
#[derive(Debug, Clone, Default)]
pub struct AcmgEvidence {
    /// Applied criteria.
    pub criteria: Vec<AcmgCriterion>,
}

impl AcmgEvidence {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a pathogenic very-strong criterion (PVS1).
    pub fn pvs1(mut self, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            "PVS1", true, EvidenceStrength::VeryStrong, desc,
        ));
        self
    }

    /// Add a pathogenic strong criterion.
    pub fn strong_pathogenic(mut self, code: &str, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            code, true, EvidenceStrength::Strong, desc,
        ));
        self
    }

    /// Add a pathogenic moderate criterion.
    pub fn moderate_pathogenic(mut self, code: &str, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            code, true, EvidenceStrength::Moderate, desc,
        ));
        self
    }

    /// Add a pathogenic supporting criterion.
    pub fn supporting_pathogenic(mut self, code: &str, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            code, true, EvidenceStrength::Supporting, desc,
        ));
        self
    }

    /// Add a benign stand-alone criterion (BA1).
    pub fn ba1(mut self, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            "BA1", false, EvidenceStrength::VeryStrong, desc,
        ));
        self
    }

    /// Add a benign strong criterion.
    pub fn strong_benign(mut self, code: &str, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            code, false, EvidenceStrength::Strong, desc,
        ));
        self
    }

    /// Add a benign supporting criterion.
    pub fn supporting_benign(mut self, code: &str, desc: &str) -> Self {
        self.criteria.push(AcmgCriterion::new(
            code, false, EvidenceStrength::Supporting, desc,
        ));
        self
    }

    /// Classify the variant based on collected evidence.
    ///
    /// Implements the ACMG/AMP 2015 combining rules (Richards et al., 2015):
    ///
    /// **Pathogenic**: PVS1 + ≥1 PS, or PVS1 + ≥2 PM, or PVS1 + 1 PM + 1 PP,
    /// or PVS1 + ≥2 PP, or ≥2 PS, or 1 PS + ≥3 PM, or 1 PS + 2 PM + ≥2 PP,
    /// or 1 PS + 1 PM + ≥4 PP.
    ///
    /// **Likely Pathogenic**: 1 PVS1 + 1 PM, or 1 PS + 1-2 PM,
    /// or 1 PS + ≥2 PP, or ≥3 PM, or 2 PM + ≥2 PP, or 1 PM + ≥4 PP.
    ///
    /// **Benign**: 1 BA1, or ≥2 BS.
    ///
    /// **Likely Benign**: 1 BS + 1 BP, or ≥2 BP.
    pub fn classify(&self) -> AcmgClassification {
        let mut pvs = 0usize;
        let mut ps = 0usize;
        let mut pm = 0usize;
        let mut pp = 0usize;
        let mut ba = 0usize;
        let mut bs = 0usize;
        let mut bp = 0usize;

        for criterion in &self.criteria {
            if criterion.is_pathogenic {
                match criterion.strength {
                    EvidenceStrength::VeryStrong => pvs += 1,
                    EvidenceStrength::Strong => ps += 1,
                    EvidenceStrength::Moderate => pm += 1,
                    EvidenceStrength::Supporting => pp += 1,
                }
            } else {
                match criterion.strength {
                    EvidenceStrength::VeryStrong => ba += 1,
                    EvidenceStrength::Strong => bs += 1,
                    EvidenceStrength::Moderate => {} // no benign moderate in standard rules
                    EvidenceStrength::Supporting => bp += 1,
                }
            }
        }

        let classification = classify_from_counts(pvs, ps, pm, pp, ba, bs, bp);

        AcmgClassification {
            classification,
            criteria: self.criteria.clone(),
            pvs_count: pvs,
            ps_count: ps,
            pm_count: pm,
            pp_count: pp,
            ba_count: ba,
            bs_count: bs,
            bp_count: bp,
        }
    }
}

/// Apply the ACMG combining rules from evidence counts.
fn classify_from_counts(
    pvs: usize, ps: usize, pm: usize, pp: usize,
    ba: usize, bs: usize, bp: usize,
) -> AcmgClass {
    // Benign rules (checked first — BA1 is stand-alone)
    if ba >= 1 {
        return AcmgClass::Benign;
    }
    if bs >= 2 {
        return AcmgClass::Benign;
    }

    // Likely Benign
    if (bs >= 1 && bp >= 1) || bp >= 2 {
        return AcmgClass::LikelyBenign;
    }

    // Pathogenic rules
    let pathogenic = (pvs >= 1 && ps >= 1)
        || (pvs >= 1 && pm >= 2)
        || (pvs >= 1 && pm >= 1 && pp >= 1)
        || (pvs >= 1 && pp >= 2)
        || (ps >= 2)
        || (ps >= 1 && pm >= 3)
        || (ps >= 1 && pm >= 2 && pp >= 2)
        || (ps >= 1 && pm >= 1 && pp >= 4);

    if pathogenic {
        return AcmgClass::Pathogenic;
    }

    // Likely Pathogenic
    let likely_path = (pvs >= 1 && pm >= 1)
        || (ps >= 1 && pm >= 1)
        || (ps >= 1 && pp >= 2)
        || (pm >= 3)
        || (pm >= 2 && pp >= 2)
        || (pm >= 1 && pp >= 4);

    if likely_path {
        return AcmgClass::LikelyPathogenic;
    }

    AcmgClass::Vus
}

/// Auto-assign basic ACMG criteria from variant properties.
///
/// This is a simplified auto-classifier. For clinical use, manual review
/// and additional database lookups are required.
///
/// Checks:
/// - PVS1: null variant in gene with known LOF mechanism
/// - PM2: absent from population databases (allele_freq == None or < 0.0001)
/// - PP3: in silico predictions (if provided)
/// - BA1: allele frequency > 5% in population
/// - BS1: allele frequency > expected for disorder
pub fn auto_evidence(
    variant: &Variant,
    allele_freq: Option<f64>,
    is_lof_gene: bool,
    in_silico_pathogenic: Option<bool>,
) -> AcmgEvidence {
    let mut ev = AcmgEvidence::new();

    let is_null = matches!(
        variant.variant_type(),
        VariantType::Insertion | VariantType::Deletion
    ) && variant.ref_allele.len() != variant.alt_alleles[0].len()
        && (variant.ref_allele.len().abs_diff(variant.alt_alleles[0].len()) % 3 != 0);

    // PVS1: null variant in LOF gene
    if is_null && is_lof_gene {
        ev = ev.pvs1("Null variant (frameshift) in gene where LOF is a known mechanism");
    }

    // Population frequency criteria
    if let Some(freq) = allele_freq {
        if freq > 0.05 {
            ev = ev.ba1(&format!("Allele frequency {:.4} > 5% in population", freq));
        } else if freq > 0.01 {
            ev = ev.strong_benign("BS1", &format!("Allele frequency {:.4} > 1%", freq));
        } else if freq < 0.0001 {
            ev = ev.moderate_pathogenic("PM2", &format!("Absent/rare in population (AF={:.6})", freq));
        }
    } else {
        ev = ev.moderate_pathogenic("PM2", "Absent from population databases");
    }

    // In silico predictions
    if let Some(pathogenic) = in_silico_pathogenic {
        if pathogenic {
            ev = ev.supporting_pathogenic("PP3", "In silico tools predict damaging effect");
        } else {
            ev = ev.supporting_benign("BP4", "In silico tools predict benign effect");
        }
    }

    ev
}

/// ClinVar-style clinical significance.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ClinVarAnnotation {
    /// ClinVar variation ID.
    pub variation_id: Option<String>,
    /// Clinical significance string (e.g., "Pathogenic", "Likely benign").
    pub significance: String,
    /// Review status (e.g., "criteria provided, single submitter").
    pub review_status: String,
    /// Associated conditions/diseases.
    pub conditions: Vec<String>,
    /// Number of submitters.
    pub submitter_count: usize,
    /// Star rating (0–4).
    pub star_rating: u8,
}

/// Match a variant against a ClinVar-style annotation database.
///
/// Annotations are keyed by (chrom, position, ref, alt).
pub fn match_clinvar(
    variant: &Variant,
    annotations: &[(String, u64, Vec<u8>, Vec<u8>, ClinVarAnnotation)],
) -> Option<ClinVarAnnotation> {
    annotations.iter().find_map(|(chrom, pos, ref_a, alt_a, ann)| {
        if variant.chrom == *chrom
            && variant.position == *pos
            && variant.ref_allele == *ref_a
            && variant.alt_alleles[0] == *alt_a
        {
            Some(ann.clone())
        } else {
            None
        }
    })
}

/// Parse simple ClinVar TSV format.
///
/// Expected columns: chrom, pos, ref, alt, significance, review_status, conditions, submitter_count, star_rating
pub fn parse_clinvar_tsv(content: &str) -> Result<Vec<(String, u64, Vec<u8>, Vec<u8>, ClinVarAnnotation)>> {
    let mut results = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("chrom") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chrom = fields[0].to_string();
        let pos: u64 = fields[1].parse().unwrap_or(0);
        let ref_a = fields[2].as_bytes().to_vec();
        let alt_a = fields[3].as_bytes().to_vec();
        let significance = fields[4].to_string();
        let review_status = fields[5].to_string();
        let conditions: Vec<String> = fields[6].split(';').map(|s| s.to_string()).collect();
        let submitter_count: usize = fields[7].parse().unwrap_or(0);
        let star_rating: u8 = fields[8].parse().unwrap_or(0);

        results.push((
            chrom,
            pos,
            ref_a,
            alt_a,
            ClinVarAnnotation {
                variation_id: None,
                significance,
                review_status,
                conditions,
                submitter_count,
                star_rating,
            },
        ));
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_variant(chrom: &str, pos: u64, ref_a: &[u8], alt_a: &[u8]) -> Variant {
        Variant::new(chrom, pos, ref_a.to_vec(), vec![alt_a.to_vec()]).unwrap()
    }

    #[test]
    fn test_acmg_class_ordering() {
        assert!(AcmgClass::Benign < AcmgClass::LikelyBenign);
        assert!(AcmgClass::LikelyBenign < AcmgClass::Vus);
        assert!(AcmgClass::Vus < AcmgClass::LikelyPathogenic);
        assert!(AcmgClass::LikelyPathogenic < AcmgClass::Pathogenic);
    }

    #[test]
    fn test_pathogenic_pvs1_ps1() {
        let ev = AcmgEvidence::new()
            .pvs1("Frameshift in BRCA1")
            .strong_pathogenic("PS1", "Same amino acid change as known pathogenic");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Pathogenic);
    }

    #[test]
    fn test_pathogenic_two_strong() {
        let ev = AcmgEvidence::new()
            .strong_pathogenic("PS1", "Same AA change")
            .strong_pathogenic("PS3", "Functional studies");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Pathogenic);
    }

    #[test]
    fn test_likely_pathogenic_pvs1_pm() {
        let ev = AcmgEvidence::new()
            .pvs1("Frameshift")
            .moderate_pathogenic("PM2", "Absent from population");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::LikelyPathogenic);
    }

    #[test]
    fn test_likely_pathogenic_three_moderate() {
        let ev = AcmgEvidence::new()
            .moderate_pathogenic("PM1", "In functional domain")
            .moderate_pathogenic("PM2", "Absent from population")
            .moderate_pathogenic("PM5", "Novel missense at same position");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::LikelyPathogenic);
    }

    #[test]
    fn test_benign_ba1() {
        let ev = AcmgEvidence::new()
            .ba1("AF > 5% in gnomAD");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Benign);
    }

    #[test]
    fn test_benign_two_strong() {
        let ev = AcmgEvidence::new()
            .strong_benign("BS1", "AF > 1%")
            .strong_benign("BS2", "Observed in healthy adults");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Benign);
    }

    #[test]
    fn test_likely_benign() {
        let ev = AcmgEvidence::new()
            .strong_benign("BS1", "AF > 1%")
            .supporting_benign("BP4", "In silico benign");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::LikelyBenign);
    }

    #[test]
    fn test_vus_no_evidence() {
        let ev = AcmgEvidence::new();
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Vus);
    }

    #[test]
    fn test_vus_conflicting() {
        // One pathogenic supporting + one benign supporting = VUS
        let ev = AcmgEvidence::new()
            .supporting_pathogenic("PP3", "In silico pathogenic")
            .supporting_benign("BP1", "Missense in gene where truncation causes disease");
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Vus);
    }

    #[test]
    fn test_auto_evidence_rare_lof() {
        let v = make_variant("chr17", 43092919, b"ACGTG", b"A"); // 4bp deletion → frameshift
        let ev = auto_evidence(&v, None, true, Some(true));
        let result = ev.classify();
        // PVS1 + PM2 + PP3 → likely pathogenic or pathogenic
        assert!(result.classification >= AcmgClass::LikelyPathogenic);
    }

    #[test]
    fn test_auto_evidence_common() {
        let v = make_variant("chr1", 100, b"A", b"G");
        let ev = auto_evidence(&v, Some(0.15), false, Some(false));
        let result = ev.classify();
        assert_eq!(result.classification, AcmgClass::Benign); // BA1
    }

    #[test]
    fn test_clinvar_match() {
        let v = make_variant("chr17", 43092919, b"A", b"G");
        let db = vec![(
            "chr17".to_string(),
            43092919u64,
            b"A".to_vec(),
            b"G".to_vec(),
            ClinVarAnnotation {
                variation_id: Some("12345".into()),
                significance: "Pathogenic".into(),
                review_status: "reviewed by expert panel".into(),
                conditions: vec!["Breast cancer".into()],
                submitter_count: 5,
                star_rating: 3,
            },
        )];
        let result = match_clinvar(&v, &db);
        assert!(result.is_some());
        assert_eq!(result.unwrap().significance, "Pathogenic");
    }

    #[test]
    fn test_clinvar_no_match() {
        let v = make_variant("chr1", 100, b"A", b"G");
        let db: Vec<(String, u64, Vec<u8>, Vec<u8>, ClinVarAnnotation)> = vec![];
        assert!(match_clinvar(&v, &db).is_none());
    }

    #[test]
    fn test_parse_clinvar_tsv() {
        let tsv = "chrom\tpos\tref\talt\tsignificance\treview_status\tconditions\tsubmitter_count\tstar_rating\n\
                   chr17\t43092919\tA\tG\tPathogenic\texpert panel\tBreast cancer\t5\t3\n";
        let db = parse_clinvar_tsv(tsv).unwrap();
        assert_eq!(db.len(), 1);
        assert_eq!(db[0].4.significance, "Pathogenic");
    }

    #[test]
    fn test_acmg_criterion_counts() {
        let ev = AcmgEvidence::new()
            .pvs1("Test")
            .strong_pathogenic("PS1", "Test")
            .moderate_pathogenic("PM2", "Test")
            .supporting_pathogenic("PP3", "Test");
        let result = ev.classify();
        assert_eq!(result.pvs_count, 1);
        assert_eq!(result.ps_count, 1);
        assert_eq!(result.pm_count, 1);
        assert_eq!(result.pp_count, 1);
    }
}
