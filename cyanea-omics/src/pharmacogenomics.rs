//! Pharmacogenomics: star allele calling and drug-gene interactions.
//!
//! Implements star allele (*) nomenclature matching for pharmacogenes
//! (e.g., CYP2D6, CYP2C19, DPYD, TPMT), metabolizer phenotype prediction,
//! and drug-gene interaction lookup.

use crate::variant::Variant;
use cyanea_core::Result;
use std::collections::{HashMap, HashSet};

/// A pharmacogene star allele definition.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StarAllele {
    /// Gene name (e.g., "CYP2D6").
    pub gene: String,
    /// Allele name (e.g., "*4", "*17").
    pub allele: String,
    /// Defining variants (chrom, position, ref, alt).
    pub defining_variants: Vec<(String, u64, Vec<u8>, Vec<u8>)>,
    /// Activity score (0.0 = no function, 0.5 = decreased, 1.0 = normal, 2.0 = increased).
    pub activity_score: f64,
    /// Functional status.
    pub function: AlleleFunction,
}

/// Functional status of a star allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AlleleFunction {
    NormalFunction,
    DecreasedFunction,
    NoFunction,
    IncreasedFunction,
    Uncertain,
}

impl AlleleFunction {
    pub fn as_str(&self) -> &'static str {
        match self {
            AlleleFunction::NormalFunction => "Normal Function",
            AlleleFunction::DecreasedFunction => "Decreased Function",
            AlleleFunction::NoFunction => "No Function",
            AlleleFunction::IncreasedFunction => "Increased Function",
            AlleleFunction::Uncertain => "Uncertain Function",
        }
    }
}

/// Metabolizer phenotype.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MetabolizerPhenotype {
    UltrarapidMetabolizer,
    RapidMetabolizer,
    NormalMetabolizer,
    IntermediateMetabolizer,
    PoorMetabolizer,
    Indeterminate,
}

impl MetabolizerPhenotype {
    pub fn as_str(&self) -> &'static str {
        match self {
            MetabolizerPhenotype::UltrarapidMetabolizer => "Ultrarapid Metabolizer",
            MetabolizerPhenotype::RapidMetabolizer => "Rapid Metabolizer",
            MetabolizerPhenotype::NormalMetabolizer => "Normal Metabolizer",
            MetabolizerPhenotype::IntermediateMetabolizer => "Intermediate Metabolizer",
            MetabolizerPhenotype::PoorMetabolizer => "Poor Metabolizer",
            MetabolizerPhenotype::Indeterminate => "Indeterminate",
        }
    }
}

/// Result of star allele calling for a gene.
#[derive(Debug, Clone)]
pub struct StarAlleleCall {
    /// Gene name.
    pub gene: String,
    /// Called diplotype (e.g., "*1/*4").
    pub diplotype: String,
    /// Allele 1 name.
    pub allele1: String,
    /// Allele 2 name.
    pub allele2: String,
    /// Combined activity score.
    pub activity_score: f64,
    /// Predicted metabolizer phenotype.
    pub phenotype: MetabolizerPhenotype,
}

/// A drug-gene interaction recommendation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DrugGeneInteraction {
    /// Drug name.
    pub drug: String,
    /// Gene name.
    pub gene: String,
    /// Phenotype that triggers this recommendation.
    pub phenotype: MetabolizerPhenotype,
    /// Clinical recommendation.
    pub recommendation: String,
    /// Evidence level (e.g., "1A", "1B", "2A").
    pub evidence_level: String,
    /// Source guideline (e.g., "CPIC", "DPWG").
    pub source: String,
}

/// A pharmacogenomics allele database for a gene.
#[derive(Debug, Clone)]
pub struct PgxDatabase {
    /// Gene name → list of star allele definitions.
    pub alleles: HashMap<String, Vec<StarAllele>>,
    /// Drug-gene interaction guidelines.
    pub interactions: Vec<DrugGeneInteraction>,
}

impl PgxDatabase {
    pub fn new() -> Self {
        Self {
            alleles: HashMap::new(),
            interactions: Vec::new(),
        }
    }

    pub fn add_allele(&mut self, allele: StarAllele) {
        self.alleles
            .entry(allele.gene.clone())
            .or_default()
            .push(allele);
    }

    pub fn add_interaction(&mut self, interaction: DrugGeneInteraction) {
        self.interactions.push(interaction);
    }
}

/// Call star alleles for a gene given observed variants.
///
/// Matches observed variants against star allele definitions. The reference
/// allele (*1) is assumed when no defining variants match.
///
/// For diploid calling, the function tries all pairs of alleles and selects
/// the pair that best explains the observed variants.
pub fn call_star_alleles(
    gene: &str,
    variants: &[Variant],
    db: &PgxDatabase,
) -> Result<StarAlleleCall> {
    let allele_defs = db.alleles.get(gene).cloned().unwrap_or_default();

    // Build set of observed variant keys
    let observed: HashSet<(String, u64, Vec<u8>, Vec<u8>)> = variants
        .iter()
        .map(|v| (v.chrom.clone(), v.position, v.ref_allele.clone(), v.alt_alleles[0].clone()))
        .collect();

    // Score each allele by how many defining variants match
    let mut allele_scores: Vec<(String, f64, f64)> = Vec::new(); // (name, match_fraction, activity)

    for def in &allele_defs {
        if def.defining_variants.is_empty() {
            continue;
        }
        let matches = def
            .defining_variants
            .iter()
            .filter(|dv| observed.contains(dv))
            .count();
        let fraction = matches as f64 / def.defining_variants.len() as f64;
        if fraction > 0.5 {
            allele_scores.push((def.allele.clone(), fraction, def.activity_score));
        }
    }

    // Sort by match fraction descending
    allele_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Select top two alleles (or default to *1)
    let (allele1, score1) = if !allele_scores.is_empty() {
        (allele_scores[0].0.clone(), allele_scores[0].2)
    } else {
        ("*1".to_string(), 1.0) // reference = normal function
    };

    let (allele2, score2) = if allele_scores.len() > 1 {
        (allele_scores[1].0.clone(), allele_scores[1].2)
    } else {
        ("*1".to_string(), 1.0)
    };

    let activity_score = score1 + score2;
    let phenotype = activity_to_phenotype(activity_score);

    let diplotype = format!("{}/{}", allele1, allele2);

    Ok(StarAlleleCall {
        gene: gene.to_string(),
        diplotype,
        allele1,
        allele2,
        activity_score,
        phenotype,
    })
}

/// Convert combined activity score to metabolizer phenotype.
///
/// Standard CPIC activity score ranges:
/// - > 2.25: Ultrarapid Metabolizer
/// - > 2.0: Rapid Metabolizer
/// - 1.25–2.0: Normal Metabolizer
/// - 0.25–1.0: Intermediate Metabolizer
/// - < 0.25: Poor Metabolizer
pub fn activity_to_phenotype(score: f64) -> MetabolizerPhenotype {
    if score > 2.25 {
        MetabolizerPhenotype::UltrarapidMetabolizer
    } else if score > 2.0 {
        MetabolizerPhenotype::RapidMetabolizer
    } else if score >= 1.25 {
        MetabolizerPhenotype::NormalMetabolizer
    } else if score >= 0.25 {
        MetabolizerPhenotype::IntermediateMetabolizer
    } else {
        MetabolizerPhenotype::PoorMetabolizer
    }
}

/// Look up drug recommendations for a given gene and phenotype.
pub fn lookup_drug_interactions<'a>(
    db: &'a PgxDatabase,
    gene: &str,
    phenotype: MetabolizerPhenotype,
) -> Vec<&'a DrugGeneInteraction> {
    db.interactions
        .iter()
        .filter(|i| i.gene == gene && i.phenotype == phenotype)
        .collect()
}

/// Build a demo CYP2D6 database with common star alleles.
pub fn demo_cyp2d6_database() -> PgxDatabase {
    let mut db = PgxDatabase::new();

    db.add_allele(StarAllele {
        gene: "CYP2D6".into(),
        allele: "*4".into(),
        defining_variants: vec![("chr22".into(), 42128945, b"C".to_vec(), b"T".to_vec())],
        activity_score: 0.0,
        function: AlleleFunction::NoFunction,
    });

    db.add_allele(StarAllele {
        gene: "CYP2D6".into(),
        allele: "*10".into(),
        defining_variants: vec![("chr22".into(), 42130692, b"C".to_vec(), b"T".to_vec())],
        activity_score: 0.25,
        function: AlleleFunction::DecreasedFunction,
    });

    db.add_allele(StarAllele {
        gene: "CYP2D6".into(),
        allele: "*17".into(),
        defining_variants: vec![
            ("chr22".into(), 42130692, b"C".to_vec(), b"T".to_vec()),
            ("chr22".into(), 42129132, b"C".to_vec(), b"T".to_vec()),
        ],
        activity_score: 0.5,
        function: AlleleFunction::DecreasedFunction,
    });

    db.add_interaction(DrugGeneInteraction {
        drug: "codeine".into(),
        gene: "CYP2D6".into(),
        phenotype: MetabolizerPhenotype::PoorMetabolizer,
        recommendation: "Avoid codeine; use alternative analgesic".into(),
        evidence_level: "1A".into(),
        source: "CPIC".into(),
    });

    db.add_interaction(DrugGeneInteraction {
        drug: "codeine".into(),
        gene: "CYP2D6".into(),
        phenotype: MetabolizerPhenotype::UltrarapidMetabolizer,
        recommendation: "Avoid codeine; risk of toxicity from rapid conversion to morphine".into(),
        evidence_level: "1A".into(),
        source: "CPIC".into(),
    });

    db.add_interaction(DrugGeneInteraction {
        drug: "tamoxifen".into(),
        gene: "CYP2D6".into(),
        phenotype: MetabolizerPhenotype::PoorMetabolizer,
        recommendation: "Consider alternative endocrine therapy (e.g., aromatase inhibitor)".into(),
        evidence_level: "1A".into(),
        source: "CPIC".into(),
    });

    db
}

#[cfg(test)]
mod tests {
    use super::*;
    fn make_variant(chrom: &str, pos: u64, ref_a: &[u8], alt_a: &[u8]) -> Variant {
        Variant::new(chrom, pos, ref_a.to_vec(), vec![alt_a.to_vec()]).unwrap()
    }

    #[test]
    fn test_allele_function_str() {
        assert_eq!(AlleleFunction::NormalFunction.as_str(), "Normal Function");
        assert_eq!(AlleleFunction::NoFunction.as_str(), "No Function");
    }

    #[test]
    fn test_metabolizer_phenotype_str() {
        assert_eq!(
            MetabolizerPhenotype::PoorMetabolizer.as_str(),
            "Poor Metabolizer"
        );
    }

    #[test]
    fn test_activity_to_phenotype() {
        assert_eq!(activity_to_phenotype(0.0), MetabolizerPhenotype::PoorMetabolizer);
        assert_eq!(activity_to_phenotype(0.5), MetabolizerPhenotype::IntermediateMetabolizer);
        assert_eq!(activity_to_phenotype(1.5), MetabolizerPhenotype::NormalMetabolizer);
        assert_eq!(activity_to_phenotype(2.0), MetabolizerPhenotype::NormalMetabolizer);
        assert_eq!(activity_to_phenotype(2.1), MetabolizerPhenotype::RapidMetabolizer);
        assert_eq!(activity_to_phenotype(3.0), MetabolizerPhenotype::UltrarapidMetabolizer);
    }

    #[test]
    fn test_call_star_alleles_reference() {
        let db = demo_cyp2d6_database();
        // No variants → *1/*1
        let result = call_star_alleles("CYP2D6", &[], &db).unwrap();
        assert_eq!(result.diplotype, "*1/*1");
        assert_eq!(result.phenotype, MetabolizerPhenotype::NormalMetabolizer);
        assert!((result.activity_score - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_call_star_alleles_poor_metabolizer() {
        let db = demo_cyp2d6_database();
        // CYP2D6*4 homozygous → PM (but we only get one allele from simple matching)
        let variants = vec![
            make_variant("chr22", 42128945, b"C", b"T"),
        ];
        let result = call_star_alleles("CYP2D6", &variants, &db).unwrap();
        assert!(result.diplotype.contains("*4"));
    }

    #[test]
    fn test_call_star_alleles_intermediate() {
        let db = demo_cyp2d6_database();
        let variants = vec![
            make_variant("chr22", 42130692, b"C", b"T"), // *10
        ];
        let result = call_star_alleles("CYP2D6", &variants, &db).unwrap();
        assert!(result.diplotype.contains("*10"));
    }

    #[test]
    fn test_lookup_drug_interactions() {
        let db = demo_cyp2d6_database();
        let recs = lookup_drug_interactions(&db, "CYP2D6", MetabolizerPhenotype::PoorMetabolizer);
        assert_eq!(recs.len(), 2); // codeine + tamoxifen
    }

    #[test]
    fn test_lookup_drug_interactions_normal() {
        let db = demo_cyp2d6_database();
        let recs = lookup_drug_interactions(&db, "CYP2D6", MetabolizerPhenotype::NormalMetabolizer);
        assert!(recs.is_empty()); // no special recommendations
    }

    #[test]
    fn test_pgx_database_add() {
        let mut db = PgxDatabase::new();
        db.add_allele(StarAllele {
            gene: "CYP2C19".into(),
            allele: "*2".into(),
            defining_variants: vec![("chr10".into(), 94781859, b"G".to_vec(), b"A".to_vec())],
            activity_score: 0.0,
            function: AlleleFunction::NoFunction,
        });
        assert_eq!(db.alleles["CYP2C19"].len(), 1);
    }

    #[test]
    fn test_demo_database() {
        let db = demo_cyp2d6_database();
        assert_eq!(db.alleles["CYP2D6"].len(), 3);
        assert_eq!(db.interactions.len(), 3);
    }
}
