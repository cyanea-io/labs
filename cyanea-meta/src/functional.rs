//! Functional annotation and metabolic analysis.
//!
//! Maps taxa to functional categories (KEGG/COG), aggregates to pathways,
//! and computes functional diversity.

use std::collections::HashMap;
use crate::error::{MetaError, Result};

/// Functional category ID and its abundance.
#[derive(Debug, Clone)]
pub struct FunctionalProfile {
    /// Function ID (KEGG K number, COG ID, etc.) → abundance.
    pub abundances: HashMap<String, f64>,
}

impl FunctionalProfile {
    /// Create a new empty functional profile.
    pub fn new() -> Self {
        Self {
            abundances: HashMap::new(),
        }
    }

    /// Insert or update a function's abundance.
    pub fn insert(&mut self, function_id: &str, abundance: f64) {
        self.abundances.insert(function_id.to_string(), abundance);
    }

    /// Get abundance for a function.
    pub fn get(&self, function_id: &str) -> Option<f64> {
        self.abundances.get(function_id).copied()
    }

    /// Total abundance across all functions.
    pub fn total_abundance(&self) -> f64 {
        self.abundances.values().sum()
    }

    /// Number of functions in the profile.
    pub fn len(&self) -> usize {
        self.abundances.len()
    }

    /// Check if profile is empty.
    pub fn is_empty(&self) -> bool {
        self.abundances.is_empty()
    }

    /// Get functions sorted by abundance (descending).
    pub fn top_functions(&self, n: usize) -> Vec<(String, f64)> {
        let mut vec: Vec<_> = self
            .abundances
            .iter()
            .map(|(id, &abund)| (id.clone(), abund))
            .collect();
        vec.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        vec.into_iter().take(n).collect()
    }
}

impl Default for FunctionalProfile {
    fn default() -> Self {
        Self::new()
    }
}

/// Map taxa to functional categories.
///
/// Given a taxonomic profile and a taxon → functions mapping, computes
/// the functional profile by aggregating abundances.
///
/// # Arguments
///
/// * `taxon_abundances` — taxon ID → abundance
/// * `taxon_functions` — taxon ID → list of function IDs
///
/// # Errors
///
/// Returns an error if abundances or mapping is empty.
pub fn map_to_functions(
    taxon_abundances: &HashMap<u32, f64>,
    taxon_functions: &HashMap<u32, Vec<String>>,
) -> Result<FunctionalProfile> {
    if taxon_abundances.is_empty() {
        return Err(MetaError::Functional(
            "taxon abundances map is empty".into(),
        ));
    }

    let mut profile = FunctionalProfile::new();

    for (taxid, abundance) in taxon_abundances {
        if let Some(functions) = taxon_functions.get(taxid) {
            // Distribute abundance equally among functions for this taxon
            let per_function = abundance / functions.len() as f64;
            for func_id in functions {
                let current = profile.get(func_id).unwrap_or(0.0);
                profile.insert(func_id, current + per_function);
            }
        }
    }

    if profile.is_empty() {
        return Err(MetaError::Functional(
            "no functions mapped for any taxon".into(),
        ));
    }

    Ok(profile)
}

/// Pathway abundance: aggregate gene families to pathways.
///
/// Given a functional profile and a function → pathways mapping,
/// computes pathway-level abundance by summing constituent functions.
///
/// # Arguments
///
/// * `functional_profile` — function ID → abundance
/// * `function_pathways` — function ID → list of pathway IDs
///
/// # Errors
///
/// Returns an error if profile is empty.
pub fn pathway_abundance(
    functional_profile: &FunctionalProfile,
    function_pathways: &HashMap<String, Vec<String>>,
) -> Result<HashMap<String, f64>> {
    if functional_profile.is_empty() {
        return Err(MetaError::Functional(
            "functional profile is empty".into(),
        ));
    }

    let mut pathways: HashMap<String, f64> = HashMap::new();

    for (func_id, abundance) in &functional_profile.abundances {
        if let Some(pathway_ids) = function_pathways.get(func_id) {
            for pathway_id in pathway_ids {
                *pathways.entry(pathway_id.clone()).or_insert(0.0) += abundance;
            }
        }
    }

    if pathways.is_empty() {
        return Err(MetaError::Functional(
            "no pathways mapped for any function".into(),
        ));
    }

    Ok(pathways)
}

/// Functional diversity: Shannon entropy of functional abundance distribution.
///
/// Measures how evenly abundant different functions are.
///
/// # Errors
///
/// Returns an error if profile is empty.
pub fn functional_diversity(profile: &FunctionalProfile) -> Result<f64> {
    if profile.is_empty() {
        return Err(MetaError::Functional(
            "profile is empty".into(),
        ));
    }

    let total = profile.total_abundance();
    if total <= 0.0 {
        return Err(MetaError::Functional(
            "total abundance must be positive".into(),
        ));
    }

    let mut diversity = 0.0;
    for &abundance in profile.abundances.values() {
        if abundance > 0.0 {
            let p = abundance / total;
            diversity -= p * p.ln();
        }
    }

    Ok(diversity)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn functional_profile_basic() {
        let mut profile = FunctionalProfile::new();
        profile.insert("K001", 100.0);
        profile.insert("K002", 50.0);
        assert_eq!(profile.len(), 2);
        assert!((profile.total_abundance() - 150.0).abs() < 1e-10);
    }

    #[test]
    fn map_to_functions_basic() {
        let mut abundances = HashMap::new();
        abundances.insert(1u32, 100.0);
        abundances.insert(2u32, 50.0);

        let mut mapping = HashMap::new();
        mapping.insert(1u32, vec!["K001".to_string(), "K002".to_string()]);
        mapping.insert(2u32, vec!["K002".to_string(), "K003".to_string()]);

        let profile = map_to_functions(&abundances, &mapping).unwrap();
        assert_eq!(profile.len(), 3);
        // K001: 100 / 2 = 50
        // K002: 50 + 25 = 75
        // K003: 25
        assert!((profile.get("K001").unwrap() - 50.0).abs() < 1e-10);
        assert!((profile.get("K002").unwrap() - 75.0).abs() < 1e-10);
    }

    #[test]
    fn pathway_abundance_basic() {
        let mut profile = FunctionalProfile::new();
        profile.insert("K001", 50.0);
        profile.insert("K002", 75.0);

        let mut function_pathways = HashMap::new();
        function_pathways.insert("K001".to_string(), vec!["path1".to_string()]);
        function_pathways.insert("K002".to_string(), vec!["path1".to_string(), "path2".to_string()]);

        let pathways = pathway_abundance(&profile, &function_pathways).unwrap();
        // path1: 50 + 75 = 125
        // path2: 75
        assert!((pathways.get("path1").unwrap() - 125.0).abs() < 1e-10);
        assert!((pathways.get("path2").unwrap() - 75.0).abs() < 1e-10);
    }

    #[test]
    fn functional_diversity_uniform() {
        let mut profile = FunctionalProfile::new();
        profile.insert("K001", 25.0);
        profile.insert("K002", 25.0);
        profile.insert("K003", 25.0);
        profile.insert("K004", 25.0);

        let diversity = functional_diversity(&profile).unwrap();
        // Uniform: H = ln(4)
        let expected = (4.0f64).ln();
        assert!((diversity - expected).abs() < 1e-10);
    }

    #[test]
    fn empty_profile_error() {
        let profile = FunctionalProfile::new();
        assert!(functional_diversity(&profile).is_err());
    }
}
