//! Effect size measures for quantifying the magnitude of statistical effects.
//!
//! Provides common effect sizes used in biostatistics:
//! - [`cohens_d`] — standardized mean difference between two groups
//! - [`eta_squared`] — proportion of variance explained (ANOVA)
//! - [`odds_ratio`] — association strength for 2×2 tables
//! - [`relative_risk`] — risk ratio for 2×2 tables

use cyanea_core::{CyaneaError, Result};

use crate::descriptive;

/// Cohen's d: standardized difference between two group means.
///
/// Uses pooled standard deviation: `d = (mean1 - mean2) / s_pooled`.
///
/// Each group needs at least 2 observations.
pub fn cohens_d(group1: &[f64], group2: &[f64]) -> Result<f64> {
    if group1.len() < 2 || group2.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "cohens_d: each group needs at least 2 observations".into(),
        ));
    }

    let mean1 = descriptive::mean(group1)?;
    let mean2 = descriptive::mean(group2)?;
    let var1 = descriptive::variance(group1, 1)?;
    let var2 = descriptive::variance(group2, 1)?;
    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;

    let pooled_var = ((n1 - 1.0) * var1 + (n2 - 1.0) * var2) / (n1 + n2 - 2.0);
    let s_pooled = pooled_var.sqrt();

    if s_pooled == 0.0 {
        return Ok(0.0);
    }

    Ok((mean1 - mean2) / s_pooled)
}

/// Eta-squared (η²): proportion of total variance explained by group membership.
///
/// Computed from ANOVA sums of squares: `η² = SS_between / SS_total`.
///
/// Returns a value in [0, 1]. Values near 0 indicate no effect, near 1
/// indicate the grouping explains most of the variance.
pub fn eta_squared(groups: &[&[f64]]) -> Result<f64> {
    if groups.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "eta_squared: need at least 2 groups".into(),
        ));
    }
    for (i, g) in groups.iter().enumerate() {
        if g.is_empty() {
            return Err(CyaneaError::InvalidInput(
                format!("eta_squared: group {} is empty", i),
            ));
        }
    }

    let n_total: usize = groups.iter().map(|g| g.len()).sum();
    let grand_sum: f64 = groups.iter().flat_map(|g| g.iter()).sum();
    let grand_mean = grand_sum / n_total as f64;

    let ss_between: f64 = groups
        .iter()
        .map(|g| {
            let group_mean: f64 = g.iter().sum::<f64>() / g.len() as f64;
            g.len() as f64 * (group_mean - grand_mean).powi(2)
        })
        .sum();

    let ss_total: f64 = groups
        .iter()
        .flat_map(|g| g.iter())
        .map(|&x| (x - grand_mean).powi(2))
        .sum();

    if ss_total == 0.0 {
        return Ok(0.0);
    }

    Ok(ss_between / ss_total)
}

/// Odds ratio for a 2×2 contingency table.
///
/// The table is `[[a, b], [c, d]]`:
///
/// ```text
///           Exposed   Unexposed
/// Cases        a          b
/// Controls     c          d
/// ```
///
/// OR = (a × d) / (b × c). A value of 1.0 means no association.
pub fn odds_ratio(table: &[[usize; 2]; 2]) -> Result<f64> {
    let a = table[0][0] as f64;
    let b = table[0][1] as f64;
    let c = table[1][0] as f64;
    let d = table[1][1] as f64;

    let denom = b * c;
    if denom == 0.0 {
        return Err(CyaneaError::InvalidInput(
            "odds_ratio: zero in denominator cell (b×c = 0)".into(),
        ));
    }

    Ok((a * d) / denom)
}

/// Relative risk (risk ratio) for a 2×2 contingency table.
///
/// The table is `[[a, b], [c, d]]`:
///
/// ```text
///              Event   No Event
/// Exposed        a        b
/// Unexposed      c        d
/// ```
///
/// RR = [a/(a+b)] / [c/(c+d)]. A value of 1.0 means no difference in risk.
pub fn relative_risk(table: &[[usize; 2]; 2]) -> Result<f64> {
    let a = table[0][0] as f64;
    let b = table[0][1] as f64;
    let c = table[1][0] as f64;
    let d = table[1][1] as f64;

    let risk_exposed = a + b;
    let risk_unexposed = c + d;

    if risk_exposed == 0.0 || risk_unexposed == 0.0 {
        return Err(CyaneaError::InvalidInput(
            "relative_risk: row total is zero".into(),
        ));
    }

    let rr_exposed = a / risk_exposed;
    let rr_unexposed = c / risk_unexposed;

    if rr_unexposed == 0.0 {
        return Err(CyaneaError::InvalidInput(
            "relative_risk: unexposed risk is zero".into(),
        ));
    }

    Ok(rr_exposed / rr_unexposed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cohens_d_zero_difference() {
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let d = cohens_d(&g1, &g2).unwrap();
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn cohens_d_large_effect() {
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [10.0, 11.0, 12.0, 13.0, 14.0];
        let d = cohens_d(&g1, &g2).unwrap();
        // Large negative effect (g1 < g2)
        assert!(d < -3.0, "d={}", d);
    }

    #[test]
    fn cohens_d_small_effect() {
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [1.5, 2.5, 3.5, 4.5, 5.5];
        let d = cohens_d(&g1, &g2).unwrap();
        // d ≈ -0.316 (small effect)
        assert!(d.abs() < 0.5, "d={}", d);
    }

    #[test]
    fn cohens_d_too_few() {
        assert!(cohens_d(&[1.0], &[2.0, 3.0]).is_err());
    }

    #[test]
    fn eta_squared_no_effect() {
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let eta = eta_squared(&[&g1, &g2]).unwrap();
        assert!(eta < 0.01, "eta={}", eta);
    }

    #[test]
    fn eta_squared_strong_effect() {
        let g1 = [1.0, 2.0, 3.0];
        let g2 = [100.0, 101.0, 102.0];
        let g3 = [200.0, 201.0, 202.0];
        let eta = eta_squared(&[&g1, &g2, &g3]).unwrap();
        assert!(eta > 0.99, "eta={}", eta);
    }

    #[test]
    fn eta_squared_range() {
        let g1 = [1.0, 2.0, 3.0, 4.0];
        let g2 = [3.0, 4.0, 5.0, 6.0];
        let eta = eta_squared(&[&g1, &g2]).unwrap();
        assert!(eta >= 0.0 && eta <= 1.0, "eta={}", eta);
    }

    #[test]
    fn odds_ratio_no_association() {
        let table = [[10, 10], [10, 10]];
        let or = odds_ratio(&table).unwrap();
        assert!((or - 1.0).abs() < 1e-10);
    }

    #[test]
    fn odds_ratio_strong_positive() {
        let table = [[9, 1], [1, 9]];
        let or = odds_ratio(&table).unwrap();
        assert!((or - 81.0).abs() < 1e-10);
    }

    #[test]
    fn odds_ratio_zero_denominator() {
        let table = [[5, 0], [3, 2]];
        assert!(odds_ratio(&table).is_err());
    }

    #[test]
    fn relative_risk_no_difference() {
        let table = [[10, 90], [10, 90]];
        let rr = relative_risk(&table).unwrap();
        assert!((rr - 1.0).abs() < 1e-10);
    }

    #[test]
    fn relative_risk_doubled() {
        // Exposed: 20/100, Unexposed: 10/100
        let table = [[20, 80], [10, 90]];
        let rr = relative_risk(&table).unwrap();
        assert!((rr - 2.0).abs() < 1e-10);
    }

    #[test]
    fn relative_risk_zero_row() {
        let table = [[0, 0], [5, 5]];
        assert!(relative_risk(&table).is_err());
    }
}
