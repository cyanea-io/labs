//! VCF variant manipulation operations.
//!
//! Provides normalization, splitting/joining multi-allelic variants, merging,
//! filtering with expression parsing, set operations, and detailed statistics.

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::variant::{Variant, VariantFilter};

// ===========================================================================
// Normalization
// ===========================================================================

/// Split a multi-allelic variant into biallelic records.
///
/// Each resulting variant has exactly one alternate allele.
/// A variant with a single alternate allele is returned unchanged (in a Vec of one).
pub fn split_multiallelic(variant: &Variant) -> Vec<Variant> {
    if variant.alt_alleles.len() <= 1 {
        return vec![variant.clone()];
    }

    variant
        .alt_alleles
        .iter()
        .map(|alt| Variant {
            chrom: variant.chrom.clone(),
            position: variant.position,
            id: variant.id.clone(),
            ref_allele: variant.ref_allele.clone(),
            alt_alleles: vec![alt.clone()],
            quality: variant.quality,
            filter: variant.filter.clone(),
        })
        .collect()
}

/// Join multiple biallelic variants at the same position into a single multi-allelic variant.
///
/// Returns `None` if variants have different chrom/position/ref_allele.
pub fn join_biallelic(variants: &[Variant]) -> Option<Variant> {
    if variants.is_empty() {
        return None;
    }

    let first = &variants[0];

    // Verify all have same chrom, position, ref_allele
    for v in &variants[1..] {
        if v.chrom != first.chrom || v.position != first.position || v.ref_allele != first.ref_allele
        {
            return None;
        }
    }

    let mut alt_alleles: Vec<Vec<u8>> = Vec::new();
    for v in variants {
        for alt in &v.alt_alleles {
            if !alt_alleles.contains(alt) {
                alt_alleles.push(alt.clone());
            }
        }
    }

    // Use quality from the first variant that has one
    let quality = variants.iter().find_map(|v| v.quality);

    Some(Variant {
        chrom: first.chrom.clone(),
        position: first.position,
        id: first.id.clone(),
        ref_allele: first.ref_allele.clone(),
        alt_alleles,
        quality,
        filter: first.filter.clone(),
    })
}

/// Left-align and trim a variant's alleles.
///
/// Removes common suffix, then common prefix (preserving at least 1 base in each
/// allele), adjusting position accordingly.
pub fn normalize_variant(variant: &Variant) -> Variant {
    if variant.alt_alleles.is_empty() {
        return variant.clone();
    }

    let mut ref_allele = variant.ref_allele.clone();
    let mut alt_alleles = variant.alt_alleles.clone();
    let mut position = variant.position;

    // Trim common suffix
    loop {
        if ref_allele.len() <= 1 || alt_alleles.iter().any(|a| a.len() <= 1) {
            break;
        }
        let last_ref = *ref_allele.last().unwrap();
        if alt_alleles.iter().all(|a| *a.last().unwrap() == last_ref) {
            ref_allele.pop();
            for alt in &mut alt_alleles {
                alt.pop();
            }
        } else {
            break;
        }
    }

    // Trim common prefix
    loop {
        if ref_allele.len() <= 1 || alt_alleles.iter().any(|a| a.len() <= 1) {
            break;
        }
        let first_ref = ref_allele[0];
        if alt_alleles.iter().all(|a| a[0] == first_ref) {
            ref_allele.remove(0);
            for alt in &mut alt_alleles {
                alt.remove(0);
            }
            position += 1;
        } else {
            break;
        }
    }

    Variant {
        chrom: variant.chrom.clone(),
        position,
        id: variant.id.clone(),
        ref_allele,
        alt_alleles,
        quality: variant.quality,
        filter: variant.filter.clone(),
    }
}

// ===========================================================================
// Merging
// ===========================================================================

/// Merge multiple sorted variant sets into one sorted set.
///
/// Variants at the same position are combined by merging alternate alleles.
pub fn merge_variants(inputs: &[&[Variant]]) -> Vec<Variant> {
    use std::collections::BTreeMap;

    // Group by (chrom, position, ref_allele)
    let mut groups: BTreeMap<(String, u64, Vec<u8>), Vec<&Variant>> = BTreeMap::new();

    for &input in inputs {
        for v in input {
            let key = (v.chrom.clone(), v.position, v.ref_allele.clone());
            groups.entry(key).or_default().push(v);
        }
    }

    let mut result = Vec::new();
    for ((_chrom, _pos, _ref_allele), variants) in groups {
        let owned: Vec<Variant> = variants.into_iter().cloned().collect();
        if let Some(merged) = join_biallelic(&owned) {
            result.push(merged);
        } else {
            // If join fails (shouldn't, since they share key), add individually
            result.extend(owned);
        }
    }

    result
}

// ===========================================================================
// Filter expressions
// ===========================================================================

/// A parsed filter expression AST node.
#[derive(Debug, Clone)]
pub enum FilterExpr {
    /// Compare quality: QUAL > 30, QUAL >= 30, QUAL < 30, etc.
    QualCmp {
        op: CmpOp,
        value: f64,
    },
    /// Compare variant type: TYPE == SNV, TYPE == INDEL, etc.
    TypeEq(String),
    /// Compare chromosome: CHROM == chr1.
    ChromEq(String),
    /// FILTER == PASS.
    FilterPass,
    /// Logical AND.
    And(Box<FilterExpr>, Box<FilterExpr>),
    /// Logical OR.
    Or(Box<FilterExpr>, Box<FilterExpr>),
    /// Logical NOT.
    Not(Box<FilterExpr>),
}

/// Comparison operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmpOp {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
    Ne,
}

/// Parse a filter expression string.
///
/// Supports:
/// - `QUAL>30`, `QUAL>=30`, `QUAL<10`, `QUAL<=10`, `QUAL==30`, `QUAL!=0`
/// - `TYPE==SNV`, `TYPE==INDEL`
/// - `CHROM==chr1`
/// - `FILTER==PASS`
/// - `&&` (AND), `||` (OR), `!` (NOT)
/// - Parentheses for grouping
pub fn parse_filter(expr: &str) -> Result<FilterExpr> {
    let tokens = tokenize(expr)?;
    let (node, rest) = parse_or(&tokens)?;
    if !rest.is_empty() {
        return Err(CyaneaError::Parse(format!(
            "unexpected tokens in filter expression: {:?}",
            rest
        )));
    }
    Ok(node)
}

/// Apply a filter expression to a variant, returning true if it passes.
pub fn eval_filter(variant: &Variant, expr: &FilterExpr) -> bool {
    match expr {
        FilterExpr::QualCmp { op, value } => {
            let q = variant.quality.unwrap_or(0.0);
            match op {
                CmpOp::Gt => q > *value,
                CmpOp::Ge => q >= *value,
                CmpOp::Lt => q < *value,
                CmpOp::Le => q <= *value,
                CmpOp::Eq => (q - value).abs() < f64::EPSILON,
                CmpOp::Ne => (q - value).abs() >= f64::EPSILON,
            }
        }
        FilterExpr::TypeEq(t) => {
            let t_upper = t.to_uppercase();
            match t_upper.as_str() {
                "SNV" | "SNP" => variant.is_snv(),
                "INDEL" => variant.is_indel(),
                _ => false,
            }
        }
        FilterExpr::ChromEq(c) => variant.chrom == *c,
        FilterExpr::FilterPass => variant.filter == VariantFilter::Pass,
        FilterExpr::And(a, b) => eval_filter(variant, a) && eval_filter(variant, b),
        FilterExpr::Or(a, b) => eval_filter(variant, a) || eval_filter(variant, b),
        FilterExpr::Not(a) => !eval_filter(variant, a),
    }
}

/// Filter variants using a parsed filter expression.
pub fn filter_variants(variants: &[Variant], expr: &FilterExpr) -> Vec<Variant> {
    variants
        .iter()
        .filter(|v| eval_filter(v, expr))
        .cloned()
        .collect()
}

// --- Tokenizer and recursive-descent parser ---

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Ident(String),
    Number(f64),
    Op(String), // ==, !=, >, >=, <, <=
    And,        // &&
    Or,         // ||
    Not,        // !
    LParen,
    RParen,
}

fn tokenize(expr: &str) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let bytes = expr.as_bytes();
    let mut i = 0;

    while i < bytes.len() {
        match bytes[i] {
            b' ' | b'\t' => i += 1,
            b'(' => {
                tokens.push(Token::LParen);
                i += 1;
            }
            b')' => {
                tokens.push(Token::RParen);
                i += 1;
            }
            b'&' if i + 1 < bytes.len() && bytes[i + 1] == b'&' => {
                tokens.push(Token::And);
                i += 2;
            }
            b'|' if i + 1 < bytes.len() && bytes[i + 1] == b'|' => {
                tokens.push(Token::Or);
                i += 2;
            }
            b'!' if i + 1 < bytes.len() && bytes[i + 1] == b'=' => {
                tokens.push(Token::Op("!=".to_string()));
                i += 2;
            }
            b'!' => {
                tokens.push(Token::Not);
                i += 1;
            }
            b'=' if i + 1 < bytes.len() && bytes[i + 1] == b'=' => {
                tokens.push(Token::Op("==".to_string()));
                i += 2;
            }
            b'>' if i + 1 < bytes.len() && bytes[i + 1] == b'=' => {
                tokens.push(Token::Op(">=".to_string()));
                i += 2;
            }
            b'<' if i + 1 < bytes.len() && bytes[i + 1] == b'=' => {
                tokens.push(Token::Op("<=".to_string()));
                i += 2;
            }
            b'>' => {
                tokens.push(Token::Op(">".to_string()));
                i += 1;
            }
            b'<' => {
                tokens.push(Token::Op("<".to_string()));
                i += 1;
            }
            b'0'..=b'9' | b'.' if {
                // Check it's a number, not just a dot in an identifier
                bytes[i] != b'.' || (i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit())
            } =>
            {
                let start = i;
                while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                    i += 1;
                }
                let num_str = &expr[start..i];
                let num: f64 = num_str.parse().map_err(|_| {
                    CyaneaError::Parse(format!("invalid number in filter: '{}'", num_str))
                })?;
                tokens.push(Token::Number(num));
            }
            _ if bytes[i].is_ascii_alphabetic() || bytes[i] == b'_' => {
                let start = i;
                while i < bytes.len()
                    && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_' || bytes[i] == b'.')
                {
                    i += 1;
                }
                tokens.push(Token::Ident(expr[start..i].to_string()));
            }
            _ => {
                return Err(CyaneaError::Parse(format!(
                    "unexpected character '{}' in filter expression",
                    bytes[i] as char
                )));
            }
        }
    }

    Ok(tokens)
}

// Recursive-descent: or -> and -> not -> primary
fn parse_or<'a>(tokens: &'a [Token]) -> Result<(FilterExpr, &'a [Token])> {
    let (mut left, mut rest) = parse_and(tokens)?;
    while !rest.is_empty() && rest[0] == Token::Or {
        let (right, r) = parse_and(&rest[1..])?;
        left = FilterExpr::Or(Box::new(left), Box::new(right));
        rest = r;
    }
    Ok((left, rest))
}

fn parse_and<'a>(tokens: &'a [Token]) -> Result<(FilterExpr, &'a [Token])> {
    let (mut left, mut rest) = parse_not(tokens)?;
    while !rest.is_empty() && rest[0] == Token::And {
        let (right, r) = parse_not(&rest[1..])?;
        left = FilterExpr::And(Box::new(left), Box::new(right));
        rest = r;
    }
    Ok((left, rest))
}

fn parse_not<'a>(tokens: &'a [Token]) -> Result<(FilterExpr, &'a [Token])> {
    if !tokens.is_empty() && tokens[0] == Token::Not {
        let (inner, rest) = parse_not(&tokens[1..])?;
        Ok((FilterExpr::Not(Box::new(inner)), rest))
    } else {
        parse_primary(tokens)
    }
}

fn parse_primary<'a>(tokens: &'a [Token]) -> Result<(FilterExpr, &'a [Token])> {
    if tokens.is_empty() {
        return Err(CyaneaError::Parse(
            "unexpected end of filter expression".into(),
        ));
    }

    // Parenthesized expression
    if tokens[0] == Token::LParen {
        let (inner, rest) = parse_or(&tokens[1..])?;
        if rest.is_empty() || rest[0] != Token::RParen {
            return Err(CyaneaError::Parse(
                "missing closing parenthesis in filter".into(),
            ));
        }
        return Ok((inner, &rest[1..]));
    }

    // IDENT op VALUE
    if let Token::Ident(name) = &tokens[0] {
        if tokens.len() < 3 {
            return Err(CyaneaError::Parse(format!(
                "incomplete filter expression after '{}'",
                name
            )));
        }

        if let Token::Op(op) = &tokens[1] {
            match name.as_str() {
                "QUAL" => {
                    if let Token::Number(val) = &tokens[2] {
                        let cmp = match op.as_str() {
                            ">" => CmpOp::Gt,
                            ">=" => CmpOp::Ge,
                            "<" => CmpOp::Lt,
                            "<=" => CmpOp::Le,
                            "==" => CmpOp::Eq,
                            "!=" => CmpOp::Ne,
                            _ => {
                                return Err(CyaneaError::Parse(format!(
                                    "unsupported operator '{}'",
                                    op
                                )))
                            }
                        };
                        return Ok((FilterExpr::QualCmp { op: cmp, value: *val }, &tokens[3..]));
                    }
                }
                "TYPE" => {
                    if op == "==" {
                        if let Token::Ident(val) = &tokens[2] {
                            return Ok((FilterExpr::TypeEq(val.clone()), &tokens[3..]));
                        }
                    }
                }
                "CHROM" => {
                    if op == "==" {
                        if let Token::Ident(val) = &tokens[2] {
                            return Ok((FilterExpr::ChromEq(val.clone()), &tokens[3..]));
                        }
                    }
                }
                "FILTER" => {
                    if op == "==" {
                        if let Token::Ident(val) = &tokens[2] {
                            if val == "PASS" {
                                return Ok((FilterExpr::FilterPass, &tokens[3..]));
                            }
                        }
                    }
                }
                _ => {}
            }
        }

        return Err(CyaneaError::Parse(format!(
            "unsupported filter field or expression: '{}'",
            name
        )));
    }

    Err(CyaneaError::Parse(format!(
        "unexpected token in filter: {:?}",
        tokens[0]
    )))
}

// ===========================================================================
// Set operations
// ===========================================================================

/// Key for matching variants: (chrom, position, ref_allele, sorted alt_alleles).
fn variant_key(v: &Variant) -> (String, u64, Vec<u8>, Vec<Vec<u8>>) {
    let mut alts = v.alt_alleles.clone();
    alts.sort();
    (v.chrom.clone(), v.position, v.ref_allele.clone(), alts)
}

/// Variants present in both `a` and `b` (by position + alleles).
pub fn variant_intersection(a: &[Variant], b: &[Variant]) -> Vec<Variant> {
    let b_keys: std::collections::HashSet<_> = b.iter().map(|v| variant_key(v)).collect();
    a.iter()
        .filter(|v| b_keys.contains(&variant_key(v)))
        .cloned()
        .collect()
}

/// Variants in `a` but not in `b` (by position + alleles).
pub fn variant_complement(a: &[Variant], b: &[Variant]) -> Vec<Variant> {
    let b_keys: std::collections::HashSet<_> = b.iter().map(|v| variant_key(v)).collect();
    a.iter()
        .filter(|v| !b_keys.contains(&variant_key(v)))
        .cloned()
        .collect()
}

/// Concordance statistics between two variant sets.
#[derive(Debug, Clone, Default)]
pub struct ConcordanceStats {
    /// Variants shared between both sets.
    pub shared: usize,
    /// Variants only in set A.
    pub a_only: usize,
    /// Variants only in set B.
    pub b_only: usize,
}

/// Compute concordance between two variant sets.
pub fn variant_concordance(a: &[Variant], b: &[Variant]) -> ConcordanceStats {
    let a_keys: std::collections::HashSet<_> = a.iter().map(|v| variant_key(v)).collect();
    let b_keys: std::collections::HashSet<_> = b.iter().map(|v| variant_key(v)).collect();

    let shared = a_keys.intersection(&b_keys).count();
    let a_only = a_keys.difference(&b_keys).count();
    let b_only = b_keys.difference(&a_keys).count();

    ConcordanceStats {
        shared,
        a_only,
        b_only,
    }
}

// ===========================================================================
// Detailed statistics
// ===========================================================================

/// Detailed variant-level statistics.
#[derive(Debug, Clone)]
pub struct DetailedVcfStats {
    /// Total number of variants.
    pub total: u64,
    /// Number of SNVs.
    pub snv_count: u64,
    /// Number of insertions.
    pub insertion_count: u64,
    /// Number of deletions.
    pub deletion_count: u64,
    /// Transition/transversion ratio.
    pub ti_tv_ratio: f64,
    /// Variants that pass all filters.
    pub pass_count: u64,
    /// Variants with more than one alternate allele.
    pub multiallelic_count: u64,
    /// Quality score histogram (bin_start, count).
    pub quality_histogram: Vec<(f64, u64)>,
}

/// Compute detailed statistics over a variant set.
pub fn detailed_vcf_stats(variants: &[Variant]) -> DetailedVcfStats {
    let mut total: u64 = 0;
    let mut snv_count: u64 = 0;
    let mut insertion_count: u64 = 0;
    let mut deletion_count: u64 = 0;
    let mut transitions: u64 = 0;
    let mut transversions: u64 = 0;
    let mut pass_count: u64 = 0;
    let mut multiallelic_count: u64 = 0;

    // Quality histogram bins: 0-10, 10-20, ..., 90-100, 100+
    let bin_count = 11;
    let mut quality_bins = vec![0u64; bin_count];

    for v in variants {
        total += 1;

        if v.is_snv() {
            snv_count += 1;
            if v.is_transition() {
                transitions += 1;
            } else if v.is_transversion() {
                transversions += 1;
            }
        }

        // Check each alt allele for insertion/deletion
        for alt in &v.alt_alleles {
            if alt.len() > v.ref_allele.len() {
                insertion_count += 1;
            } else if alt.len() < v.ref_allele.len() {
                deletion_count += 1;
            }
        }

        if v.filter == VariantFilter::Pass {
            pass_count += 1;
        }

        if v.alt_alleles.len() > 1 {
            multiallelic_count += 1;
        }

        if let Some(q) = v.quality {
            let bin = (q / 10.0).floor() as usize;
            let bin = bin.min(bin_count - 1);
            quality_bins[bin] += 1;
        }
    }

    let ti_tv_ratio = if transversions > 0 {
        transitions as f64 / transversions as f64
    } else if transitions > 0 {
        f64::INFINITY
    } else {
        0.0
    };

    let quality_histogram: Vec<(f64, u64)> = quality_bins
        .into_iter()
        .enumerate()
        .filter(|(_, count)| *count > 0)
        .map(|(i, count)| (i as f64 * 10.0, count))
        .collect();

    DetailedVcfStats {
        total,
        snv_count,
        insertion_count,
        deletion_count,
        ti_tv_ratio,
        pass_count,
        multiallelic_count,
        quality_histogram,
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn snv(chrom: &str, pos: u64, ref_a: &str, alt: &str, qual: Option<f64>) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            position: pos,
            id: None,
            ref_allele: ref_a.as_bytes().to_vec(),
            alt_alleles: vec![alt.as_bytes().to_vec()],
            quality: qual,
            filter: VariantFilter::Pass,
        }
    }

    fn variant(
        chrom: &str,
        pos: u64,
        ref_a: &str,
        alts: &[&str],
        qual: Option<f64>,
        filter: VariantFilter,
    ) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            position: pos,
            id: None,
            ref_allele: ref_a.as_bytes().to_vec(),
            alt_alleles: alts.iter().map(|a| a.as_bytes().to_vec()).collect(),
            quality: qual,
            filter,
        }
    }

    // -----------------------------------------------------------------------
    // Split / Join multi-allelic
    // -----------------------------------------------------------------------

    #[test]
    fn split_multiallelic_single() {
        let v = snv("chr1", 100, "A", "G", Some(30.0));
        let result = split_multiallelic(&v);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].alt_alleles.len(), 1);
    }

    #[test]
    fn split_multiallelic_multi() {
        let v = variant("chr1", 100, "A", &["G", "T", "C"], Some(30.0), VariantFilter::Pass);
        let result = split_multiallelic(&v);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(result[1].alt_alleles, vec![b"T".to_vec()]);
        assert_eq!(result[2].alt_alleles, vec![b"C".to_vec()]);
        // All share same position and ref
        for r in &result {
            assert_eq!(r.position, 100);
            assert_eq!(r.ref_allele, b"A");
        }
    }

    #[test]
    fn join_biallelic_basic() {
        let v1 = snv("chr1", 100, "A", "G", Some(30.0));
        let v2 = snv("chr1", 100, "A", "T", Some(40.0));
        let joined = join_biallelic(&[v1, v2]).unwrap();
        assert_eq!(joined.alt_alleles.len(), 2);
        assert!(joined.alt_alleles.contains(&b"G".to_vec()));
        assert!(joined.alt_alleles.contains(&b"T".to_vec()));
    }

    #[test]
    fn join_biallelic_different_pos() {
        let v1 = snv("chr1", 100, "A", "G", None);
        let v2 = snv("chr1", 200, "A", "T", None);
        assert!(join_biallelic(&[v1, v2]).is_none());
    }

    #[test]
    fn join_biallelic_empty() {
        assert!(join_biallelic(&[]).is_none());
    }

    #[test]
    fn split_join_roundtrip() {
        let v = variant("chr1", 100, "A", &["G", "T"], Some(30.0), VariantFilter::Pass);
        let split = split_multiallelic(&v);
        let joined = join_biallelic(&split).unwrap();
        assert_eq!(joined.alt_alleles.len(), 2);
    }

    // -----------------------------------------------------------------------
    // Normalization
    // -----------------------------------------------------------------------

    #[test]
    fn normalize_trim_suffix() {
        // ACGT -> AGT should become AC -> A
        let v = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: None,
            ref_allele: b"ACGT".to_vec(),
            alt_alleles: vec![b"AGT".to_vec()],
            quality: None,
            filter: VariantFilter::Pass,
        };
        let n = normalize_variant(&v);
        assert_eq!(n.ref_allele, b"AC");
        assert_eq!(n.alt_alleles, vec![b"A".to_vec()]);
        assert_eq!(n.position, 100);
    }

    #[test]
    fn normalize_trim_prefix() {
        // TAC -> TA should become C -> (empty but we keep 1) -> stays as TAC->TA then trim suffix: TC->T, pos+0
        // Actually: common suffix first: TAC vs TA -> no common suffix (C vs A)
        // Common prefix: T is common -> remove T, pos becomes 101, AC vs A -> common suffix: C is not A, stop
        // Result: pos=101, ref=AC, alt=A
        let v = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: None,
            ref_allele: b"TAC".to_vec(),
            alt_alleles: vec![b"TA".to_vec()],
            quality: None,
            filter: VariantFilter::Pass,
        };
        let n = normalize_variant(&v);
        assert_eq!(n.ref_allele, b"AC");
        assert_eq!(n.alt_alleles, vec![b"A".to_vec()]);
        assert_eq!(n.position, 101);
    }

    #[test]
    fn normalize_no_change_for_snv() {
        let v = snv("chr1", 100, "A", "G", None);
        let n = normalize_variant(&v);
        assert_eq!(n.ref_allele, b"A");
        assert_eq!(n.alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(n.position, 100);
    }

    // -----------------------------------------------------------------------
    // Merging
    // -----------------------------------------------------------------------

    #[test]
    fn merge_overlapping() {
        let set1 = vec![snv("chr1", 100, "A", "G", Some(30.0))];
        let set2 = vec![snv("chr1", 100, "A", "T", Some(40.0))];
        let merged = merge_variants(&[&set1, &set2]);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].alt_alleles.len(), 2);
    }

    #[test]
    fn merge_non_overlapping() {
        let set1 = vec![snv("chr1", 100, "A", "G", Some(30.0))];
        let set2 = vec![snv("chr1", 200, "C", "T", Some(40.0))];
        let merged = merge_variants(&[&set1, &set2]);
        assert_eq!(merged.len(), 2);
    }

    // -----------------------------------------------------------------------
    // Filter expressions
    // -----------------------------------------------------------------------

    #[test]
    fn filter_qual_gt() {
        let expr = parse_filter("QUAL>30").unwrap();
        let v1 = snv("chr1", 100, "A", "G", Some(50.0));
        let v2 = snv("chr1", 200, "C", "T", Some(20.0));
        assert!(eval_filter(&v1, &expr));
        assert!(!eval_filter(&v2, &expr));
    }

    #[test]
    fn filter_qual_ge() {
        let expr = parse_filter("QUAL>=30").unwrap();
        let v = snv("chr1", 100, "A", "G", Some(30.0));
        assert!(eval_filter(&v, &expr));
    }

    #[test]
    fn filter_type_snv() {
        let expr = parse_filter("TYPE==SNV").unwrap();
        let snv_v = snv("chr1", 100, "A", "G", None);
        let indel = variant("chr1", 200, "AC", &["A"], None, VariantFilter::Pass);
        assert!(eval_filter(&snv_v, &expr));
        assert!(!eval_filter(&indel, &expr));
    }

    #[test]
    fn filter_type_indel() {
        let expr = parse_filter("TYPE==INDEL").unwrap();
        let indel = variant("chr1", 200, "AC", &["A"], None, VariantFilter::Pass);
        assert!(eval_filter(&indel, &expr));
    }

    #[test]
    fn filter_chrom() {
        let expr = parse_filter("CHROM==chr1").unwrap();
        let v1 = snv("chr1", 100, "A", "G", None);
        let v2 = snv("chr2", 100, "A", "G", None);
        assert!(eval_filter(&v1, &expr));
        assert!(!eval_filter(&v2, &expr));
    }

    #[test]
    fn filter_pass() {
        let expr = parse_filter("FILTER==PASS").unwrap();
        let pass_v = snv("chr1", 100, "A", "G", None);
        let fail_v = variant("chr1", 100, "A", &["G"], None, VariantFilter::Fail(vec!["LQ".into()]));
        assert!(eval_filter(&pass_v, &expr));
        assert!(!eval_filter(&fail_v, &expr));
    }

    #[test]
    fn filter_and() {
        let expr = parse_filter("QUAL>30 && CHROM==chr1").unwrap();
        let v1 = snv("chr1", 100, "A", "G", Some(50.0));
        let v2 = snv("chr2", 100, "A", "G", Some(50.0));
        let v3 = snv("chr1", 100, "A", "G", Some(10.0));
        assert!(eval_filter(&v1, &expr));
        assert!(!eval_filter(&v2, &expr));
        assert!(!eval_filter(&v3, &expr));
    }

    #[test]
    fn filter_or() {
        let expr = parse_filter("CHROM==chr1 || CHROM==chr2").unwrap();
        let v1 = snv("chr1", 100, "A", "G", None);
        let v2 = snv("chr2", 100, "A", "G", None);
        let v3 = snv("chr3", 100, "A", "G", None);
        assert!(eval_filter(&v1, &expr));
        assert!(eval_filter(&v2, &expr));
        assert!(!eval_filter(&v3, &expr));
    }

    #[test]
    fn filter_not() {
        let expr = parse_filter("!CHROM==chr1").unwrap();
        let v1 = snv("chr1", 100, "A", "G", None);
        let v2 = snv("chr2", 100, "A", "G", None);
        assert!(!eval_filter(&v1, &expr));
        assert!(eval_filter(&v2, &expr));
    }

    #[test]
    fn filter_parentheses() {
        let expr = parse_filter("(QUAL>30 || QUAL<10) && CHROM==chr1").unwrap();
        let v1 = snv("chr1", 100, "A", "G", Some(50.0));
        let v2 = snv("chr1", 100, "A", "G", Some(5.0));
        let v3 = snv("chr1", 100, "A", "G", Some(20.0));
        assert!(eval_filter(&v1, &expr));
        assert!(eval_filter(&v2, &expr));
        assert!(!eval_filter(&v3, &expr));
    }

    #[test]
    fn filter_variants_vec() {
        let expr = parse_filter("QUAL>30").unwrap();
        let variants = vec![
            snv("chr1", 100, "A", "G", Some(50.0)),
            snv("chr1", 200, "C", "T", Some(20.0)),
            snv("chr1", 300, "G", "A", Some(40.0)),
        ];
        let filtered = filter_variants(&variants, &expr);
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn filter_invalid_expr() {
        assert!(parse_filter("").is_err());
        assert!(parse_filter("UNKNOWN>5").is_err());
    }

    // -----------------------------------------------------------------------
    // Set operations
    // -----------------------------------------------------------------------

    #[test]
    fn intersection_basic() {
        let a = vec![
            snv("chr1", 100, "A", "G", None),
            snv("chr1", 200, "C", "T", None),
        ];
        let b = vec![
            snv("chr1", 100, "A", "G", None),
            snv("chr1", 300, "G", "A", None),
        ];
        let result = variant_intersection(&a, &b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].position, 100);
    }

    #[test]
    fn complement_basic() {
        let a = vec![
            snv("chr1", 100, "A", "G", None),
            snv("chr1", 200, "C", "T", None),
        ];
        let b = vec![snv("chr1", 100, "A", "G", None)];
        let result = variant_complement(&a, &b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].position, 200);
    }

    #[test]
    fn concordance_basic() {
        let a = vec![
            snv("chr1", 100, "A", "G", None),
            snv("chr1", 200, "C", "T", None),
        ];
        let b = vec![
            snv("chr1", 100, "A", "G", None),
            snv("chr1", 300, "G", "A", None),
        ];
        let stats = variant_concordance(&a, &b);
        assert_eq!(stats.shared, 1);
        assert_eq!(stats.a_only, 1);
        assert_eq!(stats.b_only, 1);
    }

    #[test]
    fn concordance_identical() {
        let a = vec![snv("chr1", 100, "A", "G", None)];
        let b = vec![snv("chr1", 100, "A", "G", None)];
        let stats = variant_concordance(&a, &b);
        assert_eq!(stats.shared, 1);
        assert_eq!(stats.a_only, 0);
        assert_eq!(stats.b_only, 0);
    }

    #[test]
    fn concordance_disjoint() {
        let a = vec![snv("chr1", 100, "A", "G", None)];
        let b = vec![snv("chr1", 200, "C", "T", None)];
        let stats = variant_concordance(&a, &b);
        assert_eq!(stats.shared, 0);
        assert_eq!(stats.a_only, 1);
        assert_eq!(stats.b_only, 1);
    }

    // -----------------------------------------------------------------------
    // Detailed stats
    // -----------------------------------------------------------------------

    #[test]
    fn detailed_stats_basic() {
        let variants = vec![
            snv("chr1", 100, "A", "G", Some(30.0)),        // transition
            snv("chr1", 200, "C", "T", Some(40.0)),        // transition
            snv("chr1", 300, "A", "C", Some(50.0)),        // transversion
            variant("chr1", 400, "AC", &["A"], Some(20.0), VariantFilter::Pass), // deletion
            variant("chr1", 500, "A", &["G", "T"], Some(60.0), VariantFilter::Pass), // multi-allelic SNV
        ];

        let stats = detailed_vcf_stats(&variants);
        assert_eq!(stats.total, 5);
        assert_eq!(stats.snv_count, 4); // positions 100, 200, 300, 500
        assert_eq!(stats.deletion_count, 1);
        assert_eq!(stats.pass_count, 5);
        assert_eq!(stats.multiallelic_count, 1);
        // Ti/Tv: 3 transitions (100: A>G, 200: C>T, 500: A>G), 1 transversion (300: A>C)
        assert!((stats.ti_tv_ratio - 3.0).abs() < 0.01);
    }

    #[test]
    fn detailed_stats_empty() {
        let stats = detailed_vcf_stats(&[]);
        assert_eq!(stats.total, 0);
        assert_eq!(stats.ti_tv_ratio, 0.0);
    }

    #[test]
    fn detailed_stats_quality_histogram() {
        let variants = vec![
            snv("chr1", 100, "A", "G", Some(5.0)),
            snv("chr1", 200, "C", "T", Some(15.0)),
            snv("chr1", 300, "G", "A", Some(15.0)),
            snv("chr1", 400, "T", "C", Some(105.0)),
        ];
        let stats = detailed_vcf_stats(&variants);
        // Bins: [0,10): 1, [10,20): 2, [100+]: 1
        assert_eq!(stats.quality_histogram.len(), 3);
    }

    #[test]
    fn detailed_stats_insertions() {
        let variants = vec![
            variant("chr1", 100, "A", &["AT"], None, VariantFilter::Pass),
            variant("chr1", 200, "A", &["ATG"], None, VariantFilter::Pass),
        ];
        let stats = detailed_vcf_stats(&variants);
        assert_eq!(stats.insertion_count, 2);
        assert_eq!(stats.deletion_count, 0);
    }
}
