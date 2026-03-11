//! mzTab output format writer.
//!
//! Produces mzTab 1.0 text output for peptide and protein identifications.
//! mzTab is a tab-delimited format standardized by the HUPO-PSI for
//! reporting proteomics results.

use crate::error::Result;
use crate::fdr::FdrResult;
use crate::protein::ProteinGroup;
use crate::quantification::ProteinQuant;

/// Write protein-level results in mzTab format.
///
/// Produces the PRH (protein header) and PRT (protein rows) sections.
pub fn write_mztab_proteins(
    groups: &[ProteinGroup],
    quants: Option<&[ProteinQuant]>,
) -> String {
    let mut out = String::new();

    // Metadata section
    out.push_str("MTD\tmzTab-version\t1.0.0\n");
    out.push_str("MTD\tmzTab-mode\tSummary\n");
    out.push_str("MTD\tmzTab-type\tQuantification\n");
    out.push_str("MTD\tdescription\tCyanea proteomics analysis\n");
    out.push_str("MTD\tprotein_search_engine_score[1]\t[MS, MS:1001171, Mascot:score, ]\n");
    out.push('\n');

    // Protein header
    out.push_str("PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tnum_psms\t");
    out.push_str("coverage\tbest_search_engine_score[1]\tnum_peptides_distinct\t");
    out.push_str("num_peptides_unique\n");

    // Protein rows
    for (i, group) in groups.iter().enumerate() {
        let accession = group.accessions.join(",");
        let num_peptides = group.unique_peptides.len() + group.shared_peptides.len();
        let unique_count = group.unique_peptides.len();

        // Get quant value if available
        let _quant_value = quants.and_then(|q| q.get(i)).map(|q| q.normalized_value);

        out.push_str(&format!(
            "PRT\t{}\tnull\tnull\tnull\tnull\t{}\t{:.4}\t{:.4}\t{}\t{}\n",
            accession,
            group.psm_count,
            group.coverage,
            group.best_score,
            num_peptides,
            unique_count,
        ));
    }

    out
}

/// Write PSM-level results in mzTab format.
///
/// Produces the PSH (PSM header) and PSM (PSM rows) sections.
pub fn write_mztab_psms(results: &[FdrResult]) -> String {
    let mut out = String::new();

    // PSM header
    out.push_str("PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\t");
    out.push_str("search_engine_score[1]\tspectra_ref\tcharge\texp_mass_to_charge\t");
    out.push_str("opt_q_value\topt_is_decoy\n");

    // PSM rows
    for (i, r) in results.iter().enumerate() {
        out.push_str(&format!(
            "PSM\t{}\t{}\tnull\tnull\tnull\t{:.4}\t{}\tnull\tnull\t{:.6}\t{}\n",
            r.peptide_sequence,
            i + 1,
            r.score,
            r.spectrum_id,
            r.q_value,
            if r.is_decoy { 1 } else { 0 },
        ));
    }

    out
}

/// Parse a simple mzTab protein section.
///
/// Returns protein accessions and PSM counts from PRT rows.
pub fn parse_mztab_proteins(text: &str) -> Result<Vec<(String, usize)>> {
    let mut proteins = Vec::new();

    for line in text.lines() {
        if !line.starts_with("PRT\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 7 {
            let accession = fields[1].to_string();
            let psm_count = fields[6].parse::<usize>().unwrap_or(0);
            proteins.push((accession, psm_count));
        }
    }

    Ok(proteins)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_mztab_proteins() {
        let groups = vec![
            ProteinGroup {
                accessions: vec!["P12345".to_string()],
                unique_peptides: vec!["AAAK".to_string(), "BBBK".to_string()],
                shared_peptides: vec!["CCCK".to_string()],
                psm_count: 5,
                best_score: 42.5,
                coverage: 0.35,
                is_decoy: false,
            },
        ];

        let output = write_mztab_proteins(&groups, None);
        assert!(output.contains("MTD\tmzTab-version\t1.0.0"));
        assert!(output.contains("PRH\taccession"));
        assert!(output.contains("PRT\tP12345"));
        assert!(output.contains("\t5\t")); // psm_count
        assert!(output.contains("0.3500")); // coverage
    }

    #[test]
    fn test_write_mztab_psms() {
        let results = vec![
            FdrResult {
                spectrum_id: "scan_1".to_string(),
                peptide_sequence: "PEPTIDEK".to_string(),
                score: 42.5,
                q_value: 0.001,
                is_decoy: false,
                passes: true,
            },
            FdrResult {
                spectrum_id: "scan_2".to_string(),
                peptide_sequence: "REVAAAK".to_string(),
                score: 15.0,
                q_value: 0.5,
                is_decoy: true,
                passes: false,
            },
        ];

        let output = write_mztab_psms(&results);
        assert!(output.contains("PSH\tsequence"));
        assert!(output.contains("PSM\tPEPTIDEK"));
        assert!(output.contains("\t0\n")); // not decoy
        assert!(output.contains("\t1\n")); // is decoy
    }

    #[test]
    fn test_roundtrip_proteins() {
        let groups = vec![
            ProteinGroup {
                accessions: vec!["P1".to_string()],
                unique_peptides: vec!["AAK".to_string()],
                shared_peptides: vec![],
                psm_count: 3,
                best_score: 20.0,
                coverage: 0.5,
                is_decoy: false,
            },
            ProteinGroup {
                accessions: vec!["P2".to_string()],
                unique_peptides: vec!["BBK".to_string()],
                shared_peptides: vec![],
                psm_count: 1,
                best_score: 10.0,
                coverage: 0.2,
                is_decoy: false,
            },
        ];

        let output = write_mztab_proteins(&groups, None);
        let parsed = parse_mztab_proteins(&output).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].0, "P1");
        assert_eq!(parsed[0].1, 3);
        assert_eq!(parsed[1].0, "P2");
        assert_eq!(parsed[1].1, 1);
    }
}
