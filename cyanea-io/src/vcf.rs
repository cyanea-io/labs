//! VCF (Variant Call Format) parser and writer.
//!
//! Parses VCF files into [`Variant`] records and writes variants to VCF format.
//! Supports VCF 4.x with the standard 8 mandatory columns
//! (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO).

use std::fs::File;
use std::io::{BufRead, BufReader, Write as IoWrite};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::variant::{Variant, VariantFilter};

/// Parse a VCF file and return all variant records.
///
/// Header lines (starting with `#`) are skipped. Each data line is parsed
/// into a [`Variant`] using VCF's 1-based coordinate system.
pub fn parse_vcf(path: impl AsRef<Path>) -> Result<Vec<Variant>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    // Read all data lines (skip headers)
    let mut data_lines: Vec<(usize, String)> = Vec::new();
    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let trimmed = line.trim().to_string();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        data_lines.push((line_num + 1, trimmed));
    }

    // Parse data lines (optionally in parallel)
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data_lines
            .par_iter()
            .map(|(line_num, line)| parse_vcf_line(line, *line_num, path))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    data_lines
        .iter()
        .map(|(line_num, line)| parse_vcf_line(line, *line_num, path))
        .collect()
}

/// Summary statistics for a VCF file.
#[derive(Debug, Clone)]
pub struct VcfStats {
    pub variant_count: u64,
    pub snv_count: u64,
    pub indel_count: u64,
    pub pass_count: u64,
    pub chromosomes: Vec<String>,
}

/// Parse a VCF file and return summary statistics without storing all records.
pub fn vcf_stats(path: impl AsRef<Path>) -> Result<VcfStats> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut variant_count: u64 = 0;
    let mut snv_count: u64 = 0;
    let mut indel_count: u64 = 0;
    let mut pass_count: u64 = 0;
    let mut chroms = Vec::new();
    let mut seen_chroms = std::collections::HashSet::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(CyaneaError::Io)?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let variant = parse_vcf_line(line, line_num + 1, path)?;
        variant_count += 1;

        if variant.is_snv() {
            snv_count += 1;
        }
        if variant.is_indel() {
            indel_count += 1;
        }
        if variant.filter == VariantFilter::Pass {
            pass_count += 1;
        }
        if seen_chroms.insert(variant.chrom.clone()) {
            chroms.push(variant.chrom);
        }
    }

    Ok(VcfStats {
        variant_count,
        snv_count,
        indel_count,
        pass_count,
        chromosomes: chroms,
    })
}

/// Parse a single VCF data line into a Variant.
fn parse_vcf_line(line: &str, line_num: usize, path: &Path) -> Result<Variant> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected at least 8 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let chrom = fields[0].to_string();

    let position: u64 = fields[1].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid position '{}'",
            path.display(),
            line_num,
            fields[1]
        ))
    })?;

    let id = match fields[2] {
        "." => None,
        s => Some(s.to_string()),
    };

    let ref_allele = fields[3].as_bytes().to_vec();

    let alt_alleles: Vec<Vec<u8>> = fields[4]
        .split(',')
        .filter(|a| *a != ".")
        .map(|a| a.as_bytes().to_vec())
        .collect();

    let quality: Option<f64> = match fields[5] {
        "." => None,
        s => Some(s.parse().map_err(|_| {
            CyaneaError::Parse(format!(
                "{}: line {}: invalid quality '{}'",
                path.display(),
                line_num,
                s
            ))
        })?),
    };

    let filter = match fields[6] {
        "PASS" => VariantFilter::Pass,
        "." => VariantFilter::Missing,
        s => VariantFilter::Fail(s.split(';').map(|f| f.to_string()).collect()),
    };

    Ok(Variant {
        chrom,
        position,
        id,
        ref_allele,
        alt_alleles,
        quality,
        filter,
    })
}

// ---------------------------------------------------------------------------
// VCF writing
// ---------------------------------------------------------------------------

/// Format a FILTER field for VCF output.
fn format_filter(filter: &VariantFilter) -> String {
    match filter {
        VariantFilter::Pass => "PASS".to_string(),
        VariantFilter::Missing => ".".to_string(),
        VariantFilter::Fail(reasons) => reasons.join(";"),
    }
}

/// Format variants as a VCF 4.3 string (8 mandatory columns, no genotype).
///
/// Includes file format header, column header, and one data line per variant.
pub fn write_vcf_string(variants: &[Variant]) -> String {
    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.3\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    for v in variants {
        let id = v.id.as_deref().unwrap_or(".");
        let ref_str = std::str::from_utf8(&v.ref_allele).unwrap_or(".");
        let alt_str: Vec<&str> = v
            .alt_alleles
            .iter()
            .map(|a| std::str::from_utf8(a).unwrap_or("."))
            .collect();
        let alt = if alt_str.is_empty() {
            ".".to_string()
        } else {
            alt_str.join(",")
        };
        let qual = v
            .quality
            .map(|q| format!("{:.1}", q))
            .unwrap_or_else(|| ".".to_string());
        let filter = format_filter(&v.filter);

        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\n",
            v.chrom, v.position, id, ref_str, alt, qual, filter,
        ));
    }

    out
}

/// Write variants to a VCF file (8 mandatory columns, no genotype).
pub fn write_vcf(variants: &[Variant], path: impl AsRef<Path>) -> Result<()> {
    let path = path.as_ref();
    let content = write_vcf_string(variants);
    let mut file = File::create(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    file.write_all(content.as_bytes()).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    Ok(())
}

/// Format called variants as a VCF 4.3 string with genotype columns.
///
/// FORMAT: `GT:DP:AD:GQ:PL`. INFO: `DP=N;AF=0.XX;SB=X.XXX`.
#[cfg(feature = "variant-calling")]
pub fn write_called_vcf_string(
    variants: &[crate::variant_call::CalledVariant],
    sample_name: &str,
) -> String {
    use crate::variant_call::Genotype;

    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.3\n");
    out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
    out.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
    out.push_str(
        "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias Fisher p-value\">\n",
    );
    out.push_str(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    );
    out.push_str(
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
    );
    out.push_str(
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depths\">\n",
    );
    out.push_str(
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n",
    );
    out.push_str(
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled Likelihoods\">\n",
    );
    out.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        sample_name
    ));

    for cv in variants {
        let v = &cv.variant;
        let id = v.id.as_deref().unwrap_or(".");
        let ref_str = std::str::from_utf8(&v.ref_allele).unwrap_or(".");
        let alt_str: Vec<&str> = v
            .alt_alleles
            .iter()
            .map(|a| std::str::from_utf8(a).unwrap_or("."))
            .collect();
        let alt = if alt_str.is_empty() {
            ".".to_string()
        } else {
            alt_str.join(",")
        };
        let qual = v
            .quality
            .map(|q| format!("{:.1}", q))
            .unwrap_or_else(|| ".".to_string());
        let filter = format_filter(&v.filter);

        let af = if cv.depth > 0 {
            cv.allele_depth[1] as f64 / cv.depth as f64
        } else {
            0.0
        };
        let info = format!("DP={};AF={:.3};SB={:.3}", cv.depth, af, cv.strand_bias_p);

        let gt = match cv.genotype {
            Genotype::HomRef => "0/0",
            Genotype::Het => "0/1",
            Genotype::HomAlt => "1/1",
        };
        let gq = cv.genotype_quality.round() as u32;
        let sample = format!(
            "{}:{}:{},{}:{}:{},{},{}",
            gt,
            cv.depth,
            cv.allele_depth[0],
            cv.allele_depth[1],
            gq,
            cv.genotype_likelihoods[0],
            cv.genotype_likelihoods[1],
            cv.genotype_likelihoods[2],
        );

        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT:DP:AD:GQ:PL\t{}\n",
            v.chrom, v.position, id, ref_str, alt, qual, filter, info, sample,
        ));
    }

    out
}

/// Write called variants to a VCF file with genotype columns.
#[cfg(feature = "variant-calling")]
pub fn write_called_vcf(
    variants: &[crate::variant_call::CalledVariant],
    sample_name: &str,
    path: impl AsRef<Path>,
) -> Result<()> {
    let path = path.as_ref();
    let content = write_called_vcf_string(variants, sample_name);
    let mut file = File::create(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    file.write_all(content.as_bytes()).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_vcf(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".vcf").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_parse_basic_vcf() {
        let file = write_vcf(
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\trs123\tA\tG\t30.0\tPASS\t.\n\
             chr1\t200\t.\tAC\tA\t.\t.\t.\n\
             chr2\t300\t.\tT\tTA,TG\t50.5\tLowQual\t.\n",
        );

        let variants = parse_vcf(file.path()).unwrap();
        assert_eq!(variants.len(), 3);

        // SNV
        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].position, 100);
        assert_eq!(variants[0].id, Some("rs123".to_string()));
        assert_eq!(variants[0].ref_allele, b"A");
        assert_eq!(variants[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(variants[0].quality, Some(30.0));
        assert_eq!(variants[0].filter, VariantFilter::Pass);
        assert!(variants[0].is_snv());

        // Deletion
        assert_eq!(variants[1].ref_allele, b"AC");
        assert_eq!(variants[1].alt_alleles, vec![b"A".to_vec()]);
        assert!(variants[1].is_indel());
        assert_eq!(variants[1].quality, None);
        assert_eq!(variants[1].filter, VariantFilter::Missing);

        // Multi-allelic
        assert_eq!(variants[2].alt_alleles.len(), 2);
        assert_eq!(variants[2].filter, VariantFilter::Fail(vec!["LowQual".to_string()]));
    }

    #[test]
    fn test_vcf_stats() {
        let file = write_vcf(
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG\t30\tPASS\t.\n\
             chr1\t200\t.\tAC\tA\t.\tPASS\t.\n\
             chr2\t300\t.\tT\tC\t50\tLowQual\t.\n",
        );

        let stats = vcf_stats(file.path()).unwrap();
        assert_eq!(stats.variant_count, 3);
        assert_eq!(stats.snv_count, 2);
        assert_eq!(stats.indel_count, 1);
        assert_eq!(stats.pass_count, 2);
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn test_vcf_empty_file() {
        let file = write_vcf("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        let variants = parse_vcf(file.path()).unwrap();
        assert!(variants.is_empty());
    }

    #[test]
    fn test_vcf_missing_columns() {
        let file = write_vcf("#CHROM\tPOS\n\nchr1\t100\n");
        let result = parse_vcf(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_vcf_file_not_found() {
        let result = parse_vcf("/nonexistent/file.vcf");
        assert!(result.is_err());
    }

    // -- VCF writing tests --

    #[test]
    fn test_write_vcf_string_basic() {
        let v1 = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: Some("rs123".to_string()),
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            quality: Some(30.0),
            filter: VariantFilter::Pass,
        };
        let v2 = Variant {
            chrom: "chr1".to_string(),
            position: 200,
            id: None,
            ref_allele: b"AC".to_vec(),
            alt_alleles: vec![b"A".to_vec()],
            quality: None,
            filter: VariantFilter::Missing,
        };

        let output = write_vcf_string(&[v1, v2]);

        assert!(output.starts_with("##fileformat=VCFv4.3\n"));
        assert!(output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
        assert!(output.contains("chr1\t100\trs123\tA\tG\t30.0\tPASS\t.\n"));
        assert!(output.contains("chr1\t200\t.\tAC\tA\t.\t.\t.\n"));
    }

    #[test]
    fn test_write_vcf_header() {
        let output = write_vcf_string(&[]);
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.3");
        assert_eq!(lines[1], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    }

    #[test]
    fn test_write_vcf_file_roundtrip() {
        let v = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: None,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            quality: Some(42.0),
            filter: VariantFilter::Pass,
        };

        let tmp = NamedTempFile::with_suffix(".vcf").unwrap();
        super::write_vcf(&[v], tmp.path()).unwrap();

        let parsed = parse_vcf(tmp.path()).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].chrom, "chr1");
        assert_eq!(parsed[0].position, 100);
        assert_eq!(parsed[0].ref_allele, b"A");
        assert_eq!(parsed[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(parsed[0].filter, VariantFilter::Pass);
    }

    #[test]
    fn test_write_vcf_multi_allelic() {
        let v = Variant {
            chrom: "chr2".to_string(),
            position: 300,
            id: None,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"TA".to_vec(), b"TG".to_vec()],
            quality: Some(50.5),
            filter: VariantFilter::Fail(vec!["LowQual".to_string()]),
        };

        let output = write_vcf_string(&[v]);
        assert!(output.contains("chr2\t300\t.\tT\tTA,TG\t50.5\tLowQual\t.\n"));
    }

    #[cfg(feature = "variant-calling")]
    #[test]
    fn test_write_called_vcf_string() {
        use crate::variant_call::{CalledVariant, Genotype};

        let variant = Variant {
            chrom: "chr1".to_string(),
            position: 100,
            id: None,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            quality: Some(45.0),
            filter: VariantFilter::Pass,
        };

        let cv = CalledVariant {
            variant,
            genotype: Genotype::Het,
            genotype_quality: 35.0,
            depth: 20,
            allele_depth: [10, 10],
            allele_depth_fwd: [5, 5],
            allele_depth_rev: [5, 5],
            strand_bias_p: 1.0,
            genotype_likelihoods: [150, 0, 200],
        };

        let output = write_called_vcf_string(&[cv], "SAMPLE1");

        assert!(output.contains("##fileformat=VCFv4.3\n"));
        assert!(output.contains("##FORMAT=<ID=GT"));
        assert!(output.contains("SAMPLE1\n"));
        assert!(output.contains("GT:DP:AD:GQ:PL\t0/1:20:10,10:35:150,0,200\n"));
        assert!(output.contains("DP=20;AF=0.500;SB=1.000"));
    }
}
