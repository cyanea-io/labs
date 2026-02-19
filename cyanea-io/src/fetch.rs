//! Database URL builders and response parsers for common bioinformatics APIs.
//!
//! This module provides URL construction helpers for NCBI Entrez, UniProt,
//! KEGG, htsget, and refget APIs, along with lightweight response parsers.
//! No HTTP client is included — these are pure URL builders and text parsers.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Percent-encoding helper
// ---------------------------------------------------------------------------

/// Percent-encode common special characters in a query string value.
fn percent_encode(input: &str) -> String {
    let mut out = String::with_capacity(input.len() * 2);
    for b in input.bytes() {
        match b {
            b' ' => out.push_str("%20"),
            b'&' => out.push_str("%26"),
            b'=' => out.push_str("%3D"),
            b'+' => out.push_str("%2B"),
            b'#' => out.push_str("%23"),
            _ => out.push(b as char),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// NCBI Entrez
// ---------------------------------------------------------------------------

/// URL builder for NCBI Entrez E-utilities.
pub struct EntrezUrl;

impl EntrezUrl {
    const BASE: &'static str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

    /// Build an esearch URL.
    ///
    /// ```
    /// use cyanea_io::fetch::EntrezUrl;
    /// let url = EntrezUrl::esearch("pubmed", "cancer", 10);
    /// assert!(url.starts_with("https://eutils.ncbi.nlm.nih.gov"));
    /// ```
    pub fn esearch(db: &str, term: &str, retmax: usize) -> String {
        format!(
            "{}/esearch.fcgi?db={}&term={}&retmax={}&retmode=xml",
            Self::BASE,
            percent_encode(db),
            percent_encode(term),
            retmax,
        )
    }

    /// Build an efetch URL with one or more IDs.
    pub fn efetch(db: &str, ids: &[&str], rettype: &str, retmode: &str) -> String {
        let id_list = ids.join(",");
        format!(
            "{}/efetch.fcgi?db={}&id={}&rettype={}&retmode={}",
            Self::BASE,
            percent_encode(db),
            percent_encode(&id_list),
            percent_encode(rettype),
            percent_encode(retmode),
        )
    }

    /// Build an elink URL to connect IDs between databases.
    pub fn elink(db_from: &str, db_to: &str, ids: &[&str]) -> String {
        let id_list = ids.join(",");
        format!(
            "{}/elink.fcgi?dbfrom={}&db={}&id={}",
            Self::BASE,
            percent_encode(db_from),
            percent_encode(db_to),
            percent_encode(&id_list),
        )
    }
}

/// Parse an NCBI esearch XML response, extracting IDs from `<IdList>`.
///
/// Finds all `<Id>...</Id>` tags and returns their text content.
pub fn parse_esearch_response(xml: &str) -> Result<Vec<String>> {
    let mut ids = Vec::new();
    let mut rest = xml;
    while let Some(start_pos) = rest.find("<Id>") {
        let after_tag = &rest[start_pos + 4..];
        let end_pos = after_tag.find("</Id>").ok_or_else(|| {
            CyaneaError::Parse("esearch XML: found <Id> without matching </Id>".into())
        })?;
        let id_text = after_tag[..end_pos].trim();
        if !id_text.is_empty() {
            ids.push(id_text.to_string());
        }
        rest = &after_tag[end_pos + 5..];
    }
    Ok(ids)
}

/// Parse an NCBI efetch FASTA response into (header, sequence) pairs.
///
/// Each entry starts with `>`, the first line is the header, and
/// subsequent lines (until the next `>` or EOF) are concatenated as the
/// sequence bytes.
pub fn parse_efetch_fasta(response: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let mut records = Vec::new();
    for entry in response.split('>') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let mut lines = entry.lines();
        let header = lines
            .next()
            .ok_or_else(|| CyaneaError::Parse("FASTA entry has no header line".into()))?
            .to_string();
        let mut seq = Vec::new();
        for line in lines {
            seq.extend_from_slice(line.trim().as_bytes());
        }
        records.push((header, seq));
    }
    if records.is_empty() {
        return Err(CyaneaError::Parse("no FASTA records found in response".into()));
    }
    Ok(records)
}

// ---------------------------------------------------------------------------
// UniProt
// ---------------------------------------------------------------------------

/// URL builder for the UniProt REST API.
pub struct UniprotUrl;

impl UniprotUrl {
    const BASE: &'static str = "https://rest.uniprot.org/uniprotkb";

    /// Build an entry retrieval URL.
    ///
    /// ```
    /// use cyanea_io::fetch::UniprotUrl;
    /// let url = UniprotUrl::entry("P12345", "fasta");
    /// assert!(url.contains("P12345.fasta"));
    /// ```
    pub fn entry(accession: &str, format: &str) -> String {
        format!("{}/{}.{}", Self::BASE, accession, format)
    }

    /// Build a search URL.
    pub fn search(query: &str, format: &str, size: usize) -> String {
        format!(
            "{}/search?query={}&format={}&size={}",
            Self::BASE,
            percent_encode(query),
            percent_encode(format),
            size,
        )
    }
}

// ---------------------------------------------------------------------------
// KEGG
// ---------------------------------------------------------------------------

/// URL builder for the KEGG REST API.
pub struct KeggUrl;

impl KeggUrl {
    const BASE: &'static str = "https://rest.kegg.jp";

    /// Build a URL to retrieve a pathway entry.
    pub fn get_pathway(pathway_id: &str) -> String {
        format!("{}/get/{}", Self::BASE, percent_encode(pathway_id))
    }

    /// Build a search/find URL.
    pub fn find(database: &str, query: &str) -> String {
        format!(
            "{}/find/{}/{}",
            Self::BASE,
            percent_encode(database),
            percent_encode(query),
        )
    }

    /// Build a URL to retrieve a specific database entry.
    pub fn get(database: &str, entry_id: &str) -> String {
        format!(
            "{}/get/{}:{}",
            Self::BASE,
            percent_encode(database),
            percent_encode(entry_id),
        )
    }
}

// ---------------------------------------------------------------------------
// htsget
// ---------------------------------------------------------------------------

/// URL builder for the htsget protocol (GA4GH).
pub struct HtsgetUrl;

impl HtsgetUrl {
    /// Build a reads URL with optional region and format filters.
    ///
    /// Only query parameters whose corresponding `Option` is `Some` are
    /// included in the URL.
    pub fn reads(
        base_url: &str,
        id: &str,
        reference_name: Option<&str>,
        start: Option<u64>,
        end: Option<u64>,
        format: Option<&str>,
    ) -> String {
        let base = format!("{}/reads/{}", base_url.trim_end_matches('/'), id);
        let mut params = Vec::new();
        if let Some(r) = reference_name {
            params.push(format!("referenceName={}", percent_encode(r)));
        }
        if let Some(s) = start {
            params.push(format!("start={}", s));
        }
        if let Some(e) = end {
            params.push(format!("end={}", e));
        }
        if let Some(f) = format {
            params.push(format!("format={}", percent_encode(f)));
        }
        if params.is_empty() {
            base
        } else {
            format!("{}?{}", base, params.join("&"))
        }
    }

    /// Build a variants URL with optional region filters.
    pub fn variants(
        base_url: &str,
        id: &str,
        reference_name: Option<&str>,
        start: Option<u64>,
        end: Option<u64>,
    ) -> String {
        let base = format!("{}/variants/{}", base_url.trim_end_matches('/'), id);
        let mut params = Vec::new();
        if let Some(r) = reference_name {
            params.push(format!("referenceName={}", percent_encode(r)));
        }
        if let Some(s) = start {
            params.push(format!("start={}", s));
        }
        if let Some(e) = end {
            params.push(format!("end={}", e));
        }
        if params.is_empty() {
            base
        } else {
            format!("{}?{}", base, params.join("&"))
        }
    }
}

// ---------------------------------------------------------------------------
// refget
// ---------------------------------------------------------------------------

/// URL builder for the refget protocol (GA4GH).
pub struct RefgetUrl;

impl RefgetUrl {
    /// Build a sequence retrieval URL with optional byte-range parameters.
    pub fn sequence(
        base_url: &str,
        sequence_id: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> String {
        let base = format!(
            "{}/sequence/{}",
            base_url.trim_end_matches('/'),
            sequence_id,
        );
        let mut params = Vec::new();
        if let Some(s) = start {
            params.push(format!("start={}", s));
        }
        if let Some(e) = end {
            params.push(format!("end={}", e));
        }
        if params.is_empty() {
            base
        } else {
            format!("{}?{}", base, params.join("&"))
        }
    }

    /// Build a sequence metadata URL.
    pub fn metadata(base_url: &str, sequence_id: &str) -> String {
        format!(
            "{}/sequence/{}/metadata",
            base_url.trim_end_matches('/'),
            sequence_id,
        )
    }
}

/// Metadata for a refget sequence.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RefgetMetadata {
    /// MD5 checksum of the sequence.
    pub md5: String,
    /// Length of the sequence in bases.
    pub length: u64,
    /// Known aliases for this sequence.
    pub aliases: Vec<RefgetAlias>,
}

/// A single alias entry in refget metadata.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RefgetAlias {
    /// The naming authority (e.g., "insdc", "ensembl").
    pub naming_authority: String,
    /// The alias identifier.
    pub alias: String,
}

/// Parse a refget metadata JSON response into a [`RefgetMetadata`].
///
/// Performs simple manual JSON parsing (no serde dependency at runtime).
/// Expected input shape:
/// ```json
/// {"metadata":{"md5":"abc","length":1000,"aliases":[{"naming_authority":"insdc","alias":"chr1"}]}}
/// ```
pub fn parse_refget_metadata(json: &str) -> Result<RefgetMetadata> {
    // Helper: extract a quoted string value for the given key from a substring.
    fn extract_string<'a>(src: &'a str, key: &str) -> std::result::Result<&'a str, CyaneaError> {
        let needle = format!("\"{}\":\"", key);
        let pos = src.find(&needle).ok_or_else(|| {
            CyaneaError::Parse(format!("refget metadata: missing key \"{}\"", key))
        })?;
        let after = &src[pos + needle.len()..];
        let end = after.find('"').ok_or_else(|| {
            CyaneaError::Parse(format!(
                "refget metadata: unterminated string for key \"{}\"",
                key
            ))
        })?;
        Ok(&after[..end])
    }

    // Helper: extract a numeric value for the given key.
    fn extract_u64(src: &str, key: &str) -> std::result::Result<u64, CyaneaError> {
        let needle = format!("\"{}\":", key);
        let pos = src.find(&needle).ok_or_else(|| {
            CyaneaError::Parse(format!("refget metadata: missing key \"{}\"", key))
        })?;
        let after = &src[pos + needle.len()..];
        let after = after.trim_start();
        let end = after
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(after.len());
        after[..end].parse::<u64>().map_err(|e| {
            CyaneaError::Parse(format!(
                "refget metadata: invalid number for \"{}\": {}",
                key, e
            ))
        })
    }

    // Find the metadata object.
    let meta_start = json.find("\"metadata\"").ok_or_else(|| {
        CyaneaError::Parse("refget metadata: missing \"metadata\" key".into())
    })?;
    let meta_body = &json[meta_start..];

    let md5 = extract_string(meta_body, "md5")?.to_string();
    let length = extract_u64(meta_body, "length")?;

    // Parse aliases array.
    let mut aliases = Vec::new();
    if let Some(arr_start) = meta_body.find("\"aliases\"") {
        let after_key = &meta_body[arr_start..];
        if let Some(bracket) = after_key.find('[') {
            let arr_body = &after_key[bracket..];
            // Find each object in the array.
            let mut rest = arr_body;
            while let Some(obj_start) = rest.find('{') {
                let obj_body = &rest[obj_start..];
                let obj_end = obj_body.find('}').ok_or_else(|| {
                    CyaneaError::Parse(
                        "refget metadata: unterminated alias object".into(),
                    )
                })?;
                let obj_str = &obj_body[..=obj_end];
                let na = extract_string(obj_str, "naming_authority")?.to_string();
                let al = extract_string(obj_str, "alias")?.to_string();
                aliases.push(RefgetAlias {
                    naming_authority: na,
                    alias: al,
                });
                rest = &obj_body[obj_end + 1..];
            }
        }
    }

    Ok(RefgetMetadata {
        md5,
        length,
        aliases,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_esearch_url() {
        let url = EntrezUrl::esearch("pubmed", "BRCA1 cancer", 20);
        assert!(url.starts_with("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"));
        assert!(url.contains("db=pubmed"));
        assert!(url.contains("term=BRCA1%20cancer"));
        assert!(url.contains("retmax=20"));
        assert!(url.contains("retmode=xml"));
    }

    #[test]
    fn test_efetch_multi_id() {
        let url = EntrezUrl::efetch("nucleotide", &["12345", "67890", "11111"], "fasta", "text");
        assert!(url.contains("db=nucleotide"));
        assert!(url.contains("id=12345,67890,11111"));
        assert!(url.contains("rettype=fasta"));
        assert!(url.contains("retmode=text"));
    }

    #[test]
    fn test_elink_url() {
        let url = EntrezUrl::elink("pubmed", "pmc", &["9001", "9002"]);
        assert!(url.starts_with("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"));
        assert!(url.contains("dbfrom=pubmed"));
        assert!(url.contains("db=pmc"));
        assert!(url.contains("id=9001,9002"));
    }

    #[test]
    fn test_parse_esearch_response() {
        let xml = r#"<?xml version="1.0"?>
<eSearchResult>
  <Count>3</Count>
  <IdList>
    <Id>39854103</Id>
    <Id>39851825</Id>
    <Id>39437063</Id>
  </IdList>
</eSearchResult>"#;
        let ids = parse_esearch_response(xml).unwrap();
        assert_eq!(ids, vec!["39854103", "39851825", "39437063"]);
    }

    #[test]
    fn test_parse_efetch_fasta() {
        let fasta = ">seq1 human gene\nATCGATCG\nGCTAGCTA\n>seq2 mouse gene\nTTTTAAAA\n";
        let records = parse_efetch_fasta(fasta).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "seq1 human gene");
        assert_eq!(records[0].1, b"ATCGATCGGCTAGCTA");
        assert_eq!(records[1].0, "seq2 mouse gene");
        assert_eq!(records[1].1, b"TTTTAAAA");
    }

    #[test]
    fn test_uniprot_entry_url() {
        let url = UniprotUrl::entry("P12345", "fasta");
        assert_eq!(url, "https://rest.uniprot.org/uniprotkb/P12345.fasta");
    }

    #[test]
    fn test_uniprot_search_url() {
        let url = UniprotUrl::search("organism_id:9606 AND gene:TP53", "json", 25);
        assert!(url.starts_with("https://rest.uniprot.org/uniprotkb/search"));
        assert!(url.contains("query=organism_id:9606%20AND%20gene:TP53"));
        assert!(url.contains("format=json"));
        assert!(url.contains("size=25"));
    }

    #[test]
    fn test_kegg_pathway_url() {
        let url = KeggUrl::get_pathway("hsa00010");
        assert_eq!(url, "https://rest.kegg.jp/get/hsa00010");
    }

    #[test]
    fn test_htsget_reads_with_region() {
        let url = HtsgetUrl::reads(
            "https://htsget.example.org",
            "NA12878",
            Some("chr1"),
            Some(10000),
            Some(20000),
            Some("BAM"),
        );
        assert!(url.starts_with("https://htsget.example.org/reads/NA12878"));
        assert!(url.contains("referenceName=chr1"));
        assert!(url.contains("start=10000"));
        assert!(url.contains("end=20000"));
        assert!(url.contains("format=BAM"));
    }

    #[test]
    fn test_htsget_reads_no_region() {
        let url = HtsgetUrl::reads(
            "https://htsget.example.org",
            "NA12878",
            None,
            None,
            None,
            None,
        );
        assert_eq!(url, "https://htsget.example.org/reads/NA12878");
        assert!(!url.contains('?'));
    }

    #[test]
    fn test_refget_sequence_with_range() {
        let url = RefgetUrl::sequence(
            "https://refget.example.org",
            "abc123def456",
            Some(100),
            Some(200),
        );
        assert!(url.starts_with("https://refget.example.org/sequence/abc123def456"));
        assert!(url.contains("start=100"));
        assert!(url.contains("end=200"));
    }

    #[test]
    fn test_refget_metadata_url() {
        let url = RefgetUrl::metadata("https://refget.example.org", "abc123def456");
        assert_eq!(
            url,
            "https://refget.example.org/sequence/abc123def456/metadata"
        );
    }
}
