//! BLAST XML output parser.
//!
//! Parses BLAST XML output (the `-outfmt 5` format produced by `blastn`,
//! `blastp`, `blastx`, etc.). Uses a simple state-machine approach with
//! tag extraction — no external XML library dependency required.

use cyanea_core::{CyaneaError, Result};

/// Top-level result from a BLAST XML output file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BlastXmlResult {
    /// BLAST program name (e.g. "blastn").
    pub program: String,
    /// Program version string.
    pub version: String,
    /// Database searched.
    pub db: String,
    /// Query sequence identifier.
    pub query_id: String,
    /// Query sequence length.
    pub query_len: u64,
    /// Search iterations (one per query in multi-query runs).
    pub iterations: Vec<BlastXmlIteration>,
}

/// A single search iteration within BLAST XML output.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BlastXmlIteration {
    /// Iteration number (1-based).
    pub iteration_num: u32,
    /// Hits found in this iteration.
    pub hits: Vec<BlastXmlHit>,
}

/// A single database hit.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BlastXmlHit {
    /// Hit identifier.
    pub hit_id: String,
    /// Hit definition/description.
    pub hit_def: String,
    /// Hit accession number.
    pub hit_accession: String,
    /// Length of the hit sequence.
    pub hit_len: u64,
    /// High-scoring segment pairs for this hit.
    pub hsps: Vec<BlastXmlHsp>,
}

/// A single high-scoring segment pair (HSP).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BlastXmlHsp {
    /// Bit score.
    pub bit_score: f64,
    /// Expect value.
    pub evalue: f64,
    /// Start of alignment in query (1-based).
    pub query_from: u64,
    /// End of alignment in query (1-based).
    pub query_to: u64,
    /// Start of alignment in hit (1-based).
    pub hit_from: u64,
    /// End of alignment in hit (1-based).
    pub hit_to: u64,
    /// Number of identical positions.
    pub identity: u64,
    /// Number of gap positions.
    pub gaps: u64,
    /// Alignment length.
    pub align_len: u64,
    /// Aligned query sequence.
    pub qseq: String,
    /// Aligned hit sequence.
    pub hseq: String,
    /// Midline (match/mismatch indicators).
    pub midline: String,
}

/// Extract the text content between `<tag>` and `</tag>`.
///
/// Returns `None` if the tag is not found.
fn extract_tag_content(xml: &str, tag: &str) -> Option<String> {
    let open = format!("<{}>", tag);
    let close = format!("</{}>", tag);
    let start = xml.find(&open)?;
    let content_start = start + open.len();
    let end = xml[content_start..].find(&close)?;
    Some(xml[content_start..content_start + end].to_string())
}

/// Extract all blocks delimited by `<tag>...</tag>`.
///
/// Returns the inner content of each block (excluding the outer tags).
fn extract_all_blocks(xml: &str, tag: &str) -> Vec<String> {
    let open = format!("<{}>", tag);
    let close = format!("</{}>", tag);
    let mut blocks = Vec::new();
    let mut search_from = 0;

    while search_from < xml.len() {
        let start = match xml[search_from..].find(&open) {
            Some(pos) => search_from + pos,
            None => break,
        };
        let content_start = start + open.len();
        let end = match xml[content_start..].find(&close) {
            Some(pos) => content_start + pos,
            None => break,
        };
        blocks.push(xml[content_start..end].to_string());
        search_from = end + close.len();
    }

    blocks
}

/// Parse a BLAST XML string into a [`BlastXmlResult`].
///
/// Expects the full `<BlastOutput>` XML document as produced by BLAST
/// with `-outfmt 5`. Uses simple tag extraction rather than a full XML
/// parser, so the input must be well-formed BLAST XML.
///
/// # Errors
///
/// Returns an error if required fields are missing or numeric values
/// cannot be parsed.
pub fn parse_blast_xml_str(xml: &str) -> Result<BlastXmlResult> {
    // Verify this looks like BLAST XML
    if !xml.contains("<BlastOutput>") {
        return Err(CyaneaError::Parse(
            "not a valid BLAST XML document: missing <BlastOutput>".to_string(),
        ));
    }

    let program = extract_tag_content(xml, "BlastOutput_program").unwrap_or_default();
    let version = extract_tag_content(xml, "BlastOutput_version").unwrap_or_default();
    let db = extract_tag_content(xml, "BlastOutput_db").unwrap_or_default();
    let query_id = extract_tag_content(xml, "BlastOutput_query-ID").unwrap_or_default();

    let query_len_str = extract_tag_content(xml, "BlastOutput_query-len").ok_or_else(|| {
        CyaneaError::Parse("missing <BlastOutput_query-len>".to_string())
    })?;
    let query_len: u64 = query_len_str.trim().parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "invalid query length: '{}'",
            query_len_str.trim()
        ))
    })?;

    // Parse iterations
    let iteration_blocks = extract_all_blocks(xml, "Iteration");
    let mut iterations = Vec::with_capacity(iteration_blocks.len());

    for iter_xml in &iteration_blocks {
        let iter_num_str =
            extract_tag_content(iter_xml, "Iteration_iter-num").unwrap_or_default();
        let iteration_num: u32 = iter_num_str.trim().parse().map_err(|_| {
            CyaneaError::Parse(format!(
                "invalid iteration number: '{}'",
                iter_num_str.trim()
            ))
        })?;

        // Parse hits within this iteration
        let hit_blocks = extract_all_blocks(iter_xml, "Hit");
        let mut hits = Vec::with_capacity(hit_blocks.len());

        for hit_xml in &hit_blocks {
            let hit = parse_hit(hit_xml)?;
            hits.push(hit);
        }

        iterations.push(BlastXmlIteration {
            iteration_num,
            hits,
        });
    }

    Ok(BlastXmlResult {
        program,
        version,
        db,
        query_id,
        query_len,
        iterations,
    })
}

/// Parse a single `<Hit>` block.
fn parse_hit(hit_xml: &str) -> Result<BlastXmlHit> {
    let hit_id = extract_tag_content(hit_xml, "Hit_id").unwrap_or_default();
    let hit_def = extract_tag_content(hit_xml, "Hit_def").unwrap_or_default();
    let hit_accession = extract_tag_content(hit_xml, "Hit_accession").unwrap_or_default();

    let hit_len_str = extract_tag_content(hit_xml, "Hit_len").ok_or_else(|| {
        CyaneaError::Parse("missing <Hit_len>".to_string())
    })?;
    let hit_len: u64 = hit_len_str.trim().parse().map_err(|_| {
        CyaneaError::Parse(format!("invalid hit length: '{}'", hit_len_str.trim()))
    })?;

    let hsp_blocks = extract_all_blocks(hit_xml, "Hsp");
    let mut hsps = Vec::with_capacity(hsp_blocks.len());

    for hsp_xml in &hsp_blocks {
        let hsp = parse_hsp(hsp_xml)?;
        hsps.push(hsp);
    }

    Ok(BlastXmlHit {
        hit_id,
        hit_def,
        hit_accession,
        hit_len,
        hsps,
    })
}

/// Parse a single `<Hsp>` block.
fn parse_hsp(hsp_xml: &str) -> Result<BlastXmlHsp> {
    let parse_f64 = |tag: &str| -> Result<f64> {
        let s = extract_tag_content(hsp_xml, tag)
            .ok_or_else(|| CyaneaError::Parse(format!("missing <{}>", tag)))?;
        s.trim()
            .parse::<f64>()
            .map_err(|_| CyaneaError::Parse(format!("invalid {}: '{}'", tag, s.trim())))
    };

    let parse_u64 = |tag: &str| -> Result<u64> {
        let s = extract_tag_content(hsp_xml, tag)
            .ok_or_else(|| CyaneaError::Parse(format!("missing <{}>", tag)))?;
        s.trim()
            .parse::<u64>()
            .map_err(|_| CyaneaError::Parse(format!("invalid {}: '{}'", tag, s.trim())))
    };

    Ok(BlastXmlHsp {
        bit_score: parse_f64("Hsp_bit-score")?,
        evalue: parse_f64("Hsp_evalue")?,
        query_from: parse_u64("Hsp_query-from")?,
        query_to: parse_u64("Hsp_query-to")?,
        hit_from: parse_u64("Hsp_hit-from")?,
        hit_to: parse_u64("Hsp_hit-to")?,
        identity: parse_u64("Hsp_identity")?,
        gaps: parse_u64("Hsp_gaps")?,
        align_len: parse_u64("Hsp_align-len")?,
        qseq: extract_tag_content(hsp_xml, "Hsp_qseq").unwrap_or_default(),
        hseq: extract_tag_content(hsp_xml, "Hsp_hseq").unwrap_or_default(),
        midline: extract_tag_content(hsp_xml, "Hsp_midline").unwrap_or_default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_BLAST_XML: &str = r#"<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.12.0+</BlastOutput_version>
  <BlastOutput_db>nt</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_hits>
        <Hit>
          <Hit_id>ref|NM_001</Hit_id>
          <Hit_def>Example gene</Hit_def>
          <Hit_accession>NM_001</Hit_accession>
          <Hit_len>500</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>150.0</Hsp_bit-score>
              <Hsp_evalue>1e-40</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>50</Hsp_hit-from>
              <Hsp_hit-to>149</Hsp_hit-to>
              <Hsp_identity>95</Hsp_identity>
              <Hsp_gaps>2</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>ATGCATGC</Hsp_qseq>
              <Hsp_hseq>ATGNATGC</Hsp_hseq>
              <Hsp_midline>||| ||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"#;

    #[test]
    fn blast_xml_minimal_one_hit() {
        let result = parse_blast_xml_str(MINIMAL_BLAST_XML).unwrap();
        assert_eq!(result.program, "blastn");
        assert_eq!(result.version, "BLASTN 2.12.0+");
        assert_eq!(result.db, "nt");
        assert_eq!(result.query_id, "Query_1");
        assert_eq!(result.query_len, 100);
        assert_eq!(result.iterations.len(), 1);
        assert_eq!(result.iterations[0].iteration_num, 1);
        assert_eq!(result.iterations[0].hits.len(), 1);

        let hit = &result.iterations[0].hits[0];
        assert_eq!(hit.hit_id, "ref|NM_001");
        assert_eq!(hit.hit_def, "Example gene");
        assert_eq!(hit.hit_accession, "NM_001");
        assert_eq!(hit.hit_len, 500);
        assert_eq!(hit.hsps.len(), 1);
    }

    #[test]
    fn blast_xml_multi_hit() {
        let xml = r#"<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.12.0+</BlastOutput_version>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_2</BlastOutput_query-ID>
  <BlastOutput_query-len>200</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_hits>
        <Hit>
          <Hit_id>sp|P12345</Hit_id>
          <Hit_def>Protein alpha</Hit_def>
          <Hit_accession>P12345</Hit_accession>
          <Hit_len>300</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>200.5</Hsp_bit-score>
              <Hsp_evalue>1e-55</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>200</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>209</Hsp_hit-to>
              <Hsp_identity>180</Hsp_identity>
              <Hsp_gaps>5</Hsp_gaps>
              <Hsp_align-len>200</Hsp_align-len>
              <Hsp_qseq>MVLK</Hsp_qseq>
              <Hsp_hseq>MVLR</Hsp_hseq>
              <Hsp_midline>MVL </Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_id>sp|Q67890</Hit_id>
          <Hit_def>Protein beta</Hit_def>
          <Hit_accession>Q67890</Hit_accession>
          <Hit_len>250</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>120.3</Hsp_bit-score>
              <Hsp_evalue>1e-30</Hsp_evalue>
              <Hsp_query-from>5</Hsp_query-from>
              <Hsp_query-to>180</Hsp_query-to>
              <Hsp_hit-from>20</Hsp_hit-from>
              <Hsp_hit-to>195</Hsp_hit-to>
              <Hsp_identity>140</Hsp_identity>
              <Hsp_gaps>3</Hsp_gaps>
              <Hsp_align-len>176</Hsp_align-len>
              <Hsp_qseq>ARND</Hsp_qseq>
              <Hsp_hseq>ARSD</Hsp_hseq>
              <Hsp_midline>AR D</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"#;

        let result = parse_blast_xml_str(xml).unwrap();
        assert_eq!(result.iterations[0].hits.len(), 2);
        assert_eq!(result.iterations[0].hits[0].hit_id, "sp|P12345");
        assert_eq!(result.iterations[0].hits[0].hit_def, "Protein alpha");
        assert_eq!(result.iterations[0].hits[1].hit_id, "sp|Q67890");
        assert_eq!(result.iterations[0].hits[1].hit_def, "Protein beta");
        assert_eq!(result.iterations[0].hits[1].hit_accession, "Q67890");
        assert_eq!(result.iterations[0].hits[1].hit_len, 250);
    }

    #[test]
    fn blast_xml_hsp_fields() {
        let result = parse_blast_xml_str(MINIMAL_BLAST_XML).unwrap();
        let hsp = &result.iterations[0].hits[0].hsps[0];

        assert!((hsp.bit_score - 150.0).abs() < f64::EPSILON);
        assert!((hsp.evalue - 1e-40).abs() < 1e-50);
        assert_eq!(hsp.query_from, 1);
        assert_eq!(hsp.query_to, 100);
        assert_eq!(hsp.hit_from, 50);
        assert_eq!(hsp.hit_to, 149);
        assert_eq!(hsp.identity, 95);
        assert_eq!(hsp.gaps, 2);
        assert_eq!(hsp.align_len, 100);
    }

    #[test]
    fn blast_xml_sequences() {
        let result = parse_blast_xml_str(MINIMAL_BLAST_XML).unwrap();
        let hsp = &result.iterations[0].hits[0].hsps[0];

        assert_eq!(hsp.qseq, "ATGCATGC");
        assert_eq!(hsp.hseq, "ATGNATGC");
        assert_eq!(hsp.midline, "||| ||||");
    }

    #[test]
    fn blast_xml_malformed_error() {
        let bad_xml = "<not_blast>some garbage</not_blast>";
        let result = parse_blast_xml_str(bad_xml);
        assert!(result.is_err());

        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("BlastOutput"));
    }
}
