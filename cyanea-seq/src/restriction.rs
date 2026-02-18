//! Restriction enzyme recognition, cut-site finding, and in-silico digestion.
//!
//! Provides a database of common restriction enzymes and functions to
//! find cut sites, digest sequences, and compute fragment sizes.

/// A restriction enzyme with its recognition site and cut offsets.
#[derive(Debug, Clone)]
pub struct RestrictionEnzyme {
    /// Enzyme name (e.g., "EcoRI").
    pub name: String,
    /// Recognition site in IUPAC DNA (uppercase).
    pub recognition_site: Vec<u8>,
    /// Cut offset on the forward (5'→3') strand, from the start of the recognition site.
    pub cut_forward: isize,
    /// Cut offset on the reverse (3'→5') strand, from the start of the recognition site.
    pub cut_reverse: isize,
    /// Whether the recognition site is palindromic.
    pub is_palindromic: bool,
}

/// A located cut site on a sequence.
#[derive(Debug, Clone)]
pub struct CutSite {
    /// Position in the original sequence where the forward strand is cut.
    pub position: usize,
    /// Name of the enzyme that cuts here.
    pub enzyme: String,
    /// Type of overhang produced.
    pub overhang: Overhang,
}

/// Overhang type produced by a restriction enzyme cut.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Overhang {
    /// 5' overhang (sticky end).
    FivePrime(Vec<u8>),
    /// 3' overhang (sticky end).
    ThreePrime(Vec<u8>),
    /// Blunt end (no overhang).
    Blunt,
}

/// A fragment produced by restriction digestion.
#[derive(Debug, Clone)]
pub struct Fragment {
    /// Start position in the original sequence (inclusive).
    pub start: usize,
    /// End position in the original sequence (exclusive).
    pub end: usize,
    /// Fragment length.
    pub length: usize,
}

/// Return a curated set of ~20 common restriction enzymes.
pub fn common_enzymes() -> Vec<RestrictionEnzyme> {
    vec![
        enzyme("EcoRI", b"GAATTC", 1, 5, true),
        enzyme("BamHI", b"GGATCC", 1, 5, true),
        enzyme("HindIII", b"AAGCTT", 1, 5, true),
        enzyme("NotI", b"GCGGCCGC", 2, 6, true),
        enzyme("XhoI", b"CTCGAG", 1, 5, true),
        enzyme("SalI", b"GTCGAC", 1, 5, true),
        enzyme("BglII", b"AGATCT", 1, 5, true),
        enzyme("NcoI", b"CCATGG", 1, 5, true),
        enzyme("NdeI", b"CATATG", 2, 4, true),
        enzyme("XbaI", b"TCTAGA", 1, 5, true),
        enzyme("SpeI", b"ACTAGT", 1, 5, true),
        enzyme("KpnI", b"GGTACC", 5, 1, true),
        enzyme("SacI", b"GAGCTC", 5, 1, true),
        enzyme("PstI", b"CTGCAG", 5, 1, true),
        enzyme("SphI", b"GCATGC", 5, 1, true),
        enzyme("ApaI", b"GGGCCC", 5, 1, true),
        enzyme("EcoRV", b"GATATC", 3, 3, true),
        enzyme("SmaI", b"CCCGGG", 3, 3, true),
        enzyme("HpaI", b"GTTAAC", 3, 3, true),
        enzyme("ScaI", b"AGTACT", 3, 3, true),
    ]
}

fn enzyme(name: &str, site: &[u8], fwd: isize, rev: isize, pal: bool) -> RestrictionEnzyme {
    RestrictionEnzyme {
        name: name.to_string(),
        recognition_site: site.to_vec(),
        cut_forward: fwd,
        cut_reverse: rev,
        is_palindromic: pal,
    }
}

/// Check if a single base matches an IUPAC degenerate code.
fn iupac_matches(code: u8, base: u8) -> bool {
    let base_upper = base.to_ascii_uppercase();
    match code.to_ascii_uppercase() {
        b'A' => base_upper == b'A',
        b'C' => base_upper == b'C',
        b'G' => base_upper == b'G',
        b'T' => base_upper == b'T',
        b'R' => matches!(base_upper, b'A' | b'G'),
        b'Y' => matches!(base_upper, b'C' | b'T'),
        b'M' => matches!(base_upper, b'A' | b'C'),
        b'K' => matches!(base_upper, b'G' | b'T'),
        b'S' => matches!(base_upper, b'G' | b'C'),
        b'W' => matches!(base_upper, b'A' | b'T'),
        b'H' => matches!(base_upper, b'A' | b'C' | b'T'),
        b'B' => matches!(base_upper, b'C' | b'G' | b'T'),
        b'V' => matches!(base_upper, b'A' | b'C' | b'G'),
        b'D' => matches!(base_upper, b'A' | b'G' | b'T'),
        b'N' => matches!(base_upper, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

/// Find all cut sites for a single enzyme on the given sequence.
///
/// Scans the sequence for IUPAC-aware matches of the recognition site.
pub fn find_cut_sites(seq: &[u8], enzyme: &RestrictionEnzyme) -> Vec<CutSite> {
    let site_len = enzyme.recognition_site.len();
    if seq.len() < site_len {
        return Vec::new();
    }

    let mut sites = Vec::new();
    for i in 0..=seq.len() - site_len {
        let matches = enzyme
            .recognition_site
            .iter()
            .zip(&seq[i..i + site_len])
            .all(|(&code, &base)| iupac_matches(code, base));

        if matches {
            let cut_pos = i as isize + enzyme.cut_forward;
            let cut_pos = cut_pos.max(0) as usize;

            let overhang = if enzyme.cut_forward == enzyme.cut_reverse {
                Overhang::Blunt
            } else if enzyme.cut_forward < enzyme.cut_reverse {
                // 5' overhang
                let start = enzyme.cut_forward.max(0) as usize;
                let end = (enzyme.cut_reverse as usize).min(site_len);
                let oh = enzyme.recognition_site[start..end].to_vec();
                Overhang::FivePrime(oh)
            } else {
                // 3' overhang
                let start = enzyme.cut_reverse.max(0) as usize;
                let end = (enzyme.cut_forward as usize).min(site_len);
                let oh = enzyme.recognition_site[start..end].to_vec();
                Overhang::ThreePrime(oh)
            };

            sites.push(CutSite {
                position: cut_pos,
                enzyme: enzyme.name.clone(),
                overhang,
            });
        }
    }
    sites
}

/// Digest a sequence with one or more enzymes, returning the resulting fragments.
///
/// Fragments are returned sorted by position. The first fragment starts at
/// position 0 and the last ends at the sequence length.
pub fn digest(seq: &[u8], enzymes: &[&RestrictionEnzyme]) -> Vec<Fragment> {
    let mut cut_positions: Vec<usize> = Vec::new();
    for enz in enzymes {
        for site in find_cut_sites(seq, enz) {
            cut_positions.push(site.position);
        }
    }
    cut_positions.sort_unstable();
    cut_positions.dedup();

    // Build fragments between consecutive cut positions.
    let mut fragments = Vec::new();
    let mut prev = 0usize;
    for &pos in &cut_positions {
        if pos > prev && pos <= seq.len() {
            fragments.push(Fragment {
                start: prev,
                end: pos,
                length: pos - prev,
            });
            prev = pos;
        }
    }
    // Final fragment from last cut to end.
    if prev < seq.len() {
        fragments.push(Fragment {
            start: prev,
            end: seq.len(),
            length: seq.len() - prev,
        });
    }
    fragments
}

/// Compute fragment sizes from digestion with one or more enzymes.
pub fn fragment_sizes(seq: &[u8], enzymes: &[&RestrictionEnzyme]) -> Vec<usize> {
    digest(seq, enzymes)
        .into_iter()
        .map(|f| f.length)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ecori_cuts_correctly() {
        // EcoRI recognizes GAATTC, cuts between G and AATTC (offset 1).
        let seq = b"AAAGAATTCAAA";
        let ecori = &common_enzymes().into_iter().find(|e| e.name == "EcoRI").unwrap();
        let sites = find_cut_sites(seq, ecori);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].position, 4); // 3 + 1
        assert!(matches!(sites[0].overhang, Overhang::FivePrime(_)));
    }

    #[test]
    fn blunt_end_enzyme() {
        // EcoRV: GATATC, cuts at 3/3 → blunt.
        let seq = b"AAAGATATCAAA";
        let ecorv = &common_enzymes().into_iter().find(|e| e.name == "EcoRV").unwrap();
        let sites = find_cut_sites(seq, ecorv);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].overhang, Overhang::Blunt);
    }

    #[test]
    fn no_sites_single_fragment() {
        let seq = b"AAAAAAAAAA";
        let ecori = &common_enzymes().into_iter().find(|e| e.name == "EcoRI").unwrap();
        let fragments = digest(seq, &[ecori]);
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].length, 10);
    }

    #[test]
    fn iupac_degenerate_match() {
        // Create an enzyme with degenerate site: GAANTC (N matches any).
        let enz = RestrictionEnzyme {
            name: "TestEnz".into(),
            recognition_site: b"GAANTC".to_vec(),
            cut_forward: 1,
            cut_reverse: 5,
            is_palindromic: true,
        };
        // GAATTC should match (N=T), GAACTC should also match (N=C).
        assert_eq!(find_cut_sites(b"GAATTC", &enz).len(), 1);
        assert_eq!(find_cut_sites(b"GAACTC", &enz).len(), 1);
        assert_eq!(find_cut_sites(b"GAAATC", &enz).len(), 1);
        assert_eq!(find_cut_sites(b"GAAGTC", &enz).len(), 1);
    }

    #[test]
    fn common_enzymes_nonempty() {
        let enzymes = common_enzymes();
        assert_eq!(enzymes.len(), 20);
        for e in &enzymes {
            assert!(!e.name.is_empty());
            assert!(!e.recognition_site.is_empty());
        }
    }
}
