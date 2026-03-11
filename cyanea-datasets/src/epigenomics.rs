//! Epigenomics demo datasets: ChIP-seq peaks, methylation sites.

/// Demo ChIP-seq narrow peaks (H3K27ac, promoter marks).
///
/// 10 peaks on chr1 with typical ChIP-seq signal characteristics.
pub fn demo_chipseq_peaks() -> Vec<DemoPeak> {
    vec![
        DemoPeak { chrom: "chr1", start: 9873, end: 10472, name: "peak1", score: 245.0, signal: 12.5, pvalue: 4.2e-8, qvalue: 1.1e-5, summit: 10150 },
        DemoPeak { chrom: "chr1", start: 15523, end: 15903, name: "peak2", score: 180.0, signal: 9.3, pvalue: 2.1e-6, qvalue: 3.4e-4, summit: 15710 },
        DemoPeak { chrom: "chr1", start: 29553, end: 30039, name: "peak3", score: 312.0, signal: 15.8, pvalue: 1.3e-10, qvalue: 8.7e-8, summit: 29800 },
        DemoPeak { chrom: "chr1", start: 69090, end: 70008, name: "peak4", score: 520.0, signal: 25.1, pvalue: 6.5e-15, qvalue: 2.2e-12, summit: 69540 },
        DemoPeak { chrom: "chr1", start: 103801, end: 104285, name: "peak5", score: 198.0, signal: 10.2, pvalue: 8.9e-7, qvalue: 1.5e-4, summit: 104050 },
        DemoPeak { chrom: "chr1", start: 137511, end: 137960, name: "peak6", score: 156.0, signal: 8.1, pvalue: 5.4e-5, qvalue: 6.2e-3, summit: 137730 },
        DemoPeak { chrom: "chr1", start: 200001, end: 200510, name: "peak7", score: 410.0, signal: 20.3, pvalue: 3.2e-12, qvalue: 1.8e-9, summit: 200250 },
        DemoPeak { chrom: "chr1", start: 235230, end: 235890, name: "peak8", score: 275.0, signal: 13.7, pvalue: 7.8e-9, qvalue: 3.1e-6, summit: 235560 },
        DemoPeak { chrom: "chr1", start: 567890, end: 568350, name: "peak9", score: 340.0, signal: 17.2, pvalue: 4.5e-11, qvalue: 5.3e-8, summit: 568120 },
        DemoPeak { chrom: "chr1", start: 778210, end: 778690, name: "peak10", score: 290.0, signal: 14.5, pvalue: 2.3e-9, qvalue: 1.7e-6, summit: 778450 },
    ]
}

/// A ChIP-seq narrow peak.
#[derive(Debug, Clone)]
pub struct DemoPeak {
    pub chrom: &'static str,
    pub start: u64,
    pub end: u64,
    pub name: &'static str,
    pub score: f64,
    pub signal: f64,
    pub pvalue: f64,
    pub qvalue: f64,
    pub summit: u64,
}

impl DemoPeak {
    /// Format as narrowPeak (BED6+4) line.
    pub fn to_narrowpeak(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t.\t{:.1}\t{:.2e}\t{:.2e}\t{}",
            self.chrom,
            self.start,
            self.end,
            self.name,
            self.score as u32,
            self.signal,
            self.pvalue,
            self.qvalue,
            self.summit - self.start
        )
    }

    /// Format as BED line.
    pub fn to_bed(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t.",
            self.chrom, self.start, self.end, self.name, self.score as u32
        )
    }
}

/// Generate narrowPeak format output from demo peaks.
pub fn demo_narrowpeak() -> String {
    demo_chipseq_peaks()
        .iter()
        .map(|p| p.to_narrowpeak())
        .collect::<Vec<_>>()
        .join("\n")
}

/// Demo CpG methylation data (bisulfite sequencing).
///
/// 12 CpG sites in a promoter region with methylation levels.
pub fn demo_cpg_methylation() -> Vec<DemoCpgSite> {
    vec![
        DemoCpgSite { chrom: "chr1", pos: 10468, methylated: 2, total: 15, beta: 0.133 },
        DemoCpgSite { chrom: "chr1", pos: 10470, methylated: 1, total: 12, beta: 0.083 },
        DemoCpgSite { chrom: "chr1", pos: 10483, methylated: 0, total: 18, beta: 0.0 },
        DemoCpgSite { chrom: "chr1", pos: 10488, methylated: 8, total: 20, beta: 0.4 },
        DemoCpgSite { chrom: "chr1", pos: 10492, methylated: 14, total: 22, beta: 0.636 },
        DemoCpgSite { chrom: "chr1", pos: 10496, methylated: 18, total: 25, beta: 0.72 },
        DemoCpgSite { chrom: "chr1", pos: 10524, methylated: 20, total: 21, beta: 0.952 },
        DemoCpgSite { chrom: "chr1", pos: 10541, methylated: 19, total: 20, beta: 0.95 },
        DemoCpgSite { chrom: "chr1", pos: 10562, methylated: 17, total: 19, beta: 0.895 },
        DemoCpgSite { chrom: "chr1", pos: 10571, methylated: 3, total: 16, beta: 0.188 },
        DemoCpgSite { chrom: "chr1", pos: 10576, methylated: 1, total: 14, beta: 0.071 },
        DemoCpgSite { chrom: "chr1", pos: 10589, methylated: 0, total: 17, beta: 0.0 },
    ]
}

/// A CpG methylation site.
#[derive(Debug, Clone)]
pub struct DemoCpgSite {
    pub chrom: &'static str,
    pub pos: u64,
    pub methylated: u32,
    pub total: u32,
    pub beta: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chipseq_peaks() {
        let peaks = demo_chipseq_peaks();
        assert_eq!(peaks.len(), 10);
        assert!(peaks.iter().all(|p| p.start < p.end));
        assert!(peaks.iter().all(|p| p.summit >= p.start && p.summit <= p.end));
    }

    #[test]
    fn test_narrowpeak_format() {
        let peak = &demo_chipseq_peaks()[0];
        let line = peak.to_narrowpeak();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 10);
        assert_eq!(fields[0], "chr1");
    }

    #[test]
    fn test_bed_format() {
        let peak = &demo_chipseq_peaks()[0];
        let line = peak.to_bed();
        assert!(line.starts_with("chr1\t"));
    }

    #[test]
    fn test_narrowpeak_output() {
        let np = demo_narrowpeak();
        let lines: Vec<&str> = np.lines().collect();
        assert_eq!(lines.len(), 10);
    }

    #[test]
    fn test_cpg_methylation() {
        let sites = demo_cpg_methylation();
        assert_eq!(sites.len(), 12);
        assert!(sites.iter().all(|s| s.beta >= 0.0 && s.beta <= 1.0));
        assert!(sites.iter().all(|s| s.methylated <= s.total));
    }
}
