//! Nanopore-specific analysis: signal metadata, methylation calling, and read QC.
//!
//! Provides types and functions for Oxford Nanopore-specific data including
//! raw signal metadata (FAST5/POD5/BLOW5), CpG methylation detection from
//! modified base calls, and Nanopore-specific quality metrics.

use cyanea_core::{CyaneaError, Result};

/// Nanopore signal metadata (extracted from FAST5/POD5/BLOW5 headers).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SignalMetadata {
    /// Read identifier (UUID).
    pub read_id: String,
    /// Run identifier.
    pub run_id: String,
    /// Flow cell ID.
    pub flow_cell_id: String,
    /// Experiment name.
    pub experiment: String,
    /// Sample name.
    pub sample_id: String,
    /// Channel number.
    pub channel: u32,
    /// Start time (seconds from run start).
    pub start_time: f64,
    /// Duration of the read (seconds).
    pub duration: f64,
    /// Sampling rate (Hz).
    pub sampling_rate: f64,
    /// Number of raw signal samples.
    pub num_samples: u64,
    /// Median raw signal level (pA).
    pub median_signal: f64,
    /// Signal-to-noise ratio estimate.
    pub signal_to_noise: Option<f64>,
    /// Basecaller used.
    pub basecaller: Option<String>,
    /// Basecaller version.
    pub basecaller_version: Option<String>,
}

/// Parse a FAST5/POD5/BLOW5 metadata summary line.
///
/// Accepts tab-separated lines with key=value pairs as exported by
/// `pod5 inspect` or `slow5tools stats`.
///
/// # Format
///
/// ```text
/// read_id=uuid run_id=run channel=42 start_time=100.5 duration=5.2 ...
/// ```
pub fn parse_signal_metadata(line: &str) -> Result<SignalMetadata> {
    let mut meta = SignalMetadata {
        read_id: String::new(),
        run_id: String::new(),
        flow_cell_id: String::new(),
        experiment: String::new(),
        sample_id: String::new(),
        channel: 0,
        start_time: 0.0,
        duration: 0.0,
        sampling_rate: 4000.0,
        num_samples: 0,
        median_signal: 0.0,
        signal_to_noise: None,
        basecaller: None,
        basecaller_version: None,
    };

    for part in line.split_whitespace() {
        if let Some((key, val)) = part.split_once('=') {
            match key {
                "read_id" => meta.read_id = val.to_string(),
                "run_id" => meta.run_id = val.to_string(),
                "flow_cell_id" => meta.flow_cell_id = val.to_string(),
                "experiment" => meta.experiment = val.to_string(),
                "sample_id" => meta.sample_id = val.to_string(),
                "channel" => meta.channel = val.parse().unwrap_or(0),
                "start_time" => meta.start_time = val.parse().unwrap_or(0.0),
                "duration" => meta.duration = val.parse().unwrap_or(0.0),
                "sampling_rate" => meta.sampling_rate = val.parse().unwrap_or(4000.0),
                "num_samples" => meta.num_samples = val.parse().unwrap_or(0),
                "median_signal" => meta.median_signal = val.parse().unwrap_or(0.0),
                "signal_to_noise" => meta.signal_to_noise = val.parse().ok(),
                "basecaller" => meta.basecaller = Some(val.to_string()),
                "basecaller_version" => meta.basecaller_version = Some(val.to_string()),
                _ => {}
            }
        }
    }

    if meta.read_id.is_empty() {
        return Err(CyaneaError::Parse("missing read_id in signal metadata".into()));
    }

    Ok(meta)
}

/// Methylation call from modified base detection.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MethylationCall {
    /// Chromosome / contig.
    pub chrom: String,
    /// Position (0-based) on reference.
    pub position: u64,
    /// Strand ('+' or '-').
    pub strand: char,
    /// Modification type.
    pub mod_type: ModificationType,
    /// Probability of modification (0.0 - 1.0).
    pub probability: f64,
    /// Read name that produced this call.
    pub read_id: String,
}

/// Type of base modification detectable by Nanopore.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ModificationType {
    /// 5-methylcytosine (5mC) — most common.
    FiveMC,
    /// 5-hydroxymethylcytosine (5hmC).
    FiveHMC,
    /// 6-methyladenine (6mA).
    SixMA,
    /// 4-methylcytosine (4mC).
    FourMC,
}

impl ModificationType {
    /// Standard abbreviation.
    pub fn as_str(&self) -> &'static str {
        match self {
            ModificationType::FiveMC => "5mC",
            ModificationType::FiveHMC => "5hmC",
            ModificationType::SixMA => "6mA",
            ModificationType::FourMC => "4mC",
        }
    }
}

/// Parse methylation calls from a modBAM-style format.
///
/// Each line: `read_id\tchrom\tposition\tstrand\tmod_type\tprobability`
pub fn parse_methylation_calls(text: &str) -> Result<Vec<MethylationCall>> {
    let mut calls = Vec::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            return Err(CyaneaError::Parse(format!("expected 6 tab-separated fields, got {}", fields.len())));
        }

        let mod_type = match fields[4] {
            "5mC" | "m" => ModificationType::FiveMC,
            "5hmC" | "h" => ModificationType::FiveHMC,
            "6mA" | "a" => ModificationType::SixMA,
            "4mC" => ModificationType::FourMC,
            other => return Err(CyaneaError::Parse(format!("unknown modification: {}", other))),
        };

        calls.push(MethylationCall {
            read_id: fields[0].to_string(),
            chrom: fields[1].to_string(),
            position: fields[2].parse().map_err(|_| CyaneaError::Parse("invalid position".into()))?,
            strand: fields[3].chars().next().unwrap_or('+'),
            mod_type,
            probability: fields[5].parse().map_err(|_| CyaneaError::Parse("invalid probability".into()))?,
        });
    }

    Ok(calls)
}

/// Per-site methylation summary across reads.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MethylationSite {
    /// Chromosome.
    pub chrom: String,
    /// Position.
    pub position: u64,
    /// Strand.
    pub strand: char,
    /// Modification type.
    pub mod_type: ModificationType,
    /// Number of reads covering this site.
    pub coverage: usize,
    /// Number of reads showing modification.
    pub modified_count: usize,
    /// Mean modification probability.
    pub mean_probability: f64,
    /// Methylation frequency (modified_count / coverage).
    pub frequency: f64,
}

/// Aggregate per-read methylation calls into per-site summaries.
///
/// Groups calls by (chrom, position, strand, mod_type) and computes
/// coverage, modified count (probability >= threshold), and mean probability.
pub fn aggregate_methylation(
    calls: &[MethylationCall],
    threshold: f64,
) -> Vec<MethylationSite> {
    let mut sites: std::collections::HashMap<(String, u64, char, &str), Vec<f64>> =
        std::collections::HashMap::new();

    for call in calls {
        let key = (call.chrom.clone(), call.position, call.strand, call.mod_type.as_str());
        sites.entry(key).or_default().push(call.probability);
    }

    let mut results: Vec<MethylationSite> = sites.into_iter().map(|((chrom, position, strand, mod_str), probs)| {
        let coverage = probs.len();
        let modified_count = probs.iter().filter(|&&p| p >= threshold).count();
        let mean_probability = probs.iter().sum::<f64>() / coverage as f64;
        let frequency = modified_count as f64 / coverage as f64;

        let mod_type = match mod_str {
            "5mC" => ModificationType::FiveMC,
            "5hmC" => ModificationType::FiveHMC,
            "6mA" => ModificationType::SixMA,
            "4mC" => ModificationType::FourMC,
            _ => ModificationType::FiveMC,
        };

        MethylationSite {
            chrom, position, strand, mod_type,
            coverage, modified_count, mean_probability, frequency,
        }
    }).collect();

    results.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.position.cmp(&b.position)));
    results
}

/// Nanopore-specific QC metrics for a sequencing run.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct NanoporeQC {
    /// Total reads.
    pub total_reads: usize,
    /// Total bases.
    pub total_bases: u64,
    /// N50 read length.
    pub n50: usize,
    /// Reads passing Q-score filter.
    pub pass_reads: usize,
    /// Reads failing Q-score filter.
    pub fail_reads: usize,
    /// Mean quality score.
    pub mean_quality: f64,
    /// Estimated pore occupancy (fraction of channels producing reads).
    pub pore_occupancy: f64,
    /// Reads per channel (mean).
    pub reads_per_channel: f64,
    /// Translocation speed (bases/second, estimated).
    pub translocation_speed: f64,
}

/// Compute Nanopore-specific QC from signal metadata and read stats.
pub fn nanopore_qc(
    metadata: &[SignalMetadata],
    read_lengths: &[usize],
    read_qualities: &[f64],
    total_channels: u32,
    quality_threshold: f64,
) -> Result<NanoporeQC> {
    if metadata.is_empty() || read_lengths.is_empty() {
        return Err(CyaneaError::InvalidInput("no data for QC".into()));
    }

    let total_reads = read_lengths.len();
    let total_bases: u64 = read_lengths.iter().map(|&l| l as u64).sum();

    // N50
    let mut sorted_lengths = read_lengths.to_vec();
    sorted_lengths.sort_unstable();
    let half = total_bases / 2;
    let mut cumulative = 0u64;
    let mut n50 = 0usize;
    for &l in sorted_lengths.iter().rev() {
        cumulative += l as u64;
        if cumulative >= half {
            n50 = l;
            break;
        }
    }

    let pass_reads = read_qualities.iter().filter(|&&q| q >= quality_threshold).count();
    let fail_reads = total_reads - pass_reads;
    let mean_quality = read_qualities.iter().sum::<f64>() / total_reads as f64;

    // Pore occupancy: fraction of channels that produced at least one read
    let active_channels: std::collections::HashSet<u32> = metadata.iter()
        .map(|m| m.channel)
        .collect();
    let pore_occupancy = active_channels.len() as f64 / total_channels.max(1) as f64;
    let reads_per_channel = total_reads as f64 / active_channels.len().max(1) as f64;

    // Translocation speed: bases / duration
    let total_duration: f64 = metadata.iter().map(|m| m.duration).sum();
    let translocation_speed = if total_duration > 0.0 {
        total_bases as f64 / total_duration
    } else {
        0.0
    };

    Ok(NanoporeQC {
        total_reads,
        total_bases,
        n50,
        pass_reads,
        fail_reads,
        mean_quality,
        pore_occupancy,
        reads_per_channel,
        translocation_speed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_signal_metadata() {
        let line = "read_id=abc123 run_id=run1 flow_cell_id=FC001 channel=42 start_time=100.5 duration=5.2 sampling_rate=4000 num_samples=20800 median_signal=85.3";
        let meta = parse_signal_metadata(line).unwrap();
        assert_eq!(meta.read_id, "abc123");
        assert_eq!(meta.run_id, "run1");
        assert_eq!(meta.flow_cell_id, "FC001");
        assert_eq!(meta.channel, 42);
        assert!((meta.start_time - 100.5).abs() < 1e-10);
        assert!((meta.duration - 5.2).abs() < 1e-10);
        assert_eq!(meta.num_samples, 20800);
    }

    #[test]
    fn test_parse_signal_metadata_missing_read_id() {
        let line = "run_id=run1 channel=42";
        assert!(parse_signal_metadata(line).is_err());
    }

    #[test]
    fn test_parse_methylation_calls() {
        let text = "\
read1\tchr1\t1000\t+\t5mC\t0.95
read2\tchr1\t1000\t+\t5mC\t0.85
read3\tchr1\t2000\t-\t5mC\t0.10
";
        let calls = parse_methylation_calls(text).unwrap();
        assert_eq!(calls.len(), 3);
        assert_eq!(calls[0].chrom, "chr1");
        assert_eq!(calls[0].position, 1000);
        assert_eq!(calls[0].mod_type, ModificationType::FiveMC);
        assert!((calls[0].probability - 0.95).abs() < 1e-10);
    }

    #[test]
    fn test_parse_methylation_shorthand() {
        let text = "read1\tchr1\t500\t+\tm\t0.9\n";
        let calls = parse_methylation_calls(text).unwrap();
        assert_eq!(calls[0].mod_type, ModificationType::FiveMC);
    }

    #[test]
    fn test_parse_methylation_6ma() {
        let text = "read1\tchr1\t500\t+\t6mA\t0.8\n";
        let calls = parse_methylation_calls(text).unwrap();
        assert_eq!(calls[0].mod_type, ModificationType::SixMA);
    }

    #[test]
    fn test_parse_methylation_invalid() {
        let text = "read1\tchr1\t500";
        assert!(parse_methylation_calls(text).is_err());
    }

    #[test]
    fn test_aggregate_methylation() {
        let calls = vec![
            MethylationCall { read_id: "r1".into(), chrom: "chr1".into(), position: 1000, strand: '+', mod_type: ModificationType::FiveMC, probability: 0.95 },
            MethylationCall { read_id: "r2".into(), chrom: "chr1".into(), position: 1000, strand: '+', mod_type: ModificationType::FiveMC, probability: 0.85 },
            MethylationCall { read_id: "r3".into(), chrom: "chr1".into(), position: 1000, strand: '+', mod_type: ModificationType::FiveMC, probability: 0.10 },
        ];

        let sites = aggregate_methylation(&calls, 0.5);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].coverage, 3);
        assert_eq!(sites[0].modified_count, 2); // 0.95 and 0.85 pass 0.5 threshold
        assert!((sites[0].frequency - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_aggregate_multiple_sites() {
        let calls = vec![
            MethylationCall { read_id: "r1".into(), chrom: "chr1".into(), position: 100, strand: '+', mod_type: ModificationType::FiveMC, probability: 0.9 },
            MethylationCall { read_id: "r2".into(), chrom: "chr1".into(), position: 200, strand: '+', mod_type: ModificationType::FiveMC, probability: 0.8 },
            MethylationCall { read_id: "r3".into(), chrom: "chr2".into(), position: 100, strand: '-', mod_type: ModificationType::SixMA, probability: 0.7 },
        ];

        let sites = aggregate_methylation(&calls, 0.5);
        assert_eq!(sites.len(), 3); // three different sites
    }

    #[test]
    fn test_modification_type_str() {
        assert_eq!(ModificationType::FiveMC.as_str(), "5mC");
        assert_eq!(ModificationType::FiveHMC.as_str(), "5hmC");
        assert_eq!(ModificationType::SixMA.as_str(), "6mA");
        assert_eq!(ModificationType::FourMC.as_str(), "4mC");
    }

    #[test]
    fn test_nanopore_qc() {
        let metadata: Vec<SignalMetadata> = (0..100).map(|i| SignalMetadata {
            read_id: format!("read_{}", i),
            run_id: "run1".into(),
            flow_cell_id: "FC001".into(),
            experiment: "exp1".into(),
            sample_id: "sample1".into(),
            channel: (i % 50) as u32,
            start_time: i as f64 * 10.0,
            duration: 2.0,
            sampling_rate: 4000.0,
            num_samples: 8000,
            median_signal: 80.0,
            signal_to_noise: Some(10.0),
            basecaller: Some("dorado".into()),
            basecaller_version: Some("0.5.0".into()),
        }).collect();

        let lengths: Vec<usize> = (0..100).map(|i| 5000 + i * 100).collect();
        let qualities: Vec<f64> = (0..100).map(|i| 10.0 + (i as f64) * 0.2).collect();

        let qc = nanopore_qc(&metadata, &lengths, &qualities, 512, 10.0).unwrap();
        assert_eq!(qc.total_reads, 100);
        assert!(qc.total_bases > 0);
        assert!(qc.n50 > 0);
        assert!(qc.pore_occupancy > 0.0 && qc.pore_occupancy <= 1.0);
        assert!(qc.translocation_speed > 0.0);
        assert!(qc.pass_reads + qc.fail_reads == 100);
    }

    #[test]
    fn test_nanopore_qc_empty() {
        assert!(nanopore_qc(&[], &[], &[], 512, 10.0).is_err());
    }
}
