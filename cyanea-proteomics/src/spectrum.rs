//! Mass spectrum data types and operations.
//!
//! Core types for representing mass spectrometry data: peaks, spectra,
//! precursor information, and spectrum processing operations.

use crate::error::{ProteomicsError, Result};

/// A single peak in a mass spectrum.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Peak {
    /// Mass-to-charge ratio (m/z).
    pub mz: f64,
    /// Signal intensity.
    pub intensity: f64,
}

/// MS level (MS1 = survey scan, MS2 = fragmentation, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MsLevel {
    Ms1,
    Ms2,
    Ms3,
    MsN(u8),
}

impl MsLevel {
    pub fn from_u8(level: u8) -> Self {
        match level {
            1 => MsLevel::Ms1,
            2 => MsLevel::Ms2,
            3 => MsLevel::Ms3,
            n => MsLevel::MsN(n),
        }
    }

    pub fn as_u8(&self) -> u8 {
        match self {
            MsLevel::Ms1 => 1,
            MsLevel::Ms2 => 2,
            MsLevel::Ms3 => 3,
            MsLevel::MsN(n) => *n,
        }
    }
}

/// Fragmentation method used to generate MS2+ spectra.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum FragmentationMethod {
    /// Collision-induced dissociation.
    CID,
    /// Higher-energy collisional dissociation.
    HCD,
    /// Electron-transfer dissociation.
    ETD,
    /// Electron-transfer/higher-energy collision dissociation.
    EThcD,
    /// Unknown or unspecified method.
    Unknown,
}

/// Precursor ion information for MS2+ scans.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Precursor {
    /// Precursor m/z.
    pub mz: f64,
    /// Charge state (if known).
    pub charge: Option<i32>,
    /// Precursor intensity (if known).
    pub intensity: Option<f64>,
    /// Isolation window width (if known).
    pub isolation_width: Option<f64>,
    /// Fragmentation method.
    pub fragmentation: FragmentationMethod,
}

/// A mass spectrum with metadata.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MassSpectrum {
    /// Spectrum identifier (scan number or native ID).
    pub id: String,
    /// MS level.
    pub ms_level: MsLevel,
    /// Retention time in seconds.
    pub retention_time: f64,
    /// Peaks sorted by m/z.
    pub peaks: Vec<Peak>,
    /// Precursor information (for MS2+ scans).
    pub precursor: Option<Precursor>,
    /// Total ion current.
    pub tic: f64,
    /// Base peak intensity.
    pub base_peak_intensity: f64,
    /// Base peak m/z.
    pub base_peak_mz: f64,
}

impl MassSpectrum {
    /// Create a new mass spectrum from peaks.
    ///
    /// Peaks are sorted by m/z. TIC, base peak intensity, and base peak m/z
    /// are computed automatically.
    pub fn new(id: impl Into<String>, ms_level: MsLevel, retention_time: f64, mut peaks: Vec<Peak>) -> Result<Self> {
        if peaks.is_empty() {
            return Err(ProteomicsError::Spectrum("spectrum must have at least one peak".into()));
        }
        for p in &peaks {
            if p.mz < 0.0 || p.intensity < 0.0 {
                return Err(ProteomicsError::Spectrum("negative m/z or intensity".into()));
            }
        }

        peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));

        let tic: f64 = peaks.iter().map(|p| p.intensity).sum();
        let base = peaks.iter().max_by(|a, b| a.intensity.partial_cmp(&b.intensity).unwrap()).unwrap();

        Ok(Self {
            id: id.into(),
            ms_level,
            retention_time,
            tic,
            base_peak_intensity: base.intensity,
            base_peak_mz: base.mz,
            peaks,
            precursor: None,
        })
    }

    /// Set precursor information (for MS2+ spectra).
    pub fn with_precursor(mut self, precursor: Precursor) -> Self {
        self.precursor = Some(precursor);
        self
    }

    /// Number of peaks.
    pub fn num_peaks(&self) -> usize {
        self.peaks.len()
    }

    /// Find the closest peak to a target m/z within a tolerance (in Da).
    pub fn find_peak(&self, target_mz: f64, tolerance: f64) -> Option<&Peak> {
        let mut best: Option<&Peak> = None;
        let mut best_dist = tolerance;

        for p in &self.peaks {
            let dist = (p.mz - target_mz).abs();
            if dist < best_dist {
                best_dist = dist;
                best = Some(p);
            }
            // Since peaks are sorted, we can stop once we pass the window
            if p.mz > target_mz + tolerance {
                break;
            }
        }
        best
    }

    /// Filter peaks below a relative intensity threshold (0.0-1.0 of base peak).
    pub fn filter_by_relative_intensity(&mut self, min_relative: f64) {
        let threshold = self.base_peak_intensity * min_relative;
        self.peaks.retain(|p| p.intensity >= threshold);
        self.recalculate_stats();
    }

    /// Keep only the top N most intense peaks.
    pub fn top_n(&mut self, n: usize) {
        if self.peaks.len() <= n {
            return;
        }
        let mut by_intensity: Vec<usize> = (0..self.peaks.len()).collect();
        by_intensity.sort_by(|&a, &b| {
            self.peaks[b].intensity.partial_cmp(&self.peaks[a].intensity).unwrap()
        });
        by_intensity.truncate(n);
        by_intensity.sort();
        self.peaks = by_intensity.into_iter().map(|i| self.peaks[i].clone()).collect();
        self.recalculate_stats();
    }

    /// Bin peaks within a tolerance (Da), keeping the most intense per bin.
    pub fn deisotope(&mut self, tolerance: f64) {
        if self.peaks.len() < 2 {
            return;
        }
        let mut deisotoped = Vec::new();
        let mut i = 0;
        while i < self.peaks.len() {
            let mut best = self.peaks[i].clone();
            let mut j = i + 1;
            while j < self.peaks.len() && (self.peaks[j].mz - best.mz) < tolerance {
                if self.peaks[j].intensity > best.intensity {
                    best = self.peaks[j].clone();
                }
                j += 1;
            }
            deisotoped.push(best);
            i = j;
        }
        self.peaks = deisotoped;
        self.recalculate_stats();
    }

    /// Normalize intensities so the base peak has intensity 100.
    pub fn normalize(&mut self) {
        if self.base_peak_intensity == 0.0 {
            return;
        }
        let factor = 100.0 / self.base_peak_intensity;
        for p in &mut self.peaks {
            p.intensity *= factor;
        }
        self.recalculate_stats();
    }

    /// Sqrt-transform intensities (variance-stabilizing).
    pub fn sqrt_transform(&mut self) {
        for p in &mut self.peaks {
            p.intensity = p.intensity.sqrt();
        }
        self.recalculate_stats();
    }

    fn recalculate_stats(&mut self) {
        self.tic = self.peaks.iter().map(|p| p.intensity).sum();
        if let Some(base) = self.peaks.iter().max_by(|a, b| a.intensity.partial_cmp(&b.intensity).unwrap()) {
            self.base_peak_intensity = base.intensity;
            self.base_peak_mz = base.mz;
        }
    }
}

/// Summary statistics for a collection of spectra (e.g., an entire run).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RunStats {
    /// Total number of spectra.
    pub total_spectra: usize,
    /// Number of MS1 spectra.
    pub ms1_count: usize,
    /// Number of MS2 spectra.
    pub ms2_count: usize,
    /// Number of MS3+ spectra.
    pub msn_count: usize,
    /// Minimum retention time (seconds).
    pub min_rt: f64,
    /// Maximum retention time (seconds).
    pub max_rt: f64,
    /// Mean peaks per MS2.
    pub mean_peaks_per_ms2: f64,
    /// Median peaks per MS2.
    pub median_peaks_per_ms2: f64,
    /// Mean TIC across all spectra.
    pub mean_tic: f64,
}

/// Compute summary statistics for a collection of spectra.
pub fn run_stats(spectra: &[MassSpectrum]) -> Result<RunStats> {
    if spectra.is_empty() {
        return Err(ProteomicsError::Spectrum("no spectra to summarize".into()));
    }

    let mut ms1_count = 0usize;
    let mut ms2_count = 0usize;
    let mut msn_count = 0usize;
    let mut min_rt = f64::MAX;
    let mut max_rt = f64::MIN;
    let mut total_tic = 0.0f64;
    let mut ms2_peaks: Vec<usize> = Vec::new();

    for s in spectra {
        match s.ms_level {
            MsLevel::Ms1 => ms1_count += 1,
            MsLevel::Ms2 => {
                ms2_count += 1;
                ms2_peaks.push(s.num_peaks());
            }
            _ => {
                msn_count += 1;
            }
        }
        if s.retention_time < min_rt {
            min_rt = s.retention_time;
        }
        if s.retention_time > max_rt {
            max_rt = s.retention_time;
        }
        total_tic += s.tic;
    }

    let mean_peaks = if ms2_peaks.is_empty() {
        0.0
    } else {
        ms2_peaks.iter().sum::<usize>() as f64 / ms2_peaks.len() as f64
    };

    let median_peaks = if ms2_peaks.is_empty() {
        0.0
    } else {
        ms2_peaks.sort();
        let mid = ms2_peaks.len() / 2;
        if ms2_peaks.len() % 2 == 0 {
            (ms2_peaks[mid - 1] + ms2_peaks[mid]) as f64 / 2.0
        } else {
            ms2_peaks[mid] as f64
        }
    };

    Ok(RunStats {
        total_spectra: spectra.len(),
        ms1_count,
        ms2_count,
        msn_count,
        min_rt,
        max_rt,
        mean_peaks_per_ms2: mean_peaks,
        median_peaks_per_ms2: median_peaks,
        mean_tic: total_tic / spectra.len() as f64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_peaks(data: &[(f64, f64)]) -> Vec<Peak> {
        data.iter().map(|&(mz, intensity)| Peak { mz, intensity }).collect()
    }

    #[test]
    fn test_spectrum_creation() {
        let peaks = make_peaks(&[(100.0, 1000.0), (200.0, 500.0), (150.0, 2000.0)]);
        let spec = MassSpectrum::new("scan1", MsLevel::Ms1, 60.0, peaks).unwrap();
        assert_eq!(spec.num_peaks(), 3);
        assert_eq!(spec.base_peak_mz, 150.0);
        assert_eq!(spec.base_peak_intensity, 2000.0);
        assert_eq!(spec.tic, 3500.0);
        // Peaks should be sorted by m/z
        assert_eq!(spec.peaks[0].mz, 100.0);
        assert_eq!(spec.peaks[1].mz, 150.0);
        assert_eq!(spec.peaks[2].mz, 200.0);
    }

    #[test]
    fn test_empty_spectrum_error() {
        let result = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_negative_mz_error() {
        let peaks = make_peaks(&[(-1.0, 100.0)]);
        let result = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, peaks);
        assert!(result.is_err());
    }

    #[test]
    fn test_find_peak() {
        let peaks = make_peaks(&[(100.0, 1000.0), (200.0, 500.0), (300.0, 750.0)]);
        let spec = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, peaks).unwrap();

        let found = spec.find_peak(199.5, 1.0);
        assert!(found.is_some());
        assert_eq!(found.unwrap().mz, 200.0);

        let not_found = spec.find_peak(250.0, 1.0);
        assert!(not_found.is_none());
    }

    #[test]
    fn test_filter_by_relative_intensity() {
        let peaks = make_peaks(&[(100.0, 1000.0), (200.0, 50.0), (300.0, 500.0)]);
        let mut spec = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, peaks).unwrap();
        spec.filter_by_relative_intensity(0.1); // keep >= 10% of base peak (1000)
        assert_eq!(spec.num_peaks(), 2);
        assert_eq!(spec.peaks[0].mz, 100.0);
        assert_eq!(spec.peaks[1].mz, 300.0);
    }

    #[test]
    fn test_top_n() {
        let peaks = make_peaks(&[(100.0, 300.0), (200.0, 100.0), (300.0, 500.0), (400.0, 200.0)]);
        let mut spec = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, peaks).unwrap();
        spec.top_n(2);
        assert_eq!(spec.num_peaks(), 2);
        assert_eq!(spec.peaks[0].mz, 100.0);
        assert_eq!(spec.peaks[1].mz, 300.0);
    }

    #[test]
    fn test_normalize() {
        let peaks = make_peaks(&[(100.0, 1000.0), (200.0, 500.0)]);
        let mut spec = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, peaks).unwrap();
        spec.normalize();
        assert_eq!(spec.peaks[0].intensity, 100.0);
        assert_eq!(spec.peaks[1].intensity, 50.0);
    }

    #[test]
    fn test_sqrt_transform() {
        let peaks = make_peaks(&[(100.0, 100.0), (200.0, 400.0)]);
        let mut spec = MassSpectrum::new("scan1", MsLevel::Ms1, 0.0, peaks).unwrap();
        spec.sqrt_transform();
        assert!((spec.peaks[0].intensity - 10.0).abs() < 1e-10);
        assert!((spec.peaks[1].intensity - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_deisotope() {
        let peaks = make_peaks(&[(100.0, 1000.0), (100.5, 300.0), (101.0, 100.0), (200.0, 500.0)]);
        let mut spec = MassSpectrum::new("scan1", MsLevel::Ms2, 0.0, peaks).unwrap();
        spec.deisotope(1.1);
        assert_eq!(spec.num_peaks(), 2);
        assert_eq!(spec.peaks[0].mz, 100.0); // most intense in first cluster
        assert_eq!(spec.peaks[1].mz, 200.0);
    }

    #[test]
    fn test_ms_level() {
        assert_eq!(MsLevel::from_u8(1), MsLevel::Ms1);
        assert_eq!(MsLevel::from_u8(2), MsLevel::Ms2);
        assert_eq!(MsLevel::from_u8(3), MsLevel::Ms3);
        assert_eq!(MsLevel::from_u8(5), MsLevel::MsN(5));
        assert_eq!(MsLevel::Ms2.as_u8(), 2);
    }

    #[test]
    fn test_with_precursor() {
        let peaks = make_peaks(&[(100.0, 1000.0)]);
        let spec = MassSpectrum::new("scan1", MsLevel::Ms2, 30.0, peaks)
            .unwrap()
            .with_precursor(Precursor {
                mz: 500.0,
                charge: Some(2),
                intensity: Some(5000.0),
                isolation_width: Some(2.0),
                fragmentation: FragmentationMethod::HCD,
            });
        assert!(spec.precursor.is_some());
        let pre = spec.precursor.unwrap();
        assert_eq!(pre.mz, 500.0);
        assert_eq!(pre.charge, Some(2));
        assert_eq!(pre.fragmentation, FragmentationMethod::HCD);
    }

    #[test]
    fn test_run_stats() {
        let ms1 = MassSpectrum::new("s1", MsLevel::Ms1, 10.0, make_peaks(&[(100.0, 1000.0)])).unwrap();
        let ms2a = MassSpectrum::new("s2", MsLevel::Ms2, 15.0, make_peaks(&[(100.0, 500.0), (200.0, 300.0)])).unwrap();
        let ms2b = MassSpectrum::new("s3", MsLevel::Ms2, 20.0, make_peaks(&[(100.0, 200.0), (200.0, 300.0), (300.0, 100.0), (400.0, 50.0)])).unwrap();

        let stats = run_stats(&[ms1, ms2a, ms2b]).unwrap();
        assert_eq!(stats.total_spectra, 3);
        assert_eq!(stats.ms1_count, 1);
        assert_eq!(stats.ms2_count, 2);
        assert_eq!(stats.min_rt, 10.0);
        assert_eq!(stats.max_rt, 20.0);
        assert_eq!(stats.mean_peaks_per_ms2, 3.0); // (2+4)/2
    }

    #[test]
    fn test_run_stats_empty() {
        assert!(run_stats(&[]).is_err());
    }
}
