use cyanea_proteomics::spectrum::*;
use cyanea_proteomics::mgf::{parse_mgf, write_mgf, mgf_stats};
use cyanea_proteomics::mzml::{parse_mzml_text, mzml_stats};

#[test]
fn test_spectrum_processing_pipeline() {
    let peaks: Vec<Peak> = (0..100)
        .map(|i| Peak {
            mz: 100.0 + i as f64 * 2.0,
            intensity: if i % 10 == 0 { 10000.0 } else { 100.0 },
        })
        .collect();

    let mut spec = MassSpectrum::new("pipeline_test", MsLevel::Ms2, 120.0, peaks).unwrap();
    assert_eq!(spec.num_peaks(), 100);

    // Filter low-intensity
    spec.filter_by_relative_intensity(0.05);
    assert!(spec.num_peaks() < 100);

    // Top N
    spec.top_n(10);
    assert!(spec.num_peaks() <= 10);

    // Normalize
    spec.normalize();
    assert!((spec.base_peak_intensity - 100.0).abs() < 1e-10);
}

#[test]
fn test_mgf_multiple_spectra() {
    let mgf = "\
BEGIN IONS
TITLE=scan1
PEPMASS=400.0 2000
CHARGE=2+
RTINSECONDS=60
100 500
200 1000
300 750
END IONS

BEGIN IONS
TITLE=scan2
PEPMASS=500.0 3000
CHARGE=3+
RTINSECONDS=120
150 800
250 600
350 400
450 200
END IONS

BEGIN IONS
TITLE=scan3
PEPMASS=600.0
CHARGE=2+
RTINSECONDS=180
200 300
END IONS
";

    let spectra = parse_mgf(mgf).unwrap();
    assert_eq!(spectra.len(), 3);

    let stats = mgf_stats(&spectra).unwrap();
    assert_eq!(stats.num_spectra, 3);
    assert_eq!(stats.with_charge, 3);
    assert!((stats.min_rt - 60.0).abs() < 1e-10);
    assert!((stats.max_rt - 180.0).abs() < 1e-10);

    // Roundtrip
    let written = write_mgf(&spectra);
    let reparsed = parse_mgf(&written).unwrap();
    assert_eq!(reparsed.len(), 3);
    for (orig, re) in spectra.iter().zip(reparsed.iter()) {
        assert_eq!(orig.id, re.id);
        assert_eq!(orig.num_peaks(), re.num_peaks());
    }
}

#[test]
fn test_mzml_mixed_levels() {
    let mzml = "\
SPECTRUM id=ms1_1 ms_level=1 rt=10.0
PEAKS
200.0 50000
400.0 30000
600.0 20000
END

SPECTRUM id=ms2_1 ms_level=2 rt=12.0
PRECURSOR mz=400.25 charge=2 fragmentation=HCD
PEAKS
100.0 1000
200.0 800
300.0 500
END

SPECTRUM id=ms2_2 ms_level=2 rt=14.0
PRECURSOR mz=600.5 charge=3 fragmentation=CID
PEAKS
150.0 600
250.0 400
END
";

    let spectra = parse_mzml_text(mzml).unwrap();
    assert_eq!(spectra.len(), 3);

    let stats = mzml_stats(&spectra).unwrap();
    assert_eq!(stats.ms1_count, 1);
    assert_eq!(stats.ms2_count, 2);
    assert!((stats.min_rt - 10.0).abs() < 1e-10);
    assert!((stats.max_rt - 14.0).abs() < 1e-10);
}

#[test]
fn test_run_stats_comprehensive() {
    let mut spectra = Vec::new();

    // Generate 10 MS1 and 20 MS2 spectra
    for i in 0..10 {
        let peaks = vec![
            Peak { mz: 100.0 + i as f64, intensity: 1000.0 * (i + 1) as f64 },
        ];
        spectra.push(MassSpectrum::new(
            format!("ms1_{}", i), MsLevel::Ms1, i as f64 * 10.0, peaks,
        ).unwrap());
    }

    for i in 0..20 {
        let num_peaks = (i % 5) + 2;
        let peaks: Vec<Peak> = (0..num_peaks)
            .map(|j| Peak { mz: 100.0 + j as f64 * 50.0, intensity: 500.0 })
            .collect();
        spectra.push(MassSpectrum::new(
            format!("ms2_{}", i), MsLevel::Ms2, 5.0 + i as f64 * 5.0, peaks,
        ).unwrap());
    }

    let stats = run_stats(&spectra).unwrap();
    assert_eq!(stats.total_spectra, 30);
    assert_eq!(stats.ms1_count, 10);
    assert_eq!(stats.ms2_count, 20);
    assert!(stats.mean_peaks_per_ms2 > 0.0);
}
