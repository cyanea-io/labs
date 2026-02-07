//! Python bindings for the Cyanea bioinformatics ecosystem.
//!
//! Submodules:
//!
//! - `cyanea.seq` — Sequence types (DNA, RNA, Protein) and FASTA/FASTQ parsing
//! - `cyanea.align` — Pairwise sequence alignment
//! - `cyanea.stats` — Descriptive statistics and hypothesis testing
//! - `cyanea.core` — SHA-256 hashing and zstd compression

mod align;
mod core_utils;
mod error;
mod seq;
mod stats;

use pyo3::prelude::*;

#[pymodule]
fn cyanea(m: &Bound<'_, PyModule>) -> PyResult<()> {
    seq::register(m)?;
    align::register(m)?;
    stats::register(m)?;
    core_utils::register(m)?;

    // Register submodules in sys.modules so `from cyanea.seq import DnaSequence` works.
    let py = m.py();
    let sys = py.import("sys")?;
    let modules = sys.getattr("modules")?;
    for name in &["cyanea.seq", "cyanea.align", "cyanea.stats", "cyanea.core"] {
        let submod = m.getattr(name.rsplit('.').next().unwrap())?;
        modules.set_item(name, submod)?;
    }

    Ok(())
}
