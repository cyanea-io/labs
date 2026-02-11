//! Python bindings for the Cyanea bioinformatics ecosystem.
//!
//! Submodules:
//!
//! - `cyanea.seq` — Sequence types (DNA, RNA, Protein) and FASTA/FASTQ parsing
//! - `cyanea.align` — Pairwise sequence alignment
//! - `cyanea.stats` — Descriptive statistics and hypothesis testing
//! - `cyanea.core` — SHA-256 hashing and zstd compression
//! - `cyanea.ml` — Distance metrics
//! - `cyanea.chem` — Chemistry and small-molecule analysis
//! - `cyanea.struct_bio` — Protein 3D structure analysis
//! - `cyanea.phylo` — Phylogenetic trees and distances
//! - `cyanea.io` — Bioinformatics file format parsing

mod align;
mod chem;
mod core_utils;
mod error;
mod io;
mod ml;
mod phylo;
mod seq;
mod stats;
mod struct_bio;

use pyo3::prelude::*;

#[pymodule]
fn cyanea(m: &Bound<'_, PyModule>) -> PyResult<()> {
    seq::register(m)?;
    align::register(m)?;
    stats::register(m)?;
    core_utils::register(m)?;
    ml::register(m)?;
    chem::register(m)?;
    struct_bio::register(m)?;
    phylo::register(m)?;
    io::register(m)?;

    // Register submodules in sys.modules so `from cyanea.seq import DnaSequence` works.
    let py = m.py();
    let sys = py.import("sys")?;
    let modules = sys.getattr("modules")?;
    for name in &[
        "cyanea.seq",
        "cyanea.align",
        "cyanea.stats",
        "cyanea.core",
        "cyanea.ml",
        "cyanea.chem",
        "cyanea.struct_bio",
        "cyanea.phylo",
        "cyanea.io",
    ] {
        let submod = m.getattr(name.rsplit('.').next().unwrap())?;
        modules.set_item(name, submod)?;
    }

    Ok(())
}
