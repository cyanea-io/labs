//! Structured protocol templates for wet-lab and dry-lab workflows.
//!
//! Each template provides title, category, description, materials/software,
//! numbered steps with timing and tips, and expected outputs. Designed for
//! rendering in the Cyanea notebook platform as "open and run" experiences.

/// A protocol template with structured steps and metadata.
#[derive(Debug, Clone)]
pub struct Protocol {
    pub title: &'static str,
    pub slug: &'static str,
    pub category: ProtocolCategory,
    pub description: &'static str,
    pub estimated_time: &'static str,
    pub difficulty: Difficulty,
    pub requirements: Vec<&'static str>,
    pub steps: Vec<ProtocolStep>,
    pub expected_outputs: Vec<&'static str>,
    pub references: Vec<&'static str>,
}

/// Protocol category (wet lab vs dry lab).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProtocolCategory {
    WetLab,
    DryLab,
}

/// Protocol difficulty level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Difficulty {
    Beginner,
    Intermediate,
    Advanced,
}

/// A single step in a protocol.
#[derive(Debug, Clone)]
pub struct ProtocolStep {
    pub number: u32,
    pub title: &'static str,
    pub description: &'static str,
    pub duration: Option<&'static str>,
    pub tips: Vec<&'static str>,
    pub caution: Option<&'static str>,
}

impl Protocol {
    /// Render protocol as markdown.
    pub fn to_markdown(&self) -> String {
        let mut md = String::new();
        md.push_str(&format!("# {}\n\n", self.title));
        md.push_str(&format!("**Category:** {}\n", match self.category {
            ProtocolCategory::WetLab => "Wet Lab",
            ProtocolCategory::DryLab => "Dry Lab (Computational)",
        }));
        md.push_str(&format!("**Difficulty:** {}\n", match self.difficulty {
            Difficulty::Beginner => "Beginner",
            Difficulty::Intermediate => "Intermediate",
            Difficulty::Advanced => "Advanced",
        }));
        md.push_str(&format!("**Estimated time:** {}\n\n", self.estimated_time));
        md.push_str(&format!("{}\n\n", self.description));

        md.push_str("## Requirements\n\n");
        for r in &self.requirements {
            md.push_str(&format!("- {}\n", r));
        }
        md.push('\n');

        md.push_str("## Steps\n\n");
        for step in &self.steps {
            md.push_str(&format!("### Step {}: {}\n\n", step.number, step.title));
            if let Some(dur) = step.duration {
                md.push_str(&format!("*Duration: {}*\n\n", dur));
            }
            md.push_str(&format!("{}\n\n", step.description));
            if let Some(caution) = step.caution {
                md.push_str(&format!("> **Caution:** {}\n\n", caution));
            }
            for tip in &step.tips {
                md.push_str(&format!("- **Tip:** {}\n", tip));
            }
            if !step.tips.is_empty() {
                md.push('\n');
            }
        }

        md.push_str("## Expected Outputs\n\n");
        for o in &self.expected_outputs {
            md.push_str(&format!("- {}\n", o));
        }
        md.push('\n');

        if !self.references.is_empty() {
            md.push_str("## References\n\n");
            for (i, r) in self.references.iter().enumerate() {
                md.push_str(&format!("{}. {}\n", i + 1, r));
            }
        }

        md
    }
}

/// List all available protocol templates.
pub fn all_protocols() -> Vec<Protocol> {
    vec![
        elisa(),
        qpcr(),
        immunofluorescence(),
        cell_viability_mtt(),
        site_directed_mutagenesis(),
        bacterial_transformation(),
        coimmunoprecipitation(),
        chipseq_library_prep(),
        atacseq_library_prep(),
        hic_library_prep(),
        gwas_pipeline(),
        longread_genome_assembly(),
        proteomics_dia_search(),
        spatial_transcriptomics_analysis(),
        hic_analysis(),
        alphafold_structure_prediction(),
    ]
}

/// List only wet-lab protocols.
pub fn wet_lab_protocols() -> Vec<Protocol> {
    all_protocols().into_iter().filter(|p| p.category == ProtocolCategory::WetLab).collect()
}

/// List only dry-lab (computational) protocols.
pub fn dry_lab_protocols() -> Vec<Protocol> {
    all_protocols().into_iter().filter(|p| p.category == ProtocolCategory::DryLab).collect()
}

// ---------------------------------------------------------------------------
// Wet lab protocols (10)
// ---------------------------------------------------------------------------

/// ELISA (Enzyme-Linked Immunosorbent Assay) protocol.
pub fn elisa() -> Protocol {
    Protocol {
        title: "ELISA (Enzyme-Linked Immunosorbent Assay)",
        slug: "elisa",
        category: ProtocolCategory::WetLab,
        description: "Sandwich ELISA for quantitative detection of a target protein in cell lysates, serum, or culture supernatants. Uses a capture antibody, detection antibody, and HRP-conjugated secondary for colorimetric readout.",
        estimated_time: "2 days (overnight coating + 6h assay)",
        difficulty: Difficulty::Beginner,
        requirements: vec![
            "96-well ELISA plate (high-binding)",
            "Capture antibody (target-specific)",
            "Detection antibody (biotinylated)",
            "Streptavidin-HRP conjugate",
            "TMB substrate solution",
            "Stop solution (2N H2SO4)",
            "Wash buffer (PBS + 0.05% Tween-20)",
            "Blocking buffer (PBS + 1% BSA)",
            "Plate reader (450 nm)",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Coat plate",
                description: "Dilute capture antibody to 1-10 ug/mL in coating buffer (0.1M sodium bicarbonate, pH 9.6). Add 100 uL/well. Seal and incubate overnight at 4C.",
                duration: Some("Overnight (12-16h)"),
                tips: vec!["Use multichannel pipette for consistency"],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Block",
                description: "Wash plate 3x with wash buffer (300 uL/well). Add 200 uL/well blocking buffer. Incubate 1h at room temperature.",
                duration: Some("1h"),
                tips: vec!["Tap plate firmly on paper towels between washes"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Add standards and samples",
                description: "Prepare serial dilutions of recombinant protein standard (7-point, 2-fold). Add 100 uL/well of standards and diluted samples in duplicate. Incubate 2h at room temperature.",
                duration: Some("2h"),
                tips: vec![
                    "Include blank wells (diluent only)",
                    "Pre-dilute samples if concentration is expected to be high",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Detection antibody",
                description: "Wash 3x. Add 100 uL/well biotinylated detection antibody (per manufacturer). Incubate 1h at room temperature.",
                duration: Some("1h"),
                tips: vec![],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Streptavidin-HRP",
                description: "Wash 3x. Add 100 uL/well streptavidin-HRP (1:200 dilution). Incubate 30 min at room temperature.",
                duration: Some("30 min"),
                tips: vec!["Protect from light after this step"],
                caution: None,
            },
            ProtocolStep {
                number: 6,
                title: "Develop and read",
                description: "Wash 5x. Add 100 uL/well TMB substrate. Incubate 15-30 min in the dark until blue color develops. Add 50 uL/well stop solution. Read absorbance at 450 nm within 30 min.",
                duration: Some("30-60 min"),
                tips: vec!["Add stop solution in same order as TMB for consistent timing"],
                caution: Some("H2SO4 is corrosive — wear gloves and eye protection"),
            },
        ],
        expected_outputs: vec![
            "Standard curve (OD450 vs concentration)",
            "Sample concentrations interpolated from standard curve",
            "CV% < 10% between duplicates",
        ],
        references: vec![
            "Engvall E, Perlmann P. Immunochemistry. 1971;8(9):871-874.",
        ],
    }
}

/// qPCR / RT-qPCR protocol.
pub fn qpcr() -> Protocol {
    Protocol {
        title: "qPCR / RT-qPCR (Quantitative PCR)",
        slug: "qpcr-rt-qpcr",
        category: ProtocolCategory::WetLab,
        description: "Real-time quantitative PCR for gene expression analysis using SYBR Green or TaqMan chemistry. Includes reverse transcription for RNA targets.",
        estimated_time: "4-6 hours",
        difficulty: Difficulty::Beginner,
        requirements: vec![
            "RNA samples (50-500 ng total RNA per reaction)",
            "Reverse transcriptase kit (for RT-qPCR)",
            "SYBR Green or TaqMan master mix",
            "Forward and reverse primers (10 uM stocks)",
            "96-well qPCR plates or strips",
            "Real-time PCR instrument",
            "Nuclease-free water",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Reverse transcription (RT-qPCR only)",
                description: "Combine 500 ng RNA with oligo(dT) or random hexamer primers. Heat 65C 5 min, cool on ice. Add RT enzyme, buffer, dNTPs. Incubate: 25C 10 min, 42C 60 min, 70C 15 min. Dilute cDNA 1:5.",
                duration: Some("90 min"),
                tips: vec![
                    "Include no-RT control to detect genomic DNA contamination",
                    "Use same RNA input across all samples",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Primer validation",
                description: "Run standard curve with 5-point serial dilution (1:4) of pooled cDNA. Calculate efficiency: E = 10^(-1/slope) - 1. Acceptable range: 90-110%. Check melt curve for single peak (SYBR).",
                duration: Some("2h"),
                tips: vec!["Optimal primer Tm: 58-62C, amplicon 80-200 bp"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Set up qPCR reactions",
                description: "Per 20 uL reaction: 10 uL 2x master mix, 0.5 uL each primer (10 uM), 2 uL cDNA template, 7 uL water. Prepare master mix for all replicates. Load in triplicate.",
                duration: Some("30 min"),
                tips: vec![
                    "Use a multichannel pipette and reagent reservoirs",
                    "Include no-template control (NTC) for each primer pair",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Run thermal cycling",
                description: "Program: 95C 10 min (hot start), then 40 cycles of 95C 15s / 60C 60s. Add melt curve: 65C-95C ramp at 0.5C/5s (SYBR only).",
                duration: Some("90 min"),
                tips: vec!["Seal plate firmly to prevent evaporation"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Analyze results",
                description: "Set threshold in exponential phase. Export Ct values. Calculate relative expression using delta-delta Ct method: fold change = 2^(-ddCt). Normalize to reference gene (GAPDH, ACTB, or 18S).",
                duration: Some("30 min"),
                tips: vec![
                    "Ct SD among triplicates should be < 0.5",
                    "NTC Ct should be > 35 or undetected",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Ct values for target and reference genes",
            "Amplification curves and melt curves",
            "Fold-change table (delta-delta Ct)",
        ],
        references: vec![
            "Livak KJ, Schmittgen TD. Methods. 2001;25(4):402-408.",
        ],
    }
}

/// Immunohistochemistry / Immunofluorescence protocol.
pub fn immunofluorescence() -> Protocol {
    Protocol {
        title: "Immunofluorescence (IF) Staining",
        slug: "immunofluorescence",
        category: ProtocolCategory::WetLab,
        description: "Fluorescent antibody staining of fixed cells or tissue sections for protein localization. Supports direct (fluorophore-conjugated primary) or indirect (primary + fluorescent secondary) methods.",
        estimated_time: "2 days (overnight primary antibody)",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "Fixed cells on coverslips or tissue sections on slides",
            "Primary antibody (target-specific)",
            "Fluorophore-conjugated secondary antibody (e.g., Alexa Fluor 488/594)",
            "DAPI nuclear counterstain",
            "Blocking buffer (5% normal serum + 0.3% Triton X-100 in PBS)",
            "Mounting medium (anti-fade)",
            "Fluorescence or confocal microscope",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Fixation",
                description: "Fix cells in 4% paraformaldehyde (PFA) in PBS for 15 min at room temperature. Wash 3x5 min with PBS.",
                duration: Some("30 min"),
                tips: vec![],
                caution: Some("PFA is a fixative and irritant — use in fume hood with gloves"),
            },
            ProtocolStep {
                number: 2,
                title: "Permeabilization and blocking",
                description: "Incubate in blocking buffer (5% normal serum from secondary antibody host species + 0.3% Triton X-100 in PBS) for 1h at room temperature.",
                duration: Some("1h"),
                tips: vec!["Use serum matching the secondary antibody host to reduce background"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Primary antibody",
                description: "Dilute primary antibody in antibody dilution buffer (1% BSA, 0.3% Triton X-100 in PBS). Apply to samples and incubate overnight at 4C in a humidified chamber.",
                duration: Some("Overnight"),
                tips: vec![
                    "Optimal dilution varies (typically 1:100 to 1:1000) — titrate",
                    "Include a no-primary control",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Secondary antibody",
                description: "Wash 3x5 min with PBS. Apply fluorophore-conjugated secondary antibody (1:500-1:1000) in antibody dilution buffer. Incubate 1-2h at room temperature in the dark.",
                duration: Some("1-2h"),
                tips: vec!["Protect from light from this step onward"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Counterstain and mount",
                description: "Wash 3x5 min with PBS. Incubate with DAPI (1 ug/mL) for 5 min. Wash 2x5 min. Mount with anti-fade mounting medium. Seal edges with nail polish. Allow to cure 24h at 4C.",
                duration: Some("30 min"),
                tips: vec!["Store slides at 4C protected from light, image within 1 week"],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Fluorescence images showing protein localization",
            "Merged images with nuclear (DAPI) and protein channels",
            "Quantification: % positive cells, mean fluorescence intensity",
        ],
        references: vec![
            "Im K, et al. Methods Mol Biol. 2019;1897:299-311.",
        ],
    }
}

/// MTT / Cell viability assay protocol.
pub fn cell_viability_mtt() -> Protocol {
    Protocol {
        title: "MTT Cell Viability Assay",
        slug: "mtt-cell-viability",
        category: ProtocolCategory::WetLab,
        description: "Colorimetric assay measuring mitochondrial metabolic activity as a proxy for cell viability and proliferation. MTT tetrazolium is reduced to purple formazan by live cells.",
        estimated_time: "4-5 hours (+ overnight cell seeding)",
        difficulty: Difficulty::Beginner,
        requirements: vec![
            "Cells seeded in 96-well plate (5,000-10,000 cells/well)",
            "MTT stock solution (5 mg/mL in PBS, sterile-filtered)",
            "Solubilization solution (DMSO or 10% SDS in 0.01N HCl)",
            "Plate reader (570 nm, reference 690 nm)",
            "Multichannel pipette",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Seed cells",
                description: "Seed cells at 5,000-10,000 cells/well in 100 uL complete medium in a 96-well plate. Include blank wells (medium only). Incubate overnight at 37C, 5% CO2.",
                duration: Some("Overnight"),
                tips: vec!["Use inner 60 wells; fill outer wells with PBS to reduce edge effects"],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Treat cells",
                description: "Add treatment compounds at desired concentrations. Include vehicle control. Use at least 6 replicates per condition. Incubate for desired time (24-72h).",
                duration: Some("24-72h"),
                tips: vec!["Prepare serial dilutions of drug for dose-response curves"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Add MTT",
                description: "Add 10 uL of 5 mg/mL MTT stock per well (final 0.5 mg/mL). Incubate 2-4h at 37C. Purple formazan crystals should be visible under microscope.",
                duration: Some("2-4h"),
                tips: vec!["MTT is light-sensitive — wrap in foil"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Solubilize and read",
                description: "Remove medium carefully. Add 100 uL DMSO per well. Pipette up and down to dissolve crystals. Incubate 10 min at room temperature with gentle shaking. Read absorbance at 570 nm (reference 690 nm).",
                duration: Some("30 min"),
                tips: vec!["Read within 1h of adding DMSO"],
                caution: Some("DMSO penetrates skin — wear nitrile gloves"),
            },
        ],
        expected_outputs: vec![
            "Absorbance values (570-690 nm) per well",
            "% viability = (OD_treated / OD_control) x 100",
            "Dose-response curve and IC50 value",
        ],
        references: vec![
            "Mosmann T. J Immunol Methods. 1983;65(1-2):55-63.",
        ],
    }
}

/// Site-directed mutagenesis protocol.
pub fn site_directed_mutagenesis() -> Protocol {
    Protocol {
        title: "Site-Directed Mutagenesis",
        slug: "site-directed-mutagenesis",
        category: ProtocolCategory::WetLab,
        description: "PCR-based introduction of specific point mutations, insertions, or deletions into a plasmid. Uses overlapping mutagenic primers and DpnI digestion of parental methylated DNA.",
        estimated_time: "2 days",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "Template plasmid (50-100 ng, from dam+ E. coli)",
            "Mutagenic primers (forward and reverse, 25-45 nt with mutation centered)",
            "High-fidelity DNA polymerase (e.g., Q5, Pfu Turbo)",
            "dNTP mix (10 mM each)",
            "DpnI restriction enzyme",
            "Competent E. coli cells (DH5a or similar)",
            "LB agar plates with appropriate antibiotic",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Design mutagenic primers",
                description: "Design complementary primers (25-45 bp) with the desired mutation in the center, flanked by >= 10 bp of matching sequence on each side. Tm >= 78C recommended. Primers should overlap by >= 20 bp.",
                duration: Some("30 min"),
                tips: vec![
                    "GC clamp at 3' end improves annealing",
                    "For deletions, omit the target sequence from the primer",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "PCR amplification",
                description: "25 uL reaction: 50 ng template, 125 ng each primer, 200 uM dNTPs, 1x buffer, polymerase. Cycling: 95C 30s, then 18 cycles of 95C 30s / 55C 60s / 68C 1 min/kb plasmid.",
                duration: Some("2-3h"),
                tips: vec!["Use 12-18 cycles to minimize secondary mutations"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "DpnI digestion",
                description: "Add 1 uL DpnI (20 U) directly to the PCR product. Incubate 1h at 37C to digest methylated parental DNA.",
                duration: Some("1h"),
                tips: vec!["DpnI only cuts dam-methylated (parental) DNA"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Transform and screen",
                description: "Transform 2-5 uL of DpnI-treated product into 50 uL competent cells. Heat shock 42C 45s. Recover in SOC 1h at 37C. Plate on selective agar. Pick 3-4 colonies for sequencing.",
                duration: Some("Overnight + 1 day sequencing"),
                tips: vec!["Expect 50-80% mutation efficiency with this method"],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Colonies on selective plate (typically 10-100)",
            "Sanger sequencing confirming desired mutation",
            "Confirmed mutant plasmid prep for downstream use",
        ],
        references: vec![
            "Zheng L, et al. Nucleic Acids Res. 2004;32(14):e115.",
        ],
    }
}

/// Bacterial transformation protocol.
pub fn bacterial_transformation() -> Protocol {
    Protocol {
        title: "Bacterial Transformation (Heat Shock)",
        slug: "bacterial-transformation",
        category: ProtocolCategory::WetLab,
        description: "Introduction of plasmid DNA into chemically competent E. coli by heat shock. Standard cloning workflow step for amplifying and maintaining recombinant DNA.",
        estimated_time: "2 hours + overnight growth",
        difficulty: Difficulty::Beginner,
        requirements: vec![
            "Chemically competent E. coli (DH5a, TOP10, or BL21 for expression)",
            "Plasmid DNA (1-100 ng) or ligation mix",
            "SOC medium (pre-warmed to 37C)",
            "LB agar plates with selective antibiotic",
            "Ice bucket, 42C water bath, 37C shaking incubator",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Thaw competent cells",
                description: "Thaw one aliquot (50-100 uL) of competent cells on ice for 10-20 min. Do not vortex.",
                duration: Some("20 min"),
                tips: vec!["Keep cells on ice at all times until heat shock"],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Add DNA",
                description: "Add 1-5 uL plasmid DNA (1-100 ng) or ligation product. Gently flick tube to mix. Incubate on ice for 30 min.",
                duration: Some("30 min"),
                tips: vec![
                    "Less is more — excess DNA decreases efficiency",
                    "Include a positive control (known plasmid) and negative control (no DNA)",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Heat shock",
                description: "Heat shock at exactly 42C for 45 seconds. Immediately return to ice for 2 min.",
                duration: Some("3 min"),
                tips: vec!["Timing is critical — use a timer"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Recovery",
                description: "Add 250-950 uL pre-warmed SOC medium. Incubate at 37C with shaking (225 rpm) for 1h.",
                duration: Some("1h"),
                tips: vec!["Recovery is essential for antibiotic resistance expression"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Plate and grow",
                description: "Plate 50-200 uL on pre-warmed LB-antibiotic agar plates. For low-efficiency transformations, pellet cells and plate entire volume. Incubate overnight at 37C.",
                duration: Some("Overnight"),
                tips: vec!["Spread evenly using glass beads or a sterile spreader"],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Colonies on selective plate (10-1000+ depending on DNA amount)",
            "Transformation efficiency: CFU/ug DNA",
            "No colonies on negative control plate",
        ],
        references: vec![
            "Hanahan D. J Mol Biol. 1983;166(4):557-580.",
        ],
    }
}

/// Co-immunoprecipitation (Co-IP) protocol.
pub fn coimmunoprecipitation() -> Protocol {
    Protocol {
        title: "Co-Immunoprecipitation (Co-IP)",
        slug: "co-immunoprecipitation",
        category: ProtocolCategory::WetLab,
        description: "Identification of protein-protein interactions by immunoprecipitating a bait protein and detecting co-precipitated binding partners by western blot.",
        estimated_time: "2 days",
        difficulty: Difficulty::Advanced,
        requirements: vec![
            "Cell lysate (1-5 mg total protein)",
            "IP antibody (bait protein-specific, 2-5 ug)",
            "Protein A/G magnetic beads (25-50 uL slurry)",
            "Lysis buffer (50 mM Tris pH 7.4, 150 mM NaCl, 1% NP-40, protease inhibitors)",
            "Wash buffer (lysis buffer with 0.1% NP-40)",
            "SDS-PAGE and western blot reagents",
            "Detection antibodies for bait and prey proteins",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Prepare cell lysate",
                description: "Lyse 1-5 x 10^7 cells in ice-cold lysis buffer (1 mL per 10^7 cells). Add protease inhibitor cocktail fresh. Incubate on ice 30 min with occasional vortexing. Clarify by centrifugation at 14,000g, 15 min, 4C.",
                duration: Some("1h"),
                tips: vec!["Save 50 uL as 'input' control (10% of total)"],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Pre-clear lysate",
                description: "Add 25 uL Protein A/G beads to lysate. Rotate 1h at 4C. Remove beads by magnetic separation. This reduces non-specific binding.",
                duration: Some("1h"),
                tips: vec!["Pre-clearing is optional but reduces background"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Immunoprecipitation",
                description: "Add IP antibody (2-5 ug) to pre-cleared lysate. Rotate overnight at 4C. Include IgG isotype control at same concentration.",
                duration: Some("Overnight"),
                tips: vec!["Cross-linking antibody to beads (DMP) reduces IgG heavy chain in blot"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Capture and wash",
                description: "Add 50 uL Protein A/G beads. Rotate 2h at 4C. Wash beads 4x with wash buffer (5 min each at 4C). Keep samples cold.",
                duration: Some("2.5h"),
                tips: vec!["Increasing NaCl to 300 mM in wash reduces background at the cost of weaker interactions"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Elute and detect",
                description: "Elute by adding 30 uL 2x SDS sample buffer and heating at 95C for 5 min. Load eluates and input on SDS-PAGE gel. Transfer and blot for bait and prey proteins.",
                duration: Some("4-6h"),
                tips: vec!["Run IP and IgG control side by side for comparison"],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Western blot showing bait protein in IP lane (confirms pulldown)",
            "Prey protein band in IP lane but not in IgG control",
            "Input lane showing both proteins in starting material",
        ],
        references: vec![
            "Bhatt DK, et al. Expert Rev Proteomics. 2020;17(4):301-317.",
        ],
    }
}

/// ChIP-seq library preparation protocol.
pub fn chipseq_library_prep() -> Protocol {
    Protocol {
        title: "ChIP-seq Library Preparation",
        slug: "chipseq-library-prep",
        category: ProtocolCategory::WetLab,
        description: "Chromatin immunoprecipitation followed by sequencing to map genome-wide protein-DNA interactions. Covers cross-linking, sonication, immunoprecipitation, and library construction.",
        estimated_time: "3-4 days",
        difficulty: Difficulty::Advanced,
        requirements: vec![
            "10-20 million cells per IP",
            "1% formaldehyde (methanol-free) for cross-linking",
            "Glycine (1.25 M stock) to quench",
            "Sonication device (Bioruptor or Covaris)",
            "ChIP-grade antibody (2-5 ug per IP)",
            "Protein A/G magnetic beads",
            "Proteinase K, RNase A",
            "DNA purification columns or beads",
            "Library prep kit (NEBNext Ultra II or similar)",
            "AMPure XP beads for size selection",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Cross-link chromatin",
                description: "Add formaldehyde to 1% final concentration to cells in medium. Incubate 10 min at room temperature with rotation. Quench with 125 mM glycine for 5 min. Wash 2x with cold PBS. Pellet cells.",
                duration: Some("30 min"),
                tips: vec!["Cross-linking time is critical — over-fixation reduces sonication efficiency"],
                caution: Some("Formaldehyde is toxic and carcinogenic — use in fume hood"),
            },
            ProtocolStep {
                number: 2,
                title: "Sonicate chromatin",
                description: "Resuspend pellet in lysis buffer with protease inhibitors. Sonicate to generate 200-500 bp fragments. Verify size on agarose gel or Bioanalyzer. Clarify by centrifugation at 16,000g, 10 min, 4C.",
                duration: Some("1-2h"),
                tips: vec![
                    "Optimize sonication cycles for your cell type",
                    "Target 200-300 bp for TF ChIP, 200-500 bp for histone ChIP",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Immunoprecipitation",
                description: "Save 5% of sonicated chromatin as input control. Pre-block beads with BSA. Add antibody to chromatin and rotate overnight at 4C. Add beads and rotate 2h. Wash sequentially: low salt, high salt, LiCl, TE (2x each).",
                duration: Some("Overnight + 3h"),
                tips: vec!["Include IgG control for specificity assessment"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Reverse cross-links and purify",
                description: "Elute DNA from beads (1% SDS, 100 mM NaHCO3, 15 min at 65C, 2x). Add NaCl to 200 mM. Incubate overnight at 65C to reverse cross-links. Treat with RNase A (30 min, 37C) and Proteinase K (2h, 45C). Purify DNA.",
                duration: Some("Overnight + 3h"),
                tips: vec!["Process input sample in parallel from this step"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Build sequencing library",
                description: "Use 1-10 ng ChIP DNA. End repair, A-tailing, adapter ligation per library prep kit. PCR amplify (10-15 cycles). Size-select 200-400 bp with AMPure XP beads (0.7x-1.0x double-sided). QC on Bioanalyzer and quantify by qPCR.",
                duration: Some("4-5h"),
                tips: vec![
                    "Minimize PCR cycles to reduce duplicates",
                    "Target 20-30 million reads for TF, 50-100M for broad histone marks",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Sequencing library (200-400 bp inserts)",
            "Bioanalyzer trace showing peak at expected size",
            "Enrichment validated by qPCR at known target loci",
        ],
        references: vec![
            "Park PJ. Nat Rev Genet. 2009;10(10):669-680.",
        ],
    }
}

/// ATAC-seq library preparation protocol.
pub fn atacseq_library_prep() -> Protocol {
    Protocol {
        title: "ATAC-seq Library Preparation",
        slug: "atacseq-library-prep",
        category: ProtocolCategory::WetLab,
        description: "Assay for Transposase-Accessible Chromatin with sequencing. Uses Tn5 transposase to simultaneously fragment and tag open chromatin regions, enabling fast library prep from low cell input.",
        estimated_time: "6-8 hours",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "50,000-100,000 viable cells (fresh, not fixed)",
            "Tn5 transposase (Illumina Tagment DNA Enzyme, or homemade)",
            "Tagmentation buffer (2x TD buffer)",
            "Cell lysis buffer (10 mM Tris pH 7.4, 10 mM NaCl, 3 mM MgCl2, 0.1% NP-40)",
            "Nextera i5/i7 index primers",
            "NEBNext High-Fidelity PCR master mix",
            "MinElute or AMPure XP for cleanup",
            "Bioanalyzer or TapeStation",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Cell lysis",
                description: "Pellet 50,000 cells at 500g, 5 min, 4C. Wash once with cold PBS. Resuspend in 50 uL cold lysis buffer. Incubate on ice 3 min. Centrifuge at 500g, 10 min, 4C. Discard supernatant.",
                duration: Some("20 min"),
                tips: vec![
                    "Cell viability > 90% is critical — check with trypan blue",
                    "Do not over-lyse (avoid >3 min) to preserve nuclear integrity",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Tagmentation",
                description: "Resuspend nuclei in 50 uL tagmentation reaction: 25 uL 2x TD buffer, 2.5 uL Tn5, 22.5 uL water. Incubate at 37C for 30 min with gentle mixing (300 rpm on thermomixer).",
                duration: Some("30 min"),
                tips: vec!["Tagmentation time and Tn5 concentration may need optimization per cell type"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Purify tagmented DNA",
                description: "Immediately purify using MinElute or Zymo columns. Elute in 10 uL elution buffer.",
                duration: Some("15 min"),
                tips: vec!["Do not let tagmented DNA sit — purify immediately after incubation"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "PCR amplification",
                description: "Set up 50 uL PCR: 10 uL tagmented DNA, 25 uL NEBNext HiFi master mix, 2.5 uL each i5/i7 primer (25 uM). Initial extension: 72C 5 min, 98C 30s. Then 5 cycles of 98C 10s / 63C 30s / 72C 1 min. Run qPCR side-reaction to determine additional cycles needed (typically 5-7 more).",
                duration: Some("2h"),
                tips: vec![
                    "Use qPCR to determine optimal cycle number — stop at 1/3 maximum fluorescence",
                    "Total cycles should be 9-12 for 50K cells",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Size selection and QC",
                description: "Double-sided AMPure XP cleanup: 0.5x to remove large fragments, then 1.5x to collect library. Elute in 20 uL. Run Bioanalyzer — expect nucleosomal laddering (sub-nucleosomal ~200 bp, mono ~400 bp, di ~600 bp).",
                duration: Some("1h"),
                tips: vec!["Sub-nucleosomal (<200 bp) fraction should be the dominant peak"],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Sequencing library with nucleosomal laddering pattern on Bioanalyzer",
            "Library concentration > 1 nM",
            "Recommended sequencing: paired-end 50 bp, 50-100M read pairs",
        ],
        references: vec![
            "Buenrostro JD, et al. Nat Methods. 2013;10(12):1213-1218.",
            "Corces MR, et al. Nat Methods. 2017;14(10):959-962.",
        ],
    }
}

/// Hi-C library preparation protocol.
pub fn hic_library_prep() -> Protocol {
    Protocol {
        title: "Hi-C Library Preparation",
        slug: "hic-library-prep",
        category: ProtocolCategory::WetLab,
        description: "Chromosome conformation capture with sequencing to map 3D genome organization. Cross-linked chromatin is digested, proximity-ligated, and sequenced to identify long-range chromatin interactions, TADs, and compartments.",
        estimated_time: "4-5 days",
        difficulty: Difficulty::Advanced,
        requirements: vec![
            "5-10 million cells",
            "1% formaldehyde (methanol-free)",
            "Restriction enzyme (DpnII or MboI for 4-cutter; HindIII for 6-cutter)",
            "Biotin-14-dCTP for fill-in",
            "T4 DNA ligase",
            "Streptavidin beads for biotin pulldown",
            "Library prep kit (NEBNext Ultra II or equivalent)",
            "Proteinase K, RNase A",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Cross-link and lyse",
                description: "Fix cells with 1% formaldehyde for 10 min at room temperature. Quench with 125 mM glycine. Lyse cells in Hi-C lysis buffer (10 mM Tris pH 8.0, 10 mM NaCl, 0.2% Igepal CA-630 + protease inhibitors) on ice for 15 min.",
                duration: Some("45 min"),
                tips: vec![],
                caution: Some("Formaldehyde — use in fume hood"),
            },
            ProtocolStep {
                number: 2,
                title: "Restriction digest",
                description: "Resuspend nuclei in 1x restriction buffer. Add 0.1% SDS, incubate 10 min at 65C. Quench SDS with 1% Triton X-100. Add 400 U restriction enzyme (MboI or DpnII). Digest overnight at 37C with rotation.",
                duration: Some("Overnight"),
                tips: vec![
                    "4-cutter (MboI/DpnII) gives higher resolution than 6-cutter (HindIII)",
                    "Check digestion efficiency by reversing cross-links on a small aliquot and running on gel",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Biotin fill-in and ligation",
                description: "Fill in restriction overhangs with biotin-14-dCTP, dATP, dGTP, dTTP using Klenow (37C, 45 min). Ligate in dilute conditions (large volume, 800 uL) with T4 DNA ligase for 4h at 16C to favor intramolecular (proximity) ligation.",
                duration: Some("5h"),
                tips: vec!["Dilute ligation favors proximity (intramolecular) over random (intermolecular) events"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Reverse cross-links and shear",
                description: "Reverse cross-links overnight at 65C with Proteinase K. Purify DNA. Shear to 300-500 bp using Covaris. Size-select with AMPure XP beads.",
                duration: Some("Overnight + 2h"),
                tips: vec![],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Biotin pulldown and library prep",
                description: "Pull down biotin-labeled ligation junctions with streptavidin beads. Wash stringently. Perform on-bead library prep: end repair, A-tailing, adapter ligation, PCR (6-10 cycles). Size-select 300-500 bp.",
                duration: Some("6h"),
                tips: vec![
                    "Streptavidin pulldown enriches for true Hi-C ligation junctions",
                    "Target 200-400M read pairs for TAD-resolution maps",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Sequencing library (300-500 bp, paired-end)",
            "Expected valid pairs: 40-60% of total reads",
            "Contact matrix showing checkerboard pattern (compartments) and TAD triangles",
        ],
        references: vec![
            "Lieberman-Aiden E, et al. Science. 2009;326(5950):289-293.",
            "Rao SS, et al. Cell. 2014;159(7):1665-1680.",
        ],
    }
}

// ---------------------------------------------------------------------------
// Dry lab protocols (6)
// ---------------------------------------------------------------------------

/// GWAS analysis pipeline protocol.
pub fn gwas_pipeline() -> Protocol {
    Protocol {
        title: "GWAS Analysis Pipeline",
        slug: "gwas-pipeline",
        category: ProtocolCategory::DryLab,
        description: "Genome-wide association study pipeline from genotype data to Manhattan plot. Covers QC, population stratification, association testing, and result visualization.",
        estimated_time: "1-3 days (compute-dependent)",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "Genotype data (PLINK .bed/.bim/.fam or VCF)",
            "Phenotype file (case/control or quantitative trait)",
            "PLINK 1.9 / PLINK 2.0",
            "R with qqman, ggplot2 packages",
            "Reference panel for PCA (e.g., 1000 Genomes)",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Sample QC",
                description: "Filter samples: call rate > 98% (--mind 0.02), heterozygosity within 3 SD of mean, sex check concordance (--check-sex), relatedness (--genome, remove one of each pair with pi-hat > 0.2). Remove ancestry outliers after PCA.",
                duration: Some("1-2h"),
                tips: vec![
                    "Plot heterozygosity vs call rate to identify problematic samples",
                    "Use IBD estimation to detect unexpected duplicates or relatives",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Variant QC",
                description: "Filter variants: call rate > 98% (--geno 0.02), MAF > 1% (--maf 0.01), Hardy-Weinberg equilibrium p > 1e-6 in controls (--hwe 1e-6). Remove strand-ambiguous A/T and C/G SNPs if merging datasets.",
                duration: Some("30 min"),
                tips: vec!["Apply HWE filter in controls only — deviations in cases may be true associations"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Population stratification",
                description: "Merge with reference panel (1000 Genomes). Prune LD (--indep-pairwise 50 5 0.2). Run PCA (--pca 20). Plot PC1 vs PC2 colored by population. Select ethnically matched samples or include PCs as covariates.",
                duration: Some("1-2h"),
                tips: vec!["Include top 10-20 PCs as covariates in association testing to control for stratification"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Association testing",
                description: "For binary traits: logistic regression (--logistic) with PCs and sex as covariates. For quantitative traits: linear regression (--linear). Apply genomic control (lambda_GC should be ~1.0). Genome-wide significance: p < 5e-8.",
                duration: Some("1-4h"),
                tips: vec![
                    "Lambda_GC > 1.1 suggests residual population stratification or polygenicity",
                    "Use SAIGE or BOLT-LMM for biobank-scale (>100K samples) data",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Visualization and interpretation",
                description: "Generate Manhattan plot (genome-wide p-values by position) and QQ plot (observed vs expected -log10 p). Annotate significant loci with nearest genes. Perform conditional analysis to identify independent signals. Check LD with known GWAS hits.",
                duration: Some("1-2h"),
                tips: vec!["Use LocusZoom for detailed regional visualization of significant hits"],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Manhattan plot with genome-wide significant loci highlighted",
            "QQ plot with lambda_GC",
            "Table of significant SNPs (p < 5e-8) with effect sizes, nearest genes, and LD info",
        ],
        references: vec![
            "Marees AT, et al. Int J Methods Psychiatr Res. 2018;27(2):e1608.",
            "Uffelmann E, et al. Nat Rev Methods Primers. 2021;1:59.",
        ],
    }
}

/// Long-read genome assembly protocol.
pub fn longread_genome_assembly() -> Protocol {
    Protocol {
        title: "Long-Read Genome Assembly",
        slug: "longread-genome-assembly",
        category: ProtocolCategory::DryLab,
        description: "De novo genome assembly from PacBio HiFi or Oxford Nanopore long reads. Covers read QC, assembly, polishing, scaffolding, and quality assessment.",
        estimated_time: "2-7 days (genome size-dependent)",
        difficulty: Difficulty::Advanced,
        requirements: vec![
            "Long-read sequencing data (HiFi FASTQ or Nanopore FASTQ/FAST5)",
            "30-50x genome coverage recommended",
            "HiFiasm, Flye, or Hifiasm (assembler)",
            "Short reads for polishing (optional but recommended for ONT)",
            "BUSCO for completeness assessment",
            "QUAST for assembly statistics",
            "Compute: 32+ cores, 64+ GB RAM (genome-size dependent)",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Read QC and filtering",
                description: "Assess read length distribution and quality. For HiFi: filter reads Q < 20 (if not already filtered). For Nanopore: filter by length (> 1 kb) and quality (Q > 10). Target N50 read length > 10 kb. Calculate genome coverage.",
                duration: Some("1-2h"),
                tips: vec![
                    "HiFi reads (Q30+) generally don't need filtering",
                    "For ONT: simplex Q20+ or duplex reads give best assemblies",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Assembly",
                description: "For HiFi: run hifiasm (hifiasm -o asm -t 32 reads.fq.gz). For ONT: run Flye (flye --nano-hq reads.fq -o assembly -t 32 --genome-size 3g). Hifiasm outputs GFA; convert to FASTA with awk. Flye outputs FASTA directly.",
                duration: Some("4h-3 days"),
                tips: vec![
                    "Hifiasm can produce phased assemblies with Hi-C or trio data",
                    "Use --nano-hq for Q20+ ONT reads, --nano-raw for older data",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Polishing (ONT assemblies)",
                description: "For ONT assemblies, polish with Medaka (1 round with ONT reads) then Pilon or NextPolish (2-3 rounds with short reads). For HiFi assemblies, polishing is usually unnecessary.",
                duration: Some("4-12h"),
                tips: vec!["Each polishing round has diminishing returns — 2-3 rounds is typically sufficient"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Scaffolding (optional)",
                description: "If Hi-C data is available, scaffold contigs with SALSA2 or YaHS. Align Hi-C reads, filter for valid pairs, and scaffold. Visualize contact map to validate scaffold correctness.",
                duration: Some("4-8h"),
                tips: vec!["Hi-C scaffolding can achieve chromosome-level assemblies"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Quality assessment",
                description: "Run QUAST for assembly statistics (N50, L50, total length, # contigs). Run BUSCO for gene completeness against lineage database. Run Merqury for k-mer-based QV estimate if short reads available. Expected: BUSCO > 95% complete, QV > 40 for HiFi.",
                duration: Some("1-2h"),
                tips: vec![
                    "Compare assembly size to expected genome size",
                    "Check for contamination with BlobToolKit if sample source is unclear",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Assembly FASTA (contigs or scaffolds)",
            "QUAST report: N50, total length, # contigs",
            "BUSCO completeness score (> 95% target)",
            "Assembly QV (> 40 for HiFi, > 30 for polished ONT)",
        ],
        references: vec![
            "Cheng H, et al. Nat Methods. 2021;18(2):170-175.",
            "Kolmogorov M, et al. Nat Biotechnol. 2019;37(5):540-546.",
        ],
    }
}

/// DIA proteomics search protocol.
pub fn proteomics_dia_search() -> Protocol {
    Protocol {
        title: "DIA Proteomics Analysis",
        slug: "proteomics-dia-search",
        category: ProtocolCategory::DryLab,
        description: "Data-independent acquisition (DIA) mass spectrometry analysis pipeline. Covers spectral library generation, DIA search, protein quantification, and statistical analysis.",
        estimated_time: "1-2 days",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "DIA raw files (.raw, .mzML, or .d)",
            "Protein FASTA database (UniProt, species-specific)",
            "DIA-NN, Spectronaut, or MaxDIA",
            "R or Python for statistical analysis",
            "Optional: DDA runs for spectral library building",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Spectral library preparation",
                description: "Option A (library-free): Use DIA-NN in library-free mode with predicted spectra from a FASTA database. Option B (DDA library): Search DDA files with MSFragger/MaxQuant, build library with Spectronaut or EasyPQP. Library should contain RT, fragment ions, and charge states.",
                duration: Some("2-6h"),
                tips: vec![
                    "Library-free DIA-NN performs comparably to library-based for most experiments",
                    "Include iRT peptides for retention time calibration",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "DIA database search",
                description: "Search DIA files against the spectral library. Key parameters: precursor mass tolerance 10 ppm, fragment mass tolerance 20 ppm, 1% FDR at peptide and protein level. Enable match-between-runs (MBR) for improved quantification coverage.",
                duration: Some("2-8h"),
                tips: vec!["Use target-decoy approach for FDR control"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Protein quantification",
                description: "Extract peptide-level quantities (MS2 fragment ion areas or MS1 precursor areas). Aggregate to protein level using MaxLFQ algorithm or Top3 method. Apply median normalization across runs.",
                duration: Some("30 min"),
                tips: vec!["Filter: require >= 2 unique peptides per protein identification"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Quality control",
                description: "Check: protein ID counts across runs (CV < 10%), total ion chromatograms, RT alignment quality, missing value patterns. Plot PCA of protein quantities colored by condition/batch. Identify and flag outlier runs.",
                duration: Some("1h"),
                tips: vec!["Missing values in DIA are often MNAR (missing not at random) — imputation strategy matters"],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Statistical analysis",
                description: "Impute missing values (MinProb, KNN, or QRILC). Log2-transform quantities. Run limma or MSstats for differential abundance testing. Apply BH correction. Threshold: adjusted p < 0.05, |log2FC| > 1. Generate volcano plot.",
                duration: Some("1-2h"),
                tips: vec![
                    "MSstats handles missing values and experimental design explicitly",
                    "For TMT DIA, use MSstatsTMT",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Protein quantification matrix (proteins x samples)",
            "List of differentially abundant proteins (adj. p < 0.05, |FC| > 2)",
            "Volcano plot and PCA plot",
            "Pathway enrichment results",
        ],
        references: vec![
            "Demichev V, et al. Nat Methods. 2020;17(1):41-44.",
        ],
    }
}

/// Spatial transcriptomics analysis protocol.
pub fn spatial_transcriptomics_analysis() -> Protocol {
    Protocol {
        title: "Spatial Transcriptomics Analysis",
        slug: "spatial-transcriptomics-analysis",
        category: ProtocolCategory::DryLab,
        description: "Analysis pipeline for spatially-resolved gene expression data (10x Visium, MERFISH, Slide-seq). Covers preprocessing, clustering, spatial domain detection, cell-cell communication, and deconvolution.",
        estimated_time: "1-2 days",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "Spatial transcriptomics data (Space Ranger output, MERFISH coordinates, or Slide-seq DGE)",
            "Tissue image (H&E for Visium)",
            "Reference single-cell dataset (for deconvolution)",
            "Python (scanpy, squidpy) or R (Seurat, STUtility)",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Data loading and QC",
                description: "Load spatial data with coordinates and expression matrix. For Visium: filter spots by tissue coverage and total UMI counts (> 200 genes, < 20% mitochondrial). For MERFISH: filter cells by minimum transcript count and estimated FPR < 0.2.",
                duration: Some("30 min"),
                tips: vec![
                    "Overlay QC metrics on tissue image to check spatial patterns",
                    "Mitochondrial genes may show spatial bias in tissue sections",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Normalization and feature selection",
                description: "Normalize per-spot/cell library size (SCTransform for Visium, or total-count normalization). Identify highly variable genes (HVGs). For MERFISH, gene panel is pre-selected so HVG step may be skipped.",
                duration: Some("30 min"),
                tips: vec!["SCTransform regularizes variance and often works better than log-normalization for spatial data"],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Clustering and spatial domains",
                description: "Run PCA and build spatial neighbor graph (KNN + spatial coordinates). Cluster with Leiden algorithm. Alternatively, detect spatial domains with BayesSpace or SpaGCN that integrate expression and spatial information. Visualize clusters on tissue image.",
                duration: Some("1-2h"),
                tips: vec![
                    "Spatial-aware clustering (BayesSpace, SpaGCN) gives more spatially coherent domains than standard Leiden",
                    "Compare domain boundaries with histological annotations",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Spatially variable genes and communication",
                description: "Identify spatially variable genes (SVGs) using Moran's I or SpatialDE. Run ligand-receptor analysis with CellChat or COMMOT to infer cell-cell communication. Visualize communication patterns on tissue.",
                duration: Some("2-4h"),
                tips: vec![
                    "Moran's I > 0.3 with p < 0.01 indicates strong spatial autocorrelation",
                    "L-R analysis requires cell type annotations — use deconvolution results",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Deconvolution and integration",
                description: "Deconvolve multi-cellular spots (Visium) into cell type proportions using RCTD, Cell2location, or NNLS with a scRNA-seq reference. Validate by checking consistency with known tissue architecture. Integrate with other spatial datasets or modalities if available.",
                duration: Some("2-4h"),
                tips: vec![
                    "Reference quality matters — use well-annotated scRNA-seq from the same tissue",
                    "RCTD works well for Visium; Cell2location for larger panels",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Spatial domain map overlaid on tissue image",
            "Spatially variable gene list ranked by Moran's I",
            "Cell type proportion maps per spot",
            "Ligand-receptor interaction network and spatial communication patterns",
        ],
        references: vec![
            "Palla G, et al. Nat Methods. 2022;19(2):171-178.",
            "Cable DM, et al. Nat Biotechnol. 2022;40(4):517-526.",
        ],
    }
}

/// Hi-C data analysis protocol.
pub fn hic_analysis() -> Protocol {
    Protocol {
        title: "Hi-C Data Analysis Pipeline",
        slug: "hic-analysis",
        category: ProtocolCategory::DryLab,
        description: "Computational analysis of Hi-C data for 3D genome organization. Covers read alignment, contact matrix generation, normalization, TAD calling, A/B compartment identification, and loop detection.",
        estimated_time: "1-3 days",
        difficulty: Difficulty::Advanced,
        requirements: vec![
            "Hi-C paired-end FASTQ files (200-400M read pairs for TAD resolution)",
            "Reference genome (BWA index)",
            "BWA-MEM2 or Chromap (aligner)",
            "pairtools / HiC-Pro (contact processing)",
            "cooler (contact matrix storage)",
            "cooltools / FAN-C (analysis)",
            "Compute: 16+ cores, 32+ GB RAM",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Alignment and contact parsing",
                description: "Align each read end independently to the reference genome (BWA-MEM2 -SP5M). Parse contacts with pairtools: assign strands, classify pair types (UU, UR, RU, etc.), remove duplicates. Filter for valid cis contacts. Expected valid pairs: 40-60% of total.",
                duration: Some("4-12h"),
                tips: vec![
                    "Align R1 and R2 separately then pair — chimeric reads span ligation junctions",
                    "Use pairtools stats to monitor library quality metrics",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "Contact matrix generation",
                description: "Bin valid pairs into contact matrices at multiple resolutions (1 kb, 5 kb, 10 kb, 25 kb, 100 kb, 500 kb, 1 Mb) using cooler. Store as .mcool file. Verify sufficient coverage at target resolution (> 80% non-zero bins).",
                duration: Some("1-2h"),
                tips: vec![
                    "1 kb resolution needs ~1 billion valid pairs for human genome",
                    "10-25 kb is typical for TAD analysis; 100 kb-1 Mb for compartments",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Normalization",
                description: "Apply iterative correction (ICE) or Knight-Ruiz (KR) balancing to correct for GC content, mappability, and restriction site biases. Use cooltools balance. Check that balanced matrix diagonal is uniform.",
                duration: Some("30-60 min"),
                tips: vec!["KR balancing is faster; ICE is more robust to sparse matrices"],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "TADs and compartments",
                description: "Call TADs: compute insulation scores at 10-25 kb resolution (cooltools insulation, window 500 kb). Boundaries are local minima. Call A/B compartments: compute observed/expected matrix at 100 kb, correlate eigenvector with GC content to orient. A = active/gene-rich, B = repressed.",
                duration: Some("1-2h"),
                tips: vec![
                    "Insulation window size affects TAD granularity — try 250 kb, 500 kb, 1 Mb",
                    "Compartment switching between conditions indicates epigenetic reprogramming",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Loop detection and visualization",
                description: "Call loops with HiCCUPS (Juicer) or cooltools dots at 5-10 kb resolution. Loops appear as focal enrichments off-diagonal. Require >= 2-fold enrichment over donut background, FDR < 0.1. Visualize with HiGlass or cooler show for interactive exploration.",
                duration: Some("2-4h"),
                tips: vec![
                    "Loop calling requires high resolution (5-10 kb) — ensure sufficient read depth",
                    "CTCF ChIP-seq at loop anchors validates Hi-C loops",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            ".mcool file with multi-resolution contact matrices",
            "TAD boundary list with insulation scores",
            "A/B compartment assignments per genomic bin",
            "Loop list with coordinates and enrichment scores",
            "Hi-C contact map visualizations",
        ],
        references: vec![
            "Abdennur N, Mirny LA. Bioinformatics. 2020;36(11):3525-3528.",
            "Open2C, et al. bioRxiv. 2022. doi:10.1101/2022.10.31.514564.",
        ],
    }
}

/// AlphaFold structure prediction protocol.
pub fn alphafold_structure_prediction() -> Protocol {
    Protocol {
        title: "AlphaFold Structure Prediction",
        slug: "alphafold-structure-prediction",
        category: ProtocolCategory::DryLab,
        description: "Protein structure prediction using AlphaFold2 or AlphaFold3. Covers sequence preparation, MSA generation, structure prediction, confidence assessment, and model validation.",
        estimated_time: "1-24 hours (protein size-dependent)",
        difficulty: Difficulty::Intermediate,
        requirements: vec![
            "Target protein sequence(s) in FASTA format",
            "AlphaFold2 installation with model weights (or ColabFold for cloud)",
            "Sequence databases (UniRef90, MGnify, PDB70, BFD/UniClust30)",
            "GPU with >= 16 GB VRAM (A100 recommended for large proteins)",
            "PyMOL or ChimeraX for visualization",
        ],
        steps: vec![
            ProtocolStep {
                number: 1,
                title: "Sequence preparation",
                description: "Prepare target sequence in FASTA format. For multimers, provide all chain sequences. Remove signal peptides, tags, and disordered regions if not needed. Check sequence against UniProt for known structures or close homologs in PDB.",
                duration: Some("15 min"),
                tips: vec![
                    "For complexes, AlphaFold-Multimer predicts inter-chain contacts",
                    "Sequences > 2500 residues may need to be split into domains",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 2,
                title: "MSA generation",
                description: "Run jackhmmer against UniRef90 and MGnify. Run HHblits against BFD/UniClust30. Search PDB70 for structural templates. ColabFold uses MMseqs2 for faster MSA generation. Assess MSA depth — > 100 effective sequences recommended.",
                duration: Some("30 min - 4h"),
                tips: vec![
                    "Shallow MSA (< 30 sequences) reduces prediction confidence",
                    "ColabFold's MMseqs2 is ~100x faster than jackhmmer with similar accuracy",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 3,
                title: "Structure prediction",
                description: "Run AlphaFold2 with 5 model presets (model_1 through model_5). Each produces a predicted structure with per-residue confidence (pLDDT) and predicted aligned error (PAE). For multimers, use multimer mode. GPU time: minutes (small proteins) to hours (large complexes).",
                duration: Some("30 min - 12h"),
                tips: vec![
                    "Enable templates for better accuracy on proteins with PDB homologs",
                    "Recycle 3-20 times (default 3) for difficult targets",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 4,
                title: "Confidence assessment",
                description: "Evaluate predictions using: pLDDT (per-residue, > 90 very high, 70-90 confident, 50-70 low, < 50 disordered). PAE for domain-domain orientation confidence. ipTM for multimer interface quality (> 0.7 is good). Rank models by ipTM*0.8 + pTM*0.2 (multimer) or pLDDT mean (monomer).",
                duration: Some("30 min"),
                tips: vec![
                    "Low pLDDT regions are often intrinsically disordered — don't interpret their structure",
                    "High PAE between domains means their relative orientation is uncertain",
                ],
                caution: None,
            },
            ProtocolStep {
                number: 5,
                title: "Validation and analysis",
                description: "Validate with MolProbity (Ramachandran, clashscore, rotamer outliers). Compare to experimental structures if available (RMSD, TM-score). Analyze predicted contacts, binding sites, and functional regions. Visualize in PyMOL/ChimeraX colored by pLDDT.",
                duration: Some("1h"),
                tips: vec![
                    "TM-score > 0.5 indicates same fold; > 0.7 indicates close structural match",
                    "Deposit models in ModelArchive for community access",
                ],
                caution: None,
            },
        ],
        expected_outputs: vec![
            "Ranked predicted structures in PDB format (5 models)",
            "Per-residue pLDDT confidence plot",
            "PAE heatmap for domain orientation confidence",
            "MolProbity validation report",
        ],
        references: vec![
            "Jumper J, et al. Nature. 2021;596(7873):583-589.",
            "Mirdita M, et al. Nat Methods. 2022;19(6):679-682.",
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_protocols_count() {
        let all = all_protocols();
        assert_eq!(all.len(), 16);
    }

    #[test]
    fn test_wet_lab_count() {
        let wet = wet_lab_protocols();
        assert_eq!(wet.len(), 10);
        assert!(wet.iter().all(|p| p.category == ProtocolCategory::WetLab));
    }

    #[test]
    fn test_dry_lab_count() {
        let dry = dry_lab_protocols();
        assert_eq!(dry.len(), 6);
        assert!(dry.iter().all(|p| p.category == ProtocolCategory::DryLab));
    }

    #[test]
    fn test_unique_slugs() {
        let all = all_protocols();
        let mut slugs: Vec<&str> = all.iter().map(|p| p.slug).collect();
        let count = slugs.len();
        slugs.sort();
        slugs.dedup();
        assert_eq!(slugs.len(), count, "duplicate slugs found");
    }

    #[test]
    fn test_all_have_steps() {
        let all = all_protocols();
        for p in &all {
            assert!(!p.steps.is_empty(), "{} has no steps", p.title);
            assert!(!p.requirements.is_empty(), "{} has no requirements", p.title);
            assert!(!p.expected_outputs.is_empty(), "{} has no expected outputs", p.title);
        }
    }

    #[test]
    fn test_step_numbering() {
        let all = all_protocols();
        for p in &all {
            for (i, step) in p.steps.iter().enumerate() {
                assert_eq!(step.number as usize, i + 1, "{}: step {} misnumbered", p.title, i + 1);
            }
        }
    }

    #[test]
    fn test_markdown_rendering() {
        let p = elisa();
        let md = p.to_markdown();
        assert!(md.contains("# ELISA"));
        assert!(md.contains("## Requirements"));
        assert!(md.contains("## Steps"));
        assert!(md.contains("### Step 1:"));
        assert!(md.contains("## Expected Outputs"));
        assert!(md.contains("## References"));
    }

    #[test]
    fn test_markdown_contains_tips() {
        let p = qpcr();
        let md = p.to_markdown();
        assert!(md.contains("**Tip:**"));
    }

    #[test]
    fn test_markdown_contains_caution() {
        let p = elisa();
        let md = p.to_markdown();
        assert!(md.contains("**Caution:**"));
    }

    #[test]
    fn test_protocol_difficulty_levels() {
        let all = all_protocols();
        let beginner: Vec<_> = all.iter().filter(|p| p.difficulty == Difficulty::Beginner).collect();
        let intermediate: Vec<_> = all.iter().filter(|p| p.difficulty == Difficulty::Intermediate).collect();
        let advanced: Vec<_> = all.iter().filter(|p| p.difficulty == Difficulty::Advanced).collect();
        assert!(!beginner.is_empty());
        assert!(!intermediate.is_empty());
        assert!(!advanced.is_empty());
    }

    #[test]
    fn test_protocol_references() {
        let all = all_protocols();
        for p in &all {
            assert!(!p.references.is_empty(), "{} has no references", p.title);
        }
    }

    #[test]
    fn test_specific_protocols_exist() {
        let all = all_protocols();
        let slugs: Vec<&str> = all.iter().map(|p| p.slug).collect();
        assert!(slugs.contains(&"elisa"));
        assert!(slugs.contains(&"qpcr-rt-qpcr"));
        assert!(slugs.contains(&"chipseq-library-prep"));
        assert!(slugs.contains(&"atacseq-library-prep"));
        assert!(slugs.contains(&"hic-library-prep"));
        assert!(slugs.contains(&"gwas-pipeline"));
        assert!(slugs.contains(&"longread-genome-assembly"));
        assert!(slugs.contains(&"alphafold-structure-prediction"));
    }
}
