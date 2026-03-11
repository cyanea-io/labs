//! Metagenomics demo datasets: OTU tables, taxonomy profiles.

/// Demo OTU table from a gut microbiome study.
///
/// 8 taxa × 6 samples (3 healthy, 3 IBD).
pub fn demo_otu_table() -> DemoOtuTable {
    DemoOtuTable {
        taxa: vec![
            "Bacteroides_fragilis",
            "Faecalibacterium_prausnitzii",
            "Escherichia_coli",
            "Lactobacillus_acidophilus",
            "Ruminococcus_bromii",
            "Akkermansia_muciniphila",
            "Clostridium_difficile",
            "Bifidobacterium_longum",
        ],
        samples: vec![
            "Healthy_1", "Healthy_2", "Healthy_3",
            "IBD_1", "IBD_2", "IBD_3",
        ],
        counts: vec![
            // Bacteroides_fragilis — common in both
            vec![800, 780, 820, 250, 220, 260],
            // Faecalibacterium — depleted in IBD
            vec![750, 700, 780, 80, 60, 90],
            // E. coli — elevated in IBD (dominant in IBD)
            vec![150, 120, 180, 1800, 2020, 1950],
            // Lactobacillus — moderate healthy, low IBD
            vec![600, 650, 580, 50, 40, 60],
            // Ruminococcus — depleted in IBD
            vec![680, 750, 620, 30, 25, 35],
            // Akkermansia — moderate healthy, low IBD
            vec![500, 480, 520, 20, 15, 25],
            // C. difficile — elevated in IBD
            vec![20, 15, 25, 900, 1050, 980],
            // Bifidobacterium — moderate
            vec![550, 620, 480, 60, 50, 70],
        ],
        metadata: vec![
            ("Healthy_1", "healthy"),
            ("Healthy_2", "healthy"),
            ("Healthy_3", "healthy"),
            ("IBD_1", "ibd"),
            ("IBD_2", "ibd"),
            ("IBD_3", "ibd"),
        ],
    }
}

/// A demo OTU table.
#[derive(Debug, Clone)]
pub struct DemoOtuTable {
    pub taxa: Vec<&'static str>,
    pub samples: Vec<&'static str>,
    /// Counts matrix (taxa × samples).
    pub counts: Vec<Vec<u64>>,
    /// Sample metadata: (sample_name, group).
    pub metadata: Vec<(&'static str, &'static str)>,
}

impl DemoOtuTable {
    /// Total counts per sample.
    pub fn sample_totals(&self) -> Vec<u64> {
        (0..self.samples.len())
            .map(|s| self.counts.iter().map(|row| row[s]).sum())
            .collect()
    }

    /// Relative abundance matrix.
    pub fn relative_abundance(&self) -> Vec<Vec<f64>> {
        let totals = self.sample_totals();
        self.counts
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(s, &count)| count as f64 / totals[s] as f64)
                    .collect()
            })
            .collect()
    }

    /// Shannon diversity per sample.
    pub fn shannon_diversity(&self) -> Vec<f64> {
        let rel = self.relative_abundance();
        let n_samples = self.samples.len();
        (0..n_samples)
            .map(|s| {
                let mut h = 0.0f64;
                for row in &rel {
                    let p = row[s];
                    if p > 0.0 {
                        h -= p * p.ln();
                    }
                }
                h
            })
            .collect()
    }
}

/// Demo taxonomy lineage for the OTU table taxa.
pub fn demo_taxonomy() -> Vec<DemoTaxonomy> {
    vec![
        DemoTaxonomy { species: "Bacteroides_fragilis", genus: "Bacteroides", family: "Bacteroidaceae", order: "Bacteroidales", class: "Bacteroidia", phylum: "Bacteroidetes" },
        DemoTaxonomy { species: "Faecalibacterium_prausnitzii", genus: "Faecalibacterium", family: "Ruminococcaceae", order: "Clostridiales", class: "Clostridia", phylum: "Firmicutes" },
        DemoTaxonomy { species: "Escherichia_coli", genus: "Escherichia", family: "Enterobacteriaceae", order: "Enterobacterales", class: "Gammaproteobacteria", phylum: "Proteobacteria" },
        DemoTaxonomy { species: "Lactobacillus_acidophilus", genus: "Lactobacillus", family: "Lactobacillaceae", order: "Lactobacillales", class: "Bacilli", phylum: "Firmicutes" },
        DemoTaxonomy { species: "Ruminococcus_bromii", genus: "Ruminococcus", family: "Ruminococcaceae", order: "Clostridiales", class: "Clostridia", phylum: "Firmicutes" },
        DemoTaxonomy { species: "Akkermansia_muciniphila", genus: "Akkermansia", family: "Akkermansiaceae", order: "Verrucomicrobiales", class: "Verrucomicrobiae", phylum: "Verrucomicrobia" },
        DemoTaxonomy { species: "Clostridium_difficile", genus: "Clostridioides", family: "Peptostreptococcaceae", order: "Clostridiales", class: "Clostridia", phylum: "Firmicutes" },
        DemoTaxonomy { species: "Bifidobacterium_longum", genus: "Bifidobacterium", family: "Bifidobacteriaceae", order: "Bifidobacteriales", class: "Actinobacteria", phylum: "Actinobacteria" },
    ]
}

/// Taxonomy lineage for a species.
#[derive(Debug, Clone)]
pub struct DemoTaxonomy {
    pub species: &'static str,
    pub genus: &'static str,
    pub family: &'static str,
    pub order: &'static str,
    pub class: &'static str,
    pub phylum: &'static str,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_otu_table() {
        let otu = demo_otu_table();
        assert_eq!(otu.taxa.len(), 8);
        assert_eq!(otu.samples.len(), 6);
        assert!(otu.counts.iter().all(|row| row.len() == 6));
    }

    #[test]
    fn test_sample_totals() {
        let otu = demo_otu_table();
        let totals = otu.sample_totals();
        assert_eq!(totals.len(), 6);
        assert!(totals.iter().all(|&t| t > 1000));
    }

    #[test]
    fn test_relative_abundance() {
        let otu = demo_otu_table();
        let rel = otu.relative_abundance();
        // Each column should sum to ~1.0
        for s in 0..6 {
            let col_sum: f64 = rel.iter().map(|row| row[s]).sum();
            assert!((col_sum - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_shannon_diversity() {
        let otu = demo_otu_table();
        let div = otu.shannon_diversity();
        assert_eq!(div.len(), 6);
        // Healthy samples should have higher diversity than IBD
        let healthy_mean: f64 = div[0..3].iter().sum::<f64>() / 3.0;
        let ibd_mean: f64 = div[3..6].iter().sum::<f64>() / 3.0;
        assert!(healthy_mean > ibd_mean);
    }

    #[test]
    fn test_taxonomy() {
        let taxa = demo_taxonomy();
        assert_eq!(taxa.len(), 8);
        // Check Firmicutes are correctly assigned
        let firmicutes: Vec<&str> = taxa.iter()
            .filter(|t| t.phylum == "Firmicutes")
            .map(|t| t.species)
            .collect();
        assert!(firmicutes.contains(&"Faecalibacterium_prausnitzii"));
    }
}
