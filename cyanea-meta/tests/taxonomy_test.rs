/// Integration tests for taxonomy module.

use cyanea_meta::taxonomy::{TaxonomyDB, TaxonRank};

#[test]
fn build_taxonomy_tree() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();
    db.add_node(202, 2, TaxonRank::Phylum, "Firmicutes").unwrap();
    db.add_node(2011, 201, TaxonRank::Species, "E. coli").unwrap();
    db.add_node(2012, 201, TaxonRank::Species, "Salmonella").unwrap();
    db.add_node(2021, 202, TaxonRank::Species, "B. subtilis").unwrap();

    assert_eq!(db.len(), 7);
}

#[test]
fn get_lineage() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();
    db.add_node(2011, 201, TaxonRank::Species, "E. coli").unwrap();

    let lineage = db.get_lineage(2011).unwrap();
    assert_eq!(lineage, vec![2011, 201, 2, 1]);
}

#[test]
fn lca_sibling_species() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();
    db.add_node(2011, 201, TaxonRank::Species, "E. coli").unwrap();
    db.add_node(2012, 201, TaxonRank::Species, "Salmonella").unwrap();

    // LCA(E. coli, Salmonella) should be Proteobacteria
    let lca = db.lca(&[2011, 2012]).unwrap();
    assert_eq!(lca, 201);
}

#[test]
fn lca_cross_phylum() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();
    db.add_node(202, 2, TaxonRank::Phylum, "Firmicutes").unwrap();
    db.add_node(2011, 201, TaxonRank::Species, "E. coli").unwrap();
    db.add_node(2021, 202, TaxonRank::Species, "B. subtilis").unwrap();

    // LCA(E. coli, B. subtilis) should be Bacteria
    let lca = db.lca(&[2011, 2021]).unwrap();
    assert_eq!(lca, 2);
}

#[test]
fn lca_single_taxon() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();

    let lca = db.lca(&[201]).unwrap();
    assert_eq!(lca, 201);
}

#[test]
fn classify_sequence_kmer() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();

    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp, long enough for k=31
    db.add_reference(seq, 201).unwrap();

    let result = db.classify_sequence(seq).unwrap();
    assert_eq!(result, Some(201));
}

#[test]
fn classify_no_match() {
    let db = TaxonomyDB::new(1);
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let result = db.classify_sequence(seq).unwrap();
    assert_eq!(result, None);
}

#[test]
fn get_rank() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Phylum, "Bacteria").unwrap();
    db.add_node(201, 2, TaxonRank::Species, "E. coli").unwrap();

    let rank = db.get_rank(201).unwrap();
    assert_eq!(rank, TaxonRank::Species);
}

#[test]
fn get_name() {
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();

    let name = db.get_name(2).unwrap();
    assert_eq!(name, "Bacteria");
}

#[test]
fn taxon_not_found_error() {
    let db = TaxonomyDB::new(1);
    assert!(db.get_rank(9999).is_err());
    assert!(db.get_name(9999).is_err());
    assert!(db.get_lineage(9999).is_err());
}

#[test]
fn invalid_parent_error() {
    let mut db = TaxonomyDB::new(1);
    let result = db.add_node(2, 9999, TaxonRank::Phylum, "Bacteria");
    assert!(result.is_err());
}
