import os
from pathlib import Path

import pytest

from pathway_census.genomes import load_genome_table, parse_genus

DATA = Path(__file__).parent / "data"


# seqdb layout per locustag: (has proteins file, has scaffolds file); absent key = not in genomes.tsv
SEQDB = {
    "AAAAAAA1": (True, True),
    "AAAAAAA2": (True, True),
    "AAAAAAA3": (True, True),
    "AAAAAAA4": (True, False),
    "AAAAAAA5": (False, False),
    "AAAAAAA7": (True, True),
}


@pytest.fixture(scope="module")
def df(tmp_path_factory):
    root = tmp_path_factory.mktemp("bfd")
    tables = root / "tables"
    tables.mkdir()
    for real, link in [("species_mini", "species"), ("busco_mini", "busco_genome"), ("asm_stats_mini", "asm_stats")]:
        os.symlink(DATA / f"{real}.parquet", tables / f"{link}.parquet")
    (root / "seqdb").mkdir()
    lines = []
    for lt, (prot, scaf) in SEQDB.items():
        d = root / "ga" / lt / "predict_results"
        d.mkdir(parents=True)
        pp = d / f"{lt}.proteins.fa.gz"
        if prot:
            pp.write_text("x")
        if scaf:
            (d / f"{lt}.scaffolds.fa.gz").write_text("x")
        lines.append(f"{lt}\t{lt}\t{pp}\t{d / (lt + '.cds-transcripts.fa.gz')}")
    (root / "seqdb" / "genomes.tsv").write_text("\n".join(lines) + "\n")
    return load_genome_table(tables, DATA / "prior_status_mini.tsv").set_index("locustag")


def test_empty_genus_is_parsed_from_species(df):
    for k in ("AAAAAAA2", "AAAAAAA3"):  # None and empty string
        assert df.loc[k, "genus"] == "Antarctolithica"
        assert df.loc[k, "genus_source"] == "parsed_species"
    assert df.loc["AAAAAAA1", "genus_source"] == "genus_field"


def test_parse_genus_rules():
    assert parse_genus("Ascomycota sp. Antarctolithica obscura") == "Antarctolithica"
    assert parse_genus("Candida mesenterica") == "Candida"
    assert parse_genus("Candida sp.") == "Candida"
    assert parse_genus("Saccharomycetaceae sp.") is None
    assert parse_genus("Ascomycota sp.") is None
    assert parse_genus("fungal sp.") is None
    assert parse_genus("aff. Porphyrellus sp. UT-M0000030") == "Porphyrellus"
    assert parse_genus("Microsporidia sp. MB") is None  # phylum, strain code
    assert parse_genus("Tulasnella sp. JGI-2019a") == "Tulasnella"


def test_unresolvable_genus_gets_label(df):
    assert df.loc["AAAAAAA6", "genus"] == "Saccharomycetaceae sp."
    assert df.loc["AAAAAAA6", "genus_source"] == "unresolved"
    assert df.loc["AAAAAAA5", "genus"] == "Mucor"  # GENUS field present
    assert df.loc["AAAAAAA8", "genus_source"] == "genus_field"


def test_genus_never_null(df):
    assert df["genus"].notna().all()
    assert (df["genus"].str.len() > 0).all()
    assert set(df["genus_source"]) <= {"genus_field", "parsed_species", "unresolved"}


def test_status_mapping(df):
    assert df.loc["AAAAAAA1", "prior_status"] == "baseline_ok"
    assert df.loc["AAAAAAA2", "prior_status"] == "no_prior_annotation_present"
    assert df.loc["AAAAAAA3", "prior_status"] == "no_prior_annotation_not_downloaded"
    assert df.loc["AAAAAAA4", "prior_status"] == "annotated_folder_missing_protein"
    assert df.loc["AAAAAAA5", "prior_status"] == "annotated_NOT_downloaded"
    assert df.loc["AAAAAAA1", "annotation_provider"] == "NCBI"


def test_downloaded_means_bfd_genome_fasta_exists(df):
    # prior state "not downloaded" must not force downloaded=False
    assert df.loc["AAAAAAA3", "prior_folder_state"] == "not_downloaded"
    assert df.loc["AAAAAAA3", "downloaded"]
    assert df.loc["AAAAAAA1", "downloaded"]
    assert not df.loc["AAAAAAA4", "downloaded"]  # proteins but no scaffolds file
    assert not df.loc["AAAAAAA5", "downloaded"]
    assert not df.loc["AAAAAAA6", "downloaded"]  # not in seqdb


def test_has_proteins_means_bfd_proteins_exist(df):
    # no prior annotation, but BFD re-annotated: has_proteins True, deposited False
    assert df.loc["AAAAAAA2", "has_proteins"]
    assert not df.loc["AAAAAAA2", "has_deposited_proteins"]
    assert df.loc["AAAAAAA1", "has_proteins"] and df.loc["AAAAAAA1", "has_deposited_proteins"]
    assert not df.loc["AAAAAAA5", "has_proteins"]  # listed in seqdb, file missing
    assert not df.loc["AAAAAAA8", "has_proteins"]  # not in seqdb
    assert df.loc["AAAAAAA4", "has_proteins"] and not df.loc["AAAAAAA4", "has_deposited_proteins"]


def test_row_count_equals_species_rows(df):
    assert len(df) == 8
    # row 8 has no prior-status, busco or asm_stats row: kept, with nulls
    assert df.loc["AAAAAAA8", "prior_status"] == "unknown"
    assert df.loc["AAAAAAA8", "busco_complete"] != df.loc["AAAAAAA8", "busco_complete"]  # NaN


def test_columns_and_values(df):
    assert df.loc["AAAAAAA1", "busco_complete"] == 98.5
    assert df.loc["AAAAAAA1", "n50"] == 100000
    assert df.loc["AAAAAAA1", "class_"] == "Dothideomycetes"
    assert df.loc["AAAAAAA1", "phylum"] == "Ascomycota"


def test_asmid_and_folder_are_carried(df):
    # asmid from species.parquet; folder from the prior-annotation TSV (it can differ from ASMID)
    assert df.loc["AAAAAAA1", "asmid"] == "GCA_000000001.1_A"
    assert df.loc["AAAAAAA1", "folder"] == "GCA_000000001.1_A"
    assert df.loc["AAAAAAA7", "asmid"] == "GCA_000000007.1_G"
    assert df.loc["AAAAAAA7", "folder"] == "GCA_000000007.1_G_alt"
    assert df.loc["AAAAAAA8", "folder"] != df.loc["AAAAAAA8", "folder"] or df.loc["AAAAAAA8", "folder"] is None  # no TSV row
    assert list(df.reset_index().columns) == [
        "locustag", "asmid", "folder", "species", "phylum", "class_", "genus", "genus_source", "prior_status",
        "prior_folder_state", "annotation_provider", "busco_complete", "n50", "downloaded", "has_proteins",
        "has_deposited_proteins"]
