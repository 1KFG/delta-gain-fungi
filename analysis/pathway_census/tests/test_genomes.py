from pathlib import Path

import pytest

from pathway_census.genomes import load_genome_table, parse_genus

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def df():
    # fixture dir needs species.parquet etc.; build a tables dir of symlinks
    import tempfile, os
    d = Path(tempfile.mkdtemp())
    for real, link in [("species_mini", "species"), ("busco_mini", "busco_genome"), ("asm_stats_mini", "asm_stats")]:
        os.symlink(DATA / f"{real}.parquet", d / f"{link}.parquet")
    return load_genome_table(d, DATA / "prior_status_mini.tsv").set_index("locustag")


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


def test_unresolvable_genus_is_none(df):
    assert df.loc["AAAAAAA6", "genus"] is None
    assert df.loc["AAAAAAA6", "genus_source"] is None


def test_status_mapping(df):
    assert df.loc["AAAAAAA1", "prior_status"] == "baseline_ok"
    assert df.loc["AAAAAAA2", "prior_status"] == "no_prior_annotation_present"
    assert df.loc["AAAAAAA3", "prior_status"] == "no_prior_annotation_not_downloaded"
    assert df.loc["AAAAAAA4", "prior_status"] == "annotated_folder_missing_protein"
    assert df.loc["AAAAAAA5", "prior_status"] == "annotated_NOT_downloaded"
    assert df.loc["AAAAAAA1", "annotation_provider"] == "NCBI"


def test_not_downloaded_has_downloaded_false(df):
    assert not df.loc["AAAAAAA3", "downloaded"]
    assert not df.loc["AAAAAAA5", "downloaded"]
    assert df.loc["AAAAAAA1", "downloaded"]
    assert df.loc["AAAAAAA4", "downloaded"]


def test_no_protein_row_has_has_proteins_false(df):
    assert not df.loc["AAAAAAA4", "has_proteins"]
    assert not df.loc["AAAAAAA2", "has_proteins"]
    assert df.loc["AAAAAAA1", "has_proteins"]


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
