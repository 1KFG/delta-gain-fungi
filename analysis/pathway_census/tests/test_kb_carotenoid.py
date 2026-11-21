"""Tests for the carotenoid flat knowledge tables (stage 1, task S3)."""
from pathlib import Path

import pytest

from pathway_census.defs import load_pathway
from pathway_census.kb import (EVIDENCE_TIERS, formula_mass, load_kb, path_possible,
                               validate_kb)

ROOT = Path(__file__).resolve().parent.parent
KB = ROOT / "kb" / "carotenoid"
YAML = ROOT / "pathways" / "carotenoid.yaml"
# Steps that exist only in the tables: hypothesis-level (T4) steps with no reference protein.
T4_ONLY_STEPS = {"torularhodin_formation"}
ASTAXANTHIN = "CHEBI:40968"


@pytest.fixture(scope="module")
def kb():
    return load_kb(KB)


def _parts(kb, role):
    out = {}
    for _, r in kb["reaction_participant"].iterrows():
        if r["role"] == role:
            out.setdefault(r["reaction_id"], set()).add(r["compound"])
    return out


def test_tables_validate(kb):
    assert validate_kb(kb) == []


def test_tables_not_empty(kb):
    for name, df in kb.items():
        assert len(df) > 0, name


def test_every_t1_t2_link_has_pmid_and_quote(kb):
    sr = kb["step_reaction"]
    ev = sr[sr["tier"].isin(EVIDENCE_TIERS)]
    assert len(ev) > 0
    for _, r in ev.iterrows():
        assert r["pmid"].strip().isdigit(), (r["step_id"], r["reaction_id"])
        assert r["quote"].strip(), (r["step_id"], r["reaction_id"])
        assert r["assay"].strip(), (r["step_id"], r["reaction_id"])


def test_one_quoted_sentence_per_paper(kb):
    sr = kb["step_reaction"]
    sr = sr[sr["quote"].str.strip() != ""]
    for pmid, g in sr.groupby("pmid"):
        assert g["quote"].nunique() == 1, pmid


def test_every_link_names_organism_gene_and_accession_of_a_yaml_ref(kb):
    refs = {}
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            refs[(st.id, r.accession)] = r
    for _, r in kb["step_reaction"].iterrows():
        if r["step_id"] in T4_ONLY_STEPS:
            assert r["tier"] == "T4"
            continue
        ref = refs[(r["step_id"], r["protein_accession"] or None)]
        assert r["organism"] == ref.species, (r["step_id"], r["protein_accession"])
        assert r["gene_name"] == ref.gene_name, (r["step_id"], r["protein_accession"])


def test_kb_steps_match_yaml_steps(kb):
    yaml_steps = {s.id for s in load_pathway(YAML).steps}
    assert set(kb["step"]["step_id"]) == yaml_steps | T4_ONLY_STEPS
    assert not (yaml_steps & T4_ONLY_STEPS)


def test_yaml_step_reactions_are_the_kb_reactions(kb):
    sr = kb["step_reaction"]
    for st in load_pathway(YAML).steps:
        assert set(st.reactions) == set(sr[sr["step_id"] == st.id]["reaction_id"]), st.id


def test_masses_are_computed_from_formulas(kb):
    for _, r in kb["compound"].iterrows():
        assert r["neutral_exact_mass"] == f"{formula_mass(r['formula']):.4f}", r["chebi_id"]
        assert r["inchikey"].strip(), r["chebi_id"]


def test_phytoene_desaturase_two_reactions_with_organism_bound_tiers(kb):
    """Each tier is bound to one protein and one paper. The five-step reaction is T1a for
    Neurospora AL-1 only. The four-step reaction is T1a for Phycomyces carB; the abstracts read
    on 2026-11-21 also show lycopene for X. dendrorhous crtI and as a major product of AL-1."""
    sr = kb["step_reaction"]
    pds = sr[sr["step_id"] == "phytoene_desaturase"]
    assert set(pds["reaction_id"]) == {"RHEA:15585", "RHEA:30975"}
    five = pds[(pds["reaction_id"] == "RHEA:30975") & (pds["tier"] == "T1a")]
    assert list(five["protein_accession"]) == ["P21334"]
    assert list(five["pmid"]) == ["11017770"]
    four = pds[(pds["reaction_id"] == "RHEA:15585") & (pds["tier"] == "T1a")]
    assert dict(zip(four["protein_accession"], four["pmid"])) == {
        "P54982": "9079885", "O13506": "10397832", "P21334": "11017770"}
    # every desaturase link names its own protein and paper; no tier is shared by orthology
    assert not pds.duplicated(["reaction_id", "protein_accession"]).any()
    assert (pds["pmid"] != "").all()
    # the Phycomyces enzyme has no five-step link
    assert "P54982" not in set(pds[pds["reaction_id"] == "RHEA:30975"]["protein_accession"])
    rx = kb["reaction"].set_index("reaction_id")
    assert rx.loc["RHEA:15585", "ec"] == "1.3.99.31"
    assert rx.loc["RHEA:30975", "ec"] == "1.3.99.30"


def test_rhea_56772_not_used(kb):
    for col in ("reaction_id", "rhea_id"):
        assert not kb["reaction"][col].str.contains("56772").any()


def test_astaxanthin_route_is_a_grid(kb):
    sr = kb["step_reaction"]
    ast = sr[sr["step_id"] == "astaxanthin_synthase"]
    rids = set(ast["reaction_id"])
    subs, prods = _parts(kb, "substrate"), _parts(kb, "product")
    # one compound (phoenicoxanthin) is made by two reactions with different substrates
    made_by = {}
    for rid in rids:
        for c in prods.get(rid, ()):
            made_by.setdefault(c, set()).add(frozenset(subs.get(rid, ())))
    grid_nodes = [c for c, s in made_by.items() if len(s) >= 2]
    assert len(grid_nodes) >= 1
    # astaxanthin is made from a grid node
    to_ast = [rid for rid in rids if ASTAXANTHIN in prods.get(rid, ())]
    assert to_ast and all(subs[rid] & set(grid_nodes) for rid in to_ast)
    # the two Rhea ketolation reactions are T1a; no link of this step is T1a outside them
    t1a = set(ast[ast["tier"] == "T1a"]["reaction_id"])
    assert t1a == {"RHEA:55660", "RHEA:55664"}
    assert "T1b" in set(ast["tier"])
    assert set(ast["protein_accession"]) == {"Q3HR17"}


def test_crtr_links_are_separate(kb):
    sr = kb["step_reaction"]
    assert len(sr[sr["step_id"] == "crtr"]) >= 1
    assert set(sr[sr["step_id"] == "crtr"]["pmid"]) == {"18837978"}


def test_curated_reaction_ids(kb):
    for _, r in kb["reaction"].iterrows():
        if r["curated"] == "true":
            assert r["reaction_id"].startswith("CUR:carotenoid:"), r["reaction_id"]
            assert r["rhea_id"] == ""
        else:
            assert r["reaction_id"] == r["rhea_id"] and r["rhea_id"].startswith("RHEA:")


def test_torularhodin_is_t4_and_not_reported_possible(kb):
    sr = kb["step_reaction"]
    t4 = sr[sr["step_id"] == "torularhodin_formation"]
    assert len(t4) >= 1 and set(t4["tier"]) == {"T4"}
    assert path_possible(kb, ["torularhodin_formation"]).empty


def test_neurosporaxanthin_path_is_possible_at_t1a(kb):
    got = path_possible(kb, ["cart", "ylo1"])
    assert list(got["min_tier"]) == ["T1a"]


def test_no_kegg_pubchem_column(kb):
    assert list(kb["compound"].columns) == ["chebi_id", "inchikey", "inchikey_block1", "smiles",
                                            "formula", "neutral_exact_mass", "kegg_id",
                                            "np_class"]
