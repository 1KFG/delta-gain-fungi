"""Tests for the DHN-melanin pathway definition and its reference provenance (Task 4)."""
import csv
import subprocess
import sys
from pathlib import Path

from pathway_census.defs import load_pathway, validate

ROOT = Path(__file__).resolve().parent.parent
YAML = ROOT / "pathways" / "melanin_dhn.yaml"
PROV = ROOT / "refs" / "melanin_dhn" / "provenance.tsv"
PROV_COLS = ["step_id", "accession", "species", "gene_name", "region", "pmid",
             "source_database", "retrieved_on", "verified"]
DECOY_DIR = PROV.parent / "decoys"
DECOY_PROV = DECOY_DIR / "provenance.tsv"
DECOY_COLS = ["step_id", "accession", "species", "gene_name", "pmid", "reason",
              "source_database", "retrieved_on", "verified"]
UNVERIFIED = PROV.parent / "UNVERIFIED.md"
PKS_STEPS = {"pks_thn": "thn_synthase", "pks_hexaketide": "hexaketide_pks",
             "pks_heptaketide": "heptaketide_pks"}
STEPS = set(PKS_STEPS.values()) | {"ayg1", "t4hn_reductase", "t3hn_reductase",
                                   "scytalone_dehydratase", "laccase_abr2", "mco_abr1"}
# Steps for which the brief names a decoy family (requirement 6): 3 to 8 decoys each.
DECOY_STEPS = set(PKS_STEPS.values()) | {"t4hn_reductase", "t3hn_reductase", "laccase_abr2",
                                         "mco_abr1"}
TAIL = ["t4hn_reductase", "scytalone_dehydratase", "t3hn_reductase", "scytalone_dehydratase"]


def _fasta(path):
    out, name = {}, None
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = ""
        elif name is not None:
            out[name] += ln.strip()
    return out


def _steps():
    return {s.id: s for s in load_pathway(YAML).steps}


def _routes():
    return {r.id: r for r in load_pathway(YAML).routes}


def _rows(path, cols):
    with open(path, encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == cols
        return list(rdr)


def test_yaml_validates():
    assert validate(load_pathway(YAML)) == []


def test_pathway_not_verified():
    p = load_pathway(YAML)
    assert p.id == "melanin_dhn" and p.verified is False


def test_step_ids():
    assert set(_steps()) == STEPS


# ---- requirement 1: three PKS routes; Ayg1/Yg1 only with the hexaketide and heptaketide PKS ----

def test_three_pks_routes_and_one_polymerisation_route():
    assert set(_routes()) == set(PKS_STEPS) | {"polymerisation"}


def test_each_pks_route_has_exactly_its_own_pks_step_first():
    routes = _routes()
    for rid, pks in PKS_STEPS.items():
        ids = routes[rid].step_ids
        assert ids[0] == pks, rid
        assert set(ids) & set(PKS_STEPS.values()) == {pks}, rid


def test_ayg1_is_only_in_the_hexaketide_and_heptaketide_routes():
    routes = _routes()
    with_ayg1 = {rid for rid, r in routes.items() if "ayg1" in r.step_ids}
    assert with_ayg1 == {"pks_hexaketide", "pks_heptaketide"}
    assert routes["pks_thn"].step_ids == ["thn_synthase"] + TAIL
    assert routes["pks_hexaketide"].step_ids == ["hexaketide_pks", "ayg1"] + TAIL
    assert routes["pks_heptaketide"].step_ids == ["heptaketide_pks", "ayg1"] + TAIL


def test_pks_alternatives_are_never_consecutive_steps():
    for r in _routes().values():
        assert len(set(r.step_ids) & set(PKS_STEPS.values())) <= 1, r.id


def test_pks_references_are_the_named_enzymes():
    s = _steps()
    assert {r.accession for r in s["thn_synthase"].refs} == {"P79068"}
    assert {r.accession for r in s["hexaketide_pks"].refs} == {"Q9Y7A7"}
    assert "Q4WZA8" in {r.accession for r in s["heptaketide_pks"].refs}
    ayg = {r.accession for r in s["ayg1"].refs}
    assert {"Q4WZB3", "Q6BCB5"} <= ayg


def test_pks_refs_do_not_overlap_between_routes():
    s = _steps()
    sets = [{r.accession for r in s[x].refs} for x in PKS_STEPS.values()]
    for i in range(3):
        for j in range(i + 1, 3):
            assert not (sets[i] & sets[j])


# ---- requirement 2: two reductase types ----

def test_two_reductase_steps_with_distinct_reference_sets():
    s = _steps()
    a = {r.accession for r in s["t4hn_reductase"].refs}
    b = {r.accession for r in s["t3hn_reductase"].refs}
    assert a and b and not (a & b)
    # 4HNR/Arp2 type
    assert {"Q9HFV6", "E9QUT3"} <= a
    # 3HNR/THR1 type
    assert {"Q12634", "P87025"} <= b
    assert s["t4hn_reductase"].ec == "1.1.1.252"
    assert s["t3hn_reductase"].ec is None


def test_q12634_is_the_3hnr_not_the_4hnr():
    """The pathway review named Q12634 as 4HNR. The records fetched on 2026-11-21 show that
    Q12634 is the trihydroxynaphthalene reductase (it cites PMID 8112349 and the 3HNR crystal
    structures) and that Q9HFV6 is the 4HNR of PMID 10956664."""
    s = _steps()
    assert "Q12634" not in {r.accession for r in s["t4hn_reductase"].refs}
    assert "Q9HFV6" not in {r.accession for r in s["t3hn_reductase"].refs}


def test_reductase_lineage_notes_warn_about_paralogues():
    s = _steps()
    for sid in ("t4hn_reductase", "t3hn_reductase"):
        note = s[sid].lineage_note
        assert "paralogue" in note and "does not mean" in note, sid


# ---- requirement 3: scytalone dehydratase ----

def test_scytalone_dehydratase_references_and_reactions():
    st = _steps()["scytalone_dehydratase"]
    assert {"Q00455", "O14434", "P56221"} <= {r.accession for r in st.refs}
    assert st.ec == "4.2.1.94"
    assert st.reactions == ["RHEA:24396", "CUR:melanin_dhn:vermelone_dehydration"]


# ---- requirement 4: two oxidase families ----

def test_oxidases_are_two_steps_with_distinct_references():
    s = _steps()
    lac = {r.accession for r in s["laccase_abr2"].refs}
    mco = {r.accession for r in s["mco_abr1"].refs}
    assert "E9RBR0" in lac and "Q4WZB4" in mco and not (lac & mco)
    assert _routes()["polymerisation"].step_ids == ["laccase_abr2", "mco_abr1"]


# ---- requirement 5: lineage restrictions ----

def test_lineage_restrictions_are_recorded():
    notes = " ".join(st.lineage_note or "" for st in load_pathway(YAML).steps)
    assert "Aspergillus terreus" in notes and "Asp-melanin" in notes and "27133313" in notes
    assert "Botrytis" in notes and "Leotiomycetes" in notes and "26514268" in notes


def test_same_gene_name_for_different_enzymes_is_noted():
    """PKS1 names the THN synthase of C. lagenarium and the hexaketide synthase of
    Exophiala (Wangiella) dermatitidis."""
    s = _steps()
    for sid in ("thn_synthase", "hexaketide_pks"):
        assert {r.gene_name for r in s[sid].refs} == {"PKS1"}
        assert "PKS1" in s[sid].lineage_note and "same gene name" in s[sid].lineage_note


# ---- requirement 6: decoys ----

def test_wa_is_a_decoy_of_the_thn_and_hexaketide_steps_only():
    """A. nidulans wA (Q03149) makes YWA1 (RHEA:62652), the reference reaction of
    heptaketide_pks. A decoy must not do the reference reaction of its step."""
    s = _steps()
    for sid in ("thn_synthase", "hexaketide_pks"):
        assert "Q03149" in {d.accession for d in s[sid].decoys}, sid
    assert "Q03149" not in {d.accession for d in s["heptaketide_pks"].decoys}
    assert "RHEA:62652" in s["heptaketide_pks"].reactions


def test_decoys_in_provenance_and_not_references():
    prov = {(r["step_id"], r["accession"]): r for r in _rows(DECOY_PROV, DECOY_COLS)}
    all_refs = {r.accession for st in load_pathway(YAML).steps for r in st.refs}
    yaml_decoys = set()
    for st in load_pathway(YAML).steps:
        for d in st.decoys:
            assert d.accession, st.id
            assert d.accession not in all_refs, (st.id, d.accession)
            row = prov[(st.id, d.accession)]
            assert row["reason"].strip(), (st.id, d.accession)
            assert row["verified"] == str(d.verified).lower()
            assert row["species"] == d.species and row["gene_name"] == d.gene_name
            yaml_decoys.add((st.id, d.accession))
    assert yaml_decoys == set(prov)


def test_decoy_counts_and_fasta():
    s = _steps()
    for sid in DECOY_STEPS:
        assert 3 <= len(s[sid].decoys) <= 8, (sid, len(s[sid].decoys))
    for sid, st in s.items():
        if not st.decoys:
            continue
        seqs = _fasta(DECOY_DIR / f"{sid}.faa")
        assert set(seqs) == {d.accession for d in st.decoys}, sid
        assert all(seqs.values()), sid


def test_steps_without_decoys_are_explained():
    text = UNVERIFIED.read_text(encoding="utf-8")
    for sid, st in _steps().items():
        if not st.decoys:
            assert f"`{sid}`: no decoy" in text, sid


# ---- provenance ----

def test_every_step_has_a_ref():
    for sid, st in _steps().items():
        assert st.refs, sid
        assert all(r.accession for r in st.refs), sid


def test_provenance_matches_yaml_refs():
    prov = {(r["step_id"], r["accession"]): r for r in _rows(PROV, PROV_COLS)}
    yaml_refs = set()
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            row = prov[(st.id, r.accession)]
            assert row["verified"] == str(r.verified).lower(), (st.id, r.accession)
            assert row["pmid"] == (r.pmid or ""), (st.id, r.accession)
            assert row["gene_name"] == r.gene_name and row["species"] == r.species
            assert row["retrieved_on"] == "2026-11-21"
            yaml_refs.add((st.id, r.accession))
    assert yaml_refs == set(prov)


def test_ref_fasta_ids_equal_yaml_accessions():
    for sid, st in _steps().items():
        seqs = _fasta(PROV.parent / f"{sid}.faa")
        assert set(seqs) == {r.accession for r in st.refs}, sid
        assert all(seqs.values()), sid


def test_every_ref_and_decoy_has_gene_name_species_and_pmid_rule():
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            assert r.gene_name and r.species and r.pmid, (st.id, r.accession)
        for d in st.decoys:
            assert d.gene_name and d.species, (st.id, d.accession)


def test_unverified_refs_are_listed_in_unverified_md():
    text = UNVERIFIED.read_text(encoding="utf-8")
    for st in load_pathway(YAML).steps:
        for r in list(st.refs) + list(st.decoys):
            if not r.verified:
                assert r.accession in text, (st.id, r.accession)


def test_references_with_a_gap_are_listed_in_unverified_md():
    """Q9HFV6 has no gene name in UniProt; Q12634 and the A. niger entries have only an ORF
    name. The gap is written down for the user."""
    text = UNVERIFIED.read_text(encoding="utf-8")
    for acc in ("Q9HFV6", "Q12634", "A2R3S8", "A2R3T0", "P79068"):
        assert acc in text, acc


def test_list_unverified_refs_runs():
    out = subprocess.run([sys.executable, str(ROOT / "bin" / "list_unverified_refs.py"),
                          str(YAML)], capture_output=True, text=True, check=True).stdout
    assert "melanin_dhn" in out


def test_curation_scripts_are_committed_with_a_readme():
    d = ROOT / "curation_scripts" / "melanin_dhn"
    assert (d / "README.md").exists() and (d / "build.py").exists() and (d / "verify.py").exists()
    text = (d / "README.md").read_text(encoding="utf-8")
    assert "not run by the tests" in text and "Task 7" in text
    for f in d.glob("*.py"):
        assert "/bigdata/" not in f.read_text(encoding="utf-8"), f.name
