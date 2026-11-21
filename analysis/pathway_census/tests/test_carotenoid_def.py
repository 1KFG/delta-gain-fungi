"""Tests for the carotenoid pathway definition and its reference provenance (Task 3)."""
import csv
import subprocess
import sys
from pathlib import Path

from pathway_census.defs import load_pathway, validate

ROOT = Path(__file__).resolve().parent.parent
YAML = ROOT / "pathways" / "carotenoid.yaml"
PROV = ROOT / "refs" / "carotenoid" / "provenance.tsv"
PROV_COLS = ["step_id", "accession", "species", "gene_name", "region", "pmid",
             "source_database", "retrieved_on", "verified"]
DECOY_DIR = PROV.parent / "decoys"
DECOY_PROV = DECOY_DIR / "provenance.tsv"
DECOY_COLS = ["step_id", "accession", "species", "gene_name", "pmid", "reason",
              "source_database", "retrieved_on", "verified"]
UNVERIFIED = PROV.parent / "UNVERIFIED.md"
REQUIRED_STEPS = {"ggpp_synthase", "phytoene_synthase", "phytoene_desaturase",
                  "lycopene_cyclase", "carx", "cart", "ylo1", "astaxanthin_synthase",
                  "crtr", "cars", "acaa", "tsp3", "tsp2", "tdh"}
# Steps for which the brief names a decoy family (requirement 6): 3 to 8 decoys each.
DECOY_STEPS = {"ggpp_synthase", "phytoene_synthase", "phytoene_desaturase", "carx", "cart",
               "ylo1", "astaxanthin_synthase"}
FUSED = ("phytoene_synthase", "lycopene_cyclase")


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


def _prov_rows():
    with open(PROV, encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == PROV_COLS
        return list(rdr)


def test_carotenoid_yaml_validates():
    assert validate(load_pathway(YAML)) == []


def test_pathway_not_verified():
    assert load_pathway(YAML).verified is False


def test_step_ids_cover_spec_list():
    assert REQUIRED_STEPS <= set(_steps())


def test_fused_psy_cyclase_declared():
    s = _steps()
    assert "lycopene_cyclase" in s["phytoene_synthase"].fused_with
    assert "phytoene_synthase" in s["lycopene_cyclase"].fused_with


def test_every_required_step_has_a_ref():
    s = _steps()
    for sid in REQUIRED_STEPS:
        assert s[sid].refs, sid


def test_provenance_matches_yaml_refs():
    prov = {(r["step_id"], r["accession"]) for r in _prov_rows()}
    for s in load_pathway(YAML).steps:
        for r in s.refs:
            if r.accession:
                assert (s.id, r.accession) in prov, (s.id, r.accession)


def test_provenance_verified_agrees_with_yaml():
    prov = {(r["step_id"], r["accession"]): r for r in _prov_rows()}
    for s in load_pathway(YAML).steps:
        for r in s.refs:
            if r.accession:
                row = prov[(s.id, r.accession)]
                assert row["verified"] == str(r.verified).lower(), (s.id, r.accession)
                assert row["pmid"] == (r.pmid or ""), (s.id, r.accession)


def test_ref_fasta_ids_are_yaml_accessions():
    """Each step FASTA holds only accessions declared for that step in the YAML."""
    s = _steps()
    fdir = PROV.parent
    for sid, step in s.items():
        accs = {r.accession for r in step.refs if r.accession}
        f = fdir / f"{sid}.faa"
        if not accs:
            continue
        assert f.exists(), f
        ids = {ln[1:].split()[0] for ln in f.read_text(encoding="utf-8").splitlines()
               if ln.startswith(">")}
        assert ids <= accs, (sid, ids - accs)


def test_list_unverified_refs_runs():
    out = subprocess.run([sys.executable, str(ROOT / "bin" / "list_unverified_refs.py"),
                          str(YAML)], capture_output=True, text=True, check=True).stdout
    assert "carotenoid" in out


# ---- Stage 1 (S3) curation requirements ----

def test_phytoene_desaturase_has_null_ec_and_two_reactions():
    st = _steps()["phytoene_desaturase"]
    assert st.ec is None
    assert sorted(st.reactions) == ["RHEA:15585", "RHEA:30975"]


def test_astaxanthin_synthase_ref_is_q3hr17_and_ec_null():
    st = _steps()["astaxanthin_synthase"]
    assert st.ec is None
    assert "Q3HR17" in {r.accession for r in st.refs}
    assert all(r.accession for r in st.refs)
    seqs = _fasta(PROV.parent / "astaxanthin_synthase.faa")
    assert len(seqs["Q3HR17"]) == 557


def test_crtr_is_a_separate_auxiliary_step():
    s = _steps()
    assert s["crtr"].role == "auxiliary"
    assert {r.pmid for r in s["crtr"].refs} == {"18837978"}
    assert not ({r.accession for r in s["crtr"].refs}
                & {r.accession for r in s["astaxanthin_synthase"].refs})


def test_carx_has_ustilago_reference():
    refs = _steps()["carx"].refs
    assert any(r.species.startswith(("Ustilago maydis", "Mycosarcoma maydis"))
               and r.pmid == "19584000" for r in refs)
    assert any(r.accession == "S0ENQ2" for r in refs)


def test_kept_references_still_present():
    s = _steps()
    assert {"A1KQY3", "A7UXI1"} <= {r.accession for r in s["cart"].refs}
    assert {"Q1K615", "F6IBC7"} <= {r.accession for r in s["ylo1"].refs}


def test_every_ref_and_decoy_has_gene_name_and_species():
    """Gene names are keyed by organism (requirement 8)."""
    for st in load_pathway(YAML).steps:
        for r in list(st.refs) + list(st.decoys):
            assert r.gene_name, (st.id, r.accession)
            assert r.species, (st.id, r.accession)


def test_unverified_refs_are_listed_in_unverified_md():
    text = UNVERIFIED.read_text(encoding="utf-8")
    for st in load_pathway(YAML).steps:
        for r in list(st.refs) + list(st.decoys):
            if not r.verified:
                assert (r.accession or r.pmid) in text, (st.id, r.accession, r.pmid)


def test_routes_declared():
    p = load_pathway(YAML)
    ids = {r.id for r in p.routes}
    assert {"neurosporaxanthin", "astaxanthin", "retinal", "trisporoid"} <= ids
    ast = next(r for r in p.routes if r.id == "astaxanthin")
    assert {"astaxanthin_synthase", "crtr"} <= set(ast.step_ids)


def _full_regions():
    """(step, accession) -> (start, end) on the full UniProt sequence, from provenance.tsv."""
    out = {}
    for r in _prov_rows():
        if r["region"]:
            a, b = r["region"].split("-")
            out[(r["step_id"], r["accession"])] = (int(a), int(b))
    return out


def test_split_fasta_files_differ_and_regions_lie_inside():
    """`defs.validate` checks a ref region against the record in <step>.faa. The records of the
    fused proteins are cut to the domain, so the YAML region is in record coordinates and covers
    the whole record. The coordinates on the full protein are in provenance.tsv."""
    s = _steps()
    psy = _fasta(PROV.parent / "phytoene_synthase.faa")
    cyc = _fasta(PROV.parent / "lycopene_cyclase.faa")
    assert psy and cyc
    assert (PROV.parent / "phytoene_synthase.faa").read_bytes() != \
        (PROV.parent / "lycopene_cyclase.faa").read_bytes()
    assert not (set(psy.values()) & set(cyc.values()))
    full_reg = _full_regions()
    for sid, seqs in (("phytoene_synthase", psy), ("lycopene_cyclase", cyc)):
        for r in s[sid].refs:
            if r.accession not in seqs:
                assert r.verified is False, (sid, r.accession)  # dropped from the split file
                continue
            assert r.region == (1, len(seqs[r.accession])), (sid, r.accession)
            start, end = full_reg[(sid, r.accession)]
            assert end - start + 1 == len(seqs[r.accession]), (sid, r.accession)
    # the two regions of one fused protein do not overlap on the full protein
    for acc in psy.keys() & cyc.keys():
        a, b = full_reg[("phytoene_synthase", acc)], full_reg[("lycopene_cyclase", acc)]
        assert a[1] < b[0] or b[1] < a[0], acc


def test_split_regions_are_substrings_of_the_full_protein():
    """The cut record equals full_sequence[start-1:end]. The full sequences are kept in
    `full_length/fused.faa` (a subdirectory, so no tool reads it as a step FASTA)."""
    full = _fasta(PROV.parent / "full_length" / "fused.faa")
    full_reg = _full_regions()
    assert len(full) == 5
    for sid in FUSED:
        seqs = _fasta(PROV.parent / f"{sid}.faa")
        assert set(seqs) <= set(full)
        for acc, seq in seqs.items():
            start, end = full_reg[(sid, acc)]
            assert 1 <= start <= end <= len(full[acc]), (sid, acc)
            assert seq == full[acc][start - 1:end], (sid, acc)


def _decoy_rows():
    with open(DECOY_PROV, encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == DECOY_COLS
        return list(rdr)


def test_decoys_in_provenance_and_not_references():
    prov = {(r["step_id"], r["accession"]): r for r in _decoy_rows()}
    yaml_decoys = set()
    for st in load_pathway(YAML).steps:
        ref_acc = {r.accession for r in st.refs}
        for d in st.decoys:
            assert d.accession, st.id
            assert d.accession not in ref_acc, (st.id, d.accession)
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


def test_cao1_is_a_decoy_of_carx_and_cart():
    s = _steps()
    for sid in ("carx", "cart"):
        assert any(d.accession == "Q7S860" and d.pmid == "23893079" for d in s[sid].decoys), sid


def test_steps_without_decoys_are_explained():
    text = UNVERIFIED.read_text(encoding="utf-8")
    for sid, st in _steps().items():
        if not st.decoys:
            assert f"`{sid}`: no decoy" in text, sid


def test_provenance_gene_and_species_match_yaml():
    prov = {(r["step_id"], r["accession"]): r for r in _prov_rows()}
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            if r.accession:
                row = prov[(st.id, r.accession)]
                assert bool(row["region"]) == (r.region is not None), (st.id, r.accession)
                assert row["gene_name"] == r.gene_name
                assert row["species"] == r.species
