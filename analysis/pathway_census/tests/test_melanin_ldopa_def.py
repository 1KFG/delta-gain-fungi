"""Tests for the L-DOPA melanin / other melanin routes definition and its provenance (Task 5)."""
import csv
import subprocess
import sys
from pathlib import Path

from pathway_census.defs import load_pathway, validate

ROOT = Path(__file__).resolve().parent.parent
YAML = ROOT / "pathways" / "melanin_ldopa.yaml"
PROV = ROOT / "refs" / "melanin_ldopa" / "provenance.tsv"
PROV_COLS = ["step_id", "accession", "species", "gene_name", "region", "pmid",
             "source_database", "retrieved_on", "verified"]
DECOY_DIR = PROV.parent / "decoys"
DECOY_PROV = DECOY_DIR / "provenance.tsv"
DECOY_COLS = ["step_id", "accession", "species", "gene_name", "pmid", "reason",
              "source_database", "retrieved_on", "verified"]
UNVERIFIED = PROV.parent / "UNVERIFIED.md"
STEPS = {"tyrosinase", "laccase_diphenol", "hppd", "hmga", "aspulvinone_e_synthase",
         "asp_melanin_tyrosinase"}
ROUTES = {
    "tyrosinase_dopa": ["tyrosinase"],
    "laccase_exogenous": ["laccase_diphenol"],
    "pyomelanin": ["hppd"],
    "tyrosine_degradation": ["hppd", "hmga"],
    "asp_melanin": ["aspulvinone_e_synthase", "asp_melanin_tyrosinase"],
}
# Steps for which the brief asks for decoys (requirement 4), and the second tyrosinase-family step.
DECOY_STEPS = {"tyrosinase", "laccase_diphenol", "asp_melanin_tyrosinase"}
FLAG = "exogenous_substrate: true"
TYR_RHEA = {"RHEA:18117", "RHEA:34287"}


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
    assert p.id == "melanin_ldopa" and p.verified is False


def test_step_ids_and_routes():
    assert set(_steps()) == STEPS
    assert {rid: r.step_ids for rid, r in _routes().items()} == ROUTES


# ---- requirement 1: tyrosinase and laccase are separate steps ----

def test_tyrosinase_and_laccase_are_distinct_steps():
    s = _steps()
    tyr, lac = s["tyrosinase"], s["laccase_diphenol"]
    assert tyr.ec == "1.14.18.1" and lac.ec == "1.10.3.2"
    a = {r.accession for r in tyr.refs}
    b = {r.accession for r in lac.refs}
    assert a and b and not (a & b)
    assert not (set(tyr.domains_any) & set(lac.domains_any))
    assert tyr.domains_any == ["PF00264"]
    assert lac.domains_any == ["PF07732", "PF00394", "PF07731"]
    assert not (set(tyr.reactions) & set(lac.reactions))


def test_tyrosinase_and_laccase_are_never_in_one_route():
    """Alternatives are separate routes, not consecutive steps of one route."""
    for r in _routes().values():
        assert not {"tyrosinase", "laccase_diphenol"} <= set(r.step_ids), r.id
        assert not {"tyrosinase", "asp_melanin_tyrosinase"} <= set(r.step_ids), r.id


def test_tyrosinase_references_are_the_named_enzymes():
    refs = {r.accession: r for r in _steps()["tyrosinase"].refs}
    # Neurospora crassa T, Aspergillus oryzae melB (ORF name only in UniProt), Agaricus bisporus
    assert refs["P00440"].gene_name == "T" and refs["P00440"].species.startswith("Neurospora crassa")
    # melB: A0A1S9DK56 has the sequence of the assayed protein (PDB 3W6W). Q2UP46, the RIB40
    # gene model of the locus, is not a reference (fix round 1, I1).
    assert refs["A0A1S9DK56"].gene_name == "OAory_01107480"
    assert refs["A0A1S9DK56"].species == "Aspergillus oryzae" and refs["A0A1S9DK56"].verified
    assert "Q2UP46" not in refs
    text = UNVERIFIED.read_text(encoding="utf-8")
    for word in ("3W6W", "615 of 616", "Q2UP46", "not the assayed sequence"):
        assert word in text, word
    assert {"C7FF04", "C7FF05"} <= set(refs)
    assert refs["C7FF05"].gene_name == "PPO4" and refs["C7FF05"].species == "Agaricus bisporus"
    assert _steps()["tyrosinase"].reactions == ["RHEA:18117", "RHEA:34287"]


def test_laccase_references_are_cryptococcus_lac1_and_lac2():
    refs = {r.accession: r for r in _steps()["laccase_diphenol"].refs}
    assert set(refs) == {"Q55P57", "J9VY90", "J9VQZ4"}
    assert refs["Q55P57"].gene_name == "LAC1" and refs["Q55P57"].pmid == "8300520"
    assert refs["Q55P57"].species == "Cryptococcus deneoformans (strain B-3501A)"
    assert refs["J9VY90"].gene_name == "LAC1" and refs["J9VQZ4"].gene_name == "LAC2"
    for acc in ("J9VY90", "J9VQZ4"):
        assert refs[acc].species.startswith("Cryptococcus neoformans (strain H99")


def test_laccase_step_carries_no_tyrosinase_reaction():
    lac = _steps()["laccase_diphenol"]
    assert not (set(lac.reactions) & TYR_RHEA)
    assert all(x.startswith("CUR:melanin_ldopa:") for x in lac.reactions)
    note = lac.lineage_note
    assert "does not hydroxylate" in note and "tyrosine" in note


def test_annotation_that_separates_tyrosinase_and_laccase_is_stated():
    s = _steps()
    for sid in ("tyrosinase", "laccase_diphenol"):
        note = s[sid].lineage_note
        for word in ("EC 1.14.18.1", "EC 1.10.3.2", "RHEA:18117", "RHEA:11276", "PF00264",
                     "PF07732"):
            assert word in note, (sid, word)


# ---- requirement 2: the Cryptococcus route is substrate-dependent ----

def test_exogenous_substrate_flag_is_on_the_cryptococcus_route_only():
    """The schema has no field for an exogenous substrate. The flag is the text
    'exogenous_substrate: true' at the start of the route note."""
    routes = _routes()
    flagged = [rid for rid, r in routes.items() if FLAG in (r.note or "")]
    assert flagged == ["laccase_exogenous"]
    note = routes["laccase_exogenous"].note
    assert note.startswith(FLAG)
    assert "Cryptococcus" in routes["laccase_exogenous"].name
    for word in ("17098915", "only from exogenous substrates", "L-DOPA", "dopamine",
                 "homogentisic acid"):
        assert word in note, word
    for rid, r in routes.items():
        if rid != "laccase_exogenous":
            assert "exogenous_substrate: false" in r.note, rid


def test_exogenous_substrate_is_in_the_laccase_lineage_note():
    note = _steps()["laccase_diphenol"].lineage_note
    assert FLAG in note
    for word in ("17098915", "only from exogenous substrates", "8300520", "Cryptococcus"):
        assert word in note, word
    # no other step is flagged
    for sid, st in _steps().items():
        if sid != "laccase_diphenol":
            assert FLAG not in (st.lineage_note or ""), sid


# ---- requirement 3: other melanin routes ----

def test_pyomelanin_steps():
    s = _steps()
    assert s["hppd"].ec == "1.13.11.27" and s["hmga"].ec == "1.13.11.5"
    assert {r.accession for r in s["hppd"].refs} == {"Q4WHU1"}
    assert "Q4WHT9" in {r.accession for r in s["hmga"].refs}
    assert s["hppd"].reactions == ["RHEA:16189"] and s["hmga"].reactions == ["RHEA:15449"]
    assert "19028908" in s["hppd"].lineage_note and "19028908" in s["hmga"].lineage_note
    # the census cannot tell pigment production from tyrosine catabolism
    assert "catabolism" in s["hppd"].lineage_note
    assert "removes" in s["hmga"].lineage_note
    assert "not a pigment route" in _routes()["tyrosine_degradation"].note


def test_asp_melanin_steps():
    s = _steps()
    mel, tyrp = s["aspulvinone_e_synthase"], s["asp_melanin_tyrosinase"]
    assert {r.accession for r in mel.refs} == {"Q0CRX1"}
    assert {r.accession for r in tyrp.refs} == {"Q0CRX0"}
    assert mel.reactions == ["RHEA:63824"]
    assert tyrp.reactions == ["RHEA:74195", "RHEA:74199"]
    assert mel.ec is None and tyrp.ec is None
    assert "27133313" in mel.lineage_note and "Aspergillus terreus" in mel.lineage_note
    # apvA does the same reaction for another product class
    assert "apvA" in mel.lineage_note and "Q0CWD0" in mel.lineage_note


def test_same_gene_name_for_different_enzymes_is_noted():
    """melB names a tyrosinase of Aspergillus oryzae and is a UniProt synonym of tyrP of
    Aspergillus terreus. LAC1 names different oxidases in Cryptococcus and Verticillium."""
    s = _steps()
    for sid in ("tyrosinase", "asp_melanin_tyrosinase"):
        assert "melB" in s[sid].lineage_note and "only with its organism" in s[sid].lineage_note
    assert "LAC1" in s["laccase_diphenol"].lineage_note
    assert "only with its organism" in s["laccase_diphenol"].lineage_note


def test_dopachrome_onward_is_not_a_yaml_step():
    """Dopaquinone to dopachrome and onward is non-enzymatic in the table. The T4 step is in
    kb/melanin_ldopa/step.tsv only (same convention as carotenoid torularhodin_formation)."""
    assert "dopachrome_onward" not in _steps()
    notes = _steps()["tyrosinase"].lineage_note + _routes()["tyrosinase_dopa"].note
    assert "non-enzymatic" in notes and "T4" in notes


# ---- requirement 4: decoys ----

def test_decoys_are_not_references_of_any_step_and_are_in_provenance():
    prov = {(r["step_id"], r["accession"]): r for r in _rows(DECOY_PROV, DECOY_COLS)}
    all_refs = {r.accession for st in load_pathway(YAML).steps for r in st.refs}
    yaml_decoys = set()
    for st in load_pathway(YAML).steps:
        own = {r.accession for r in st.refs}
        for d in st.decoys:
            assert d.accession, st.id
            assert d.accession not in own, (st.id, d.accession)
            assert d.accession not in all_refs, (st.id, d.accession)
            row = prov[(st.id, d.accession)]
            assert row["reason"].strip(), (st.id, d.accession)
            assert row["verified"] == str(d.verified).lower()
            assert row["species"] == d.species and row["gene_name"] == d.gene_name
            assert row["pmid"] == (d.pmid or "")
            yaml_decoys.add((st.id, d.accession))
    assert yaml_decoys == set(prov)


def test_decoy_counts_and_fasta():
    s = _steps()
    for sid in DECOY_STEPS:
        assert 3 <= len(s[sid].decoys) <= 8, (sid, len(s[sid].decoys))
    for sid, st in s.items():
        if not st.decoys:
            assert not (DECOY_DIR / f"{sid}.faa").exists(), sid
            continue
        seqs = _fasta(DECOY_DIR / f"{sid}.faa")
        assert set(seqs) == {d.accession for d in st.decoys}, sid
        assert all(seqs.values()), sid


def test_decoy_reasons_state_the_separating_annotation():
    """A decoy must not do the reference reaction. The reason names the EC or Rhea annotation
    (or the assay result) that separates the decoy from the step."""
    rows = _rows(DECOY_PROV, DECOY_COLS)
    for r in rows:
        if r["step_id"] == "tyrosinase":
            assert "no EC 1.14.18.1" in r["reason"], r["accession"]
            assert "RHEA:18117" in r["reason"] and "RHEA:34287" in r["reason"], r["accession"]
        elif r["step_id"] == "laccase_diphenol":
            assert "no EC 1.10.3.2" in r["reason"] and "RHEA:11276" in r["reason"], r["accession"]
        elif r["step_id"] == "asp_melanin_tyrosinase":
            assert "RHEA:74195" in r["reason"] and "RHEA:74199" in r["reason"], r["accession"]
    tyr = {r["accession"]: r for r in rows if r["step_id"] == "tyrosinase"}
    # GriF is excluded from the monophenolase reaction only; abr1 is not shown to lack
    # diphenol oxidation (fix round 1, I4)
    assert "monophenolase reaction only" in tyr["B1VTI5"]["reason"]
    assert "not shown to lack diphenol oxidation" in tyr["B1VTI5"]["reason"]
    for rid in ("RHEA:41420", "RHEA:41424", "RHEA:40963"):
        assert rid in tyr["B1VTI5"]["reason"], rid
    abr1 = [r for r in rows if r["accession"] == "Q4WZB4"][0]
    assert "14970241" in abr1["reason"] and "not shown to lack diphenol oxidation" in abr1["reason"]
    text = UNVERIFIED.read_text(encoding="utf-8")
    assert "B1VTI5" in text and "Q4WZB4" in text and "14970241" in text
    # two decoys have an assay that excludes the reference substrates (abstracts read)
    assert tyr["Q5BGU9"]["pmid"] == "2126551" and "not tyrosine" in tyr["Q5BGU9"]["reason"]
    assert "19798497" in tyr["Q2UNF9"]["reason"] and "no activity" in tyr["Q2UNF9"]["reason"]
    lac = {r["accession"] for r in rows if r["step_id"] == "laccase_diphenol"}
    assert "P38993" in lac


def test_no_laccase_and_no_tyrosinase_is_a_decoy_of_its_own_class():
    """An entry annotated EC 1.14.18.1 is not a tyrosinase decoy; an entry annotated
    EC 1.10.3.2 is not a laccase decoy (ustQ, yA and the lignin laccases are excluded)."""
    s = _steps()
    assert not {"B8NM74", "A0A8E5HUS5"} & {d.accession for d in s["tyrosinase"].decoys}
    assert not {"P17489", "P06811", "Q70KY3"} & {d.accession for d in s["laccase_diphenol"].decoys}
    assert "ustQ" in s["tyrosinase"].lineage_note
    text = UNVERIFIED.read_text(encoding="utf-8")
    assert "B8NM74" in text and "P17489" in text


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


def test_q2tl94_is_not_verified_because_uniprot_gives_no_gene_name():
    refs = {r.accession: r for r in _steps()["tyrosinase"].refs}
    assert refs["Q2TL94"].verified is False
    assert sum(not r.verified for st in load_pathway(YAML).steps for r in st.refs) == 1
    assert all(d.verified for st in load_pathway(YAML).steps for d in st.decoys)
    row = [r for r in _rows(PROV, PROV_COLS) if r["accession"] == "Q2TL94"][0]
    assert row["verified"] == "false"
    text = UNVERIFIED.read_text(encoding="utf-8")
    table = text.split("## Refs with verified: false")[1].split("## Refs missing")[0]
    for word in ("Q2TL94", "no gene name", "16151802", "USER DECISION"):
        assert word in table, word
    assert "not a UniProt gene name" in YAML.read_text(encoding="utf-8").split(
        "accession: Q2TL94")[0][-500:]


def test_references_with_a_gap_are_listed_in_unverified_md():
    """A0A1S9DK56 (melB) has an ORF name only; the name melB is not in the entry. P00440,
    C7FF05, J9VY90 and J9VQZ4 do not cite the paper of their link."""
    text = UNVERIFIED.read_text(encoding="utf-8")
    for acc in ("Q2UP46", "A0A1S9DK56", "P00440", "C7FF05", "J9VY90", "J9VQZ4", "Q2TL94"):
        assert acc in text, acc


def test_user_decisions_are_written_down():
    text = UNVERIFIED.read_text(encoding="utf-8")
    assert text.count("USER DECISION") >= 4
    # product-less route: the TyrP products have charged ChEBI entries only
    assert "CHEBI:193114" in text and "CHEBI:193115" in text and "asp_melanin" in text
    # the exogenous-substrate flag is text, not a schema field
    assert "exogenous_substrate" in text and "schema" in text
    # binding rule 1 is applied as written; the substrate-class rule is withdrawn
    assert "applied as written" in text and "withdrawn" in text
    # the dopamine sentence does not name the enzyme
    assert "8300520" in text and "does not name" in text


def test_list_unverified_refs_runs():
    out = subprocess.run([sys.executable, str(ROOT / "bin" / "list_unverified_refs.py"),
                          str(YAML)], capture_output=True, text=True, check=True).stdout
    assert "melanin_ldopa" in out and "Q2TL94" in out


def test_curation_scripts_are_committed_with_a_readme():
    d = ROOT / "curation_scripts" / "melanin_ldopa"
    assert (d / "README.md").exists() and (d / "build.py").exists() and (d / "verify.py").exists()
    text = (d / "README.md").read_text(encoding="utf-8")
    assert "not run by the tests" in text and "Task 7" in text
    for f in d.glob("*.py"):
        src = f.read_text(encoding="utf-8")
        assert "/bigdata/" not in src and "/rhome/" not in src, f.name
