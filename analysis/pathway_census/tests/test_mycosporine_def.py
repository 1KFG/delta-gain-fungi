"""Tests for the mycosporine / mycosporine-like amino acid definition and its provenance (Task 6)."""
import csv
import subprocess
import sys
from pathlib import Path

from pathway_census.defs import load_pathway, validate

ROOT = Path(__file__).resolve().parent.parent
YAML = ROOT / "pathways" / "mycosporine.yaml"
REFS = ROOT / "refs" / "mycosporine"
PROV = REFS / "provenance.tsv"
PROV_COLS = ["step_id", "accession", "species", "gene_name", "region", "pmid",
             "source_database", "retrieved_on", "verified"]
DECOY_DIR = REFS / "decoys"
DECOY_PROV = DECOY_DIR / "provenance.tsv"
DECOY_COLS = ["step_id", "accession", "species", "gene_name", "pmid", "reason",
              "source_database", "retrieved_on", "verified"]
UNVERIFIED = REFS / "UNVERIFIED.md"
LINEAGES = REFS / "fungal_lineages.tsv"
LINEAGE_COLS = ["organism", "group", "compound", "pmid", "evidence", "note"]

ENZYME_STEPS = {"ddgs", "omt", "atp_grasp_ligase"}          # steps with a reference protein
T4_STEPS = {"mg_glucosyltransferase", "glutamicol_formation"}  # no gene identified
ROUTES = {
    "mycosporine_glutaminol": ["ddgs", "omt", "atp_grasp_ligase"],
    "mycosporine_glutaminol_glucoside": ["ddgs", "omt", "atp_grasp_ligase",
                                         "mg_glucosyltransferase"],
    "mycosporine_glutamicol": ["ddgs", "omt", "glutamicol_formation"],
}
# Phaffia rhodozyma UCD 67-385: GenBank proteins of the cluster OQ547787 (no UniProtKB entry).
FUNGAL = {"ddgs": "WEV88679.1", "omt": "WEV88680.1", "atp_grasp_ligase": "WEV88681.1"}
FUNGAL_PMID = "36983003"
# Cyanobacterial references (Nostoc punctiforme, Trichormus variabilis). Never fungal evidence.
BACTERIAL = {"ddgs": {"B2J6X9", "Q3M6C3"}, "omt": {"B2J6X8"}, "atp_grasp_ligase": {"B2J6X7"}}


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


def _is_fungal(species):
    return species.startswith("Phaffia rhodozyma")


def test_yaml_validates():
    assert validate(load_pathway(YAML)) == []


def test_pathway_not_verified():
    p = load_pathway(YAML)
    assert p.id == "mycosporine" and p.verified is False


def test_step_ids_and_routes():
    assert set(_steps()) == ENZYME_STEPS | T4_STEPS
    assert {rid: r.step_ids for rid, r in _routes().items()} == ROUTES


# ---- the step set starts from the fungal literature ----

def test_every_enzyme_step_has_the_phaffia_reference_and_its_paper():
    s = _steps()
    for sid, acc in FUNGAL.items():
        refs = {r.accession: r for r in s[sid].refs}
        assert acc in refs, sid
        assert refs[acc].species == "Phaffia rhodozyma (strain UCD 67-385)", sid
        assert refs[acc].pmid == FUNGAL_PMID, sid
    assert {sid: {r.accession: r.gene_name for r in s[sid].refs}[acc]
            for sid, acc in FUNGAL.items()} == {"ddgs": "DDGS", "omt": "OMT",
                                                "atp_grasp_ligase": "ATPG"}


def test_step_set_is_not_a_copy_of_the_cyanobacterial_gene_list():
    """Cyanobacteria: MysA (DDGS), MysB (O-methyltransferase), MysC (ATP-grasp, glycine), then
    MysD (D-Ala-D-Ala ligase-like) or MysE (NRPS-like) for the second amino acid. The fungal set
    has no second-amine step, its ATP-grasp step makes mycosporine-glutaminol, and it has a
    glucosylation step. No step rests on bacterial references alone."""
    s = _steps()
    for sid in s:
        assert not sid.lower().startswith("mys"), sid
        for word in ("nrps", "second_amin", "d_ala", "shinorine"):
            assert word not in sid.lower(), sid
    for sid, st in s.items():
        if not st.refs:
            assert sid in T4_STEPS, sid
            continue
        fungal = [r for r in st.refs if _is_fungal(r.species)]
        assert fungal, sid                      # at least one fungal reference
        assert {r.accession for r in st.refs if not _is_fungal(r.species)} == BACTERIAL[sid]
    all_refs = {r.accession for st in s.values() for r in st.refs}
    # MysD of Nostoc punctiforme and the NRPS-like enzyme of Anabaena variabilis are not references
    assert not {"B2J6X6", "Q3M6C6"} & all_refs
    note = s["atp_grasp_ligase"].lineage_note
    for word in ("glutaminol", "glycine", "MysC", "MysD", "NRPS-like", "not a step"):
        assert word in note, word
    text = UNVERIFIED.read_text(encoding="utf-8")
    assert "not a copy of the cyanobacterial gene list" in text


def test_bacterial_references_are_labelled_as_bacterial():
    s = _steps()
    for sid, accs in BACTERIAL.items():
        note = s[sid].lineage_note
        assert "T2" in note and "cyanobacteri" in note.lower(), sid
        for r in s[sid].refs:
            if r.accession in accs:
                assert r.species.startswith(("Nostoc punctiforme", "Trichormus variabilis")), sid
    text = YAML.read_text(encoding="utf-8")
    assert "never fungal evidence" in text


def test_t4_steps_have_no_reference_and_say_that_the_gene_is_unknown():
    s = _steps()
    for sid in T4_STEPS:
        st = s[sid]
        assert st.refs == [] and st.decoys == [] and st.optional is True, sid
        assert "T4" in st.lineage_note and "not identified" in st.lineage_note, sid
        assert not (REFS / f"{sid}.faa").exists(), sid
    assert "36983003" in s["mg_glucosyltransferase"].lineage_note
    # the enzyme activity of Ascochyta fabae is recorded; whether mycosporin-2 is the aglycone is not
    for word in ("24258137", "Ascochyta fabae", "UDP-glucose", "no sequence is known", "UNVERIFIED"):
        assert word in s["mg_glucosyltransferase"].lineage_note, word
    assert "not part of the" in s["mg_glucosyltransferase"].lineage_note
    text = UNVERIFIED.read_text(encoding="utf-8")
    for sid in T4_STEPS:
        assert f"`{sid}`" in text, sid
    assert "The fungal genes of these steps are unknown" in text


def test_routes_say_what_is_shown_in_fungi():
    r = _routes()
    assert r["mycosporine_glutaminol"].note.startswith("fungal_evidence: true")
    for rid in ("mycosporine_glutaminol_glucoside", "mycosporine_glutamicol"):
        assert r[rid].note.startswith("fungal_evidence: partial"), rid
        assert "T4" in r[rid].note, rid
    for rid in r:
        assert "4-deoxygadusol has no ChEBI" in r[rid].note, rid
    assert "Phaffia rhodozyma" in r["mycosporine_glutaminol"].name


def test_steps_carry_ec_reactions_and_domains():
    s = _steps()
    assert s["ddgs"].ec == "4.2.3.154" and s["ddgs"].reactions == ["RHEA:49560"]
    assert s["omt"].ec is None
    assert s["omt"].reactions == ["CUR:mycosporine:ddg_o_methylation"]
    assert s["atp_grasp_ligase"].reactions == ["CUR:mycosporine:mycosporine_glutaminol_formation",
                                               "CUR:mycosporine:mycosporine_glycine_formation"]
    assert s["mg_glucosyltransferase"].reactions == ["CUR:mycosporine:mg_glucosylation"]
    assert s["glutamicol_formation"].reactions == [
        "CUR:mycosporine:mycosporine_glutamicol_formation"]
    assert s["ddgs"].domains_any == ["PF01761", "PF24621"]
    assert s["omt"].domains_any == ["PF13578", "PF01596"]
    # the Phaffia ATP-grasp ligase has no Pfam family at the gathering cutoff (hmmscan, MEASURED)
    assert s["atp_grasp_ligase"].domains_any == []
    assert "gathering cutoff" in s["atp_grasp_ligase"].lineage_note
    assert "PS50975" in s["atp_grasp_ligase"].lineage_note


def test_gene_names_are_keyed_by_organism():
    for sid in ENZYME_STEPS:
        assert "only with its organism" in _steps()[sid].lineage_note, sid


# ---- decoys ----

DECOYS = {
    "ddgs": {"H2K887", "Q9ZAE9", "C6WFL3", "P07639", "P07547", "Q9S5E2"},
    "omt": {"O42898", "B8NI24", "A0A397HK53", "P21964"},
    "atp_grasp_ligase": {"Q4WD45", "C4JF81", "A0A5N6TQK6", "P07862", "P04425", "B2J6X6"},
}


def test_decoys_are_the_expected_set_and_not_references():
    s = _steps()
    prov = {(r["step_id"], r["accession"]): r for r in _rows(DECOY_PROV, DECOY_COLS)}
    all_refs = {r.accession for st in s.values() for r in st.refs}
    seen = set()
    for sid, st in s.items():
        assert {d.accession for d in st.decoys} == DECOYS.get(sid, set()), sid
        for d in st.decoys:
            assert d.accession not in all_refs, (sid, d.accession)
            row = prov[(sid, d.accession)]
            assert row["reason"].strip() and row["verified"] == str(d.verified).lower()
            assert row["species"] == d.species and row["gene_name"] == d.gene_name
            assert row["pmid"] == (d.pmid or "")
            assert d.gene_name and d.species, (sid, d.accession)
            seen.add((sid, d.accession))
        if st.decoys:
            seqs = _fasta(DECOY_DIR / f"{sid}.faa")
            assert set(seqs) == {d.accession for d in st.decoys} and all(seqs.values()), sid
    assert seen == set(prov)


def test_decoy_reasons_state_the_annotation_and_what_is_not_shown():
    rows = _rows(DECOY_PROV, DECOY_COLS)
    by = {(r["step_id"], r["accession"]): r["reason"] for r in rows}
    for (sid, acc), reason in by.items():
        if sid == "ddgs":
            assert "no EC 4.2.3.154" in reason and "no RHEA:49560" in reason, acc
            assert "PF01761" in reason, acc
        elif sid == "omt":
            assert "PF01596" in reason and "not assayed with demethyl-4-deoxygadusol" in reason
        else:
            assert "PS50975" in reason and "not assayed with 4-deoxygadusol" in reason, acc
        assert ("annotation only" in reason) != ("assay read" in reason), acc
    # sugar phosphate cyclases with another product
    assert "RHEA:44184" in by[("ddgs", "H2K887")] and "EC 4.2.3.152" in by[("ddgs", "H2K887")]
    assert "RHEA:49564" in by[("ddgs", "C6WFL3")] and "22741921" in by[("ddgs", "C6WFL3")]
    assert "RHEA:21968" in by[("ddgs", "P07547")] and "AROM" in by[("ddgs", "P07547")]
    # MysD does another reaction; it is the one decoy that is not verified
    assert "21890703" in by[("atp_grasp_ligase", "B2J6X6")]
    # I4: the probable catechol O-methyltransferase is not excluded as a true member
    o = by[("omt", "O42898")]
    for word in ("not excluded as a true member", "no evidence code", "closer than MysB", "tree placement"):
        assert word in o, word
    assert "USER DECISION 7" in UNVERIFIED.read_text(encoding="utf-8")
    d = {x.accession: x for x in _steps()["atp_grasp_ligase"].decoys}
    assert d["B2J6X6"].verified is False
    assert sum(not x.verified for st in _steps().values() for x in st.decoys) == 1


# ---- provenance ----

def test_provenance_matches_yaml_refs():
    prov = {(r["step_id"], r["accession"]): r for r in _rows(PROV, PROV_COLS)}
    yaml_refs = set()
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            row = prov[(st.id, r.accession)]
            assert row["verified"] == str(r.verified).lower(), (st.id, r.accession)
            assert row["pmid"] == (r.pmid or "") and r.pmid, (st.id, r.accession)
            assert row["gene_name"] == r.gene_name and row["species"] == r.species
            assert r.gene_name and r.species
            assert row["retrieved_on"] == "2026-11-21"
            yaml_refs.add((st.id, r.accession))
    assert yaml_refs == set(prov)
    src = {r["accession"]: r["source_database"] for r in _rows(PROV, PROV_COLS)}
    for acc in FUNGAL.values():
        assert src[acc].startswith("GenBank"), acc


def test_ref_fasta_ids_equal_yaml_accessions():
    for sid, st in _steps().items():
        if not st.refs:
            continue
        seqs = _fasta(REFS / f"{sid}.faa")
        assert set(seqs) == {r.accession for r in st.refs}, sid
        assert all(seqs.values()), sid
    assert len(_fasta(REFS / "ddgs.faa")["WEV88679.1"]) == 477
    assert len(_fasta(REFS / "omt.faa")["WEV88680.1"]) == 284
    assert len(_fasta(REFS / "atp_grasp_ligase.faa")["WEV88681.1"]) == 500


def test_verified_flags_follow_binding_rule_3():
    """verified: true only for a UniProt entry with a gene name that cites the paper and carries
    the EC and the Rhea reaction of the step. The Phaffia proteins have no UniProtKB entry."""
    flags = {r.accession: r.verified for st in _steps().values() for r in st.refs}
    assert flags == {"WEV88679.1": False, "B2J6X9": True, "Q3M6C3": True,
                     "WEV88680.1": False, "B2J6X8": False,
                     "WEV88681.1": False, "A0A9Y1GGR7": False, "B2J6X7": False}
    text = UNVERIFIED.read_text(encoding="utf-8")
    table = text.split("## Refs with verified: false")[1].split("\n## ")[0]
    for acc, ok in flags.items():
        if not ok:
            assert acc in table, acc
    assert "no UniProtKB entry" in table and "OQ547787" in table
    for st in _steps().values():
        for d in st.decoys:
            if not d.verified:
                assert d.accession in text


def test_user_decisions_and_context_are_written_down():
    text = UNVERIFIED.read_text(encoding="utf-8")
    assert text.count("USER DECISION") >= 5
    for n in range(1, 8):
        assert f"USER DECISION {n}" in text, n
    # I2: the OMT and DDGS tiers rest on one sentence and one experiment
    assert "decide the pair together" in text and "no strain with DDGS alone" in text
    # I5: the decoy with an ordered locus name is part of the decision
    assert "C6WFL3 (Amir_2000)" in text
    # M1: the sugar-transporter paper is not a paper about the step
    assert "sugar transporter" in text
    # compounds without a ChEBI id are never invented
    for name in ("4-deoxygadusol", "mycosporine-glycine", "mycosporine-glutaminol-glucoside",
                 "mycosporine-glutamicol-glucoside", "asterina-330"):
        assert name in text, name
    assert "no ChEBI" in text
    # the EBC report is context only
    for word in ("36/87", "16.9", "0 of 82", "176 corpus accessions", "32 of 88", "14.5",
                 "REPORTED", "context only", "not evidence"):
        assert word in text, word
    # what was read
    assert "full text" in text and "abstract only" in text
    # ordered locus names
    assert "ordered locus name" in text


def test_fungal_lineage_table():
    rows = _rows(LINEAGES, LINEAGE_COLS)
    assert len(rows) == 29
    for r in rows:
        assert r["organism"].strip() and r["compound"].strip() and r["pmid"].isdigit(), r
        assert r["evidence"] in ("abstract", "full text"), r
    got = {(r["organism"], r["compound"], r["pmid"]) for r in rows}
    for want in [("Phaffia rhodozyma", "mycosporine-glutaminol-glucoside", "20955516"),
                 ("Phaffia rhodozyma", "mycosporine-glutaminol", "36983003"),
                 ("Ascochyta fabae", "mycosporin-2 (acceptor of a spore glucosyltransferase)", "24258137"),
                 ("Naganishia friedmannii", "mycosporine-glutaminol-glucoside", "40036870"),
                 ("Rhodotorula minuta", "mycosporine-glutaminol-glucoside", "15484286"),
                 ("Cystobasidium keelungensis", "mycosporine-glutamicol-glucoside", "30182254"),
                 ("Hortaea werneckii", "mycosporine-glutaminol-glucoside", "18048939"),
                 ("Colletotrichum graminicola", "mycosporine-glutamine", "31550526")]:
        assert want in got, want
    # negative reports are kept, and nothing is generalised to the kingdom
    assert any(r["compound"] == "none detected" for r in rows)
    text = UNVERIFIED.read_text(encoding="utf-8")
    assert "not generalised to the kingdom" in text


def test_list_unverified_refs_runs():
    out = subprocess.run([sys.executable, str(ROOT / "bin" / "list_unverified_refs.py"),
                          str(YAML)], capture_output=True, text=True, check=True).stdout
    assert "mycosporine" in out and "WEV88679.1" in out and "mg_glucosyltransferase" in out


def test_curation_scripts_are_committed_with_a_readme():
    d = ROOT / "curation_scripts" / "mycosporine"
    assert (d / "README.md").exists() and (d / "build.py").exists() and (d / "verify.py").exists()
    text = (d / "README.md").read_text(encoding="utf-8")
    assert "not run by the tests" in text and "Task 7" in text
    names = {f.name for f in d.glob("*.py")}
    assert not [n for n in names if len(n) <= 5], names      # no a1.py, b1.py, h.py
    for f in d.glob("*.py"):
        assert f.name in text, f.name
        src = f.read_text(encoding="utf-8")
        assert "/bigdata/" not in src and "/rhome/" not in src, f.name


def test_reference_pmids_and_lineage_evidence_after_fix_round_1():
    refs = {r.accession: r for st in _steps().values() for r in st.refs}
    # M1: PMID 23463784 (a sugar-transporter paper) is not the pmid of any reference
    assert "23463784" not in {r.pmid for r in refs.values()}
    assert refs["B2J6X8"].pmid == refs["B2J6X7"].pmid == "21890703"
    src = {r["accession"]: r["source_database"] for r in _rows(PROV, PROV_COLS)}
    for acc in ("B2J6X8", "B2J6X7"):
        assert "not cited by the entry" in src[acc] and "23463784" in src[acc], acc
    # M4: mycosporine-glutaminol in wild-type Phaffia is from the full text
    row = [r for r in _rows(LINEAGES, LINEAGE_COLS)
           if (r["organism"], r["compound"], r["pmid"]) == (
               "Phaffia rhodozyma", "mycosporine-glutaminol", "36983003")]
    assert len(row) == 1 and row[0]["evidence"] == "full text"
    text = UNVERIFIED.read_text(encoding="utf-8")
    for word in ("24258137", "the amide form is the best acceptor", "It is UNVERIFIED whether"):
        assert word in text, word
