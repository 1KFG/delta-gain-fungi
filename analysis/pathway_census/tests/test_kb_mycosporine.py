"""Tests for the mycosporine flat knowledge tables (stage 1, task S6)."""
import hashlib
import re
from pathlib import Path

import pandas as pd
import pytest

from pathway_census.defs import load_pathway
from pathway_census.kb import (EVIDENCE_TIERS, SPECIES_ALIASES, formula_mass, load_kb,
                               path_possible, species_key, validate_kb)

ROOT = Path(__file__).resolve().parent.parent
KB = ROOT / "kb" / "mycosporine"
YAML = ROOT / "pathways" / "mycosporine.yaml"
QUOTES = ROOT / "tests" / "data" / "kb" / "mycosporine_quotes.tsv"
TAXA = ROOT / "tests" / "data" / "kb" / "ncbi_taxonomy_names.tsv"
C = "CUR:mycosporine:"
SH7P, DDG, MG, MGLUOH = "CHEBI:15721", "CHEBI:132137", "CHEBI:144275", "CHEBI:144274"
MSER, GAD, SHI, P334, PAL = ("CHEBI:7039", "CHEBI:172445", "CHEBI:131712", "CHEBI:35671",
                             "CHEBI:182637")
T4_STEPS = {"mg_glucosyltransferase", "glutamicol_formation"}
NO_ENZYME = "not assigned (no enzyme identified)"
PHAFFIA = "xanthophyllomyces dendrorhous"      # species key of Phaffia rhodozyma (taxid 264483)
FUNGAL_PMID = "36983003"
S_DEL = ("Homozygous deletion mutants of the entire cluster, single-gene mutants, and the "
         "Δddgs-/-;Δomt-/- and Δomt-/-;Δatpg-/- double-gene mutants did not produce "
         "mycosporines.")
S_HET = ("Heterologous expression of the DDGS and OMT or DDGS, OMT, and ATPG cDNAs in "
         "Saccharomyces cerevisiae led to 4-deoxygadusol or MG production, respectively.")
S_MYS = ("Heterologous expression in Escherichia coli of three contiguous N. punctiforme genes "
         "(NpR5600, NpR5599, and NpR5598, here named mysA, mysB, and mysC, respectively) led to "
         "the production of mycosporine-glycine, an oxomycosporine.")


@pytest.fixture(scope="module")
def kb():
    return load_kb(KB)


def _parts(kb, rid, role):
    p = kb["reaction_participant"]
    return set(p[(p["reaction_id"] == rid) & (p["role"] == role)]["compound"])


def _routes():
    return {r.id: r.step_ids for r in load_pathway(YAML).routes}


def _links(kb, step):
    sr = kb["step_reaction"]
    sr = sr[sr["step_id"] == step]
    return {(r, a, p): t for r, a, p, t in zip(sr["reaction_id"], sr["protein_accession"],
                                              sr["pmid"], sr["tier"])}


def test_tables_validate(kb):
    assert validate_kb(kb) == []
    for name, df in kb.items():
        assert len(df) > 0, name


def test_every_t1_t2_link_has_pmid_quote_and_assay(kb):
    sr = kb["step_reaction"]
    ev = sr[sr["tier"].isin(EVIDENCE_TIERS)]
    assert len(ev) == 5
    for _, r in ev.iterrows():
        assert r["pmid"].strip().isdigit() and r["quote"].strip() and r["assay"].strip()
    assert set(sr["tier"]) == {"T1a", "T1c", "T2", "T3", "T4"}


def test_one_quoted_sentence_per_link(kb):
    sr = kb["step_reaction"]
    assert not sr.duplicated(["step_id", "reaction_id", "protein_accession", "pmid"]).any()
    for _, r in sr.iterrows():
        q = r["quote"]
        if q:
            assert r["pmid"], (r["step_id"], r["protein_accession"])
            assert len(q) <= 400 and q == q.strip() and q.endswith(".")
            assert not re.search(r"[.!?] [A-Z]", q[:-1].replace("N. punctiforme", "")), q


def test_links_are_the_expected_set(kb):
    assert _links(kb, "ddgs") == {
        ("RHEA:49560", "WEV88679.1", FUNGAL_PMID): "T1c",
        ("RHEA:49560", "B2J6X9", "21890703"): "T2",
        ("RHEA:49560", "Q3M6C3", "22741921"): "T3"}
    assert _links(kb, "omt") == {
        (C + "ddg_o_methylation", "WEV88680.1", FUNGAL_PMID): "T1a",
        (C + "ddg_o_methylation", "B2J6X8", "21890703"): "T3"}   # fix round 1 (I3)
    assert _links(kb, "atp_grasp_ligase") == {
        (C + "mycosporine_glutaminol_formation", "WEV88681.1", FUNGAL_PMID): "T1a",
        (C + "mycosporine_glutaminol_formation", "A0A9Y1GGR7", ""): "T3",
        (C + "mycosporine_glycine_formation", "B2J6X7", "21890703"): "T2"}
    assert _links(kb, "mg_glucosyltransferase") == {(C + "mg_glucosylation", "", ""): "T4"}
    assert _links(kb, "glutamicol_formation") == {
        (C + "mycosporine_glutamicol_formation", "", ""): "T4"}


def test_no_t1_link_without_a_fungal_protein_and_a_fungal_paper(kb):
    """Binding rule: a cyanobacterial or actinomycete enzyme is T2 at most. A T1 tier needs the
    fungal protein and the fungal paper."""
    sr = kb["step_reaction"]
    t1 = sr[sr["tier"].isin(["T1a", "T1b", "T1c"])]
    assert len(t1) == 3
    for _, r in t1.iterrows():
        assert r["organism"] == "Phaffia rhodozyma (strain UCD 67-385)", r["protein_accession"]
        assert species_key(r["organism"]) == PHAFFIA
        assert r["pmid"] == FUNGAL_PMID
        assert r["protein_accession"] in {"WEV88679.1", "WEV88680.1", "WEV88681.1"}
    t2 = sr[sr["tier"] == "T2"]
    assert set(t2["protein_accession"]) == {"B2J6X9", "B2J6X7"}
    for _, r in t2.iterrows():
        assert species_key(r["organism"]) == "nostoc punctiforme"
        assert r["pmid"] == "21890703" and r["quote"] == S_MYS
    bact = sr[sr["organism"].str.startswith(("Nostoc", "Trichormus"))]
    assert dict(bact["tier"].value_counts()) == {"T2": 2, "T3": 2}


def test_no_step_has_only_bacterial_t2_references_presented_as_fungal(kb):
    sr = kb["step_reaction"]
    for sid in ("ddgs", "omt", "atp_grasp_ligase"):
        own = sr[sr["step_id"] == sid]
        fungal = own[own["organism"].str.startswith("Phaffia rhodozyma")]
        assert set(fungal["tier"]) & {"T1a", "T1b", "T1c"}, sid
    # the bacterial ATP-grasp link is on another reaction (glycine, not glutaminol)
    gly = sr[sr["reaction_id"] == C + "mycosporine_glycine_formation"]
    assert set(gly["protein_accession"]) == {"B2J6X7"}
    glu = sr[sr["reaction_id"] == C + "mycosporine_glutaminol_formation"]
    assert set(glu["organism"].map(species_key)) == {PHAFFIA}
    # the second-amine reactions of the bacterial route have no step link at all
    for name in ("shinorine_formation", "porphyra334_formation"):
        assert (sr["reaction_id"] == C + name).sum() == 0, name
        note = kb["reaction"].set_index("reaction_id").loc[C + name, "note"]
        assert "no fungal enzyme" in note and "no step" in note, name


# Binding rule 1: a T1 or T2 link needs a quoted sentence that names the gene of its protein.
NAMED = {"WEV88679.1": ("Δddgs", S_DEL), "WEV88680.1": ("OMT", S_HET),
         "WEV88681.1": ("ATPG", S_HET), "B2J6X9": ("NpR5600", S_MYS),
         "B2J6X7": ("NpR5598", S_MYS)}


def test_every_t1_t2_quote_names_the_protein_of_its_link(kb):
    sr = kb["step_reaction"]
    ev = sr[sr["tier"].isin(EVIDENCE_TIERS)]
    assert set(ev["protein_accession"]) == set(NAMED)
    for _, r in ev.iterrows():
        name, sentence = NAMED[r["protein_accession"]]
        assert r["quote"] == sentence and name in r["quote"], r["protein_accession"]
    t1a = sr[sr["tier"] == "T1a"].set_index("protein_accession")
    # T1a: the sentence names the protein and the product of its reaction
    assert "4-deoxygadusol" in t1a.loc["WEV88680.1", "quote"]
    assert "MG production" in t1a.loc["WEV88681.1", "quote"]
    for a in t1a["assay"]:
        assert "retention time and absorption spectr" in a and "co-expression" in a
    # DDGS: the sentence states no product of its own reaction, so the link is not T1a
    ddgs = sr[sr["protein_accession"] == "WEV88679.1"].iloc[0]
    assert ddgs["tier"] == "T1c" and "demethyl-4-deoxygadusol is not named" in ddgs["assay"]
    assert "USER DECISION" in ddgs["assay"]


def test_held_links_say_why(kb):
    sr = kb["step_reaction"]
    t3 = sr[sr["tier"] == "T3"].set_index("protein_accession")
    assert set(t3.index) == {"Q3M6C3", "A0A9Y1GGR7", "B2J6X8"}
    # mysB: the gene is named, the product of the methylation is not, and no record supports it
    assert t3.loc["B2J6X8", "assay"].startswith("Held at T3.") and t3.loc["B2J6X8", "quote"] == S_MYS
    assert "no database record supports the activity" in t3.loc["B2J6X8", "assay"]
    # mysA stays T2 on the experimental evidence of the reviewed UniProt entry
    mysa = sr[sr["protein_accession"] == "B2J6X9"].iloc[0]
    assert mysa["tier"] == "T2" and "ECO:0000269" in mysa["assay"] and "20813918" in mysa["assay"]
    # OMT T1a and DDGS T1c rest on one experiment (I2)
    omt = sr[sr["protein_accession"] == "WEV88680.1"].iloc[0]
    assert "no strain with DDGS alone" in omt["assay"] and "USER DECISION 3" in omt["assay"]
    assert t3.loc["Q3M6C3", "assay"].startswith("Held at T3.")
    assert "does not name" in t3.loc["Q3M6C3", "assay"]
    assert "desmethyl-4-deoxygadusol synthases (DDGS)" in t3.loc["Q3M6C3", "quote"]
    assert t3.loc["A0A9Y1GGR7", "assay"].startswith("Inferred.")
    assert t3.loc["A0A9Y1GGR7", "quote"] == ""
    t4 = sr[sr["tier"] == "T4"]
    assert set(t4["step_id"]) == T4_STEPS
    assert set(t4["organism"]) == {NO_ENZYME} and set(t4["gene_name"]) == {"unidentified"}
    for a in t4["assay"]:
        assert a.startswith("hypothesis")
    # the glucosyltransferase activity of Ascochyta fabae is recorded; the link stays T4 (no sequence)
    gt = t4[t4["step_id"] == "mg_glucosyltransferase"].iloc[0]["assay"]
    for word in ("24258137", "Ascochyta fabae", "UDP-glucose", "no sequence is known", "UNVERIFIED"):
        assert word in gt, word


def test_every_link_names_organism_gene_and_accession_of_a_yaml_ref(kb):
    refs = {(st.id, r.accession): r for st in load_pathway(YAML).steps for r in st.refs}
    sr = kb["step_reaction"]
    real = sr[~sr["step_id"].isin(T4_STEPS)]
    for _, r in real.iterrows():
        ref = refs[(r["step_id"], r["protein_accession"])]
        assert r["organism"] == ref.species and r["gene_name"] == ref.gene_name
    assert {(s, a) for s, a in zip(real["step_id"], real["protein_accession"])} == set(refs)
    for st in load_pathway(YAML).steps:
        assert set(st.reactions) == set(sr[sr["step_id"] == st.id]["reaction_id"]), st.id
    assert set(kb["step"]["step_id"]) == {s.id for s in load_pathway(YAML).steps}


# ---- compounds and reactions ----

def test_compounds_and_kegg_cross_references(kb):
    c = kb["compound"].set_index("chebi_id")
    want = {SH7P: ("C7H15O10P", "C05382"), DDG: ("C7H10O5", "C21281"),
            MG: ("C13H22N2O6", ""), MGLUOH: ("C13H21NO7", ""), MSER: ("C11H19NO6", "C10607"),
            GAD: ("C8H12O6", "C22846"), SHI: ("C13H20N2O8", ""), P334: ("C14H22N2O8", ""),
            PAL: ("C10H16N2O5", "")}
    assert set(c.index) == set(want)
    for cid, (formula, kegg) in want.items():
        assert c.loc[cid, "formula"] == formula and c.loc[cid, "kegg_id"] == kegg, cid
    for _, r in kb["compound"].iterrows():
        assert r["neutral_exact_mass"] == f"{formula_mass(r['formula']):.4f}", r["chebi_id"]
        assert r["inchikey_block1"] == r["inchikey"][:14] and r["inchikey"].strip()
        assert not re.search(r"[+-]", r["formula"])
        assert "+" not in r["smiles"] and "-]" not in r["smiles"], r["chebi_id"]
    # the charged Rhea species and the anion-only entries are not in the table
    for charged in ("57483", "132139", "131711", "228248", "228249"):
        assert f"CHEBI:{charged}" not in c.index
    assert (kb["reaction"]["kegg_rn"] == "").all()
    ko = dict(zip(kb["step"]["step_id"], kb["step"]["ko"]))
    assert ko == {"ddgs": "K21342", "omt": "", "atp_grasp_ligase": "",
                  "mg_glucosyltransferase": "", "glutamicol_formation": ""}


def test_reactions_and_participants(kb):
    rx = kb["reaction"].set_index("reaction_id")
    assert set(rx.index) == {"RHEA:49560", C + "ddg_o_methylation",
                             C + "mycosporine_glutaminol_formation",
                             C + "mycosporine_glycine_formation", C + "mg_glucosylation",
                             C + "mycosporine_glutamicol_formation", C + "shinorine_formation",
                             C + "porphyra334_formation"}
    for rid, r in rx.iterrows():
        if r["curated"] == "true":
            assert rid.startswith(C) and r["rhea_id"] == ""
        else:
            assert rid == r["rhea_id"]
    assert rx.loc["RHEA:49560", "ec"] == "4.2.3.154"
    assert _parts(kb, "RHEA:49560", "substrate") == {SH7P}
    assert _parts(kb, "RHEA:49560", "product") == {DDG}
    # 4-deoxygadusol has no ChEBI entry: it is not listed, and the notes say so
    assert _parts(kb, C + "ddg_o_methylation", "substrate") == {DDG}
    assert _parts(kb, C + "ddg_o_methylation", "product") == set()
    assert _parts(kb, C + "mycosporine_glutaminol_formation", "substrate") == set()
    assert _parts(kb, C + "mycosporine_glutaminol_formation", "product") == {MG}
    assert _parts(kb, C + "mg_glucosylation", "substrate") == {MG}
    assert _parts(kb, C + "mg_glucosylation", "product") == set()
    assert _parts(kb, C + "mycosporine_glutamicol_formation", "product") == {MGLUOH}
    assert _parts(kb, C + "shinorine_formation", "product") == {SHI}
    assert _parts(kb, C + "porphyra334_formation", "product") == {P334}
    for name in ("ddg_o_methylation", "mycosporine_glutaminol_formation",
                 "mycosporine_glycine_formation"):
        assert "4-deoxygadusol has no ChEBI" in rx.loc[C + name, "note"], name
    assert "Rhea has no reaction" in rx.loc[C + "ddg_o_methylation", "note"]
    # context compounds that no reaction of the table makes
    p = kb["reaction_participant"]
    assert not {MSER, GAD, PAL} & set(p["compound"])


# ---- routes ----

def test_declared_routes_give_no_product_because_4_deoxygadusol_has_no_chebi_id(kb):
    for rid, steps in _routes().items():
        assert path_possible(kb, steps).empty, rid


def _one(got):
    assert len(got) == 1
    r = got.iloc[0]
    return (r["compound"], r["min_tier"], r["organisms"], bool(r["single_organism"]),
            int(r["n_equal_paths"]))


def test_route_prefixes_and_single_steps(kb):
    """Exact product, tier, organisms, single_organism and n_equal_paths."""
    assert _one(path_possible(kb, ["ddgs"])) == (DDG, "T1c", PHAFFIA, True, 1)
    assert path_possible(kb, ["ddgs", "omt"]).empty
    # one step alone: the first step of a call applies every link with no substrate requirement
    assert _one(path_possible(kb, ["atp_grasp_ligase"])) == (MG, "T1a", PHAFFIA, True, 1)
    assert path_possible(kb, ["mg_glucosyltransferase"]).empty      # T4 links are never used
    assert path_possible(kb, ["glutamicol_formation"]).empty


def test_routes_one_organism_at_a_time(kb):
    assert _one(path_possible(kb, ["ddgs"], organism=PHAFFIA)) == (DDG, "T1c", PHAFFIA, True, 1)
    assert _one(path_possible(kb, ["ddgs"], organism="Nostoc punctiforme")) == (
        DDG, "T2", "nostoc punctiforme", True, 1)
    assert _one(path_possible(kb, ["ddgs"], organism="Trichormus variabilis")) == (
        DDG, "T3", "trichormus variabilis", True, 1)
    assert path_possible(kb, ["atp_grasp_ligase"], organism="Nostoc punctiforme").empty
    orgs = {species_key(o) for o in kb["step_reaction"]["organism"] if o != NO_ENZYME}
    assert orgs == {PHAFFIA, "nostoc punctiforme", "trichormus variabilis"}


def test_species_keys_are_checked_against_ncbi_taxids():
    taxa = pd.read_csv(TAXA, sep="\t", dtype=str, keep_default_na=False).set_index("name")
    assert taxa.loc["Phaffia rhodozyma", "taxid"] == "264483"
    assert species_key("Phaffia rhodozyma (strain UCD 67-385)") == PHAFFIA
    assert SPECIES_ALIASES["phaffia rhodozyma"][1] == 264483
    assert taxa.loc["Nostoc punctiforme", "taxid"] == "272131"
    # one taxid for both names: the paper name Anabaena variabilis is an alias of the UniProt name
    assert taxa.loc["Trichormus variabilis", "taxid"] == taxa.loc["Anabaena variabilis", "taxid"]
    assert SPECIES_ALIASES["anabaena variabilis"][:2] == ("trichormus variabilis", 264691)
    assert species_key("Anabaena variabilis ATCC 29413") == "trichormus variabilis"
    assert species_key("Trichormus variabilis (strain ATCC 29413 / PCC 7937)") == (
        "trichormus variabilis")


# ---- change detector and README ----

def test_quote_fixture_equals_the_table(kb):
    fx = pd.read_csv(QUOTES, sep="\t", dtype=str, keep_default_na=False)
    assert list(fx.columns) == ["step_id", "reaction_id", "protein_accession", "pmid", "quote",
                                "tier", "assay"]
    sr = kb["step_reaction"]
    got = sr[["step_id", "reaction_id", "protein_accession", "pmid", "quote", "tier", "assay"]]
    assert len(fx) == len(sr) == 10
    pd.testing.assert_frame_equal(fx.reset_index(drop=True), got.reset_index(drop=True))


def test_readme_documents_the_mycosporine_conventions():
    text = (ROOT / "kb" / "README.md").read_text(encoding="utf-8")
    assert "`mycosporine`" in text and "kb/mycosporine/" in text
    sec = text.split("## Conventions used by the mycosporine tables")[1]
    for word in ("4-deoxygadusol", "no ChEBI", "T2", "cyanobacterial", "GenBank", "WEV88679.1",
                 "ordered locus name", "fungal_evidence", "mycosporine_quotes.tsv",
                 "fungal_lineages.tsv", "refs: []"):
        assert word in sec, word


def test_fasta_sequences_equal_the_fetched_records():
    """tests/data/kb/mycosporine_sequences.tsv holds the length and SHA-256 of each sequence as
    fetched from UniProt or GenBank by build.py. Every FASTA record must equal it."""
    fx = pd.read_csv(ROOT / "tests" / "data" / "kb" / "mycosporine_sequences.tsv", sep="\t",
                     dtype=str, keep_default_na=False).set_index("accession")
    assert len(fx) == 24 and set(fx["source"]) == {"UniProtKB", "GenBank"}
    seen = set()
    refs = ROOT / "refs" / "mycosporine"
    for path in sorted(refs.glob("*.faa")) + sorted((refs / "decoys").glob("*.faa")):
        name, seqs = None, {}
        for ln in path.read_text(encoding="utf-8").splitlines():
            if ln.startswith(">"):
                name = ln[1:].split()[0]
                seqs[name] = ""
            else:
                seqs[name] += ln.strip()
        for acc, seq in seqs.items():
            assert str(len(seq)) == fx.loc[acc, "length"], (path.name, acc)
            assert hashlib.sha256(seq.encode()).hexdigest() == fx.loc[acc, "sha256"], (path.name, acc)
            seen.add(acc)
    assert seen == set(fx.index)
