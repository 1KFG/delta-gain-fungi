"""Tests for the DHN-melanin flat knowledge tables (stage 1, task S4)."""
import re
from pathlib import Path

import pandas as pd
import pytest

from pathway_census.defs import load_pathway
from pathway_census.kb import (EVIDENCE_TIERS, formula_mass, load_kb, path_possible,
                               species_key, validate_kb)

ROOT = Path(__file__).resolve().parent.parent
KB = ROOT / "kb" / "melanin_dhn"
YAML = ROOT / "pathways" / "melanin_dhn.yaml"
QUOTES = ROOT / "tests" / "data" / "kb" / "melanin_dhn_quotes.tsv"
C = "CUR:melanin_dhn:"
SCYTALONE, T4HN, T3HN = "CHEBI:16945", "CHEBI:18365", "CHEBI:18393"
VERMELONE, DHN, FLAVIOLIN, YWA1 = "CHEBI:9957", "CHEBI:167604", "CHEBI:42646", "CHEBI:133763"


@pytest.fixture(scope="module")
def kb():
    return load_kb(KB)


def _parts(kb, rid, role):
    p = kb["reaction_participant"]
    return set(p[(p["reaction_id"] == rid) & (p["role"] == role)]["compound"])


def _routes():
    return {r.id: r.step_ids for r in load_pathway(YAML).routes}


def test_tables_validate(kb):
    assert validate_kb(kb) == []


def test_tables_not_empty(kb):
    for name, df in kb.items():
        assert len(df) > 0, name


def test_every_t1_t2_link_has_pmid_quote_and_assay(kb):
    sr = kb["step_reaction"]
    ev = sr[sr["tier"].isin(EVIDENCE_TIERS)]
    assert len(ev) > 0
    for _, r in ev.iterrows():
        assert r["pmid"].strip().isdigit(), (r["step_id"], r["reaction_id"])
        assert r["quote"].strip(), (r["step_id"], r["reaction_id"])
        assert r["assay"].strip(), (r["step_id"], r["reaction_id"])


def test_no_t4_and_no_t2_links(kb):
    assert set(kb["step_reaction"]["tier"]) <= {"T1a", "T1b", "T1c", "T3"}


def test_one_quoted_sentence_per_link(kb):
    """The plan's rule: at most one quoted sentence per link. A paper may give different
    sentences to different links when the links rest on different statements."""
    sr = kb["step_reaction"]
    assert not sr.duplicated(["step_id", "reaction_id", "protein_accession", "pmid"]).any()
    for _, r in sr.iterrows():
        q = r["quote"]
        if q:
            assert r["pmid"], (r["step_id"], r["protein_accession"])
            assert len(q) <= 400 and q.endswith(".")
            assert not re.search(r"[.!?] [A-Z]", q[:-1]), (r["pmid"], q)
    # PMID 10515939 gives two sentences: intermediates (arp1, arp2) and colour (abr1)
    by = sr[sr["pmid"] == "10515939"]
    quotes = dict(zip(by["protein_accession"], by["quote"]))
    assert quotes["E9QUT3"] == quotes["O14434"] != quotes["Q4WZB4"]
    assert "scytalone and flaviolin in the arp1 deletant" in quotes["O14434"]
    assert quotes["Q4WZB4"] == (
        "Involvement of the six genes in conidial pigmentation was confirmed by the altered "
        "conidial color phenotypes that resulted from disruption of each gene in A. fumigatus.")


def test_every_link_names_organism_gene_and_accession_of_a_yaml_ref(kb):
    refs = {}
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            refs[(st.id, r.accession)] = r
    for _, r in kb["step_reaction"].iterrows():
        ref = refs[(r["step_id"], r["protein_accession"])]
        assert r["organism"] == ref.species, (r["step_id"], r["protein_accession"])
        assert r["gene_name"] == ref.gene_name, (r["step_id"], r["protein_accession"])
    # every reference has at least one link
    assert {(s, a) for s, a in zip(kb["step_reaction"]["step_id"],
                                   kb["step_reaction"]["protein_accession"])} == set(refs)


def test_kb_steps_match_yaml_steps(kb):
    assert set(kb["step"]["step_id"]) == {s.id for s in load_pathway(YAML).steps}


def test_yaml_step_reactions_are_the_kb_reactions(kb):
    sr = kb["step_reaction"]
    for st in load_pathway(YAML).steps:
        assert set(st.reactions) == set(sr[sr["step_id"] == st.id]["reaction_id"]), st.id


def test_compounds_of_requirement_7_with_kegg_cross_references(kb):
    c = kb["compound"].set_index("chebi_id")
    want = {SCYTALONE: "C00779", T4HN: "C04033", T3HN: "C01173", VERMELONE: "C01624",
            DHN: "", FLAVIOLIN: "C18012", YWA1: ""}
    for cid, kegg in want.items():
        assert c.loc[cid, "kegg_id"] == kegg, cid
    assert c.loc[SCYTALONE, "formula"] == "C10H10O4"
    assert c.loc[T4HN, "formula"] == "C10H8O4"
    assert c.loc[YWA1, "formula"] == "C14H12O6"


def test_masses_are_computed_from_formulas(kb):
    for _, r in kb["compound"].iterrows():
        assert r["neutral_exact_mass"] == f"{formula_mass(r['formula']):.4f}", r["chebi_id"]
        assert r["inchikey"].strip() and r["inchikey_block1"] == r["inchikey"][:14]


def test_no_kegg_reaction_content_and_no_bacterial_rppa(kb):
    assert set(kb["step"]["ko"]) <= {"", "K17739", "K17740"}
    ko = dict(zip(kb["step"]["step_id"], kb["step"]["ko"]))
    assert ko["t4hn_reductase"] == "K17739" and ko["scytalone_dehydratase"] == "K17740"
    assert (kb["reaction"]["kegg_rn"] == "").all()
    for col in ("reaction_id", "rhea_id", "ec"):
        assert not kb["reaction"][col].str.contains("41524|2.3.1.233").any()
    for name in ("step", "reaction"):
        for col in kb[name].columns:
            if col != "note":
                assert not kb[name][col].str.contains("K19580").any(), (name, col)


def test_curated_reaction_ids(kb):
    for _, r in kb["reaction"].iterrows():
        if r["curated"] == "true":
            assert r["reaction_id"].startswith(C), r["reaction_id"]
            assert r["rhea_id"] == ""
        else:
            assert r["reaction_id"] == r["rhea_id"] and r["rhea_id"].startswith("RHEA:")
    assert {"RHEA:62652", "RHEA:62700", "RHEA:21908", "RHEA:24396"} <= set(
        kb["reaction"]["reaction_id"])


def test_reductase_reactions_are_stored_in_the_direction_of_the_step(kb):
    assert _parts(kb, "RHEA:21908", "substrate") == {T4HN}
    assert _parts(kb, "RHEA:21908", "product") == {SCYTALONE}
    assert _parts(kb, C + "t3hn_reduction", "substrate") == {T3HN}
    assert _parts(kb, C + "t3hn_reduction", "product") == {VERMELONE}
    rx = kb["reaction"].set_index("reaction_id")
    assert rx.loc["RHEA:21908", "ec"] == "1.1.1.252"
    assert rx.loc[C + "t3hn_reduction", "ec"] == ""


def test_two_reductase_types_have_organism_bound_tiers(kb):
    sr = kb["step_reaction"]
    main4 = sr[(sr["step_id"] == "t4hn_reductase") & (sr["reaction_id"] == "RHEA:21908")]
    assert dict(zip(main4["protein_accession"], main4["tier"])) == {
        "Q9HFV6": "T1a", "E9QUT3": "T1b", "W3XAA0": "T1a", "F8WTI4": "T3", "A7DWM4": "T3"}
    main3 = sr[(sr["step_id"] == "t3hn_reductase") & (sr["reaction_id"] == C + "t3hn_reduction")]
    assert dict(zip(main3["protein_accession"], main3["tier"])) == {
        "Q12634": "T1a", "P87025": "T1b", "A7DWM6": "T3"}
    # Arp2 is not shown to reduce T3HN in the abstracts that were read
    arp2 = sr[(sr["protein_accession"] == "E9QUT3") & (sr["reaction_id"] == C + "t3hn_reduction")]
    assert list(arp2["tier"]) == ["T3"]
    # substrate overlap: the 3HNR type also reduces T4HN (two organisms)
    side = sr[(sr["step_id"] == "t3hn_reductase") & (sr["reaction_id"] == "RHEA:21908")]
    assert dict(zip(side["protein_accession"], side["pmid"])) == {
        "Q12634": "10956664", "P87025": "12967213"}


def test_vermelone_dehydration_is_a_separate_curated_reaction_with_its_own_evidence(kb):
    rid = C + "vermelone_dehydration"
    rx = kb["reaction"].set_index("reaction_id")
    assert rx.loc[rid, "curated"] == "true" and rx.loc[rid, "rhea_id"] == ""
    assert "unsettled" in rx.loc[rid, "note"] or "not settled" in rx.loc[rid, "note"]
    assert _parts(kb, rid, "substrate") == {VERMELONE}
    assert _parts(kb, rid, "product") == {DHN}
    assert _parts(kb, "RHEA:24396", "substrate") == {SCYTALONE}
    assert _parts(kb, "RHEA:24396", "product") == {T3HN}
    sr = kb["step_reaction"]
    links = sr[sr["reaction_id"] == rid]
    assert set(links["step_id"]) == {"scytalone_dehydratase"}
    # tier T3, or a paper whose abstract shows the reaction (with its quote)
    for _, r in links.iterrows():
        if r["tier"] != "T3":
            assert r["tier"] in EVIDENCE_TIERS and r["pmid"] and (
                "vermelone" in r["quote"].lower() or "veremelone" in r["quote"].lower()), \
                r["protein_accession"]
    tiers = dict(zip(links["protein_accession"], links["tier"]))
    assert tiers == {"Q00455": "T3", "O14434": "T3", "P56221": "T1a", "A7DWM5": "T3"}
    row = links[links["protein_accession"] == "P56221"].iloc[0]
    assert row["pmid"] == "10320327"
    # the evidence record is not a copy of the scytalone edge
    scy = sr[(sr["reaction_id"] == "RHEA:24396") & (sr["protein_accession"] == "P56221")].iloc[0]
    assert (scy["pmid"], scy["quote"]) != (row["pmid"], row["quote"])


def test_scytalone_dehydration_tiers(kb):
    sr = kb["step_reaction"]
    links = sr[sr["reaction_id"] == "RHEA:24396"]
    assert dict(zip(links["protein_accession"], links["tier"])) == {
        "Q00455": "T1c", "O14434": "T1b", "P56221": "T1a", "A7DWM5": "T3"}


def test_oxidase_links_are_phenotype_only(kb):
    sr = kb["step_reaction"]
    ox = sr[sr["step_id"].isin(["laccase_abr2", "mco_abr1"])]
    assert set(ox["reaction_id"]) == {C + "dhn_polymerisation"}
    # The quoted sentence of PMID 10515939 speaks of "the six genes" and does not name abr1,
    # so the abr1 link is held at T3 with that sentence as context (fix round 1).
    assert dict(zip(ox["protein_accession"], ox["tier"])) == {
        "E9RBR0": "T1c", "Q4WZB4": "T3", "A2R3T0": "T1c"}
    abr1 = ox[ox["protein_accession"] == "Q4WZB4"].iloc[0]
    assert "abr1" not in abr1["quote"] and "disruption of each gene" in abr1["quote"]
    assert "does not name abr1" in abr1["assay"]
    assert dict(zip(ox["protein_accession"], ox["pmid"])) == {
        "E9RBR0": "16988817", "Q4WZB4": "10515939", "A2R3T0": "21277986"}
    assert _parts(kb, C + "dhn_polymerisation", "substrate") == {DHN}
    assert _parts(kb, C + "dhn_polymerisation", "product") == set()


def test_pks_and_ayg1_tiers(kb):
    sr = kb["step_reaction"]
    got = {(s, r, a, p): t for s, r, a, p, t in zip(
        sr["step_id"], sr["reaction_id"], sr["protein_accession"], sr["pmid"], sr["tier"])
        if s in ("thn_synthase", "hexaketide_pks", "heptaketide_pks", "ayg1")}
    assert got == {
        ("thn_synthase", C + "t4hn_synthesis", "P79068", "10501004"): "T1a",
        ("thn_synthase", C + "t4hn_synthesis", "P79068", "23261597"): "T1a",
        ("hexaketide_pks", C + "at4hn_synthesis", "Q9Y7A7", "11179356"): "T1c",
        ("heptaketide_pks", "RHEA:62652", "Q4WZA8", "11040426"): "T1a",
        ("heptaketide_pks", "RHEA:62652", "A2QUI2", "21277986"): "T1c",
        ("ayg1", "RHEA:62700", "Q4WZB3", "15310761"): "T1a",
        ("ayg1", "RHEA:62700", "Q4WZB3", "11350964"): "T1a",
        ("ayg1", "RHEA:62700", "A2R3S8", "21277986"): "T1c",
        ("ayg1", C + "at4hn_deacetylation", "Q6BCB5", "18676950"): "T1b",
    }
    assert _parts(kb, C + "t4hn_synthesis", "product") == {T4HN}
    assert _parts(kb, "RHEA:62652", "product") == {YWA1}
    assert _parts(kb, "RHEA:62700", "substrate") == {YWA1}
    assert T4HN in _parts(kb, "RHEA:62700", "product")
    # 2-acetyl-1,3,6,8-tetrahydroxynaphthalene has no ChEBI entry: no participant for it
    assert _parts(kb, C + "at4hn_synthesis", "product") == set()
    assert _parts(kb, C + "at4hn_deacetylation", "substrate") == set()


def test_flaviolin_is_a_shunt_product_with_no_enzyme_link(kb):
    rid = C + "t4hn_autooxidation"
    assert _parts(kb, rid, "substrate") == {T4HN} and _parts(kb, rid, "product") == {FLAVIOLIN}
    assert rid not in set(kb["step_reaction"]["reaction_id"])


# Expected final products of every declared route (path_possible over the route's step order):
# compound -> (weakest tier along the best path, organisms of the links on that path, number of
# organism lists with that tier and that size; the list shown is the alphabetically first one).
# The tier of a route does NOT hold for one organism unless the organism list has one name.
# An empty dict means: the route cannot be followed in the tables (UNVERIFIED.md, section D).
AF, PF = "aspergillus fumigatus", "pestalotiopsis fici"
PO, PG = "pyricularia oryzae", "pyricularia grisea"
GL, CO = "glomerella lagenarium", "colletotrichum orbiculare"
ROUTE_PRODUCTS = {
    # 1,8-DHN. The second equal list has Pyricularia grisea in place of Pestalotiopsis fici.
    "pks_thn": {DHN: ("T1a", [GL, PF, PO], 2)},
    "pks_heptaketide": {DHN: ("T1a", [AF, PF, PO], 2)},
    # the hexaketide 2-acetyl-1,3,6,8-tetrahydroxynaphthalene has no ChEBI entry
    "pks_hexaketide": {},
    # the polymer is not a compound of the table
    "polymerisation": {},
}


def _have(got):
    return {c: (t, o.split(","), single, n) for c, t, o, single, n in
            zip(got["compound"], got["min_tier"], got["organisms"], got["single_organism"],
                got["n_equal_paths"])}


def test_every_declared_route_gives_its_expected_products(kb):
    routes = _routes()
    assert set(routes) == set(ROUTE_PRODUCTS)
    for rid, want in ROUTE_PRODUCTS.items():
        got = path_possible(kb, routes[rid])
        assert _have(got) == {c: (t, o, len(o) == 1, n) for c, (t, o, n) in want.items()}, rid
        # no route tier of the unfiltered call holds for one organism
        assert not got["single_organism"].any(), rid


def test_no_organism_completes_a_route_alone(kb):
    routes = _routes()
    organisms = sorted({species_key(o) for o in kb["step_reaction"]["organism"]})
    assert organisms == sorted([AF, "aspergillus niger", CO, "exophiala dermatitidis", GL, PF,
                                PG, PO, "sordaria macrospora"])
    for rid, steps in routes.items():
        for org in organisms:
            assert path_possible(kb, steps, organism=org).empty, (rid, org)


def test_how_far_single_organisms_get(kb):
    """Prefixes of the declared routes that one organism completes with its own links."""
    hepta = _routes()["pks_heptaketide"]
    # Aspergillus fumigatus: alb1, ayg1, arp2, arp1 reach 1,3,8-trihydroxynaphthalene at T1b.
    got = path_possible(kb, hepta[:4], organism="Aspergillus fumigatus")
    assert _have(got) == {T3HN: ("T1b", [AF], True, 1)}
    # The next step has no A. fumigatus link (Arp2 on T3HN is a T3 link of t4hn_reductase).
    assert path_possible(kb, hepta[:5], organism="Aspergillus fumigatus").empty
    # Tail from T4HN: Colletotrichum orbiculare (T4HR1, SCD1, THR1) and Sordaria macrospora
    # (teh, sdh, tih) reach 1,8-DHN at T3. Their PKS is not in the same species key.
    tail = hepta[2:]
    for org in (CO, "sordaria macrospora"):
        assert _have(path_possible(kb, tail, organism=org)) == {DHN: ("T3", [org], True, 1)}


def test_tail_of_the_hexaketide_route_reaches_dhn(kb):
    """Route pks_hexaketide stops at the hexaketide because it has no ChEBI entry. From Yg1
    onwards the steps reach 1,8-DHN (first step: no requirement on substrates)."""
    got = path_possible(kb, _routes()["pks_hexaketide"][1:])
    assert set(got["compound"]) == {DHN}


def test_organism_filter_accepts_the_names_used_in_the_papers(kb):
    assert species_key("Wangiella dermatitidis") == species_key("Exophiala dermatitidis")
    assert species_key("Magnaporthe oryzae") == species_key(
        "Pyricularia oryzae (strain 70-15 / ATCC MYA-4617 / FGSC 8958)") == PO
    assert species_key("Magnaporthe grisea") == PG
    # two NCBI taxa: never merged
    assert species_key("Pyricularia grisea") != species_key("Pyricularia oryzae")
    assert species_key("Glomerella lagenarium") != species_key("Colletotrichum orbiculare")


def test_quotes_are_single_sentences_of_limited_length(kb):
    for _, r in kb["step_reaction"].iterrows():
        q = r["quote"]
        if not q:
            continue
        assert len(q) <= 400, (r["pmid"], len(q))
        assert q == q.strip() and q.endswith("."), r["pmid"]
        assert not re.search(r"[.!?] [A-Z]", q[:-1]), r["pmid"]


def test_quotes_equal_the_committed_fixture(kb):
    """Change detector: any edit of a quote (or of the links that carry one) fails here.
    The fixture is not a check against PubMed; Task 7 does that."""
    cols = ["step_id", "reaction_id", "protein_accession", "tier", "assay", "pmid", "quote"]
    fix = pd.read_csv(QUOTES, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8",
                      quoting=3)
    assert list(fix.columns) == cols
    cur = kb["step_reaction"][cols].reset_index(drop=True)
    assert len(cur) > 0 and (fix["quote"] != "").sum() >= 20
    assert cur.to_dict("records") == fix.to_dict("records")


def test_no_kegg_pubchem_column(kb):
    assert list(kb["compound"].columns) == ["chebi_id", "inchikey", "inchikey_block1", "smiles",
                                            "formula", "neutral_exact_mass", "kegg_id",
                                            "np_class"]
