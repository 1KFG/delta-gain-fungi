"""Tests for bin/make_review_packet.py (stage 1 review packet).

The packet is generated from the YAMLs, the kb tables, refs/<pathway>/ and the committed
verification report. These tests check that it is deterministic, complete, and that every number
in its generated text has a source. No network.
"""
import importlib.util
import re
from pathlib import Path

import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
spec = importlib.util.spec_from_file_location("make_review_packet", ROOT / "bin" / "make_review_packet.py")
mrp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mrp)

PATHWAYS = ["carotenoid", "melanin_dhn", "melanin_ldopa", "mycosporine"]
MARKERS = ("USER DECISION", "DECISION FOR THE USER")
SPEC = REPO / "docs" / "superpowers" / "specs" / "2026-11-21-pathway-knowledge-graph-design.md"
REPORT = ROOT / "kb" / "VERIFICATION_REPORT.md"
PACKET = REPO / "docs" / "superpowers" / "plans" / "2026-11-21-stage1-review-packet.md"


@pytest.fixture(scope="module")
def built():
    return mrp.build_packet(ROOT)


@pytest.fixture(scope="module")
def packet(built):
    return built[0]


def test_deterministic_byte_for_byte(built):
    again = mrp.build_packet(ROOT)
    assert again[0].encode("utf-8") == built[0].encode("utf-8")
    assert again[1] == built[1]


def test_committed_packet_is_current(packet):
    assert PACKET.read_text(encoding="utf-8") == packet


@pytest.mark.parametrize("name", PATHWAYS)
def test_every_step_id_is_in_the_packet(packet, name):
    d = yaml.safe_load((ROOT / "pathways" / f"{name}.yaml").read_text(encoding="utf-8"))
    for st in d["steps"]:
        assert f"`{st['id']}`" in packet, (name, st["id"])
    kb_steps = pd.read_csv(ROOT / "kb" / name / "step.tsv", sep="\t", dtype=str)["step_id"]
    for sid in kb_steps:
        assert f"`{sid}`" in packet, (name, sid)
    for r in d["routes"]:
        assert f"`{r['id']}`" in packet, (name, r["id"])


def test_verification_counts_equal_report(packet):
    rep = REPORT.read_text(encoding="utf-8")
    sec = rep.split("## Counts per pathway")[1].split("## Counts per check")[0]
    rows = [ln for ln in sec.splitlines() if re.match(r"\| (carotenoid|melanin_dhn|melanin_ldopa|mycosporine|all) \|", ln)]
    assert len(rows) == 5
    for ln in rows:
        assert ln in packet.splitlines(), ln
    sec = rep.split("## Counts per check")[1].split("## Gene-name passes by kind")[0]
    chk = [ln for ln in sec.splitlines() if re.match(r"\| [a-z_0-9]+ \| \d+ \| \d+ \| \d+ \|$", ln)]
    assert len(chk) == 23
    for ln in chk:
        assert ln in packet.splitlines(), ln


def test_every_failure_row_of_the_report_is_in_the_packet(packet):
    rep = REPORT.read_text(encoding="utf-8")
    sec = rep.split("## Failures")[1].split("## Unchecked")[0]
    rows = [ln for ln in sec.splitlines() if re.match(r"\| (carotenoid|melanin_dhn|melanin_ldopa|mycosporine) \|", ln)]
    assert len(rows) == 6
    for ln in rows:
        assert ln in packet.splitlines(), ln


def test_unchecked_summary_counts_sum_to_the_report(built, packet):
    rep = REPORT.read_text(encoding="utf-8")
    sec = rep.split("## Unchecked")[1]
    n = len([ln for ln in sec.splitlines() if re.match(r"\| (carotenoid|melanin_dhn|melanin_ldopa|mycosporine) \|", ln)])
    assert n == 24
    rows = [ln for ln in packet.split("### Unchecked, by check and reason")[1].split("\n## ")[0].splitlines()
            if re.match(r"\| (carotenoid|melanin_dhn|melanin_ldopa|mycosporine) \|", ln)]
    total = sum(int(ln.strip("|").split("|")[-1]) for ln in rows)
    assert total == n


@pytest.mark.parametrize("name", PATHWAYS)
def test_link_counts_per_tier(packet, name):
    sr = pd.read_csv(ROOT / "kb" / name / "step_reaction.tsv", sep="\t", dtype=str)
    for tier, n in sr["tier"].value_counts().items():
        assert f"| {name} | {tier} | {n} |" in packet, (name, tier, n)


@pytest.mark.parametrize("name", PATHWAYS)
def test_route_rows_equal_path_possible(packet, name):
    from pathway_census.kb import load_kb, path_possible
    kb = load_kb(ROOT / "kb" / name)
    d = yaml.safe_load((ROOT / "pathways" / f"{name}.yaml").read_text(encoding="utf-8"))
    for r in d["routes"]:
        try:
            df = path_possible(kb, r["step_ids"])
        except ValueError:
            continue
        for row in df.itertuples():
            cell = f"| `{r['id']}` | {row.compound} |"
            line = next(ln for ln in packet.splitlines() if ln.startswith(cell))
            cols = [c.strip() for c in line.strip("|").split("|")]
            assert row.min_tier == cols[3]
            assert row.organisms == cols[4].strip("`")
            assert str(bool(row.single_organism)) == cols[5]
            assert str(row.n_equal_paths) == cols[6]


def test_tier_definitions_match_spec(packet):
    text = SPEC.read_text(encoding="utf-8")
    sec = text.split("### 6.2 ")[1].split("### 6.3")[0].split("\n", 1)[1]
    lines = [ln for ln in sec.splitlines() if ln.strip()]
    assert len(lines) >= 7
    for ln in lines:
        assert ln in packet, ln


@pytest.mark.parametrize("name", PATHWAYS)
def test_every_decision_item_is_verbatim_in_the_packet(packet, name):
    text = (ROOT / "refs" / name / "UNVERIFIED.md").read_text(encoding="utf-8")
    body = "\n## " + text.split("\n## ", 1)[1]  # drop the header bullets before the first section
    lines = [ln for ln in body.splitlines() if any(m in ln for m in MARKERS)]
    assert lines, name
    for ln in lines:
        assert ln in packet, (name, ln[:80])
    ids = re.findall(rf"\*\*D{mrp.LETTER[name]}(\d+)\*\*", packet)
    assert ids and ids == [str(i) for i in range(1, len(ids) + 1)], (name, ids)


def test_decision_blocks_are_in_verbatim_divs(packet):
    divs = re.findall(mrp.VERBATIM_RE, packet, flags=re.S | re.M)
    assert divs
    text = (ROOT / "refs" / "melanin_dhn" / "UNVERIFIED.md").read_text(encoding="utf-8")
    sec = text.split("### G. USER DECISION")[1]
    assert sec.split("\n", 1)[1].strip() in packet


def test_rule1_reading_is_at_the_top(packet):
    top = packet.split("## How to read this packet")[0]
    assert "strict" in top and "context" in top
    assert "not recomputed" in top
    assert "carotenoid" in top and "melanin_dhn" in top


def test_numbers_in_text_are_sourced(built):
    text, numbers = built
    plain = mrp.strip_verbatim(text)
    plain = re.sub(r"\A---\n.*?\n---\n", "", plain, flags=re.S)           # front matter
    plain = re.sub(r"\b(PMID|PMIDs)\s+\d+(,?\s+\d+)*", "", plain)         # PMIDs
    plain = re.sub(r"\b(CHEBI|RHEA|CID|taxid):\s?\d+", "", plain)         # ChEBI, Rhea ids
    plain = re.sub(r"\b\d{4}-\d{2}-\d{2}\b", "", plain)                   # dates
    plain = re.sub(r"\*\*D[CDLM]\d+\*\*", "", plain)                      # decision numbering
    plain = re.sub(r"\bD[CDLM]\d+\b", "", plain)                          # decision references
    toks = set()
    for tok in re.split(r"[\s|`(),;:\[\]*_/=<>]+", plain):
        tok = tok.strip(".'\"")
        if re.fullmatch(r"\d+(\.\d+)?", tok):
            toks.add(tok)
    assert toks, "no numbers found: the test is not testing anything"
    missing = sorted(toks - set(numbers))
    assert not missing, missing
    assert len(numbers) > 50
