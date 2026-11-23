"""Tests for reports/carotenoid/figures/bin/fig_pathway_diagram.py."""
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "fig_pathway_diagram", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_pathway_diagram.py")
fp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fp)


def write_kb(d):
    d.mkdir(parents=True, exist_ok=True)
    (d / "reaction.tsv").write_text("reaction_id\tec\tkegg_rn\tdirection\tcurated\tnote\trhea_id\n"
                                    "R1\t2.5.1.29\t\tforward\tfalse\t\t\nR2\t1.1.1.1\t\tforward\tfalse\t\t\nR3\t\t\tforward\tfalse\t\t\n")
    (d / "reaction_participant.tsv").write_text(
        "reaction_id\tcompound\trole\tstoich\n"
        "R1\tCHEBI:16584\tsubstrate\t1\nR1\tCHEBI:17407\tsubstrate\t1\nR1\tCHEBI:48861\tproduct\t1\n"
        "R2\tCHEBI:48861\tsubstrate\t1\nR2\tCHEBI:27787\tproduct\t1\nR2\tCHEBI:15825\tproduct\t1\n")
    (d / "step_reaction.tsv").write_text(
        "step_id\treaction_id\torganism\tgene_name\tprotein_accession\ttier\tassay\tpmid\tquote\n"
        "s1\tR1\to\tg\tP1\tT3\t\t\t\ns2\tR2\to\tg\tP2\tT3\t\t\t\ns2\tR2\to\tg\tP3\tT1a\t\t\t\ns3\tR3\to\tg\tP4\tT3\t\t\t\n")


def test_best_tier_picks_the_strongest_and_ignores_unknown():
    assert fp.best_tier(["T3", "T1b", "T4"]) == "T1b"
    assert fp.best_tier(["X"]) is None and fp.best_tier([]) is None


def test_build_graph_edges_feeder_node_side_product_and_undrawn(tmp_path):
    write_kb(tmp_path / "kb")
    edges, undrawn = fp.build_graph(tmp_path / "kb")
    got = {(e["src"], e["dst"], e["tier"], tuple(e["steps"])) for e in edges}
    assert ("FEED", "CHEBI:48861", "T3", ("s1",)) in got               # IPP and FPP merged into one grey node
    assert ("CHEBI:48861", "CHEBI:27787", "T1a", ("s2",)) in got       # best tier of the reaction, not the weakest
    assert not any(e["dst"] == "CHEBI:15825" for e in edges)           # side product is not an edge end
    assert undrawn == [("R3", ["s3"], "no substrate or no product compound in the tables")]


def test_compound_name_strips_html_and_falls_back_to_id(tmp_path):
    c = tmp_path / "chebi"
    c.mkdir()
    (c / "27787.json").write_text(json.dumps({"data": {"name": "15-<i>cis</i>-phytoene"}}))
    assert fp.compound_name("CHEBI:27787", tmp_path) == "15-cis-phytoene"
    assert fp.compound_name("CHEBI:1", tmp_path) == "CHEBI:1"


def test_dot_text_has_edge_labels_styles_legend_and_note(tmp_path):
    write_kb(tmp_path / "kb")
    edges, undrawn = fp.build_graph(tmp_path / "kb")
    txt = fp.dot_text(edges, undrawn, tmp_path, "T")
    assert "digraph G" in txt and 's2\\nT1a' in txt and 'style=solid' in txt and "legend" in txt and "Steps with no drawable edge" in txt


@pytest.mark.skipif(shutil.which("dot") is None, reason="graphviz dot not on PATH")
def test_main_writes_png_pdf_svg(tmp_path):
    write_kb(tmp_path / "kb")
    out = tmp_path / "out"
    assert fp.main(["--kb-dir", str(tmp_path / "kb"), "--cache-dir", str(tmp_path), "--out-dir", str(out), "--name", "x"]) == 0
    assert (out / "x.pdf").read_bytes().startswith(b"%PDF")
    assert (out / "x.png").stat().st_size > 500 and (out / "x.svg").read_text().lstrip().startswith("<?xml")


def test_short_name_strips_stereo_prefix_and_names_ggpp():
    assert fp.short_name("2-trans,6-trans,10-trans-geranylgeranyl diphosphate") == "geranylgeranyl diphosphate (GGPP)"
    assert fp.short_name("all-trans-retinal") == "all-trans-retinal"


def test_dot_text_rankdir_is_set(tmp_path):
    write_kb(tmp_path / "kb")
    edges, undrawn = fp.build_graph(tmp_path / "kb")
    assert "rankdir=TB" in fp.dot_text(edges, undrawn, tmp_path, "T") and "rankdir=LR" in fp.dot_text(edges, undrawn, tmp_path, "T", "LR")
