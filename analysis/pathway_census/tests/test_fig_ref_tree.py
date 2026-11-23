"""Tests for reports/carotenoid/figures/bin/fig_ref_tree.py."""
import importlib.util
import io
import subprocess
from pathlib import Path

import pytest
from Bio import Phylo

SPEC = importlib.util.spec_from_file_location(
    "fig_ref_tree", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_ref_tree.py")
ft = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ft)

NWK = "((A:1,B:2):1,(C:1,D:3):2);"


def test_layout_x_is_cumulative_branch_length_and_y_is_leaf_order():
    tree = Phylo.read(io.StringIO(NWK), "newick")
    xy, edges, leaves = ft.layout(tree)
    assert [l.name for l in leaves] == ["A", "B", "C", "D"]
    pos = {l.name: xy[id(l)] for l in leaves}
    assert pos["A"] == (2.0, 0) and pos["B"] == (3.0, 1) and pos["C"] == (3.0, 2) and pos["D"] == (5.0, 3)
    assert len(edges) == 6


def test_layout_internal_node_y_is_mean_of_children():
    tree = Phylo.read(io.StringIO(NWK), "newick")
    xy, _, _ = ft.layout(tree)
    first = tree.root.clades[0]
    assert xy[id(first)][1] == 0.5


def test_read_cluster_sizes_counts_members_per_representative(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("r1\tr1\nr1\tm2\nr1\tm3\nr2\tr2\n\n")
    assert ft.read_cluster_sizes(p) == {"r1": 3, "r2": 1}


def test_read_roles(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("label\tour_clade\tfig1_clade\trole\nP1\ttarget\t\treference\nX\tdecoy\t\tdecoy\n")
    assert ft.read_roles(p) == {"P1": ("target", "reference"), "X": ("decoy", "decoy")}


@pytest.mark.skipif(subprocess.run(["which", "zstd"], capture_output=True).returncode != 0, reason="zstd not on PATH")
def test_main_writes_png_and_pdf(tmp_path):
    d = tmp_path / "s5"
    d.mkdir()
    (d / "ref.rooted.nwk").write_text(NWK)
    (d / "clade_members.tsv").write_text(
        "label\tour_clade\tfig1_clade\trole\nA\ttarget\t\treference\nB\ttarget\t\trepresentative\nC\tdecoy\t\tdecoy\nD\ttarget\t\trepresentative\n")
    (d / "clusters.tsv").write_text("B\tB\nB\tx\nD\tD\n")
    subprocess.run(["zstd", "-q", "-o", str(d / "clusters.tsv.zst"), str(d / "clusters.tsv")], check=True)
    out = tmp_path / "out"
    assert ft.main(["--step5-dir", str(d), "--out-dir", str(out)]) == 0
    assert (out / "fig4_reference_tree.pdf").read_bytes().startswith(b"%PDF")
    assert (out / "fig4_reference_tree.png").stat().st_size > 1000


def test_read_subclades_returns_names_and_blanks(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("label\tour_clade\tfig1_clade\trole\nA\ttarget\tEurotiomycetes\tnamed_member\nB\ttarget\t\trepresentative\n")
    assert ft.read_subclades(p) == {"A": "Eurotiomycetes", "B": ""}


@pytest.mark.skipif(subprocess.run(["which", "zstd"], capture_output=True).returncode != 0, reason="zstd not on PATH")
def test_main_draws_named_members_outgroups_and_subclades(tmp_path):
    d = tmp_path / "s5"
    d.mkdir()
    (d / "ref.rooted.nwk").write_text("((A:1,B:2):1,((C:1,D:3):2,E:1):1);")
    (d / "clade_members.tsv").write_text(
        "label\tour_clade\tfig1_clade\trole\nA\ttarget\tClassX\tnamed_member\nB\ttarget\tClassX\treference\n"
        "C\tdecoy\t\ttree_outgroup\nD\tdecoy\t\tdecoy\nE\ttarget\t\trepresentative\n")
    (d / "clusters.tsv").write_text("E\tE\n")
    subprocess.run(["zstd", "-q", "-o", str(d / "clusters.tsv.zst"), str(d / "clusters.tsv")], check=True)
    out = tmp_path / "out"
    assert ft.main(["--step5-dir", str(d), "--out-dir", str(out), "--name", "figx", "--title", "T"]) == 0
    assert (out / "figx.pdf").read_bytes().startswith(b"%PDF")
