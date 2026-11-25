"""Tests for bin/apply_prefix_aliases.py."""
import importlib.util
import subprocess
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("apply_prefix_aliases", Path(__file__).parents[1] / "bin" / "apply_prefix_aliases.py")
ap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ap)

G = [dict(locustag="A", species="x"), dict(locustag="B", species="y")]


def test_adds_column_and_keeps_rest():
    rows, amap = ap.apply_aliases(G, [dict(locustag="A", alias_prefix="Z")])
    assert amap == {"A": "Z"}
    assert [r["alias_prefix"] for r in rows] == ["Z", ""]
    assert [r["species"] for r in rows] == ["x", "y"]
    assert "alias_prefix" not in G[0]  # input rows are not changed


@pytest.mark.parametrize("al", [
    [dict(locustag="Q", alias_prefix="Z")],                                   # unknown locustag
    [dict(locustag="A", alias_prefix="B")],                                   # alias is another table locustag
    [dict(locustag="A", alias_prefix="Z"), dict(locustag="A", alias_prefix="Y")],
    [dict(locustag="A", alias_prefix="Z"), dict(locustag="B", alias_prefix="Z")],
])
def test_bad_aliases_raise(al):
    with pytest.raises(ValueError):
        ap.apply_aliases(G, al)


def test_cli_zst_roundtrip(tmp_path):
    g = tmp_path / "g.tsv"
    g.write_text("locustag\tspecies\nA\tx\nB\ty\n")
    subprocess.run(["zstd", "-q", "-f", str(g), "-o", str(g) + ".zst"], check=True)
    al = tmp_path / "al.tsv"
    al.write_text("locustag\talias_prefix\tnote\nA\tZ\tn\n")
    out = tmp_path / "o.tsv.zst"
    assert ap.main(["--genomes", str(g) + ".zst", "--aliases", str(al), "--out", str(out)]) == 0
    txt = subprocess.run(["zstd", "-dc", str(out)], capture_output=True, check=True).stdout.decode()
    assert txt == "locustag\tspecies\talias_prefix\nA\tx\tZ\nB\ty\t\n"
