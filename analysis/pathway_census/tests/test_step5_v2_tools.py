import importlib.util
from pathlib import Path

BIN = Path(__file__).parents[1] / "bin"


def load(name):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cc = load("step5_compare_calls")
idd = load("step5_identity_distribution")


def test_compare_counts_transitions_and_direction():
    old = {"a": "target", "b": "decoy", "c": "decoy", "d": "target", "x": "target"}
    new = {"a": "target", "b": "target", "c": "unassigned", "d": "target", "y": "decoy"}
    r = cc.compare(old, new)
    assert r["n_both"] == 4 and r["only_old"] == 1 and r["only_new"] == 1
    assert r["transitions"][("target", "target")] == 2
    assert r["changed"] == {("decoy", "target"): 1, ("decoy", "unassigned"): 1}


def test_read_calls_plain_and_long_pendant(tmp_path):
    p = tmp_path / "d.tsv"
    p.write_text("protein_id\tclade_call\tlong_pendant\nA\ttarget\tTrue\nB\tdecoy\tFalse\nC\ttarget\tTrue\n")
    calls, lp = cc.read_calls(p)
    assert calls == {"A": "target", "B": "decoy", "C": "target"} and lp == 2


def test_cli_writes_transition_table(tmp_path):
    o, n = tmp_path / "o.tsv", tmp_path / "n.tsv"
    o.write_text("protein_id\tclade_call\tlong_pendant\nA\tdecoy\tFalse\n")
    n.write_text("protein_id\tclade_call\tlong_pendant\nA\ttarget\tFalse\n")
    import subprocess, sys
    out = tmp_path / "t.tsv"
    r = subprocess.run([sys.executable, str(BIN / "step5_compare_calls.py"), "--old", str(o), "--new", str(n), "--out", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "decoy>target: 1" in r.stdout
    assert "decoy\ttarget\t1" in out.read_text()


def test_bin_of_edges():
    assert idd.bin_of(49.9) == "<50" and idd.bin_of(50) == "50-70" and idd.bin_of(94.99) == "90-95"
    assert idd.bin_of(99.5) == "99-<100" and idd.bin_of(100.0) == "100"


def test_duplicate_summary():
    d = idd.duplicate_summary({"a": "MK", "b": "MK", "c": "MR"})
    assert d == {"sequences": 3, "distinct": 2, "largest_group": 2}


def test_nearest_and_pairs_ignore_self_and_unknown():
    rows = [("a", "a", 100.0), ("a", "b", 90.0), ("b", "a", 91.0), ("a", "c", 70.0), ("a", "zz", 99.0)]
    nearest, pairs = idd.nearest_and_pairs(rows, ["a", "b", "c"])
    assert nearest == {"a": 90.0, "b": 91.0}
    assert pairs == {("a", "b"): 91.0, ("a", "c"): 70.0}


def test_read_fasta_joins_wrapped_lines(tmp_path):
    p = tmp_path / "s.faa"
    p.write_text(">a x\nMK\nLV\n>b\nMR\n")
    assert idd.read_fasta(p) == {"a": "MKLV", "b": "MR"}
