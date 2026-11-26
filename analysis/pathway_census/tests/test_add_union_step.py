import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("add_union_step", Path(__file__).resolve().parents[1] / "bin" / "add_union_step.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_union_value_rules():
    assert mod.union_value(["0", "0", "1"]) == "1"
    assert mod.union_value(["0", "0", "0"]) == "0"
    assert mod.union_value(["NA", "NA", "NA"]) == "NA"
    assert mod.union_value(["NA", "0", "0"]) == "NA"
    assert mod.union_value(["NA", "1", "0"]) == "1"


def test_add_union_columns_and_source():
    rows = [
        {"locustag": "A", "a": "1", "b": "0", "a_source": "placed", "b_source": "placed"},
        {"locustag": "B", "a": "0", "b": "0", "a_source": "placed", "b_source": "placed"},
        {"locustag": "C", "a": "NA", "b": "NA", "a_source": "placed", "b_source": "placed"},
    ]
    out, extra = mod.add_union(rows, "any_ab", ["a", "b"])
    assert extra == ["any_ab", "any_ab_source"]
    assert [r["any_ab"] for r in out] == ["1", "0", "NA"]
    assert out[0]["any_ab_source"] == "placed"
    assert "any_ab" not in rows[0]  # input rows are not changed


def test_mixed_sources_give_hits():
    rows = [{"locustag": "A", "a": "1", "b": "0", "a_source": "placed", "b_source": "hits"}]
    out, _ = mod.add_union(rows, "u", ["a", "b"])
    assert out[0]["u_source"] == "hits"


def test_missing_member_raises():
    import pytest
    with pytest.raises(KeyError):
        mod.add_union([{"locustag": "A", "a": "1"}], "u", ["a", "zz"])


def test_drop_members_cli(tmp_path):
    src = tmp_path / "c.tsv"
    src.write_text("locustag\ta\tb\tc\ta_source\tb_source\tc_source\nX\t1\t0\t1\tplaced\tplaced\tplaced\nY\t0\t0\t0\tplaced\tplaced\tplaced\n")
    out = tmp_path / "o.tsv"
    mod.main(["--calls", str(src), "--name", "any_ab", "--members", "a", "b", "--out", str(out), "--drop-members"])
    lines = [l.split("\t") for l in out.read_text().splitlines()]
    assert lines[0] == ["locustag", "c", "c_source", "any_ab", "any_ab_source"]
    assert lines[1] == ["X", "1", "placed", "1", "placed"]
    assert lines[2] == ["Y", "0", "placed", "0", "placed"]
