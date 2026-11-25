import importlib.util
import json
from pathlib import Path

import pytest

BIN = Path(__file__).parents[1] / "bin"
_SPEC = importlib.util.spec_from_file_location("step5_multi_ref_clades", BIN / "step5_multi_ref_clades.py")
mr = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(mr)
sp = mr.sp
needs_upo = pytest.mark.skipif(not sp.UPO_BIN.exists(), reason="enzymes_UPO tools not found")

# refs r1 r2 form a clade with rep1; named member m1 and outgroup og are outside
NW = "((((r1:1,r2:1):1,rep1:1):1,rep2:1):1,((m1:1,rep3:1):1,(og:1,rep4:1):1):1);"


def files(tmp_path, nw=NW, refs=">r1\nMK\n>r2\nMK\n"):
    (tmp_path / "refs.faa").write_text(refs)
    (tmp_path / "og.faa").write_text(">og\nMK\n")
    (tmp_path / "mem.faa").write_text(">m1\nMK\n")
    (tmp_path / "reps").write_text("rep1\nrep2\nrep3\nrep4\n")
    (tmp_path / "t.nwk").write_text(nw)


def call(tmp_path):
    return mr.run(tmp_path / "t.nwk", tmp_path / "refs.faa", None, tmp_path / "og.faa", tmp_path / "mem.faa",
                  tmp_path / "reps", tmp_path, "s")


@needs_upo
def test_two_references_largest_side(tmp_path):
    files(tmp_path)
    assert call(tmp_path) == 0
    rep = json.loads((tmp_path / "clade_report.json").read_text())
    assert rep["refs_form_clade_free_of_decoys"] is True
    assert rep["smallest_target_size"] == 2            # default rule: the two references only
    assert rep["target_size"] == 4                      # r1 r2 rep1 rep2: stops at the labelled leaves
    assert rep["target_refs"] == 2 and rep["target_representatives"] == 2
    rows = (tmp_path / "clade_members.tsv").read_text().splitlines()
    assert "rep2\ttarget\t\trepresentative" in rows and "rep3\tdecoy\t\trepresentative" in rows


@needs_upo
def test_references_split_by_a_labelled_leaf_gives_exit_3(tmp_path):
    files(tmp_path, nw="((r1:1,m1:1):1,(r2:1,rep1:1):1,(og:1,rep2:1):1);")
    (tmp_path / "reps").write_text("rep1\nrep2\n")
    rc = mr.run(tmp_path / "t.nwk", tmp_path / "refs.faa", None, tmp_path / "og.faa", tmp_path / "mem.faa",
                tmp_path / "reps", tmp_path, "s")
    assert rc == 3
    rep = json.loads((tmp_path / "clade_report.json").read_text())
    assert rep["refs_form_clade_free_of_decoys"] is False and "target_size" not in rep
    assert not (tmp_path / "clade_members.tsv").exists()


@needs_upo
def test_no_labelled_leaf_exits_4(tmp_path):
    files(tmp_path)
    assert mr.run(tmp_path / "t.nwk", tmp_path / "refs.faa", None, None, None, tmp_path / "reps", tmp_path, "s") == 4


@needs_upo
def test_keeps_default_files_once(tmp_path):
    files(tmp_path)
    (tmp_path / "clade_members.tsv").write_text("old\n")
    (tmp_path / "clade_report.json").write_text("{}")
    call(tmp_path)
    assert (tmp_path / "clade_members.smallest.tsv").read_text() == "old\n"
    call(tmp_path)                                       # second run must not overwrite the kept copy
    assert (tmp_path / "clade_members.smallest.tsv").read_text() == "old\n"


def test_member_that_is_reference_stops(tmp_path):
    files(tmp_path)
    (tmp_path / "mem.faa").write_text(">r1\nMK\n")
    with pytest.raises(SystemExit):
        call(tmp_path)
