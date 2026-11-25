import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

BIN = Path(__file__).parents[1] / "bin"
_SPEC = importlib.util.spec_from_file_location("step5_run_report", BIN / "step5_run_report.py")
rr = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rr)
needs_upo = pytest.mark.skipif(not rr.sp.UPO_BIN.exists(), reason="enzymes_UPO tools not found")


def test_genome_of():
    assert rr.genome_of("F85C5601_003317-T1") == "F85C5601"


def test_class_table_merges_small_and_counts():
    searched = ["a1", "a2", "a3", "b1", "c1"]
    label = {"a1": "A", "a2": "A", "a3": "A", "b1": "B"}
    rows = rr.class_table(searched, label, {"a1", "a2", "b1"}, min_n=2)
    d = {r[0]: r for r in rows}
    assert d["A"] == ("A", 3, 2, 2 / 3)
    assert d["(other, small groups)"][1:3] == (2, 1)      # b1 (assigned) and c1 (no label)


def test_nearest_summary_median_and_candidate_sum():
    near = [("r1", "x", 0.2, "reference"), ("r2", "x", 0.4, "reference"), ("r3", "y", 0.1, "tree_outgroup")]
    rows = rr.nearest_summary(near, {"r1": 10, "r2": 5})
    assert rows == [("reference", 2, 15, 0.3), ("tree_outgroup", 1, 1, 0.1)]


@needs_upo
def test_run_end_to_end(tmp_path):
    (tmp_path / "ref.raxml.bestTree").write_text("(((r1:1,rep1:1):1,rep2:1):1,((m1:1,rep3:1):1,og:1):1);")
    rows = ["label\tour_clade\tfig1_clade\trole", "r1\ttarget\t\treference", "rep1\ttarget\t\trepresentative",
            "rep2\ttarget\t\trepresentative", "m1\t\t\tnamed_member", "rep3\tdecoy\t\trepresentative", "og\toutgroup\t\ttree_outgroup"]
    (tmp_path / "clade_members.tsv").write_text("\n".join(rows) + "\n")
    (tmp_path / "clade_report.json").write_text(json.dumps({"step_id": "s", "refs_form_clade_free_of_decoys": True, "target_size": 3,
                                                            "smallest_target_size": 1, "target_refs": 1, "target_representatives": 2}))
    subprocess.run(["zstd", "-q", "-f", "-o", str(tmp_path / "clusters.tsv.zst")],
                   input=b"rep1\trep1\nrep1\tG1_1\nrep2\trep2\nrep3\trep3\nrep3\tG2_1\nrep3\tG2_2\n", check=True)
    (tmp_path / "placements_detail.tsv").write_text(
        "protein_id\tstep_id\tclade_call\tbest_clade\tsupport\tstem_lwr\tpendant_length\tbest_edge_lwr\tlong_pendant\tsubclade_call\tsubclade_support\n"
        "G1_1\ts\ttarget\ttarget\t0.9\t0\t0.1\t0.5\tFalse\t\t\n" "G2_1\ts\tdecoy\tdecoy\t0.8\t0\t0.1\t0.5\tFalse\t\t\n"
        "G2_2\ts\tunassigned\tdecoy\t0.4\t0\t0.1\t0.5\tFalse\t\t\n")
    (tmp_path / "candidates_step.tsv").write_text("protein_id\tstep_id\nG1_1\ts\nG2_1\ts\nG2_2\ts\n")
    (tmp_path / "assigned.PROVISIONAL.tsv").write_text("protein_id\tlabel\nG1_1\tassigned\nG2_1\tunassigned\n")
    (tmp_path / "genomes.tsv").write_text("locustag\tphylum\tclass_\tgenus\thas_proteins\nG1\tAsco\tC1\tg\tTrue\nG2\tAsco\tC1\tg\tTrue\nG3\tBasi\tC2\tg\tFalse\n")
    s = rr.run(tmp_path, tmp_path / "genomes.tsv", min_n=1)
    assert s["candidates_equal_placements"] is True and s["n_placements"] == 3
    assert s["calls"] == {"target": 1, "decoy": 1, "unassigned": 1}
    assert s["assigned_proteins"] == 1 and s["genomes_with_assigned"] == 1 and s["genomes_searched"] == 2
    assert s["leaves_by_role"]["representative"] == 3
    near = (tmp_path / "report" / "nearest_labelled.tsv").read_text().splitlines()
    assert near[0].startswith("representative") and len(near) == 4
    cls = (tmp_path / "report" / "assigned_by_class.tsv").read_text()
    assert "class\tC1\t2\t1\t0.5000" in cls
