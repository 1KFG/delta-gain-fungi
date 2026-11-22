"""Tests for bin/step5_place.py (census step 5, step-generic). Synthetic inputs only."""
import csv
import gzip
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

from pathway_census import hits

ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location("step5_place", ROOT / "bin" / "step5_place.py")
s5 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(s5)

needs_samtools = pytest.mark.skipif(not hits.tool_available("samtools"),
                                    reason="samtools not on PATH (module load samtools/1.19.2)")
needs_bgzip = pytest.mark.skipif(not hits.tool_available("bgzip"),
                                 reason="bgzip not on PATH (module load htslib or samtools)")
needs_zstd = pytest.mark.skipif(not hits.tool_available("zstd"), reason="zstd not on PATH")
needs_upo = pytest.mark.skipif(not s5.UPO_BIN.exists(), reason="enzymes_UPO upo_place not found")


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


# ---- text input ----

def test_open_text_reads_plain_and_gz(tmp_path):
    (tmp_path / "a.txt").write_text("x\n")
    (tmp_path / "b.txt.gz").write_bytes(gzip.compress(b"y\n"))
    assert s5.open_text(tmp_path / "a.txt").read() == "x\n"
    assert s5.open_text(tmp_path / "b.txt.gz").read() == "y\n"


@needs_zstd
def test_open_text_reads_zst(tmp_path):
    (tmp_path / "c.txt").write_text("z\n")
    subprocess.run(["zstd", "-q", str(tmp_path / "c.txt")], check=True)
    assert s5.open_text(tmp_path / "c.txt.zst").read() == "z\n"


def test_read_tsv_skips_comment_lines(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("# PROVISIONAL\na\tb\n1\t2\n")
    assert s5.read_tsv(p) == [{"a": "1", "b": "2"}]


def test_read_fasta_text_uses_first_word():
    assert s5.read_fasta_text(">x desc\nAB\nCD\n>y\nEF\n") == [("x", "ABCD"), ("y", "EF")]


# ---- extract ----

def test_locustag_and_path():
    assert s5.locustag_of("F8027CAE_004645-T1") == "F8027CAE"
    p = s5.protein_path(Path("/db"), "F01BDC31")
    assert p == Path("/db/proteins/00/F01BDC31.proteins.fa.gz")  # sha1("F01BDC31") starts 00


def test_cut_domain_pads_and_clamps():
    seq = "".join(chr(65 + i % 26) for i in range(100))
    assert s5.cut_domain(seq, 30, 60, pad=20)[:2] == (10, 80)
    assert s5.cut_domain(seq, 5, 95, pad=20)[:2] == (1, 100)
    s, e, sub = s5.cut_domain(seq, 30, 60, pad=0)
    assert sub == seq[29:60] and len(sub) == 31


def test_cut_domain_rejects_bad_coordinates():
    with pytest.raises(ValueError):
        s5.cut_domain("ABCDE", 4, 2)
    with pytest.raises(ValueError):
        s5.cut_domain("ABCDE", 0, 2)
    with pytest.raises(ValueError):
        s5.cut_domain("ABCDE", 50, 60, pad=5)  # window lies past the protein end


def test_best_rows_filters_step_and_keeps_best_score():
    rows = [{"protein_id": "p", "step_id": "A", "score": "10"},
            {"protein_id": "p", "step_id": "A", "score": "30"},
            {"protein_id": "p", "step_id": "B", "score": "99"},
            {"protein_id": "q", "step_id": "A", "score": "5"}]
    got, n_dup = s5.best_rows(rows, "A")
    assert {r["protein_id"]: r["score"] for r in got} == {"p": "30", "q": "5"}
    assert n_dup == 1


@needs_samtools
@needs_bgzip
def test_extract_end_to_end_with_failures(tmp_path):
    seqdb = tmp_path / "seqdb"
    tag = "F01BDC31"
    f = s5.protein_path(seqdb, tag)
    f.parent.mkdir(parents=True)
    plain = tmp_path / "p.fa"
    plain.write_text(f">{tag}_1-T1\n" + "A" * 50 + "\n" + f">{tag}_2-T1\n" + "C" * 30 + "\n")
    subprocess.run(["bgzip", "-c", str(plain)], stdout=open(f, "wb"), check=True)
    subprocess.run(["samtools", "faidx", str(f)], check=True)
    cand = tmp_path / "cand.tsv"
    cols = ["protein_id", "hmm_name", "score", "evalue", "ali_from", "ali_to", "hmm_from", "hmm_to", "step_id"]
    rows = [[f"{tag}_1-T1", "h", "50", "1e-9", "10", "40", "1", 30, "S"],
            [f"{tag}_2-T1", "h", "50", "1e-9", "5", "20", "1", 15, "S"],
            [f"{tag}_3-T1", "h", "50", "1e-9", "5", "20", "1", 15, "S"],      # id not in file
            ["FFFFFFFF_1-T1", "h", "50", "1e-9", "5", "20", "1", 15, "S"],    # genome file missing
            [f"{tag}_1-T1", "h", "50", "1e-9", "5", "20", "1", 15, "OTHER"]]  # other step, ignored
    write_tsv(cand, cols, rows)
    r = s5.extract([cand], "S", seqdb, tmp_path / "o.faa", tmp_path / "c.tsv", tmp_path / "f.tsv",
                   pad=5, workers=2)
    assert r["domains_extracted"] == 2 and r["failures"] == 2 and r["rows_in_step"] == 4
    seqs = dict(s5.read_fasta_text((tmp_path / "o.faa").read_text()))
    assert len(seqs[f"{tag}_1-T1"]) == 41  # 5..45 = 41 residues
    assert len(seqs[f"{tag}_2-T1"]) == 25  # 1..25 (clamped at the start)
    fails = {r["protein_id"]: r["reason"] for r in s5.read_tsv(tmp_path / "f.tsv")}
    assert fails[f"{tag}_3-T1"] == "id_not_in_protein_file"
    assert fails["FFFFFFFF_1-T1"] == "protein_file_missing"


# ---- representatives ----

def test_parse_clusters(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("a\ta\na\tb\nc\tc\n")
    assert s5.parse_clusters(p) == {"a": ["a", "b"], "c": ["c"]}


def test_stratified_quota_proportional_with_minimum_one():
    q = s5.stratified_quota({"A": 900, "B": 90, "C": 10}, 100)
    assert sum(q.values()) == 100
    assert q["C"] >= 1 and q["A"] > q["B"] > q["C"]


def test_stratified_quota_small_set_is_kept_whole():
    assert s5.stratified_quota({"A": 3, "B": 2}, 10) == {"A": 3, "B": 2}


def test_stratified_quota_never_exceeds_stratum_size():
    q = s5.stratified_quota({"A": 1000, "B": 2}, 10)
    assert q["B"] <= 2 and sum(q.values()) == 10


def test_stratified_quota_more_strata_than_n_gives_one_each():
    q = s5.stratified_quota({"A": 5, "B": 5, "C": 5}, 2)
    assert all(v == 1 for v in q.values())


def test_stratified_draw_is_seeded_and_exact():
    by = {"A": [f"a{i}" for i in range(50)], "B": [f"b{i}" for i in range(10)]}
    one, q = s5.stratified_draw(by, 12, seed=1)
    two, _ = s5.stratified_draw(by, 12, seed=1)
    other, _ = s5.stratified_draw(by, 12, seed=2)
    assert one == two and one != other and len(one) == 12
    assert sum(1 for x in one if x.startswith("b")) == q["B"]


@needs_zstd
def test_select_reps_reports_counts(tmp_path):
    clu = tmp_path / "clu.tsv"
    clu.write_text("AAAAAAAA_1\tAAAAAAAA_1\nAAAAAAAA_1\tBBBBBBBB_1\nBBBBBBBB_2\tBBBBBBBB_2\n"
                   "CCCCCCCC_1\tCCCCCCCC_1\n")
    gen = tmp_path / "g.tsv"
    write_tsv(gen, ["locustag", "phylum"], [["AAAAAAAA", "Asco"], ["BBBBBBBB", "Asco"], ["CCCCCCCC", "Basidio"]])
    r = s5.select_reps(clu, gen, 2, 1, tmp_path / "ids", tmp_path / "rep.tsv")
    assert r == {"clusters": 3, "representatives_total": 3, "representatives_selected": 2,
                 "phyla": 2, "capped": True}
    ids = (tmp_path / "ids").read_text().split()
    assert "CCCCCCCC_1" in ids and len(ids) == 2  # the single Basidio cluster keeps its place


# ---- clades ----

def bip_of(newick):
    upo = s5.load_upo()
    top = upo.parse_newick(newick)
    for i, n in enumerate(s5._nodes(top)[1:]):
        n.edge = i
    return upo.edge_bipartitions(top)


NW = "(((r1:1,r2:1):1,(rep1:1,r3:1):1):1,((d1:1,d2:1):1,drep:1):1,og:1);"


@needs_upo
def test_define_clades_target_and_decoy():
    c = s5.define_clades(bip_of(NW), {"r1", "r2", "r3"}, {"d1", "d2"}, {"og"})
    assert c["target"] == frozenset({"r1", "r2", "r3", "rep1"})
    assert c["decoy"] == frozenset({"d1", "d2", "drep"})
    assert c["decoys_share_side_without_refs"] is True


@needs_upo
def test_define_clades_refs_mixed_with_decoy_gives_none():
    nw = "((r1:1,d1:1):1,(r2:1,r3:1):1,(d2:1,og:1):1);"
    c = s5.define_clades(bip_of(nw), {"r1", "r2", "r3"}, {"d1", "d2"}, {"og"})
    assert c["target"] is None


@needs_upo
def test_define_clades_decoy_is_everything_outside_target_and_outgroup():
    nw = "((r1:1,r2:1):1,(d1:1,x1:1):1,(d2:1,og:1):1);"
    c = s5.define_clades(bip_of(nw), {"r1", "r2"}, {"d1", "d2"}, {"og"})
    assert c["target"] == frozenset({"r1", "r2"})
    assert c["decoy"] == frozenset({"d1", "x1", "d2"})


@needs_upo
def test_define_clades_refs_split_falls_back_to_decoy_leaves():
    nw = "((r1:1,d1:1):1,(r2:1,r3:1):1,(d2:1,og:1):1);"
    c = s5.define_clades(bip_of(nw), {"r1", "r2", "r3"}, {"d1", "d2"}, {"og"})
    assert c["target"] is None and c["decoy"] == frozenset({"d1", "d2"})


@needs_upo
def test_define_clades_without_outgroup():
    nw = "((r1:1,r2:1):1,(d1:1,d2:1):1,x:1);"
    c = s5.define_clades(bip_of(nw), {"r1", "r2"}, {"d1", "d2"}, set())
    assert c["target"] == frozenset({"r1", "r2"}) and c["outgroup"] == frozenset()


def test_smallest_side_picks_the_smallest():
    sides = [frozenset("abc"), frozenset("ab"), frozenset("abd")]
    assert s5.smallest_side(sides, "a", "d") == frozenset("ab")
    assert s5.smallest_side(sides, "a", "ab") is None


def test_leaf_roles_priority():
    r = s5.leaf_roles({"a", "b", "c", "d", "e"}, {"a"}, {"b"}, {"c"}, {"d"})
    assert r == {"a": "reference", "b": "decoy", "c": "outgroup", "d": "representative", "e": "other"}


def test_write_clade_members_one_label_per_leaf(tmp_path):
    clades = {"target": frozenset({"a", "d"}), "decoy": frozenset({"b"}), "outgroup": frozenset({"c"})}
    roles = {"a": "reference", "b": "decoy", "c": "outgroup", "d": "representative"}
    s5.write_clade_members(tmp_path / "m.tsv", clades, roles)
    rows = {r["label"]: r for r in s5.read_tsv(tmp_path / "m.tsv")}
    assert rows["d"]["our_clade"] == "target" and rows["d"]["role"] == "representative"
    assert rows["c"]["our_clade"] == "outgroup" and rows["a"]["fig1_clade"] == ""


@needs_upo
def test_run_clades_reports_and_roots(tmp_path):
    (tmp_path / "t.nwk").write_text(NW)
    (tmp_path / "refs.faa").write_text(">r1\nA\n>r2\nA\n>r3\nA\n")
    (tmp_path / "dec.faa").write_text(">d1\nA\n>d2\nA\n")
    (tmp_path / "og.faa").write_text(">og\nA\n")
    (tmp_path / "reps").write_text("rep1\ndrep\n")
    rep = s5.run_clades(tmp_path / "t.nwk", tmp_path / "refs.faa", tmp_path / "dec.faa", tmp_path / "og.faa",
                        tmp_path / "reps", tmp_path, "S")
    assert rep["refs_form_clade_free_of_decoys"] is True
    assert rep["target_refs"] == 3 and rep["target_representatives"] == 1
    assert rep["decoy_representatives"] == 1
    assert (tmp_path / "ref.rooted.nwk").read_text().count("og") == 1
    assert (tmp_path / "clade_members.tsv").exists()


@needs_upo
def test_run_clades_without_outgroup_roots_on_decoys(tmp_path):
    (tmp_path / "t.nwk").write_text("((r1:1,r2:1):1,(d1:1,d2:1):1,x:1);")
    (tmp_path / "refs.faa").write_text(">r1\nA\n>r2\nA\n")
    (tmp_path / "dec.faa").write_text(">d1\nA\n>d2\nA\n")
    (tmp_path / "reps").write_text("")
    rep = s5.run_clades(tmp_path / "t.nwk", tmp_path / "refs.faa", tmp_path / "dec.faa", None,
                        tmp_path / "reps", tmp_path, "S")
    assert rep["rooting_leaves"] == ["d1", "d2"] and "decoys" in rep["rooting_source"]


@needs_upo
def test_run_clades_missing_leaf_stops(tmp_path):
    (tmp_path / "t.nwk").write_text("((r1:1,r2:1):1,(d1:1,d2:1):1,x:1);")
    (tmp_path / "refs.faa").write_text(">r1\nA\n>zzz\nA\n")
    (tmp_path / "dec.faa").write_text(">d1\nA\n")
    (tmp_path / "reps").write_text("")
    with pytest.raises(SystemExit):
        s5.run_clades(tmp_path / "t.nwk", tmp_path / "refs.faa", tmp_path / "dec.faa", None,
                      tmp_path / "reps", tmp_path, "S")


# ---- placements ----

def arow(q, call, best, lwr, flags=""):
    return {"query": q, "our_clade": call, "our_best": best, "our_lwr": lwr, "our_stem_lwr": "0.000",
            "pendant_length": "0.2", "best_edge_lwr": "0.9", "flags": flags}


def test_to_placements_only_target_gets_a_clade():
    rows = [arow("p1", "target", "target", "0.900"), arow("p2", "decoy", "decoy", "0.800"),
            arow("p3", "unassigned", "target", "0.400", "long_pendant")]
    plc, det = s5.to_placements(rows, "S")
    assert [(r["protein_id"], r["clade"], r["support"]) for r in plc] == [
        ("p1", "target", "0.900"), ("p2", "", "0.800"), ("p3", "", "0.400")]
    assert [r["clade_call"] for r in det] == ["target", "decoy", "unassigned"]
    assert [r["long_pendant"] for r in det] == [False, False, True]
    assert all(r["step_id"] == "S" for r in plc)


def test_run_placements_writes_both_files(tmp_path):
    write_tsv(tmp_path / "a.tsv", ["query", "our_clade", "our_best", "our_lwr", "our_stem_lwr",
                                   "pendant_length", "best_edge_lwr", "flags"],
              [["p1", "target", "target", "0.9", "0", "0.1", "0.9", ""]])
    r = s5.run_placements(tmp_path / "a.tsv", "S", tmp_path / "p.tsv", tmp_path / "d.tsv")
    assert r == {"queries": 1}
    assert s5.read_tsv(tmp_path / "p.tsv")[0]["clade"] == "target"


def test_placements_feed_assign_steps_decoy_is_unassigned():
    import pandas as pd
    from pathway_census.assign import assign_steps
    from pathway_census.defs import load_pathway
    pw = load_pathway(ROOT / "pathways" / "carotenoid.yaml")
    plc, _ = s5.to_placements([arow("p1", "target", "target", "0.9"), arow("p2", "decoy", "decoy", "0.9")],
                              "phytoene_synthase")
    cand = pd.DataFrame({"protein_id": ["p1", "p2"], "hmm_name": "h", "score": 100.0, "evalue": 1e-9,
                         "ali_from": 1, "ali_to": 50, "hmm_from": 1, "hmm_to": 50,
                         "step_id": "phytoene_synthase"})
    res = assign_steps(cand, pd.DataFrame(plc).astype({"support": float}), pw, 0.5)
    assert dict(zip(res.protein_id, res.label)) == {"p1": "assigned", "p2": "unassigned"}


# ---- summary ----

def test_support_distribution():
    d = s5.support_distribution(["0.5", "0.8", "0.95", "1.0"])
    assert d["n"] == 4 and d["n_ge_0.5"] == 4 and d["n_ge_0.9"] == 2 and d["n_ge_0.95"] == 2
    assert d["min"] == 0.5 and d["max"] == 1.0


def test_support_distribution_empty():
    d = s5.support_distribution([])
    assert d["n"] == 0 and d["median"] == 0.0 and d["n_ge_0.5"] == 0


def test_genomes_per_assigned():
    rows = [{"protein_id": "AAAAAAAA_1-T1", "label": "assigned"},
            {"protein_id": "AAAAAAAA_2-T1", "label": "assigned"},
            {"protein_id": "BBBBBBBB_1-T1", "label": "assigned"},
            {"protein_id": "CCCCCCCC_1-T1", "label": "unassigned"}]
    g = s5.genomes_per_assigned(rows)
    assert g["assigned_proteins"] == 3 and g["genomes_with_assigned"] == 2
    assert g["copy_distribution"] == {1: 1, 2: 1}


def test_run_summarize_collects_files(tmp_path):
    (tmp_path / "extract_summary.json").write_text(json.dumps({"domains_extracted": 5}))
    write_tsv(tmp_path / "placements_detail.tsv",
              ["protein_id", "clade_call", "support", "long_pendant"],
              [["p1", "target", "0.9", "False"], ["p2", "decoy", "0.7", "True"]])
    write_tsv(tmp_path / "assigned.tsv", ["protein_id", "step_id", "support", "label", "fusion_candidate"],
              [["AAAAAAAA_1-T1", "S", "0.9", "assigned", "False"]])
    n = s5.run_summarize(tmp_path, tmp_path / "assigned.tsv", tmp_path / "SUMMARY.tsv")
    rows = {r["key"]: r for r in s5.read_tsv(tmp_path / "SUMMARY.tsv")}
    assert n == len(rows)
    assert rows["extract.domains_extracted"]["value"] == "5"
    assert rows["placement.call.decoy"]["value"] == "1"
    assert rows["placement.long_pendant"]["value"] == "1"
    assert rows["assign_steps.label.assigned"]["value"] == "1"
    assert all(r["status"] == "MEASURED" for r in rows.values())
