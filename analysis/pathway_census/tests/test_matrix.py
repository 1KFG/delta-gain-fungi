import pandas as pd
import pytest

from pathway_census.defs import Pathway, Step
import subprocess
import sys
from pathlib import Path

from pathway_census.matrix import build_matrix, parse_bool_col, phylum_summary

A, B, C, D = "AAAAAAAA", "BBBBBBBB", "CCCCCCCC", "DDDDDDDD"


def _step(i):
    return Step(id=i, name=i, ec=None, refs=[], domains=[], copy_expect="single", role="catalytic")


PW = Pathway(id="t", verified=False, steps=[_step("s1"), _step("s2"), _step("s3")])


def _genomes():
    return pd.DataFrame({
        "locustag": [A, B, C, D],
        "phylum": ["Asco", "Asco", "Basidio", "Basidio"],
        "has_proteins": [True, True, True, False]})


def _asg(rows):
    return pd.DataFrame([dict(protein_id=p, step_id=s, support=0.9, label=l, fusion_candidate=False)
                         for p, s, l in rows])


def _flags(rows):
    """rows: (protein, step, truncated, contam, hap)"""
    t = lambda x: pd.array(x, dtype="boolean")
    r = list(zip(*rows))
    return pd.DataFrame({"protein_id": r[0], "step_id": r[1], "domain_ok": t([pd.NA] * len(rows)),
                         "motif_ok": t([pd.NA] * len(rows)), "truncated": t(list(r[2])),
                         "contam_suspect": t(list(r[3])), "hap_dup_suspect": t(list(r[4]))})


def _m(asg, flags=None, **kw):
    if flags is None:
        flags = _flags([(r.protein_id, r.step_id, False, False, False) for r in asg.itertuples()]
                       ) if len(asg) else _flags([("x", "s1", False, False, False)]).iloc[:0]
    return build_matrix(asg, flags, _genomes(), PW, **kw)


def _get(m, lt, step, col):
    r = m[(m.locustag == lt) & (m.step_id == step)]
    assert len(r) == 1
    return r.iloc[0][col]


def test_every_genome_step_pair_present():
    m = _m(_asg([(f"{A}_000001-T1", "s1", "assigned")]))
    assert len(m) == 4 * 3
    assert not m.duplicated(["locustag", "step_id"]).any()
    assert _get(m, A, "s1", "copy_raw") == 1
    assert _get(m, A, "s2", "copy_raw") == 0
    assert _get(m, B, "s1", "copy_raw") == 0  # has_proteins True, no assigned rows (zero, in n_genomes)
    assert not m[["copy_raw", "copy_filtered", "copy_unknown"]].isna().any().any()
    for lt in (A, B, C, D):
        assert set(m[m.locustag == lt].step_id) == {"s1", "s2", "s3"}


def test_zero_protein_genome_marked_not_zero_evidence():
    m = _m(_asg([(f"{A}_000001-T1", "s1", "assigned")]))
    d = m[m.locustag == D]
    assert (d.has_proteins == False).all() and (d.copy_raw == 0).all()  # noqa: E712
    s = phylum_summary(m, _genomes())
    basid = s[s.phylum == "Basidio"].set_index("step_id")
    assert basid.loc["s1", "n_genomes"] == 1  # C only, D excluded
    asco = s[s.phylum == "Asco"].set_index("step_id")
    assert asco.loc["s1", "n_genomes"] == 2


def test_zero_protein_genome_with_assigned_rows_is_an_error():
    with pytest.raises(ValueError, match="has_proteins"):
        _m(_asg([(f"{D}_000001-T1", "s1", "assigned")]))


def test_fused_protein_counts_once_per_step():
    p = f"{A}_000001-T1"
    m = _m(_asg([(p, "s1", "assigned"), (p, "s2", "assigned")]))
    assert _get(m, A, "s1", "copy_raw") == 1
    assert _get(m, A, "s2", "copy_raw") == 1
    assert (m[m.locustag == A].n_proteins_distinct == 1).all()


def test_filtered_excludes_truncated_and_contam_and_invariant():
    ps = [f"{A}_00000{i}-T1" for i in range(1, 6)]
    asg = _asg([(p, "s1", "assigned") for p in ps])
    fl = _flags([(ps[0], "s1", False, False, False),   # clean
                 (ps[1], "s1", True, False, False),    # truncated
                 (ps[2], "s1", False, True, False),    # contam
                 (ps[3], "s1", pd.NA, False, False),   # unknown truncated
                 (ps[4], "s1", False, False, True)])   # hap suspect only
    m = _m(asg, fl)
    assert _get(m, A, "s1", "copy_raw") == 5
    assert _get(m, A, "s1", "copy_filtered") == 2  # ps0 and ps4 (hap not required)
    assert _get(m, A, "s1", "copy_unknown") == 1
    assert _get(m, A, "s1", "n_hap_suspect") == 1
    assert _get(m, A, "s1", "n_hap_unknown") == 0
    assert _get(m, A, "s1", "copy_flagged") == 2
    assert _get(m, A, "s1", "genome_flags") == "unknown_flags"
    m2 = _m(asg, fl, require=("truncated", "contam_suspect", "hap_dup_suspect"))
    assert _get(m2, A, "s1", "copy_filtered") == 1


def test_missing_flag_row_is_unknown_not_clean():
    p = f"{A}_000001-T1"
    m = _m(_asg([(p, "s1", "assigned")]), _flags([("other", "s1", False, False, False)]).iloc[:0])
    assert _get(m, A, "s1", "copy_raw") == 1
    assert _get(m, A, "s1", "copy_filtered") == 0
    assert _get(m, A, "s1", "copy_unknown") == 1


def test_ubiquitous_family_alone_does_not_count():
    # bare HMM hit: label unassigned. Never counted toward any step.
    asg = _asg([(f"{A}_000001-T1", "s1", "unassigned"), (f"{A}_000002-T1", "s1", "assigned")])
    m = _m(asg)
    assert _get(m, A, "s1", "copy_raw") == 1
    assert (m[m.locustag == A].n_proteins_distinct == 1).all()
    m0 = _m(_asg([(f"{A}_000001-T1", "s1", "unassigned")]))
    assert (m0.copy_raw == 0).all() and (m0.copy_filtered == 0).all() and (m0.copy_unknown == 0).all()


def test_unknown_protein_locustag_raises():
    with pytest.raises(ValueError, match="LOCUSTAG"):
        _m(_asg([("ZZZZZZZZ_000001-T1", "s1", "assigned")]))


def test_phylum_summary_values():
    asg = _asg([(f"{A}_000001-T1", "s1", "assigned"), (f"{A}_000002-T1", "s1", "assigned"),
                (f"{A}_000003-T1", "s1", "assigned"), (f"{B}_000001-T1", "s1", "assigned")])
    fl = _flags([(f"{A}_000001-T1", "s1", False, False, False),
                 (f"{A}_000002-T1", "s1", False, False, False),
                 (f"{A}_000003-T1", "s1", True, False, False),
                 (f"{B}_000001-T1", "s1", pd.NA, pd.NA, False)])
    m = _m(asg, fl)
    s = phylum_summary(m, _genomes())
    assert list(s.columns) == ["phylum", "step_id", "n_genomes", "n_with_step", "median_copy",
                               "max_copy", "frac_with_step", "n_unknown", "copy_raw_total",
                               "copy_filtered_total", "copy_unknown_total", "copy_flagged_total"]
    assert len(s) == 2 * 3
    r = s[(s.phylum == "Asco") & (s.step_id == "s1")].iloc[0]
    # filtered: A=2, B=0 (unknown) -> median 1.0, max 2, with_step 1 of 2
    assert (r.n_genomes, r.n_with_step, r.median_copy, r.max_copy, r.frac_with_step, r.n_unknown) == (
        2, 1, 1.0, 2, 0.5, 1)
    r = phylum_summary(m, _genomes(), count_col="copy_raw")
    r = r[(r.phylum == "Asco") & (r.step_id == "s1")].iloc[0]
    assert (r.n_with_step, r.median_copy, r.max_copy, r.frac_with_step) == (2, 2.0, 3, 1.0)
    z = s[(s.phylum == "Basidio") & (s.step_id == "s2")].iloc[0]
    assert (z.n_genomes, z.n_with_step, z.median_copy, z.max_copy, z.frac_with_step) == (1, 0, 0.0, 0, 0.0)


def test_genome_flags_column():
    m = _m(_asg([(f"{A}_000001-T1", "s1", "assigned")]))
    assert "no_proteins" in _get(m, D, "s1", "genome_flags")
    assert _get(m, A, "s1", "genome_flags") == ""


def test_buckets_are_mutually_exclusive_for_every_row():
    ps = [f"{A}_00000{i}-T1" for i in range(1, 8)]
    asg = _asg([(p, "s1", "assigned") for p in ps])
    fl = _flags([(ps[0], "s1", False, False, False),   # all clean
                 (ps[1], "s1", True, pd.NA, pd.NA),    # mixed: True + NA -> flagged
                 (ps[2], "s1", pd.NA, False, False),   # NA + False -> unknown
                 (ps[3], "s1", True, True, False),     # both True -> flagged
                 (ps[4], "s1", pd.NA, pd.NA, pd.NA)])  # both NA -> unknown
    # ps[5] and ps[6]: missing flags rows -> unknown
    m = _m(asg, fl)
    assert (m.copy_filtered + m.copy_unknown + m.copy_flagged == m.copy_raw).all()
    r = m[(m.locustag == A) & (m.step_id == "s1")].iloc[0]
    assert (r.copy_raw, r.copy_filtered, r.copy_unknown, r.copy_flagged) == (7, 1, 4, 2)
    assert r.n_hap_unknown == 4  # ps1, ps4, ps5, ps6
    assert r.genome_flags == "unknown_flags"
    # mixed row alone
    m1 = _m(_asg([(ps[1], "s1", "assigned")]), _flags([(ps[1], "s1", True, pd.NA, pd.NA)]))
    r = m1[(m1.locustag == A) & (m1.step_id == "s1")].iloc[0]
    assert (r.copy_filtered, r.copy_unknown, r.copy_flagged, r.genome_flags) == (0, 0, 1, "")
    # other steps and genomes: all zero, flag-free
    z = m[(m.locustag == B)]
    assert (z.copy_flagged == 0).all() and (z.genome_flags == "").all()


def test_hap_unknown_and_missing_flags_row_counts():
    p1, p2 = f"{A}_000001-T1", f"{A}_000002-T1"
    m = _m(_asg([(p1, "s1", "assigned"), (p2, "s1", "assigned")]),
           _flags([(p1, "s1", False, False, pd.NA)]))
    r = m[(m.locustag == A) & (m.step_id == "s1")].iloc[0]
    assert (r.copy_filtered, r.copy_unknown, r.n_hap_unknown, r.n_hap_suspect) == (1, 1, 2, 0)
    assert r.genome_flags == "unknown_flags"


def test_parse_bool_col():
    s = parse_bool_col(pd.Series(["True", "False", "", None]), "f", allow_na=True)
    assert s.tolist()[:2] == [True, False] and s.isna().tolist() == [False, False, True, True]
    assert str(s.dtype) == "boolean"
    assert parse_bool_col(pd.Series([True, False]), "f", allow_na=False).tolist() == [True, False]
    assert parse_bool_col(pd.Series([True, None], dtype="object"), "f", allow_na=True).isna().tolist() == [False, True]
    for bad in ("1", "yes", "nan", "FALSE ", "true"):
        with pytest.raises(ValueError, match=r"flagcol.*" + repr(bad).replace(" ", r"\ ")):
            parse_bool_col(pd.Series(["True", bad]), "flagcol", allow_na=True)
    with pytest.raises(ValueError, match="has_proteins"):
        parse_bool_col(pd.Series(["True", ""]), "has_proteins", allow_na=False)
    with pytest.raises(ValueError, match="has_proteins"):
        parse_bool_col(pd.Series([True, float("nan")]), "has_proteins", allow_na=False)


CLI = Path(__file__).resolve().parents[1] / "bin" / "build_matrix.py"
PW_YAML = Path(__file__).resolve().parents[1] / "pathways" / "carotenoid.yaml"


def _cli(tmp_path, *extra, flags_text=None):
    (tmp_path / "a.tsv").write_text(
        "protein_id\tstep_id\tsupport\tlabel\tfusion_candidate\n"
        f"{A}_000001-T1\tggpp_synthase\t0.9\tassigned\tFalse\n"
        f"{A}_000002-T1\tggpp_synthase\t0.9\tassigned\tFalse\n")
    (tmp_path / "f.tsv").write_text(flags_text or (
        "protein_id\tstep_id\tdomain_ok\tmotif_ok\ttruncated\tcontam_suspect\thap_dup_suspect\n"
        f"{A}_000001-T1\tggpp_synthase\t\t\tFalse\tFalse\tFalse\n"
        f"{A}_000002-T1\tggpp_synthase\t\t\tTrue\t\t\n"))
    (tmp_path / "g.tsv").write_text(
        f"locustag\tphylum\thas_proteins\n{A}\tAsco\tTrue\n{B}\tAsco\tFalse\n")
    return subprocess.run(
        ["/usr/bin/python3.12", str(CLI), "--pathway", str(PW_YAML),
         "--assigned", str(tmp_path / "a.tsv"), "--flags", str(tmp_path / "f.tsv"),
         "--genomes", str(tmp_path / "g.tsv"), "--out-dir", str(tmp_path / "o"), *extra],
        capture_output=True, text=True)


def test_cli_refuses_unverified_without_provisional(tmp_path):
    r = _cli(tmp_path)
    assert r.returncode != 0 and "verified: false" in r.stderr
    assert not (tmp_path / "o").exists()


def test_cli_provisional_writes_outputs(tmp_path):
    r = _cli(tmp_path, "--provisional")
    assert r.returncode == 0, r.stderr
    mp = tmp_path / "o" / "carotenoid_matrix.PROVISIONAL.tsv.zst"
    sp = tmp_path / "o" / "carotenoid_phylum_summary.PROVISIONAL.tsv"
    assert mp.exists() and sp.exists()
    import io
    m = pd.read_csv(io.BytesIO(subprocess.run(["zstd", "-dc", str(mp)], capture_output=True, check=True).stdout),
                    sep="\t")
    assert len(m) == 2 * 8
    row = m[(m.locustag == A) & (m.step_id == "ggpp_synthase")].iloc[0]
    assert (row.copy_raw, row.copy_filtered, row.copy_flagged, row.copy_unknown) == (2, 1, 1, 0)
    s = pd.read_csv(sp, sep="\t", comment="#")
    g = s[(s.phylum == "Asco") & (s.step_id == "ggpp_synthase")].iloc[0]
    assert (g.n_genomes, g.n_with_step, g.max_copy) == (1, 1, 1)


def test_cli_rejects_bad_flag_value(tmp_path):
    r = _cli(tmp_path, "--provisional", flags_text=(
        "protein_id\tstep_id\tdomain_ok\tmotif_ok\ttruncated\tcontam_suspect\thap_dup_suspect\n"
        f"{A}_000001-T1\tggpp_synthase\t\t\tyes\tFalse\tFalse\n"))
    assert r.returncode != 0 and "truncated" in r.stderr and "yes" in r.stderr


# ---- I3: has_proteins given as strings must parse strictly ----
def _genomes_str():
    g = _genomes()
    g["has_proteins"] = g.has_proteins.map({True: "True", False: "False"})
    return g


def test_build_matrix_parses_string_has_proteins_strictly():
    asg = _asg([(f"{A}_000001-T1", "s1", "assigned")])
    fl = _flags([(f"{A}_000001-T1", "s1", False, False, False)])
    real = build_matrix(asg, fl, _genomes(), PW)
    strs = build_matrix(asg, fl, _genomes_str(), PW)
    pd.testing.assert_frame_equal(real, strs)
    assert (strs[strs.locustag == D].has_proteins == False).all()  # noqa: E712
    assert "no_proteins" in _get(strs, D, "s1", "genome_flags")
    bad = _genomes()
    bad["has_proteins"] = ["True", "yes", "True", "False"]
    with pytest.raises(ValueError, match="has_proteins"):
        build_matrix(asg, fl, bad, PW)


def test_phylum_summary_parses_string_has_proteins_strictly():
    m = _m(_asg([(f"{A}_000001-T1", "s1", "assigned")]))
    real = phylum_summary(m, _genomes())
    ms = m.copy()
    ms["has_proteins"] = ms.has_proteins.map({True: "True", False: "False"})
    pd.testing.assert_frame_equal(real, phylum_summary(ms, _genomes()))
    basid = phylum_summary(ms, _genomes()).set_index(["phylum", "step_id"])
    assert basid.loc[("Basidio", "s1"), "n_genomes"] == 1  # D (no proteins) excluded
    ms["has_proteins"] = "1"
    with pytest.raises(ValueError, match="has_proteins"):
        phylum_summary(ms, _genomes())


def test_phylum_summary_copy_totals():
    # final-review I9(b): bucket totals per (phylum, step) over genomes with has_proteins True.
    asg = _asg([(f"{A}_000001-T1", "s1", "assigned"), (f"{A}_000002-T1", "s1", "assigned"),
                (f"{A}_000003-T1", "s1", "assigned"), (f"{B}_000001-T1", "s1", "assigned"),
                (f"{C}_000001-T1", "s2", "assigned")])
    fl = _flags([(f"{A}_000001-T1", "s1", False, False, False),
                 (f"{A}_000002-T1", "s1", False, False, False),
                 (f"{A}_000003-T1", "s1", True, False, False),
                 (f"{B}_000001-T1", "s1", pd.NA, pd.NA, False),
                 (f"{C}_000001-T1", "s2", False, pd.NA, False)])
    m = _m(asg, fl)
    for col in ("copy_filtered", "copy_raw"):  # totals do not depend on count_col
        s = phylum_summary(m, _genomes(), count_col=col).set_index(["phylum", "step_id"])
        tot = ["copy_raw_total", "copy_filtered_total", "copy_unknown_total", "copy_flagged_total"]
        assert s.loc[("Asco", "s1"), tot].tolist() == [4, 2, 1, 1]
        assert s.loc[("Basidio", "s2"), tot].tolist() == [1, 0, 1, 0]
        assert s.loc[("Basidio", "s1"), tot].tolist() == [0, 0, 0, 0]
        # invariant carried to the totals
        assert (s.copy_filtered_total + s.copy_unknown_total + s.copy_flagged_total == s.copy_raw_total).all()


def test_no_contamination_table_gives_filtered_zero():
    # final-review I9(a): with contam_suspect NA for every copy (no contamination table), every copy
    # is copy_unknown and copy_filtered is 0. Reviewer MEASURED copy_raw 1, copy_filtered 0,
    # frac_with_step 0.0 for this case.
    asg = _asg([(f"{A}_000001-T1", "s1", "assigned")])
    fl = _flags([(f"{A}_000001-T1", "s1", False, pd.NA, pd.NA)])
    m = _m(asg, fl)
    assert (_get(m, A, "s1", "copy_raw"), _get(m, A, "s1", "copy_filtered"),
            _get(m, A, "s1", "copy_unknown")) == (1, 0, 1)
    s = phylum_summary(m, _genomes()).set_index(["phylum", "step_id"])
    r = s.loc[("Asco", "s1")]
    assert (r.copy_raw_total, r.copy_filtered_total, r.frac_with_step) == (1, 0, 0.0)
