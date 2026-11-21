import pandas as pd
import pytest

from pathway_census.defs import Pathway, Step
from pathway_census.matrix import build_matrix, phylum_summary

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
    flagged_true = 2  # truncated, contam
    assert 2 + 1 + flagged_true == 5
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
                               "max_copy", "frac_with_step", "n_unknown"]
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
