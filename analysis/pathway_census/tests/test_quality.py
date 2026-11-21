import pandas as pd
import pytest

from pathway_census.quality import clean_mask, quality_flags, unknown_counts

REFS_LEN = {"psy": (300, 500), "ggpp": (300, 400)}
EXPECTED = {"psy": ["PF00494"], "ggpp": ["PF00348"]}
LT = "AAAAAAAA"


def _genomes():
    return pd.DataFrame({"locustag": [LT]})


def _gi(rows):
    """rows: (gene number, chrom, start, length_aa)"""
    return pd.DataFrame([dict(gene_id=f"{LT}_{n:06d}", LOCUSTAG=LT, chrom=c, start=s, end=s + 900,
                              strand=1, length=L) for n, c, s, L in rows])


def _assigned(pairs):
    return pd.DataFrame([dict(protein_id=f"{LT}_{n:06d}-T1", step_id=s, support=0.9,
                              label="assigned", fusion_candidate=False) for n, s in pairs])


def _dom(rows):
    return pd.DataFrame([dict(protein_id=f"{LT}_{n:06d}-T1", hmm_name=h, score=100.0, evalue=ev,
                              ali_from=1, ali_to=200, hmm_from=1, hmm_to=200) for n, h, ev in rows])


def _flags(assigned, domains, gi, **kw):
    out = quality_flags(assigned, domains, gi, REFS_LEN, _genomes(), expected_domains=EXPECTED, **kw)
    return out.set_index(["protein_id", "step_id"])


def test_truncated_when_below_length_floor():
    gi = _gi([(1, "s1", 1000, 209), (2, "s1", 5000, 210), (3, "s1", 9000, 400)])
    a = _assigned([(1, "psy"), (2, "psy"), (3, "ggpp")])
    d = _dom([(1, "PF00494", 1e-30), (2, "PF00494", 1e-30), (3, "PF00348", 1e-30)])
    f = _flags(a, d, gi, length_floor=0.7)  # floor for psy = 0.7 * 300 = 210
    assert bool(f.loc[(f"{LT}_000001-T1", "psy"), "truncated"]) is True   # 209 < 210
    assert bool(f.loc[(f"{LT}_000002-T1", "psy"), "truncated"]) is False  # 210 is not below 210
    assert bool(f.loc[(f"{LT}_000003-T1", "ggpp"), "truncated"]) is False


def test_missing_expected_domain_sets_domain_ok_false():
    gi = _gi([(1, "s1", 1000, 400), (2, "s1", 5000, 400), (3, "s1", 9000, 400)])
    a = _assigned([(1, "psy"), (2, "psy"), (3, "psy")])
    d = _dom([(1, "PF00494", 1e-30),      # present
              (2, "PF00348", 1e-30),      # wrong domain only
              (3, "PF00494", 1e-2)])      # present but evalue above 1e-5
    f = _flags(a, d, gi, domain_evalue=1e-5)
    assert bool(f.loc[(f"{LT}_000001-T1", "psy"), "domain_ok"]) is True
    assert bool(f.loc[(f"{LT}_000002-T1", "psy"), "domain_ok"]) is False
    assert bool(f.loc[(f"{LT}_000003-T1", "psy"), "domain_ok"]) is False


def test_fused_gene_domains_are_checked_per_step():
    gi = _gi([(1, "s1", 1000, 700)])
    a = _assigned([(1, "psy"), (1, "ggpp")])
    d = _dom([(1, "PF00494", 1e-30)])  # one domain of the fusion only
    f = _flags(a, d, gi)
    assert bool(f.loc[(f"{LT}_000001-T1", "psy"), "domain_ok"]) is True
    assert bool(f.loc[(f"{LT}_000001-T1", "ggpp"), "domain_ok"]) is False


def test_motif_ok_is_na_without_motif_data_and_set_with_it():
    gi = _gi([(1, "s1", 1000, 400), (2, "s1", 5000, 400)])
    a = _assigned([(1, "psy"), (2, "psy")])
    d = _dom([(1, "PF00494", 1e-30), (2, "PF00494", 1e-30)])
    f = _flags(a, d, gi)
    assert f.motif_ok.isna().all()
    m = pd.DataFrame([dict(protein_id=f"{LT}_000001-T1", step_id="psy", motif_present=True),
                      dict(protein_id=f"{LT}_000002-T1", step_id="psy", motif_present=False)])
    f = _flags(a, d, gi, motif_hits=m)
    assert bool(f.loc[(f"{LT}_000001-T1", "psy"), "motif_ok"]) is True
    assert bool(f.loc[(f"{LT}_000002-T1", "psy"), "motif_ok"]) is False


def _hits(classes):
    return pd.DataFrame([dict(protein_id=f"{LT}_{n:06d}-T1", taxon_class=c) for n, c in classes])


def test_contam_suspect_when_neighbours_have_non_fungal_best_hits():
    # scaffold bact: focal gene 5 with neighbours 1-4,6-9 (4 bacterial of 6 within window 3)
    rows = [(n, "bact", 1000 * n, 400) for n in range(1, 10)]
    # scaffold fung: focal gene 20, neighbours fungal
    rows += [(n, "fung", 1000 * (n - 10), 400) for n in range(11, 20)] + [(20, "fung", 10000, 400)]
    gi = _gi(rows)
    a = _assigned([(5, "psy"), (15, "psy")])
    d = _dom([(5, "PF00494", 1e-30), (15, "PF00494", 1e-30)])
    h = _hits([(2, "non_fungal"), (3, "non_fungal"), (4, "non_fungal"), (6, "non_fungal"),
               (7, "fungal"), (8, "fungal"),
               (12, "fungal"), (13, "fungal"), (14, "fungal"), (16, "fungal"), (17, "fungal"),
               (18, "non_fungal")])
    f = _flags(a, d, gi, best_hits=h, neighbour_window=3, min_neighbours_with_hits=3,
               non_fungal_fraction=0.5)
    assert bool(f.loc[(f"{LT}_000005-T1", "psy"), "contam_suspect"]) is True    # 4/6 non-fungal in window 2-8
    assert bool(f.loc[(f"{LT}_000015-T1", "psy"), "contam_suspect"]) is False   # 1/6 non-fungal


def test_contam_not_judged_with_too_few_neighbour_hits_and_without_table():
    rows = [(n, "s", 1000 * n, 400) for n in range(1, 8)]
    gi = _gi(rows)
    a = _assigned([(4, "psy")])
    d = _dom([(4, "PF00494", 1e-30)])
    h = _hits([(3, "non_fungal"), (5, "non_fungal")])  # only 2 neighbours with hits
    f = _flags(a, d, gi, best_hits=h, neighbour_window=3, min_neighbours_with_hits=3)
    assert pd.isna(f.iloc[0].contam_suspect)  # insufficient evidence: unknown, not False
    f = _flags(a, d, gi)  # no best-hit table: not evaluated
    assert pd.isna(f.iloc[0].contam_suspect)


def _pairs(rows):
    return pd.DataFrame([dict(protein_a=f"{LT}_{x:06d}-T1", protein_b=f"{LT}_{y:06d}-T1",
                              pident=p, qcov=c) for x, y, p, c in rows])


def _hap_fixture():
    # scaffold A: genes 1-5, focal 3. scaffold B: genes 11-15, focal 13. scaffold C: genes 21-25, focal 23.
    rows = [(n, "A", 1000 * n, 400) for n in range(1, 6)]
    rows += [(n, "B", 1000 * (n - 10), 400) for n in range(11, 16)]
    rows += [(n, "C", 1000 * (n - 20), 400) for n in range(21, 26)]
    gi = _gi(rows)
    a = _assigned([(3, "psy"), (13, "psy"), (23, "psy")])
    d = _dom([(3, "PF00494", 1e-30), (13, "PF00494", 1e-30), (23, "PF00494", 1e-30)])
    return gi, a, d


def test_hap_dup_suspect_when_pair_is_near_identical_and_adjacent_scaffolds_duplicate():
    gi, a, d = _hap_fixture()
    # A and B are duplicated: focal pair 3-13 plus neighbour pairs 1-11, 2-12, 4-14 (all >= 0.98).
    # C has only a near-identical focal copy (23 vs 3) with no duplicated neighbourhood.
    p = _pairs([(3, 13, 0.995, 0.98), (1, 11, 0.99, 0.95), (2, 12, 0.99, 0.95), (4, 14, 0.99, 0.95),
                (3, 23, 0.995, 0.98)])
    f = _flags(a, d, gi, pair_identity=p, neighbour_window=2, min_shared_neighbours=3,
               hap_identity=0.98, hap_coverage=0.9)
    assert bool(f.loc[(f"{LT}_000003-T1", "psy"), "hap_dup_suspect"]) is True
    assert bool(f.loc[(f"{LT}_000013-T1", "psy"), "hap_dup_suspect"]) is True
    assert bool(f.loc[(f"{LT}_000023-T1", "psy"), "hap_dup_suspect"]) is False


def test_hap_dup_not_suspect_when_pair_is_diverged_or_same_scaffold():
    gi, a, d = _hap_fixture()
    # same neighbourhood pairs, but the focal pair is 90% identical: a real paralog
    p = _pairs([(3, 13, 0.90, 0.98), (1, 11, 0.99, 0.95), (2, 12, 0.99, 0.95), (4, 14, 0.99, 0.95)])
    f = _flags(a, d, gi, pair_identity=p, neighbour_window=2, min_shared_neighbours=3)
    assert not f.hap_dup_suspect.any()
    # near-identical pair on one scaffold (tandem copy): not judged a haplotype
    rows = [(n, "A", 1000 * n, 400) for n in range(1, 8)]
    gi2 = _gi(rows)
    a2 = _assigned([(3, "psy"), (5, "psy")])
    d2 = _dom([(3, "PF00494", 1e-30), (5, "PF00494", 1e-30)])
    p2 = _pairs([(3, 5, 0.999, 1.0), (2, 4, 0.99, 1.0), (1, 6, 0.99, 1.0), (4, 7, 0.99, 1.0)])
    f2 = _flags(a2, d2, gi2, pair_identity=p2, neighbour_window=3, min_shared_neighbours=1)
    assert not f2.hap_dup_suspect.any()


def test_unknown_locustag_raises():
    gi, a, d = _hap_fixture()
    with pytest.raises(ValueError):
        quality_flags(a, d, gi, REFS_LEN, pd.DataFrame({"locustag": ["ZZZZZZZZ"]}), expected_domains=EXPECTED)


def test_clean_mask_and_unknown_counts_exact():
    fl = pd.DataFrame({
        "protein_id": list("abcde"), "step_id": ["s"] * 5,
        "domain_ok": pd.array([True, True, True, True, True], dtype="boolean"),
        "motif_ok": pd.array([pd.NA] * 5, dtype="boolean"),
        "truncated": pd.array([False, False, True, False, False], dtype="boolean"),
        "contam_suspect": pd.array([False, pd.NA, False, False, False], dtype="boolean"),
        "hap_dup_suspect": pd.array([False, False, False, True, pd.NA], dtype="boolean")})
    assert clean_mask(fl).tolist() == [True, False, False, False, False]
    assert clean_mask(fl, require=("truncated",)).tolist() == [True, True, False, True, True]
    assert clean_mask(fl, require=("contam_suspect", "hap_dup_suspect")).tolist() == \
        [True, False, True, False, False]
    u = unknown_counts(fl)
    assert u.to_dict() == {"domain_ok": 0, "motif_ok": 5, "truncated": 0, "contam_suspect": 1,
                           "hap_dup_suspect": 1}


def test_na_when_inputs_missing_for_truncated_domain_and_hap():
    gi = _gi([(1, "s1", 1000, 400), (2, "s1", 5000, pd.NA)])
    a = _assigned([(1, "nostep"), (2, "psy")])
    d = _dom([(1, "PF00494", 1e-30)])
    out = quality_flags(a, d, gi, REFS_LEN, _genomes(), expected_domains=EXPECTED).set_index("protein_id")
    r1, r2 = out.loc[f"{LT}_000001-T1"], out.loc[f"{LT}_000002-T1"]
    assert pd.isna(r1.truncated)      # step has no refs_len entry
    assert pd.isna(r1.domain_ok)      # step absent from expected_domains
    assert pd.isna(r2.truncated)      # length missing
    assert bool(r2.domain_ok) is False
    assert out.hap_dup_suspect.isna().all()  # no pair_identity table


def test_hap_na_when_gene_missing_from_gene_info():
    gi = _gi([(1, "s1", 1000, 400)])
    a = _assigned([(1, "psy"), (9, "psy")])  # gene 9 not in gene_info
    d = _dom([(1, "PF00494", 1e-30)])
    p = _pairs([(1, 9, 0.99, 1.0)])
    f = _flags(a, d, gi, pair_identity=p)
    assert bool(f.loc[(f"{LT}_000001-T1", "psy"), "hap_dup_suspect"]) is False
    assert pd.isna(f.loc[(f"{LT}_000009-T1", "psy"), "hap_dup_suspect"])
