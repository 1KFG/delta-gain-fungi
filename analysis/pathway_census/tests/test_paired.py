import io
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pathway_census.paired import (accession_from_asmid, count_pairs, deposited_protein_path, paired_table,
                                   summarize_paired, summarize_paired_by_busco)

STEPS = ["s1", "s2"]


def _genomes(rows):
    """rows: (locustag, provider, busco, prior_status, has_bfd, has_dep)"""
    return pd.DataFrame([dict(
        locustag=l, asmid=f"GCA_{i:09d}.1_x", folder=f"GCA_{i:09d}.1_x", annotation_provider=p, busco_complete=b, prior_status=ps,
        has_proteins=hb, has_deposited_proteins=hd, species="sp", phylum="Asco", class_="c",
        genus="g", genus_source="genus_field", prior_folder_state="present", n50=1, downloaded=True)
        for i, (l, p, b, ps, hb, hd) in enumerate(rows)])


def _matrix(counts, has=None, unknown=None):
    """counts: {(locustag, step): copy_raw}; one row per (locustag, step)."""
    tags = sorted({k[0] for k in counts})
    rows = []
    for t in tags:
        for s in STEPS:
            c = counts.get((t, s), 0)
            u = (unknown or {}).get((t, s), 0)
            rows.append(dict(locustag=t, step_id=s, copy_raw=c, copy_filtered=c - u, copy_unknown=u,
                             copy_flagged=0, n_hap_suspect=0, n_hap_unknown=0, n_proteins_distinct=c,
                             has_proteins=(has or {}).get(t, True), genome_flags=""))
    return pd.DataFrame(rows)


def _pairs(deltas, provider="P1", busco=99.0):
    """Build (dep, bfd, genomes) with one genome per entry in deltas for step s1 (BFD - deposited)."""
    g, dep, bfd = [], {}, {}
    for i, d in enumerate(deltas):
        t = f"G{i:07d}"
        g.append((t, provider, busco, "baseline_ok", True, True))
        dep[(t, "s1")] = 0
        bfd[(t, "s1")] = d
        dep[(t, "s2")] = 0
        bfd[(t, "s2")] = 0
    return _matrix(dep), _matrix(bfd), _genomes(g)


def test_delta_sign():
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True)])
    t = paired_table(_matrix({("A0000001", "s1"): 2}), _matrix({("A0000001", "s1"): 3}), g)
    r = t[t.step_id == "s1"].iloc[0]
    assert (r.copy_deposited, r.copy_bfd, r.delta) == (2, 3, 1)
    assert list(t.columns[:7]) == ["locustag", "step_id", "copy_deposited", "copy_bfd", "delta",
                                   "annotation_provider", "busco_complete"]


def test_same_protein_id_in_both_sources_not_merged():
    # The same protein id exists in both sources; counts are per source and never merged by id.
    # Matrices carry counts only, so a shared id shows up as 1 copy in each source, not 1 total.
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True)])
    t = paired_table(_matrix({("A0000001", "s1"): 1}), _matrix({("A0000001", "s1"): 1}), g)
    r = t[t.step_id == "s1"].iloc[0]
    assert (r.copy_deposited, r.copy_bfd, r.delta) == (1, 1, 0)
    assert len(t) == 2


def test_pair_requires_both_sources():
    g = _genomes([
        ("A0000001", "P", 99.0, "baseline_ok", True, True),      # pair
        ("B0000002", "P", 99.0, "baseline_ok", True, True),      # missing from BFD matrix
        ("C0000003", "P", 99.0, "no_prior_annotation_present", True, False),
        ("D0000004", "P", 99.0, "baseline_ok", False, True),     # no BFD proteins
        ("E0000005", "P", 99.0, "baseline_ok", True, False),     # no deposited proteins
        ("F0000006", "P", 99.0, "baseline_ok", True, True)])     # deposited matrix has_proteins False
    dep = _matrix({("A0000001", "s1"): 1, ("B0000002", "s1"): 1, ("F0000006", "s1"): 0},
                  has={"F0000006": False})
    bfd = _matrix({("A0000001", "s1"): 1, ("C0000003", "s1"): 1, ("D0000004", "s1"): 0, ("E0000005", "s1"): 1,
                   ("F0000006", "s1"): 1}, has={"D0000004": False})
    t, n = paired_table(dep, bfd, g, return_unpaired=True)
    assert set(t.locustag) == {"A0000001"}
    assert n == {"not_baseline": 1, "no_bfd_proteins": 2, "no_deposited_proteins": 1,
                 "deposited_file_missing": 1, "unknown_flag": 0}
    assert t.attrs["n_unpaired"] == n


def test_matrix_locustag_not_in_genomes_raises():
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True)])
    with pytest.raises(ValueError, match="not in genome"):
        paired_table(_matrix({("Z0000009", "s1"): 1}), _matrix({("A0000001", "s1"): 1}), g)


def test_delta_filtered_na_when_unknown():
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True),
                  ("B0000002", "P", 99.0, "baseline_ok", True, True)])
    dep = _matrix({("A0000001", "s1"): 2, ("B0000002", "s1"): 2}, unknown={("B0000002", "s1"): 1})
    bfd = _matrix({("A0000001", "s1"): 3, ("B0000002", "s1"): 3})
    t = paired_table(dep, bfd, g)
    a = t[(t.locustag == "A0000001") & (t.step_id == "s1")].iloc[0]
    b = t[(t.locustag == "B0000002") & (t.step_id == "s1")].iloc[0]
    assert a.delta_filtered == 1
    assert pd.isna(b.delta_filtered)
    assert (b.unknown_deposited, b.unknown_bfd, b.delta) == (1, 0, 1)
    assert b.copy_filtered_deposited == 1 and b.copy_filtered_bfd == 3


def test_summary_by_provider_drops_providers_below_min_pairs():
    d1, b1, g1 = _pairs([1, 0, 2, 1], provider="BIG")
    d2, b2, g2 = _pairs([1], provider="SMALL")
    g2["locustag"] = ["Z0000000"]
    d2["locustag"] = "Z0000000"
    b2["locustag"] = "Z0000000"
    t = paired_table(pd.concat([d1, d2]), pd.concat([b1, b2]), pd.concat([g1, g2]))
    s, dropped = summarize_paired(t, min_pairs=3, n_boot=50, seed=1)
    assert set(s.annotation_provider) == {"BIG", "ALL"}
    assert dropped.to_dict("records") == [{"provider": "SMALL", "n_pairs": 1}]
    allrow = s[(s.annotation_provider == "ALL") & (s.step_id == "s1")].iloc[0]
    assert allrow.n_pairs == 5 and allrow.mean_delta == pytest.approx(1.0)
    assert list(s.columns) == ["step_id", "annotation_provider", "n_pairs", "mean_delta", "ci_low",
                               "ci_high", "frac_bfd_more", "frac_bfd_less"]
    big = s[(s.annotation_provider == "BIG") & (s.step_id == "s1")].iloc[0]
    assert (big.frac_bfd_more, big.frac_bfd_less) == (0.75, 0.0)


def test_ci_is_built_from_genome_resampling_and_reflects_outlier():
    # One genome of 20 has delta +50. This test shows the interval is built by resampling genomes
    # with a documented scheme (matches independent numpy code) and that the upper bound reflects the
    # outlier. It does not prove that a protein-level resample would behave differently.
    deltas = [0] * 19 + [50]
    d, b, g = _pairs(deltas)
    t = paired_table(d, b, g)
    s, _ = summarize_paired(t, min_pairs=1, n_boot=2000, seed=7)
    r = s[(s.annotation_provider == "ALL") & (s.step_id == "s1")].iloc[0]
    assert r.mean_delta == pytest.approx(2.5)
    assert r.ci_low == 0.0
    assert r.ci_high >= 7.5  # k copies of the outlier in 20 draws gives mean 2.5k; P(k<=2) < 0.975
    # independent numpy reference with the documented scheme
    rng = np.random.default_rng(7)
    x = np.array(deltas, dtype=float)
    means = x[rng.integers(0, len(x), size=(2000, len(x)))].mean(axis=1)
    assert r.ci_low == pytest.approx(np.percentile(means, 2.5))
    assert r.ci_high == pytest.approx(np.percentile(means, 97.5))
    # same seed, same answer
    s2, _ = summarize_paired(t, min_pairs=1, n_boot=2000, seed=7)
    pd.testing.assert_frame_equal(s, s2)


def test_seed_and_n_boot_are_honoured():
    d, b, g = _pairs([0, 0, 1, 5, 0, 2, 0, 9, 1, 0])
    t = paired_table(d, b, g)
    s1, _ = summarize_paired(t, min_pairs=1, n_boot=300, seed=1)
    s2, _ = summarize_paired(t, min_pairs=1, n_boot=300, seed=2)
    s3, _ = summarize_paired(t, min_pairs=1, n_boot=301, seed=1)
    row = lambda x: x[(x.annotation_provider == "ALL") & (x.step_id == "s1")].iloc[0]
    assert (row(s1).ci_low, row(s1).ci_high) != (row(s2).ci_low, row(s2).ci_high)
    assert (row(s1).ci_low, row(s1).ci_high) != (row(s3).ci_low, row(s3).ci_high)
    x = np.array([0, 0, 1, 5, 0, 2, 0, 9, 1, 0], dtype=float)
    for ns, sd, r in ((300, 1, row(s1)), (301, 1, row(s3))):
        means = x[np.random.default_rng(sd).integers(0, 10, size=(ns, 10))].mean(axis=1)
        assert r.ci_high == pytest.approx(np.percentile(means, 97.5))


@pytest.mark.parametrize("name", ["ALL", "(missing)"])
def test_reserved_provider_names_raise(name):
    d, b, g = _pairs([1, 1], provider=name)
    with pytest.raises(ValueError, match="reserved"):
        summarize_paired(paired_table(d, b, g), min_pairs=1, n_boot=10, seed=0)


def test_unknown_flag_counted_separately():
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True),
                  ("B0000002", "P", 99.0, "baseline_ok", True, pd.NA),
                  ("C0000003", "P", 99.0, "baseline_ok", True, False)])
    g["has_deposited_proteins"] = g.has_deposited_proteins.astype("boolean")
    m = _matrix({("A0000001", "s1"): 1, ("B0000002", "s1"): 1, ("C0000003", "s1"): 1})
    t, n = paired_table(m, m, g, return_unpaired=True)
    assert n["no_deposited_proteins"] == 2 and n["unknown_flag"] == 1


def test_all_zero_step_still_reported():
    d, b, g = _pairs([0, 0, 0])
    s, _ = summarize_paired(paired_table(d, b, g), min_pairs=1, n_boot=20, seed=0)
    r = s[(s.annotation_provider == "ALL") & (s.step_id == "s2")].iloc[0]
    assert (r.n_pairs, r.mean_delta, r.ci_low, r.ci_high) == (3, 0, 0, 0)


def test_null_provider_is_labelled_not_dropped():
    d, b, g = _pairs([1, 1])
    g["annotation_provider"] = None
    s, dropped = summarize_paired(paired_table(d, b, g), min_pairs=1, n_boot=10, seed=0)
    assert "(missing)" in set(s.annotation_provider) and dropped.empty


def test_busco_bins_with_unknown_bin():
    g, dep, bfd = [], {}, {}
    for i, (bu, d) in enumerate([(50.0, 3), (79.9, 1), (80.0, 0), (96.0, 2), (100.0, 2), (None, 5)]):
        t = f"G{i:07d}"
        g.append((t, "P", bu, "baseline_ok", True, True))
        dep[(t, "s1")] = 0
        bfd[(t, "s1")] = d
    t = paired_table(_matrix(dep), _matrix(bfd), _genomes(g))
    r = summarize_paired_by_busco(t, bins=(0, 80, 95, 100))
    s1 = r[r.step_id == "s1"].set_index("busco_bin")
    assert s1.loc["[0,80)", "n_pairs"] == 2 and s1.loc["[0,80)", "mean_delta"] == 2.0
    assert s1.loc["[80,95)", "n_pairs"] == 1 and s1.loc["[80,95)", "mean_delta"] == 0.0
    assert s1.loc["[95,100]", "n_pairs"] == 2 and s1.loc["[95,100]", "mean_delta"] == 2.0
    assert s1.loc["unknown", "n_pairs"] == 1 and s1.loc["unknown", "mean_delta"] == 5.0


def test_busco_out_of_range_raises():
    d, b, g = _pairs([1], busco=101.0)
    with pytest.raises(ValueError, match="busco"):
        summarize_paired_by_busco(paired_table(d, b, g))


def test_accession_from_asmid():
    assert accession_from_asmid("GCF_010015735.1_Aaoar1") == "GCF_010015735.1"
    with pytest.raises(ValueError):
        accession_from_asmid("weird")


def test_count_pairs_on_fixtures(tmp_path):
    root = tmp_path / "ASM"
    g = _genomes([
        ("A0000001", "P", 99.0, "baseline_ok", True, True),    # pair, file present
        ("B0000002", "P", 99.0, "baseline_ok", True, True),    # file missing
        ("C0000003", "P", 99.0, "annotated_folder_missing_protein", True, False),
        ("D0000004", "P", 99.0, "baseline_ok", False, True),
        ("E0000005", "P", 99.0, "baseline_ok", True, False),
        ("F0000006", "P", 99.0, "no_prior_annotation_present", True, True)])
    f0 = g.folder[0]
    (root / f0).mkdir(parents=True)
    (root / f0 / f"{f0}_protein.faa.gz").write_bytes(b"x")
    # B: the old spec layout <accession>/protein.faa is NOT the deposited layout and is not found
    accb = accession_from_asmid(g.asmid[1])
    (root / accb).mkdir(parents=True)
    (root / accb / "protein.faa").write_text(">x\nM\n")
    r = count_pairs(g, root)
    assert r["n_pairs"] == 1 and r["n_genomes"] == 6
    assert r["n_unpaired"] == {"not_baseline": 2, "no_bfd_proteins": 1, "no_deposited_proteins": 1,
                               "deposited_file_missing": 1, "unknown_flag": 0}
    assert r["pair_locustags"] == ["A0000001"]


def test_count_pairs_needs_folder(tmp_path):
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True)]).drop(columns="folder")
    with pytest.raises(ValueError, match="folder"):
        count_pairs(g, tmp_path)


def test_count_pairs_null_folder_is_file_missing(tmp_path):
    g = _genomes([("A0000001", "P", 99.0, "baseline_ok", True, True)])
    g["folder"] = None
    r = count_pairs(g, tmp_path)
    assert r["n_pairs"] == 0 and r["n_unpaired"]["deposited_file_missing"] == 1


def test_deposited_protein_path_order(tmp_path):
    fo = "GCF_000001985.1_JCVI-PMFA1-2.0"
    d = tmp_path / fo
    gz, plain = d / f"{fo}_protein.faa.gz", d / f"{fo}_protein.faa"
    # neither exists: the primary (.faa.gz) path is returned
    assert deposited_protein_path(tmp_path, fo) == gz
    d.mkdir()
    plain.write_text(">x\nM\n")
    assert deposited_protein_path(tmp_path, fo) == plain  # fallback only because it exists
    gz.write_bytes(b"x")
    assert deposited_protein_path(tmp_path, fo) == gz     # .faa.gz wins when both exist
    with pytest.raises(ValueError, match="folder"):
        deposited_protein_path(tmp_path, None)


# ---- CLI ----
CLI = Path(__file__).resolve().parents[1] / "bin" / "paired_compare.py"
PW_YAML = Path(__file__).resolve().parents[1] / "pathways" / "carotenoid.yaml"
CSTEPS = ["ggpp_synthase"]


def _cli(tmp_path, *extra, gtext=None):
    cols = "locustag\tstep_id\tcopy_raw\tcopy_filtered\tcopy_unknown\tcopy_flagged\tn_hap_suspect\t" \
           "n_hap_unknown\tn_proteins_distinct\thas_proteins\tgenome_flags\n"
    (tmp_path / "dep.tsv").write_text(cols + "AAAAAAAA\tggpp_synthase\t1\t1\t0\t0\t0\t0\t1\tTrue\t\n"
                                      "BBBBBBBB\tggpp_synthase\t0\t0\t0\t0\t0\t0\t0\tTrue\t\n")
    (tmp_path / "bfd.tsv").write_text(cols + "AAAAAAAA\tggpp_synthase\t3\t3\t0\t0\t0\t0\t3\tTrue\t\n"
                                      "BBBBBBBB\tggpp_synthase\t0\t0\t0\t0\t0\t0\t0\tTrue\t\n")
    (tmp_path / "g.tsv").write_text(gtext or (
        "locustag\tasmid\tannotation_provider\tbusco_complete\tprior_status\thas_proteins\thas_deposited_proteins\n"
        "AAAAAAAA\tGCA_1.1_a\tNCBI\t99\tbaseline_ok\tTrue\tTrue\n"
        "BBBBBBBB\tGCA_2.1_b\tNCBI\t\tbaseline_ok\tTrue\tTrue\n"))
    return subprocess.run(
        ["/usr/bin/python3.12", str(CLI), "--pathway", str(PW_YAML),
         "--deposited-matrix", str(tmp_path / "dep.tsv"), "--bfd-matrix", str(tmp_path / "bfd.tsv"),
         "--genomes", str(tmp_path / "g.tsv"), "--out-dir", str(tmp_path / "o"), *extra],
        capture_output=True, text=True)


def test_cli_refuses_unverified_without_provisional(tmp_path):
    r = _cli(tmp_path)
    assert r.returncode != 0 and "verified: false" in r.stderr
    assert not (tmp_path / "o").exists()


def test_cli_provisional_writes_outputs(tmp_path):
    r = _cli(tmp_path, "--provisional", "--min-pairs", "1", "--n-boot", "50", "--seed", "3")
    assert r.returncode == 0, r.stderr
    o = tmp_path / "o"
    for n in ("carotenoid_paired_table.PROVISIONAL.tsv.zst", "carotenoid_paired_summary.PROVISIONAL.tsv",
              "carotenoid_paired_dropped_providers.PROVISIONAL.tsv",
              "carotenoid_paired_by_busco.PROVISIONAL.tsv", "carotenoid_paired_unpaired.PROVISIONAL.json"):
        assert (o / n).exists(), n
    t = pd.read_csv(io.BytesIO(subprocess.run(
        ["zstd", "-dc", str(o / "carotenoid_paired_table.PROVISIONAL.tsv.zst")], capture_output=True,
        check=True).stdout), sep="\t", comment="#")
    assert t.loc[t.locustag == "AAAAAAAA", "delta"].iloc[0] == 2
    s = pd.read_csv(o / "carotenoid_paired_summary.PROVISIONAL.tsv", sep="\t", comment="#")
    a = s[s.annotation_provider == "ALL"].iloc[0]
    assert a.n_pairs == 2 and a.mean_delta == 1.0
    js = json.loads((o / "carotenoid_paired_unpaired.PROVISIONAL.json").read_text())
    assert js["n_pairs"] == 2
    assert js["params"] == {"min_pairs": 1, "n_boot": 50, "seed": 3, "busco_bins": [0, 80, 95, 100]}
    for n in ("carotenoid_paired_summary.PROVISIONAL.tsv", "carotenoid_paired_dropped_providers.PROVISIONAL.tsv",
              "carotenoid_paired_by_busco.PROVISIONAL.tsv"):
        assert "# params: min_pairs=1 n_boot=50 seed=3 busco_bins=0,80,95,100" in (o / n).read_text(), n
    zt = subprocess.run(["zstd", "-dc", str(o / "carotenoid_paired_table.PROVISIONAL.tsv.zst")],
                        capture_output=True, text=True, check=True).stdout
    assert "# params: min_pairs=1 n_boot=50 seed=3 busco_bins=0,80,95,100" in zt


def test_cli_rejects_bad_boolean(tmp_path):
    r = _cli(tmp_path, "--provisional", gtext=(
        "locustag\tasmid\tannotation_provider\tbusco_complete\tprior_status\thas_proteins\thas_deposited_proteins\n"
        "AAAAAAAA\tGCA_1.1_a\tNCBI\t99\tbaseline_ok\tyes\tTrue\n"))
    assert r.returncode != 0 and "has_proteins" in r.stderr and "yes" in r.stderr
