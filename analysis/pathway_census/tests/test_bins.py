import subprocess
from pathlib import Path

import pandas as pd
import pytest

from pathway_census.bins import triage

B = pd.NA


def _row(pid, step, chrom="s1", start=1000, end=2000, **kw):
    r = dict(protein_id=pid, step_id=step, label="assigned", truncated=False, contam_suspect=False,
             domain_ok=True, motif_ok=True, hap_dup_suspect=False, ref_identity=0.9,
             copy_filtered=1, chrom=chrom, start=start, end=end, length=400,
             signal_peptide=False, tm_helices=0)
    r.update(kw)
    return r


def _df(rows):
    d = pd.DataFrame(rows)
    for c in ("truncated", "contam_suspect", "domain_ok", "motif_ok", "hap_dup_suspect"):
        d[c] = d[c].astype("boolean")
    return d


def _get(out, pid, step="a"):
    return out[(out.protein_id == pid) & (out.step_id == step)].iloc[0]


def test_truncated_copy_is_excluded():
    out = triage(_df([_row("AAAAAAAA_1", "a", truncated=True)]))
    r = out.iloc[0]
    assert r.bin == "excluded" and r.reasons == "truncated"


def test_neighbour_within_window_promotes():
    d = _df([_row("AAAAAAAA_1", "a", start=1000, end=2000),
             _row("AAAAAAAA_2", "b", start=12000, end=13000)])
    out = triage(d, window_bp=20000)
    assert out.in_cluster.all()


def test_neighbour_on_other_scaffold_does_not_promote():
    d = _df([_row("AAAAAAAA_1", "a"), _row("AAAAAAAA_2", "b", chrom="s2")])
    assert not triage(d, window_bp=20000).in_cluster.any()


def test_neighbour_beyond_window_does_not_promote():
    d = _df([_row("AAAAAAAA_1", "a", start=1, end=100), _row("AAAAAAAA_2", "b", start=50000, end=51000)])
    assert not triage(d, window_bp=20000).in_cluster.any()


def test_neighbour_of_same_step_does_not_count():
    d = _df([_row("AAAAAAAA_1", "a"), _row("AAAAAAAA_2", "a", start=3000, end=4000)])
    assert not triage(d, window_bp=20000).in_cluster.any()


def test_same_scaffold_name_in_other_genome_does_not_count():
    d = _df([_row("AAAAAAAA_1", "a"), _row("BBBBBBBB_1", "b")])
    assert not triage(d, window_bp=20000).in_cluster.any()


def test_reasons_lists_every_term_that_fired():
    d = _df([_row("AAAAAAAA_1", "a", ref_identity=0.3, copy_filtered=3),
             _row("AAAAAAAA_2", "b", start=3000, end=4000)])
    r = _get(triage(d, window_bp=20000), "AAAAAAAA_1")
    assert r.reasons == "divergent;lineage_expansion;in_cluster;expressible"
    assert r.bin == "tier1"


def test_same_input_same_bin():
    d = _df([_row("AAAAAAAA_1", "a", ref_identity=0.3), _row("AAAAAAAA_2", "b", start=3000, end=4000),
             _row("BBBBBBBB_1", "a", truncated=True)])
    pd.testing.assert_frame_equal(triage(d, window_bp=5000), triage(d.copy(), window_bp=5000))
    pd.testing.assert_frame_equal(triage(d), triage(d.iloc[::-1].iloc[::-1]))


def test_na_required_flag_gives_unknown():
    out = triage(_df([_row("AAAAAAAA_1", "a", motif_ok=B)]))
    r = out.iloc[0]
    assert r.bin == "unknown" and r.reasons == "unevaluated:motif_ok"


def test_known_bad_beats_na():
    out = triage(_df([_row("AAAAAAAA_1", "a", motif_ok=B, truncated=True)]))
    assert out.iloc[0].bin == "excluded"


def test_domain_ok_false_excluded():
    r = triage(_df([_row("AAAAAAAA_1", "a", domain_ok=False)])).iloc[0]
    assert r.bin == "excluded" and r.reasons == "domain_not_ok"


def test_motif_ok_false_excluded_and_require_parameter():
    d = _df([_row("AAAAAAAA_1", "a", motif_ok=False)])
    assert triage(d).iloc[0].reasons == "motif_not_ok"
    assert triage(d, require=("truncated",)).iloc[0].bin != "excluded"


def test_hap_dup_true_demotes_tier1():
    d = _df([_row("AAAAAAAA_1", "a", ref_identity=0.3, copy_filtered=3, hap_dup_suspect=True)])
    r = triage(d).iloc[0]
    assert r.bin == "tier2"
    assert r.reasons == "divergent;lineage_expansion;expressible;hap_dup_suspect_demoted"


def test_hap_dup_na_does_nothing_but_is_shown():
    d = _df([_row("AAAAAAAA_1", "a", ref_identity=0.3, copy_filtered=3, hap_dup_suspect=B)])
    r = triage(d).iloc[0]
    assert r.bin == "tier1" and r.reasons.endswith("hap_dup_unevaluated")


def test_unassigned_rows_never_appear():
    d = _df([_row("AAAAAAAA_1", "a"), _row("AAAAAAAA_2", "b", label="unassigned")])
    out = triage(d)
    assert list(out.protein_id) == ["AAAAAAAA_1"]


def test_zero_terms_is_tier3():
    d = _df([_row("AAAAAAAA_1", "a", length=5)])  # not expressible, identity high, one copy
    r = triage(d).iloc[0]
    assert r.bin == "tier3" and r.reasons == ""


def test_na_identity_does_not_promote():
    r = triage(_df([_row("AAAAAAAA_1", "a", ref_identity=B)])).iloc[0]
    assert not r.divergent


def test_soluble_step_needs_no_signal_or_tm():
    d = _df([_row("AAAAAAAA_1", "a", signal_peptide=True), _row("AAAAAAAA_2", "b", tm_helices=5)])
    out = triage(d, soluble_steps=("a", "b"))
    assert not out.expressible.any()
    assert triage(d).expressible.all()


CLI = Path(__file__).resolve().parents[1] / "bin" / "triage.py"
PW = Path(__file__).parent / "data" / "pathway_ok.yaml"


def _run(tmp_path, *extra, pathway=None):
    d = _df([_row("AAAAAAAA_1", "a", ref_identity=0.3, copy_filtered=3), _row("AAAAAAAA_2", "b", start=3000, end=4000)])
    inp = tmp_path / "copies.tsv"
    d.astype(object).where(d.notna(), "").to_csv(inp, sep="\t", index=False)
    return subprocess.run(["/usr/bin/python3.12", str(CLI), "--pathway", str(pathway or PW),
                           "--copies", str(inp), "--out", str(tmp_path / "o.tsv"), *extra],
                          capture_output=True, text=True)


def _verified(tmp_path):
    p = tmp_path / "pw.yaml"
    p.write_text(PW.read_text().replace("verified: false", "verified: true", 1))
    return p


def test_cli_refuses_unverified_pathway(tmp_path):
    txt = PW.read_text()
    if "verified: true" in txt.split("steps:")[0]:
        pytest.skip("fixture pathway is verified")
    r = _run(tmp_path)
    assert r.returncode != 0 and "verified: false" in r.stderr
    assert not (tmp_path / "o.tsv").exists()


def test_cli_runs_with_provisional_and_writes_params_line(tmp_path):
    r = _run(tmp_path, "--provisional")
    assert r.returncode == 0, r.stderr
    out = (tmp_path / "o.PROVISIONAL.tsv")
    first = out.read_text().splitlines()[0]
    assert first.startswith("# params:") and "window_bp=20000" in first
    t = pd.read_csv(out, sep="\t", comment="#")
    assert set(t.bin) <= {"tier1", "tier2", "tier3", "excluded", "unknown"} and len(t) == 2


def test_cli_rejects_bad_boolean(tmp_path):
    d = _df([_row("AAAAAAAA_1", "a")]).astype(object)
    d["truncated"] = "yes"
    inp = tmp_path / "c.tsv"
    d.to_csv(inp, sep="\t", index=False)
    r = subprocess.run(["/usr/bin/python3.12", str(CLI), "--pathway", str(PW), "--copies", str(inp),
                        "--out", str(tmp_path / "o.tsv"), "--provisional"], capture_output=True, text=True)
    assert r.returncode != 0 and "truncated" in r.stderr
