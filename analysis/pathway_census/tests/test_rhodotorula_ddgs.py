import importlib.util
from pathlib import Path

import pandas as pd

p = Path(__file__).resolve().parents[1] / "curation_scripts/mycosporine/rhodotorula_ddgs_check.py"
spec = importlib.util.spec_from_file_location("rdc", p)
rdc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rdc)


def test_bins_and_sensitivity():
    s = pd.Series([None, 10, 50, 149.9, 150, 299, 300, 463.4, 463.5, 900])
    b = rdc.bin_scores(s, 463.5)
    assert b == {"no_hit_E10": 1, "<50": 1, "50-150": 2, "150-300": 2, "300-463.5": 2, ">=463.5": 2}
    sens = rdc.sensitivity(s, 463.5)
    assert [r["n_genomes_at_or_above"] for r in sens] == [2, 3, 3, 5]
    assert sens[1]["cutoff"] == 417.15


def test_best_hit_coverage():
    d = pd.DataFrame([("a", "h", 10.0, 1e-3, 1, 5, 1, 5), ("b", "h", 100.0, 1e-30, 1, 200, 11, 210)],
                     columns=["protein_id", "hmm_name", "score", "evalue", "ali_from", "ali_to", "hmm_from", "hmm_to"])
    r = rdc.best_hit(d, 400)
    assert r["best_protein"] == "b" and r["hmm_coverage"] == 0.5 and r["n_proteins_hit"] == 2
    assert rdc.best_hit(d.iloc[0:0], 400)["best_score"] == ""


def test_select_genomes_cap_and_proteins():
    g = pd.DataFrame(dict(locustag=["C", "A", "B", "D"], genus=["X", "X", "X", "Y"], species=["s"] * 4,
                          has_proteins=["True", "True", "False", "True"]))
    m = pd.DataFrame(dict(locustag=["A", "B", "C", "D"], proteins_path=["a", "b", "c", "d"]))
    r = rdc.select_genomes(g, m, {"X": 1, "Y": None})
    assert list(r.locustag) == ["A", "D"]
