"""Tests for reports/melanin/bin/build_genome_calls.py (rules) and summarise_calls.py (four classes)."""
import importlib.util
from pathlib import Path

import pandas as pd

BIN = Path(__file__).parents[1] / "reports" / "melanin" / "bin"


def load(name):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


bgc = load("build_genome_calls")
sc = load("summarise_calls")


def frame(rows):
    steps = bgc.DHN_STEPS + bgc.LDOPA_STEPS
    cols = steps + [f"hit_{s}" for s in steps]
    df = pd.DataFrame(0, index=range(len(rows)), columns=cols)
    for i, on in enumerate(rows):
        for s in on:
            df.loc[i, s] = 1
            df.loc[i, f"hit_{s}"] = 1
    return df


def test_route_needs_every_step_of_one_route():
    full_a = ["thn_synthase", "t4hn_reductase", "scytalone_dehydratase", "t3hn_reductase"]
    hepta_no_t3hn = ["heptaketide_pks", "ayg1", "t4hn_reductase", "scytalone_dehydratase"]
    df = bgc.apply_rules(frame([full_a, hepta_no_t3hn, ["thn_synthase"], []]))
    assert list(df.D_route) == [1, 0, 0, 0]
    assert list(df.D_pks) == [1, 1, 1, 0]
    assert list(df.D_any) == [1, 1, 1, 0]


def test_ldopa_rules_do_not_use_dhn_laccase_and_need_both_asp_steps():
    df = bgc.apply_rules(frame([["laccase_abr2"], ["mco_abr1"], ["aspulvinone_e_synthase"],
                                ["aspulvinone_e_synthase", "asp_melanin_tyrosinase"], ["tyrosinase"], ["laccase_diphenol"]]))
    assert list(df.L_route) == [0, 0, 0, 1, 1, 1]
    assert list(df.L_nontyr) == [0, 0, 0, 1, 0, 1]
    assert list(df.L_any) == [0, 0, 1, 1, 1, 1]
    assert list(df.L_tyr) == [0, 0, 0, 0, 1, 0]


def test_family_present_steps_are_not_counted():
    df = frame([["hppd"], ["hmga"]])
    df["hppd"] = [1, 0]
    df["hmga"] = [0, 1]
    out = bgc.apply_rules(df)
    assert out.L_any.sum() == 0 and out.D_any.sum() == 0


def test_four_class_counts_sum_to_n():
    d = pd.DataFrame({"D": [1, 1, 0, 0, 0], "L": [1, 0, 1, 0, 0]})
    r = sc.four_class(d, "D", "L")
    assert (r.DHN_only, r.LDOPA_only, r.both, r.neither, r.n) == (1, 1, 1, 2, 5)


def test_alias_maps_new_prefix_to_table_locustag(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text("locustag\talias_prefix\nOLD\tNEW\n")
    al = bgc.load_aliases(f)
    assert bgc.genome_of("NEW_0001-T1", al) == "OLD" and bgc.genome_of("X_1", al) == "X"
