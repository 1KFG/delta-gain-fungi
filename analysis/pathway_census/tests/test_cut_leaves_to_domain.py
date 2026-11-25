import importlib.util
from pathlib import Path

import pandas as pd
import pytest

_SPEC = importlib.util.spec_from_file_location("cut_leaves_to_domain", Path(__file__).resolve().parents[1] / "bin" / "cut_leaves_to_domain.py")
cl = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cl)


def test_cut_bounds_pad_and_clamp():
    assert cl.cut_bounds(300, 100, 200) == (80, 220)
    assert cl.cut_bounds(300, 5, 290) == (1, 300)
    assert cl.cut_bounds(100, 10, 90, pad=0) == (10, 90)


def test_cut_bounds_rejects_bad_coordinates():
    with pytest.raises(ValueError):
        cl.cut_bounds(100, 0, 10)
    with pytest.raises(ValueError):
        cl.cut_bounds(100, 50, 40)
    with pytest.raises(ValueError):
        cl.cut_bounds(100, 150, 160, pad=0)


def test_best_domain_coords_takes_highest_score():
    df = pd.DataFrame({"protein_id": ["p", "p", "q"], "score": [30.0, 90.0, 5.0],
                       "ali_from": [1, 40, 7], "ali_to": [20, 120, 30]})
    assert cl.best_domain_coords(df) == {"p": (40, 120), "q": (7, 30)}
    assert cl.best_domain_coords(pd.DataFrame(columns=["protein_id", "score", "ali_from", "ali_to"])) == {}


def test_cut_records_cuts_and_keeps_unhit():
    seq = "".join(chr(65 + i % 20) for i in range(100))
    recs = [("p", "p desc", seq), ("n", "n desc", "MKV")]
    out, table = cl.cut_records(recs, {"p": (40, 60)}, pad=5)
    assert out[0] == ("p", "p desc", seq[34:65])
    assert out[1] == ("n", "n desc", "MKV")
    assert table[0] == ("p", 100, 40, 60, 35, 65, 31)
    assert table[1][2] == "" and table[1][6] == 3
