"""Tests for bin/drop_phylum_calls.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("drop_phylum_calls", Path(__file__).parents[1] / "bin" / "drop_phylum_calls.py")
dp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dp)


def rows():
    mk = lambda g, ph, a, b: dict(locustag=g, PHYLUM=ph, searched="1" if a != "NA" else "0", s1=a, s2=b, s1_source="placed", s2_source="hits")
    return [mk("A", "Microsporidia", "0", "0"), mk("B", "Ascomycota", "1", "0"), mk("C", "Microsporidia", "NA", "NA"), mk("D", "Rozellomycota", "0", "0"),
            mk("E", "", "1", "1")]


def test_only_named_phylum_is_dropped_and_other_tips_untouched():
    new, dropped = dp.drop_phyla(rows(), ["Microsporidia"])
    by = {r["locustag"]: r for r in new}
    assert dropped == ["A", "C"]
    assert (by["A"]["s1"], by["A"]["s2"], by["A"]["dropped_phylum"], by["A"]["searched"]) == ("NA", "NA", "1", "0")
    assert by["B"]["s1"] == "1" and by["B"]["dropped_phylum"] == "0"
    assert by["D"]["s1"] == "0"  # Rozellomycota kept
    assert by["E"]["s2"] == "1"  # empty phylum kept


def test_two_phyla_and_no_match():
    _, d2 = dp.drop_phyla(rows(), ["Microsporidia", "Rozellomycota"])
    assert d2 == ["A", "C", "D"]
    _, d0 = dp.drop_phyla(rows(), ["Nonexistent"])
    assert d0 == []
