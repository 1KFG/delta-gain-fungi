"""Offline tests for shortlist.py. Run: python3 -m unittest discover -s candidates/test"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import shortlist as sl  # noqa: E402


def row(p, clade="Clade IV", lwr="0.99", e="E", flags="", length="260", lineage="Sordariomycetes",
        species="Fusarium x"):
    return {"protein_id": p, "place_our_clade": clade, "place_our_lwr": lwr,
            "place_fig1_clade": "", "place_fig1_lwr": "0", "place_flags": flags,
            "p39_Cys_proximal": "C", "p109_E_of_EHD": "E", "p111_D_of_EHD": "D", "p180_Acid": e,
            "length": length, "LINEAGE_GROUP": lineage, "SPECIES": species}


def qc(partial="false", start="TRUE", stop="TRUE"):
    return {"is_partial": partial, "has_start_codon": start, "has_stop_codon": stop}


PRED = lambda r: r["place_our_clade"] == "Clade IV"  # noqa: E731


class T(unittest.TestCase):
    def test_eligibility_reasons(self):
        rows = [row("ok1"), row("pend", flags="long_pendant"), row("part"), row("noglu", e="A"),
                row("low", lwr="0.6"), row("other", clade="Clade V")]
        q = {p: qc() for p in ("ok1", "pend", "noglu", "low", "other")}
        q["part"] = qc(partial="true")
        kept, why, _ = sl.eligible(rows, PRED, "place_our_lwr", True, q)
        self.assertEqual([r["protein_id"] for r in kept], ["ok1"])
        self.assertEqual((why["in pool"], why["long_pendant"], why["incomplete gene model"],
                          why["catalytic residues"], why["LWR < 0.95"]), (5, 1, 1, 1, 1))

    def test_pool_c_ignores_glu180(self):
        kept, _, _ = sl.eligible([row("noglu", e="A")], PRED, "place_our_lwr", False,
                                 {"noglu": qc()})
        self.assertEqual(len(kept), 1)

    def test_select_cluster_genus_lineage_and_refseq(self):
        pool = [row("a1", species="Fusarium a"), row("a2", species="Botrytis b"),   # cluster A
                row("b1", species="Fusarium c"),                                    # cluster B, same genus as a1
                row("c1", lineage="Eurotiomycetes", species="Aspergillus d"),       # cluster C
                row("d1", lineage="Dothideomycetes", species="Alternaria e")]       # cluster D
        cl = {"a1": "A", "a2": "A", "b1": "B", "c1": "C", "d1": "D"}
        size = {"A": 50, "B": 40, "C": 10, "D": 5}
        gcf = lambda r: r["protein_id"] == "a2"  # noqa: E731
        got = [r["protein_id"] for _, _, r in sl.select(pool, 3, "place_our_lwr", cl, size, gcf)]
        # cluster A represented by its RefSeq member a2; lineages round-robin; b1 is not needed
        self.assertEqual(got[0], "a2")
        self.assertEqual(set(got), {"a2", "c1", "d1"})

    def test_select_one_per_genus(self):
        pool = [row("x1", species="Fusarium a"), row("x2", species="Fusarium b")]
        got = sl.select(pool, 5, "place_our_lwr", {}, {}, lambda r: False)
        self.assertEqual(len(got), 1)

    def test_percentile(self):
        self.assertEqual(sl.percentile([1, 2, 3, 4, 5], 0.5), 3)


if __name__ == "__main__":
    unittest.main()
