"""Offline tests for upo_place.py. Run: python3 -m unittest discover -s upo_place/test"""
import csv
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import upo_place as up  # noqa: E402
from tree_outline import parse_newick  # noqa: E402

TREE = "((A:1{0},B:1{1}):1{2},(C:1{3},D:1{4}):1{5},(E:1{6},F:1{7}):1{8}){9};"
CLADES = ("label\tfig1_clade\tour_clade\n"
          "A\tBig\tClade X\nB\tBig\tClade X\nC\tBig\t\nD\tBig\t\nE\t\t\nF\t\t\n")
FIELDS = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]


def jplace(placements, tree=TREE):
    return {"tree": tree, "fields": FIELDS, "version": 3, "metadata": {},
            "placements": [{"p": p, "n": [n]} for n, p in placements]}


class T(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.clades = os.path.join(self.d, "c.tsv")
        open(self.clades, "w").write(CLADES)

    def run_assign(self, placements, tree=TREE, min_lwr=0.5):
        p = os.path.join(self.d, "q.jplace")
        json.dump(jplace(placements, tree), open(p, "w"))
        return {r["query"]: r for r in up.assign(p, self.clades, min_lwr)}

    def test_match_columns(self):
        self.assertEqual(up.match_columns("abC-D..eF"), "C-DF")

    def test_inside_stem_and_sum(self):
        r = self.run_assign([
            ("q1", [[0, -1, 0.6, 0, 0.1], [1, -1, 0.3, 0, 0.1], [3, -1, 0.1, 0, 0.1]]),
            ("q2", [[2, -1, 0.9, 0, 0.2], [5, -1, 0.1, 0, 0.2]]),   # stem of X; inside Big
            ("q3", [[6, -1, 1.0, 0, 0.5]]),                          # outside everything
        ])
        self.assertEqual(r["q1"]["our_clade"], "Clade X")
        self.assertEqual(r["q1"]["our_lwr"], "0.900")
        self.assertEqual(r["q1"]["fig1_clade"], "Big")
        self.assertEqual(r["q1"]["fig1_lwr"], "1.000")
        self.assertEqual(r["q2"]["our_clade"], "unassigned")
        self.assertEqual(r["q2"]["our_stem_lwr"], "0.900")
        self.assertEqual(r["q2"]["fig1_clade"], "Big")          # stem of X lies inside Big
        self.assertEqual(r["q3"]["our_clade"], "unassigned")
        self.assertEqual(r["q3"]["fig1_clade"], "unassigned")
        self.assertEqual(r["q1"]["flags"], "")

    def test_long_pendant_flag(self):
        r = self.run_assign([("q", [[0, -1, 1.0, 0, 1.5]])])
        self.assertEqual(r["q"]["our_clade"], "Clade X")          # flagged, not rejected
        self.assertEqual(r["q"]["flags"], "long_pendant")

    def test_loo_stem_with_leaf_absent(self):
        # A pruned: the B terminal edge is now the stem of Clade X (members present: {B})
        tree = "(B:2{0},(C:1{1},D:1{2}):1{3},(E:1{4},F:1{5}):1{6}){7};"
        r = self.run_assign([("A", [[0, -1, 1.0, 0, 0.1]])], tree=tree)
        self.assertEqual(r["A"]["our_clade"], "unassigned")
        self.assertEqual(r["A"]["our_stem_lwr"], "1.000")

    def test_prune_roundtrip(self):
        root = up.prune(parse_newick("((A:1,B:1)90:1,(C:1,D:1)80:2,E:1);"), "A")
        nw = up.to_newick(root)
        self.assertEqual(sorted(l.name for l in parse_newick(nw).leaves()), ["B", "C", "D", "E"])
        self.assertIn("B:2", nw)                                   # 1 + 1 merged

    def test_prune_root_child_stays_unrooted(self):
        root = up.prune(parse_newick("((A:1,B:1)90:1,(C:1,D:1)80:2,E:1);"), "E")
        self.assertEqual(len(root.children), 3)
        self.assertEqual(sorted(l.name for l in root.leaves()), ["A", "B", "C", "D"])

    def test_prune_many_and_matched_cases(self):
        import matched_identity as mi
        d = self.d
        tree, msa = os.path.join(d, "t.nwk"), os.path.join(d, "a.fa")
        open(tree, "w").write("((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);")
        # A~B identical, C/D distinct; clades: X = {A, B}
        open(msa, "w").write(">A\nACDEFGHIKL\n>B\nACDEFGHIKL\n>C\nMNPQRSTVWY\n"
                             ">D\nMNPQRSTVWA\n>E\nACDEFMNPQR\n>F\nACDEFMNPQA\n")
        up.main(["prune", "--tree", tree, "--msa", msa, "--leaf", "A", "--remove", "B,C",
                 "--out-tree", os.path.join(d, "o.nwk"), "--out-msa", os.path.join(d, "o.fa"),
                 "--out-query", os.path.join(d, "q.fa")])
        left = sorted(l.name for l in parse_newick(open(os.path.join(d, "o.nwk")).read()).leaves())
        self.assertEqual(left, ["D", "E", "F"])
        cl = os.path.join(d, "cl.tsv")
        open(cl, "w").write("label\tour_clade\nA\tClade X\nB\tClade X\nC\t\nD\t\nE\t\nF\t\n")
        rows = {r["case"]: r for r in mi.cases(msa, cl, [None, 70.0])}
        self.assertEqual(rows["A@loo"]["possible"], "yes")        # B survives
        self.assertEqual(rows["A@70.0"]["possible"], "no")        # B (100 %) removed: impossible
        self.assertEqual(rows["A@70.0"]["nearest_surviving"], "E")
        self.assertEqual(rows["C@70.0"]["possible"], "yes")       # outside leaf: always possible

    def test_loo_summary(self):
        out = os.path.join(self.d, "a.tsv")
        up.write_rows(out, [{"query": "A", "our_clade": "Clade X", "fig1_clade": "Big"},
                            {"query": "E", "our_clade": "Clade X", "fig1_clade": "unassigned"}])
        rows = {r["query"]: r for r in up.loo([out], self.clades)}
        self.assertEqual(rows["A"]["our_ok"], "yes")
        self.assertEqual(rows["E"]["our_ok"], "no")
        self.assertEqual(rows["E"]["fig1_ok"], "yes")


    def test_compare_bfd(self):
        import compare_bfd as cb
        placed = os.path.join(self.d, "chunk_0001.tsv")
        up.write_rows(placed, [
            {"query": "p1", "our_clade": "Clade I", "fig1_clade": "Fam II", "flags": ""},
            {"query": "p2", "our_clade": "unassigned", "fig1_clade": "Basidio 1st", "flags": ""}])
        old = os.path.join(self.d, "old.tsv")
        open(old, "w").write(
            "# comment\nsample\tprotein_id\tclade\tclade_confidence\tp180_Acid\tPHYLUM\n"
            "c\tp1 x\tI\tvalidated_range\tE\tBasidiomycota\n"
            "c\tp2\tI\textrapolated_low_margin\tN\tBasidiomycota\n"
            "c\tp3\tV\textrapolated\tE\tAscomycota\n")
        merged = os.path.join(self.d, "m.tsv")
        stats, mp, mo = cb.compare(cb.read_placed([placed]), old, merged)
        self.assertEqual((mp, mo), (1, 0))                       # p3 not placed: counted
        self.assertEqual(stats["old_vs_new"][("Clade I / extrapolated", "unassigned")], 1)
        self.assertEqual(stats["oldI_fig1"][("Basidio 1st", "notE")], 1)
        rows = list(csv.DictReader(open(merged), delimiter="\t"))
        self.assertEqual(rows[0]["place_our_clade"], "Clade I")


    def test_summarize_place(self):
        import summarize_place as sp
        base = {"clade_confidence": "extrapolated", "pocket_score": "3", "catalytic_ok": "yes",
                "p39_Cys_proximal": "C", "p109_E_of_EHD": "E", "p111_D_of_EHD": "D",
                "p110_Base_EHD": "H", "PHYLUM": "Ascomycota", "LINEAGE_GROUP": "Sordariomycetes",
                "SPECIES": "X y", "place_flags": "", "place_our_stem_lwr": "0.000"}
        rows = [
            dict(base, protein_id="G1_1", clade="I", p180_Acid="A", place_our_clade="Clade I",
                 place_our_lwr="0.990", place_fig1_clade="Basidiomycetous 3rd (Family II, long)"),
            dict(base, protein_id="G1_2", clade="V", p180_Acid="E", place_our_clade="unassigned",
                 place_our_lwr="0.000", place_our_stem_lwr="0.900",
                 place_fig1_clade="Basidiomycetous 1st", place_flags="long_pendant"),
        ]
        md = sp.build(rows, "t")
        self.assertIn("| Clade I | 1 | 1 | 0 |", md)                 # placement / old / old-val
        self.assertIn("| unassigned: on the stem of a clade | 1 |", md)
        self.assertIn("**1** proteins", md)                            # long_pendant count
        self.assertIn("| Ascomycota | 1 | 0.0 % | A 1 |", md)
        self.assertIn("**1/2 (50.0 %)** proteins", md)                 # catalytic core


if __name__ == "__main__":
    unittest.main()
