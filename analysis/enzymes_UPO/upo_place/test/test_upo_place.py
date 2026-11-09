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

    def test_loo_summary(self):
        out = os.path.join(self.d, "a.tsv")
        up.write_rows(out, [{"query": "A", "our_clade": "Clade X", "fig1_clade": "Big"},
                            {"query": "E", "our_clade": "Clade X", "fig1_clade": "unassigned"}])
        rows = {r["query"]: r for r in up.loo([out], self.clades)}
        self.assertEqual(rows["A"]["our_ok"], "yes")
        self.assertEqual(rows["E"]["our_ok"], "no")
        self.assertEqual(rows["E"]["fig1_ok"], "yes")


if __name__ == "__main__":
    unittest.main()
