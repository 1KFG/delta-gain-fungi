"""Offline tests for tree_outline.py and name_clades.py on a toy tree."""
import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import name_clades as nc  # noqa: E402
import tree_outline as to  # noqa: E402

# Unrooted IQ-TREE style; the Mucoro pair is the outgroup. Supports are SH-aLRT/UFBoot.
TREE = ("((Z__m1__Mucoro:1,Z__m2__Mucoro:1)100/100:1,"
        "((B__a__cladeV:1,Z__x__Asco:1)99/99:1,"
        "(B__b__cladeIV:1,(Z__y__Asco:1,Z__z__Asco:1)90/100:1)90/96:1)80/80:1,"
        "Z__w__Asco:1);")
ANCHORS = ("# comment\nclade\tanchors\tpaper_mean_length\tevidence\n"
           "root\tZ__m1__Mucoro,Z__m2__Mucoro\t\t\n"
           "Muc\tZ__m1__Mucoro,Z__m2__Mucoro\t100\tt\n"
           "FamI\tZ__x__Asco,Z__y__Asco\t200\tt\n")
LENGTHS = {"Z__m1__Mucoro": 100, "Z__m2__Mucoro": 100, "B__a__cladeV": 200, "Z__x__Asco": 200,
           "B__b__cladeIV": 200, "Z__y__Asco": 200, "Z__z__Asco": 200, "Z__w__Asco": 300}


class TreeTools(unittest.TestCase):
    def test_parse_and_support(self):
        t = to.parse_newick(TREE)
        self.assertEqual(len(t.leaves()), 8)
        self.assertEqual(len(t.children), 3)          # unrooted trifurcation
        self.assertEqual(to.ufb(t.children[0]), 100)

    def test_reroot_on_outgroup(self):
        r = to.reroot(to.parse_newick(TREE), {"Z__m1__Mucoro", "Z__m2__Mucoro"})
        sides = sorted(sorted(l.name for l in c.leaves()) for c in r.children)
        self.assertIn(["Z__m1__Mucoro", "Z__m2__Mucoro"], sides)
        self.assertEqual(sum(len(s) for s in sides), 8)

    def test_ladder_order(self):
        r = to.reroot(to.parse_newick(TREE), {"Z__m1__Mucoro", "Z__m2__Mucoro"})
        lines = to.ladder(r)
        self.assertIn("Mucoro 2", lines[0])
        self.assertIn("Asco 1", lines[1])             # Z__w splits off next

    def test_name_clades(self):
        d = tempfile.mkdtemp()
        p = {k: os.path.join(d, k) for k in ("tree", "anchors", "dedup", "members", "report")}
        open(p["tree"], "w").write(TREE)
        open(p["anchors"], "w").write(ANCHORS)
        with open(p["dedup"], "w") as fh:
            fh.write("label\tlength\n" + "".join(f"{k}\t{v}\n" for k, v in LENGTHS.items()))
        nc.main(["--tree", p["tree"], "--anchors", p["anchors"], "--dedup", p["dedup"],
                 "--out-members", p["members"], "--out-report", p["report"]])
        m = {r["name"]: r for r in csv.DictReader(open(p["members"]), delimiter="\t")}
        rep = {r["clade"]: r for r in csv.DictReader(open(p["report"]), delimiter="\t")}
        # Clade V: member a, extended to (a,x) at UFB 99; stops below the UFB-80 node.
        self.assertEqual(m["x"]["our_clade"], "Clade V")
        # Clade IV: member b, extended to (b,(y,z)) at UFB 96.
        self.assertEqual(m["z"]["our_clade"], "Clade IV")
        # FamI = MRCA(x, y) = the UFB-80 node; contains both backbone clades.
        self.assertEqual(rep["FamI"]["n_leaves"], "5")
        self.assertEqual(m["b"]["fig1_clade"], "FamI")
        self.assertEqual(rep["FamI"]["status"], "ok")
        # w sits on the backbone: no Figure 1 clade, no backbone clade.
        self.assertEqual(m["w"]["fig1_clade"], "")
        self.assertEqual(m["w"]["our_clade"], "")
        self.assertEqual(rep["Muc"]["mean_length"], "100")
        self.assertEqual(rep["Muc"]["n_leaves"], "2")  # the outgroup, named from the root side

    def test_grow_skips_unsupported_intermediate(self):
        # a's MRCA with p is unsupported (50) but the clade above it is supported (99)
        tree = ("((Z__m1__Mucoro:1,Z__m2__Mucoro:1)100:1,"
                "(((B__a__cladeI:1,Z__p__Asco:1)50:1,Z__q__Asco:1)99:1,B__c__cladeV:1)70:1,"
                "Z__w__Asco:1);")
        r = to.reroot(to.parse_newick(tree), {"Z__m1__Mucoro", "Z__m2__Mucoro"})
        leaves, below = nc.index(r)
        bb = nc.backbone_clades(leaves, below, 95)
        self.assertEqual(sorted(below[id(bb["I"][0])]),
                         ["B__a__cladeI", "Z__p__Asco", "Z__q__Asco"])
        self.assertEqual(sorted(below[id(bb["V"][0])]), ["B__c__cladeV"])


if __name__ == "__main__":
    unittest.main()
