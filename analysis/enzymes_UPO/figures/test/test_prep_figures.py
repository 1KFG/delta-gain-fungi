"""Offline tests for prep_figures.py. Run: python3 -m unittest discover -s figures/test"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import prep_figures as pf  # noqa: E402

MD = """Cases: 4

### characterized, matched identity (ceilings 70/55/45/35)

| identity to nearest surviving reference | n | correct | declined | wrong | precision when committing |
|---|---|---|---|---|---|
| 45-55 | 14 | 13 | 1 | 0 | 100 % |
| <35 | 16 | 14 | 2 | 0 | 100 % |
| **all** | 30 | 27 | 3 | 0 | 100 % |

### characterized, plain leave-one-out

| identity to nearest surviving reference | n | correct | declined | wrong | precision when committing |
|---|---|---|---|---|---|
| 70-100 | 10 | 10 | 0 | 0 | 100 % |
"""
SKY = ("query\ttrue_clade\tthreshold\tnearest_surviving_identity\tcall\tweight\toutcome\n"
       "A\tV\t95.0\t73.5\tV\t1.0\tcorrect\n"          # threshold 95: excluded
       "A\tV\t55.0\t50.0\tIII\t0.9\tWRONG\n"
       "B\tI\t45.0\t47.0\tI\t1.0\tcorrect\n"
       "C\tIV\t45.0\t40.0\t\t\timpossible_no_same_clade_ref\n")
S8 = ("# comment\nenzyme\tclade\tperox_pH3.5\tperox_pH5.5\tperox_pH8.0\n"
      "X\tI\t0\t0\t3\nY\tI\t51\t53\t44\nZ\tIV\t0\t3\t5\n")


def bfd_row(clade, conf, new, e="E", lwr="0.99", lineage="Sordariomycetes"):
    return {"clade": clade, "clade_confidence": conf, "place_our_clade": new,
            "place_fig1_clade": "", "place_our_lwr": lwr, "p180_Acid": e, "p39_Cys_proximal": "C",
            "p109_E_of_EHD": "E", "p111_D_of_EHD": "D", "length": "300", "LINEAGE_GROUP": lineage}


class T(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def f(self, name, txt):
        p = os.path.join(self.d, name)
        with open(p, "w") as fh:
            fh.write(txt)
        return p

    def test_bins(self):
        self.assertEqual([pf.bin_of(x) for x in (70, 69.9, 55, 45, 35, 34.9)],
                         ["70-100", "55-70", "55-70", "45-55", "35-45", "<35"])

    def test_parse_matched_md_takes_only_named_section(self):
        rows = pf.parse_matched_md(self.f("m.md", MD))
        self.assertEqual([r["bin"] for r in rows], ["45-55", "<35"])
        self.assertEqual(rows[0]["correct"], 13)

    def test_fig4_skylark_filter_and_wrong(self):
        rows = pf.fig4(self.f("m.md", MD), self.f("s.tsv", SKY))
        sky = {r["bin"]: r for r in rows if r["caller"].startswith("19")}
        self.assertEqual(set(sky), {"45-55"})                   # 95 and impossible dropped
        self.assertEqual((sky["45-55"]["correct"], sky["45-55"]["wrong"]), (1, 1))

    def test_t1_activity_is_max_over_ph_then_range(self):
        act = pf.t1_activity(self.f("s8.tsv", S8))
        self.assertEqual(act["Clade I"], "3-53 (n=2)")
        self.assertEqual(act["Clade IV"], "5-5 (n=1)")

    def test_fig3_flows_and_tiers(self):
        flows = pf.fig3([bfd_row("V", "extrapolated_low_margin", "Clade IV"),
                         bfd_row("V", "extrapolated", "Clade IV"),
                         bfd_row("", "below_identity_floor", "unassigned")])
        d = {(r["old_clade"], r["old_tier"], r["new_clade"]): r["n"] for r in flows}
        self.assertEqual(d[("Clade V", "extrapolated", "Clade IV")], 2)
        self.assertEqual(d[("unassigned", "below floor", "unassigned")], 1)

    def test_table_t1_rows(self):
        bfd = [bfd_row("I", "validated_range", "Clade I", e="A"),
               bfd_row("I", "validated_range", "Clade I"),
               bfd_row("V", "extrapolated", "unassigned", lwr="0")]
        members = [{"label": "B__x__cladeI", "set": "backbone", "our_clade": "Clade I",
                    "fig1_clade": ""}]
        sup = [{"system": "our", "clade": "Clade I", "sh_alrt": "98.5", "ufboot": "97"}]
        loo = self.f("loo.tsv", "query\tour_expected\tour_ok\nB__x__cladeI\tClade I\tyes\n")
        rows = {r["clade"]: r for r in pf.table_t1(bfd, members, sup, loo, self.f("s8.tsv", S8))}
        c1 = rows["Clade I"]
        self.assertEqual((c1["bfd_proteins"], c1["catalytic_core"], c1["glu180_pct"]), (2, 1, "50.0"))
        self.assertEqual((c1["support_sh_ufb"], c1["loo_correct"]), ("98.5/97", "1/1"))
        self.assertEqual(c1["t1_peroxygenase_pct"], "3-53 (n=2)")
        self.assertIn("outside Clades I-V", rows)


if __name__ == "__main__":
    unittest.main()
