"""Offline tests for build_refset.py (combine + select). Run:
    python3 -m unittest discover -s analysis/enzymes_UPO/refset/test
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import build_refset as br  # noqa: E402

BACKBONE = """>HspUPO|cladeV|Short|OTA57433.1
MKLLAAACDEFGHIK
>MroUPO|cladeIII|Short|X
MPPPQQQRRRSSS
"""
# Zm1 = HspUPO with a signal peptide (contains it) -> duplicate of the backbone
# Zm2 = new; Zm3 = identical to Zm2 -> duplicate of Zm2; Zm4 = no PF01328 hit in sites
ZAM = """>Zm1
MSSSMKLLAAACDEFGHIK
>Zm2
MWWWYYYVVV
>Zm3
MWWWYYYVVV*
>Zm4
MNNNNNNN
"""
TABLE = ("abbrev\tgiven_id\tid_type\tresolved_id\troute\torganism\ttaxon\tlength_aa\tpdb\tnote\n"
         "Zm1\tA\tu\tA\tu\to\tAscomycota \t19\t\t\n"
         "Zm2\tB\tu\tB\tu\to\tAmobeozoa\t12\t\t\n"
         "Zm3\tC\tu\tC\tu\to\tMucoromycotina\t10\t\t\n"
         "Zm4\tD\tu\tD\tu\to\tBasidiomycota\t8\t\t\n")


def sites_tsv(rows):
    head = ("# comment\nsample\tprotein_id\thmm_evalue\tlength\tclade\tclade_identity\t"
            "clade_margin\tclade_confidence\tnearest_ref\tcatalytic_ok\tcatalytic_detail\n")
    return head + "".join(
        f"s\t{pid}\t1e-50\t10\t{cl}\t60.0\t5\tok\tHspUPO\t{ok}\t39C\n" for pid, cl, ok in rows)


class T(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.p = {}
        for k, txt in {"bb": BACKBONE, "zam": ZAM, "tab": TABLE}.items():
            self.p[k] = os.path.join(self.d, k)
            with open(self.p[k], "w") as fh:
                fh.write(txt)

    def test_taxon_variants(self):
        self.assertEqual(br.taxon_group("Ascomycota "), "Asco")
        self.assertEqual(br.taxon_group("Amobeozoa"), "Amoebozoa")
        self.assertEqual(br.taxon_group("Zygomycota_Mucoromycotina"), "Mucoro")
        with self.assertRaises(ValueError):
            br.taxon_group("Plantae")

    def test_combine_dedup(self):
        recs, rows = br.combine(self.p["zam"], self.p["tab"], self.p["bb"])
        labels = [r[0] for r in recs]
        self.assertEqual(labels, ["B__HspUPO__cladeV", "B__MroUPO__cladeIII",
                                  "Z__Zm2__Amoebozoa", "Z__Zm4__Basidio"])
        by = {r["name"]: r for r in rows}
        self.assertEqual(by["Zm1"]["dup_of"], "B__HspUPO__cladeV")   # containment
        self.assertEqual(by["HspUPO"]["dup_of"], "Z__Zm1__Asco")    # reverse link
        self.assertEqual(by["Zm3"]["dup_of"], "Z__Zm2__Amoebozoa")  # stop codon stripped
        self.assertEqual(by["Zm2"]["dup_of"], "")

    def test_select(self):
        recs, rows = br.combine(self.p["zam"], self.p["tab"], self.p["bb"])
        comb, dd, st = (os.path.join(self.d, x) for x in ("c.faa", "d.tsv", "s.tsv"))
        br.write_fasta(comb, recs)
        br.write_tsv(dd, rows)
        with open(st, "w") as fh:
            fh.write(sites_tsv([("B__HspUPO__cladeV", "V", "yes"),
                                ("B__MroUPO__cladeIII", "III", "yes"),
                                ("Z__Zm2__Amoebozoa", "II", "no")]))
        qc, tree = br.select(comb, dd, st, self.p["tab"])
        self.assertEqual([t[0] for t in tree],
                         ["B__HspUPO__cladeV", "B__MroUPO__cladeIII", "Z__Zm2__Amoebozoa"])
        by = {r["name"]: r for r in qc}
        self.assertEqual(by["Zm4"]["reason"], "no PF01328 hit")
        self.assertEqual(by["Zm4"]["in_tree"], "no")
        self.assertTrue(by["Zm3"]["reason"].startswith("duplicate"))
        self.assertIn("catalytic_residues_differ", by["Zm2"]["reason"])
        self.assertIn("length_10_vs_tableS1_12", by["Zm2"]["reason"])


if __name__ == "__main__":
    unittest.main()
