"""Offline tests for build_secretome.py and extract_proteins.py.
Run: python3 -m unittest discover -s mucor_secretome/test"""
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import build_secretome as bs  # noqa: E402
import extract_proteins as ep  # noqa: E402


class Classify(unittest.TestCase):
    def test_topology_parsing(self):
        self.assertEqual(bs.parse_topology("i7-29o44-66i"), [(7, 29), (44, 66)])
        self.assertEqual(bs.parse_topology(""), [])

    def test_signal_peptide_helix_is_not_a_membrane_helix(self):
        # the only helix sits inside the signal peptide (cleaved at 22)
        self.assertEqual(bs.membrane_helices("n5-20o", 1, 22), 0)
        # a helix well past the cleavage site is a real TM
        self.assertEqual(bs.membrane_helices("i60-82o", 1, 22), 1)

    def test_membrane_helices_without_topology_falls_back_to_count(self):
        self.assertEqual(bs.membrane_helices("", 3, 22), 2)
        self.assertEqual(bs.membrane_helices(None, 1, 0), 1)

    def test_classify(self):
        self.assertEqual(bs.classify(0.99, "n5-20o", 1, 22, False, 0.0), "secreted")
        self.assertEqual(bs.classify(0.6, "", 0, 20, False, 0.0), "weak_signal")
        self.assertEqual(bs.classify(0.99, "i60-82o", 1, 22, False, 0.0), "membrane")
        self.assertEqual(bs.classify(0.99, "", 0, 22, True, 0.0), "gpi_anchored")
        self.assertEqual(bs.classify(0.99, "", 0, 22, False, 0.9), "mitochondrial")
        # a membrane helix outranks a GPI anchor
        self.assertEqual(bs.classify(0.99, "i60-82o", 1, 22, True, 0.0), "membrane")


class Tiers(unittest.TestCase):
    def test_chitin_beats_plant(self):
        t, why = bs.tier_of(["CE4(145-269)"], ["chitin"], [], 300, 5, 0.02)
        self.assertEqual(t, "T1")
        self.assertIn("CE4", why)

    def test_lysm_pfam_is_t1(self):
        self.assertEqual(bs.tier_of([], [], ["PF01476"], 200, 0, 0.0)[0], "T1")

    def test_plant_substrate_is_t2(self):
        t, why = bs.tier_of(["GH28(20-300)"], ["pectin, polygalacturonan"], [], 350, 2, 0.01)
        self.assertEqual(t, "T2")
        self.assertIn("pectin", why)

    def test_upo_pfam_is_t2(self):
        self.assertEqual(bs.tier_of([], [], ["PF01328"], 280, 6, 0.03)[0], "T2")

    def test_effector_like_needs_no_domain(self):
        self.assertEqual(bs.tier_of([], [], [], 180, 8, 0.045)[0], "T3")
        # same protein with a Pfam domain is not "unknown"
        self.assertEqual(bs.tier_of([], [], ["PF99999"], 180, 8, 0.045)[0], "T4")
        # too long
        self.assertEqual(bs.tier_of([], [], [], 400, 8, 0.02)[0], "T4")
        # too few cysteines
        self.assertEqual(bs.tier_of([], [], [], 180, 3, 0.05)[0], "T4")

    def test_cazy_root(self):
        self.assertEqual(bs.cazy_root("CE4(145-269)"), "CE4")
        self.assertEqual(bs.cazy_root("GH18_3"), "GH18")
        self.assertEqual(bs.cazy_root(""), "")

    def test_cys_stats_ignores_trailing_stop(self):
        n, frac, length = bs.cys_stats("MCCAA*")
        self.assertEqual((n, length), (2, 5))
        self.assertAlmostEqual(frac, 0.4)


class Extract(unittest.TestCase):
    def test_read_fasta_filters_and_keeps_last_record(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "p.fa.gz")
            with gzip.open(p, "wt") as fh:
                fh.write(">A_1-T1 gene1\nMKV\nAA\n>A_2-T1 gene2\nMMM\n")
            self.assertEqual(dict(ep.read_fasta(p)), {"A_1-T1": "MKVAA", "A_2-T1": "MMM"})
            self.assertEqual(dict(ep.read_fasta(p, {"A_2-T1"})), {"A_2-T1": "MMM"})

    def test_read_index(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "genomes.tsv")
            with open(p, "w") as fh:
                fh.write("LT1\tout1\t/x/p.fa.gz\t/x/c.fa.gz\n\nLT2\tout2\t/y/p.fa.gz\n")
            self.assertEqual(ep.read_index(p), {"LT1": "/x/p.fa.gz", "LT2": "/y/p.fa.gz"})


if __name__ == "__main__":
    unittest.main()
