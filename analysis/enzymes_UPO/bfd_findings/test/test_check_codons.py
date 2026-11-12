"""Offline tests for check_codons.py. Run: python3 -m unittest discover -s bfd_findings/test"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import check_codons as cc  # noqa: E402


class T(unittest.TestCase):
    # contig: 2 bp pad, ATG AAA | intron GTAAG | TTT TAA, then 3 bp
    CONTIG = "CC" + "ATGAAA" + "GTAAG" + "TTTTAA" + "GGG"

    def test_plus_strand_splice_with_stop_in_cds(self):
        cds, down = cc.splice(self.CONTIG, [(14, 19), (3, 8)], "1")
        self.assertEqual(cds, "ATGAAATTTTAA")
        self.assertEqual(down, "GGG")
        r = cc.assess(cds, down, "MKF")
        self.assertTrue(r["seq_start"] and r["seq_stop"] and r["stop_in_cds"])
        self.assertEqual(r["internal_stops"], 0)
        self.assertTrue(r["matches_bfd_protein"])

    def test_stop_just_downstream_counts(self):
        cds, down = cc.splice(self.CONTIG, [(3, 8), (14, 16)], "1")  # CDS stops before TAA
        self.assertEqual((cds, down), ("ATGAAATTT", "TAA"))
        r = cc.assess(cds, down)
        self.assertTrue(r["seq_stop"])
        self.assertFalse(r["stop_in_cds"])

    def test_minus_strand(self):
        rc = cc.revcomp(self.CONTIG)
        n = len(rc)
        # map plus-strand exons (3..8, 14..19) onto rc coordinates
        segs = [(n - 19 + 1, n - 14 + 1), (n - 8 + 1, n - 3 + 1)]
        cds, down = cc.splice(rc, segs, "-1")
        self.assertEqual((cds, down), ("ATGAAATTTTAA", "GGG"))

    def test_judge(self):
        self.assertEqual(cc.judge("TRUE", True), "agree")
        self.assertEqual(cc.judge("FALSE", True), "flag_says_missing_but_present")
        self.assertEqual(cc.judge("TRUE", False), "flag_says_present_but_missing")
        self.assertEqual(cc.judge(None, True), "flag_null")


if __name__ == "__main__":
    unittest.main()
