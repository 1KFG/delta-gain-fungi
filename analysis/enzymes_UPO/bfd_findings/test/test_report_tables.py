"""Offline tests for report_tables.py. Run: python3 -m unittest discover -s bfd_findings/test"""
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import report_tables as rt  # noqa: E402


class T(unittest.TestCase):
    def test_category(self):
        self.assertEqual(rt.category("TRUE", "TRUE"), "complete")
        self.assertEqual(rt.category("FALSE", "TRUE"), "no_start")
        self.assertEqual(rt.category("true", "false"), "no_stop")
        self.assertEqual(rt.category("FALSE", "FALSE"), "neither")
        self.assertEqual(rt.category("", "TRUE"), "flags_null")
        self.assertEqual(rt.category(None, None), "flags_null")

    def test_end_distance(self):
        # 1-based inclusive: a model at 1..300 touches the left end
        self.assertEqual(rt.end_distance(1, 300, 10000), 0)
        self.assertEqual(rt.end_distance(5001, 6000, 10000), 4000)
        self.assertEqual(rt.end_distance(9001, 9950, 10000), 50)
        self.assertIsNone(rt.end_distance(1, 2, None))

    def test_contig_lengths_only_wanted(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "g.fa.gz")
            with gzip.open(p, "wt") as fh:
                fh.write(">c1 desc\nACGT\nAC\n>c2\nAAAAAAAAAA\n>c3\nA\n")
            path, got = rt.contig_lengths((p, {"c1", "c3"}))
            self.assertEqual(got, {"c1": 6, "c3": 1})

    def test_quartiles_and_pct(self):
        self.assertEqual(rt.pct(1, 4), "25.0")
        self.assertEqual(rt.pct(1, 0), "NA")
        self.assertEqual(rt.quartiles([]), ("NA", "NA", "NA"))
        self.assertEqual(rt.quartiles([1.0, 2.0, 3.0]), ("1.000", "2.000", "3.000"))


if __name__ == "__main__":
    unittest.main()
