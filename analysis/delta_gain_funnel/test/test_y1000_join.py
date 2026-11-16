"""Tests for the Y1000+ deposition join.

The three things that would silently give a wrong answer:
  * BIOPROJECT holds several accessions joined by ';'. Matching by equality misses every
    genome that belongs to an umbrella project as well, which would undercount Y1000+ and
    make the deposition effect look smaller than it is.
  * A cluster counts as novel when NO member is baseline. The "if Y1000+ were known" column
    must use the same all-members rule, not drop Y1000+ rows and re-ask.
  * The family set comes from the null run's filenames. A family with a null result but no
    regions file must be skipped, not counted as zero clusters.
"""
import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from y1000_join import family_names, family_rows, summarize, y1000_asmids  # noqa: E402

try:
    import duckdb
except ImportError:  # pragma: no cover
    duckdb = None


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


class TestFamilyNames(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_basenames_without_suffix(self):
        for n in ("pfam_PF00001", "pfam_PF00002"):
            open(os.path.join(self.tmp.name, f"{n}.ordnull.tsv"), "w").close()
        self.assertEqual(family_names(self.tmp.name), ["pfam_PF00001", "pfam_PF00002"])

    def test_prefix_selects_only_its_own(self):
        for n in ("pfam_PF1", "cazy_family_GH18"):
            open(os.path.join(self.tmp.name, f"{n}.ordnull.tsv"), "w").close()
        self.assertEqual(family_names(self.tmp.name, "cazy_"), ["cazy_family_GH18"])

    def test_empty_directory(self):
        self.assertEqual(family_names(self.tmp.name), [])


@unittest.skipIf(duckdb is None, "duckdb not installed")
class TestY1000Asmids(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.con = duckdb.connect()
        self.pq = os.path.join(self.tmp.name, "samples.parquet")
        self.con.execute(
            "create table t (ASMID varchar, BIOPROJECT varchar)")
        self.con.executemany("insert into t values (?, ?)", [
            ("G_alone", "PRJNA736342"),
            ("G_first", "PRJNA736342;PRJEB65317"),
            ("G_last", "PRJEB65317;PRJNA736342"),
            ("G_middle", "PRJEB1;PRJNA736342;PRJEB2"),
            ("G_other", "PRJNA999999"),
            ("G_substring", "PRJNA7363421"),   # must NOT match PRJNA736342
        ])
        self.con.execute(f"copy t to '{self.pq}' (format parquet)")

    def tearDown(self):
        self.tmp.cleanup()

    def test_matches_in_every_position_of_a_multi_project_field(self):
        got = y1000_asmids(self.con, self.pq, ["PRJNA736342"])
        self.assertEqual(got, {"G_alone", "G_first", "G_last", "G_middle"})

    def test_does_not_match_a_longer_accession_with_the_same_prefix(self):
        self.assertNotIn("G_substring", y1000_asmids(self.con, self.pq, ["PRJNA736342"]))

    def test_several_bioprojects_are_unioned(self):
        got = y1000_asmids(self.con, self.pq, ["PRJNA736342", "PRJNA999999"])
        self.assertIn("G_other", got)

    def test_no_bioprojects_means_empty(self):
        self.assertEqual(y1000_asmids(self.con, self.pq, []), set())


@unittest.skipIf(duckdb is None, "duckdb not installed")
class TestFamilyRows(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.reg = os.path.join(self.tmp.name, "regions")
        self.clu = os.path.join(self.tmp.name, "clusters")
        os.makedirs(self.reg)
        os.makedirs(self.clu)
        self.con = duckdb.connect()
        self.con.execute("create table y (ASMID varchar)")
        self.con.execute("insert into y values ('GY1'), ('GY2')")

    def tearDown(self):
        self.tmp.cleanup()

    def make(self, fam, members, clusters):
        write_tsv(os.path.join(self.reg, f"{fam}.members.tsv"),
                  ["region_id", "asmid", "class"], members)
        write_tsv(os.path.join(self.clu, f"{fam}.clusters.tsv"),
                  ["region_id", "c90", "c70"], clusters)

    def rows(self, fam):
        return {r["tier"]: r for r in family_rows(self.con, fam, self.reg, self.clu)}

    def test_cluster_with_a_baseline_member_is_not_novel(self):
        self.make("f", [["r1", "GB", "baseline"], ["r2", "GN", "new"]],
                  [["r1", "r1", "r1"], ["r2", "r1", "r1"]])
        r = self.rows("f")["c70"]
        self.assertEqual((r["clusters"], r["novel_now"]), (1, 0))

    def test_new_only_cluster_is_novel_and_unaffected_when_no_y1000_member(self):
        self.make("f", [["r1", "GN", "new"]], [["r1", "r1", "r1"]])
        r = self.rows("f")["c70"]
        self.assertEqual((r["novel_now"], r["novel_if_y_known"], r["y_only"]), (1, 1, 0))

    def test_cluster_seen_only_in_y1000_genomes_stops_being_novel(self):
        self.make("f", [["r1", "GY1", "new"], ["r2", "GY2", "new"]],
                  [["r1", "r1", "r1"], ["r2", "r1", "r1"]])
        r = self.rows("f")["c70"]
        self.assertEqual((r["novel_now"], r["novel_if_y_known"], r["y_only"]), (1, 0, 1))

    def test_one_y1000_member_is_enough_to_make_a_mixed_cluster_known(self):
        # The rule has to be the same one the funnel applies to baseline: a cluster is novel
        # only if NO member is already public. A cluster holding one Y1000+ protein and one
        # protein from an unrelated new genome is a cluster whose sequence has been public
        # since Y1000+ released it, so it stops being novel - it does not survive on the
        # strength of its other member. Asserting the opposite here was the one wrong
        # expectation in this file, and it would have understated the deposition effect.
        self.make("f", [["r1", "GY1", "new"], ["r2", "GN", "new"]],
                  [["r1", "r1", "r1"], ["r2", "r1", "r1"]])
        r = self.rows("f")["c70"]
        self.assertEqual((r["novel_now"], r["novel_if_y_known"], r["y_only"]), (1, 0, 1))

    def test_tiers_are_reported_separately(self):
        self.make("f", [["r1", "GN", "new"], ["r2", "GN", "new"]],
                  [["r1", "r1", "r1"], ["r2", "r2", "r1"]])
        r = self.rows("f")
        self.assertEqual(r["c90"]["clusters"], 2)
        self.assertEqual(r["c70"]["clusters"], 1)

    def test_absent_tier_column_is_skipped_not_invented(self):
        self.make("f", [["r1", "GN", "new"]], [["r1", "r1", "r1"]])
        self.assertEqual(set(self.rows("f")), {"c90", "c70"})

    def test_missing_family_files_give_no_rows(self):
        self.assertEqual(family_rows(self.con, "absent", self.reg, self.clu), [])


class TestSummarize(unittest.TestCase):
    def test_pools_one_tier_only(self):
        rows = [
            {"tier": "c70", "clusters": 10, "novel_now": 6, "novel_if_y_known": 4},
            {"tier": "c70", "clusters": 10, "novel_now": 4, "novel_if_y_known": 4},
            {"tier": "c50", "clusters": 99, "novel_now": 99, "novel_if_y_known": 0},
        ]
        s = summarize(rows, "c70")
        self.assertEqual((s["families"], s["clusters"], s["novel_now"]), (2, 20, 10))
        self.assertEqual(s["y_only"], 2)
        self.assertAlmostEqual(s["pct_of_novel_lost"], 20.0)

    def test_no_novel_clusters_does_not_divide_by_zero(self):
        s = summarize([{"tier": "c70", "clusters": 5, "novel_now": 0,
                        "novel_if_y_known": 0}], "c70")
        self.assertEqual(s["pct_of_novel_lost"], 0.0)

    def test_absent_tier_is_empty(self):
        self.assertEqual(summarize([], "c70")["families"], 0)


if __name__ == "__main__":
    unittest.main()
