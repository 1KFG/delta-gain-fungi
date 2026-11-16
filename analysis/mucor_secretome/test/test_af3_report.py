"""Tests for which predictions af3_report.py counts as results.

Two failures these pin down, both from the first completed chitin run:
  * AF3 writes a summary for the prediction AND one inside every seed-N_sample-M/ subdirectory.
    Walking the tree gave six rows per job, five labelled with the sample directory name instead
    of the job - 1,010 of 1,212 rows - which then sorted into the "best" list unattributable.
  * A prediction directory from an earlier panel is not a result. Six CTRL_4EBZ_1__* directories
    survived the control rename and were reported beside their renamed twins with identical
    scores.
"""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from af3_report import find_summaries, job_names, rows_from  # noqa: E402

SUMMARY = {"ptm": 0.8, "iptm": 0.7, "ranking_score": 0.9,
           "fraction_disordered": 0.1, "has_clash": 0.0,
           "chain_pair_iptm": [[0.8, 0.7], [0.7, 0.8]]}


def make_pred(root, job, with_samples=True):
    d = os.path.join(root, job)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{job}_summary_confidences.json"), "w") as fh:
        json.dump(SUMMARY, fh)
    if with_samples:
        for i in range(5):
            sd = os.path.join(d, f"seed-1_sample-{i}")
            os.makedirs(sd, exist_ok=True)
            with open(os.path.join(sd, f"seed-1_sample-{i}_summary_confidences.json"), "w") as fh:
                json.dump(SUMMARY, fh)
    return d


def make_jobs(root, names):
    d = os.path.join(root, "pairs")
    os.makedirs(d, exist_ok=True)
    for n in names:
        open(os.path.join(d, f"{n}.json"), "w").write("{}")
    return d


class TestFindSummaries(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.pred = os.path.join(self.tmp.name, "pred")
        os.makedirs(self.pred)

    def tearDown(self):
        self.tmp.cleanup()

    def test_one_row_per_job_not_one_per_sample(self):
        make_pred(self.pred, "A__B")
        got = find_summaries(self.pred)
        self.assertEqual([n for n, _ in got], ["A__B"])

    def test_sample_directories_never_become_job_names(self):
        make_pred(self.pred, "A__B")
        names = [n for n, _ in find_summaries(self.pred)]
        self.assertFalse([n for n in names if "sample" in n])

    def test_job_list_excludes_a_stale_directory(self):
        make_pred(self.pred, "CTRL_4EBZ_03a7bf__X")
        make_pred(self.pred, "CTRL_4EBZ_1__X")          # left over from the old naming
        jobs = {"CTRL_4EBZ_03a7bf__X"}
        got = [n for n, _ in find_summaries(self.pred, jobs)]
        self.assertEqual(got, ["CTRL_4EBZ_03a7bf__X"])

    def test_without_a_job_list_everything_top_level_is_kept(self):
        make_pred(self.pred, "A__B")
        make_pred(self.pred, "C__D")
        self.assertEqual(len(find_summaries(self.pred)), 2)

    def test_lowercase_summary_filename_is_accepted(self):
        d = os.path.join(self.pred, "A__B")
        os.makedirs(d)
        with open(os.path.join(d, "a__b_summary_confidences.json"), "w") as fh:
            json.dump(SUMMARY, fh)
        self.assertEqual([n for n, _ in find_summaries(self.pred)], ["A__B"])

    def test_directory_without_a_summary_is_skipped(self):
        os.makedirs(os.path.join(self.pred, "UNFINISHED"))
        self.assertEqual(find_summaries(self.pred), [])


class TestJobNames(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_reads_json_basenames(self):
        d = make_jobs(self.tmp.name, ["A__B", "C__CO4"])
        self.assertEqual(job_names(d), {"A__B", "C__CO4"})

    def test_missing_directory_means_no_filtering(self):
        self.assertIsNone(job_names(os.path.join(self.tmp.name, "nope")))

    def test_none_means_no_filtering(self):
        self.assertIsNone(job_names(None))


class TestRows(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.pred = os.path.join(self.tmp.name, "pred")
        os.makedirs(self.pred)

    def tearDown(self):
        self.tmp.cleanup()

    def test_row_count_matches_job_count(self):
        make_pred(self.pred, "A__B")
        make_pred(self.pred, "C__D")
        self.assertEqual(len(rows_from(self.pred)), 2)

    def test_control_role_is_carried(self):
        make_pred(self.pred, "A__B")
        rows = rows_from(self.pred, {"A__B": "positive"})
        self.assertEqual(rows[0]["kind"], "control")
        self.assertEqual(rows[0]["control_role"], "positive")

    def test_stale_rows_are_dropped_when_a_job_list_is_given(self):
        make_pred(self.pred, "A__B")
        make_pred(self.pred, "OLD__B")
        rows = rows_from(self.pred, None, {"A__B"})
        self.assertEqual([r["job"] for r in rows], ["A__B"])


if __name__ == "__main__":
    unittest.main()
