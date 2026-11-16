"""Tests for the MSA/sequence pairing in af3_input.py.

These cover the failure that cost job 32069698 all 26 of its control predictions and job
32069699 all 426 of its monomers: an MSA attached to the wrong sequence, and no MSA attached
at all. Both produce a job JSON that AF3 accepts and then dies on in featurising, with nothing
in the error that points back here - so they are worth a test apiece.
"""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from af3_input import extract_msas, find_msa, msa_for  # noqa: E402


def write_data_json(root, name, seq, unpaired="A\nB\n"):
    """Imitate an AF3 data-pipeline output for one chain."""
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    doc = {"name": name, "sequences": [{"protein": {
        "id": "A", "sequence": seq, "unpairedMsa": unpaired, "pairedMsa": "", "templates": []}}]}
    with open(os.path.join(d, f"{name}_data.json"), "w") as fh:
        json.dump(doc, fh)
    return d


class TestMsaPairing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_matching_sequence_returns_the_msa(self):
        write_data_json(self.root, "CHAIN_X", "MKVAAA")
        got = msa_for(self.root, "CHAIN_X", "MKVAAA")
        self.assertIsNotNone(got)
        self.assertEqual(got["unpairedMsa"], "A\nB\n")

    def test_mismatched_sequence_is_refused(self):
        # The 32069698 failure: the name resolves, but the MSA was built for another chain.
        write_data_json(self.root, "CTRL_6Q76_1", "M" * 75)
        self.assertIsNone(msa_for(self.root, "CTRL_6Q76_1", "M" * 68))

    def test_no_check_when_sequence_not_supplied(self):
        # Stage 1 has no panel sequence to check against; it must still work.
        write_data_json(self.root, "CHAIN_X", "MKVAAA")
        self.assertIsNotNone(msa_for(self.root, "CHAIN_X"))

    def test_private_keys_are_not_emitted_into_the_job(self):
        # _query_seq is for the check only; AF3 would reject an unknown protein key.
        write_data_json(self.root, "CHAIN_X", "MKVAAA")
        got = msa_for(self.root, "CHAIN_X", "MKVAAA")
        self.assertFalse([k for k in got if k.startswith("_")])
        self.assertEqual(set(got), {"unpairedMsa", "pairedMsa", "templates"})

    def test_missing_directory_returns_none(self):
        self.assertIsNone(msa_for(self.root, "ABSENT", "MKV"))

    def test_find_msa_accepts_lowercase_dir(self):
        # AF3 lowercases the job directory name; the file inside may keep either spelling.
        d = os.path.join(self.root, "chain_y")
        os.makedirs(d)
        open(os.path.join(d, "chain_y_data.json"), "w").write("{}")
        self.assertIsNotNone(find_msa(self.root, "CHAIN_Y"))

    def test_extract_msas_keeps_the_query_sequence(self):
        write_data_json(self.root, "CHAIN_X", "MKVAAA")
        m = extract_msas(os.path.join(self.root, "CHAIN_X", "CHAIN_X_data.json"))
        self.assertEqual(m["A"]["_query_seq"], "MKVAAA")


class TestSingleWithMsaDir(unittest.TestCase):
    """`single` must attach MSAs when asked, and must still work without them."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.faa = os.path.join(self.root, "panel.faa")
        with open(self.faa, "w") as fh:
            fh.write(">P1\nMKVAAA\n>P2\nMMMMMM\n")
        self.msas = os.path.join(self.root, "msas")
        os.makedirs(self.msas)

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, args):
        from af3_input import main
        return main(args)

    def test_without_msa_dir_writes_bare_sequences(self):
        # Stage 1 behaviour: the data pipeline will build the MSA, so none is attached.
        out = os.path.join(self.root, "jobs1")
        self.assertEqual(self._run(["single", "--fasta", self.faa, "--out", out]), 0)
        doc = json.load(open(os.path.join(out, "P1.json")))
        prot = doc["sequences"][0]["protein"]
        self.assertNotIn("unpairedMsa", prot)

    def test_with_msa_dir_attaches_the_msa(self):
        # Stage 2 behaviour: --norun_data_pipeline, so the MSA must be in the JSON.
        write_data_json(self.msas, "P1", "MKVAAA", unpaired=">q\nMKVAAA\n")
        write_data_json(self.msas, "P2", "MMMMMM", unpaired=">q\nMMMMMM\n")
        out = os.path.join(self.root, "jobs2")
        self.assertEqual(self._run(
            ["single", "--fasta", self.faa, "--msa-dir", self.msas, "--out", out]), 0)
        prot = json.load(open(os.path.join(out, "P1.json")))["sequences"][0]["protein"]
        self.assertEqual(prot["unpairedMsa"], ">q\nMKVAAA\n")

    def test_a_chain_with_a_stale_msa_is_skipped_not_written(self):
        # Better an absent job than one that fails in featurising for an opaque reason.
        write_data_json(self.msas, "P1", "MKVAAA", unpaired=">q\nMKVAAA\n")
        write_data_json(self.msas, "P2", "DIFFERENT", unpaired=">q\nDIFFERENT\n")
        out = os.path.join(self.root, "jobs3")
        self._run(["single", "--fasta", self.faa, "--msa-dir", self.msas, "--out", out])
        self.assertTrue(os.path.exists(os.path.join(out, "P1.json")))
        self.assertFalse(os.path.exists(os.path.join(out, "P2.json")))


if __name__ == "__main__":
    unittest.main()
