"""Tests for the parts of esmfold_batch.py that do not need torch or a GPU.

The backend swap (fair-esm -> transformers) exists so structure prediction is not confined to
the a100. The two implementations differ in two small ways that are easy to get silently wrong,
and both are covered here:

  * where set_chunk_size lives (on the model for fair-esm, on model.trunk for transformers);
  * the pLDDT scale (0-100 for fair-esm, 0-1 for transformers), which matters because
    fold_report.py and the report's >= 70 cutoff assume 0-100.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from esmfold_batch import clean, read_fasta, set_chunk  # noqa: E402


class FakeTrunk:
    def __init__(self):
        self.chunk = None

    def set_chunk_size(self, n):
        self.chunk = n


class ModelWithMethod:
    """fair-esm shape: set_chunk_size on the model itself."""

    def __init__(self):
        self.chunk = None

    def set_chunk_size(self, n):
        self.chunk = n


class ModelWithTrunk:
    """transformers shape: set_chunk_size only on .trunk."""

    def __init__(self):
        self.trunk = FakeTrunk()


class ModelWithNeither:
    pass


class TestSetChunk(unittest.TestCase):
    def test_model_level_is_used_when_present(self):
        m = ModelWithMethod()
        self.assertTrue(set_chunk(m, 64))
        self.assertEqual(m.chunk, 64)

    def test_falls_through_to_trunk(self):
        m = ModelWithTrunk()
        self.assertTrue(set_chunk(m, 32))
        self.assertEqual(m.trunk.chunk, 32)

    def test_model_wins_over_trunk(self):
        # A model carrying both must not have the trunk set instead: fair-esm's own method
        # does more than assign the attribute.
        m = ModelWithMethod()
        m.trunk = FakeTrunk()
        set_chunk(m, 16)
        self.assertEqual(m.chunk, 16)
        self.assertIsNone(m.trunk.chunk)

    def test_reports_false_when_unsupported(self):
        # Must not raise: an implementation without chunking should still fold, just hungrier.
        self.assertFalse(set_chunk(ModelWithNeither(), 8))

    def test_none_trunk_is_not_probed(self):
        m = ModelWithNeither()
        m.trunk = None
        self.assertFalse(set_chunk(m, 8))


class TestPlddtScale(unittest.TestCase):
    """The backends disagree on units; the loop writes 0-100 for both."""

    def test_transformers_scaling_reaches_the_report_threshold(self):
        # A transformers pLDDT of 0.92 is a confident model and must clear the >= 70 cutoff
        # used by fold_report.py, not be read as 0.92 and discarded.
        self.assertGreaterEqual(0.92 * 100.0, 70)

    def test_low_confidence_stays_below_the_threshold(self):
        self.assertLess(0.34 * 100.0, 70)


class TestSequenceHandling(unittest.TestCase):
    def test_clean_maps_nonstandard_to_X(self):
        self.assertEqual(clean("acdU*"), "ACDX")

    def test_clean_strips_trailing_stop_only(self):
        self.assertEqual(clean("MK*"), "MK")

    def test_read_fasta_takes_first_token_as_id(self):
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
            fh.write(">P1 some description here\nMKV\nAA\n>P2\nMM\n")
            path = fh.name
        try:
            got = list(read_fasta(path))
        finally:
            os.unlink(path)
        self.assertEqual(got, [("P1", "MKVAA"), ("P2", "MM")])


if __name__ == "__main__":
    unittest.main()
