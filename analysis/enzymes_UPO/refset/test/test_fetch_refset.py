"""Offline tests for fetch_refset.py. No network: every HTTP response comes from a fake fetcher.

Fixtures are trimmed from real responses (see the comment at the top of each).
Run:  /usr/bin/python3.12 -m unittest discover -s analysis/enzymes_UPO/refset/test -v
"""
import datetime as dt
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bin"))
import fetch_refset as fr  # noqa: E402

FIX = HERE / "fixtures"
PB_PAGE = (FIX / "peroxibase_10810.html").read_text()            # entry page with inline sequence
PB_FASTA = (FIX / "peroxibase_get_fasta_12616.html").read_text()  # archived FASTA endpoint
RCSB = (FIX / "rcsb_5FUJ.fasta").read_text()


class FakeFetcher:
    """Answers by URL substring (first match in insertion order); records every URL."""

    def __init__(self, routes):
        self.routes, self.calls = dict(routes), []

    def get(self, url):
        self.calls.append(url)
        for key, body in self.routes.items():
            if key in url:
                if isinstance(body, Exception):
                    raise body
                return body
        raise RuntimeError(f"unexpected URL {url}")


def row(abbrev, route, resolved, length, given=None):
    return {"abbrev": abbrev, "given_id": given or resolved, "route": route,
            "resolved_id": resolved, "length_aa": str(length)}


def cdx(*pairs):
    return json.dumps([["timestamp", "original"]] + [list(p) for p in pairs])


def enc(s):
    return s.replace("/", "%2F")


GF = "tools/get_fasta/{}/PEP"
VP = "display_perox/view_perox/{}"
OLD, NEW = fr.PB_HOSTS


class Parsers(unittest.TestCase):
    def test_peroxibase_entry_page(self):
        info = fr.parse_peroxibase_page(PB_PAGE)
        self.assertEqual(len(info["seq"]), 371)                 # Table S1: AbHTP5, 371 aa
        self.assertTrue(info["seq"].startswith("MNPFLKLAV"))
        self.assertNotIn("ATGAATCC", info["seq"])               # not the GenSeq (DNA) field
        self.assertEqual(info["pb_entry"], "AbHalPrx05_bisporus")
        self.assertEqual(info["pb_class"], "Haloperoxidase (haem)")
        self.assertEqual(info["pb_orthogroup"], "HalPrx001")

    def test_peroxibase_entry_page_without_sequence(self):
        self.assertEqual(fr.parse_peroxibase_page("<html>Entry information : X </html>")["seq"], "")

    def test_peroxibase_fasta_endpoint(self):
        info = fr.parse_peroxibase_fasta(PB_FASTA)
        self.assertEqual(len(info["seq"]), 288)                 # Table S1: SbreHTP1, 288 aa
        self.assertTrue(info["seq"].startswith("MGFFSSVSRLF"))
        self.assertTrue(info["seq"].endswith("ATST"))           # trailing stop '*' removed
        self.assertEqual(info["pb_entry"], "SbreHalPrx01")
        self.assertEqual(info["pb_class"], "Haloperoxidase (haem)")
        self.assertEqual(fr.parse_peroxibase_fasta("<html>nothing</html>")["seq"], "")

    def test_pdb_chain(self):
        self.assertTrue(fr.pick_pdb_chain(RCSB, "A").startswith("SAHPWKAPGPND"))
        self.assertEqual(fr.pick_pdb_chain(RCSB, "Z"), "")

    def test_fasta_and_uniprot_header(self):
        recs = fr.parse_fasta(">sp|P04963|PRXC_LEPFU x\nmkl*\nAC\n>tr|Q9|Y z\nGG\n")
        self.assertEqual(recs, [("sp|P04963|PRXC_LEPFU x", "MKLAC"), ("tr|Q9|Y z", "GG")])
        self.assertEqual(fr.uniprot_accession(recs[0][0]), "P04963")

    def test_closest_capture_to_access_date(self):
        caps = [("20240807005644", "a"), ("20220518034748", "b"), ("20210119214828", "c")]
        self.assertEqual([c[1] for c in fr.closest_first(caps, dt.datetime(2022, 4, 7))], ["b", "c", "a"])


class Run(unittest.TestCase):
    def setUp(self):
        self.fake = FakeFetcher([
            # UniProtKB batch: P9 (deleted) is absent, like the real deleted entries
            ("uniprotkb/accessions", ">sp|P1|A\nMAAAA\n>tr|P2|B\nMCC\n>tr|P3|C\nMCC\n>tr|P4|D\nMEEEEEE\n"),
            ("uniprotkb/P9.json", json.dumps({
                "entryType": "Inactive", "primaryAccession": "P9",
                "inactiveReason": {"inactiveReasonType": "DELETED",
                                   "deletedReason": "Not part of a reference proteome"},
                "extraAttributes": {"uniParcId": "UPI0000000009"}})),
            ("uniparc/UPI0000000009.fasta", ">UPI0000000009 status=active\nMKKKK\n"),
            ("efetch.fcgi", ">XP_1.1 something\nMDDDDD\n"),
            # PeroxiBase 12616: only the old domain's FASTA endpoint was captured
            (enc(f"{OLD}/" + GF.format(12616)),
             cdx(("20230604012120", "https://x/12616"), ("20220518034748", "https://x/12616"))),
            ("20220518034748id_/https://x/12616", PB_FASTA),
            # PeroxiBase 10810: no FASTA capture; the entry page has the sequence
            (enc(f"{OLD}/" + VP.format(10810)), cdx(("20240525203944", "https://y/10810"))),
            ("20240525203944id_/https://y/10810", PB_PAGE),
            # PeroxiBase 4242 (xref row): FASTA endpoint captured, 288 aa
            (enc(f"{OLD}/" + GF.format(4242)), cdx(("20210119214828", "https://z/4242"))),
            ("20210119214828id_/https://z/4242", PB_FASTA),
            ("cdx/search/cdx", ""),                              # every other lookup: no captures
            ("rcsb.org/fasta/entry/5FUJ", RCSB),
        ])

    def run_rows(self, rows):
        return fr.run(rows, self.fake, log=io.StringIO())

    def test_every_route(self):
        rows = [
            row("U1", "uniprot_direct", "P1", 5),
            row("X2", "uniprot_xref", "P2;P3", 3, given="4207"),
            row("N1", "ncbi_efetch", "XP_1", 7),                  # length disagrees on purpose
            row("W1", "wayback_view_perox", "peroxibase:12616", 288, given="12616"),
            row("W2", "wayback_view_perox", "peroxibase:10810", 371, given="10810"),
            row("W3", "wayback_view_perox", "peroxibase:4233", 451, given="4233"),
            row("I1", "inferred_unique_PF01328_len_match", "P1", 5),
            row("D1", "pdb (MroUPO)", "5FUJ_A", 236),
        ]
        recs, man = self.run_rows(rows)
        by = {m["abbrev"]: m for m in man}
        self.assertEqual(len(man), len(rows))                      # nothing dropped
        self.assertEqual([n for n, _ in recs], ["U1", "X2", "N1", "W1", "W2", "I1", "D1"])
        self.assertEqual((by["U1"]["source"], by["U1"]["len_match"]), ("uniprotkb", "yes"))
        self.assertIn("identical", by["X2"]["note"])
        self.assertEqual(by["N1"]["len_match"], "no")
        self.assertEqual((by["W1"]["source"], by["W1"]["fetched_len"]), ("peroxibase_wayback", "288"))
        self.assertIn("20220518034748id_", by["W1"]["source_url"])   # capture closest to 2022-04-07
        self.assertEqual(by["W1"]["pb_class"], "Haloperoxidase (haem)")
        self.assertEqual((by["W2"]["fetched_len"], by["W2"]["pb_orthogroup"]), ("371", "HalPrx001"))
        self.assertEqual(by["W3"]["status"], "missing")
        self.assertIn("0 captures", by["W3"]["note"])
        self.assertIn("inferred", by["I1"]["note"])
        self.assertEqual(by["D1"]["source"], "rcsb")
        up_batch = [c for c in self.fake.calls if "uniprotkb/accessions" in c]
        self.assertEqual(len(up_batch), 1)                         # one batched request
        self.assertIn("P1%2CP2%2CP3", up_batch[0])

    def test_deleted_uniprot_entry_falls_back_to_uniparc(self):
        _, man = self.run_rows([row("G1", "uniprot_direct", "P9", 5)])
        self.assertEqual((man[0]["status"], man[0]["source"]), ("ok", "uniparc"))
        self.assertIn("DELETED: Not part of a reference proteome", man[0]["note"])
        self.assertIn("UPI0000000009", man[0]["note"])

    def test_xref_length_mismatch_prefers_peroxibase_match(self):
        # UniProt P4 is 7 aa; Table S1 says 288; the archived PeroxiBase FASTA is 288 -> keep it
        _, man = self.run_rows([row("A2", "uniprot_xref", "P4", 288, given="4242")])
        self.assertEqual((man[0]["source"], man[0]["len_match"]), ("peroxibase_wayback", "yes"))
        self.assertIn("uniprotkb=7, peroxibase_wayback=288; kept peroxibase_wayback", man[0]["note"])

    def test_xref_length_match_skips_wayback(self):
        self.run_rows([row("U1", "uniprot_xref", "P1", 5, given="4242")])
        self.assertFalse(any("cdx" in c for c in self.fake.calls))

    def test_no_source_matches_keeps_first_and_flags(self):
        _, man = self.run_rows([row("A3", "uniprot_xref", "P4", 999, given="4242")])
        self.assertEqual((man[0]["source"], man[0]["len_match"]), ("uniprotkb", "no"))

    def test_differing_duplicate_accessions_are_flagged(self):
        self.fake.routes["uniprotkb/accessions"] = ">tr|P2|B\nMCC\n>tr|P3|C\nMCD\n"
        _, man = self.run_rows([row("X2", "uniprot_xref", "P2;P3", 3)])
        self.assertIn("DIFFER", man[0]["note"])

    def test_errors_are_recorded_not_raised(self):
        self.fake.routes = {"uniprotkb/accessions": "", "rcsb.org": RuntimeError("boom")}
        _, man = self.run_rows([row("D1", "pdb x", "9ZZZ_A", 100)])
        self.assertEqual(man[0]["status"], "missing")
        self.assertIn("error: boom", man[0]["note"])

    def test_write_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            fr.write_outputs(Path(d), [("A", "M" * 130)], [{c: "" for c in fr.MANIFEST_COLS}])
            faa = (Path(d) / "zamocky2022_refs.faa").read_text().splitlines()
            self.assertEqual(faa[0], ">A")
            self.assertEqual([len(x) for x in faa[1:]], [60, 60, 10])
            hdr = (Path(d) / "fetch_manifest.tsv").read_text().splitlines()[0]
            self.assertEqual(hdr.split("\t"), fr.MANIFEST_COLS)


class FetcherCache(unittest.TestCase):
    def test_cache_hit_needs_no_network(self):
        with tempfile.TemporaryDirectory() as d:
            f = fr.Fetcher(cache=Path(d))
            url = "https://example.invalid/x"
            key = Path(d) / (fr.hashlib.sha1(url.encode()).hexdigest() + ".txt")
            key.write_text("cached")
            self.assertEqual(f.get(url), "cached")             # served from cache, no network


if __name__ == "__main__":
    unittest.main()
