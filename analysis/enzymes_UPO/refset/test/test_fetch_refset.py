"""Offline tests for fetch_refset.py. No network: every HTTP response comes from a fake fetcher.

Run:  /usr/bin/python3.12 -m unittest discover -s analysis/enzymes_UPO/refset/test -v
"""
import io
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bin"))
import fetch_refset as fr  # noqa: E402

FIX = HERE / "fixtures"
PB_PAGE = (FIX / "peroxibase_10810.html").read_text()
RCSB = (FIX / "rcsb_5FUJ.fasta").read_text()


class FakeFetcher:
    """Answers by URL substring; records every URL asked for."""

    def __init__(self, routes):
        self.routes, self.calls = routes, []

    def get(self, url):
        self.calls.append(url)
        for key, body in self.routes.items():
            if key in url:
                return body
        raise RuntimeError(f"unexpected URL {url}")


def row(abbrev, route, resolved, length, given=None):
    return {"abbrev": abbrev, "given_id": given or resolved, "route": route,
            "resolved_id": resolved, "length_aa": str(length)}


AVAIL_OK = json.dumps({"archived_snapshots": {"closest": {
    "available": True, "status": "200", "timestamp": "20240525203944",
    "url": "http://web.archive.org/web/20240525203944/https://peroxibase.toulouse.inra.fr/"
           "display_perox/view_perox/10810"}}})
AVAIL_NONE = json.dumps({"archived_snapshots": {}})


class Parsers(unittest.TestCase):
    def test_peroxibase_page(self):
        info = fr.parse_peroxibase_page(PB_PAGE)
        self.assertEqual(len(info["seq"]), 371)          # Table S1 length for AbHTP5
        self.assertTrue(info["seq"].startswith("MNPFLKLAV"))
        self.assertNotIn("ATGAATCC", info["seq"])        # not the GenSeq (DNA) field
        self.assertEqual(info["pb_entry"], "AbHalPrx05_bisporus")
        self.assertEqual(info["pb_class"], "Haloperoxidase (haem)")
        self.assertEqual(info["pb_orthogroup"], "HalPrx001")

    def test_pdb_chain(self):
        seq = fr.pick_pdb_chain(RCSB, "A")
        self.assertTrue(seq.startswith("SAHPWKAPGPND"))
        self.assertEqual(fr.pick_pdb_chain(RCSB, "Z"), "")

    def test_fasta_and_uniprot_header(self):
        recs = fr.parse_fasta(">sp|P04963|PRXC_LEPFU x\nmkl*\nAC\n>tr|Q9|Y z\nGG\n")
        self.assertEqual(recs, [("sp|P04963|PRXC_LEPFU x", "MKLAC"), ("tr|Q9|Y z", "GG")])
        self.assertEqual(fr.uniprot_accession(recs[0][0]), "P04963")


class Run(unittest.TestCase):
    def setUp(self):
        self.fake = FakeFetcher({
            "rest.uniprot.org": ">sp|P1|A\nMAAAA\n>tr|P2|B\nMCC\n>tr|P3|C\nMCC\n",
            "efetch.fcgi": ">XP_1.1 something\nMDDDDD\n",
            "wayback/available?url=peroxibase.toulouse.inra.fr%2Fdisplay_perox%2Fview_perox%2F10810": AVAIL_OK,
            "wayback/available": AVAIL_NONE,
            "20240525203944id_/https://peroxibase.toulouse.inra.fr/display_perox/view_perox/10810": PB_PAGE,
            "rcsb.org/fasta/entry/5FUJ": RCSB,
        })

    def test_every_route_and_manifest(self):
        rows = [
            row("U1", "uniprot_direct", "P1", 5),
            row("X2", "uniprot_xref", "P2;P3", 3, given="4207"),
            row("N1", "ncbi_efetch", "XP_1", 7),                    # length disagrees on purpose
            row("W1", "wayback_view_perox", "peroxibase:10810", 371, given="10810"),
            row("W2", "wayback_view_perox", "peroxibase:4233", 451, given="4233"),
            row("I1", "inferred_unique_PF01328_len_match", "P1", 5),
            row("D1", "pdb (MroUPO)", "5FUJ_A", 236),
        ]
        recs, man = fr.run(rows, self.fake, log=io.StringIO())
        by = {m["abbrev"]: m for m in man}

        self.assertEqual(len(man), len(rows))                     # nothing dropped
        self.assertEqual([n for n, _ in recs], ["U1", "X2", "N1", "W1", "I1", "D1"])
        self.assertEqual(by["U1"]["len_match"], "yes")
        self.assertIn("identical", by["X2"]["note"])
        self.assertEqual(by["N1"]["len_match"], "no")
        self.assertEqual(by["W1"]["fetched_len"], "371")
        self.assertEqual(by["W1"]["pb_orthogroup"], "HalPrx001")
        self.assertTrue(by["W1"]["source_url"].startswith("https://web.archive.org/web/20240525203944id_/"))
        self.assertEqual(by["W2"]["status"], "missing")
        self.assertIn("no Wayback capture", by["W2"]["note"])
        self.assertIn("inferred", by["I1"]["note"])
        self.assertEqual(by["D1"]["status"], "ok")
        # UniProt asked once for the de-duplicated accession set
        up_calls = [c for c in self.fake.calls if "rest.uniprot.org" in c]
        self.assertEqual(len(up_calls), 1)
        self.assertIn("P1%2CP2%2CP3", up_calls[0])

    def test_differing_duplicate_accessions_are_flagged(self):
        self.fake.routes["rest.uniprot.org"] = ">tr|P2|B\nMCC\n>tr|P3|C\nMCD\n"
        _, man = fr.run([row("X2", "uniprot_xref", "P2;P3", 3)], self.fake, log=io.StringIO())
        self.assertIn("DIFFER", man[0]["note"])

    def test_errors_are_recorded_not_raised(self):
        rows = [row("D1", "pdb x", "9ZZZ_A", 100)]
        _, man = fr.run(rows, self.fake, log=io.StringIO())
        self.assertEqual(man[0]["status"], "missing")
        self.assertIn("error:", man[0]["note"])

    def test_write_outputs(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            fr.write_outputs(Path(d), [("A", "M" * 130)], [{c: "" for c in fr.MANIFEST_COLS}])
            faa = (Path(d) / "zamocky2022_refs.faa").read_text().splitlines()
            self.assertEqual(faa[0], ">A")
            self.assertEqual([len(x) for x in faa[1:]], [60, 60, 10])
            hdr = (Path(d) / "fetch_manifest.tsv").read_text().splitlines()[0]
            self.assertEqual(hdr.split("\t"), fr.MANIFEST_COLS)


if __name__ == "__main__":
    unittest.main()
