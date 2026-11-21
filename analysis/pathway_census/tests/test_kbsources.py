import json
from pathlib import Path

import pytest

from pathway_census import kbsources as ks

KB = Path(__file__).parent / "data" / "kb"


class Fake:
    """Injected fetch: maps a URL substring to a stored fixture file; counts calls."""

    def __init__(self, mapping):
        self.mapping = mapping
        self.urls = []

    def __call__(self, url):
        self.urls.append(url)
        for key, fn in self.mapping.items():
            if key in url:
                return (KB / fn).read_bytes()
        raise AssertionError(f"unexpected url {url}")


def test_chebi_astaxanthin(tmp_path):
    f = Fake({"compound/40968": "chebi_40968.json"})
    r = ks.chebi_compound("CHEBI:40968", fetch=f, cache_dir=tmp_path)
    assert r["inchikey"] == "MQZIGYBFDRPAKN-UWFIBFSHSA-N"
    assert r["name"] == "astaxanthin"
    assert r["formula"] == "C40H52O4"
    assert r["monoisotopic_mass"] == pytest.approx(596.38656)
    assert r["smiles"].startswith("CC1=C(")
    assert r["chebi_id"] == "CHEBI:40968"
    assert "backend/api/public/compound/40968/" in f.urls[0]


def test_chebi_bare_number_same_cache(tmp_path):
    f = Fake({"compound/40968": "chebi_40968.json"})
    ks.chebi_compound("40968", fetch=f, cache_dir=tmp_path)
    ks.chebi_compound("CHEBI:40968", fetch=f, cache_dir=tmp_path)
    assert len(f.urls) == 1


def test_pubchem_by_inchikey(tmp_path):
    f = Fake({"inchikey/MQZIGYBFDRPAKN-UWFIBFSHSA-N": "pubchem_MQZIGYBFDRPAKN.json"})
    r = ks.pubchem_by_inchikey("MQZIGYBFDRPAKN-UWFIBFSHSA-N", fetch=f, cache_dir=tmp_path)
    assert r["cid"] == 5281224
    assert r["formula"] == "C40H52O4"
    assert r["exact_mass"] == pytest.approx(596.38656014)
    assert r["inchikey"] == "MQZIGYBFDRPAKN-UWFIBFSHSA-N"


def test_uniprot_q3hr17(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    r = ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    assert r["length"] == 557
    assert len(r["sequence"]) == 557
    assert "Ast" in r["gene_names"]
    assert r["organism"] == "Phaffia rhodozyma"
    assert r["reviewed"] is False
    assert "25504221" in r["pmids"]
    assert isinstance(r["ec"], list)


def test_uniprot_p21334(tmp_path):
    f = Fake({"P21334.json": "uniprot_P21334.json"})
    r = ks.uniprot_entry("P21334", fetch=f, cache_dir=tmp_path)
    assert r["accession"] == "P21334"
    assert r["length"] == len(r["sequence"])
    assert "Neurospora crassa" in r["organism"]
    assert isinstance(r["reviewed"], bool)


def test_pubmed_abstract(tmp_path):
    f = Fake({"id=16416328": "pubmed_16416328.xml"})
    r = ks.pubmed_abstract("16416328", fetch=f, cache_dir=tmp_path)
    assert r["title"].startswith("Cloning of the astaxanthin synthase gene")
    assert r["abstract"].startswith("A gene has been cloned from Xanthophyllomyces dendrorhous")
    assert r["year"] == 2006
    assert r["pmid"] == "16416328"


def test_cache_hit_makes_no_second_call(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    a = ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    b = ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    assert a == b
    assert len(f.urls) == 1
    assert list(tmp_path.rglob("*.json"))


def test_cache_survives_new_call_with_failing_fetch(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)

    def boom(url):
        raise AssertionError("network used")

    assert ks.uniprot_entry("Q3HR17", fetch=boom, cache_dir=tmp_path)["length"] == 557


def test_user_agent_has_no_personal_info():
    assert ks.USER_AGENT == "pathway-census/0.1 (BFD-HPCC-master agent)"


def test_rate_limiter_max_three_per_second():
    now = [0.0]
    slept = []

    def clock():
        return now[0]

    def sleep(s):
        slept.append(s)
        now[0] += s

    rl = ks.RateLimiter(3, clock=clock, sleep=sleep)
    for _ in range(4):
        rl.wait()
    assert now[0] == pytest.approx(3 / 3)  # 4 requests need 3 gaps of 1/3 s
    assert all(s <= 1 / 3 + 1e-9 for s in slept)


def test_default_fetch_sets_user_agent(monkeypatch):
    seen = {}

    class Resp:
        def read(self):
            return b"{}"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        seen["ua"] = req.get_header("User-agent")
        return Resp()

    monkeypatch.setattr(ks.urllib.request, "urlopen", fake_urlopen)
    ks.RateLimiter(1000)  # constructing must not hit the network
    assert ks.default_fetch("https://example.invalid/x") == b"{}"
    assert seen["ua"] == ks.USER_AGENT
