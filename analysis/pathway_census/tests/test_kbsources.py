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


# ---- fix round 1 -------------------------------------------------------------------------
import os  # noqa: E402
import urllib.error  # noqa: E402


def _bytes(s):
    return s.encode()


def test_truncated_cache_is_a_miss(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    p = tmp_path / "uniprot" / "Q3HR17.json"
    p.parent.mkdir(parents=True)
    p.write_text('{"_v": 1, "da')
    r = ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    assert r["length"] == 557
    assert len(f.urls) == 1
    assert json.loads(p.read_text())["data"]["length"] == 557  # rewritten whole


def test_wrong_shape_cache_is_a_miss(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    p = tmp_path / "uniprot" / "Q3HR17.json"
    p.parent.mkdir(parents=True)
    p.write_text("[1, 2]")
    assert ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)["length"] == 557


def test_cache_version_mismatch_is_a_miss(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    p = tmp_path / "uniprot" / "Q3HR17.json"
    e = json.loads(p.read_text())
    e["_v"] = -1
    p.write_text(json.dumps(e))
    ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    assert len(f.urls) == 2


def test_cache_url_mismatch_is_a_miss(tmp_path):
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    p = tmp_path / "uniprot" / "Q3HR17.json"
    e = json.loads(p.read_text())
    e["_url"] = "https://old.example/x"
    p.write_text(json.dumps(e))
    ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    assert len(f.urls) == 2


def test_cache_write_is_atomic(tmp_path, monkeypatch):
    calls = []
    real = os.replace
    monkeypatch.setattr(ks.os, "replace", lambda a, b: (calls.append((a, b)), real(a, b))[1])
    f = Fake({"Q3HR17.json": "uniprot_Q3HR17.json"})
    ks.uniprot_entry("Q3HR17", fetch=f, cache_dir=tmp_path)
    assert len(calls) == 1
    src, dst = calls[0]
    assert Path(src).parent == Path(dst).parent and Path(src) != Path(dst)
    assert [p.name for p in (tmp_path / "uniprot").iterdir()] == ["Q3HR17.json"]  # no temp left


def _empty_cache(tmp_path):
    return not list(tmp_path.rglob("*.json"))


def test_chebi_missing_inchikey_raises_and_not_cached(tmp_path):
    body = _bytes('{"id": 1, "name": "x", "chemical_data": {"formula": "C"}, "default_structure": {}}')
    with pytest.raises(ks.KbFetchError, match=r"CHEBI:40968.*inchikey"):
        ks.chebi_compound("40968", fetch=lambda u: body, cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_chebi_missing_formula_raises(tmp_path):
    body = _bytes('{"id": 1, "name": "x", "chemical_data": {}, "default_structure": {"standard_inchi_key": "K"}}')
    with pytest.raises(ks.KbFetchError, match=r"CHEBI:40968.*formula"):
        ks.chebi_compound("40968", fetch=lambda u: body, cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_chebi_error_json_raises(tmp_path):
    with pytest.raises(ks.KbFetchError, match="CHEBI:1"):
        ks.chebi_compound("1", fetch=lambda u: b'{"detail": "Not found"}', cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_pubchem_fault_not_cached(tmp_path):
    body = b'{"Fault": {"Code": "PUGREST.NotFound"}}'
    with pytest.raises(ks.KbFetchError, match="AAAA"):
        ks.pubchem_by_inchikey("AAAA", fetch=lambda u: body, cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_uniprot_error_not_cached(tmp_path):
    with pytest.raises(ks.KbFetchError, match="Q00000"):
        ks.uniprot_entry("Q00000", fetch=lambda u: b'{"messages": ["not found"]}', cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_uniprot_length_mismatch_not_cached(tmp_path):
    d = json.loads((KB / "uniprot_Q3HR17.json").read_text())
    d["sequence"]["length"] = 5
    with pytest.raises(ks.KbFetchError, match="Q3HR17"):
        ks.uniprot_entry("Q3HR17", fetch=lambda u: json.dumps(d), cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_pubmed_empty_set_not_cached(tmp_path):
    body = b"<?xml version='1.0'?><PubmedArticleSet></PubmedArticleSet>"
    with pytest.raises(ks.KbFetchError, match="99"):
        ks.pubmed_abstract("99", fetch=lambda u: body, cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_pubmed_garbage_xml_not_cached(tmp_path):
    with pytest.raises(ks.KbFetchError, match="99"):
        ks.pubmed_abstract("99", fetch=lambda u: b"<html>", cache_dir=tmp_path)
    assert _empty_cache(tmp_path)


def test_default_fetch_wraps_http_error(monkeypatch):
    def boom(req, timeout=None):
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(ks.urllib.request, "urlopen", boom)
    with pytest.raises(ks.KbFetchError, match="https://example.invalid/y") as ei:
        ks.default_fetch("https://example.invalid/y")
    assert isinstance(ei.value.__cause__, urllib.error.HTTPError)


def test_default_fetch_wraps_url_error(monkeypatch):
    def boom(req, timeout=None):
        raise urllib.error.URLError("refused")

    monkeypatch.setattr(ks.urllib.request, "urlopen", boom)
    with pytest.raises(ks.KbFetchError, match="https://example.invalid/z") as ei:
        ks.default_fetch("https://example.invalid/z")
    assert isinstance(ei.value.__cause__, urllib.error.URLError)


def test_rate_limiter_has_lock():
    import threading
    assert isinstance(ks.RateLimiter(3).lock, type(threading.Lock()))
    assert "per process" in ks.RateLimiter.__doc__
