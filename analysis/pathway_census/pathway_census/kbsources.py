"""Cached fetchers for public databases (ChEBI, PubChem, UniProt, PubMed).

Each fetcher takes ``fetch=`` (callable url -> bytes or str) so tests inject stored responses, and
``cache_dir=`` (default ``analysis/pathway_census/.cache/``, git-ignored). The parsed result is
cached as JSON, so a second call for the same key makes no request. The default fetch sends
``USER_AGENT`` and waits so that at most 3 requests per second go out.
"""
from __future__ import annotations

import json
import os
import re
import tempfile
import threading
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

USER_AGENT = "pathway-census/0.1 (BFD-HPCC-master agent)"
CACHE_DIR = Path(__file__).resolve().parents[1] / ".cache"
MAX_RPS = 3
CACHE_VERSION = 1  # bump when the parsed result layout changes

CHEBI_URL = "https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/"
PUBCHEM_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{ik}"
               "/property/InChIKey,MolecularFormula,ExactMass/JSON")
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
PUBMED_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"


class KbFetchError(RuntimeError):
    """A request failed or a response was incomplete. The message names the URL or the id."""


class RateLimiter:
    """Allow at most ``rps`` calls per second by spacing calls 1/rps apart.

    The limit is per process only. Several processes (for example several SLURM jobs on one
    node) each get their own limit, so the combined rate can exceed ``rps``. Threads of one
    process share the limit (a lock serialises ``wait``).
    """

    def __init__(self, rps: float = MAX_RPS, clock=time.monotonic, sleep=time.sleep):
        self.gap, self.clock, self.sleep = 1.0 / rps, clock, sleep
        self.last: float | None = None
        self.lock = threading.Lock()

    def wait(self) -> None:
        with self.lock:
            if self.last is not None:
                delay = self.last + self.gap - self.clock()
                if delay > 0:
                    self.sleep(delay)
            self.last = self.clock()


_LIMITER = RateLimiter()


def default_fetch(url: str) -> bytes:
    _LIMITER.wait()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.read()
    except (urllib.error.URLError, OSError) as e:  # HTTPError is a URLError
        raise KbFetchError(f"request failed for {url}: {e}") from e


def _read_cache(path: Path, url: str):
    """Return the cached result, or None on a miss (absent, unreadable, invalid or stale)."""
    try:
        entry = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if (not isinstance(entry, dict) or entry.get("_v") != CACHE_VERSION
            or entry.get("_url") != url or "data" not in entry):
        return None
    return entry["data"]


def _write_cache(path: Path, url: str, result) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump({"_v": CACHE_VERSION, "_url": url, "data": result}, fh, ensure_ascii=False, indent=1)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def _cached(namespace: str, key: str, url: str, parse, fetch, cache_dir):
    """Fetch, parse, validate and cache. ``parse`` raises on an incomplete response; nothing is
    written to the cache in that case."""
    cache_dir = Path(cache_dir) if cache_dir is not None else CACHE_DIR
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", key)
    path = cache_dir / namespace / f"{safe}.json"
    hit = _read_cache(path, url)
    if hit is not None:
        return hit
    raw = (fetch or default_fetch)(url)
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")
    try:
        result = parse(raw)
    except KbFetchError:
        raise
    except (ValueError, KeyError, IndexError, TypeError, AttributeError, ET.ParseError) as e:
        raise KbFetchError(f"{namespace} {key}: unusable response from {url} ({e!r})") from e
    _write_cache(path, url, result)
    return result


def chebi_compound(chebi_id: str, fetch=None, cache_dir=None) -> dict:
    """name, formula, monoisotopic mass, InChIKey, SMILES, charge from the public ChEBI API."""
    num = re.sub(r"^CHEBI:", "", str(chebi_id).strip(), flags=re.I)

    def parse(raw):
        d = json.loads(raw)
        cd, st = d.get("chemical_data") or {}, d.get("default_structure") or {}
        mono = cd.get("monoisotopic_mass")
        for field, val in (("name", d.get("name")), ("formula", cd.get("formula")),
                           ("inchikey", st.get("standard_inchi_key"))):
            if not val:
                raise KbFetchError(f"ChEBI CHEBI:{num}: response lacks {field}")
        return {"chebi_id": f"CHEBI:{num}", "name": d.get("name"), "formula": cd.get("formula"),
                "monoisotopic_mass": float(mono) if mono not in (None, "") else None,
                "charge": cd.get("charge"), "inchikey": st.get("standard_inchi_key"),
                "smiles": st.get("smiles")}

    return _cached("chebi", num, CHEBI_URL.format(num=num), parse, fetch, cache_dir)


def pubchem_by_inchikey(inchikey: str, fetch=None, cache_dir=None) -> dict:
    """CID, InChIKey, formula, exact mass from PubChem PUG REST."""
    def parse(raw):
        p = json.loads(raw)["PropertyTable"]["Properties"][0]
        if not p.get("MolecularFormula") or not p.get("InChIKey"):
            raise KbFetchError(f"PubChem {inchikey}: response lacks formula or InChIKey")
        return {"cid": p["CID"], "inchikey": p.get("InChIKey"), "formula": p.get("MolecularFormula"),
                "exact_mass": float(p["ExactMass"])}

    return _cached("pubchem", inchikey, PUBCHEM_URL.format(ik=inchikey), parse, fetch, cache_dir)


def uniprot_entry(accession: str, fetch=None, cache_dir=None) -> dict:
    """Gene names, organism, length, sequence, EC list, cited PMIDs, reviewed flag from UniProt."""
    def parse(raw):
        d = json.loads(raw)
        genes: list[str] = []
        for g in d.get("genes", []):
            if "geneName" in g:
                genes.append(g["geneName"]["value"])
            genes += [s["value"] for s in g.get("synonyms", [])]
            genes += [o["value"] for o in g.get("orfNames", [])]
        pd_ = d.get("proteinDescription", {})
        ec: list[str] = []
        names = [pd_.get("recommendedName")] + pd_.get("alternativeNames", []) + pd_.get("submissionNames", [])
        for n in names:
            for e in (n or {}).get("ecNumbers", []):
                if e["value"] not in ec:
                    ec.append(e["value"])
        pmids: list[str] = []
        for ref in d.get("references", []):
            for x in ref.get("citation", {}).get("citationCrossReferences", []) or []:
                if x.get("database") == "PubMed" and x["id"] not in pmids:
                    pmids.append(x["id"])
        if len(d["sequence"]["value"]) != d["sequence"]["length"]:
            raise KbFetchError(f"UniProt {accession}: sequence length does not match length field")
        return {"accession": d["primaryAccession"], "gene_names": genes,
                "organism": d["organism"]["scientificName"], "length": d["sequence"]["length"],
                "sequence": d["sequence"]["value"], "ec": ec, "pmids": pmids,
                "reviewed": d.get("entryType", "").startswith("UniProtKB reviewed")}

    return _cached("uniprot", accession, UNIPROT_URL.format(acc=accession), parse, fetch, cache_dir)


def pubmed_abstract(pmid: str, fetch=None, cache_dir=None) -> dict:
    """Title, abstract text and year from PubMed efetch XML."""
    def parse(raw):
        root = ET.fromstring(raw.encode("utf-8"))
        art = root.find(".//MedlineCitation/Article")
        if art is None:
            raise ValueError(f"no PubMed article for {pmid}")

        def text(el):
            return "".join(el.itertext()).strip() if el is not None else ""

        parts = []
        for at in art.findall("Abstract/AbstractText"):
            label = at.get("Label")
            parts.append((f"{label}: " if label else "") + text(at))
        year = None
        for path in ("Journal/JournalIssue/PubDate/Year", "Journal/JournalIssue/PubDate/MedlineDate",
                     "ArticleDate/Year"):
            m = re.search(r"\d{4}", text(art.find(path)))
            if m:
                year = int(m.group())
                break
        if not text(art.find("ArticleTitle")):
            raise KbFetchError(f"PubMed {pmid}: response lacks a title")
        return {"pmid": str(pmid), "title": text(art.find("ArticleTitle")),
                "abstract": " ".join(parts), "year": year}

    return _cached("pubmed", str(pmid), PUBMED_URL.format(pmid=pmid), parse, fetch, cache_dir)
