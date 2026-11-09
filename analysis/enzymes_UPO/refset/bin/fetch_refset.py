#!/usr/bin/env python3
"""Fetch the Zamocky et al. 2022 reference sequences listed in tableS1_resolved.tsv.

Sources, by the row's `route`:

  uniprot_direct, inferred_unique_PF01328_len_match
        UniProtKB by accession. Entries since deleted from UniProtKB (many were, as "not part of
        a reference proteome") are taken from UniParc, which keeps every sequence ever issued.
  uniprot_xref
        A PeroxiBase ID that UniProt cross-references. UniProtKB first; if its length disagrees
        with Table S1, the archived PeroxiBase sequence is tried, since PeroxiBase was Zamocky's
        primary source and the UniProt entry can be a different gene model.
  wayback_view_perox
        PeroxiBase ID with no UniProt link. PeroxiBase is offline, so the sequence comes from the
        Wayback Machine: first the archived FASTA endpoint (tools/get_fasta/<id>/PEP), then the
        archived entry page (display_perox/view_perox/<id>, hidden field name="Sequence"), on both
        the old (peroxibase.toulouse.inra.fr) and new (redoxibase.toulouse.inrae.fr) domains.
        The capture closest to 2022-04-07, Zamocky's stated access date, is used.
  ncbi_efetch   NCBI protein efetch.
  pdb ...       RCSB FASTA, the chain named in resolved_id (e.g. 5FUJ_A).

Table S1's length is the arbiter: when a row has more than one possible source, the first whose
length matches Table S1 wins; if none matches, the first non-empty one is kept and flagged.

Outputs, in --out:
  zamocky2022_refs.faa   one record per fetched row, ID = Table S1 abbreviation
  fetch_manifest.tsv     one row per Table S1 row (nothing is dropped): status, source, URL,
                         fetched vs. Table S1 length, sha1, PeroxiBase entry/class/orthogroup, notes

Exit status is 0 unless --strict is given and a row is missing or disagrees in length.
Python 3.9+ standard library only (on the UCR HPCC use /usr/bin/python3.12).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

USER_AGENT = "delta-gain-fungi-refset/0.2 (research; python-urllib)"
UNIPROT_BATCH = "https://rest.uniprot.org/uniprotkb/accessions?{q}"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPARC_FASTA = "https://rest.uniprot.org/uniparc/{upi}.fasta"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{q}"
WAYBACK_CDX = "https://web.archive.org/cdx/search/cdx?{q}"
WAYBACK_RAW = "https://web.archive.org/web/{ts}id_/{url}"
RCSB_FASTA = "https://www.rcsb.org/fasta/entry/{entry}"

PB_HOSTS = ("peroxibase.toulouse.inra.fr", "redoxibase.toulouse.inrae.fr")
PB_TARGET = dt.datetime(2022, 4, 7)        # Zamocky 2022, Methods 2.1: "accessed on 7 April 2022"
PB_TRIES_PER_PATH = 3                      # captures tried per URL before moving on
SLOW_HOSTS = {"web.archive.org": 3.0}      # seconds between requests; CDX throttles hard

AA = re.compile(r"[^A-Z]")
MANIFEST_COLS = ["abbrev", "given_id", "route", "resolved_id", "status", "source", "source_url",
                 "fetched_len", "table_len", "len_match", "sha1", "pb_entry", "pb_class",
                 "pb_orthogroup", "note"]


# ----------------------------------------------------------------------------- network

class NotFound(Exception):
    pass


class Fetcher:
    """HTTP GET with a per-host delay, retries with backoff (429/5xx/network), and a cache.

    404 raises NotFound immediately. Only successful responses are cached, so a rerun
    retries exactly what failed.
    """

    def __init__(self, cache: Path | None = None, delay: float = 1.0, retries: int = 6,
                 timeout: float = 60.0):
        self.cache, self.delay, self.retries, self.timeout = cache, delay, retries, timeout
        self._last: dict[str, float] = {}
        if cache:
            cache.mkdir(parents=True, exist_ok=True)

    def get(self, url: str) -> str:
        key = None
        if self.cache:
            key = self.cache / (hashlib.sha1(url.encode()).hexdigest() + ".txt")
            if key.exists():
                return key.read_text()
        host = urllib.parse.urlsplit(url).netloc
        delay = SLOW_HOSTS.get(host, self.delay)
        err: Exception | None = None
        for attempt in range(self.retries):
            wait = delay - (time.monotonic() - self._last.get(host, 0.0))
            if wait > 0:
                time.sleep(wait)
            self._last[host] = time.monotonic()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    text = r.read().decode("utf-8", "replace")
                if key and text.strip():         # an empty body may be throttling: never cache it
                    key.write_text(text)
                return text
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    raise NotFound(url) from e
                if e.code not in (429, 500, 502, 503, 504):
                    raise
                err = e
            except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
                err = e
            time.sleep(min(120, 5 * 2 ** attempt))
        raise RuntimeError(f"GET failed after {self.retries} tries: {url}: {err}")


# ----------------------------------------------------------------------------- parsers

def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Return [(header_without_>, sequence)] with sequence upper-cased, letters only."""
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, AA.sub("", "".join(buf).upper())))
            name, buf = line[1:].strip(), []
        elif name is not None:
            buf.append(line.strip())
    if name is not None:
        out.append((name, AA.sub("", "".join(buf).upper())))
    return out


def uniprot_accession(header: str) -> str:
    """'sp|P04963|PRXC_LEPFU ...' or 'tr|A0A...|...' -> accession."""
    f = header.split("|")
    return f[1] if len(f) >= 3 else header.split()[0]


def parse_peroxibase_page(page: str) -> dict:
    """PeroxiBase entry page: protein from the hidden Sequence field, plus entry/class/orthogroup.

    Some entry pages do not inline the sequence (only a FASTA link); seq is then ''.
    """
    m = re.search(r'name="Sequence"\s+value="([^"]*)"', page)
    entry = re.search(r"Entry information\s*:\s*([^\s<]+)", page)
    cls = re.search(r"/search/search_type/Class/([^\"]+)\"", page)
    og = re.search(r"view_orthogroup/([A-Za-z0-9_]+)", page)
    return {
        "seq": AA.sub("", html.unescape(m.group(1)).upper()) if m else "",
        "pb_entry": entry.group(1) if entry else "",
        "pb_class": urllib.parse.unquote(cls.group(1)) if cls else "",
        "pb_orthogroup": og.group(1) if og else "",
    }


PB_FASTA = re.compile(r">(\d+)\|([^|\n]*)\|[^|\n]*\|[^|\n]*\|([^|\n]*)\|[^\n]*\n([A-Z*\s]+)")


def parse_peroxibase_fasta(page: str) -> dict:
    """PeroxiBase tools/get_fasta/<id>/PEP page: the FASTA sits in a textarea.

    Header: >12616|SbreHalPrx01| [Suillus_brevipes]|Basidiomycota|Haloperoxidase_(haem)|complete|
    """
    m = PB_FASTA.search(html.unescape(page))
    if not m:
        return {"seq": ""}
    return {"seq": AA.sub("", m.group(4)), "pb_entry": m.group(2).strip(),
            "pb_class": m.group(3).replace("_", " ").strip(), "pb_orthogroup": ""}


def pick_pdb_chain(fasta_text: str, chain: str) -> str:
    """RCSB entry FASTA headers look like '5FUJ_1|Chains A, B|...'; return the named chain."""
    for header, seq in parse_fasta(fasta_text):
        f = header.split("|")
        if len(f) > 1:
            names = re.findall(r"\b([A-Za-z0-9]+)(?:\[auth [^\]]+\])?",
                               f[1].replace("Chains", "").replace("Chain", ""))
            if chain in names:
                return seq
    return ""


def closest_first(captures: list[tuple[str, str]], target: dt.datetime = PB_TARGET):
    """Sort Wayback (timestamp, original) pairs by distance from target."""
    def dist(c):
        try:
            return abs((dt.datetime.strptime(c[0][:14], "%Y%m%d%H%M%S") - target).total_seconds())
        except ValueError:
            return float("inf")
    return sorted(captures, key=dist)


# ----------------------------------------------------------------------------- sources

class Sources:
    """Each method returns (seq, info) where info has 'source', 'url' and optional extras."""

    def __init__(self, fetcher: Fetcher, log=sys.stderr):
        self.f, self.log = fetcher, log
        self.up: dict[str, str] = {}
        self.nc: dict[str, str] = {}

    # --- UniProt
    def preload_uniprot(self, accs: list[str]) -> None:
        for i in range(0, len(accs), 200):
            q = urllib.parse.urlencode({"accessions": ",".join(accs[i:i + 200]), "format": "fasta"})
            for h, s in parse_fasta(self.f.get(UNIPROT_BATCH.format(q=q))):
                self.up[uniprot_accession(h)] = s
        print(f"uniprotkb: {len(self.up)} of {len(accs)} accessions active", file=self.log)

    def uniprot(self, acc: str):
        if acc in self.up:
            return self.up[acc], {"source": "uniprotkb", "url": f"https://rest.uniprot.org/uniprotkb/{acc}"}
        entry = json.loads(self.f.get(UNIPROT_ENTRY.format(acc=acc)))
        upi = entry.get("extraAttributes", {}).get("uniParcId", "")
        why = entry.get("inactiveReason", {})
        reason = why.get("inactiveReasonType", "?")
        if why.get("deletedReason"):
            reason += f": {why['deletedReason']}"
        if why.get("mergeDemergeTo"):
            reason += f" -> {','.join(why['mergeDemergeTo'])}"
        if not upi:
            return "", {"source": "uniprotkb", "note": f"{acc} inactive ({reason}), no UniParc ID"}
        recs = parse_fasta(self.f.get(UNIPARC_FASTA.format(upi=upi)))
        seq = recs[0][1] if recs else ""
        return seq, {"source": "uniparc", "url": f"https://rest.uniprot.org/uniparc/{upi}",
                     "note": f"{acc} inactive in UniProtKB ({reason}); sequence from UniParc {upi}"}

    # --- NCBI
    def preload_ncbi(self, accs: list[str]) -> None:
        if not accs:
            return
        q = urllib.parse.urlencode({"db": "protein", "id": ",".join(accs), "rettype": "fasta",
                                    "retmode": "text", "tool": "fetch_refset"})
        for h, s in parse_fasta(self.f.get(NCBI_EFETCH.format(q=q))):
            acc = h.split()[0]
            self.nc[acc] = self.nc[acc.split(".")[0]] = s
        print(f"ncbi: {len(accs)} requested", file=self.log)

    def ncbi(self, acc: str):
        seq = self.nc.get(acc, self.nc.get(acc.split(".")[0], ""))
        return seq, {"source": "ncbi", "url": f"https://www.ncbi.nlm.nih.gov/protein/{acc}"}

    # --- PeroxiBase via the Wayback Machine
    def _captures(self, host_path: str) -> list[tuple[str, str]]:
        q = urllib.parse.urlencode({"url": host_path, "output": "json",
                                    "fl": "timestamp,original", "filter": "statuscode:200"})
        text = self.f.get(WAYBACK_CDX.format(q=q)).strip()
        rows = json.loads(text) if text else []
        return [(r[0], r[1]) for r in rows[1:]]          # rows[0] is the header

    def peroxibase(self, pb_id: str):
        notes = []
        for path, parse in ((f"tools/get_fasta/{pb_id}/PEP", parse_peroxibase_fasta),
                            (f"display_perox/view_perox/{pb_id}", parse_peroxibase_page)):
            caps = []
            for host in PB_HOSTS:
                caps += self._captures(f"{host}/{path}")
            for ts, orig in closest_first(caps)[:PB_TRIES_PER_PATH]:
                url = WAYBACK_RAW.format(ts=ts, url=orig)
                info = parse(self.f.get(url))
                if info["seq"]:
                    info.update(source="peroxibase_wayback", url=url)
                    return info.pop("seq"), info
            notes.append(f"{path}: {len(caps)} captures, no sequence")
        return "", {"source": "peroxibase_wayback", "note": "; ".join(notes)}

    # --- RCSB
    def pdb(self, entry_chain: str):
        entry, _, chain = entry_chain.partition("_")
        url = RCSB_FASTA.format(entry=entry)
        return pick_pdb_chain(self.f.get(url), chain or "A"), {
            "source": "rcsb", "url": url, "note": "structure chain; duplicate of a backbone enzyme"}


# ----------------------------------------------------------------------------- main logic

UNIPROT_ROUTES = ("uniprot_direct", "uniprot_xref", "inferred_unique_PF01328_len_match")


def candidates(r: dict, src: Sources):
    """Ordered zero-arg callables, each returning (seq, info)."""
    route, rid = r["route"], r["resolved_id"]
    accs = [a for a in rid.split(";") if a]
    if route in ("uniprot_direct", "inferred_unique_PF01328_len_match"):
        return [lambda: src.uniprot(accs[0])]
    if route == "uniprot_xref":
        return [lambda: src.uniprot(accs[0]), lambda: src.peroxibase(r["given_id"])]
    if route == "wayback_view_perox":
        return [lambda: src.peroxibase(rid.removeprefix("peroxibase:"))]
    if route == "ncbi_efetch":
        return [lambda: src.ncbi(rid)]
    if route.startswith("pdb"):
        return [lambda: src.pdb(rid)]
    return []


def resolve_row(r: dict, src: Sources) -> dict:
    m = {c: "" for c in MANIFEST_COLS}
    m.update({k: r.get(k, "") for k in ("abbrev", "given_id", "route", "resolved_id")})
    m["table_len"] = r.get("length_aa", "")
    want = int(m["table_len"]) if m["table_len"].isdigit() else None
    notes, tried, chosen = [], [], None
    cands = candidates(r, src)
    if not cands:
        notes.append(f"unknown route {r['route']!r}")
    for get in cands:
        try:
            seq, info = get()
        except Exception as e:                    # recorded, never fatal
            seq, info = "", {"note": f"error: {e}"}
        if info.get("note"):
            notes.append(info["note"])
        if seq:
            tried.append((seq, info))
            if want is None or len(seq) == want:
                chosen = (seq, info)
                break
    if chosen is None and tried:
        chosen = tried[0]
    if len(tried) > 1:
        notes.append("lengths by source: " + ", ".join(f"{i['source']}={len(s)}" for s, i in tried)
                     + f"; kept {chosen[1]['source']}")

    accs = [a for a in r["resolved_id"].split(";") if a]
    if r["route"] in UNIPROT_ROUTES and len(accs) > 1:
        same = len({src.up.get(a, "") for a in accs}) == 1
        notes.append(f"{len(accs)} UniProt entries for this PeroxiBase ID; first used; "
                     f"sequences {'identical' if same else 'DIFFER'}")
    if r["route"].startswith("inferred"):
        notes.append("ID inferred (unique PF01328 entry at this length in this organism)")

    if chosen:
        seq, info = chosen
        m.update(status="ok", source=info.get("source", ""), source_url=info.get("url", ""),
                 fetched_len=str(len(seq)), sha1=hashlib.sha1(seq.encode()).hexdigest(),
                 len_match="na" if want is None else ("yes" if len(seq) == want else "no"))
        for k in ("pb_entry", "pb_class", "pb_orthogroup"):
            m[k] = info.get(k, "")
        m["_seq"] = seq
    else:
        m["status"] = "missing"
    m["note"] = "; ".join(dict.fromkeys(notes))
    return m


def read_table(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(rows: list[dict], fetcher: Fetcher, log=sys.stderr):
    """Fetch every row. Returns (fasta_records, manifest_rows)."""
    src = Sources(fetcher, log)
    src.preload_uniprot(sorted({a for r in rows if r["route"] in UNIPROT_ROUTES
                                for a in r["resolved_id"].split(";") if a}))
    src.preload_ncbi([r["resolved_id"] for r in rows if r["route"] == "ncbi_efetch"])
    records, manifest = [], []
    for i, r in enumerate(rows, 1):
        m = resolve_row(r, src)
        seq = m.pop("_seq", "")
        if seq:
            records.append((r["abbrev"], seq))
        manifest.append(m)
        print(f"[{i}/{len(rows)}] {m['abbrev']} {m['status']} {m['source']} "
              f"{m['fetched_len']}/{m['table_len']}", file=log, flush=True)
    return records, manifest


def write_outputs(out: Path, records, manifest) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "zamocky2022_refs.faa").open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    with (out / "fetch_manifest.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(manifest)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--table", required=True, type=Path, help="tableS1_resolved.tsv")
    ap.add_argument("--out", required=True, type=Path, help="output directory")
    ap.add_argument("--cache", type=Path, help="cache successful HTTP responses here")
    ap.add_argument("--delay", type=float, default=1.0,
                    help="seconds between requests to one host (web.archive.org always >= 3)")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if any row is missing or disagrees with the Table S1 length")
    a = ap.parse_args(argv)

    records, manifest = run(read_table(a.table), Fetcher(cache=a.cache, delay=a.delay))
    write_outputs(a.out, records, manifest)

    missing = [m["abbrev"] for m in manifest if m["status"] == "missing"]
    bad_len = [m["abbrev"] for m in manifest if m["len_match"] == "no"]
    by_src: dict[str, int] = {}
    for m in manifest:
        if m["status"] == "ok":
            by_src[m["source"]] = by_src.get(m["source"], 0) + 1
    print(f"fetched {len(records)}/{len(manifest)} by source {by_src}; missing {len(missing)} "
          f"{missing}; length differs from Table S1: {len(bad_len)} {bad_len}", file=sys.stderr)
    return 1 if a.strict and (missing or bad_len) else 0


if __name__ == "__main__":
    sys.exit(main())
