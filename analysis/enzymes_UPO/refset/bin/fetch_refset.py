#!/usr/bin/env python3
"""Fetch the Zamocky et al. 2022 reference sequences listed in tableS1_resolved.tsv.

Each row of the resolved table names a `route` that says where its sequence comes from:

  uniprot_direct, uniprot_xref, inferred_unique_PF01328_len_match
        UniProt REST, by accession. A `resolved_id` holding two accessions separated by ';'
        (one PeroxiBase ID mapped to two UniProt entries) fetches both and keeps the first;
        the manifest says whether the two sequences were identical.
  ncbi_efetch
        NCBI E-utilities efetch, protein FASTA.
  wayback_view_perox
        PeroxiBase is offline. The Wayback Machine keeps its entry pages
        (peroxibase.toulouse.inra.fr/display_perox/view_perox/<id>); the protein is in the
        page's hidden form field name="Sequence". The PeroxiBase class and orthogroup on the
        same page are recorded too.
  pdb ...
        RCSB FASTA for the entry, the chain named in `resolved_id` (e.g. 5FUJ_A).

Outputs, in --out:
  zamocky2022_refs.faa      one record per fetched row, ID = Table S1 abbreviation
  fetch_manifest.tsv        one row per Table S1 row: status, source URL, fetched vs. Table S1
                            length, sha1 of the sequence, PeroxiBase class/orthogroup, notes

Nothing is silently dropped: a row that cannot be fetched gets status=missing in the manifest.
Exit status is 0 unless --strict is given and a row is missing or a length disagrees.

Requires only the Python 3.9+ standard library (on the UCR HPCC use /usr/bin/python3.12).
"""
from __future__ import annotations

import argparse
import csv
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

USER_AGENT = "delta-gain-fungi-refset/0.1 (research; python-urllib)"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/accessions?{q}"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{q}"
WAYBACK_AVAIL = "https://archive.org/wayback/available?{q}"
WAYBACK_RAW = "https://web.archive.org/web/{ts}id_/{url}"
PEROXIBASE_PAGE = "peroxibase.toulouse.inra.fr/display_perox/view_perox/{id}"
RCSB_FASTA = "https://www.rcsb.org/fasta/entry/{entry}"

AA = re.compile(r"[^A-Z]")
MANIFEST_COLS = ["abbrev", "given_id", "route", "resolved_id", "status", "source_url",
                 "fetched_len", "table_len", "len_match", "sha1", "pb_entry", "pb_class",
                 "pb_orthogroup", "note"]


# ----------------------------------------------------------------------------- network

class Fetcher:
    """HTTP GET with retries, a polite delay per host, and an optional on-disk cache."""

    def __init__(self, cache: Path | None = None, delay: float = 1.0, retries: int = 4,
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
        err: Exception | None = None
        for attempt in range(self.retries):
            wait = self.delay - (time.monotonic() - self._last.get(host, 0.0))
            if wait > 0:
                time.sleep(wait)
            self._last[host] = time.monotonic()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    text = r.read().decode("utf-8", "replace")
                if key:
                    key.write_text(text)
                return text
            except urllib.error.HTTPError as e:
                if e.code in (400, 404):
                    raise
                err = e
            except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
                err = e
            time.sleep(2 ** attempt * 2)
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
    """Extract protein sequence, entry name, class and orthogroup from a PeroxiBase entry page."""
    m = re.search(r'name="Sequence"\s+value="([^"]*)"', page)
    seq = AA.sub("", html.unescape(m.group(1)).upper()) if m else ""
    entry = re.search(r"Entry information\s*:\s*([^\s<]+)", page)
    cls = re.search(r"/search/search_type/Class/([^\"]+)\"", page)
    og = re.search(r"view_orthogroup/([A-Za-z0-9_]+)", page)
    return {
        "seq": seq,
        "pb_entry": entry.group(1) if entry else "",
        "pb_class": urllib.parse.unquote(cls.group(1)) if cls else "",
        "pb_orthogroup": og.group(1) if og else "",
    }


def pick_pdb_chain(fasta_text: str, chain: str) -> str:
    """RCSB entry FASTA headers look like '5FUJ_1|Chains A, B|...'; return the named chain."""
    for header, seq in parse_fasta(fasta_text):
        f = header.split("|")
        if len(f) > 1:
            chains = re.findall(r"\b([A-Za-z0-9]+)(?:\[auth [^\]]+\])?", f[1].replace("Chains", "").replace("Chain", ""))
            if chain in chains:
                return seq
    return ""


# ----------------------------------------------------------------------------- routes

def fetch_uniprot(fetcher: Fetcher, accessions: list[str]) -> dict[str, str]:
    seqs: dict[str, str] = {}
    for i in range(0, len(accessions), 200):
        chunk = accessions[i:i + 200]
        q = urllib.parse.urlencode({"accessions": ",".join(chunk), "format": "fasta"})
        for h, s in parse_fasta(fetcher.get(UNIPROT_FASTA.format(q=q))):
            seqs[uniprot_accession(h)] = s
    return seqs


def fetch_ncbi(fetcher: Fetcher, accessions: list[str]) -> dict[str, str]:
    q = urllib.parse.urlencode({"db": "protein", "id": ",".join(accessions),
                                "rettype": "fasta", "retmode": "text", "tool": "fetch_refset"})
    seqs: dict[str, str] = {}
    for h, s in parse_fasta(fetcher.get(NCBI_EFETCH.format(q=q))):
        acc = h.split()[0]
        seqs[acc] = s
        seqs[acc.split(".")[0]] = s
    return seqs


def fetch_wayback(fetcher: Fetcher, pb_id: str) -> tuple[dict, str]:
    page_url = PEROXIBASE_PAGE.format(id=pb_id)
    avail = json.loads(fetcher.get(WAYBACK_AVAIL.format(q=urllib.parse.urlencode({"url": page_url}))))
    snap = avail.get("archived_snapshots", {}).get("closest") or {}
    if not snap.get("available") or str(snap.get("status")) != "200":
        return {"seq": ""}, ""
    # snap["url"] is http://web.archive.org/web/<ts>/<orig>; rebuild it as the raw (id_) capture,
    # which serves the page as archived, without the Wayback toolbar rewriting.
    orig = re.sub(r"^https?://web\.archive\.org/web/\d+/", "", snap["url"])
    raw = WAYBACK_RAW.format(ts=snap["timestamp"], url=orig)
    return parse_peroxibase_page(fetcher.get(raw)), raw


# ----------------------------------------------------------------------------- main logic

def read_table(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(rows: list[dict], fetcher: Fetcher, log=sys.stderr) -> tuple[list[tuple[str, str]], list[dict]]:
    """Fetch every row. Returns (fasta_records, manifest_rows)."""
    up_rows = [r for r in rows if r["route"] in
               ("uniprot_direct", "uniprot_xref", "inferred_unique_PF01328_len_match")]
    up_accs = sorted({a for r in up_rows for a in r["resolved_id"].split(";") if a})
    ncbi_rows = [r for r in rows if r["route"] == "ncbi_efetch"]

    up = fetch_uniprot(fetcher, up_accs) if up_accs else {}
    print(f"uniprot: {len(up)} of {len(up_accs)} accessions returned", file=log)
    nc = fetch_ncbi(fetcher, [r["resolved_id"] for r in ncbi_rows]) if ncbi_rows else {}
    print(f"ncbi: {len(ncbi_rows)} requested", file=log)

    records, manifest = [], []
    for r in rows:
        m = {c: "" for c in MANIFEST_COLS}
        m.update({k: r.get(k, "") for k in ("abbrev", "given_id", "route", "resolved_id")})
        m["table_len"] = r.get("length_aa", "")
        notes, seq, url = [], "", ""
        route = r["route"]
        try:
            if route in ("uniprot_direct", "uniprot_xref", "inferred_unique_PF01328_len_match"):
                accs = [a for a in r["resolved_id"].split(";") if a]
                seq = up.get(accs[0], "")
                url = "https://rest.uniprot.org/uniprotkb/" + accs[0]
                if len(accs) > 1:
                    same = len({up.get(a, "") for a in accs}) == 1
                    notes.append(f"{len(accs)} UniProt entries for this PeroxiBase ID; kept {accs[0]}; "
                                 f"sequences {'identical' if same else 'DIFFER'}")
                if route.startswith("inferred"):
                    notes.append("ID inferred (unique PF01328 entry at this length in this organism)")
            elif route == "ncbi_efetch":
                acc = r["resolved_id"]
                seq = nc.get(acc, nc.get(acc.split(".")[0], ""))
                url = "https://www.ncbi.nlm.nih.gov/protein/" + acc
            elif route == "wayback_view_perox":
                pb_id = r["resolved_id"].removeprefix("peroxibase:")
                info, url = fetch_wayback(fetcher, pb_id)
                seq = info["seq"]
                for k in ("pb_entry", "pb_class", "pb_orthogroup"):
                    m[k] = info.get(k, "")
                if not url:
                    notes.append("no Wayback capture")
            elif route.startswith("pdb"):
                entry, _, chain = r["resolved_id"].partition("_")
                url = RCSB_FASTA.format(entry=entry)
                seq = pick_pdb_chain(fetcher.get(url), chain or "A")
                notes.append("structure chain; duplicate of a backbone enzyme (see route)")
            else:
                notes.append(f"unknown route {route!r}")
        except Exception as e:  # recorded, never fatal: one bad row must not lose the rest
            notes.append(f"error: {e}")
        m["source_url"] = url
        if seq:
            m["status"] = "ok"
            m["fetched_len"] = str(len(seq))
            m["sha1"] = hashlib.sha1(seq.encode()).hexdigest()
            try:
                m["len_match"] = "yes" if int(m["table_len"]) == len(seq) else "no"
            except ValueError:
                m["len_match"] = "na"
            records.append((r["abbrev"], seq))
        else:
            m["status"] = "missing"
        m["note"] = "; ".join(notes)
        manifest.append(m)
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
    ap.add_argument("--cache", type=Path, help="cache raw HTTP responses here (reruns skip the network)")
    ap.add_argument("--delay", type=float, default=1.0, help="minimum seconds between requests to one host")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if any row is missing or disagrees with the Table S1 length")
    a = ap.parse_args(argv)

    rows = read_table(a.table)
    records, manifest = run(rows, Fetcher(cache=a.cache, delay=a.delay))
    write_outputs(a.out, records, manifest)

    n_ok = sum(m["status"] == "ok" for m in manifest)
    bad_len = [m["abbrev"] for m in manifest if m["len_match"] == "no"]
    missing = [m["abbrev"] for m in manifest if m["status"] == "missing"]
    print(f"fetched {n_ok}/{len(manifest)}; missing {len(missing)} {missing}; "
          f"length differs from Table S1: {len(bad_len)} {bad_len}", file=sys.stderr)
    return 1 if a.strict and (missing or bad_len) else 0


if __name__ == "__main__":
    sys.exit(main())
