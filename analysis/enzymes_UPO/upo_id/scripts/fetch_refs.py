#!/usr/bin/env python3
"""Fetch reference protein sequences listed in a sources TSV.

Sources: UniProt REST, NCBI E-utilities (protein), Ensembl REST, RCSB PDB FASTA.
Writes one FASTA (record ID = name column) and a manifest TSV with
accession, source URL, original header, length, sha256 of the sequence and fetch date.

Polite use: one request at a time, a pause between requests, exponential
back-off on HTTP 429/5xx.

Usage:
  fetch_refs.py --sources ref/ref_sources.tsv --out-fasta ref/refs.fa --out-manifest ref/manifest.tsv
"""
import argparse
import csv
import datetime
import hashlib
import sys
import time
import urllib.error
import urllib.request

UA = "upo_id-fetch/0.1 (UCR Stajich lab; contact jason.stajich@ucr.edu)"


def get(url, tries=6, pause=1.0):
    delay = 5.0
    for i in range(tries):
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read().decode()
            time.sleep(pause)
            return data
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                ra = e.headers.get("Retry-After")
                wait = float(ra) if ra and ra.isdigit() else delay
                sys.stderr.write(f"HTTP {e.code} for {url}; wait {wait}s\n")
                time.sleep(wait)
                delay *= 2
                continue
            raise
        except urllib.error.URLError as e:
            sys.stderr.write(f"URL error {e} for {url}; wait {delay}s\n")
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(f"failed after {tries} tries: {url}")


def parse_fasta(text):
    recs, h, s = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if h is not None:
                recs.append((h, "".join(s)))
            h, s = line[1:].strip(), []
        elif line.strip():
            s.append(line.strip())
    if h is not None:
        recs.append((h, "".join(s)))
    return recs


def fetch_one(db, acc):
    if db == "uniprot":
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        recs = parse_fasta(get(url))
    elif db == "ensembl":
        url = f"https://rest.ensembl.org/sequence/id/{acc}?type=protein;content-type=text/x-fasta"
        recs = parse_fasta(get(url))
    elif db == "pdb":
        url = f"https://www.rcsb.org/fasta/entry/{acc}"
        recs = parse_fasta(get(url))[:1]  # first entity = enzyme chain
    else:
        raise ValueError(db)
    if not recs:
        raise RuntimeError(f"no sequence for {db}:{acc}")
    return url, recs[0]


def fetch_ncbi(accs):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta"
           "&retmode=text&tool=upo_id&email=jason.stajich@ucr.edu&id=" + ",".join(accs))
    recs = parse_fasta(get(url, pause=0.5))
    out = {}
    for h, s in recs:
        a = h.split()[0]
        out[a] = (h, s)
        out[a.split(".")[0]] = (h, s)
    return url.split("&id=")[0] + "&id=<acc>", out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-manifest", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.sources), delimiter="\t"))
    today = datetime.date.today().isoformat()
    ncbi_accs = [r["accession"] for r in rows if r["db"] == "ncbi"]
    ncbi = {}
    ncbi_url = ""
    for i in range(0, len(ncbi_accs), 50):
        ncbi_url, part = fetch_ncbi(ncbi_accs[i:i + 50])
        ncbi.update(part)
    fa = open(a.out_fasta, "w")
    mf = open(a.out_manifest, "w")
    cols = ["name", "accession", "db", "role", "expected_class", "expected_basis", "source_doc",
            "source_url", "original_header", "length", "sha256", "date_fetched", "status"]
    mf.write("\t".join(cols) + "\n")
    for r in rows:
        status = "ok"
        try:
            if r["db"] == "ncbi":
                key = r["accession"]
                if key not in ncbi and key.split(".")[0] not in ncbi:
                    raise RuntimeError("not returned by efetch")
                h, s = ncbi.get(key) or ncbi[key.split(".")[0]]
                url = ncbi_url.replace("<acc>", key)
            else:
                url, (h, s) = fetch_one(r["db"], r["accession"])
        except Exception as e:  # record the failure, keep going
            h, s, url, status = "", "", "", f"FAILED: {e}"
        s = s.replace("*", "").upper()
        if s:
            fa.write(f">{r['name']} {r['db']}:{r['accession']}\n")
            for i in range(0, len(s), 60):
                fa.write(s[i:i + 60] + "\n")
        sha = hashlib.sha256(s.encode()).hexdigest() if s else ""
        vals = [r[c] for c in cols[:7]] + [url, h.replace("\t", " "), str(len(s)), sha, today, status]
        mf.write("\t".join(vals) + "\n")
        sys.stderr.write(f"{r['name']}\t{len(s)}\t{status}\n")
    fa.close()
    mf.close()


if __name__ == "__main__":
    main()
