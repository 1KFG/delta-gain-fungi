#!/usr/bin/env python3.12
"""Fetch reviewed UniProt entries of one Pfam family for step 5 leaf candidates (one-off build helper, not tested code).

For each Pfam id it fetches (a) reviewed fungal entries (taxonomy 4751) and (b) reviewed non-fungal entries, at most
--max-other of the second set (UniProt order; the cap is a size choice, not a sample design). It fetches the TSV table
and the FASTA with the same query, and checks that both give the same accessions and that the FASTA length equals the
`Length` column. Nothing is taken from memory: every sequence and every field comes from rest.uniprot.org.
Output: <out>/<PF>.<fungal|other>.tsv and .faa, and <out>/fetch_log.tsv (query, date, counts, checks).
Run from analysis/pathway_census:
  /usr/bin/python3.12 curation_scripts/melanin_myco_step5/fetch_family.py --out results/step5/other_inputs/family PF13561 PF06500 ...
"""
import argparse, csv, io, sys, urllib.parse, urllib.request
from datetime import date
from pathlib import Path

UA = {"User-Agent": "pathway-census-curation"}
FIELDS = "accession,organism_name,lineage,gene_names,protein_name,length,ec,xref_pfam"
BASE = "https://rest.uniprot.org/uniprotkb/search"


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read().decode()


def fetch(query, fmt, size):
    q = urllib.parse.urlencode({"query": query, "format": fmt, "size": size, **({"fields": FIELDS} if fmt == "tsv" else {})})
    return get(f"{BASE}?{q}")


def fasta_records(text):
    name, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--max-other", type=int, default=500)
    ap.add_argument("pfams", nargs="+")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    log = []
    for pf in a.pfams:
        for kind, tail, size in (("fungal", "AND taxonomy_id:4751", 500), ("other", "AND NOT taxonomy_id:4751", a.max_other)):
            query = f"reviewed:true AND xref:pfam-{pf} {tail}"
            tsv = fetch(query, "tsv", size)
            faa = fetch(query, "fasta", size)
            rows = list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
            recs = list(fasta_records(faa))
            acc_t = [r["Entry"] for r in rows]
            acc_f = [h.split("|")[1] for h, _ in recs]
            length_ok = all(int(r["Length"]) == len(s) for r, (_, s) in zip(rows, recs)) if acc_t == acc_f else False
            ok = (acc_t == acc_f) and length_ok
            (a.out / f"{pf}.{kind}.tsv").write_text(tsv)
            (a.out / f"{pf}.{kind}.faa").write_text(faa)
            log.append([pf, kind, query, date.today().isoformat(), len(rows), len(recs), acc_t == acc_f, length_ok])
            print(pf, kind, len(rows), len(recs), "OK" if ok else "CHECK FAILED", file=sys.stderr)
            if not ok:
                sys.exit(1)
    with open(a.out / "fetch_log.tsv", "a") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if fh.tell() == 0:
            w.writerow(["pfam", "set", "query", "date", "n_table", "n_fasta", "same_accessions", "lengths_equal"])
        w.writerows(log)


if __name__ == "__main__":
    main()
