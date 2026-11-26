#!/usr/bin/env python3
"""Choose labelled non-target tree leaves (named members and outgroups) for step 5 from fetched UniProt family entries.

Inputs: the alignability screen of the family FASTA (bin/screen_outgroups.py output), the UniProt TSV tables of the same
entries (fields Entry, Organism, Taxonomic lineage, Gene Names, Protein names, Length, EC number), the FASTA, and the
step references. Rules (PROPOSED):
  1. A protein is eligible when it passes the screen (`--mode best`: best domain >= 10 bits and covering >= 50% of the
     HMM; `--mode union` uses the union coverage).
  2. A protein is excluded when it is a reference (any accession in --exclude-ids), or when it shares a complete EC number
     with a reference (--exclude-ec): it is then the same enzyme of another species and not a non-target leaf.
  2b. With --max-score F (the step floor), a protein scoring F bits or more is excluded: the search would call it a candidate,
     so it cannot be a labelled non-target leaf.
  2c. With --exclude-gene (gene-name tokens of the step references), a protein with one of those gene names is excluded: it is
     named like a reference (for example Arp1 of Metarhizium, THR1 of Sclerotinia), so it may be an ortholog, not a non-target leaf.
  3. Named members: fungal entries (lineage holds "Fungi"), best scores first, one per organism, at most --n-named.
  4. Outgroups: the other eligible entries (non-fungal first, then fungal leftovers), one per organism, taken at evenly
     spaced ranks of the score order, at most --n-outgroup. Evenly spaced ranks give near and far relatives.
Output: FASTA of named members, FASTA of outgroups (id = UniProt accession), a provenance TSV (every field from the fetched
record plus the screen numbers) and an exclusions TSV (accession, reason). Nothing is taken from memory.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_fasta(path):
    name, hdr, seq = None, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, hdr, "".join(seq)
            hdr = line[1:]
            name, seq = hdr.split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, hdr, "".join(seq)


def accession_of(token: str) -> str:
    """`sp|P12345|NAME_SPEC` or `tr|A0A..|..` -> the accession; a bare id stays as is."""
    parts = token.split("|")
    return parts[1] if len(parts) >= 3 else token


def organism_short(org: str) -> str:
    """Organism name without the strain note in parentheses: 'Aspergillus fumigatus (strain X) (Y)' -> 'Aspergillus fumigatus'."""
    return org.split(" (")[0].strip()


def eligible_rows(screen_rows, meta, mode, exclude_ids, exclude_ec, max_score=None, exclude_genes=()):
    """Return (eligible list of dicts, exclusions list of (accession, reason)). Order: best score first."""
    col = {"best": "pass_best", "union": "pass_union"}[mode]
    elig, excl = [], []
    ec_ex = {e for e in exclude_ec if e}
    for r in screen_rows:
        acc = accession_of(r["protein"])
        if r[col] != "True":
            continue
        if acc in exclude_ids:
            excl.append((acc, "reference or excluded id"))
            continue
        m = meta.get(acc)
        if m is None:
            excl.append((acc, "no UniProt table row"))
            continue
        if max_score is not None and float(r["best_score"]) >= max_score:
            excl.append((acc, f"best score {float(r['best_score']):g} is at or above the step floor {max_score:g}: it would be a candidate"))
            continue
        genes = {g.lower() for g in m.get("Gene Names", "").split()}
        if genes & {g.lower() for g in exclude_genes if g}:
            excl.append((acc, "gene name matches a reference: " + ",".join(sorted(genes & {g.lower() for g in exclude_genes}))))
            continue
        ecs = {e.strip() for e in m.get("EC number", "").split(";") if e.strip()}
        if ecs & ec_ex:
            excl.append((acc, "same EC as a reference: " + ",".join(sorted(ecs & ec_ex))))
            continue
        elig.append(dict(accession=acc, best_score=float(r["best_score"]), best_cov=float(r["best_cov"]),
                         union_cov=float(r["union_cov"]), organism=m["Organism"], gene=m.get("Gene Names", ""),
                         name=m.get("Protein names", ""), length=int(m["Length"]), ec=m.get("EC number", ""),
                         fungal="Fungi" in m.get("Taxonomic lineage", "")))
    elig.sort(key=lambda d: (-d["best_score"], d["accession"]))
    return elig, excl


def choose(elig, n_named, n_outgroup):
    """Return (named list, outgroup list). One entry per organism across both lists."""
    used, named = set(), []
    for d in elig:
        o = organism_short(d["organism"])
        if d["fungal"] and o not in used and len(named) < n_named:
            named.append(d)
            used.add(o)
    pool, seen = [], set(used)
    for d in [x for x in elig if not x["fungal"]] + [x for x in elig if x["fungal"]]:
        o = organism_short(d["organism"])
        if d["accession"] in {n["accession"] for n in named} or o in seen:
            continue
        seen.add(o)
        pool.append(d)
    # pool order: non-fungal by score, then fungal leftovers; spread evenly over the pool
    if n_outgroup <= 0 or not pool:
        return named, []
    if len(pool) <= n_outgroup:
        return named, pool
    step = (len(pool) - 1) / (n_outgroup - 1) if n_outgroup > 1 else 0
    idx = sorted({round(i * step) for i in range(n_outgroup)})
    return named, [pool[i] for i in idx]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--screen", required=True, type=Path)
    ap.add_argument("--meta", required=True, nargs="+", type=Path, help="UniProt TSV tables")
    ap.add_argument("--fasta", required=True, nargs="+", type=Path, help="FASTA files of the same entries")
    ap.add_argument("--hmm", required=True, help="HMM name in the screen table")
    ap.add_argument("--mode", choices=["best", "union"], default="best")
    ap.add_argument("--exclude-ids", default="", help="comma list of accessions (the references)")
    ap.add_argument("--exclude-ec", default="", help="comma list of complete EC numbers of the references")
    ap.add_argument("--max-score", type=float, help="step floor: exclude proteins scoring at or above it")
    ap.add_argument("--exclude-gene", default="", help="comma list of gene-name tokens of the references (case-insensitive)")
    ap.add_argument("--n-named", type=int, default=10)
    ap.add_argument("--n-outgroup", type=int, default=6)
    ap.add_argument("--out-named", required=True, type=Path)
    ap.add_argument("--out-outgroups", required=True, type=Path)
    ap.add_argument("--out-provenance", required=True, type=Path)
    ap.add_argument("--out-exclusions", required=True, type=Path)
    a = ap.parse_args(argv)
    meta = {}
    for p in a.meta:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                meta[r["Entry"]] = r
    seqs = {}
    for p in a.fasta:
        for name, hdr, seq in read_fasta(p):
            seqs[accession_of(name)] = (hdr, seq)
    with open(a.screen) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["hmm"] == a.hmm]
    ex_ids = {x for x in a.exclude_ids.split(",") if x}
    ex_ec = [x for x in a.exclude_ec.split(",") if x]
    elig, excl = eligible_rows(rows, meta, a.mode, ex_ids, ex_ec, a.max_score, [x for x in a.exclude_gene.split(",") if x])
    named, outg = choose(elig, a.n_named, a.n_outgroup)
    for path, lst in ((a.out_named, named), (a.out_outgroups, outg)):
        with open(path, "w") as fh:
            for d in lst:
                hdr, seq = seqs[d["accession"]]
                if len(seq) != d["length"]:
                    sys.exit(f"ERROR: length mismatch for {d['accession']}: FASTA {len(seq)}, table {d['length']}")
                fh.write(f">{d['accession']} {hdr.split(' ', 1)[1] if ' ' in hdr else ''}\n{seq}\n")
    cols = ["role", "accession", "organism", "gene", "name", "length", "ec", "best_score", "best_cov", "union_cov", "fungal"]
    with open(a.out_provenance, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for role, lst in (("named_member", named), ("tree_outgroup", outg)):
            for d in lst:
                w.writerow([role] + [d[c] for c in cols[1:]])
    with open(a.out_exclusions, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "reason"])
        w.writerows(excl)
    print(f"{a.hmm}: eligible {len(elig)}, excluded {len(excl)}, named {len(named)}, outgroups {len(outg)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
