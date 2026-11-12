#!/usr/bin/env python3
"""What are these proteins? Turn Foldseek hits into one fold call per effector candidate.

  fold_report.py --hits pdb=hits_pdb.tsv,cath=hits_cath.tsv,afproteome=hits_af.tsv \
      --plddt folds/plddt.tsv --candidates candidates.tsv --out fold_calls.tsv

Foldseek is run with the columns
  query target theader fident alnlen qcov tcov evalue bits alntmscore prob

The call is made on **alntmscore**, not sequence identity or E-value: TM-score >= 0.5 is the
conventional same-fold threshold and is what makes this analysis worth doing - these proteins
were selected precisely because their sequences say nothing. `prob` (Foldseek's homology
probability) is reported alongside but not used as the cut, because it is calibrated on
sequence-detectable homology.

A model with mean pLDDT below --min-plddt is reported with its hits but its fold call is held
back as `low_confidence`: a bad model can hit anything at TM 0.5. A candidate with a model but
no hit above threshold is `no_match` - which for an effector set is a result, not a failure,
since novel folds are what lineage-restricted effectors are known for.
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter

COLS = ["query", "target", "theader", "fident", "alnlen", "qcov", "tcov", "evalue", "bits",
        "alntmscore", "prob"]
TM_SAME_FOLD = 0.5
MIN_PLDDT = 70.0


def parse_hits(lines):
    """Foldseek tabular lines -> list of dicts. Rows with unparsable numbers are dropped."""
    out = []
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(COLS):
            continue
        r = dict(zip(COLS, f))
        try:
            r["alntmscore"] = float(r["alntmscore"])
            r["prob"] = float(r["prob"])
            r["fident"] = float(r["fident"])
        except ValueError:
            continue
        out.append(r)
    return out


def read_hits_from_string(text):
    return parse_hits(text.splitlines())


def read_hits(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        return parse_hits(fh)


def best_per_query(hits, min_tm=TM_SAME_FOLD):
    """-> {query: best hit by TM-score}, keeping only hits at or above min_tm."""
    best = {}
    for r in hits:
        if r["alntmscore"] < min_tm:
            continue
        cur = best.get(r["query"])
        if cur is None or r["alntmscore"] > cur["alntmscore"]:
            best[r["query"]] = r
    return best


def describe(theader):
    """Shorten a Foldseek target header to something readable in a table."""
    h = (theader or "").strip()
    for sep in (" ", "\t"):
        if sep in h:
            h = h.split(sep, 1)[1] if h.split(sep, 1)[1] else h
            break
    return h[:90]


def call_folds(per_db, plddt, min_tm=TM_SAME_FOLD, min_plddt=MIN_PLDDT):
    """-> {query: row dict}. per_db is {db_name: {query: best hit}}."""
    queries = {q for d in per_db.values() for q in d}
    queries |= set(plddt)
    rows = {}
    for q in sorted(queries):
        pl = plddt.get(q, {})
        mean_plddt = float(pl.get("plddt", 0) or 0)
        hits = {db: d[q] for db, d in per_db.items() if q in d}
        row = {"protein_id": q, "plddt": mean_plddt or "", "ptm": pl.get("ptm", ""),
               "length": pl.get("length", "")}
        if not hits:
            row["fold_call"] = "no_match" if pl else "not_folded"
            row["best_db"] = row["best_target"] = row["best_desc"] = ""
            row["best_tm"] = row["best_prob"] = row["best_fident"] = ""
        else:
            db, h = max(hits.items(), key=lambda kv: kv[1]["alntmscore"])
            row["best_db"] = db
            row["best_target"] = h["target"]
            row["best_desc"] = describe(h["theader"])
            row["best_tm"] = round(h["alntmscore"], 3)
            row["best_prob"] = round(h["prob"], 3)
            row["best_fident"] = round(h["fident"], 3)
            row["fold_call"] = ("low_confidence" if mean_plddt and mean_plddt < min_plddt
                                else "assigned")
        for db in sorted(per_db):
            h = hits.get(db)
            row[f"{db}_tm"] = round(h["alntmscore"], 3) if h else ""
            row[f"{db}_target"] = h["target"] if h else ""
        rows[q] = row
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hits", required=True, help="name=path,name=path")
    ap.add_argument("--plddt")
    ap.add_argument("--candidates", help="call_effectors.py candidates.tsv, merged in if given")
    ap.add_argument("--min-tm", type=float, default=TM_SAME_FOLD)
    ap.add_argument("--min-plddt", type=float, default=MIN_PLDDT)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    per_db = {}
    for part in a.hits.split(","):
        if not part:
            continue
        name, _, path = part.partition("=")
        per_db[name] = best_per_query(read_hits(path), a.min_tm)

    plddt = {}
    if a.plddt and os.path.exists(a.plddt):
        with open(a.plddt) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                plddt[r["protein_id"]] = r

    rows = call_folds(per_db, plddt, a.min_tm, a.min_plddt)

    extra = []
    if a.candidates and os.path.exists(a.candidates):
        with open(a.candidates) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            extra = [c for c in (rd.fieldnames or [])
                     if c not in ("protein_id", "length") and c in
                     ("LOCUSTAG", "category", "tier", "fold_reason", "n_cys", "cys_frac",
                      "pfam", "cazy", "cazy_substrate")]
            fh.seek(0)
            next(fh)
            for c in csv.DictReader(fh, delimiter="\t",
                                    fieldnames=rd.fieldnames):
                if c["protein_id"] in rows:
                    for k in extra:
                        rows[c["protein_id"]][k] = c.get(k, "")

    cols = (["protein_id", "fold_call", "best_db", "best_target", "best_desc", "best_tm",
             "best_prob", "best_fident", "plddt", "ptm", "length"] + extra
            + [f"{db}_{s}" for db in sorted(per_db) for s in ("tm", "target")])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for q in sorted(rows, key=lambda x: (-float(rows[x]["best_tm"] or 0), x)):
            w.writerow(rows[q])
    c = Counter(r["fold_call"] for r in rows.values())
    print("fold calls: " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())), file=sys.stderr)
    print(f"{len(rows)} proteins -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
