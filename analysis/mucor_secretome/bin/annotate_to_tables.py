#!/usr/bin/env python3
"""Turn SignalP 6, TMHMM 2.0 and hmmsearch output into the BFD table layout.

  annotate_to_tables.py --locustag JEM1MUCOR --signalp prediction_results.txt \
      --tmhmm tmhmm_short.txt --pfam pfam.domtbl [--dbcan dbcan.domtbl] --out <dir>

Genomes outside BFD (the MorenoLab Mucor strains, the Marchantia reference) have no precomputed
annotation, so the predictors are run here and their output reshaped into the same columns
`build_secretome.py` reads from the BFD export: signalp.tsv, tmhmm.tsv, pfam.tsv, cazy.tsv,
lengths.tsv. Rows are appended when the file already exists, so several genomes can share a
directory.

Thresholds are not applied here; they belong to build_secretome.py. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

CS_RE = re.compile(r"CS pos:\s*(\d+)-(\d+)")
KV_RE = re.compile(r"(\w+)=([^\s]+)")


def parse_signalp6(handle):
    """SignalP 6 `prediction_results.txt` -> [(id, probability, cleave_start, cleave_end)].

    Columns: ID, Prediction, OTHER, SP(Sec/SPI)[, other classes...], CS Position. The signal-peptide
    probability is 1 - OTHER, which covers every SP class the model emits (Sec/SPI, Sec/SPII,
    Tat/SPI, ...) rather than only Sec/SPI."""
    out = []
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            other = float(f[2])
        except ValueError:
            continue
        m = CS_RE.search(line)
        start, end = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
        out.append((f[0].split()[0], round(1.0 - other, 6), start, end))
    return out


def parse_tmhmm_short(handle):
    """TMHMM 2.0 `--short` -> [(id, len, ExpAA, First60, PredHel, Topology)]."""
    out = []
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t") if "\t" in line else line.split()
        kv = dict(KV_RE.findall(line))
        if "PredHel" not in kv:
            continue
        out.append((f[0].split()[0], int(kv.get("len", 0)), float(kv.get("ExpAA", 0)),
                    float(kv.get("First60", 0)), int(kv["PredHel"]), kv.get("Topology", "")))
    return out


def parse_domtbl(handle):
    """hmmsearch --domtblout -> [(target, query_name, query_acc, i_evalue, ali_from, ali_to)].

    With hmmsearch the target is the protein and the query is the HMM, which is the opposite of
    hmmscan. Columns: 0 target, 3 query name, 4 query accession, 12 i-Evalue, 17 ali_from,
    18 ali_to."""
    out = []
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 19:
            continue
        out.append((f[0], f[3], f[4], float(f[12]), int(f[17]), int(f[18])))
    return out


def fam_from_hmm(name):
    """dbCAN HMM name ('GH18.hmm', 'CBM50') -> family string."""
    return re.sub(r"\.hmm$", "", (name or "").strip())


def append(path, header, rows):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if new:
            w.writerow(header)
        w.writerows(rows)
    return len(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--locustag", required=True)
    ap.add_argument("--signalp", required=True)
    ap.add_argument("--tmhmm", required=True)
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--dbcan")
    ap.add_argument("--lengths", help="TSV of protein_id<TAB>length, from the proteome")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    D = lambda f: os.path.join(a.out, f)  # noqa: E731
    lt = a.locustag

    with open(a.signalp) as fh:
        sp = parse_signalp6(fh)
    n = append(D("signalp.tsv"),
               ["protein_id", "species_prefix", "peptide_start", "peptide_end", "probability"],
               [(p, lt, s, e, prob) for p, prob, s, e in sp])
    with open(a.tmhmm) as fh:
        tm = parse_tmhmm_short(fh)
    append(D("tmhmm.tsv"), ["protein_id", "len", "ExpAA", "First60", "PredHel", "Topology"], tm)
    with open(a.pfam) as fh:
        pf = parse_domtbl(fh)
    append(D("pfam.tsv"),
           ["protein_id", "pfam_acc", "pfam_id", "domain_i_evalue", "ali_from", "ali_to"],
           [(t, (acc.split(".")[0] if acc and acc != "-" else q), q, ev, af, at)
            for t, q, acc, ev, af, at in pf])
    if a.dbcan:
        with open(a.dbcan) as fh:
            cz = parse_domtbl(fh)
        append(D("cazy.tsv"),
               ["protein_id", "species_prefix", "EC", "cazyme_fam", "sub_fam", "substrate",
                "toolcount"],
               [(t, lt, "", fam_from_hmm(q), "", "", 1) for t, q, _a, ev, _f, _o in cz
                if ev < 1e-5])
    if a.lengths:
        with open(a.lengths) as fh:
            rows = [l.rstrip("\n").split("\t")[:2] for l in fh if l.strip()]
        append(D("lengths.tsv"), ["protein_id", "transcript_id", "length"],
               [(r[0], r[0], r[1]) for r in rows if len(r) == 2])
    print(f"{lt}: {n} signal-peptide rows, {len(tm)} tmhmm, {len(pf)} pfam -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
