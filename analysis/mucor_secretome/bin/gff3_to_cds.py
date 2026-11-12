#!/usr/bin/env python3
"""CDS segments and transcript rows from a funannotate GFF3, in the BFD table column order.

  gff3_to_cds.py --gff3 X.gff3 --asmid JEM1 --set jem --out-cds cds.tsv --out-models models.tsv

`--out-cds` has transcript_id, chrom, start, end, strand (1-based inclusive, strand as 1/-1), the
columns `check_codons.py --cds` expects. `--out-models` has set, transcript_id, ASMID,
has_start_codon, has_stop_codon, written as NA because these genomes are not in BFD and have no
flags; the point of the check is to read the codons from the sequence instead.
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attrs(field):
    out = {}
    for part in field.strip().rstrip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_gff3(handle):
    """-> {transcript_id: (chrom, strand, [(start, end)])}, CDS segments in file order."""
    cds = defaultdict(list)
    meta = {}
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = attrs(f[8])
        parent = a.get("Parent")
        if not parent:
            continue
        for t in parent.split(","):
            cds[t].append((int(f[3]), int(f[4])))
            meta[t] = (f[0], 1 if f[6] != "-" else -1)
    return {t: (meta[t][0], meta[t][1], sorted(segs)) for t, segs in cds.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--asmid", required=True)
    ap.add_argument("--set", default="jem")
    ap.add_argument("--out-cds", required=True)
    ap.add_argument("--out-models", required=True)
    a = ap.parse_args(argv)

    with opener(a.gff3) as fh:
        tx = parse_gff3(fh)
    with open(a.out_cds, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["transcript_id", "chrom", "start", "end", "strand"])
        for t, (chrom, strand, segs) in sorted(tx.items()):
            for s, e in segs:
                w.writerow([t, chrom, s, e, strand])
    with open(a.out_models, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "transcript_id", "ASMID", "has_start_codon", "has_stop_codon"])
        for t in sorted(tx):
            w.writerow([a.set, t, a.asmid, "NA", "NA"])
    print(f"{len(tx)} transcripts from {a.gff3}", file=sys.stderr)


if __name__ == "__main__":
    main()
