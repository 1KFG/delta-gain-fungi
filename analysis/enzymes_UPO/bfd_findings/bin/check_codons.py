#!/usr/bin/env python3
"""Check BFD's start/stop-codon flags against the genome sequence.

  check_codons.py --cds cds_segments.tsv --models models.tsv --genomes-dir input_clean_genomes \
      [--proteins proteins.faa] --threads 8 --out codon_check.tsv

--models: transcript_id, ASMID, set, has_start_codon, has_stop_codon (BFD flags), one row per model.
--cds:    transcript_id, chrom, start, end, strand (one row per CDS segment; 1-based inclusive).

For each model the CDS is spliced from <genomes-dir>/<ASMID>.fa.gz (reverse-complemented on the
minus strand) and the script records the first codon, the last codon, the codon immediately
downstream of the CDS (a stop there means the annotation excludes the stop from the CDS), internal
stops, and, when --proteins is given, whether the translation equals BFD's protein. Then:

  seq_start = first codon is ATG
  seq_stop  = last codon, or the codon just after the CDS, is a stop
  verdict   = agree / flag_says_missing_but_present / flag_says_present_but_missing, per end

A model BFD flags as lacking a start or stop, whose sequence has both, is a flag error, not a
broken gene model. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODE[a + b + c] = _aa[16 * i + 4 * j + k]


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(s):
    return "".join(CODE.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def minus(strand):
    return str(strand).strip() in ("-1", "-", "0")


def splice(contig_seq, segs, strand):
    """segs: [(start, end)] 1-based inclusive. Returns (cds, downstream codon)."""
    segs = sorted(segs)
    cds = "".join(contig_seq[s - 1:e] for s, e in segs)
    if minus(strand):
        lo = segs[0][0]
        down = revcomp(contig_seq[max(0, lo - 4):lo - 1])
        return revcomp(cds).upper(), down.upper()
    hi = segs[-1][1]
    return cds.upper(), contig_seq[hi:hi + 3].upper()


def judge(flag, seq_has):
    f = (flag or "").upper()
    if f not in ("TRUE", "FALSE"):
        return "flag_null"
    if (f == "TRUE") == seq_has:
        return "agree"
    return "flag_says_missing_but_present" if seq_has else "flag_says_present_but_missing"


def assess(cds, down, protein=None):
    first, last = cds[:3], cds[-3:] if len(cds) >= 3 else ""
    body = translate(cds[:-3] if last in STOPS else cds)
    r = {"cds_len": len(cds), "first_codon": first, "last_codon": last, "downstream_codon": down,
         "seq_start": first == "ATG", "seq_stop": last in STOPS or down in STOPS,
         "stop_in_cds": last in STOPS, "internal_stops": body.count("*")}
    if protein is not None:
        r["matches_bfd_protein"] = body.rstrip("*") == protein.rstrip("*")
    return r


def read_genome(job):
    asm, path, wanted, models = job
    seqs, name, buf = {}, None, []
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in wanted:
                        seqs[name] = "".join(buf)
                    name, buf = line[1:].split()[0], []
                elif name in wanted:
                    buf.append(line.strip())
        if name in wanted:
            seqs[name] = "".join(buf)
    except OSError as e:
        print(f"warn: {path}: {e}", file=sys.stderr)
    out = []
    for m in models:
        c = seqs.get(m["chrom"])
        if c is None:
            out.append({**m, "status": "contig_missing"})
            continue
        cds, down = splice(c, m["segs"], m["strand"])
        out.append({**m, "status": "ok", **assess(cds, down, m.get("protein"))})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cds", required=True)
    ap.add_argument("--models", required=True)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--proteins")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    segs, chrom, strand = defaultdict(list), {}, {}
    for r in csv.DictReader(open(a.cds), delimiter="\t"):
        t = r["transcript_id"]
        segs[t].append((int(r["start"]), int(r["end"])))
        chrom[t], strand[t] = r["chrom"], r["strand"]
    prot = {}
    if a.proteins:
        name = None
        for line in open(a.proteins):
            if line.startswith(">"):
                name = line[1:].split()[0]
                prot[name] = []
            elif name:
                prot[name].append(line.strip())
        prot = {k: "".join(v) for k, v in prot.items()}
    by_asm = defaultdict(list)
    for r in csv.DictReader(open(a.models), delimiter="\t"):
        t = r["transcript_id"]
        if t not in segs:
            continue
        by_asm[r["ASMID"]].append({"transcript_id": t, "ASMID": r["ASMID"], "set": r["set"],
                                   "has_start_codon": r["has_start_codon"],
                                   "has_stop_codon": r["has_stop_codon"], "chrom": chrom[t],
                                   "strand": strand[t], "segs": segs[t],
                                   **({"protein": prot[t]} if t in prot else {})})
    jobs = [(asm, os.path.join(a.genomes_dir, f"{asm}.fa.gz"), {m["chrom"] for m in ms}, ms)
            for asm, ms in by_asm.items()]
    rows = []
    with ProcessPoolExecutor(max_workers=a.threads) as ex:
        for out in ex.map(read_genome, jobs):
            rows.extend(out)
    cols = ["set", "transcript_id", "ASMID", "chrom", "strand", "status", "has_start_codon",
            "has_stop_codon", "cds_len", "first_codon", "last_codon", "downstream_codon",
            "seq_start", "seq_stop", "stop_in_cds", "internal_stops", "matches_bfd_protein",
            "start_verdict", "stop_verdict"]
    for r in rows:
        if r["status"] == "ok":
            r["start_verdict"] = judge(r["has_start_codon"], r["seq_start"])
            r["stop_verdict"] = judge(r["has_stop_codon"], r["seq_stop"])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["set"], r["transcript_id"])))
    print(f"checked {len(rows)} models in {len(jobs)} genomes -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
