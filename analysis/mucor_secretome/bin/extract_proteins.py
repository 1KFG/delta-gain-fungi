#!/usr/bin/env python3
"""Pull protein sequences out of the BFD seqdb for a set of genomes or protein IDs.

  extract_proteins.py --seqdb-index seqdb/genomes.tsv --locustags lt.txt \
      [--ids wanted_ids.txt] --out proteins.faa [--threads 8]

`seqdb/genomes.tsv` has no header: LOCUSTAG, out, proteins.fa.gz path, cds.fa.gz path. Headers in
those files are `>{transcript_id} {gene_id}`, so the first field is the ID used everywhere else.
Sequences are written one line per record, in the order the genomes are listed. Standard library
only.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor


def read_index(path):
    """LOCUSTAG -> protein fasta path."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0]:
                out[f[0]] = f[2]
    return out


def lines(path):
    return open(path) if not str(path).endswith(".gz") else gzip.open(path, "rt")


def read_fasta(path, wanted=None):
    """Yield (id, seq). `wanted` None keeps everything."""
    name, buf = None, []
    with lines(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        yield name, "".join(buf)


def one(job):
    lt, path, wanted = job
    if not path or not os.path.exists(path):
        return lt, None, []
    try:
        return lt, None, list(read_fasta(path, wanted))
    except OSError as e:
        return lt, str(e), []


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--seqdb-index", required=True)
    ap.add_argument("--locustags", required=True, help="file with one LOCUSTAG per line")
    ap.add_argument("--ids", help="file with one protein ID per line; default keeps every protein")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    index = read_index(a.seqdb_index)
    lts = [l.strip() for l in open(a.locustags) if l.strip()]
    wanted = None
    if a.ids:
        wanted = {l.strip() for l in open(a.ids) if l.strip()}
    by_lt = None
    if wanted is not None:
        by_lt = {}
        for p in wanted:
            by_lt.setdefault(p.split("_")[0], set()).add(p)

    jobs = [(lt, index.get(lt), None if by_lt is None else by_lt.get(lt, set())) for lt in lts]
    jobs = [j for j in jobs if j[2] is None or j[2]]
    n, missing, failed = 0, [], []
    with open(a.out, "w") as fo, ProcessPoolExecutor(max_workers=a.threads) as ex:
        for lt, err, recs in ex.map(one, jobs):
            if err:
                failed.append(f"{lt}: {err}")
            if not recs and not err:
                missing.append(lt)
            for pid, seq in recs:
                fo.write(f">{pid}\n{seq}\n")
                n += 1
    if missing:
        print(f"warn: no sequences for {len(missing)} genomes (first: {missing[:3]})",
              file=sys.stderr)
    for f in failed:
        print(f"warn: {f}", file=sys.stderr)
    print(f"wrote {n} proteins from {len(jobs)} genomes -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
