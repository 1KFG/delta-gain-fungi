#!/usr/bin/env python3
"""Make the chunk manifest for the batched BFD search.

Reads results/genomes.tsv.zst (pathway_census.genomes) and seqdb/genomes.tsv. Keeps genomes with
has_proteins == True. Shuffles them with a fixed seed so every chunk has a mixed set of genomes.
Writes TSV (zstd): chunk_id, locustag, proteins_path. Genomes without proteins are not listed here;
they are counted and written to --excluded so they can never read as absences.
"""
from __future__ import annotations

import argparse
import csv
import io
import random
import subprocess
import sys


def read_tsv(path):
    if str(path).endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout.decode()
        return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def make_chunks(genomes, seqdb_paths, size, seed):
    """genomes: dict rows with locustag, has_proteins. seqdb_paths: locustag -> proteins path."""
    if size < 1:
        raise ValueError("size must be >= 1")
    keep, excluded = [], []
    for g in genomes:
        ok = str(g["has_proteins"]).strip().lower() == "true"
        path = seqdb_paths.get(g["locustag"])
        if ok and path:
            keep.append((g["locustag"], path))
        else:
            excluded.append((g["locustag"], "has_proteins=" + str(g["has_proteins"]) + ("" if path else ";no_seqdb_path")))
    keep.sort()
    random.Random(seed).shuffle(keep)
    rows = [(i // size, lt, p) for i, (lt, p) in enumerate(keep)]
    return rows, excluded


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--seqdb", default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb/genomes.tsv")
    ap.add_argument("--size", type=int, required=True, help="genomes per chunk")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True, help="chunk manifest, .tsv.zst")
    ap.add_argument("--excluded", required=True)
    a = ap.parse_args(argv)
    genomes = read_tsv(a.genomes)
    seqdb = {}
    with open(a.seqdb) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            seqdb[r[0]] = r[2]
    rows, excluded = make_chunks(genomes, seqdb, a.size, a.seed)
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(["chunk_id", "locustag", "proteins_path"])
    w.writerows(rows)
    subprocess.run(["zstd", "-q", "-T0", "-f", "-o", a.out], input=buf.getvalue().encode(), check=True)
    with open(a.excluded, "w") as fh:
        fh.write("locustag\treason\n")
        for lt, why in excluded:
            fh.write(f"{lt}\t{why}\n")
    n_chunks = (rows[-1][0] + 1) if rows else 0
    print(f"genomes listed {len(rows)}, excluded {len(excluded)}, chunks {n_chunks} (size {a.size}, seed {a.seed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
