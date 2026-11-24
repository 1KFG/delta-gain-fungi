#!/usr/bin/env python3
"""Check that the protein id prefix of each genome matches the genome table.

For every genome in the genome table with has_proteins == True, the script opens the protein file
named in seqdb/genomes.tsv and reads the FIRST protein id. The prefix is the text before the first
underscore. The census (bin/search_chunk.py) counts a protein as `bad id prefix` when its id does not
start with `<locustag>_`. This script finds the genomes where that happens, one protein per genome.

Limits. Only the first protein of each file is read. A file that mixes two prefixes is not detected.
The mtime test (file newer than --since) is a hint about why a file changed. It is not proof.

Inputs may be plain, .gz or .zst (genome table, seqdb table; protein files: plain or .gz).
Outputs: a TSV with one row per genome (--out), a summary on stdout, and optionally (--aliases-out)
a prefix_aliases.tsv with one row per mismatching genome that has a readable first id.
Threads: at most 8 (login node rule).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

MAX_THREADS = 8
OUT_COLS = ["locustag", "proteins_path", "mtime", "newer_than_table", "first_id", "first_prefix", "status"]


def read_text(path) -> str:
    path = str(path)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return Path(path).read_text()


def read_genomes(path):
    """locustags with has_proteins true, in table order."""
    rows = csv.DictReader(io.StringIO(read_text(path)), delimiter="\t")
    return [r["locustag"] for r in rows if str(r["has_proteins"]).strip().lower() == "true"]


def read_seqdb(path):
    """locustag -> proteins path (seqdb/genomes.tsv: locustag, out, proteins, cds; no header)."""
    out = {}
    for r in csv.reader(io.StringIO(read_text(path)), delimiter="\t"):
        if len(r) >= 3:
            out[r[0]] = r[2]
    return out


def first_protein_id(path):
    """First FASTA id (text after `>` up to the first space), or None for an empty file. Raises OSError/EOFError."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0]
    return None


def check_one(locustag, path, since_ts):
    row = dict(locustag=locustag, proteins_path=path or "", mtime="", newer_than_table="", first_id="", first_prefix="", status="")
    if not path:
        row["status"] = "no_seqdb_path"
        return row
    try:
        mt = Path(path).stat().st_mtime
        row["mtime"] = datetime.fromtimestamp(mt).strftime("%Y-%m-%dT%H:%M:%S")
        row["newer_than_table"] = "" if since_ts is None else str(mt > since_ts)
        fid = first_protein_id(path)
    except (OSError, EOFError) as e:
        row["status"] = f"missing:{type(e).__name__}"
        return row
    if fid is None:
        row["status"] = "empty"
        return row
    row["first_id"] = fid
    row["first_prefix"] = fid.split("_")[0]
    row["status"] = "ok" if fid.startswith(locustag + "_") else "mismatch"
    return row


def check_all(locustags, seqdb, since_ts, threads=4):
    threads = max(1, min(int(threads), MAX_THREADS))
    with ThreadPoolExecutor(max_workers=threads) as ex:
        return list(ex.map(lambda lt: check_one(lt, seqdb.get(lt), since_ts), locustags))


def summarise(rows):
    n = len(rows)
    by = {}
    for r in rows:
        by[r["status"]] = by.get(r["status"], 0) + 1
    newer = [r for r in rows if r["newer_than_table"] == "True"]
    return dict(n_checked=n, by_status=by, n_newer=len(newer),
                n_newer_mismatch=sum(1 for r in newer if r["status"] == "mismatch"),
                n_older_mismatch=sum(1 for r in rows if r["newer_than_table"] == "False" and r["status"] == "mismatch"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True, help="genome table (genomes.tsv.zst)")
    ap.add_argument("--seqdb", default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb/genomes.tsv")
    ap.add_argument("--since", required=True, help="table build time, ISO local time, e.g. 2026-11-22T12:11:00")
    ap.add_argument("--out", required=True, help="per-genome report TSV")
    ap.add_argument("--aliases-out", default=None, help="write prefix_aliases.tsv for mismatching genomes")
    ap.add_argument("--threads", type=int, default=4, help=f"at most {MAX_THREADS}")
    a = ap.parse_args(argv)
    since_ts = datetime.fromisoformat(a.since).timestamp()
    rows = check_all(read_genomes(a.genomes), read_seqdb(a.seqdb), since_ts, a.threads)
    write_tsv(a.out, rows, OUT_COLS)
    if a.aliases_out:
        al = [dict(locustag=r["locustag"], alias_prefix=r["first_prefix"], proteins_path=r["proteins_path"], file_mtime=r["mtime"])
              for r in rows if r["status"] == "mismatch"]
        write_tsv(a.aliases_out, al, ["locustag", "alias_prefix", "proteins_path", "file_mtime"])
    s = summarise(rows)
    print(f"checked {s['n_checked']}; status {s['by_status']}")
    print(f"files newer than {a.since}: {s['n_newer']}; of these mismatch: {s['n_newer_mismatch']}; "
          f"older files that mismatch: {s['n_older_mismatch']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
