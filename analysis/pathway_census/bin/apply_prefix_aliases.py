#!/usr/bin/env python3
"""Write a corrected genome table: the original table plus one column `alias_prefix`.

alias_prefix is the protein id prefix that the BFD protein file of the genome uses now, when it differs
from `locustag`. It is empty for all other genomes. The original table is not changed and no other
column is changed (taxonomy, BUSCO and N50 stay those of the original table row).
Aliases file: TSV with columns locustag, alias_prefix (more columns are allowed and ignored).
Every aliased locustag must be in the table, and an alias prefix must not be a locustag of the table.
Inputs and output may be plain, .gz or .zst (output by extension: .zst or plain).
"""
from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from pathlib import Path


def read_text(path) -> str:
    path = str(path)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    if path.endswith(".gz"):
        import gzip
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return Path(path).read_text()


def apply_aliases(genomes, aliases):
    """genomes: list of dict rows. aliases: list of dict rows (locustag, alias_prefix). Returns new rows."""
    table = {g["locustag"] for g in genomes}
    amap = {}
    for a in aliases:
        lt, pre = a["locustag"], a["alias_prefix"]
        if lt not in table:
            raise ValueError(f"aliased locustag {lt} is not in the genome table")
        if pre in table:
            raise ValueError(f"alias prefix {pre} is already a locustag of the genome table")
        if lt in amap and amap[lt] != pre:
            raise ValueError(f"locustag {lt} has two alias prefixes")
        amap[lt] = pre
    if len(set(amap.values())) != len(amap):
        raise ValueError("two locustags share one alias prefix")
    return [{**g, "alias_prefix": amap.get(g["locustag"], "")} for g in genomes], amap


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--aliases", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    gen = list(csv.DictReader(io.StringIO(read_text(a.genomes)), delimiter="\t"))
    if gen and "alias_prefix" in gen[0]:
        raise SystemExit("input table already has an alias_prefix column")
    ali = list(csv.DictReader(io.StringIO(read_text(a.aliases)), delimiter="\t"))
    rows, amap = apply_aliases(gen, ali)
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    if a.out.endswith(".zst"):
        subprocess.run(["zstd", "-q", "-T0", "-f", "-o", a.out], input=buf.getvalue().encode(), check=True)
    else:
        Path(a.out).write_text(buf.getvalue())
    print(f"rows {len(rows)}; aliased genomes {len(amap)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
