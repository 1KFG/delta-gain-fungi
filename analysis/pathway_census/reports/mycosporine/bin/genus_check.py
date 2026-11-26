#!/usr/bin/env python3
"""Per-genus check: how many genomes of genera named in the literature have a call for each mycosporine step.

Reads genomes.tsv.zst (locustag, genus, species, has_proteins) and the two calls tables of make_myco_calls.sh.
Writes a TSV. This is a descriptive check against reported producers; it is not a sensitivity or specificity test, because
the genus lists come from reports on single strains and BFD genomes of a genus are other strains.
"""
from __future__ import annotations

import argparse
import csv
import io
import subprocess
from collections import defaultdict
from pathlib import Path

DEFAULT = ["Phaffia", "Cystofilobasidium", "Cystobasidium", "Rhodotorula", "Rhodosporidiobolus", "Sporobolomyces",
           "Naganishia", "Hortaea", "Colletotrichum", "Pyrenophora", "Parastagonospora", "Pyricularia"]
SPECIES = ["Phaffia rhodozyma", "Hortaea werneckii", "Colletotrichum graminicola", "Cystobasidium minutum",
           "Cystobasidium slooffiae", "Naganishia friedmannii"]


def read(path):
    p = str(path)
    raw = subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout.decode() if p.endswith(".zst") else Path(p).read_text()
    return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--calls-hits", required=True, type=Path)
    ap.add_argument("--calls-placed", required=True, type=Path)
    ap.add_argument("--genera", nargs="*", default=DEFAULT)
    ap.add_argument("--species", nargs="*", default=SPECIES, help="species names (prefix match on the species column), reported as rows `species: <name>`")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    H = {r["locustag"]: r for r in read(a.calls_hits)}
    P = {r["locustag"]: r for r in read(a.calls_placed)}
    agg = defaultdict(lambda: defaultdict(int))
    for r in read(a.genomes):
        if r["has_proteins"] != "True":
            continue
        keys = ([r["genus"]] if r["genus"] in a.genera else []) + \
               [f"species: {s}" for s in a.species if r["species"].startswith(s)]
        for k in keys:
            h, p, x = H[r["locustag"]], P[r["locustag"]], agg[k]
            x["n"] += 1
            x["ddgs_placed"] += p["ddgs"] == "1"
            x["omt_placed"] += p["omt"] == "1"
            x["omt_hit"] += h["omt"] == "1"
            x["atp_grasp_ligase_hit"] += h["atp_grasp_ligase"] == "1"
            x["all_three_calls"] += p["ddgs"] == "1" and p["omt"] == "1" and h["atp_grasp_ligase"] == "1"
    cols = ["genus", "n", "ddgs_placed", "omt_placed", "omt_hit", "atp_grasp_ligase_hit", "all_three_calls"]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for g in a.genera + [f"species: {s}" for s in a.species]:
            if g in agg:
                w.writerow([g] + [agg[g][c] for c in cols[1:]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
