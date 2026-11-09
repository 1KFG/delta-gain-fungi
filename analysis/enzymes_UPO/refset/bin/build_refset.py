#!/usr/bin/env python3
"""Combine the Zamocky 2022 reference sequences with the 19-enzyme clade backbone, remove
duplicates, QC each sequence, and write the input for the in-house reference tree.

Steps 2 and 3 of the reference-set expansion (see
../../data/refset_zamocky2022/README.md). Standard library only; Python >= 3.9.

  combine   Zamocky FASTA + Table S1 + backbone -> combined.faa (tree labels) + dedup.tsv
  select    combined.faa + dedup.tsv + upo_sites output -> qc.tsv + tree_input.faa

Tree labels are  B__<name>__clade<N>  for backbone enzymes and  Z__<abbrev>__<group>  for
Zamocky sequences. Characters are restricted to [A-Za-z0-9_.-] so every tree program keeps them.

Duplicates are found by sequence, not accession: two sequences are the same protein when they are
identical or one is contained in the other (precursor vs mature chain, signal peptide trimmed).
The backbone copy is kept, because its label carries our clade call; the Zamocky abbreviation is
recorded in `dup_of` so Figure 1 clade membership can still be matched.

QC (`select`) excludes a sequence from the tree only when it has no PF01328 hit. Everything else
(catalytic residues absent, length differs from Table S1) is a flag, not an exclusion: Zamocky's
set deliberately spans divergent HTPs (oomycete, amoebozoan, chytrid), for which HspUPO-numbered
catalytic positions are not expected to hold.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

# Table S1 phylum strings, normalised. The table has spelling variants (trailing spaces,
# "Amobeozoa", three Mucoromycota spellings); this maps each to one group. Unknown strings fail.
TAXON_GROUP = {
    "ascomycota": "Asco",
    "basidiomycota": "Basidio",
    "chytridiomycota": "Chytrid",
    "mucoromycota": "Mucoro",
    "mucoromycotina": "Mucoro",
    "zygomycota_mucoromycotina": "Mucoro",
    "glomeromycotina_mucoromycota": "Glomero",
    "stramenopiles_oomycota": "Oomycota",
    "amobeozoa": "Amoebozoa",
    "amoebozoa": "Amoebozoa",
}
SAFE = re.compile(r"[^A-Za-z0-9_.-]")


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            elif line:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def clean_seq(s: str) -> str:
    return s.upper().replace("*", "").replace("-", "").replace(" ", "")


def taxon_group(raw: str) -> str:
    key = raw.strip().lower().replace(" ", "_")
    if key not in TAXON_GROUP:
        raise ValueError(f"unknown Table S1 taxon {raw!r}; add it to TAXON_GROUP")
    return TAXON_GROUP[key]


def safe(s: str) -> str:
    return SAFE.sub("_", s)


def load_table(path):
    with open(path, newline="") as fh:
        return {r["abbrev"]: r for r in csv.DictReader(fh, delimiter="\t")}


def same_protein(a: str, b: str) -> bool:
    return a == b or a in b or b in a


def combine(zam_fa, table_tsv, backbone_fa):
    """Return (records, dedup_rows). records: list of (label, seq). Backbone first, so a
    Zamocky sequence that duplicates a backbone enzyme is the one dropped."""
    table = load_table(table_tsv)
    records, rows = [], []
    kept = []  # (label, seq, set, name)
    for header, seq in read_fasta(backbone_fa):
        f = header.split("|")
        name, clade = f[0], (f[1] if len(f) > 1 else "clade?")
        label = f"B__{safe(name)}__{safe(clade)}"
        kept.append((label, clean_seq(seq), "backbone", name))
        rows.append({"label": label, "set": "backbone", "name": name, "taxon_group": "",
                     "length": len(clean_seq(seq)), "dup_of": ""})
    for header, seq in read_fasta(zam_fa):
        abbrev = header.split()[0]
        if abbrev not in table:
            raise ValueError(f"{abbrev} is in the FASTA but not in Table S1")
        grp = taxon_group(table[abbrev]["taxon"])
        label = f"Z__{safe(abbrev)}__{grp}"
        s = clean_seq(seq)
        dup = next((k for k in kept if same_protein(s, k[1])), None)
        rows.append({"label": label, "set": "zamocky", "name": abbrev, "taxon_group": grp,
                     "length": len(s), "dup_of": dup[0] if dup else ""})
        if dup is None:
            kept.append((label, s, "zamocky", abbrev))
        elif dup[2] == "backbone":
            # record the reverse link so the backbone row knows its Zamocky abbreviation
            for r in rows:
                if r["label"] == dup[0]:
                    r["dup_of"] = ";".join(x for x in (r["dup_of"], label) if x)
    records = [(k[0], k[1]) for k in kept]
    return records, rows


def parse_sites(path):
    """upo_sites TSV -> {protein_id (first token): row}. Comment lines start with '#'."""
    out = {}
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        out[r["protein_id"].split()[0]] = r
    return out


def select(combined_fa, dedup_tsv, sites_tsv, table_tsv):
    table = load_table(table_tsv) if table_tsv else {}
    sites = parse_sites(sites_tsv)
    with open(dedup_tsv, newline="") as fh:
        dedup = list(csv.DictReader(fh, delimiter="\t"))
    seqs = dict(read_fasta(combined_fa))
    qc, tree = [], []
    for r in dedup:
        row = dict(r)
        s = sites.get(r["label"])
        flags = []
        if r["label"] not in seqs:
            row.update(pf01328_evalue="", catalytic_ok="", catalytic_detail="", upo_sites_clade="",
                       nearest_ref="", nearest_identity="", in_tree="no",
                       reason=f"duplicate of {r['dup_of']}")
            qc.append(row)
            continue
        if s is None:
            row.update(pf01328_evalue="", catalytic_ok="", catalytic_detail="", upo_sites_clade="",
                       nearest_ref="", nearest_identity="", in_tree="no",
                       reason="no PF01328 hit")
            qc.append(row)
            continue
        if s.get("catalytic_ok") != "yes":
            flags.append("catalytic_residues_differ")
        if r["set"] == "zamocky" and r["name"] in table:
            tl = table[r["name"]].get("length_aa", "")
            if tl.isdigit() and int(tl) != int(r["length"]):
                flags.append(f"length_{r['length']}_vs_tableS1_{tl}")
        row.update(pf01328_evalue=s.get("hmm_evalue", ""), catalytic_ok=s.get("catalytic_ok", ""),
                   catalytic_detail=s.get("catalytic_detail", ""),
                   upo_sites_clade=s.get("clade", ""), nearest_ref=s.get("nearest_ref", ""),
                   nearest_identity=s.get("clade_identity", ""), in_tree="yes",
                   reason=";".join(flags))
        qc.append(row)
        tree.append((r["label"], seqs[r["label"]]))
    return qc, tree


def write_fasta(path, records):
    with open(path, "w") as fh:
        for label, seq in records:
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_tsv(path, rows):
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("combine")
    c.add_argument("--zamocky", required=True)
    c.add_argument("--table", required=True, help="tableS1_resolved.tsv")
    c.add_argument("--backbone", required=True, help="clade_backbone/raw.fa")
    c.add_argument("--out-fasta", required=True)
    c.add_argument("--out-dedup", required=True)
    s = sub.add_parser("select")
    s.add_argument("--combined", required=True)
    s.add_argument("--dedup", required=True)
    s.add_argument("--sites", required=True, help="upo_sites output TSV (uncompressed)")
    s.add_argument("--table", help="tableS1_resolved.tsv, for the length check")
    s.add_argument("--out-qc", required=True)
    s.add_argument("--out-fasta", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "combine":
        recs, rows = combine(a.zamocky, a.table, a.backbone)
        write_fasta(a.out_fasta, recs)
        write_tsv(a.out_dedup, rows)
        nd = sum(1 for r in rows if r["set"] == "zamocky" and r["dup_of"])
        print(f"combined {len(recs)} sequences; {nd} Zamocky duplicates of earlier sequences",
              file=sys.stderr)
    else:
        qc, tree = select(a.combined, a.dedup, a.sites, a.table)
        write_tsv(a.out_qc, qc)
        write_fasta(a.out_fasta, tree)
        out = [r for r in qc if r["in_tree"] == "no" and not r["reason"].startswith("duplicate")]
        print(f"tree input {len(tree)}; excluded (no PF01328 hit) {len(out)} "
              f"{[r['name'] for r in out]}; flagged "
              f"{sum(1 for r in qc if r['in_tree'] == 'yes' and r['reason'])}", file=sys.stderr)


if __name__ == "__main__":
    main()
