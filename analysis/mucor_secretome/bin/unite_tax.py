#!/usr/bin/env python3
"""UNITE QIIME release helpers: taxonomy lives in a file of its own, not in the FASTA headers.

  unite_tax.py subset --fasta sh_refs...fasta --tax sh_taxonomy...txt --phylum Mucoromycota \
      --out mucoro.fa
  unite_tax.py rank --blast hits.tsv --tax sh_taxonomy...txt --out ranked.tsv

`sh_refs_qiime_ver10_*.fasta` headers are bare accessions (`>SH1227328.10FU_MT153946_refs`);
the lineage is in `sh_taxonomy_qiime_ver10_*.txt`, keyed by that accession
(`k__Fungi;p__Mucoromycota;...;s__Mucor_hiemalis;sh__SH...`). Filtering the FASTA on a taxon
name therefore matches nothing, which is what an earlier run of slurm/its_identify.sh did.

`rank` reads blastn -outfmt 6 with columns
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen,
where the query is a UNITE reference and the subject is a contig of the assembly being
identified, and reports the best hit per species. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import sys

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore", "qlen"]


def parse_taxonomy(handle):
    """UNITE taxonomy TSV -> {accession: lineage}. The 'Feature ID\tTaxon' header is skipped."""
    tax = {}
    for line in handle:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or f[0] == "Feature ID":
            continue
        tax[f[0]] = f[1]
    return tax


def rank_of(lineage, letter):
    """Pull one rank out of a UNITE lineage: rank_of(l, 's') -> 'Mucor_hiemalis', '' if absent."""
    for part in (lineage or "").split(";"):
        part = part.strip()
        if part.startswith(letter + "__"):
            return part[len(letter) + 2:]
    return ""


def ids_in_phylum(tax, phylum):
    return ids_at_rank(tax, "p", phylum)


def ids_at_rank(tax, letter, taxon):
    """Accessions whose lineage has `taxon` at rank `letter` ('p' phylum, 'g' genus, ...)."""
    return {acc for acc, lin in tax.items() if rank_of(lin, letter) == taxon}


def filter_fasta(src, ids, out):
    """Copy records whose header accession is in `ids`; returns the number written."""
    keep, n = False, 0
    for line in src:
        if line.startswith(">"):
            keep = line[1:].split()[0] in ids
            n += keep
        if keep:
            out.write(line)
    return n


def rank_hits(rows, tax):
    """-> [(species, lineage, pident, length, qcov, accession, bitscore)], best bitscore first.

    One row per species: the single highest-scoring reference of that species. qcov is the
    aligned fraction of the reference, so a short ITS fragment in the assembly is visible as a
    low coverage rather than hidden behind a high identity."""
    best = {}
    for r in rows:
        acc = r["qseqid"]
        lin = tax.get(acc, "")
        sp = rank_of(lin, "s") or rank_of(lin, "g") or acc
        try:
            bits, pid, ln = float(r["bitscore"]), float(r["pident"]), int(r["length"])
            qlen = int(r["qlen"]) or 1
        except (ValueError, KeyError):
            continue
        cur = best.get(sp)
        if cur is None or bits > cur[6]:
            best[sp] = (sp, lin, pid, ln, round(ln / qlen, 3), acc, bits)
    return sorted(best.values(), key=lambda t: -t[6])


def read_blast(handle):
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(BLAST_COLS):
            continue
        yield dict(zip(BLAST_COLS, f))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("subset")
    s.add_argument("--fasta", required=True)
    s.add_argument("--tax", required=True)
    s.add_argument("--phylum", default="Mucoromycota")
    s.add_argument("--rank", default="p", help="rank letter the --phylum value is matched at")
    s.add_argument("--out", required=True)
    r = sub.add_parser("rank")
    r.add_argument("--blast", required=True)
    r.add_argument("--tax", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--top", type=int, default=25)
    a = ap.parse_args(argv)

    with open(a.tax) as fh:
        tax = parse_taxonomy(fh)
    if a.cmd == "subset":
        ids = ids_at_rank(tax, a.rank, a.phylum)
        with open(a.fasta) as fh, open(a.out, "w") as fo:
            n = filter_fasta(fh, ids, fo)
        print(f"{n} {a.phylum} references of {len(ids)} in the taxonomy -> {a.out}",
              file=sys.stderr)
        return 0 if n else 1
    with open(a.blast) as fh:
        ranked = rank_hits(list(read_blast(fh)), tax)
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "pident", "aln_len", "ref_coverage", "accession", "bitscore",
                    "lineage"])
        for sp, lin, pid, ln, cov, acc, bits in ranked[:a.top]:
            w.writerow([sp, pid, ln, cov, acc, bits, lin])
    print(f"{len(ranked)} species ranked -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
