#!/usr/bin/env python3
"""Turn the BFD predictor outputs into the BFD table layout, for genomes BFD has not run.

  annotate_to_tables.py --locustag JEM1MUCOR --out <dir> \
      --signalp-gff3 output.gff3 --tmhmm tmhmm_short.txt --predgpi x.predgpi.gff3 \
      --targetp x_summary.targetp2 --pfam pfam.domtbl --cazy-overview overview.tsv \
      --lengths lengths.tsv

Genomes outside BFD (the MorenoLab Mucor strains, the Marchantia reference) have no precomputed
annotation, so the predictors are run here and their output reshaped into the same columns
`build_secretome.py` reads from the BFD export: signalp.tsv, tmhmm.tsv, gpi.tsv, targetp.tsv,
pfam.tsv, cazy.tsv, lengths.tsv. Every input is optional, so the GPU SignalP job and the CPU
job can each write the tables they produce. Rows are appended when the file already exists, so
several genomes can share a directory.

Each parser mirrors one of Fungi_BFD/nextflow/bin/merge_*.py, so a table built here means the
same thing as the BFD parquet of that name:
  signalp  - from the GFF3, like merge_signalp.py: only signal-peptide-positive proteins, and
             `probability` is the GFF3 score. parse_signalp6() reads prediction_results.txt
             instead and is the fallback for a run made without GFF3 output; it scores
             1 - OTHER, which is the same quantity for an SP-positive protein.
  tmhmm    - merge_tmhmm.py drops proteins with 0 predicted helices; they are kept here,
             because build_secretome.py reads a missing row and `PredHel=0` identically and
             the row records that TMHMM did run on that protein.
  predgpi  - merge_predgpi.py keeps every GFF3 row, and PredGPI writes a `Chain` row for every
             protein it examines. Only `GPI-anchor` rows are written here (as sql/export.sql
             also does), because that is what "GPI-anchored" means.
  targetp  - merge_targetp.py drops noTP proteins, and `probability` belongs to whichever class
             was predicted, so it is an mTP probability only when prediction is mTP.
Thresholds are not applied here; they belong to build_secretome.py. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

CS_RE = re.compile(r"CS pos:\s*(\d+)-(\d+)")
TP_CS_RE = re.compile(r"CS pos:\s+(\d+)-(\d+)\.\s+(\S+)\.\s+Pr:\s+(\S+)")
KV_RE = re.compile(r"(\w+)=([^\s]+)")


def rows_of(handle):
    """GFF3/TSV data lines as split fields, comments and blanks dropped."""
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        yield line.rstrip("\n").split("\t")


def parse_signalp_gff3(handle):
    """SignalP 6 `output.gff3` -> [(id, score, start, end)], one row per SP-positive protein."""
    out = []
    for f in rows_of(handle):
        if len(f) < 6:
            continue
        out.append((f[0].split(" ")[0], f[5], f[3], f[4]))
    return out


def parse_signalp6(handle):
    """SignalP 6 `prediction_results.txt` -> [(id, probability, cleave_start, cleave_end)].

    Columns: ID, Prediction, OTHER, SP(Sec/SPI)[, other classes...], CS Position. The signal-
    peptide probability is 1 - OTHER, which covers every SP class the model emits (Sec/SPI,
    Sec/SPII, Tat/SPI, ...) rather than only Sec/SPI."""
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


def parse_predgpi_gff3(handle):
    """PredGPI GFF3 -> [(id, feature, start, end)] for the GPI-anchor rows only.

    PredGPI emits a `Chain` row for every protein it examined and a `GPI-anchor` row only for
    the anchored ones (224 M against 3.3 M rows across BFD), so dropping the feature column
    turns "PredGPI ran on this protein" into "this protein is GPI-anchored"."""
    return [(f[0], f[2], f[3], f[4]) for f in rows_of(handle)
            if len(f) >= 5 and f[2] == "GPI-anchor"]


def parse_targetp(handle):
    """TargetP 2.0 `-format short` summary -> [(id, prediction, probability, cs_start, cs_end)].

    `probability` is the score of the predicted class, so it is an mTP probability only when
    prediction is mTP. noTP proteins are dropped, as merge_targetp.py drops them."""
    out = []
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t", 2)
        if len(parts) < 3:
            continue
        pid, prediction, rest = parts
        if prediction == "noTP":
            continue
        fields = rest.split("\t")
        sp = fields[1] if len(fields) > 1 else ""
        mtp = fields[2] if len(fields) > 2 else ""
        cs = fields[3] if len(fields) > 3 else ""
        prob = sp if prediction == "SP" else (mtp if prediction == "mTP" else "")
        m = TP_CS_RE.match(cs)
        cs_start, cs_end = (m.group(1), m.group(2)) if m else ("0", "0")
        out.append((pid.split()[0], prediction, prob, cs_start, cs_end))
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


def parse_cazy_overview(handle):
    """dbcanlight `conclude` overview.tsv -> [(id, EC, cazyme_fam, sub_fam, substrate, tools)].

    Columns: Gene_ID, EC, cazyme_fam, sub_fam, diamond_fam, Substrate, #ofTools. One output row
    per family, because build_secretome.py reads one family per row. '-' means absent."""
    out = []
    first = True
    for f in rows_of(handle):
        if first:
            first = False
            if f and f[0].strip().lower() in ("gene_id", "gene id"):
                continue
        if len(f) < 7:
            continue
        clean = [("" if x.strip() == "-" else x.strip()) for x in f]
        gene, ec, fams, subs, _dia, substrate, tools = clean[:7]
        # a protein with several families lists them as 'GH18+CBM18', and the matching sub_fam
        # field pads absent entries with '-' ('GH18_e1+-'), so each item is cleaned on its own
        split = lambda s: [("" if x.strip() == "-" else x.strip())  # noqa: E731
                           for x in re.split(r"[+,;|]", s)]
        fam_list = [x for x in split(fams) if x] or [""]
        sub_list = split(subs)
        for i, fam in enumerate(fam_list):
            out.append((gene, ec, fam, sub_list[i] if i < len(sub_list) else "", substrate,
                        tools))
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
    ap.add_argument("--out", required=True)
    ap.add_argument("--signalp-gff3", help="SignalP 6 output.gff3 (what BFD's table is built from)")
    ap.add_argument("--signalp", help="SignalP 6 prediction_results.txt (fallback)")
    ap.add_argument("--tmhmm")
    ap.add_argument("--predgpi", help="PredGPI GFF3")
    ap.add_argument("--targetp", help="TargetP 2.0 short summary")
    ap.add_argument("--pfam", help="hmmsearch --domtblout against Pfam-A")
    ap.add_argument("--cazy-overview", help="dbcanlight conclude overview.tsv")
    ap.add_argument("--dbcan", help="hmmsearch --domtblout against dbCAN (no substrates)")
    ap.add_argument("--lengths", help="TSV of protein_id<TAB>length, from the proteome")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    D = lambda f: os.path.join(a.out, f)  # noqa: E731
    lt = a.locustag
    done = []

    if a.signalp_gff3:
        with open(a.signalp_gff3) as fh:
            sp = parse_signalp_gff3(fh)
        n = append(D("signalp.tsv"),
                   ["protein_id", "species_prefix", "peptide_start", "peptide_end", "probability"],
                   [(p, lt, s, e, prob) for p, prob, s, e in sp])
        done.append(f"{n} signalp(gff3)")
    elif a.signalp:
        with open(a.signalp) as fh:
            sp = parse_signalp6(fh)
        n = append(D("signalp.tsv"),
                   ["protein_id", "species_prefix", "peptide_start", "peptide_end", "probability"],
                   [(p, lt, s, e, prob) for p, prob, s, e in sp])
        done.append(f"{n} signalp(txt)")
    if a.tmhmm:
        with open(a.tmhmm) as fh:
            tm = parse_tmhmm_short(fh)
        done.append(f"{append(D('tmhmm.tsv'), ['protein_id', 'len', 'ExpAA', 'First60', 'PredHel', 'Topology'], tm)} tmhmm")
    if a.predgpi:
        with open(a.predgpi) as fh:
            gp = parse_predgpi_gff3(fh)
        done.append(f"{append(D('gpi.tsv'), ['protein_id', 'feature', 'start', 'end'], gp)} gpi")
    if a.targetp:
        with open(a.targetp) as fh:
            tp = parse_targetp(fh)
        done.append(f"{append(D('targetp.tsv'), ['protein_id', 'prediction', 'probability', 'cleavage_position_start', 'cleavage_position_end'], tp)} targetp")
    if a.pfam:
        with open(a.pfam) as fh:
            pf = parse_domtbl(fh)
        done.append(f"{append(D('pfam.tsv'), ['protein_id', 'pfam_acc', 'pfam_id', 'domain_i_evalue', 'ali_from', 'ali_to'], [(t, (acc.split('.')[0] if acc and acc != '-' else q), q, ev, af, at) for t, q, acc, ev, af, at in pf])} pfam")
    cz_header = ["protein_id", "species_prefix", "EC", "cazyme_fam", "sub_fam", "substrate",
                 "toolcount"]
    if a.cazy_overview:
        with open(a.cazy_overview) as fh:
            cz = parse_cazy_overview(fh)
        done.append(f"{append(D('cazy.tsv'), cz_header, [(g, lt, ec, fam, sub, substrate, tools) for g, ec, fam, sub, substrate, tools in cz])} cazy")
    elif a.dbcan:
        with open(a.dbcan) as fh:
            cz = parse_domtbl(fh)
        done.append(f"{append(D('cazy.tsv'), cz_header, [(t, lt, '', fam_from_hmm(q), '', '', 1) for t, q, _a, ev, _f, _o in cz if ev < 1e-5])} cazy(hmm)")
    if a.lengths:
        with open(a.lengths) as fh:
            rows = [l.rstrip("\n").split("\t")[:2] for l in fh if l.strip()]
        done.append(f"{append(D('lengths.tsv'), ['protein_id', 'transcript_id', 'length'], [(r[0], r[0], r[1]) for r in rows if len(r) == 2])} lengths")
    print(f"{lt}: " + ", ".join(done) + f" -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
