#!/usr/bin/env python3
"""Candidate effectors from the annotated panels, and the set to fold.

  call_effectors.py --dir out/panels --locustags JEM1MUCOR,JEM2MUCOR --out out/effectors \
      [--codons out/panels/jem_codon_check.tsv]

Reads the same tables build_secretome.py reads (signalp, tmhmm, gpi, targetp, pfam, cazy,
lengths) and reuses its secretion call and tiers, so a candidate here means what it means for
the 294 BFD Mucoromycota genomes. What this adds is the fold-first selection: which secreted
proteins are worth a predicted structure because their sequence says nothing.

A protein enters the fold set if it is secreted or GPI-anchored and either
  * carries no Pfam domain and is at most --max-len aa (`unknown_small`), whether or not it is
    cysteine-rich - a fold is exactly what is missing for these, and the T3 cysteine thresholds
    would drop unknown proteins that happen to be cysteine-poor; or
  * is in tier T1, the chitin/chitosan-active set (CE4, GH18/20/46/75, LysM/CBM), where the
    structure says whether the chitin-binding site is intact.
Everything else is reported in the table but not folded.

`seq_start_ATG` is carried through when --codons is given: a model with no real ATG start may be
missing its true 5' end, so its signal peptide and its fold are both suspect. Such proteins are
flagged, not dropped, because for these genomes the flag itself has been unreliable.
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_secretome as bs  # noqa: E402
from extract_proteins import read_fasta  # noqa: E402

FOLD_MAX_LEN = 300


def load_tables(d):
    """-> (signalp, tmhmm, gpi, targetp, pfam, cazy, substrates, lengths) keyed by protein_id."""
    D = lambda f: os.path.join(d, f)  # noqa: E731
    sig = {r["protein_id"]: r for r in bs.tsv(D("signalp.tsv"))}
    tm = {r["protein_id"]: r for r in bs.tsv(D("tmhmm.tsv"), optional=True)}
    gpi = bs.gpi_set(bs.tsv(D("gpi.tsv"), optional=True))
    tp = {r["protein_id"]: r for r in bs.tsv(D("targetp.tsv"), optional=True)}
    pfam, cazy, subs = defaultdict(list), defaultdict(list), defaultdict(list)
    for r in bs.tsv(D("pfam.tsv"), optional=True):
        pfam[r["protein_id"]].append((r["pfam_acc"], r["pfam_id"]))
    for r in bs.tsv(D("cazy.tsv"), optional=True):
        if r.get("cazyme_fam"):
            cazy[r["protein_id"]].append(r["cazyme_fam"])
        if r.get("substrate"):
            subs[r["protein_id"]].append(r["substrate"])
    length = {r["protein_id"]: int(r["length"])
              for r in bs.tsv(D("lengths.tsv"), optional=True) if r.get("length")}
    return sig, tm, gpi, tp, pfam, cazy, subs, length


def mito_prob(tp_row):
    """TargetP probability, but only when the predicted class is mTP (see annotate_to_tables)."""
    if (tp_row.get("prediction") or "").lower().startswith("mt"):
        return bs.fnum(tp_row.get("probability"))
    return 0.0


def fold_reason(category, tier, pfam_accs, length, max_len):
    """-> why this protein is folded, or '' if it is not."""
    if category not in ("secreted", "gpi_anchored"):
        return ""
    if not pfam_accs and length and length <= max_len:
        return "unknown_small"
    if tier == "T1":
        return "chitin_active"
    return ""


def call(tables, seqs, locustags, codons=None, max_len=FOLD_MAX_LEN):
    """-> list of row dicts, one per signal-peptide-positive protein of the given locustags."""
    sig, tm, gpi, tp, pfam, cazy, subs, length = tables
    want = set(locustags)
    out = []
    for p, s in sig.items():
        lt = p.split("_")[0]
        if lt not in want:
            continue
        t = tm.get(p, {})
        cleave = int(bs.fnum(s.get("peptide_end")))
        cat = bs.classify(bs.fnum(s.get("probability")), t.get("Topology"), t.get("PredHel"),
                          cleave, p in gpi, mito_prob(tp.get(p, {})))
        seq = seqs.get(p, "")
        n_cys, frac, slen = bs.cys_stats(seq) if seq else (0, 0.0, length.get(p, 0))
        slen = slen or length.get(p, 0)
        accs = [a for a, _ in pfam.get(p, [])]
        tier, why = bs.tier_of(cazy.get(p, []), subs.get(p, []), accs, slen, n_cys, frac)
        row = {
            "protein_id": p, "LOCUSTAG": lt, "category": cat, "tier": tier, "tier_why": why,
            "length": slen, "sigp_prob": s.get("probability"), "sigp_cleavage": cleave,
            "n_cys": n_cys, "cys_frac": round(frac, 4),
            "pfam": ";".join(sorted({a for a in accs})),
            "cazy": ";".join(sorted({bs.cazy_root(f) for f in cazy.get(p, [])})),
            "cazy_substrate": ";".join(sorted(set(subs.get(p, [])))),
            "targetp": (tp.get(p, {}) or {}).get("prediction", ""),
            "seq_start_ATG": (codons or {}).get(p, ""),
        }
        row["fold_reason"] = fold_reason(cat, tier, accs, slen, max_len)
        out.append(row)
    return sorted(out, key=lambda r: (r["LOCUSTAG"], r["protein_id"]))


def from_secretome(path, seqs, max_len=FOLD_MAX_LEN):
    """Rows from build_secretome.py's focal_secretome.tsv, in the same shape as call().

    The BFD reference genome was classified by the earlier Mucoromycota run, whose SQL export no
    longer exists; its per-protein table carries the same category, tier, pfam and length, so the
    fold set can be selected from it without re-exporting."""
    rows = []
    for r in bs.tsv(path):
        accs = [a for a in (r.get("pfam") or "").replace(",", ";").split(";") if a]
        slen = int(bs.fnum(r.get("length")))
        row = {
            "protein_id": r["protein_id"], "LOCUSTAG": r["protein_id"].split("_")[0],
            "category": r.get("category", ""), "tier": r.get("tier", ""),
            "tier_why": r.get("tier_evidence", ""), "length": slen,
            "sigp_prob": r.get("sigp_prob", ""), "sigp_cleavage": r.get("sigp_cleavage", ""),
            "n_cys": r.get("n_cys", ""), "cys_frac": r.get("cys_pct", ""),
            "pfam": r.get("pfam", ""), "cazy": r.get("cazy", ""),
            "cazy_substrate": r.get("cazy_substrate", ""), "targetp": "",
            "seq_start_ATG": r.get("seq_start_ATG", ""),
        }
        row["fold_reason"] = fold_reason(row["category"], row["tier"], accs, slen, max_len)
        # a protein with no sequence available cannot be folded
        if row["fold_reason"] and row["protein_id"] not in seqs:
            row["fold_reason"] = ""
        rows.append(row)
    return sorted(rows, key=lambda r: (r["LOCUSTAG"], r["protein_id"]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", help="directory of the annotation tables")
    ap.add_argument("--secretome", help="build_secretome.py focal_secretome.tsv, instead of --dir")
    ap.add_argument("--fasta", help="proteome FASTA(s), comma separated", default="")
    ap.add_argument("--locustags", default="")
    ap.add_argument("--codons", help="check_codons.py output, for seq_start_ATG")
    ap.add_argument("--max-len", type=int, default=FOLD_MAX_LEN)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    lts = [x for x in a.locustags.split(",") if x]

    seqs = {}
    for f in [x for x in a.fasta.split(",") if x]:
        seqs.update(dict(read_fasta(f)))
    codons = {}
    if a.codons and os.path.exists(a.codons):
        for r in bs.tsv(a.codons):
            key = r.get("transcript_id") or r.get("protein_id")
            if key:
                codons[key] = r.get("seq_start_ATG", r.get("start_codon", ""))

    if a.secretome:
        rows = from_secretome(a.secretome, seqs, a.max_len)
    elif a.dir:
        rows = call(load_tables(a.dir), seqs, lts, codons, a.max_len)
    else:
        ap.error("one of --dir or --secretome is required")
    cols = ["protein_id", "LOCUSTAG", "category", "tier", "tier_why", "fold_reason", "length",
            "n_cys", "cys_frac", "sigp_prob", "sigp_cleavage", "pfam", "cazy", "cazy_substrate",
            "targetp", "seq_start_ATG"]
    with open(os.path.join(a.out, "candidates.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    fold = [r for r in rows if r["fold_reason"]]
    with open(os.path.join(a.out, "to_fold.faa"), "w") as fh:
        for r in fold:
            s = seqs.get(r["protein_id"], "").rstrip("*")
            if s:
                fh.write(f">{r['protein_id']}\n{s}\n")
    for name, c in (("category", Counter(r["category"] for r in rows)),
                    ("tier", Counter(r["tier"] for r in rows)),
                    ("fold_reason", Counter(r["fold_reason"] for r in fold)),
                    ("per genome", Counter(r["LOCUSTAG"] for r in fold))):
        print(f"{name}: " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())), file=sys.stderr)
    print(f"{len(rows)} candidates, {len(fold)} to fold -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
