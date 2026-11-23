#!/usr/bin/env python3
"""Build the input files for a tree-placement TEST of ggpp_synthase at a lower cutoff (results folder only; nothing under refs/).

From UniProt (reviewed, Fungi, taxonomy_id 4751; fetched at run time and recorded in provenance.tsv):
  - named members: entries of EC 2.5.1.29 whose protein name matches 'geranylgeranyl ... synthase' (tree leaves only), minus the references,
    minus any accession in --exclude (fission-yeast O59703: it scores below its own FPPS under the wider HMM);
  - bifunctional entries: the other entries of EC 2.5.1.29 (terpene synthases with a prenyltransferase domain), kept as a labelled
    second class of tree leaves, not used to build anything;
  - decoys: the step's existing decoys plus reviewed fungal farnesyl diphosphate synthases (EC 2.5.1.10) and hexaprenyl diphosphate
    synthases (EC 2.5.1.83 and 2.5.1.82), the families closest to GGPP synthase with another function.
Each non-target entry is scored against the wider HMM and the current HMM (best domain, E <= 1000) and its coverage of the HMM is
recorded, so alignability can be judged before it goes into a tree.
Nothing here is a call and nothing is evidence for a step call.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sensitivity_audit as sa  # noqa: E402
import expanded_hmm_cv as cv  # noqa: E402
import filter_decoys_for_tree as fd  # noqa: E402


def classify(entries, pattern, refs, exclude):
    """entries: (acc, header, seq). Return (named, bifunctional, excluded accessions seen)."""
    pos, oth = cv.split_classes([e for e in entries if e[0] not in refs], pattern)
    named = [e for e in pos if e[0] not in exclude]
    return named, oth, sorted({e[0] for e in pos} & set(exclude))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs-dir", required=True, type=Path)
    ap.add_argument("--current-hmm", required=True, type=Path)
    ap.add_argument("--wide-hmm", required=True, type=Path)
    ap.add_argument("--exclude", default="O59703")
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    step = "ggpp_synthase"
    own = sa.own_accessions(a.refs_dir, step)
    refs = {n for n, _, _ in sa.read_fasta((a.refs_dir / f"{step}.faa").read_text())}
    ggpps = list(sa.read_fasta(sa.fetch_fasta("2.5.1.29")))
    named, bifunc, excluded = classify(ggpps, "geranylgeranyl.*synthase", refs, set(a.exclude.split(",")))
    decoys = {n: (n, h, s) for n, h, s in sa.read_fasta((a.refs_dir / "decoys" / f"{step}.faa").read_text())}
    decoy_source = {n: "existing decoy (refs/carotenoid/decoys)" for n in decoys}
    for ec in ("2.5.1.10", "2.5.1.83", "2.5.1.82"):
        for n, h, s in sa.read_fasta(sa.fetch_fasta(ec)):
            if n not in decoys and n not in refs and n not in {e[0] for e in ggpps}:
                decoys[n] = (n, h, s)
                decoy_source[n] = f"UniProt reviewed fungal EC {ec}"
    for name, recs in (("named_members", named), ("bifunctional", bifunc), ("decoys", list(decoys.values()))):
        cv.write_fasta(a.out_dir / f"{name}.faa", recs)
    scored = {}
    for hmm_name, hmm in (("wide", a.wide_hmm), ("current", a.current_hmm)):
        for name in ("named_members", "bifunctional", "decoys"):
            sc = cv.search_scores(hmm, a.out_dir / f"{name}.faa")
            scored.setdefault(name, {})[hmm_name] = sc
    wide_len, cur_len = fd.hmm_length(a.wide_hmm), fd.hmm_length(a.current_hmm)
    cover = {}
    for hmm_name, hmm, ln in (("wide", a.wide_hmm, wide_len), ("current", a.current_hmm, cur_len)):
        for name in ("decoys", "named_members", "bifunctional"):
            bs = fd.best_scores(hmm, a.out_dir / f"{name}.faa")
            cover[(hmm_name, name)] = {k: (v[2] - v[1] + 1) / ln for k, v in bs.items()}
    with open(a.out_dir / "provenance.tsv", "w") as fh:
        fh.write("class\taccession\tname\torganism\tsource\twide_score\twide_coverage\tcurrent_score\tcurrent_coverage\n")
        for name, recs in (("named_member", named), ("bifunctional", bifunc), ("decoy", list(decoys.values()))):
            key = {"named_member": "named_members", "bifunctional": "bifunctional", "decoy": "decoys"}[name]
            for acc, hdr, _ in recs:
                w = scored[key]["wide"].get(acc, ("",))[0]
                c = scored[key]["current"].get(acc, ("",))[0]
                src = decoy_source.get(acc, "UniProt reviewed fungal EC 2.5.1.29")
                wc = cover[("wide", key)].get(acc, "")
                cc = cover[("current", key)].get(acc, "")
                wc = f"{wc:.3f}" if wc != "" else ""
                cc = f"{cc:.3f}" if cc != "" else ""
                fh.write(f"{name}\t{acc}\t{hdr.split(' OS=')[0].split(' ', 1)[-1][:70]}\t{(hdr.split(' OS=')[1].split(' OX=')[0] if ' OS=' in hdr else '')[:40]}\t{src}\t{w}\t{wc}\t{c}\t{cc}\n")
    print(f"named {len(named)}, bifunctional {len(bifunc)}, decoys {len(decoys)}, excluded {excluded}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
