#!/usr/bin/env python3
"""Does a wider reference set make a step HMM more sensitive, and where would a cutoff then sit? (leave-one-out evaluation)

Positives are reviewed fungal UniProt entries of the step's EC whose protein name matches --name-pattern (for example the entries
named geranylgeranyl pyrophosphate synthase); other entries of the same EC (for example terpene synthases) are kept as a second
class and are not used for training. Negatives are the step's decoys.

For each positive, an HMM is built from the step's own references plus every OTHER positive (leave-one-out), and the held-out entry is
scored. The lowest held-out score is the cutoff of the owner-approved rule applied to the wider set. The full wider HMM is then
scored against the decoys and the second class. The current HMM (references only) is scored the same way, with the best-domain
bit score and its i-Evalue.

Nothing here is a call and nothing under refs/ changes. The wider HMM is for comparison.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathway_census import hits  # noqa: E402
import sensitivity_audit as sa  # noqa: E402


def search_scores(hmm, fasta_path, cpus=2):
    """name -> (best domain score, i-Evalue of that domain), E <= 1000."""
    with tempfile.TemporaryDirectory() as td:
        dt = Path(td) / "d.domtbl"
        hits.run_tool(["hmmsearch", "--cpu", str(cpus), "--noali", "-o", "/dev/null", "--domtblout", str(dt),
                       "-E", "1000", "--domE", "1000", str(hmm), str(fasta_path)])
        df = hits.parse_domtbl(dt)
    if df.empty:
        return {}
    best = df.sort_values("score", ascending=False).drop_duplicates("protein_id")
    return {r.protein_id: (float(r.score), float(r.evalue)) for r in best.itertuples()}


def split_classes(entries, pattern):
    """entries: list of (acc, header, seq). Return (positives, others) by protein-name regex on the header."""
    rx = re.compile(pattern, re.I)
    pos, oth = [], []
    for e in entries:
        name = e[1].split(" OS=")[0]
        (pos if rx.search(name) else oth).append(e)
    return pos, oth


def write_fasta(path, records):
    with open(path, "w") as fh:
        for acc, _, seq in records:
            fh.write(f">{acc}\n{seq}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", required=True)
    ap.add_argument("--ec", required=True)
    ap.add_argument("--name-pattern", required=True, help="regex on the UniProt protein name that defines the positives")
    ap.add_argument("--refs-dir", required=True, type=Path)
    ap.add_argument("--current-hmm", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--hmm-out-dir", type=Path, help="where to write the full wider HMM (named <step>.PROVISIONAL.hmm)")
    a = ap.parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    own = sa.own_accessions(a.refs_dir, a.step)
    entries = [e for e in sa.read_fasta(sa.fetch_fasta(a.ec)) if e[0] not in own]
    pos, oth = split_classes(entries, a.name_pattern)
    refs = [(n, h, s) for n, h, s in sa.read_fasta((a.refs_dir / f"{a.step}.faa").read_text())]
    decoys = [(n, h, s) for n, h, s in sa.read_fasta((a.refs_dir / "decoys" / f"{a.step}.faa").read_text())]
    work = a.out_dir / "work"
    work.mkdir(exist_ok=True)
    write_fasta(work / "pos.faa", pos)
    write_fasta(work / "oth.faa", oth)
    write_fasta(work / "decoys.faa", decoys)
    cur_pos, cur_oth, cur_dec = (search_scores(a.current_hmm, work / f) for f in ("pos.faa", "oth.faa", "decoys.faa"))
    rows = []
    for acc, hdr, seq in pos:
        train = refs + [p for p in pos if p[0] != acc]
        write_fasta(work / "train.faa", train)
        hmm = work / "loo.hmm"
        hits.build_hmm(work / "train.faa", hmm)
        write_fasta(work / "one.faa", [(acc, hdr, seq)])
        sc = search_scores(hmm, work / "one.faa").get(acc)
        c = cur_pos.get(acc)
        rows.append(dict(accession=acc, name=hdr.split(" OS=")[0].split(" ", 1)[-1][:60], organism=(hdr.split(" OS=")[1].split(" OX=")[0] if " OS=" in hdr else "")[:40],
                         current_score=None if c is None else c[0], current_ievalue=None if c is None else c[1],
                         loo_score=None if sc is None else sc[0], loo_ievalue=None if sc is None else sc[1]))
    write_fasta(work / "wide_train.faa", refs + pos)
    wide = work / "wide.hmm"
    hits.build_hmm(work / "wide_train.faa", wide)
    if a.hmm_out_dir:
        a.hmm_out_dir.mkdir(parents=True, exist_ok=True)
        (a.hmm_out_dir / f"{a.step}.PROVISIONAL.hmm").write_text(wide.read_text())
    wide_dec, wide_oth = search_scores(wide, work / "decoys.faa"), search_scores(wide, work / "oth.faa")
    with open(a.out_dir / "loo_positives.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(a.out_dir / "decoys_and_others.tsv", "w") as fh:
        fh.write("class\taccession\tname\tcurrent_score\twide_score\n")
        for cls, recs, cur, wd in (("decoy", decoys, cur_dec, wide_dec), ("other_same_EC", oth, cur_oth, wide_oth)):
            for acc, hdr, _ in recs:
                fh.write(f"{cls}\t{acc}\t{hdr.split(' OS=')[0].split(' ', 1)[-1][:60]}\t{(cur.get(acc) or ('',))[0]}\t{(wd.get(acc) or ('',))[0]}\n")
    print(f"{len(pos)} positives, {len(oth)} other entries, {len(decoys)} decoys, {len(refs)} references")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
