#!/usr/bin/env python3
"""Step 2 of the census: score references and decoys against the step HMMs.

For each step HMM in --hmm-dir the script runs hmmsearch against one target set that holds
every reference and every decoy of the pathway. It reports, per step:
  - the scores of the step's own references (the PROPOSED min_score is the lowest of them);
  - every other target (decoy, or reference of another step) that scores at or above that value;
  - a leave-one-out score for each own reference when the step has two or more references
    (the HMM is rebuilt without that reference, then the held-out reference is scored).
A step with one reference has no leave-one-out result. That is reported, not hidden.

Nothing here is a biological finding. Scores are bit scores of the best domain per target.
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402


def read_fasta(path: Path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def best_per_target(domtbl: pd.DataFrame) -> dict[str, float]:
    if domtbl.empty:
        return {}
    return domtbl.groupby("protein_id").score.max().to_dict()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs-dir", required=True, type=Path)
    ap.add_argument("--hmm-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--tag", default="", help="for example PROVISIONAL; added to output names")
    ap.add_argument("--cpus", type=int, default=4)
    a = ap.parse_args()

    suffix = f".{a.tag}" if a.tag else ""
    steps = sorted(p.stem for p in a.refs_dir.glob("*.faa"))
    # target set: ids are kind|step|index|accession so that no two records collide
    targets, meta = [], {}
    for step in steps:
        for kind, d in (("ref", a.refs_dir), ("decoy", a.refs_dir / "decoys")):
            f = d / f"{step}.faa"
            if not f.exists():
                continue
            for i, (hdr, seq) in enumerate(read_fasta(f)):
                tid = f"{kind}|{step}|{i}|{hdr.split()[0]}"
                targets.append((tid, seq))
                meta[tid] = (kind, step, hdr.split()[0])
    a.out_dir.mkdir(parents=True, exist_ok=True)

    scratch = tempfile.mkdtemp(dir=None)
    tfa = Path(scratch) / "targets.faa"
    with open(tfa, "w") as fh:
        for tid, seq in targets:
            fh.write(f">{tid}\n{seq}\n")

    summary, above, loo_rows = [], [], []
    for step in steps:
        hmm = a.hmm_dir / f"{step}{suffix}.hmm"
        if not hmm.exists():
            summary.append(dict(step=step, note="no HMM"))
            continue
        dt = Path(scratch) / f"{step}.domtbl"
        hits.search_step(hmm, tfa, dt, a.cpus)
        sc = best_per_target(hits.parse_domtbl(dt))
        own = {t: sc.get(t) for t, (k, s, _) in meta.items() if k == "ref" and s == step}
        own_scores = [v for v in own.values() if v is not None]
        n_refs = len(own)
        missing = sum(v is None for v in own.values())
        min_score = min(own_scores) if own_scores and missing == 0 else None
        n_dec = n_oth = 0
        max_dec = None
        for t, v in sc.items():
            k, s, acc = meta[t]
            if k == "ref" and s == step:
                continue
            if k == "decoy" and (max_dec is None or v > max_dec):
                max_dec = v
            if min_score is not None and v >= min_score:
                above.append(dict(step=step, target_kind=k, target_step=s, accession=acc, score=v))
                n_dec += k == "decoy"
                n_oth += k == "ref"
        # leave-one-out
        loo_scores = []
        if n_refs >= 2:
            recs = list(read_fasta(a.refs_dir / f"{step}.faa"))
            for i, (hdr, seq) in enumerate(recs):
                rest = Path(scratch) / f"{step}.loo{i}.faa"
                with open(rest, "w") as fh:
                    for j, (h2, s2) in enumerate(recs):
                        if j != i:
                            fh.write(f">{h2.split()[0]}_{j}\n{s2}\n")
                lh = Path(scratch) / f"{step}.loo{i}.hmm"
                hits.build_hmm(rest, lh)
                ldt = Path(scratch) / f"{step}.loo{i}.domtbl"
                hits.search_step(lh, tfa, ldt, a.cpus)
                lsc = best_per_target(hits.parse_domtbl(ldt))
                tid = f"ref|{step}|{i}|{hdr.split()[0]}"
                loo_scores.append(lsc.get(tid))
                loo_rows.append(dict(step=step, held_out=hdr.split()[0], loo_score=lsc.get(tid),
                                     full_hmm_score=own.get(tid)))
        loo_found = [x for x in loo_scores if x is not None]
        summary.append(dict(
            step=step, n_refs=n_refs, refs_not_found_by_own_hmm=missing,
            ref_min=min(own_scores) if own_scores else None,
            ref_max=max(own_scores) if own_scores else None,
            proposed_min_score=min_score,
            n_decoys_at_or_above=n_dec, max_decoy_score=max_dec,
            n_other_step_refs_at_or_above=n_oth,
            loo_n=len(loo_scores), loo_not_found=len(loo_scores) - len(loo_found),
            loo_min=min(loo_found) if loo_found else None))

    pd.DataFrame(summary).to_csv(a.out_dir / f"step_scores{suffix}.tsv", sep="\t", index=False)
    pd.DataFrame(above, columns=["step", "target_kind", "target_step", "accession", "score"]) \
        .to_csv(a.out_dir / f"targets_at_or_above_min{suffix}.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows, columns=["step", "held_out", "loo_score", "full_hmm_score"]) \
        .to_csv(a.out_dir / f"leave_one_out{suffix}.tsv", sep="\t", index=False)
    print(f"wrote {a.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
