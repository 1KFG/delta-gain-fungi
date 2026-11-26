#!/usr/bin/env python3
"""Search one chunk of BFD genomes with every step HMM of one pathway (census step 3 and 4).

One job handles many genomes. The protein files of the chunk are joined into one FASTA on node-local
scratch. Each step HMM is searched once against that FASTA. BFD protein ids start with the locustag,
so ids are unique across genomes. The script checks that and counts violations per genome.

Outputs in --out-dir (names carry the chunk id; a pathway with verified: false needs --provisional):
  candidates.chunkNNNNN[.PROVISIONAL].tsv.zst  one row per (protein, step) at or above min_score, with step_id
  searched.chunkNNNNN.tsv                      searched-genome manifest. Written LAST: its presence means the chunk is complete.
  timing.chunkNNNNN.tsv                        seconds per stage
A genome with a missing or unreadable protein file gets status `missing` and is never an absence.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402


def read_chunk(manifest: Path, chunk_id: int):
    if str(manifest).endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", str(manifest)], check=True, capture_output=True).stdout.decode()
    else:
        raw = Path(manifest).read_text()
    rows = [r for r in csv.DictReader(io.StringIO(raw), delimiter="\t") if int(r["chunk_id"]) == chunk_id]
    if not rows:
        raise SystemExit(f"chunk {chunk_id} not in {manifest}")
    return rows


def join_proteins(rows, out_fa: Path):
    """Write all protein files of the chunk to out_fa. Return one manifest row per genome."""
    man = []
    with open(out_fa, "w") as out:
        for r in rows:
            lt, path = r["locustag"], r["proteins_path"]
            n = bad = 0
            status = "searched"
            try:
                opener = gzip.open if path.endswith(".gz") else open
                with opener(path, "rt") as fh:
                    for line in fh:
                        if line.startswith(">"):
                            n += 1
                            if not line.startswith(">" + lt + "_"):
                                bad += 1
                        out.write(line)
                    if n == 0:
                        status = "missing"
            except (OSError, EOFError) as e:  # missing, unreadable or truncated gzip
                status, n, bad = f"missing:{type(e).__name__}", 0, 0
            man.append(dict(locustag=lt, proteins_path=path, status=status, n_proteins=n, n_bad_id_prefix=bad))
    return man


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", required=True)
    ap.add_argument("--hmm-dir", required=True, type=Path)
    ap.add_argument("--chunk-manifest", required=True, type=Path)
    ap.add_argument("--chunk-id", required=True, type=int)
    ap.add_argument("--min-scores", required=True, type=Path, help="TSV step, min_score (bin/make_min_scores.py)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args(argv)

    pw = load_pathway(Path(a.pathway))
    tag = hits.provisional_tag(pw, a.provisional)
    scratch_env = os.environ.get("SCRATCH")
    if not scratch_env:
        raise SystemExit("SCRATCH is not set")
    ms_rows = list(csv.DictReader(open(a.min_scores), delimiter="\t"))
    no_floor = sorted(r["step"] for r in ms_rows if r["min_score"] == "")
    if no_floor:  # make_min_scores.py --allow-undefined wrote no value; a search must not guess one
        raise SystemExit(f"min_score is empty (no floor defined) for steps: {no_floor}")
    ms = {r["step"]: float(r["min_score"]) for r in ms_rows}
    steps = [s for s in ms if (a.hmm_dir / hits.tagged_name(s, "hmm", tag)).exists()]
    missing_hmm = sorted(set(ms) - set(steps))
    if missing_hmm:  # fail before any staging work
        raise SystemExit(f"no HMM for steps: {missing_hmm}")
    work = Path(scratch_env) / f"chunk{a.chunk_id:05d}"
    work.mkdir(parents=True, exist_ok=True)
    cid = f"chunk{a.chunk_id:05d}"
    timing = []

    t0 = time.time()
    rows = read_chunk(a.chunk_manifest, a.chunk_id)
    fa = work / "proteins.faa"
    man = join_proteins(rows, fa)
    timing.append(("stage_and_count_proteins", time.time() - t0))
    n_prot = sum(m["n_proteins"] for m in man)

    parts = []
    for step in steps:
        t1 = time.time()
        dt = work / f"{step}.domtbl"
        hits.search_step(a.hmm_dir / hits.tagged_name(step, "hmm", tag), fa, dt, a.cpus)
        cand = hits.candidates_for_step(hits.parse_domtbl(dt), ms[step])
        cand.insert(len(cand.columns), "step_id", step)
        parts.append(cand)
        dt.unlink()
        timing.append((f"hmmsearch:{step}", time.time() - t1))
    cands = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    a.out_dir.mkdir(parents=True, exist_ok=True)
    sfx = f".{tag}" if tag else ""
    cpath = a.out_dir / f"candidates.{cid}{'.' + tag if tag else ''}.tsv.zst"
    buf = io.StringIO()
    if tag:
        buf.write(f"# {tag}\n")
    cands.to_csv(buf, sep="\t", index=False)
    subprocess.run(["zstd", "-q", "-T0", "-f", "-o", str(cpath)], input=buf.getvalue().encode(), check=True)
    timing.append(("total", time.time() - t0))
    with open(a.out_dir / f"timing.{cid}.tsv", "w") as fh:
        fh.write("stage\tseconds\n")
        for k, v in timing:
            fh.write(f"{k}\t{v:.1f}\n")
        fh.write(f"n_genomes\t{len(man)}\nn_proteins\t{n_prot}\nn_candidates\t{len(cands)}\ncpus\t{a.cpus}\n")
    # manifest last
    with open(a.out_dir / f"searched.{cid}.tsv", "w") as fh:
        fh.write(f"# chunk={a.chunk_id} tag={tag or 'none'} min_scores={a.min_scores} hmm_dir={a.hmm_dir}\n")
        w = csv.DictWriter(fh, fieldnames=["locustag", "proteins_path", "status", "n_proteins", "n_bad_id_prefix"], delimiter="\t")
        w.writeheader()
        w.writerows(man)
    fa.unlink()
    print(f"{cid}: genomes {len(man)}, proteins {n_prot}, candidates {len(cands)}, total {timing[-1][1]:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
