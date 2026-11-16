#!/usr/bin/env python3
"""Build AlphaFold 3 input JSONs: single chains, protein pairs, and chito-oligosaccharides.

  af3_input.py single  --fasta in.faa --out jobs/
  af3_input.py pairs   --fasta all.faa --a fungal.txt --b host.txt --msa-dir msas/ --out jobs/
  af3_input.py ligand  --fasta all.faa --ids ids.txt --msa-dir msas/ --nag 4,7 --out jobs/

The screen is built in two stages because MSAs are the expensive part and are a property of a
chain, not of a pair: `single` writes one job per chain for AF3's data pipeline
(--norun_inference), and `pairs`/`ligand` then read those MSAs back out and embed them, so
inference runs with --norun_data_pipeline and every pair reuses work already done. With ~19
chains and ~90 pairs that is the difference between 19 and 90 jackhmmer runs.

Chito-oligosaccharides are built from CCD NAG components joined by explicit beta-1,4 bonds
(C1 of residue i to O4 of residue i+1). AF3 validates atom names, so a wrong linkage fails
loudly rather than silently modelling disconnected sugars. CO4 and CO7 are the lengths that
matter: in Marchantia paleacea the LysM pair MpaLYR/MpaCERK1 discriminates symbiosis from
immunity by oligomer length, CO4 activating symbiosis genes and CO7 triggering immunity.
Standard library only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

NAG = "NAG"


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def clean(seq):
    return seq.upper().rstrip("*").replace("*", "").replace("-", "")


def safe(name):
    """AF3 names the output directory after the job, so keep it filesystem-safe."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def protein_entry(chain_id, seq, msa=None):
    """One protein chain. `msa` is {unpairedMsa, pairedMsa, templates} from a data-pipeline run.

    Both MSA fields must be set together, and templates must be present (even empty), or AF3
    rejects the input: a chain is either fully precomputed or fully left to the pipeline."""
    e = {"id": chain_id, "sequence": seq}
    if msa is not None:
        e["unpairedMsa"] = msa.get("unpairedMsa", "")
        e["pairedMsa"] = msa.get("pairedMsa", "")
        e["templates"] = msa.get("templates", [])
    return {"protein": e}


def nag_ligand(chain_id, n):
    """-> (ligand entry, bonds) for a beta-1,4 linked (GlcNAc)n."""
    entry = {"ligand": {"id": chain_id, "ccdCodes": [NAG] * n}}
    bonds = [[[chain_id, i, "C1"], [chain_id, i + 1, "O4"]] for i in range(1, n)]
    return entry, bonds


def job_json(name, sequences, bonds=None, seed=1):
    job = {"name": name, "modelSeeds": [seed], "sequences": sequences,
           "dialect": "alphafold3", "version": 2}
    if bonds:
        job["bondedAtomPairs"] = bonds
    return job


def extract_msas(path):
    """Read an AF3 data-pipeline output JSON -> {chain_id: {unpairedMsa, pairedMsa, templates}}."""
    with open(path) as fh:
        data = json.load(fh)
    out = {}
    for s in data.get("sequences", []):
        p = s.get("protein")
        if not p:
            continue
        out[p["id"]] = {"unpairedMsa": p.get("unpairedMsa", "") or "",
                        "pairedMsa": p.get("pairedMsa", "") or "",
                        "templates": p.get("templates", []) or [],
                        # kept so the caller can check the MSA belongs to the sequence it is
                        # about to be attached to; stripped before the job JSON is written
                        "_query_seq": p.get("sequence", "") or ""}
    return out


def find_msa(msa_dir, name):
    """AF3 writes <out>/<lowercased job name>/<name>_data.json; accept either spelling."""
    n = safe(name)
    for cand in (os.path.join(msa_dir, n.lower(), f"{n.lower()}_data.json"),
                 os.path.join(msa_dir, n, f"{n}_data.json"),
                 os.path.join(msa_dir, n.lower(), f"{n}_data.json")):
        if os.path.exists(cand):
            return cand
    return None


def msa_query_seq(msa):
    """The sequence the MSA was actually built for, or None."""
    return (msa or {}).get("_query_seq")


def msa_for(msa_dir, name, seq=None):
    """-> the single-chain MSA block for `name`, or None if missing or built for another chain.

    The `seq` check is not pedantry. MSAs are looked up by chain name, and the control chains are
    named CTRL_<pdb>_<n> by their position in an RCSB fetch that af3_prepare.sh repeats on every
    run. Re-running prepare can therefore move a name onto a different chain while the stage-1
    MSAs under msas/ still belong to the old one. The job JSON then carries a sequence of one
    length with an alignment of another, and AF3 fails in featurising with nothing that points at
    the cause - which is what happened to all 26 control predictions in job 32069698 (A was 68 aa
    carrying the 75 aa chain's MSA, and B the reverse). Proteins with real locus tags are immune,
    which is why the other 170 jobs were fine.
    """
    p = find_msa(msa_dir, name)
    if not p:
        return None
    m = extract_msas(p)
    msa = m.get("A") or (list(m.values())[0] if m else None)
    if msa is None:
        return None
    if seq is not None:
        q = msa.get("_query_seq") or ""
        if q and q != seq:
            print(f"MSA for {name} was built for a {len(q)} aa sequence but the panel has "
                  f"{len(seq)} aa - stale MSA, refusing to pair them", file=sys.stderr)
            return None
    return {k: v for k, v in msa.items() if not k.startswith("_")}


def write(job, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    p = os.path.join(out_dir, f"{safe(job['name'])}.json")
    with open(p, "w") as fh:
        json.dump(job, fh, indent=2)
    return p


def read_ids(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for cmd in ("single", "pairs", "ligand"):
        s = sub.add_parser(cmd)
        s.add_argument("--fasta", required=True)
        s.add_argument("--out", required=True)
        s.add_argument("--seed", type=int, default=1)
        if cmd == "single":
            s.add_argument("--ids", help="subset of the FASTA to write")
            # Optional, and the two uses are opposites. Stage 1 (slurm/af3_msa.sh) calls `single`
            # WITHOUT it, because the whole point of that run is to let the data pipeline build
            # the MSA. Stage 2 (slurm/af3_fold.sh) runs with --norun_data_pipeline and must pass
            # it, or every job is a bare sequence with no alignment and AF3 dies in featurising -
            # which is exactly what happened to all 426 unknowns in job 32069699.
            s.add_argument("--msa-dir", help="attach the MSA each chain's stage-1 run produced; "
                                             "required when folding with --norun_data_pipeline")
        else:
            s.add_argument("--msa-dir", required=True)
        if cmd == "pairs":
            s.add_argument("--a", required=True)
            s.add_argument("--b", required=True)
        if cmd == "ligand":
            s.add_argument("--ids", required=True)
            s.add_argument("--nag", default="4,7")
    a = ap.parse_args(argv)

    seqs = {n: clean(s) for n, s in read_fasta(a.fasta)}
    n_written, missing = 0, []

    if a.cmd == "single":
        want = read_ids(a.ids) if a.ids else list(seqs)
        for name in want:
            if name not in seqs:
                missing.append(name)
                continue
            m = None
            if a.msa_dir:
                m = msa_for(a.msa_dir, name, seqs[name])
                if m is None:
                    # Writing it anyway would produce a job that fails in featurising rather
                    # than one that is obviously absent; report it instead.
                    missing.append(name)
                    continue
            write(job_json(name, [protein_entry("A", seqs[name], m)], seed=a.seed), a.out)
            n_written += 1
    elif a.cmd == "pairs":
        for x in read_ids(a.a):
            mx = msa_for(a.msa_dir, x, seqs.get(x))
            if x not in seqs or mx is None:
                missing.append(x)
                continue
            for y in read_ids(a.b):
                my = msa_for(a.msa_dir, y, seqs.get(y))
                if y not in seqs or my is None:
                    missing.append(y)
                    continue
                name = f"{x}__{y}"
                job = job_json(name, [protein_entry("A", seqs[x], mx),
                                      protein_entry("B", seqs[y], my)], seed=a.seed)
                write(job, a.out)
                n_written += 1
    else:
        lengths = [int(x) for x in a.nag.split(",") if x]
        for name in read_ids(a.ids):
            m = msa_for(a.msa_dir, name, seqs.get(name))
            if name not in seqs or m is None:
                missing.append(name)
                continue
            for n in lengths:
                lig, bonds = nag_ligand("B", n)
                job = job_json(f"{name}__CO{n}",
                               [protein_entry("A", seqs[name], m), lig], bonds, seed=a.seed)
                write(job, a.out)
                n_written += 1

    if missing:
        print(f"missing sequence or MSA for {len(set(missing))}: "
              f"{sorted(set(missing))[:8]}", file=sys.stderr)
    print(f"wrote {n_written} AF3 jobs -> {a.out}", file=sys.stderr)
    return 0 if n_written else 1


if __name__ == "__main__":
    sys.exit(main())
