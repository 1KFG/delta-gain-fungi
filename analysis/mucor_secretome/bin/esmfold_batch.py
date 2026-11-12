#!/usr/bin/env python3
"""ESMFold structures for a FASTA of candidate effectors.

  esmfold_batch.py --fasta to_fold.faa --out folds/ [--max-len 400] [--chunk 128]

ESMFold is single-sequence: no MSA, so nothing is sent anywhere and an orphan effector with no
detectable homologs - which is most of this set - folds just as well as a conserved one. That is
the whole reason it is the right tool here; AlphaFold and Boltz would spend their accuracy on an
MSA these proteins do not have.

Writes <out>/<protein_id>.pdb and appends a row per protein to <out>/plddt.tsv. Proteins whose
PDB already exists are skipped, so the job is resumable after a timeout or preemption. On CUDA
OOM the chunk size is halved and the protein retried, then skipped if it still fails.

Confidence: pLDDT below ~70 means the fold is not reliable, and for these small cysteine-rich
proteins a low pLDDT often means a genuinely disordered or ill-predicted protein rather than a
fixable run - such models are written but flagged, and bin/fold_report.py excludes them from
fold assignment by default.
"""
from __future__ import annotations

import argparse
import os
import sys
import time


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
    """ESMFold takes the 20 standard residues; anything else becomes X."""
    ok = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join(c if c in ok else "X" for c in seq.upper().rstrip("*"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=400)
    ap.add_argument("--chunk", type=int, default=128)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    import torch
    import esm

    # The GPU must be one this torch was built for, or every result is silently wrong. The
    # esm/2.0.0 module ships torch 1.13.1+cu117, built for sm_37..sm_86, so it runs correctly on
    # an a100 (sm_80) but not on the ada6000 (sm_89), h100 (sm_90) or blackwell6000 (sm_120)
    # nodes. On a Blackwell node it does not refuse: it prints a warning and then returns
    # garbage, which in job 32033882 meant 487 nonsense errors and - worse - 72 "successful"
    # models that were also nonsense. Fail loudly instead.
    if not torch.cuda.is_available():
        sys.exit("no CUDA device visible")
    cap = torch.cuda.get_device_capability()
    sm, arch = f"sm_{cap[0]}{cap[1]}", torch.cuda.get_arch_list()
    name = torch.cuda.get_device_name(0)
    print(f"device {name} ({sm}); torch {torch.__version__} built for {arch}", file=sys.stderr)
    if arch and sm not in arch:
        sys.exit(f"REFUSING TO RUN: {name} is {sm}, but torch {torch.__version__} was built for "
                 f"{arch}. Results would be silently wrong - request a supported GPU "
                 f"(--gres=gpu:a100:1).")

    model = esm.pretrained.esmfold_v1().eval().cuda()
    model.set_chunk_size(a.chunk)
    log = os.path.join(a.out, "plddt.tsv")
    new = not os.path.exists(log)
    fh = open(log, "a")
    if new:
        fh.write("protein_id\tlength\tplddt\tptm\tseconds\n")

    n_done = n_skip = n_fail = 0
    for name, raw in read_fasta(a.fasta):
        pdb_path = os.path.join(a.out, f"{name}.pdb")
        if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 0:
            n_skip += 1
            continue
        seq = clean(raw)
        if not seq or len(seq) > a.max_len:
            print(f"skip {name}: length {len(seq)}", file=sys.stderr)
            n_skip += 1
            continue
        chunk = a.chunk
        for attempt in range(3):
            try:
                t0 = time.time()
                with torch.no_grad():
                    out = model.infer(seq)
                    pdb = model.output_to_pdb(out)[0]
                dt = time.time() - t0
                plddt = float(out["mean_plddt"][0])
                ptm = float(out["ptm"][0])
                with open(pdb_path, "w") as pf:
                    pf.write(pdb)
                fh.write(f"{name}\t{len(seq)}\t{plddt:.2f}\t{ptm:.4f}\t{dt:.1f}\n")
                fh.flush()
                n_done += 1
                break
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                chunk = max(8, chunk // 2)
                model.set_chunk_size(chunk)
                print(f"OOM on {name} ({len(seq)} aa), chunk -> {chunk}", file=sys.stderr)
            except Exception as e:  # noqa: BLE001 - one bad protein must not end the run
                print(f"FAIL {name}: {type(e).__name__}: {e}", file=sys.stderr)
                n_fail += 1
                break
        else:
            print(f"FAIL {name}: still OOM at chunk {chunk}", file=sys.stderr)
            n_fail += 1
        model.set_chunk_size(a.chunk)
    fh.close()
    print(f"folded {n_done}, skipped {n_skip}, failed {n_fail} -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
