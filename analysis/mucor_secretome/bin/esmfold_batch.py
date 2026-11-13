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


def set_chunk(model, chunk):
    """Chunked attention, which is what keeps a long protein inside GPU memory.

    fair-esm exposes set_chunk_size() on the model; transformers puts it on model.trunk.
    """
    for obj in (model, getattr(model, "trunk", None)):
        fn = getattr(obj, "set_chunk_size", None)
        if callable(fn):
            fn(chunk)
            return True
    return False


def load_model(backend, device, chunk, torch):
    """Return (model, infer) where infer(seq) -> (pdb_string, plddt_0_100, ptm).

    Two implementations of the same network:

    * **fair-esm** (`esm.pretrained.esmfold_v1`) needs the `openfold` package and its compiled
      CUDA kernels. That dependency is why the cluster's esm/2.0.0 module is frozen at Python 3.7
      and torch 1.13.1+cu117, which restricts it to the a100.
    * **transformers** (`EsmForProteinFolding`) vendors openfold's structure module in pure
      torch, so it needs no openfold install and runs on whatever torch it is given - which is
      how the new env reaches the ada6000, h100 and blackwell6000 cards.

    They are the same weights but not bit-identical in output, so models from the two backends
    should not be mixed in one set without saying so.

    pLDDT is normalised to 0-100 for both. fair-esm reports `mean_plddt` on that scale already;
    transformers returns a per-atom 0-1 tensor, from which the CA column (index 1) is averaged.
    bin/fold_report.py and the report's >= 70 threshold assume the 0-100 scale.
    """
    if backend == "fair-esm":
        import esm
        model = esm.pretrained.esmfold_v1().eval()
        if device == "cuda":
            model = model.cuda()
        else:
            # esmfold_v1 keeps its ESM-2 trunk in fp16, and CPU has no half-precision LayerNorm
            # ("LayerNormKernelImpl not implemented for 'Half'"), so everything is cast to fp32.
            # That roughly doubles the memory and is the price of running without a GPU.
            model = model.cpu().float()
            if hasattr(model, "esm"):
                model.esm = model.esm.float()
        set_chunk(model, chunk)

        def infer(seq):
            out = model.infer(seq)
            return (model.output_to_pdb(out)[0],
                    float(out["mean_plddt"][0]),
                    float(out["ptm"][0]))

        return model, infer

    from transformers import AutoTokenizer, EsmForProteinFolding
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.eval()
    if device == "cuda":
        model = model.cuda()
        model.esm = model.esm.half()      # the ESM-2 trunk runs in fp16 on GPU
    else:
        model = model.cpu().float()
        model.esm = model.esm.float()
    set_chunk(model, chunk)

    def infer(seq):
        enc = tok([seq], return_tensors="pt", add_special_tokens=False)
        ids = enc["input_ids"]
        if device == "cuda":
            ids = ids.cuda()
        out = model(ids)
        pdb = model.output_to_pdb(out)[0]
        plddt = float(out["plddt"][0, :, 1].mean()) * 100.0
        ptm = float(out["ptm"])
        return pdb, plddt, ptm

    return model, infer


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=400)
    ap.add_argument("--chunk", type=int, default=128)
    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"],
                    help="cpu is the fallback when no compatible GPU is free; it needs no "
                         "architecture check, just patience")
    ap.add_argument("--limit", type=int, help="stop after this many proteins (for timing runs)")
    ap.add_argument("--backend", default="auto", choices=["auto", "transformers", "fair-esm"],
                    help="which ESMFold implementation. 'transformers' needs no openfold package "
                         "and so follows whatever torch it is given, which is what unlocks the "
                         "non-a100 GPUs; 'fair-esm' is the esm/2.0.0 module (a100 only). "
                         "'auto' prefers transformers when it is importable.")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    import torch

    # The GPU must be one this torch was built for, or every result is silently wrong. The
    # esm/2.0.0 module ships torch 1.13.1+cu117, built for sm_37..sm_86, so it runs correctly on
    # an a100 (sm_80) but not on the ada6000 (sm_89), h100 (sm_90) or blackwell6000 (sm_120)
    # nodes. On a Blackwell node it does not refuse: it prints a warning and then returns
    # garbage, which in job 32033882 meant 487 nonsense errors and - worse - 72 "successful"
    # models that were also nonsense. Fail loudly instead.
    #
    # The check stays even on the new env (torch 2.8.0+cu128, see scripts/setup_gpu_env.sh),
    # because it costs nothing and the failure mode it catches is silent.
    if a.device == "cpu":
        print(f"running on CPU, torch {torch.__version__}", file=sys.stderr)
    elif not torch.cuda.is_available():
        sys.exit("no CUDA device visible")
    if a.device == "cuda":
        cap = torch.cuda.get_device_capability()
        sm, arch = f"sm_{cap[0]}{cap[1]}", torch.cuda.get_arch_list()
        name = torch.cuda.get_device_name(0)
        print(f"device {name} ({sm}); torch {torch.__version__} built for {arch}",
              file=sys.stderr)
        if arch and sm not in arch:
            sys.exit(f"REFUSING TO RUN: {name} is {sm}, but torch {torch.__version__} was built "
                     f"for {arch}. Results would be silently wrong - request a supported GPU "
                     f"(--gres=gpu:a100:1).")

    backend = a.backend
    if backend == "auto":
        try:
            import transformers  # noqa: F401
            backend = "transformers"
        except ImportError:
            backend = "fair-esm"
    print(f"backend {backend}", file=sys.stderr)

    model, infer = load_model(backend, a.device, a.chunk, torch)
    log = os.path.join(a.out, "plddt.tsv")
    new = not os.path.exists(log)
    fh = open(log, "a")
    if new:
        fh.write("protein_id\tlength\tplddt\tptm\tseconds\n")

    n_done = n_skip = n_fail = n_tried = 0
    for name, raw in read_fasta(a.fasta):
        # attempts, not successes: a --limit that only counts successes never stops a run in
        # which everything fails, which is exactly when stopping matters
        if a.limit and n_tried >= a.limit:
            break
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
        n_tried += 1
        for attempt in range(3):
            try:
                t0 = time.time()
                with torch.no_grad():
                    pdb, plddt, ptm = infer(seq)
                dt = time.time() - t0
                with open(pdb_path, "w") as pf:
                    pf.write(pdb)
                fh.write(f"{name}\t{len(seq)}\t{plddt:.2f}\t{ptm:.4f}\t{dt:.1f}\n")
                fh.flush()
                n_done += 1
                break
            except torch.cuda.OutOfMemoryError:  # noqa: B902 - cuda only
                torch.cuda.empty_cache()
                chunk = max(8, chunk // 2)
                set_chunk(model, chunk)
                print(f"OOM on {name} ({len(seq)} aa), chunk -> {chunk}", file=sys.stderr)
            except Exception as e:  # noqa: BLE001 - one bad protein must not end the run
                print(f"FAIL {name}: {type(e).__name__}: {e}", file=sys.stderr)
                n_fail += 1
                break
        else:
            print(f"FAIL {name}: still OOM at chunk {chunk}", file=sys.stderr)
            n_fail += 1
        set_chunk(model, a.chunk)
    fh.close()
    print(f"folded {n_done}, skipped {n_skip}, failed {n_fail} -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
