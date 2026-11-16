# Bugs found and fixed, 2026-11-12/13/14

A log of defects found while moving ESMFold off the a100 and running the AF3 screen. Recorded
because **every one of them failed silently** — the jobs exited 0, or failed with an error that
pointed nowhere near the cause — and several would have produced a plausible-looking result.

The common shape, worth internalising: **derived state outliving the thing it was derived from**,
and **a check that can only ever return one answer**. Four of the eight are one or the other.

---

## 1. torch too old for the GPU — silent garbage, not a crash

**Where** the `esm/2.0.0` module (Python 3.7, torch 1.13.1+cu117, built sm_37–86).

**Symptom** job 32033882 on a Blackwell: 487 errors and, worse, **72 "successful" models that were
also nonsense**. torch does not refuse an unsupported architecture; it warns and returns rubbish.

**Fix** `bin/esmfold_batch.py` aborts when the card's `sm_XX` is absent from
`torch.cuda.get_arch_list()`. Then the real fix: a new env (`scripts/setup_gpu_env.sh`,
torch 2.8.0+cu128) and the `transformers` ESMFold backend, which needs no `openfold` and so is not
pinned to Python 3.7. 24 usable GPUs → 84.

**Caveat found later** the arch list is *sufficient, not necessary*: the ada6000 (sm_89) is absent
from torch 2.8's list and works correctly. The probe checks behaviour, not the list.

## 2. mamba stalls silently on a read-only package cache

**Where** building the env from the `miniconda3` module.

**Symptom** four `Could not open lockfile .../pkgs/cache/cache.lock` lines, then no progress ever.
Not an error exit — a hang.

**Fix** `CONDA_PKGS_DIRS` set to a directory we own. Will bite any conda env built from that module.

## 3. pLDDT averaged over the wrong atoms — a units bug dressed as a finding

**Where** `bin/esmfold_batch.py`, transformers backend.

**Symptom** the new backend scored a systematic **+3.0** mean pLDDT above fair-esm, max |diff| 6.7,
and **31 of 559 proteins crossed the pLDDT 70 line** that `fold_report.py` filters on. It looked
like a genuine difference between implementations.

**Cause** `transformers` returns `plddt` as `(batch, residues, 37)` — per *atom* — with
`atom37_atom_exists` as the mask. fair-esm's `mean_plddt` is the mask-weighted mean over all 37
slots. We averaged the CA column; backbone atoms are predicted better than side chains, hence the
offset.

**Fix** compute `(plddt * atom37_atom_exists).sum() / atom37_atom_exists.sum() * 100`. After:
mean diff **−0.0**, 1 of 559 crossing.

**Caught by** the backend-comparison step in `slurm/fold_effectors.sh`, which exists only because
the two runs were kept in separate folders. Keep it.

## 4. AF3 monomers folded with no MSA at all

**Where** `slurm/af3_fold.sh`, unknowns branch.

**Symptom** job 32069699: **0 of 426 succeeded.** AF3 died in featurising. Shards "finished" in
8 minutes instead of ~33, which is the tell — failing fast looks like succeeding fast.

**Cause** jobs built with `af3_input.py single` and **no `--msa-dir`**, then run with
`--norun_data_pipeline`. The script's own comment said "fold each with the MSA its stage-1 run
produced"; the code never passed the flag.

**Fix** `single` takes an optional `--msa-dir`; `af3_fold.sh` passes it. Stage 1 (`af3_msa.sh`)
still calls it *without* one — correct, that run exists to build the MSAs.

## 5. Control chains renamed under their own MSAs

**Where** `slurm/af3_prepare.sh`.

**Symptom** job 32069698: 170 of 196 fine, but **all 26 control predictions failed** — and the
controls are the calibration the whole panel rests on.

**Cause** controls were named `CTRL_<pdb>_<n>` **by position** in an RCSB fetch that the script
repeats every run. Rebuilding the panels re-fetched them, moving names onto different chains while
`msas/` still held the old alignments. The positive control pair went out as a 68 aa sequence
carrying the 75 aa chain's MSA, and vice versa. Proteins with real locus tags were immune, which
is exactly why only the controls broke.

**Fix** two ends. Control names are now a sequence hash (`CTRL_6Q76_265aaf`), stable across fetches
and RCSB reordering. And `msa_for()` compares the MSA's query sequence with the chain it is about
to be attached to, refusing the pair rather than emitting a job that dies later.

## 6. Stale job JSONs re-folded

**Where** `slurm/af3_fold.sh`.

**Symptom** after the rename, 23 old-name control pairs were folded again and failed again; their
empty `pred/` directories were then counted as missing and **blocked the report**.

**Cause** `af3_prepare.sh` clears `jobs/` (stage-1 inputs). `pairs/` is a different directory,
written by `af3_fold.sh`, and nothing cleared it.

**Fix** task 0 clears `pairs/*.json` before regenerating; completeness is judged against the
current job list, not against whatever `pred/` holds.

## 7. Every resume check was dead — AF3 does not lowercase its output

**Where** `slurm/af3_msa.sh` and `slurm/af3_fold.sh`.

**Symptom** reruns never skipped anything. Resubmitting a 426-protein set to finish **one** missing
prediction hit the 2 h wall (job 32080675, four shards TIMEOUT). And the completeness gate never
fired, so **no report was ever written even when all 196 and all 426 predictions existed**.

**Cause** both scripts tested `$OUT/<dir>/${n,,}/${n,,}_*.json`. AlphaFold 3.0.2 writes the job
name **as given**: `pred/CTRL_4EBZ_03a7bf__CO4/CTRL_4EBZ_03a7bf__CO4_summary_confidences.json`.
Our names are upper case, so no test ever matched. The "resumable" property the scripts advertised
had never once worked.

**Fix** a shared `af3_done()` accepting either spelling — which `bin/af3_input.py`'s `find_msa`
already did, so the inconsistency was internal.

**Cost** every rerun in this series recomputed from scratch, including a 3.5–4.2 h MSA job that
should have rebuilt three chains.

## 8. The report counted model samples as results

**Where** `bin/af3_report.py`.

**Symptom** the first complete chitin report: **1,212 rows from 196 jobs**.

**Cause** `os.walk` descended into AF3's `seed-N_sample-M/` subdirectories, each carrying its own
`*_summary_confidences.json`, and used `basename(dirpath)` as the job name. **1,010 of 1,212 rows
were per-sample duplicates labelled `seed-1_sample-N`**, unattributable to any pair — and the
report sorts by score, so they ranked alongside real results. Separately, six `CTRL_4EBZ_1__*`
directories left from before the rename were reported beside their renamed twins with identical
scores.

**Fix** look exactly one level down and take the per-job summary only; `--jobs-dir` supplies the
jobs actually requested and anything else is ignored with a note.

---

## What generalises

1. **A check that always returns the same answer is invisible.** #7 silently disabled resumption
   everywhere for the life of the scripts. A test that can only pass is worse than no test.
2. **Fast failure resembles fast success.** #4's 8-minute shards looked like good news.
3. **Derived directories must be cleared by whatever derives them** — `jobs/`, `pairs/`, `pred/`
   (#5, #6, #8). Three separate instances before the pattern was obvious.
4. **Controls break differently from data.** In #5 and #8 the real proteins were fine and only the
   controls were wrong, because the controls were the only things named positionally.
5. **The cheap comparison earns its place.** #3 was caught solely by printing two backends side by
   side. The "no report unless complete" gate (#6) is what made #7 visible at all.
6. **Numbers quoted before the pipeline was trusted need rechecking**, including in `REPORT.md`
   and anything already published.

Tests over this work: **68 → 100**.
