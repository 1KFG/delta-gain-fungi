# mucor_secretome — handoff (2026-11-13)

Branch `gpu-env-esmfold`, pushed, **not merged**. Previous branch `contam-protein-check` carries
PR #8 (protein-level *Fusarium* result), also unmerged.

## 1. What changed this session

### ESMFold is no longer limited to the a100

The cluster's `esm/2.0.0` module is a self-contained **Python 3.7** env
(`/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env`) carrying torch 1.13.1+cu117, built for
sm_37..sm_86. That restricted every fold to the a100 (sm_80). Python 3.7 is why it cannot be
upgraded in place, and the dependency holding it there is **openfold**, which fair-esm's
`esmfold_v1` needs for its compiled CUDA kernels.

`transformers.EsmForProteinFolding` vendors the same structure module in pure torch, so it needs
no openfold install and follows whatever torch it is given. That is the whole fix.

**New env, built by `scripts/setup_gpu_env.sh`:**
`/bigdata/stajichlab/jstajich/envs/esmfold-torch28` — python 3.11, **torch 2.8.0+cu128**,
transformers 4.57.6, accelerate, scipy, biopython. Weights prefetched to
`/bigdata/stajichlab/jstajich/cache/huggingface` on the **login node**, because compute nodes may
have no outbound network. Nothing was written into `/bigdata/stajichlab/shared/`.

Two traps the script now handles, both of which cost time to find:

* The miniconda3 module's package cache is under `/opt` and read-only, so mamba cannot take its
  lockfile and the solve **stalls silently** with repeated `Could not open lockfile` lines and no
  further progress. `CONDA_PKGS_DIRS` is set to a directory we own.
* `conda activate --stack` inside the esm module is the shell *function*, not the executable;
  without `eval "$(conda shell.bash hook)"` it fails with "unrecognized arguments: --stack"
  (job 32033764).

**GPUs unlocked.** Across `short_gpu` / `preempt_gpu` / `gpu` / `exfab`: 24 a100 were usable
before; 42 blackwell6000 (sm_120), 16 ada6000 (sm_89) and 2 h100 (sm_90) were not. The k80
(sm_37) and p100 (sm_60) nodes remain below torch 2.8's floor of sm_70 — **do not request a bare
`--gres=gpu:1`**, name a card type, or the job can land somewhere the guard correctly aborts.

### Verified on real silicon, not assumed

`scripts/gpu_arch_probe.sh` checks three things and only the third counts: `sm_XX` in
`get_arch_list()`, a GPU matmul agreeing with a CPU reference, and ubiquitin folding above
pLDDT 0.7 (garbage output sits at 0.2–0.3).

| card | sm | in arch list | matmul err | ubiquitin pLDDT | job |
|---|---|---|---|---|---|
| A100-SXM4-80GB | sm_80 | YES | 8.4e-05 | **0.905** | 32062542 |
| RTX 6000 Ada | sm_89 | **NO** | 5.3e-05 | **0.905** | 32062543 |
| RTX PRO 6000 Blackwell | sm_120 | YES | 7.8e-05 | **0.905** | 32062546 |
| h100 | sm_90 | — | — | — | 32062545, still queued |

**The ada6000 is the result worth remembering: sm_89 is absent from torch 2.8.0+cu128's arch list
and the card works anyway**, identically to the others. The arch list is sufficient, not
necessary. All three give exactly 0.905, so ESMFold is deterministic across these architectures.

### A pLDDT units bug, found by the comparison step and fixed

`fold_effectors.sh` prints a backend comparison whenever both runs exist. First run said:

```
mean pLDDT this 57.9 vs other 54.9; mean diff +3.0, median +3.0, max |diff| 6.7
proteins crossing the pLDDT 70 line between backends: 31 of 559
```

A uniform offset in one direction is not what two implementations of the same weights produce.
Cause: `transformers` returns `plddt` as `(batch, residues, 37)` — per **atom** — with
`atom37_atom_exists` as the mask, and fair-esm's `mean_plddt` is the mask-weighted mean over all
37 slots. The first version averaged the CA column only. Backbone atoms are predicted better than
side chains, so CA-only runs high by about 3 points — a units bug dressed as a backend
difference, and one that would have moved 31 proteins across the threshold `fold_report.py`
filters on. After the fix (job 32063774):

```
mean pLDDT this 54.9 vs other 54.9; mean diff -0.0, median +0.0, max |diff| 5.1
proteins crossing the pLDDT 70 line between backends: 1 of 559
```

**Keep that comparison step.** It is nearly free and it is the only thing that caught this.

### Folder convention — the two runs never mix

Same weights, different implementation, so the models are not bit-identical and must not be
pooled in one set.

| `ESM_BACKEND` | models | downstream |
|---|---|---|
| `fair-esm` | `out/effectors/folds/` — the original 559, a100 only | `hits_<db>.tsv`, `fold_calls.tsv` |
| `transformers` (**default** from 2026-11-13) | `out/effectors/folds_hf/` | `hits_<db>_hf.tsv`, `fold_calls_hf.tsv` |

## 2. State of the fold calls

Both backends, 559 models each, 640 attempted (81 skipped above the 400 aa cap).

| fold_call | fair-esm | transformers |
|---|---|---|
| low_confidence | 385 | 389 |
| assigned | 132 | 133 |
| no_match | 42 | 37 |

**18 of 559 change class between backends.** That is the honest cost of the switch and it belongs
in the report's methods, not in a footnote.

The headline call survives unchanged: `JEM1MUCOR_001522-T1`, no Pfam, is assigned to **cerato-platanin
1 of *M. perniciosa* (MpCP1)**, TM **0.844** in both runs. But its pLDDT is **71.12 (fair-esm) /
72.01 (transformers)** — barely over the 70 cutoff in both. That call is threshold-fragile and
should be described as such, not quoted as if it were comfortable.

## 3. Not done / open

* **`REPORT.md` and the published page are not updated** for any of this. The rendered report at
  https://dgf.fungalgenomes.org/reports/mucor_secretome.html still describes the fair-esm run
  only. The backend switch, the 18 reclassifications and the threshold-fragility of the
  cerato-platanin call all need writing in before it is re-rendered.
* **Foldseek was re-run for the `_hf` models** (`hits_*_hf.tsv` exist), but nothing downstream of
  `fold_calls_hf.tsv` has been regenerated — the AF3 panels were built from the fair-esm set.
  Decide whether to rebuild them from `_hf` or keep the fair-esm provenance for the AF3 work and
  say so explicitly.
* **AF3 MSA arrays**: 32050397 (chitin panel) **COMPLETED**; 32050398 (427 unknowns) was still
  running at last check. `slurm/af3_fold.sh` with `SET=chitin` can be submitted now; `SET=unknowns`
  waits on 32050398.
* **h100 probe 32062545 still queued.** Not blocking — the other three passed.
* **PR #8** (`contam-protein-check` → main) open, unmerged. Branch `gpu-env-esmfold` likewise.
* The repo root still holds a junk untracked file with a mangled name from a botched heredoc, plus
  an untracked `test/`. Not deleted — never asked.

## 4. How to resume

```bash
# on the HPCC
cd /bigdata/stajichlab/jstajich/projects/delta-gain-fungi
git checkout gpu-env-esmfold && git pull

# fold with the new backend (name a card type; do not use a bare gpu:1)
sbatch --gres=gpu:blackwell6000:1 analysis/mucor_secretome/slurm/fold_effectors.sh

# fold with the old one, for comparison only
sbatch --gres=gpu:a100:1 --export=ALL,ESM_BACKEND=fair-esm \
    analysis/mucor_secretome/slurm/fold_effectors.sh
```

Throughput, for planning: **640 attempted / 559 folded in ~13–20 min on one Blackwell**, against
hours on an a100 queue that was usually full. Folding is no longer the bottleneck for sets of this
size; for anything larger the MSA step is, and it is CPU-bound (the AF3 arrays ran 40 min to 5.5 h
per task).

Tests: **78 pass** (`python3 -W ignore -m unittest discover -s analysis/mucor_secretome/test`),
up from 68 — the ten new ones cover the chunk-size dispatch, which lives on the model in fair-esm
and on `model.trunk` in transformers, and the pLDDT threshold conversion in both directions.

## 5. Related, outside this directory

`HANDOFF.md` at the repo root was refreshed the same day: the delta-gain funnel finished
2026-11-06 (job 31780270) and was **closed by user decision 2026-11-08**, while the repo had
stopped tracking it on 11-05 and would have had anyone quoting v1 numbers as current. An Opus
review drafted a successor plan (orphan triage, lichen convergence, a structural census of the
secretome with effector scaffolds as labels, disorder in BUSCO orthologs). That draft is **not yet
written to `analysis/NEXT_STEPS.md`** — it is awaiting critique, and decisions D1–D3 in it are
open. Two of its verified corrections bear on structure work here: EarlGrey repeat annotation
covers **26 genomes**, not the corpus, and the corpus secretome is **18.77 M** signal-peptide
proteins out of ~249 M, roughly 25× larger than the estimate it replaced.
