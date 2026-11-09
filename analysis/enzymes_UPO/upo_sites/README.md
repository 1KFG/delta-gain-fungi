# upo_sites

Find PF01328 (heme-thiolate peroxidase) proteins in a proteome, read their active-site residues,
score the aliphatic pocket, and state whether the score is entitled to speak about them.

## One command

```bash
module load hmmer/3.4
/usr/bin/python3.12 -P upo_sites/bin/upo_scan.py \
    --proteome yourproteome.fa.gz \
    --hmm      upo_id/models/PF01328.hmm \
    --backbone data/clade_backbone/raw.fa \
    --s2       data/calib_dl/active_site_S2.tsv \
    --out      out.tsv.zst --cpu 8
```

Accepts plain, `.gz` or `.zst` FASTA, and either the raw or the aligned (`.afa`) backbone.
Several `--proteome` files can be given at once.

## Output columns added by upo_scan

| column | meaning |
|---|---|
| `t1_pocket_matches` | matches (0-5) to the active consensus `ILLAI` at HspUPO p73, p77, p81, p183, p187 |
| `t1_pocket_flag` | `clade_V_like` at >= 3 matches |
| `calib_distance` | fraction of differing residues, at 10 active-site positions the score does NOT use, to the nearest of the 19 characterized enzymes |
| `calib_domain` | `inside` at distance <= 0.50, the largest such distance within the calibration set |

## What the score is, stated plainly

Cross-validated on the 19 characterized enzymes: LOO R2 +0.496, LOO AUC 0.956, permutation
p = 0.0005. It beats a clade-only model out of sample.

**Its limit is not a formality.** 9 of the 10 active calibration enzymes are clade V. Remove
clade V and the AUC is 0.556, which is chance, on ONE remaining active. The score detects
clade-V-like aliphatic pockets. It is not shown to find actives outside clade V, and it misses
the one that exists: GmaUPO-II is 49.3 % active and carries 1 of 5 matches.

So a high score is a hypothesis worth assaying. **A low score is not evidence of inactivity.**
The substrate is 4-propylguaiacol conversion only; the score says nothing about ring cleavage.

`calib_domain` is a COVERAGE statement, not a confidence. Distance was tested against model
error and did not predict it (Spearman -0.281, p = 0.88). It tells you whether a protein's
active site resembles anything that has been measured, which is what you want when choosing
what to assay. It does not tell you the odds the score is right.

Across BFD, 69.6 % of proteins lie outside the calibration domain, as do 40.1 % of those
scoring >= 3/5.

## What upo_scan does not do

No phylogenetic placement, so there is no `clade` column and the score's central limit cannot
be checked per protein. Placement needs the reference tree and RAxML EPA and declines 78.9 % of
BFD. Run `make_placement_aln.py` then `jplace_to_clade.py` when a clade call is needed.

## Scripts

| script | role |
|---|---|
| `upo_scan.py` | the one-command entry point: sites + score + domain |
| `upo_sites.py` | site layer only: PF01328 search, hmmalign, active-site residues |
| `score_t1.py` | builds and cross-validates the score on the 19 enzymes |
| `applicability.py` | tests distance as a reliability measure, bounds the domain |
| `apply_t1_score.py` | applies the score to an existing BFD table |
| `build_assay_panel.py` | picks a 96-well panel for information, not confirmation |
| `merge_final.py`, `jplace_to_clade.py`, `make_placement_aln.py` | placement and merge stages |

## Tests

```bash
UPO_ROOT=$PWD bash upo_sites/test/test_upo_scan.sh
```

Runs the tool on 87 real BFD proteins and checks every call against what the full BFD pipeline
produced for the same proteins by a different route. Currently 87/87.

Measurements behind every claim here: `MODEL_VALIDATION.md` sections D, E and F,
`T1_SCORE_REPORT.md`, `APPLICABILITY.md`.
