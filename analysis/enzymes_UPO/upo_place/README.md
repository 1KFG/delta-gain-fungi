# `upo_place`

Clade calls for PF01328 proteins by **phylogenetic placement** on the in-house 187-leaf reference
tree (`../data/refset_zamocky2022/tree/`, 19 backbone enzymes + 168 Zamocky 2022 curated HTPs).
It replaces `upo_sites`' nearest-neighbour `clade` column, which was validated only at >= 52.7 %
identity to 19 enzymes and extrapolated below that. `upo_sites` remains the tool for residues.

## Method

1. **Reference** (`slurm/prepare_ref.sh`): `hmmbuild --symfrac 0` on the 224-column trimmed
   alignment, so every reference column is an HMM match state; raxml-ng `--evaluate` of LG+G4
   branch lengths and parameters on the ModelFinder (mfp) topology. ModelFinder chose Q.PFAM+R7,
   but epa-ng 0.3.8 rejects Q.PFAM; clade membership was identical under LG+I+G4, so the clades
   do not depend on the matrix. Output in `../data/refset_zamocky2022/placement_ref/`.
2. **Align** each query with `hmmalign` and keep match columns only (`upo_place.py match-cols`).
3. **Place** with epa-ng (tools via `pixi exec`, pinned in `slurm/env.sh`).
4. **Assign** (`upo_place.py assign`): sum the like-weight ratio (LWR) of a query's placements
   that fall *inside* each named clade; a placement on a clade's *stem* does not count for it.
   Call the clade if the sum is >= 0.5, else `unassigned`. Root-independent (bipartitions). Both
   our Clades I-V and Zamocky's Figure 1 clades are reported. `long_pendant` flags a best
   placement whose pendant branch is > 1.0 (flag, not rejection).

## Validation: leave-one-out over all 187 reference leaves (`slurm/place_loo.sh`)

Each leaf is pruned, the HMM rebuilt without it, its raw sequence re-aligned and placed (the
path a BFD protein takes). Results: `../data/refset_zamocky2022/placement_ref/loo.tsv`.

| | correct |
|---|---|
| Clades I-V | **184/187** (2 abstain on a stem: MorUPO, OpsiHTP; 1 wrong: CfumagHTP -> II, pendant 1.42) |
| Figure 1 clades | 177/187 (8 of 10 misses are leaves of 2-leaf or basal positions that land on the clade stem when held out, and are left unassigned) |

By identity to the nearest of the 19 enzymes (the old caller's regime):
>= 52.7 %: 50/51; 30-52.7 %: 102/104; < 30 %: 32/32. Of correct assigned calls, 72/73 have
LWR >= 0.95; correct pendant lengths reach 0.88 (median 0.17).

**Against the old caller on the same 168 Zamocky leaves** (neither has them as references):
placement 166/168; nearest-neighbour 79/168, having given a clade to 81 proteins that are in none
of Clades I-V.

## BFD (97,520 PF01328 proteins; jobs 31917770 + merge; `../data/bfd_place/COMPARE.md`)

| placement clade | n |
|---|---|
| Clade I (= Zamocky Family II, long) | 41,679 |
| Clade II | 3,547 |
| Clade III | 431 |
| Clade IV | 8,437 |
| Clade V | 2,425 |
| unassigned | 41,001 |

What changes against the old calls:

- **Validated-range calls are confirmed**: Clade I 3,731/3,731, II 403/403, IV 4,960/4,960,
  V 2,066/2,616 (531 -> unassigned). Clade III 361/1,079: the other 718 all sit on the stem of
  Clade III inside Family I, i.e. close relatives of MroUPO/MorUPO that branch outside the 3-leaf
  reference clade.
- **Extrapolated calls mostly do not survive outside Clade I.** Old Clade V extrapolated: 22,210
  -> 319 stay V, 2,834 -> IV, 1,316 -> II, 17,534 unassigned. Old Clade III extrapolated: 9,075 ->
  40 stay III, 7,625 unassigned. Old Clade I extrapolated: 29,633/29,693 confirmed.
- **Clade IV grows** from 4,968 validated to 8,437, mostly from old Clade V calls.
- 6,941 proteins below the old 30 % identity floor are now placed in Clade I.
- Figure 1: 6,378 proteins in Basidiomycetous 1st, 6,228 Ascomycetous 1st, 4,394
  Basidiomycetous 2nd, 784 glomeromycete, 315 Mucorales basal: lineages the 19-enzyme backbone
  had no reference for.

### Glu180 in Clade I is lineage-structured, not contamination

| Clade I, by phylum | Glu at HspUPO 180 | other |
|---|---|---|
| Ascomycota | 507 | **33,276** (Ala 29,749 overall) |
| Basidiomycota | 5,280 | 2,615 |

99.8 % of old Clade I calls (33,364/33,424) place in Family II. The non-Glu residue is almost
entirely in ascomycete long UPOs, matching the 13 curated ascomycete Family II sequences (all
A/S/T). This contradicts the reading in `../G1_BFD_SUMMARY.md` §6 that the low-identity Clade I
bin is absorbing non-UPOs. Open: whether ascomycete long UPOs lack the acid-base Glu or carry it
at a position that does not align to HspUPO 180 (check against a structure-based alignment of a
long ascomycete UPO).

Also: **Basidiomycetous 1st has no Glu180 in 6,367/6,378** BFD proteins, as in all 34 curated.

## Files

```
bin/upo_place.py      match-cols | prune | assign | loo
bin/compare_bfd.py    merge chunk calls with upo_sites output; COMPARE.md
slurm/env.sh          paths, model, pinned pixi tools
slurm/prepare_ref.sh  reference HMM + raxml-ng model
slurm/place_loo.sh    leave-one-out validation
slurm/place_bfd.sh    BFD array (reads shared chunks read-only, writes to this clone)
slurm/merge_bfd.sh    merge + compare
test/                 offline tests (python3 -m unittest discover -s upo_place/test)
```
