# Census step 5: one joint tree for the three PKS steps of melanin_dhn

Date: 2026-11-24. Branch `pathway-census-melanin-myco`. All results are PROVISIONAL (`pathways/melanin_dhn.yaml` has `verified: false`).
Labels: MEASURED = read from a file named here. PROPOSED = a choice. UNVERIFIED = not checked.
The target clades are tree regions bounded by the labelled leaves. They are not ortholog sets.

## Inputs

| Item | Value |
|---|---|
| Search | `results/search/melanin_myco_full/` (step 4). Rows are at the `min_scores` of that run. |
| Candidate folder, rows per step | `results/search/pks_joint/` (same file format; 7,286 + 8,604 + 2,200 rows) and `pks_scores.tsv` (score of each protein on each PKS HMM) |
| Candidate folder, one row per protein | `results/search/pks_joint_union/` (`step_id = pks_joint`; the row with the highest raw score over the three HMMs; 9,374 proteins). This folder is `PC_CAND_DIR` for both sbatch parts. |
| References | `inputs/refs_*.faa`: thn_synthase P79068; hexaketide_pks Q9Y7A7 (37 X characters); heptaketide_pks Q4WZA8 (alb1) and A2QUI2 (albA). Copies of `refs/melanin_dhn/<step>.faa`. |
| Decoys | `inputs/decoys_all.faa`: Q03149 (wA), Q12053 (aflC), Q12397 (stcA), Q6DQW3 (CTB1), S0DZM7 (bik1). The union of the three decoy files of `refs/melanin_dhn/decoys/`. |
| Named leaves | `inputs/named_members.faa` (18), table `inputs/named_leaves.tsv`. UniProt records fetched from rest.uniprot.org on 2026-11-24 (`inputs/fetched_records.tsv`, `inputs/fetched.faa`). For all 24 fetched records the length of the FASTA equals the length of the record and the OX value equals the organism id of the record (MEASURED). |
| Outgroups | None. The 4 outside candidates failed the alignability rule (below). |
| Alignability screen | `inputs/screen_alignability.tsv` (33 proteins x 3 HMMs, `bin/screen_outgroups.py`, rule: best domain >= 10 bits and covering >= 0.5 of the HMM; best and union mode) |

### Alignability screen (MEASURED, `inputs/screen_alignability.tsv`)

- Included as named leaves (18): best domain score 407.8 to 2,230.4 bits, best-domain coverage 0.526 to 0.998 of the HMM, `pass_best` true on all three HMMs.
- Not included, because `pass_best` fails on one HMM: F1DBA9 (mpaC, Penicillium brevicompactum; hexaketide HMM 415.0 bits, coverage 0.407) and Q5ATJ7 (ausA, Aspergillus nidulans; hexaketide HMM 374.9 bits, coverage 0.407, union coverage 0.526 passes only in union mode).
- Not included, fail both modes on all three HMMs (outside candidates): Q0C8M3 (LovB, Aspergillus terreus; best-domain coverage 0.402 to 0.411), Q03131 (EryA1, Saccharopolyspora erythraea, 3,491 aa; 0.402 to 0.418), P07149 (FAS1 yeast, 2,051 aa; under 12 bits, coverage under 0.02), P19097 (FAS2 yeast, 1,887 aa; no hit at E 1000). So no fungal reducing PKS, bacterial type I PKS or fatty acid synthase passes. The tree has no outgroup (PROPOSED: use none rather than a leaf that fails the rule). The display root is the decoys (a fallback of `step5_place.py`); the clade rule does not use the root.
- The five decoys pass (best score 1,187.0 to 3,904.0; coverage 0.598 to 1.0).
- Alb1 (Q4WZA8) is a reference of heptaketide_pks, so it is not a named leaf or outgroup.

## Method (PROPOSED unless stated)

- Domain cut: best row per protein, `ali_from..ali_to` padded by 20 (as step 5). Union rows: 9,374 cut, 0 failures (MEASURED, `extract_summary.json`). Cut length 942 to 2,250 residues, median 2,160 (MEASURED, `domain_coords.tsv`).
- Clustering: mmseqs2 `easy-cluster`, `--min-seq-id 0.7 -c 0.8 --cov-mode 0`. Cluster counts for the same 9,374 domains (MEASURED, scratch test): 0.5 gives 21; 0.6 gives 93; 0.7 gives 237; 0.8 gives 532; 0.9 gives 1,129. Identity 0.5 (the default) gives 21 representatives, too few to place 9,374 proteins; 0.7 gives 237 (under the cap of 300, so all clusters are representatives and nothing is drawn). The first attempt at 0.5 failed before clustering was used (below; `attempt1_idt0.5_failed/` keeps its extraction and cluster files). 0.7 is a choice; 0.8 (532) would exceed the cap of 300 and need a draw.
- Tree set: 4 references + 5 decoys + 18 named leaves + 237 representatives = 264 leaves (MEASURED). mafft `--auto` (7.505), clipkit smart-gap: 4,090 to 3,658 columns (10.562% removed). Gaps 44.63% (MEASURED, `raxml.log`).
- Tree: raxml-ng 2.0.4, LG+G4, 4 threads, `--tree pars{2}` (2 starts), seed 1. Best log-likelihood -630,903.14 (start 2); start 1 gave -630,909.43. alpha 0.912. One near-zero branch (warning in `raxml.log`). Two starts only: the tree is not checked for stability (UNVERIFIED).
- Clade rule (`bin/step5_joint.py clades`, new): target(S) = the LARGEST tree side that holds all references of S and no decoy, named leaf, outgroup or reference of another step. This is the rule of `bin/step5_single_ref_clades.py` for one or more references. If no such side exists, S gets no clade file (result).
- Placement: `hmmbuild --symfrac 0` of the trimmed alignment (3,658 match states), `hmmalign` of the 9,374 cut domains, epa-ng 0.3.8, `upo_place assign --min-lwr 0.5 --max-pendant 1.0`, once per step with that step's clade file; then only the proteins with a row of S in `results/search/pks_joint/`; then `assign_steps --min-support 0.5 --provisional`.
- Variants for heptaketide_pks (see Results): the default rule leaves it without a target. Two widened targets were computed afterwards from the same tree and the same placement (`variant_wA/`, `variant_YWA1/`). They are PROPOSED variants, not the rule.

## Commands

```
cd analysis/pathway_census
/usr/bin/python3.12 bin/pks_joint_candidates.py --in-dir results/search/melanin_myco_full --steps thn_synthase,hexaketide_pks,heptaketide_pks \
  --joint-id pks_joint --steps-dir results/search/pks_joint --union-dir results/search/pks_joint_union --scores-out results/search/pks_joint/pks_scores.tsv
module load hmmer/3.4; H=results/hmms/melanin_dhn
/usr/bin/python3.12 bin/screen_outgroups.py --fasta inputs/screen_input.faa --hmm $H/thn_synthase.PROVISIONAL.hmm $H/hexaketide_pks.PROVISIONAL.hmm $H/heptaketide_pks.PROVISIONAL.hmm --out inputs/screen_alignability.tsv
export PC_ROOT=$PWD PC_PATHWAY=melanin_dhn PC_STEP=pks_joint PC_CAND_DIR=$PWD/results/search/pks_joint_union PC_OUT=$PWD/results/step5/melanin_dhn/pks_joint \
  PC_REFS_ALL=$I/refs_all.faa PC_REFS_STEPS="thn_synthase=$I/refs_thn_synthase.faa hexaketide_pks=$I/refs_hexaketide_pks.faa heptaketide_pks=$I/refs_heptaketide_pks.faa" \
  PC_DECOYS_FA=$I/decoys_all.faa PC_MEMBERS_FA=$I/named_members.faa PC_RAXML_THREADS=4 PC_RAXML_STARTS='pars{2}' PC_MIN_SEQ_ID=0.7 PC_MAX_REPS=300 PC_GENOMES=$PWD/results/manifests/genomes.tsv.zst
sbatch --export=ALL slurm/step5_prepare_joint.sbatch     # part 1
export PC_STEPS="thn_synthase hexaketide_pks heptaketide_pks" PC_STEP_CAND_DIR=$PWD/results/search/pks_joint PC_MIN_SUPPORT=0.5 PC_SCORES=$PWD/results/search/pks_joint/pks_scores.tsv
sbatch --export=ALL slurm/step5_place_joint.sbatch       # part 2 (same PC_CAND_DIR)
# variants (heptaketide_pks), from the saved placement:
bin/step5_joint.py clades ... --include heptaketide_pks=Q03149 --out-dir variant_wA          # and Q03149,P9WET9,I1RF58 -> variant_YWA1
upo_place.py assign epa_result.jplace --clades variant_*/clade_members.heptaketide_pks.tsv ; step5_joint.py filter-placements ; assign_steps.py
```
$I is `results/step5/melanin_dhn/pks_joint/inputs`. `PC_ROOT` is set to this worktree in every job (the sbatch files of the carotenoid run default to another worktree; the new sbatch files require `PC_ROOT`).

## Jobs (partition `intel`, 16 CPUs)

| Job | Part | State | Elapsed | Note |
|---|---|---|---|---|
| 32434626 | prepare, identity 0.5 | FAILED, 2 min 31 s | | Extraction worked (9,374 of 9,374). Clustering gave 21 clusters. `select-reps` stopped: `step5_prepare.sbatch` reads `$PC_ROOT/results/genomes.tsv.zst`, which does not exist in this worktree (the file is `results/manifests/genomes.tsv.zst`). The new sbatch reads `PC_GENOMES`. |
| 32434736 | prepare, identity 0.7 | COMPLETED, 49 min 37 s | rc 0 | Max RSS 1.06 GB. raxml-ng 46 min 49 s for 2 starts. |
| 32435959 | place | COMPLETED, 42 min 52 s | | Max RSS 2.87 GB. Most of the time is `hmmalign` (one thread). |

CPUs: at most 16 held by this task at a time. The `short` queue holds pending jobs of the same user (reason AssocGrpCpuLimit), so all jobs ran on `intel`.

## Checks (MEASURED)

- `assign.<S>.tsv` has 9,374 queries per step = `domains.faa` (9,374) = rows of `results/search/pks_joint_union/` (9,374).
- `candidates_step.<S>.tsv` proteins = `placements.<S>.tsv` rows: thn_synthase 7,286 = 7,286; hexaketide_pks 8,604 = 8,604; heptaketide_pks variants 2,200 = 2,200 (the `PC_CAND_DIR` lesson).
- Proteins whose placement is `target` for more than one step: 0 (`joint_report/summary.json`).

## Results

See `joint_report/*.tsv`, `clade_report.json`, `target_overlap.tsv`, `variant_*/` and the note `docs/superpowers/notes/2026-11-24-melanin-myco-pks-placement.md` (numbers with denominators).
