# Census steps 3 and 4 (carotenoid): results, decisions and evidence (2026-11-22)

Status: carotenoid only. All outputs are PROVISIONAL (pathway `verified: false`). Nothing here is a biological finding.
Code: commit `32f2f90` on branch `pathway-census-realdata` (worktree `DeltaGain_Fungi_wt_realdata`), tools `bin/make_min_scores.py`,
`bin/make_chunks.py`, `bin/search_chunk.py`, `slurm/search_chunk.sbatch`. Method description: `analysis/pathway_census/README.md`, section "Batched BFD search".
Tool versions: hmmer 3.4 (`module load hmmer/3.4`), mafft 7.505, Python 3.12 (`/usr/bin/python3.12`), zstd 1.5.6. Test suite: 687 passed (2026-11-22).

## Inputs

| Item | File | Note |
|---|---|---|
| Step HMMs (14) | `analysis/pathway_census/results/hmms/carotenoid/*.PROVISIONAL.hmm` | built by `bin/build_step_hmms.py --provisional`; 7 of 14 rest on one sequence |
| Step 2 scores | `results/step2/carotenoid/step_scores.PROVISIONAL.tsv`, `leave_one_out.PROVISIONAL.tsv` | `bin/score_refs_decoys.py`, job 32375043 |
| Rule cutoffs | `results/step2/carotenoid/min_scores.PROVISIONAL.tsv` | owner rule of 2026-11-22 |
| Cutoffs used in the run | `results/step2/carotenoid/min_scores.exploratory25.PROVISIONAL.tsv` | rule values, except `astaxanthin_synthase` 1332.5 to 999.4 and `ylo1` 640.7 to 480.5 (x0.75) |
| Genome table | `results/genomes.tsv.zst` (24,363 rows) | `python -m pathway_census.genomes`, 49 s |
| Chunk manifest | `results/chunk_manifest.size4000.tsv.zst`, excluded list `results/chunk_excluded.tsv` | 23,376 genomes listed, 987 excluded (984 no seqdb path, 3 file missing), seed 1 |

All `results/` paths are git-ignored on purpose (data are not committed).

## Step 3: timing (job 32375142, 1 chunk of 200 genomes, 16 CPUs)

Job time 1 min 31 s, script time 78.9 s, peak memory 1.48 GB, 2,099,659 proteins, 1,827 candidates. Files: `results/search/timing_size200/`.

## Step 4: full carotenoid search (jobs 32375296 to 32375301, short queue, 16 CPUs, chunks of 4,000 genomes)

| Item | Value (MEASURED) | Source file |
|---|---|---|
| Jobs | 6 of 6 COMPLETED, 23 min 46 s to 28 min 0 s | `sacct`; `results/search/carotenoid_full/timing.chunk0000N.tsv` |
| Genomes searched | 23,376 of 23,376 listed, status `searched` for all, 0 genomes with an id-prefix error | `searched.chunk0000N.tsv` |
| Proteins searched | 244,451,335 | same |
| Peak memory per job | 22.5 to 22.7 GB (32 GB requested) | `sacct` MaxRSS |
| Candidate rows | 238,898 (protein by step) | `candidates.chunk0000N.PROVISIONAL.tsv.zst` |
| Estimate check | 26 min predicted from the 200-genome chunk, 24 to 28 min observed | step 3 timing |

Rows per step at the cutoffs used: acaa 26,611; astaxanthin_synthase 12; cars 28,663; cart 5,533; carx 11,528; crtr 46,698; ggpp_synthase 22,924;
lycopene_cyclase 9,530; phytoene_desaturase 9,784; phytoene_synthase 10,863; tdh 16,244; tsp2 15,976; tsp3 20,942; ylo1 13,590.
Single-reference steps with low floors (crtr 68.1, cars 8.0, acaa 12.3, tsp3 21.9) give the large counts. They need tree placement (step 5).

## Side evaluation: what a high cutoff misses (job 32375378, first 1,000 genomes of chunk 0, floor 20 bits for two steps)

Files: `results/search/side_eval_low_cutoff/` (cutoff file `min_scores.side20.tsv`, manifest `results/chunk_manifest.size1000.tsv.zst`).
The 1,000 genomes are a subset of chunk 0 of the full run, so the two runs compare directly. Every hit of the full run was also in the side run (0 missing).

| Step | Rule cutoff | Used in full run | Full run, 23,376 genomes | Side run floor 20, 1,000 genomes |
|---|---|---|---|---|
| astaxanthin_synthase | 1332.5 | 999.4 | 12 hits in 12 genomes: 8 at or above the rule value, 4 at 1216 to 1300 | 47,640 proteins in 992 genomes |
| ylo1 | 640.7 | 480.5 | 13,590 hits: 11,194 at or above the rule value, 2,396 at 480.5 to 640.6 (2,219 genomes) | 20,276 proteins in 990 genomes |

Score bands in the side run (proteins). astaxanthin_synthase: 20-50 37,255; 50-100 6,989; 100-200 2,529; 200-400 859; 400-700 6; 700-1000 0; 1000-1400 2.
ylo1: 20-50 1,208; 50-100 4,168; 100-200 13,309; 200-400 636; 400-700 659; 700-1000 214; 1000-1400 82.

Reading (a reading of counts, not a finding): both low floors return tens of proteins per genome, so they reach a large family, not a single-copy enzyme.
Whether any low-scoring hit is a true ortholog needs the step 5 trees. The identity of those proteins was not checked.

## Decision (owner, 2026-11-22): keep the lowered cutoffs

- Keep 999.4 for `astaxanthin_synthase` and 480.5 for `ylo1` in step 5 placement. Do not go lower, because lower floors flood placement with family members.
- Reason on record: the full run adds 4 genomes (astaxanthin_synthase) and 2,396 proteins (ylo1) over the rule value, while a floor of 20 adds 47,640 and 20,276 proteins per 1,000 genomes.
- The rule value stays recoverable: filter `score >= rule_min_score` (column in the cutoff file).
- Open: the rule cutoff for two-reference steps whose references are near-identical (astaxanthin_synthase, ylo1) says nothing about distant homologs.

## Species tree for step 5 and later steps (found 2026-11-22)

A genus-level phyling tree exists in the dashboard project. Not yet checked against the carotenoid candidates.

| Item | Path |
|---|---|
| ASTRAL consensus tree (tips are LOCUSTAG; branch lengths and support values present) | `/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/results/protein/tree/fungi_odb10/final_tree.nw` (link to `astral/consensus.nw`, written 2026-11-12) |
| FastTree concatenation trees, with and without support | `.../results/protein/buildtree/fungi_odb10/fasttree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fasttree.{support,nosupport}.treefile` |
| Genus representatives (1,654 rows: LOCUSTAG, GENUS, CLASS, busco_complete, pep_path) | `.../atlas_data/phylo/genus_representatives.tsv` |
| Method | `.../atlas_data/phylo/finish_atlas_tree.sh` (phyling filter, top 80 markers of `fungi_odb10`, concat, FastTree -lg with 1000 bootstraps; ASTRAL consensus) |

MEASURED by a regex scan of the newick: 1,643 tip tokens, all of the LOCUSTAG form and all in `genus_representatives.tsv` (1,654 reps). The log of the FastTree run says 1,561 taxa in the concatenation. The tip count of the ASTRAL tree and the FastTree trees was not compared. Root and outgroup were not checked.
The original Nextflow run was marked failed (PHYLING_ALIGN exit 137 after 4 h 25 min); the alignment had finished and the tree was completed by `finish_atlas_tree.sh` (log `atlas_data/logs/atlas_finish.32032053.log`).
Related plan, not built: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/docs/SPECIES_TREE_PLAN.md` (species-level tree, about 7,450 representatives).
Older fungal trees: 1KFG 2018 marker trees (`/bigdata/stajichlab/shared/projects/1KFG/2018-newanalysis/analysis_2018.v1/Phylo/run1/*.tre`), taxa named by JGI ids, not BFD locustags.
