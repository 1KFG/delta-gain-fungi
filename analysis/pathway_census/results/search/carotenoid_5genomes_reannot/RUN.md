# Carotenoid search of the 5 re-annotated genomes

Date: 2026-11-25. Branch `pathway-census-melanin-myco`. PROVISIONAL (`verified: false`). MEASURED unless marked.
Purpose: repeat the full carotenoid search for 5 genomes with the current protein files. Note: `docs/superpowers/notes/2026-11-25-carotenoid-reannotated-5genomes.md`.

## Inputs
| Item | Value |
|---|---|
| Pathway file | `pathways/carotenoid.yaml` (identical to the realdata worktree copy, `diff`) |
| HMM dir | `/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/pathway_census/results/hmms/carotenoid` (14 step HMMs, mtime 2026-11-22 11:03). sha256: `hmms_and_scores.sha256` |
| Cutoffs | `.../results/step2/carotenoid/min_scores.exploratory25.PROVISIONAL.tsv`. This is the table the full search used (`docs/superpowers/notes/2026-11-22-step4-results-and-decisions.md` of that worktree; lowest candidate score per step in the full search matches it). Not the rule table. |
| Manifest | `chunk_manifest.5genomes.tsv.zst`: 5 rows, locustag = NEW prefix, path = `proteins_path` of `prefix_aliases.tsv`. Made by hand from that file, not by `make_chunks.py`. |
| Script | `bin/search_chunk.py` of this worktree (newer than the one used for the full search: it adds a check for empty `min_score`; the search logic is unchanged, from `diff`). `slurm/search_chunk.sbatch` is identical. |

## Command
```
cd analysis/pathway_census
export PC_ROOT=$PWD PC_PATHWAY=$PWD/pathways/carotenoid.yaml PC_HMM_DIR=<HMM dir> \
  PC_CHUNKS=$PWD/results/search/carotenoid_5genomes_reannot/chunk_manifest.5genomes.tsv.zst \
  PC_MIN_SCORES=<cutoff file> PC_OUT_DIR=$PWD/results/search/carotenoid_5genomes_reannot PC_PROVISIONAL=1 PC_CHUNK_ID=0
cd $PC_OUT_DIR && sbatch --export=ALL -p intel --cpus-per-task=16 --mem=32G --time=01:00:00 --job-name=pc_car5 $PC_ROOT/slurm/search_chunk.sbatch
/usr/bin/python3.12 curation_scripts/carotenoid_reannot_compare.py <realdata>/results/search/carotenoid_full . prefix_aliases.tsv \
  /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb/proteins results/search/carotenoid_5genomes_reannot/compare
```

## Result
Job 32461202, `intel`, node i21, COMPLETED, 51 s, MaxRSS 282,896 K (`sacct`). 5 genomes, 39,700 proteins, 0 bad id prefixes, 37 candidate rows (old full search: 34 rows for the same 5 genomes). Script time 10.4 s (`timing.chunk00000.tsv`). Temp files were in `$SCRATCH` (set by `search_chunk.py`).
Comparison (by exact sequence, old sequences from the 2026-11-06 seqdb snapshot): 33 unchanged, 1 lost, 4 gained. 0 genome-step presence changes. Files: `compare/rows_per_genome_step.tsv`, `compare/changed_rows.tsv`, `compare.log`.
Candidate table `candidates.chunk00000.PROVISIONAL.tsv.zst` is not committed.
UNVERIFIED: the full search used the same HMM bytes (no checksum was recorded then).
