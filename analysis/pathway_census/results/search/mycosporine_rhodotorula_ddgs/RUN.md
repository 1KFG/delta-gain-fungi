# ddgs best hit per genome, no score cutoff

Date: 2026-11-26. Branch `pathway-census-rhodotorula-ddgs`. PROVISIONAL (mycosporine pathway has `verified: false`).
Labels: MEASURED = read from a file named here. UNVERIFIED = not checked.

## Command (MEASURED, wall time 30.8 s)

Run inside SLURM job 32339135 (4 CPUs, `$SCRATCH` temp files, hmmer/3.4 module). Script: `curation_scripts/mycosporine/rhodotorula_ddgs_check.py`.

```
cd analysis/pathway_census
R=<melanin worktree>/analysis/pathway_census/results ; T=<tables worktree>/analysis/pathway_census/results
python3.12 curation_scripts/mycosporine/rhodotorula_ddgs_check.py search \
  --hmm $R/hmms/mycosporine/ddgs.PROVISIONAL.hmm --genomes $T/genomes.corrected.tsv.zst \
  --chunk-manifest $R/manifests/chunk_manifest.size200.tsv.zst --prefix-aliases prefix_aliases.tsv \
  --group Rhodotorula --group Phaffia --group Cystobasidium --group Cystofilobasidium --group Rhodosporidiobolus \
  --group Sporobolomyces --group Naganishia --group Hortaea --group Pyricularia:100 \
  --out-dir results/search/mycosporine_rhodotorula_ddgs --cpus 4
python3.12 curation_scripts/mycosporine/rhodotorula_ddgs_check.py summarize --out-dir results/search/mycosporine_rhodotorula_ddgs --cutoff 463.5 --max-decoy 315.5
```

## Files

| File | Content |
|---|---|
| `per_genome.tsv` | 491 genomes: best ddgs domain (bit score, i-Evalue, HMM positions, coverage of the 412-position HMM, best protein id) |
| `bins.tsv` | per group: counts in score bins, genomes at or above max decoy 315.5 |
| `sensitivity.tsv` | per group: genomes with best score at or above the cutoff lowered by 0, 10, 25, 50 percent. Sensitivity only. |
| `selection_notes.txt` | HMM path, HMM length, number selected, overlap with `prefix_aliases.tsv` (none) |

Search report threshold: hmmsearch `-E 10` (`pathway_census.hits.search_step`). A genome with no reported hit would have an empty score. UNVERIFIED whether a lower threshold would report weaker hits.
Pyricularia: first 100 genomes with proteins, sorted by locustag (not random).
