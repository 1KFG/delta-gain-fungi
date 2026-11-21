# Pathway Census: Shared Tool and Carotenoid Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the shared `pathway_census` tool and use it to produce Chapter A (carotenoid): the step × phylum copy-number table, the paired deposited-versus-BFD comparison, and triage bins for functional testing.

**Architecture:** A Python package reads a pathway definition (YAML), retrieves candidates from BFD proteins with per-step HMMs, assigns each to a step by tree placement and domain checks, then builds a per-genome copy matrix with quality flags. Separate modules run the tblastn absence check, the paired comparison and the triage bins. Chapter A is a draft `.qmd` over the outputs.

**Tech Stack:** Python 3.12 (`/usr/bin/python3.12`), pandas, pyyaml, pytest, DuckDB CLI on the BFD parquet tables, HMMER, MAFFT, a tree-placement tool (reuse `analysis/enzymes_UPO/upo_place` if it fits), BLAST+ (tblastn), SLURM, Quarto.

**Spec:** `docs/superpowers/specs/2026-11-20-pathway-census-design.md` (approved 2026-11-20). Read it before starting. Review of the spec: `docs/superpowers/specs/2026-11-20-pathway-census-REVIEW_FABLE.md`.

**Out of scope here (separate plans):** (1) the species tree with `nf_phyling`, the one-per-genus list and the lineage-aware tests (spec §5.4 step 7, §5.7); (2) Chapter B, DHN-melanin and the L-DOPA family census. This plan produces the inputs both need.

## Global Constraints

- Run all Python with `/usr/bin/python3.12`. Never the default `python3` (3.9 on HPCC).
- No `BASH_SOURCE` in anything that runs under SLURM. Pass paths by environment variable or argument.
- SLURM: `-p epyc,intel --constraint="ryzen|broadwell|cascade|sapphire"`. Never `abu_dhabi` nodes (no AVX2).
- Write job temp files to `$SCRATCH` (use `${SCRATCH:?}`). Copy results to `/bigdata` before the job ends.
- Size jobs toward 1–1.5 h each.
- Compress tabular and sequence output: `zstd` for internal intermediates, `.gz` where a tool demands it. HMM files stay plain text.
- New scripts that read large text files accept plain, `.gz` and `.zst` input.
- Any step costlier than an alignment runs on dereplicated representatives, not the full family.
- Label every number as MEASURED, PROPOSED or UNVERIFIED in docs.
- A step call needs HMM plus tree-based orthology. A bare Pfam hit never counts as a pathway step.
- Pathway-specific clades only for completeness. SQS_PSY and polyprenyl synthase alone do not count.
- No activity, affinity or kinetics prediction. No fitted weights in triage.
- Genomes that are not downloaded get "unconfirmed" from the absence check, never "absent".
- All work lives in `analysis/pathway_census/`. Plans and specs stay in `docs/superpowers/`.

## Review Focus

Failure modes the spec implies that no task would otherwise test. Each has a test in the task that owns the code.

1. **Fused genes.** A protein carrying both phytoene synthase and lycopene cyclase domains must count once toward each step, and once in the raw gene count. (Task 5, Task 7.)
2. **Taxa with an empty GENUS.** The 8 "Ascomycota sp. Antarctolithica obscura" genomes have no GENUS value. They must appear in the genome table with a usable genus label, not drop out of any group-by. (Task 2.)
3. **Genomes with no proteins or no download.** `gap_noprotein_folders.tsv` lists 157 genomes (MEASURED: 158 lines including header) and `gap_notdownloaded.csv` lists 44. They must show as zero-protein or unconfirmed rows, not as steps absent. (Task 2, Task 8.)
4. **Same protein ID in both protein sets.** In the paired comparison, a deposited protein and a BFD protein can share an identifier. Counts must be kept per source. (Task 9.)
5. **A step with zero reference hits.** If an HMM finds nothing in BFD, the step row must exist with count 0, not vanish from the table. (Task 4, Task 7.)

---

## File Structure

```
analysis/pathway_census/
  README.md
  pyproject.toml                      # package metadata, pytest config
  pathways/carotenoid.yaml            # draft definition, verified: false until user review
  refs/carotenoid/                    # per-step reference FASTA, provenance.tsv
  hmms/carotenoid/                    # per-step HMMs (plain text)
  pathway_census/
    __init__.py
    defs.py                           # Pathway/Step/Ref dataclasses, load, validate
    genomes.py                        # genome table: taxonomy, status, quality
    hits.py                           # HMM build/search, domtbl parsing
    assign.py                         # step assignment, fused-gene labels
    quality.py                        # domain/motif/truncation, contamination, haplotype flags
    matrix.py                         # copy matrix, phylum summary
    absence.py                        # tblastn absence check
    paired.py                         # deposited-vs-BFD comparison
    bins.py                           # triage bins
  bin/                                # thin CLI wrappers (argparse)
  sql/                                # DuckDB queries
  slurm/                              # job scripts
  tests/                              # pytest + small fixtures in tests/data/
  chapters/A_carotenoid.qmd
```

---

### Task 1: Package scaffold and pathway definition loader

**Files:**
- Create: `analysis/pathway_census/pyproject.toml`, `analysis/pathway_census/pathway_census/__init__.py`, `analysis/pathway_census/pathway_census/defs.py`
- Test: `analysis/pathway_census/tests/test_defs.py`, fixtures `tests/data/pathway_ok.yaml`, `tests/data/pathway_bad_ref.yaml`

**Interfaces:**
- Produces:
  - `@dataclass Ref(accession: str, species: str, pmid: str | None, verified: bool)`
  - `@dataclass Step(id: str, name: str, ec: str | None, refs: list[Ref], domains: list[str], copy_expect: Literal["single","multi"], role: Literal["catalytic","auxiliary"], fused_with: list[str], optional: bool)`
  - `@dataclass Pathway(id: str, steps: list[Step], verified: bool)`
  - `load_pathway(path: Path) -> Pathway`
  - `validate(p: Pathway) -> list[str]` returns error strings, empty when valid.

- [ ] **Step 1: Write failing tests** in `test_defs.py`:
  - `test_load_ok_fixture` loads `pathway_ok.yaml` (two steps, one fused pair) and asserts `len(p.steps) == 2`, `p.steps[0].fused_with == ["crtY"]`.
  - `test_validate_ok_returns_empty`.
  - `test_ref_without_pmid_or_accession_is_error`: `pathway_bad_ref.yaml` has a ref with neither; `validate` returns one string containing the step id.
  - `test_verified_true_requires_all_refs_verified`: a pathway with `verified: true` and one ref `verified: false` yields an error.
  - `test_unknown_fused_with_target_is_error`.
- [ ] **Step 2: Run** `/usr/bin/python3.12 -m pytest analysis/pathway_census/tests/test_defs.py -v`. Expected: FAIL (module missing). If pytest or pyyaml is missing under 3.12, install into a project-local venv `analysis/pathway_census/.venv` and record that in `README.md`.
- [ ] **Step 3: Implement** `defs.py` per the interfaces. Use `yaml.safe_load`. Reject unknown keys with a clear error.
- [ ] **Step 4: Run tests.** Expected: all PASS.
- [ ] **Step 5: Commit** `feat(pathway_census): pathway definition loader and validator`.

---

### Task 2: Genome table

**Files:**
- Create: `pathway_census/genomes.py`, `sql/genomes.sql`
- Test: `tests/test_genomes.py`, fixture `tests/data/species_mini.parquet`, `tests/data/busco_mini.parquet`, `tests/data/prior_status_mini.tsv`

**Interfaces:**
- Consumes: BFD `species.parquet`, `busco_genome`, `asm_stats`, and `DeltaGain_Fungi/prior_annotation_status.tsv` (columns include `accession`, `status`, `annotation_provider`).
- Produces: `load_genome_table(tables_dir: Path, prior_status_tsv: Path) -> pd.DataFrame` with columns `locustag, species, phylum, class_, genus, genus_source, prior_status, annotation_provider, busco_complete, n50, downloaded, has_proteins`.
  - `genus_source` is `"genus_field"` or `"parsed_species"`.

- [ ] **Step 1: Discover schemas first.** Run `module load duckdb` then `duckdb -c "DESCRIBE select * from '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/species.parquet'"` and the same for `busco_genome`, `asm_stats`. Record the real column names (and the join key between `species.parquet` and the accession in `prior_annotation_status.tsv`) at the top of `genomes.py` as a comment. Do not guess names.
- [ ] **Step 2: Write failing tests.**
  - `test_empty_genus_is_parsed_from_species`: a row with GENUS empty and SPECIES "Ascomycota sp. Antarctolithica obscura" gets `genus == "Antarctolithica"`, `genus_source == "parsed_species"`.
  - `test_status_mapping`: rows map to the five status values counted in the spec (§3).
  - `test_not_downloaded_has_downloaded_false`.
  - `test_no_protein_row_has_has_proteins_false`.
  - `test_row_count_equals_species_rows`: no genome is dropped by the joins (left joins only).
- [ ] **Step 3: Run tests.** Expected: FAIL.
- [ ] **Step 4: Implement** `load_genome_table` using DuckDB via Python, left-joining everything onto `species`. For empty GENUS, parse the first capitalised token after "sp." as in the fixture. Log the count of parsed rows.
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Run on real data** and write the table to `results/genomes.tsv.zst`. Verify: row count equals the `species.parquet` row count, and the 8 Antarctolithica rows are present with `genus_source == "parsed_species"`. Record the counts per status in `README.md`, labelled MEASURED.
- [ ] **Step 7: Commit** `feat(pathway_census): genome table with taxonomy and quality`.

---

### Task 3: Carotenoid definition and reference curation (user review gate)

**Files:**
- Create: `pathways/carotenoid.yaml`, `refs/carotenoid/*.faa`, `refs/carotenoid/provenance.tsv`, `bin/list_unverified_refs.py`
- Test: `tests/test_carotenoid_def.py`

**Interfaces:**
- Consumes: `load_pathway`, `validate` from Task 1.
- Produces: a valid `carotenoid.yaml` with `verified: false`. `provenance.tsv` columns: `step_id, accession, species, pmid, source_database, retrieved_on, verified`.

- [ ] **Step 1: Write failing tests.**
  - `test_carotenoid_yaml_validates`: `validate(load_pathway("pathways/carotenoid.yaml")) == []`.
  - `test_step_ids_cover_spec_list`: the set of step ids includes GGPP synthase, phytoene synthase, phytoene desaturase, lycopene cyclase, carX, carT, astaxanthin synthase. (Spec §5.3, with carX and carT separate.)
  - `test_fused_psy_cyclase_declared`: the phytoene synthase step lists the lycopene cyclase step in `fused_with`, and vice versa.
  - `test_provenance_matches_yaml_refs`: every accession in the YAML appears in `provenance.tsv`.
- [ ] **Step 2: Run tests.** Expected: FAIL.
- [ ] **Step 3: Draft the definition and fetch references.** Use PubMed and UniProt/NCBI lookups for each reference. Mark `verified: true` only for a reference whose accession and PMID were each checked in the source. Add cao-2 and ylo-1 only if verified. Everything else stays `verified: false`.
- [ ] **Step 4: Run tests.** Expected: PASS.
- [ ] **Step 5: Run** `bin/list_unverified_refs.py pathways/carotenoid.yaml` and write its output to `refs/carotenoid/UNVERIFIED.md`.
- [ ] **Step 6: STOP for user review.** Present `carotenoid.yaml` and `UNVERIFIED.md`. The spec requires the user to review the step lists before any census runs (§5.8). Do not start Task 4 on real data until the user confirms. Task 4's tests use fixtures and may proceed.
- [ ] **Step 7: Commit** `feat(pathway_census): draft carotenoid definition and references`.

---

### Task 4: Candidate retrieval

**Files:**
- Create: `pathway_census/hits.py`, `bin/build_step_hmms.py`, `bin/search_step.py`, `slurm/search_step.sbatch`
- Test: `tests/test_hits.py`, fixtures `tests/data/tiny_refs.faa`, `tests/data/tiny_proteome.faa`, `tests/data/domtbl_sample.txt`

**Interfaces:**
- Produces:
  - `parse_domtbl(path: Path) -> pd.DataFrame` columns `protein_id, hmm_name, score, evalue, ali_from, ali_to, hmm_from, hmm_to`.
  - `build_hmm(ref_fasta: Path, out_hmm: Path) -> None` (MAFFT then `hmmbuild`).
  - `search_step(hmm: Path, proteins: Path, out_domtbl: Path, cpus: int) -> None` (wraps `hmmsearch --domtblout`).
  - `candidates_for_step(domtbl: pd.DataFrame, min_score: float) -> pd.DataFrame`; returns one row per protein, best domain score.

- [ ] **Step 1: Check tools.** Run `module avail hmmer mafft blast 2>&1 | head`. Record the exact module names and versions in `README.md`. If a tool is absent, stop and report.
- [ ] **Step 2: Write failing tests.**
  - `test_parse_domtbl_sample`: parsing `domtbl_sample.txt` yields the expected row count and the first protein's score.
  - `test_build_and_search_tiny`: build an HMM from `tiny_refs.faa`; searching `tiny_proteome.faa` finds the planted homologue and not the planted decoy.
  - `test_zero_hits_returns_empty_frame_with_columns`: searching an unrelated proteome returns an empty DataFrame that still has all columns.
  - `test_candidates_keeps_best_domain_per_protein`.
- [ ] **Step 3: Run tests.** Expected: FAIL.
- [ ] **Step 4: Implement** the four functions. `min_score` default is a parameter, not a constant. Do not choose its value here; the pilot sets it in Task 7 from the reference sets' own scores (record the rule in `README.md`).
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Build real step HMMs** from the curated refs (user-reviewed refs only, or `--provisional` flag clearly marking output). Write to `hmms/carotenoid/`.
- [ ] **Step 7: Search BFD.** Dereplicate large families first (Task 5 notes the size rule). For small families such as lycopene cyclase, search the whole BFD protein set directly with no Pfam prefilter (spec §5.4 step 2). Split into SLURM jobs of 1–1.5 h. Check the first job's real runtime, then size the rest.
- [ ] **Step 8: Commit** `feat(pathway_census): HMM build and candidate retrieval`.

---

### Task 5: Step assignment (tree placement and fused genes)

**Files:**
- Create: `pathway_census/assign.py`, `bin/assign_steps.py`
- Test: `tests/test_assign.py`, fixtures `tests/data/assign_hits.tsv`, `tests/data/assign_placements.tsv`

**Interfaces:**
- Consumes: `candidates_for_step` (Task 4); a placement result table with columns `protein_id, step_id, clade, support`.
- Produces:
  - `assign_steps(candidates: pd.DataFrame, placements: pd.DataFrame, pathway: Pathway, min_support: float) -> pd.DataFrame` columns `protein_id, step_id, support, label` where `label` is `"assigned"` or `"unassigned"`.
  - One row per (protein, step). A fused protein has two rows.

- [ ] **Step 1: Inspect `analysis/enzymes_UPO/upo_place/README.md` and `bin/`.** Decide whether `upo_place` handles step trees for these families. Write the decision, with evidence (what it takes as input, what it outputs), in `README.md`. If it does not fit, state the replacement (for example an `epa-ng` or `pplacer` run on a per-step reference tree) and why. Do not assume it scales to PKS-size families; the carotenoid families here are smaller, but record the largest family's real size.
- [ ] **Step 2: Write failing tests.**
  - `test_fused_protein_gets_two_step_rows`: a protein with domain hits to both phytoene synthase and cyclase HMMs yields two `assigned` rows.
  - `test_outside_every_clade_is_unassigned_not_dropped`.
  - `test_low_support_is_unassigned`.
  - `test_assigned_only_to_declared_steps`: a placement naming a step not in the pathway raises `ValueError`.
- [ ] **Step 3: Run tests.** Expected: FAIL.
- [ ] **Step 4: Implement** `assign_steps`. Fusion is detected from domain hits (Task 4 output) against `Step.fused_with`, not from tree placement alone.
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Run on real candidates.** Report counts per step, how many `unassigned`, and how many fused. Label MEASURED.
- [ ] **Step 7: Commit** `feat(pathway_census): step assignment with fused-gene labels`.

---

### Task 6: Quality, contamination and haplotype flags

**Files:**
- Create: `pathway_census/quality.py`
- Test: `tests/test_quality.py`, fixtures under `tests/data/`

**Interfaces:**
- Produces: `quality_flags(assigned: pd.DataFrame, domains: pd.DataFrame, gene_info: pd.DataFrame, refs_len: dict[str, tuple[int,int]], genomes: pd.DataFrame) -> pd.DataFrame` columns `protein_id, step_id, domain_ok, motif_ok, truncated, contam_suspect, hap_dup_suspect`.

- [ ] **Step 1: Write failing tests.**
  - `test_truncated_when_below_length_floor`: length below the reference range lower bound times a parameter sets `truncated`.
  - `test_missing_expected_domain_sets_domain_ok_false`.
  - `test_contam_suspect_when_neighbours_have_non_fungal_best_hits`: fixture with a scaffold whose other genes' taxonomic assignment is bacterial.
  - `test_hap_dup_suspect_when_pair_is_near_identical_and_adjacent_scaffolds_duplicate`: pins the rule as written in step 3.
- [ ] **Step 2: Run tests.** Expected: FAIL.
- [ ] **Step 3: Define the rules** in the module docstring before coding: how contamination is judged (the data source for taxonomic placement; if BFD has no per-gene taxonomy table, use best-hit against a reference set of fungal vs bacterial/algal carotenoid proteins plus the scaffold's gene neighbourhood), and how haplotype duplication is judged. State thresholds as parameters. If a data source is missing, say so in `README.md` and implement the flag only for what the data supports.
- [ ] **Step 4: Implement.** Domain-level counting uses `domains`, not whole-protein counts.
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Commit** `feat(pathway_census): quality, contamination, haplotype flags`.

---

### Task 7: Copy matrix and phylum summary

**Files:**
- Create: `pathway_census/matrix.py`, `bin/build_matrix.py`
- Test: `tests/test_matrix.py`

**Interfaces:**
- Consumes: outputs of Tasks 2, 5, 6.
- Produces:
  - `build_matrix(assigned: pd.DataFrame, flags: pd.DataFrame, genomes: pd.DataFrame, pathway: Pathway) -> pd.DataFrame` long format, one row per (genome, step), columns `locustag, step_id, copy_raw, copy_filtered, has_proteins, genome_flags`.
  - `phylum_summary(matrix: pd.DataFrame, genomes: pd.DataFrame) -> pd.DataFrame` columns `phylum, step_id, n_genomes, n_with_step, median_copy, max_copy, frac_with_step`.

- [ ] **Step 1: Write failing tests.**
  - `test_every_genome_step_pair_present`: genome count × step count rows, including zeros (Review Focus 5).
  - `test_zero_protein_genome_marked_not_zero_evidence`: `has_proteins` is false and the row is excluded from `n_genomes` in the summary.
  - `test_fused_protein_counts_once_per_step`: one fused protein gives `copy_raw == 1` for each of its two steps (Review Focus 1).
  - `test_filtered_excludes_truncated_and_contam`.
  - `test_ubiquitous_family_alone_does_not_count`: a protein assigned only by bare Pfam, with no pathway-specific clade, is not counted (spec anti-goal).
- [ ] **Step 2: Run tests.** Expected: FAIL.
- [ ] **Step 3: Implement.** Use pandas merges against the full genome × step cross product.
- [ ] **Step 4: Run tests.** Expected: PASS.
- [ ] **Step 5: Run on real data.** Write `results/carotenoid_matrix.tsv.zst` and `results/carotenoid_phylum_summary.tsv`. Check: genome rows equal the Task 2 count; no phylum is missing. Report the headline table in plain numbers, MEASURED.
- [ ] **Step 6: Commit** `feat(pathway_census): copy matrix and phylum summary`.

---

### Task 8: Absence check

**Files:**
- Create: `pathway_census/absence.py`, `bin/absence_check.py`, `slurm/absence_check.sbatch`
- Test: `tests/test_absence.py`, fixtures `tests/data/mini_genome.fna`, `tests/data/mini_step_refs.faa`

**Interfaces:**
- Produces:
  - `select_absences(matrix: pd.DataFrame, genomes: pd.DataFrame, relatives: dict[str, list[str]]) -> pd.DataFrame` columns `locustag, step_id` for genomes lacking a step that their listed relatives carry.
  - `check_absence(genome_fna: Path | None, step_refs: Path, min_bitscore: float) -> str` returns `"present_unannotated" | "fragmented" | "absent" | "unconfirmed"`.

- [ ] **Step 1: Write failing tests.**
  - `test_genome_not_downloaded_is_unconfirmed`: `genome_fna=None` returns `"unconfirmed"`.
  - `test_planted_gene_found`: a planted ORF in `mini_genome.fna` returns `"present_unannotated"`.
  - `test_no_homology_returns_absent`.
  - `test_fragmented_when_hit_spans_less_than_fraction_of_query`: pins the rule (parameter `min_query_cover`).
- [ ] **Step 2: Run tests.** Expected: FAIL.
- [ ] **Step 3: Implement** with `tblastn`. The `relatives` mapping is an input, not computed here, because the species tree is a separate plan; for the pilot use same-genus genomes from Task 2. State this in the docstring.
- [ ] **Step 4: Run tests.** Expected: PASS.
- [ ] **Step 5: Run on real absences.** Genome FASTA locations: check `seqdb/genomes.tsv` and the NCBI_ASM folders first, and record which source each genome uses. Report counts per outcome, MEASURED.
- [ ] **Step 6: Commit** `feat(pathway_census): tblastn absence check`.

---

### Task 9: Paired deposited-versus-BFD comparison

**Files:**
- Create: `pathway_census/paired.py`, `bin/paired_compare.py`
- Test: `tests/test_paired.py`

**Interfaces:**
- Consumes: Task 7 `build_matrix`, run once per protein source.
- Produces:
  - `paired_table(matrix_deposited: pd.DataFrame, matrix_bfd: pd.DataFrame, genomes: pd.DataFrame) -> pd.DataFrame` columns `locustag, step_id, copy_deposited, copy_bfd, delta, annotation_provider, busco_complete`.
  - `summarize_paired(table: pd.DataFrame) -> pd.DataFrame` columns `step_id, annotation_provider, n_pairs, mean_delta, ci_low, ci_high, frac_bfd_more, frac_bfd_less`.

- [ ] **Step 1: Count the pairs first.** Query how many `baseline_ok` genomes have both deposited proteins (`NCBI_ASM/<acc>/protein.faa` or provider equivalent) and BFD proteins. Record the count, MEASURED. (Fable counted 5,625; this has not been re-run.)
- [ ] **Step 2: Write failing tests.**
  - `test_same_protein_id_in_both_sources_not_merged`: identifiers shared across sources stay separate (Review Focus 4).
  - `test_delta_sign`: BFD 3 versus deposited 2 gives `delta == 1`.
  - `test_pair_requires_both_sources`: a genome in one source only is excluded and counted in a returned `n_unpaired`.
  - `test_summary_by_provider_drops_providers_below_min_pairs`: parameter `min_pairs`; dropped providers are listed, not silently removed.
  - `test_ci_uses_genome_resampling`: the interval is a bootstrap over genomes, not over proteins (fixed seed in the test).
- [ ] **Step 3: Run tests.** Expected: FAIL.
- [ ] **Step 4: Implement.** Run the same pipeline (Tasks 4–7) on the deposited proteins and on BFD proteins for each paired genome. Keep the HMMs and tree placement identical for both.
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Run on real data.** Report per-step mean delta with interval, overall and by provider. A difference counts as a finding only if it is consistent across providers and survives the BUSCO control (spec §5.6). Otherwise say so plainly.
- [ ] **Step 7: Commit** `feat(pathway_census): paired deposited-vs-BFD comparison`.

---

### Task 10: Triage bins

**Files:**
- Create: `pathway_census/bins.py`, `bin/triage.py`
- Test: `tests/test_bins.py`

**Interfaces:**
- Consumes: Task 6 flags, Task 5 assignments, `gene_info`, `tmhmm`, `signalp`, `wolfpsort`.
- Produces: `triage(copies: pd.DataFrame, window_bp: int) -> pd.DataFrame` columns `protein_id, step_id, bin, reasons` where `bin` is one of `"tier1" | "tier2" | "tier3" | "excluded"` and `reasons` is a semicolon-separated string.

- [ ] **Step 1: Fix the bin rules in writing before looking at results.** Put them in the module docstring: which flags exclude a copy, which terms promote it (divergence from nearest reference, lineage-specific expansion, a pathway neighbour within `window_bp` on the same scaffold, expressibility). No fitted weights. Commit the docstring alone first.
- [ ] **Step 2: Write failing tests.**
  - `test_truncated_copy_is_excluded`.
  - `test_neighbour_within_window_promotes`.
  - `test_neighbour_on_other_scaffold_does_not_promote`.
  - `test_reasons_lists_every_term_that_fired`.
  - `test_same_input_same_bin`: determinism.
- [ ] **Step 3: Run tests.** Expected: FAIL.
- [ ] **Step 4: Implement.** Pick `window_bp` and state why in the docstring (for example, from the known cluster span in the reference species; if that is unknown, say so and use a stated provisional value).
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Run on real data.** Write `results/carotenoid_triage.tsv`. Report the count per bin.
- [ ] **Step 7: Commit** `feat(pathway_census): triage bins`.

---

### Task 11: Chapter A draft and pilot review

**Files:**
- Create: `chapters/A_carotenoid.qmd`, `docs/superpowers/plans/2026-11-20-pilot-lessons.md`
- Test: render check

**Interfaces:**
- Consumes: `results/carotenoid_phylum_summary.tsv`, `carotenoid_paired_summary.tsv`, `carotenoid_triage.tsv`, absence-check summary.

- [ ] **Step 1: Write the draft** `A_carotenoid.qmd` in the style of `reports/mucor_secretome.qmd`. Sections: pathway steps table (with reference provenance), phylum × step copy-number table, paired-comparison results, absence-check summary, triage list, limitations. Every number loads from a results file. The species-tree sections get a placeholder that says the tree plan is separate.
- [ ] **Step 2: Render** with `quarto render analysis/pathway_census/chapters/A_carotenoid.qmd`. Expected: the file renders with no error and the tables are non-empty.
- [ ] **Step 3: Write the pilot lessons note.** What worked, what the data forced, which spec assumptions failed (for example orthology ambiguity, fused-gene counting, provider effects). The DHN-melanin plan starts from this note (spec §5.8).
- [ ] **Step 4: Commit** `docs(pathway_census): carotenoid chapter draft and pilot lessons`.

---

## Self-review

- **Spec coverage:** §5.4 steps 1–6, 8 (paired) and 9 map to Tasks 3, 4, 5, 6, 7, 8, 9, 10. Step 10 maps to Task 11. Steps 7 and the tree (§5.7) are deliberately in a separate plan. The Antarctolithica rule is covered in Task 2 at the genus-label level; the BUSCO inclusion rule applies to tree tips and belongs to the tree plan. Contamination, haplotype and domain-level counting are in Tasks 5–6. The "unconfirmed" label is in Task 8.
- **Open spec items:** the `phyling` marker set is in the tree plan. The paired-genome count is Task 9 step 1. The Antarctolithica BUSCO values show up in the Task 2 output.
- **Type consistency:** `Pathway`, `Step`, `Ref` (Task 1) are the types used in Tasks 3, 5, 7. `build_matrix` (Task 7) feeds Tasks 8 and 9. Column names match across tasks.
- **Known gap:** Task 6 depends on a data source for contamination that may not exist in BFD. The task says what to do if it does not.
