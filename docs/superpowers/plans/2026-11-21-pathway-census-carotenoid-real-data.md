# Pathway Census: Carotenoid Real-Data Pass Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the carotenoid pathway census on real BFD data end to end, with the paired deposited-versus-BFD comparison, the absence check and the triage bins, and render Chapter A with real numbers.

**Architecture:** Fix the definition inputs first (per-domain references, decoys, any-of domains), then run the stages the pilot built: batched HMM search, tree placement, step assignment, quality flags, copy matrix. Then add the deposited arm and the paired test, the absence check, the triage, and the chapter. Every real run is stamped PROVISIONAL until the user flips `verified: true` (R10).

**Tech Stack:** Python 3.12 (`/usr/bin/python3.12`), the `pathway_census` package from the pilot, hmmer 3.4, mafft 7.505, mmseqs2, ncbi-blast 2.17.0+, `epa-ng` via `pixi exec` (as `upo_place` does; no module exists), SLURM, DuckDB, Quarto 1.7.31.

**Specs and prior work:**
- Spec (approved): `docs/superpowers/specs/2026-11-20-pathway-census-design.md`.
- Pilot lessons and the 13-step checklist this plan implements: `docs/superpowers/plans/2026-11-20-pilot-lessons.md` (section 4).
- Pilot plan: `docs/superpowers/plans/2026-11-20-pathway-census-carotenoid-pilot.md`.
- Fable review of the spec: `docs/superpowers/specs/2026-11-20-pathway-census-REVIEW_FABLE.md`.

**Out of scope here (separate plans):** the species tree with `nf_phyling`; lineage-aware tests; DHN-melanin and L-DOPA (Chapter B).

## Decisions made for this plan (the user can reverse any of them)

- **Gate (user, 2026-11-21: "move ahead with this current plan and specs"):** the current `carotenoid.yaml` step list is accepted as reviewed. D9IVD8 is not added: it is a 518 aa fragment of the same gene as B5LEC8 (512 of 518 aligned residues identical, MEASURED 2026-11-21). `carS` (Mucorales carotene dioxygenase, UniProt G8Z9M3) is a different enzyme and stays out.
- **Provisional until R10.** The pathway has one unverified reference (PMID 16416328, no accession). All outputs from R2 to R9 carry the `.PROVISIONAL` stamp. `verified: true` is a user decision at R10.
- **Domain rule: any-of.** A step's `domains` become `domains_any` (a copy passes if it has at least one listed domain). PROPOSED; reason: `carotenoid.yaml` lists alternatives such as `[TIGR03462, PF18916]` for lycopene cyclase.
- **GGPP synthase stays out of completeness** until a carotenoid-specific clade is defined (PROPOSED).
- **Matrix and triage `require`:** keep the pilot defaults. Add raw-count outputs next to the filtered ones. With no contamination table, `copy_filtered` is 0 everywhere (MEASURED by the final reviewer on a probe); the chapter must show raw and filtered side by side.
- **A copy is a gene.** The deposited arm collapses isoforms to one protein per gene.

## Global Constraints

- Run all Python with `/usr/bin/python3.12`. Never the default `python3`.
- No `BASH_SOURCE` in anything that runs under SLURM.
- SLURM: `-p epyc,intel --constraint="ryzen|broadwell|cascade|sapphire"`. Never `abu_dhabi` nodes.
- Temp files in `${SCRATCH:?}`; copy results to `/bigdata` before the job ends.
- Size jobs toward 1-1.5 h of real runtime. Time one real job before sizing the rest. Do not submit a job per genome per step.
- Compress large tabular output with `zstd` (internal) or `.gz` (where a tool needs it). HMM files stay plain text. Do not commit data or results.
- Do not run more than one SLURM array at a time without checking the queue (`squeue -u $USER`) first.
- A step call needs HMM evidence AND a tree placement. A bare Pfam or HMM hit never counts.
- Never treat NA/unknown as clean. Never `fillna(False)` on flags. Use `clean_mask` and `unknown_counts`.
- Genomes with `has_proteins == False` or `downloaded == False` are never read as absences.
- Label every number MEASURED, PROPOSED or UNVERIFIED. Do not write any biological finding into a document before a result file supports it.
- `genus_source == "unresolved"` labels are never groups.
- Keep all work in `analysis/pathway_census/`; plans in `docs/superpowers/`. Do not commit `.superpowers/`.

## Review Focus

1. **A decoy that places inside the target clade.** A P450 or FPP-synthase decoy must place outside the step's clade, or the step gets inflated. (R1, R4.)
2. **A search that silently covers fewer genomes than it should.** The manifest must list every genome searched; a missing genome must not read as a zero. (R2, R5.)
3. **Deposited and BFD genes counted once.** The same gene must not be counted twice in the deposited arm through isoforms. (R6.)
4. **PROVISIONAL leaking into a chapter.** A provisional file must never render as a result. (R10.)
5. **Genome-id mapping for deposited proteins.** A wrong mapping would pair the wrong genomes. (R6.)

---

### Task R1: Definition inputs (references, decoys, domain rule)

**Files:**
- Modify: `analysis/pathway_census/pathway_census/defs.py` (schema: `domains_any`), `analysis/pathway_census/pathway_census/quality.py` (any-of in `domain_ok`)
- Create: `analysis/pathway_census/refs/carotenoid/decoys/<step>.faa`, `refs/carotenoid/decoys/provenance.tsv`; per-domain split files under `refs/carotenoid/` (names below)
- Modify: `analysis/pathway_census/pathways/carotenoid.yaml`
- Test: `tests/test_defs.py`, `tests/test_quality.py`, new `tests/test_carotenoid_inputs.py`

**Interfaces:**
- Produces: `Step.domains_any: list[str]` (replaces `domains`; `load_pathway` accepts the old key `domains` as an error with a clear message); `quality.quality_flags(...)` treats `domains_any` as any-of; for fused steps the YAML refs gain a `region: [start, end]` (1-based, inclusive) so per-domain FASTA can be cut; decoy FASTA per step with a `provenance.tsv` (columns `step_id, accession, species, protein_name, reason, source_database, retrieved_on`).

- [ ] **Step 1: Write failing tests.**
  - `test_domains_any_loads_and_old_key_errors`.
  - `test_domain_ok_is_any_of`: a copy with one of two listed domains is `True`; with none is `False`; a step with no listed domains is `NA`.
  - `test_fused_refs_have_regions_and_split_files_differ`: for each fused ref the region lies inside the sequence; `phytoene_synthase.faa` and `lycopene_cyclase.faa` are not byte-identical and each ref sequence equals the cut region of the fused UniProt sequence (use the stored full sequences).
  - `test_decoys_have_provenance_and_do_not_overlap_refs`: every decoy accession is in `provenance.tsv`; no decoy accession is also a reference accession of the same step.
- [ ] **Step 2: Run the tests; expect FAIL.**
- [ ] **Step 3: Implement the schema change in `defs.py` and the any-of rule in `quality.py`.** Update the carotenoid YAML to `domains_any`.
- [ ] **Step 4: Split the fused references.** From the 5 fused references (see `pathways/carotenoid.yaml`), cut the phytoene synthase region and the lycopene cyclase region using UniProt domain annotation (Pfam PF00494 for phytoene synthase, PF18916 / TIGR03462 for cyclase; check each entry on UniProt, record the coordinates in the YAML with source). Do not guess coordinates; if an entry has no usable annotation, say so in the report and drop that ref from the split files only.
- [ ] **Step 5: Curate decoys, per step, from UniProt with provenance (source_database, retrieved_on).** For each step choose 3-8 close non-pathway homologues: GGPP synthase -> farnesyl/geranyl diphosphate synthases from fungi; phytoene synthase -> squalene synthase; phytoene desaturase -> other FAD amine oxidases; lycopene cyclase -> non-carotenoid homologues if any exist, else state none; carX/carT -> other carotenoid cleavage oxygenase family members (note carX vs carT are each other's decoys); ylo1 -> other aldehyde dehydrogenases; astaxanthin synthase -> fungal cytochrome P450s (e.g. from X. dendrorhous and other fungi that are not carotenoid hydroxylase/ketolases). Verify each accession on UniProt before listing it. Never invent an accession.
- [ ] **Step 6: Run the tests; expect PASS.** Run the full suite `/usr/bin/python3.12 -m pytest analysis/pathway_census/tests -q`.
- [ ] **Step 7: Commit** `feat(pathway_census): per-domain references, decoys and any-of domain rule`.

---

### Task R2: Batched search driver with a searched-genome manifest

**Files:**
- Create: `analysis/pathway_census/pathway_census/driver.py`, `analysis/pathway_census/bin/search_genomes.py`, `analysis/pathway_census/slurm/search_genomes.sbatch`
- Modify: `analysis/pathway_census/bin/search_step.py` (reuse, not duplicate), `analysis/pathway_census/pathway_census/matrix.py` (manifest check)
- Test: `tests/test_driver.py`

**Interfaces:**
- Produces:
  - `plan_chunks(genomes: pd.DataFrame, genomes_per_chunk: int) -> list[list[str]]`: locustags with `has_proteins` True only, sorted, chunked.
  - `search_chunk(locustags: list[str], hmm_dir: Path, seqdb_root: Path, out_dir: Path, chunk_id: str, min_scores: dict[str, float], cpus: int) -> Path`: for each genome and each step HMM runs `hmmsearch`, writes `candidates.chunk<ID>.tsv.zst` (columns `hits.COLUMNS + ["step_id", "locustag"]`) and `manifest.chunk<ID>.tsv` (columns `locustag, status, n_proteins, n_hits`; status in `searched | protein_file_missing | error`).
  - `matrix.build_matrix(..., manifest: pd.DataFrame | None)`: when a manifest is given, a genome with `has_proteins` True but no `searched` row gets `genome_flags` token `not_searched` and is excluded from `n_genomes`. A genome never silently becomes a zero.
  - Protein files are located by the same logic `genomes.scaffolds_path` uses (sibling `proteins` path); read the real convention from `genomes.py`, do not guess.

- [ ] **Step 1: Failing tests** (fixtures: 3 tiny proteomes under a fake `seqdb`): `test_plan_chunks_skips_no_protein_genomes_and_is_sorted`; `test_search_chunk_writes_manifest_row_per_genome`; `test_missing_protein_file_gets_status_not_zero` (the genome appears in the manifest as `protein_file_missing` and in no candidate row); `test_candidates_carry_step_id_and_locustag`; `test_build_matrix_marks_not_searched` (genome with proteins and no manifest row has `not_searched` and is excluded from `n_genomes`); `test_chunk_outputs_do_not_overwrite` (two chunks in one directory).
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement** `driver.py` and the CLI; the sbatch passes `--chunk-id "${SLURM_ARRAY_TASK_ID:-}"`, writes to `${SCRATCH:?}`, and copies to `/bigdata`. The CLI refuses a `verified: false` pathway without `--provisional` (reuse `hits.provisional_tag`).
- [ ] **Step 4: Pilot timing.** Run ONE job on 200 genomes (random sample, fixed seed, recorded). Record wall time and peak memory (`sacct`). Compute `genomes_per_chunk` for a 1-1.5 h target and write it, with the measured numbers, in `analysis/pathway_census/README.md` (MEASURED).
- [ ] **Step 5: Run tests; expect PASS. Commit** `feat(pathway_census): batched search driver and searched-genome manifest`.
- [ ] **Step 6 (execution, user-visible):** report the chunk count and expected total node-hours, then submit the full array in a single `sbatch --array` with a concurrency cap (for example `%20`). Wait for completion; check the manifest: every `has_proteins` genome has a status; list non-`searched` ones. Do not proceed to R3 until the count of unsearched genomes is reported.

---

### Task R3: Placement wiring

**Files:**
- Create: `analysis/pathway_census/pathway_census/placement.py`, `analysis/pathway_census/bin/place_steps.py`, `analysis/pathway_census/slurm/place_steps.sbatch`
- Test: `tests/test_placement.py`

**Interfaces:**
- Produces:
  - `build_step_tree(ref_faa: Path, decoy_faa: Path, hmm: Path, out_dir: Path) -> TreeBundle` (dataclass with paths: alignment, tree, model info, clade table). Reference sequences plus decoys go into one tree. Decoy tips are labelled `decoy`; reference tips carry their step id.
  - `place_candidates(bundle: TreeBundle, candidates_faa: Path, out_jplace: Path, cpus: int) -> None` (epa-ng via `pixi exec`).
  - `placements_table(jplace: Path, bundle: TreeBundle, min_decoy_lwr: float) -> pd.DataFrame` with columns `protein_id, step_id, clade, support` (the Task 5 interface), where `clade` is the step id when the placement mass inside the step's reference clade is highest, and `decoy` when the decoy clade wins. Rows whose best clade is `decoy` are written with `clade="decoy"` so `assign_steps` leaves them `unassigned`.

- [ ] **Step 1: Check the tool.** Run `pixi exec epa-ng --version` as `upo_place/README.md` does; record it. If it fails, say so and stop with BLOCKED (do not substitute another tool without asking).
- [ ] **Step 2: Failing tests** with a small synthetic tree (real mafft, iqtree or fasttree, real epa-ng; skip with a clear reason only if the tool is missing): `test_decoy_sequence_places_in_decoy_clade`; `test_reference_like_sequence_places_in_step_clade`; `test_clade_table_has_one_row_per_protein_step`; `test_decoy_clade_rows_become_unassigned_in_assign_steps` (feed the table into `assign_steps`).
- [ ] **Step 3: Run; expect FAIL. Implement. Run; expect PASS.**
- [ ] **Step 4 (execution):** build one tree per step; place all dereplicated candidates; write `placements.<step>.tsv.zst`. Report per step: candidates, placed in step clade, placed in decoy clade, low-support.
- [ ] **Step 5: Commit** `feat(pathway_census): per-step tree placement with decoys`.

---

### Task R4: Step assignment run, contamination and haplotype tables

**Files:**
- Create: `analysis/pathway_census/pathway_census/evidence.py`, `analysis/pathway_census/bin/build_evidence.py`
- Test: `tests/test_evidence.py`

**Interfaces:**
- Produces:
  - `best_hit_table(assigned: pd.DataFrame, ref_fungal_faa: Path, ref_nonfungal_faa: Path, work_dir: Path, cpus: int) -> pd.DataFrame` columns `protein_id, taxon_class` (`fungal` | `non_fungal`), from a real `mmseqs` or BLAST search of each assigned protein against a labelled reference set; the reference set is the step references (fungal) plus bacterial and algal carotenoid proteins curated with provenance in `refs/carotenoid/nonfungal.faa` and `nonfungal_provenance.tsv`.
  - `pair_identity_table(assigned_proteins_faa: Path, gene_info: pd.DataFrame, work_dir: Path) -> pd.DataFrame` columns the `quality_flags` function expects (read its docstring; do not guess).
- Consumes: R2 candidates, R3 placements, `assign.assign_steps`.

- [ ] **Step 1:** run `assign_steps.py` on the real candidates and placements (`--provisional`); report counts per step: assigned, unassigned, fusion candidates.
- [ ] **Step 2:** curate `nonfungal.faa` (carotenoid enzymes from bacteria and algae for the same steps, UniProt, verified accessions, provenance). Never invent accessions.
- [ ] **Step 3: Failing tests** with synthetic sequences and tiny reference sets: `test_best_hit_labels_nonfungal_planted_sequence`; `test_best_hit_unknown_when_no_hit` (no hit gives no row, so `quality_flags` gets NA); `test_pair_identity_flags_near_identical_pair`.
- [ ] **Step 4: Implement, run tests, run on real assigned proteins. Commit** `feat(pathway_census): contamination and haplotype evidence tables`.

---

### Task R5: Quality flags, copy matrix and phylum summary (BFD arm)

**Files:**
- Create: `analysis/pathway_census/bin/run_quality.py` (thin glue around `quality.quality_flags`)
- Modify: `analysis/pathway_census/bin/build_matrix.py` (manifest input; raw summary output)
- Test: `tests/test_run_quality.py`

- [ ] **Step 1: Failing tests:** CLI on a synthetic fixture produces the flags table with nullable-boolean columns as `parse_bool_col` expects; `copy_filtered` is 0 and `copy_unknown` equals `copy_raw` when the evidence tables are absent (documented behaviour); a genome with no `searched` manifest row is `not_searched` and excluded from `n_genomes`.
- [ ] **Step 2: Run; expect FAIL. Implement. Run; expect PASS.**
- [ ] **Step 3 (execution):** `quality_flags` on the R4 evidence; `build_matrix` and the phylum summary on the BFD proteins. Write `results/carotenoid_matrix.PROVISIONAL.tsv.zst` and `results/carotenoid_phylum_summary.PROVISIONAL.tsv`. Report: genomes in the table, genomes searched, phyla, and per step genomes with the step (raw and filtered), MEASURED.
- [ ] **Step 4: Commit** `feat(pathway_census): quality glue and BFD copy matrix`.

---

### Task R6: Deposited arm

**Files:**
- Create: `analysis/pathway_census/pathway_census/deposited.py`, `analysis/pathway_census/bin/prepare_deposited.py`
- Modify: `analysis/pathway_census/pathway_census/matrix.py` (a `genome_of` mapping parameter)
- Test: `tests/test_deposited.py`

**Interfaces:**
- Produces:
  - `pair_genomes(genomes: pd.DataFrame, deposited_root: Path) -> pd.DataFrame`: the pair table (`paired.count_pairs` rules; folder path `<root>/<folder>/<folder>_protein.faa.gz`); returns the real pair count.
  - `collapse_isoforms(protein_faa: Path, gff: Path, out_faa: Path, locustag: str) -> dict`: one protein per gene (the longest isoform; ties by ID), rewrites each ID to `<locustag>_<gene>` so the matrix's LOCUSTAG-prefix rule holds; returns counts `{genes, proteins_in, proteins_out}`.
  - `build_matrix(..., genome_of: Callable[[str], str] | None)` as an alternative route for proteins that keep their own IDs.
- The pairing must read the deposited GFF from the same folder (check the real file name in the folder; MEASURED layout: `<folder>/<folder>_protein.faa.gz`; find the GFF name and record it).

- [ ] **Step 1 (count first):** regenerate `results/genomes.tsv.zst` (it needs `asmid` and `folder`); run `pair_genomes` and report the real number of pairs by `annotation_provider` (MEASURED). Fable counted 5,625 baseline genomes, not pairs.
- [ ] **Step 2: Failing tests:** `test_collapse_isoforms_keeps_longest_one_per_gene` (gene with 3 mRNAs gives 1 protein); `test_rewritten_ids_start_with_locustag`; `test_isoform_counts_returned`; `test_pair_count_matches_fixture`.
- [ ] **Step 3: Implement. Run tests; expect PASS.**
- [ ] **Step 4 (execution):** prepare the deposited proteomes for the pair genomes into `/bigdata` (compressed); run R2-R5 on them with the same HMMs, trees and decoys; write `results/carotenoid_matrix_deposited.PROVISIONAL.tsv.zst`. Deposited flags are NA for contamination (no coordinates), so only raw `delta` is usable; state this in the output header.
- [ ] **Step 5: Commit** `feat(pathway_census): deposited arm with isoform collapse`.

---

### Task R7: Paired comparison run

- [ ] **Step 1:** run `bin/paired_compare.py` on the R5 and R6 matrices (`--provisional`). Report per step: n pairs, mean delta (raw), 95% CI, per annotation provider, the dropped providers, and the BUSCO-bin strata.
- [ ] **Step 2: Write `analysis/pathway_census/results/PAIRED_NOTES.md`** (not committed if results-only; commit only if it contains no data tables) with the plain statement of what the paired test shows and does not show (spec §5.6). A difference is a finding only if consistent across providers and BUSCO bins. Otherwise write "not a finding".
- [ ] **Step 3: Check:** assert `n_pairs` in the summary equals the pair count from R6; assert no provider named `ALL` or `(missing)` exists in the input.

---

### Task R8: Absence check run

**Files:**
- Modify: `analysis/pathway_census/bin/absence_check.py` (per-step `--min-bitscores` from decoy scores)
- Test: `tests/test_absence.py` (new case)

- [ ] **Step 1:** for each step, run `tblastn` of the decoys against a sample of BFD genomes that carry the step (fixed seed, recorded) and of the references against the same; set `min_bitscore` to a value that separates decoys from references where a gap exists; record the decoy and reference score distributions (MEASURED) in the README. Where no gap exists, say so and use the reference-only value with the label UNVERIFIED.
- [ ] **Step 2: Failing test:** `test_per_step_min_bitscore_is_used` (two steps with different thresholds give different outcomes on one genome).
- [ ] **Step 3: Implement. Run tests. Run `select_absences` with `same_genus_relatives`, then `absence_check` in chunks (one timed job first); merge with `bin/merge_absence.py`.** Report outcome counts: present_unannotated, fragmented, absent, unconfirmed.
- [ ] **Step 4: Commit** `feat(pathway_census): per-step absence thresholds`.

---

### Task R9: Triage glue and run

**Files:**
- Create: `analysis/pathway_census/bin/triage_inputs.py`
- Test: `tests/test_triage_inputs.py`

**Interfaces:** produces the joined frame `bins.triage` expects: `signal_peptide`, `tm_helices`, `ref_identity`, `length`, `copy_filtered`, coordinates, `step_id`, `protein_id`, `genome`.

- [ ] **Step 1:** read the real column names of the BFD `tmhmm`, `signalp`, `wolfpsort` tables and `gene_info` (`duckdb -c "DESCRIBE ..."`); record them in the README (MEASURED). `signalp` holds positives only; absence means no signal peptide predicted (documented assumption).
- [ ] **Step 2: Failing tests:** each join keeps one row per `(protein_id, step_id)`; missing signalp row gives `signal_peptide=False`; missing tmhmm row gives NA `tm_helices` (not 0); `ref_identity` is the best identity of the protein to a step reference.
- [ ] **Step 3: Implement. Run tests. Run triage (`--provisional`)** with an explicit `--soluble-steps` list decided from the biology of each step (state the list and reason; if unsure leave empty and say so).
- [ ] **Step 4: Report** bin counts and the `unknown` share. Commit `feat(pathway_census): triage input glue`.

---

### Task R10: Chapter A with real numbers; user decision on `verified`

**Files:**
- Modify: `analysis/pathway_census/chapters/A_carotenoid.qmd`, `analysis/pathway_census/pathways/carotenoid.yaml` (only after the user's decision)
- Test: `tests/test_chapter_render.py`

- [ ] **Step 1: Resolve the one unverified reference.** Check UniProt entries A0A0A8KXB4 and A0A0C4MW73 and D9IVD8 for citations with a functional assay of β-carotene to astaxanthin; and the Ojima 2006 paper (PMID 16416328) for a sequence accession. Record what you find. Do not change `verified` flags without a source.
- [ ] **Step 2: Ask the user** (this is a stop point): "Set `verified: true` for `carotenoid.yaml`?" Show what is still unverified and what the results look like. If yes, change the YAML in a separate commit titled `feat: mark carotenoid definition verified (user decision <date>)`, re-run the chapter-feeding stages without `--provisional` only if the user wants it (outputs then lose the stamp). If no, the chapter stays a scaffold and the provisional files are not rendered.
- [ ] **Step 3: Render** the chapter with the final files. The existing test must still pass. Add a test: a results file carrying a `# PROVISIONAL` first line is refused even when the file name is canonical (already pinned) and a clean file renders.
- [ ] **Step 4: Update the lessons note:** tick each checklist step with the measured outcome; list what the real run changed in the spec assumptions.
- [ ] **Step 5: Commit** `docs(pathway_census): chapter A with real results`.

---

## Self-review

- **Checklist coverage:** lessons steps 1a (R1), 2-3 (R1, R2), 4 (R3), 5-7 (R3, R4, R5), 8 (R5), 9 (R6), 10 (R7), 11 (R8), 12 (R9), 13 (R10); step 1 (user gate) is recorded in the decisions block and R10.
- **Open choices left to the builder:** decoy lists and the `soluble_steps` list (stated with reasons), tree method for each step tree, chunk sizes (from the timed pilot), the `min_decoy_lwr` value (PROPOSED, record the reason).
- **Checks that need the user:** R2 step 6 (full submission, after timing), R10 step 2 (`verified`).
- **Known risks:** `epa-ng` via `pixi exec` is untested here (UNVERIFIED); deposited GFF naming is unchecked; per-step trees may be large for P450 decoys; the paired test covers only baseline genomes (spec §5.6).
