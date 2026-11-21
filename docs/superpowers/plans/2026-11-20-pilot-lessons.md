# Pathway census pilot lessons (carotenoid), 2026-11-20

Source: the controller ledger (`.superpowers/sdd/2026-11-20-pathway-census-carotenoid-pilot/progress.md`,
git-ignored, not in the repository), the final-review findings of 2026-11-20,
`analysis/pathway_census/README.md` and the module docstrings. Facts that committed files cite are
copied into this note (see "Reference curation record"). This note is the starting point for
the DHN-melanin plan (spec section 5.8).

Labels: MEASURED = read from data or a tool on 2026-11-20. PROPOSED = a builder's choice with no data
behind it. UNVERIFIED = not checked.

**No real results exist.** The census, the paired comparison, the absence check and the triage ran on
synthetic fixtures only. Nothing in this note is a biological finding.

## 1. What worked

- Rules-first order held for the triage bins (Task 10): the sort-order, duplicate and dtype rules
  were committed in docs-only commit 47529a4, before code commit 6a97602 (controller `git show` check).
- Real tools in tests, no mocked searches: hmmer 3.4, mafft 7.505, mmseqs2 17-b804f (Task 4) and
  ncbi-blast 2.17.0+ (Task 8) ran in the unit tests.
- The `verified` gate worked. Seven CLIs refuse a `verified: false` pathway without `--provisional`:
  `build_step_hmms.py`, `search_step.py`, `assign_steps.py`, `build_matrix.py`, `absence_check.py`,
  `paired_compare.py` and `triage.py`. Their provisional outputs carry a `.PROVISIONAL` name or header
  (Tasks 4, 7, 8, 9, 10; `assign_steps.py` got the gate in the final-review fix I1).
  `list_unverified_refs.py` only lists unverified references and has no gate.
- Reference curation: 23 refs, 22 verified against UniProt and PubMed (MEASURED, Task 3). A reviewer
  re-checked 6 of them (Task 3 review) (accession, PMID, sequence): all matched. Details below.

### Reference curation record (Task 3, copied from the Task 3 report)

Cited by `pathways/carotenoid.yaml` and Chapter A. MEASURED by the Task 3 implementer on 2026-11-20.
- Refs in the YAML: 23 (22 with accession, 1 PMID-only). Verified: 22. Unverified: 1 (PMID 16416328,
  Ojima 2006, astaxanthin synthase, no protein accession). Steps: 8. Steps with no verified ref: 0.
- Method for every ref with an accession: (1) fetch the UniProt entry (JSON and FASTA) on 2026-11-20 and
  check gene name, organism, length, and that the entry cites the PMID; (2) fetch the PubMed record and
  check that title and abstract name the gene and its function for the step; (3) `verified: true` only
  when 1 and 2 both pass. `verified: true` on a ref does not mean a person reviewed it.
- The Task 3 reviewer re-checked 6 refs (accession, PMID, sequence): all matched.
- The pathway flag stays `verified: false` until the user review (spec section 5.8).
- Nullable-boolean flags with explicit NA semantics (Tasks 6, 7) stopped unknown values from being
  counted as clean.
- Per-task review plus fix rounds caught Important issues in Tasks 2, 5, 6, 7, 8, 9, 10 before merge.
- The chapter scaffold (Task 11) renders with no results files, shows `NOT YET RUN: <file>` for
  each results-driven section, and is covered by `tests/test_chapter_render.py`. Quarto 1.7.31, knitr
  engine with R 4.5.2 (the Python 3.12 interpreter has no jupyter; quarto's own Python is 3.9 without jupyter).

## 2. What the data forced (by task)

- Task 1: `Ref.accession` is `str | None` because the brief's own bad-ref fixture needs it. `validate`
  requires an accession or a PMID.
- Task 2: `has_proteins` and `downloaded` are derived from BFD files on disk, not from the prior-annotation
  TSV. `genus` is never null: `genus_source` is `genus_field | parsed_species | unresolved`, and
  unresolved labels (for example "Ascomycota sp.") are not genera. MEASURED: 24,363 genomes; 23,379
  with proteins, 984 without; genus_source 24,014 / 49 / 300.
- Task 3: the loader needed type checks (list/bool fields, ref-less step under `verified: true`)
  before the curated definition could rely on it (ruling: implementer fixes minors 2, 3, 4).
  The user review gate opened here and stays open. Per ruling, Tasks 4-10 ran code and fixture tests only.
- Task 5: a step is `assigned` only with HMM evidence AND a placement with a clade and support >=
  `min_support`. Fusion only sets `fusion_candidate`. No placement tool is installed as a module
  (no epa-ng or pplacer), so the interface is a placements TSV.
- Task 6: BFD has no per-gene taxonomy table (MEASURED by the implementer). Contamination and
  haplotype flags need caller-supplied tables that do not exist. In a real run they are all NA.
- Task 7: copy buckets are exclusive with priority flagged > unknown > filtered, and
  `copy_filtered + copy_unknown + copy_flagged == copy_raw` per row. Genomes with `has_proteins`
  False get zero rows and are excluded from `n_genomes`.
- Task 8: `downloaded` False means "no BFD FASTA", not proven absence, so the outcome is `unconfirmed`.
  Same-genus relatives skip `genus_source == "unresolved"`.
- Task 9: the paired arm needs a deposited-protein matrix that is not built (see section 4). Reserved
  provider names `ALL` and `(missing)` raise ValueError.
- Task 10: a gene never neighbours itself. Output is sorted by `(protein_id, step_id)`, which changed
  the input order (documented). The signal/TM check is OFF by default (`soluble_steps` empty) and this
  is shown in params and on stderr.

## 3. Spec assumptions that failed or changed

1. **Spec section 3 genome counts were TSV-wide.** BFD-scoped counts differ: `baseline_ok` is 5,625 in
   BFD against 5,703 TSV-wide (MEASURED, Task 2). The TSV has 27,401 rows and 3,038 of them are not in
   `species.parquet` (MEASURED, README, Task 2). The spec text is not yet updated (deferred spec fix).
2. **`has_proteins` semantics.** It means "BFD `proteins.fa.gz` exists", not "the genome was annotated".
   `has_deposited_proteins` is the TSV flag (5,634 True, MEASURED). A zero copy count in a genome with
   `has_proteins` False is not evidence of absence.
3. **Flags are nullable booleans, not true/false (Task 6 ruling, ledger).** `NA` = not evaluated.
   `fillna(False)` is forbidden. Use `clean_mask` and `unknown_counts`. (No spec text was found that
   says flags are plain booleans, so no failed spec assumption is claimed here.)
4. **Fusion needs placement evidence (Task 5 ruling, ledger).** Spec section 5 line 114 says "A fused
   gene must count toward both steps", and section 5 step 3 (line 130) says "Fused genes get multiple
   step labels". Task 5 ruling: a step label is `assigned` only with HMM evidence AND a tree placement
   (clade and support >= `min_support`). Fusion only sets `fusion_candidate`. A fused protein without a
   placement for a step stays `unassigned` for it. The 50% domain-overlap threshold is PROPOSED.
5. **Fused self-neighbour.** A protein assigned to two steps could count as its own cluster neighbour in
   the triage `in_cluster` term. Fixed by the rule "different protein_id AND different step_id" (Task 10).
6. **Paired design limits.** The paired test covers only assemblies with both annotations. It says
   nothing about the genomes without a prior annotation (BFD-scoped: 18,055 `no_prior_annotation_present`
   plus 673 `no_prior_annotation_not_downloaded`, MEASURED, Task 2 report), which annotation is correct,
   or absence (spec section 5.6). The Fable review says "The 5,625 baseline genomes allow this test
   directly" (`2026-11-20-pathway-census-REVIEW_FABLE.md` line 49). 5,625 is the BFD genome-table
   `baseline_ok` count (MEASURED, Task 2), not a pair count: a pair also needs BFD proteins, deposited
   proteins and an existing deposited protein file (Task 9 pair rule). The real pair count is not run.
7. **No per-gene taxonomy in BFD.** Spec section 2 (line 61) says "The pipeline checks contamination before counting a copy". Task 6 found (MEASURED) that BFD has no per-gene taxonomy table, so the check needs a new best-hit table and a reference set labelled fungal / non-fungal. UNVERIFIED whether
   `swissprot_best.parquet` plus a taxid-to-lineage map is usable (memory notes call it "old format").
8. **Tree placement tool.** `upo_place` does not run unchanged on carotenoid steps. Only its `assign`
   subcommand is reusable (`--clades`, `--min-lwr`). Each step needs its own reference tree and
   `clade_members.tsv`.
9. **Reference purity.** B5LEC8 (astaxanthin synthase) cites PMID 18769767, a genome-organisation paper,
   not a catalysis paper (Task 3 minor). One ref (PMID 16416328) has no accession.

## 4. Still deferred: dependency-ordered checklist for the real-data pass

Each item names what it needs from earlier items.

Items marked [REQUIRED BEFORE REAL-DATA PASS] come from the final whole-branch review (2026-11-20). Step 1a was
inserted so that the numbers 1-13 and their "Needs" references stay unchanged.

1. [ ] User reviews `analysis/pathway_census/pathways/carotenoid.yaml` and `refs/carotenoid/UNVERIFIED.md`
       (spec section 5.8). Decide on `ylo1`, `cao2`, B5LEC8 and PMID 16416328. Set `verified: true` only
       after that. Needs: nothing. Gates all later items (no `--provisional` real run before this).
       The user also decides:
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] Decoy / outgroup sequences per step, with provenance (I6). With only pathway
             references in the reference tree, placement cannot reject anything, and the absence check
             calls decoy homologues `present_unannotated` (MEASURED on the test decoy: 72.8 bits at
             `min_bitscore` 50).
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] Any-of versus all-of domain rule (I10). `quality.py` `domain_ok` is all-of, but
             `carotenoid.yaml` lists alternatives `[TIGR03462, PF18916]` for lycopene cyclase.
             TIGR03462 is not Pfam, so BFD `pfam` cannot supply it.
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] Whether GGPP synthase counts toward completeness (M5). It is non-optional, with no
             carotenoid-specific clade defined.
1a. [ ] [REQUIRED BEFORE REAL-DATA PASS] Split the fused reference FASTA per domain (I5). `refs/carotenoid/phytoene_synthase.faa`
        and `lycopene_cyclase.faa` are byte-identical (`cmp`). Both steps would get the same HMM,
        `fusion_candidate` never fires, and the absence cover for a stand-alone phytoene synthase (about
        57% of the fused length, UNVERIFIED.md item 6) is `fragmented`. Needs: 1.
2. [ ] Build real step HMMs (`build_step_hmms.py`). Needs: 1 and 1a.
3. [ ] BFD search (`search_step.py`, `slurm/search_step.sbatch`) and dereplication. Set `min_score` per step
       (PROPOSED rule: lowest reference score; leave-one-out check open). Apply the dereplication rule
       (PROPOSED: >50,000 proteins at 90% identity) after counting the real family size. Needs: 2.
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] A batched driver that loops many genomes per job over all step HMMs, with the
             genome or chunk id in the output names. BFD has one protein file per genome, and
             per-genome runs of `search_step.py` overwrite each other (names carry the step only).
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] A searched-genome manifest that `build_matrix` requires. Today `has_proteins` means
             "the file exists", not "the genome was searched".
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] One timed job before sizing. 23,379 genomes x 8 steps would be 187,032 jobs at one
             job per genome-step.
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] Compress `candidates.tsv`; label the sbatch time and memory UNVERIFIED; do not
             swallow zstd errors.
       Done in the final-review wave: `search_step.py` writes a `step_id` column (I4 part).
4. [ ] Wire a placement tool and run placement for the candidates. `epa-ng` and `pplacer` are not installed
       as modules; iqtree, fasttree and raxml exist (MEASURED, Task 5). Build per-step reference trees and
       `clade_members.tsv`; output the placements TSV (`protein_id, step_id, clade, support`). Needs: 3.
5. [ ] `assign_steps.py`. Needs: 3 (candidates with `step_id`) and 4 (placements).
6. [ ] Build the contamination best-hit table (`protein_id, taxon_class`) and the haplotype
       `pair_identity` table (`protein_a, protein_b, pident, qcov`) for the ASSIGNED proteins and their
       scaffold neighbours. The fungal / non-fungal reference set does not exist yet. Needs: 5.
7. [ ] Run `quality_flags`; this produces the flags table that `build_matrix` consumes. Needs: 5 and 6
       (also `gene_info` joined with protein length, and `expected_domains`).
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] The `pfam.parquet` column mapping. Its columns are `pfam_id`, `pfam_acc`,
             `domain_i_evalue`, not `hmm_name` / `evalue`. Whether `pfam_acc` carries a version suffix
             is UNVERIFIED.
8. [ ] `build_matrix.py` and the phylum summary on BFD proteins; writes
       `results/carotenoid_phylum_summary.tsv`. Needs: 5 and 7.
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] (applies to steps 6-8) An explicit `require` decision for the matrix and the
             triage if the contamination table is not ready. `matrix.DEFAULT_REQUIRE` includes
             `contam_suspect`, so `copy_filtered` is 0 everywhere without it. `bins.DEFAULT_REQUIRE`
             includes `motif_ok`, which is always NA because the YAML has no motif, so every real copy
             is `unknown` in triage.
9. [ ] Deposited arm. Regenerate `results/genomes.tsv.zst` (the genome table now carries `asmid` and
       `folder`) and run the real pair count (`paired.count_pairs`). Decide the ID mapping: option (a)
       rewrite deposited IDs to `<locustag>_<id>`, or option (b) a `genome_of` parameter in
       `matrix.build_matrix` and `quality.quality_flags` (change surface: the `split("_")[0]` lines in
       `matrix.py` and `quality.py`, and their docstrings). Then run items 3-8 on the deposited proteins
       of the pair genomes. Needs: the designs of 1-8 (same HMMs, placement and assignment rules).
       - [x] [REQUIRED BEFORE REAL-DATA PASS] Deposited path from the TSV `folder` column:
             `<root>/<folder>/<folder>_protein.faa.gz` (done in the final-review wave, I7).
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] Collapse deposited isoforms to one protein per gene from the GFF (I8). BFD proteins
             are all `-T1` across 244,484,661 proteins. GCF_000001985.1 deposited has 10,638 proteins
             for 10,023 protein-coding genes, and 505 genes have more than one mRNA (one genome only,
             MEASURED by the reviewer). State that a copy is a gene.
       - [ ] [REQUIRED BEFORE REAL-DATA PASS] Deposited genomes have no `gene_info` equivalent. Deposited flags are then NA and
             `delta_filtered` is NA for every pair unless coordinates are parsed from the deposited
             GFF. Only raw `delta` is usable without that.
10. [ ] `paired_compare.py`. Needs: 8 (BFD matrix) and 9 (deposited matrix).
11. [ ] Absence check on SLURM (job size UNVERIFIED; target 1-1.5 h per job, from the user's HPCC guidance),
        merge chunk files into `results/carotenoid_absence_summary.tsv` (PROPOSED name) with
        `bin/merge_absence.py` (added in the final-review wave, I2). Needs: 8.
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] Per-domain references (item 1a) and a decoy-calibrated `min_bitscore` set per
              step, not only for step 8.
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] Use `bin/merge_absence.py`, not a shell recipe.
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] Consider renaming the outcome `present_unannotated` to a neutral `tblastn_hit`
              until a confirmation step exists.
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] Efficiency and correctness (M4): one decompression plus one tblastn per
              genome-step; the per-(query, subject) cover union can add distant partial hits;
              `-num_threads` is ignored with `-subject` (BLAST+ documentation, not measured); the
              sbatch has no array-task-to-locustags mapping.
12. [ ] Triage glue and `triage.py`. The glue joins signalp, tmhmm, wolfpsort, `gene_info` coordinates,
        `ref_identity`, `length` and `copy_filtered`. Check the real tmhmm / signalp / wolfpsort column
        names first (UNVERIFIED). Needs: 5, 7 and 8.
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] The `require` decision (see item 8).
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] Rename `lineage_expansion` (M3) before a testing partner sees the output. It is
              `copy_filtered >= 2` in one genome; it is not lineage-aware.
13. [ ] Re-render Chapter A with real numbers:
        `module load quarto/1.7.31; quarto render analysis/pathway_census/chapters/A_carotenoid.qmd`.
        Needs: 8, 10, 11, 12. Update spec section 3 with BFD-scoped counts.
        - [x] [REQUIRED BEFORE REAL-DATA PASS] The `#`-in-field check is DONE for SPECIES, STRAIN and `annotation_provider`
              (MEASURED by the reviewer: 0 of 24,363, 0, and 0 of 27,401 contain `#`).
        - [ ] [REQUIRED BEFORE REAL-DATA PASS] Unify the TSV readers (M1). Five readers behave differently: `build_matrix`,
              `paired_compare`, `triage` and `assign_steps` use `comment="#"`; `absence_check` strips
              only leading `#` lines; `build_matrix` does not read `.gz` despite its docstring. Use one
              shared reader before the real-data pass.

### Also recorded (final review, not checklist steps)

- M5: there is no completeness output. `Step.optional`, `copy_expect` and `role` are unused. The YAML
  has 4 non-optional steps, although the ledger says "7 required + ylo1".
- M10: the branch history holds 16 cache and result files (630,737 bytes) from commit 4b291a1. A squash
  merge removes them.
- M11: `paired_compare` does not check that the two matrices share the pathway file and `require`.

## 5. Deferred Minor findings by task (copied from the ledger)

- **Task 1 (defs.py):** silent defaults for `copy_expect` / `role`; scalar split into characters for
  `domains` / `fused_with`; `bool("false")` is True for `verified`; `verified: true` allows a ref-less
  step; missing species is `""`; no tests for duplicate id, invalid enum, positive-verified; `open()`
  has no utf-8; fixture placeholders uncommented. (Task 3 fixed minors 2, 3, 4 of these.)
- **Task 2:** minor-bare assert skipped under `python -O`; 9 `no_prior_annotation_present` rows have
  `protein_file_ok=1` and were not investigated; the `_RANK` suffix heuristic may null genuine genera;
  the tree plan and Task 8 groupby must keep `genus_source` visible.
- **Task 3:** B5LEC8 PMID 18769767 is not a catalysis paper (flag in UNVERIFIED.md or downgrade);
  the provenance test skips PMID-only refs; ylo1 refs Q1K615 / F6IBC7 not spot-checked; the reviewer
  did not re-run the full-suite count.
- **Task 4:** decoy and zero-hit tests rely on HMMER heuristic filters at `min_score=0` (no test that a
  homologue scores above a decoy); `.zst` protein input unsupported in `search_step`; `candidates.tsv`
  uncompressed; sbatch `--time=02:00:00` / `--mem=16G` unlabelled UNVERIFIED and above the 1.5 h target,
  zstd errors swallowed; `dereplicate` lacks `--threads`; HMM name truncated at the first dot; a missing
  refs file gives a bare `CalledProcessError`; `test_verified_pathway_gets_no_tag` mutates its fixture;
  empty `"${EXTRA[@]}"` needs bash >= 4.4.
- **Task 5:** fix commit 4b291a1 committed `.pytest_cache`, `__pycache__` and `results/genomes.tsv.zst`
  (Task 6 untracked them); the 50% fusion overlap threshold is PROPOSED and unvalidated; empty output
  gives an object-dtype `fusion_candidate`; `assign_steps.py` `comment="#"` may corrupt fields that
  contain `#`.
- **Task 6:** `clean_mask` raises a bare `KeyError` for a missing column (prefer `ValueError`);
  `unknown_counts` silently skips absent columns; the `by_step` lookup in quality.py is unchecked; no
  separate RED for the fix round.
- **Task 7:** test imports mid-file and unused `sys`; `len(m) == 2*8` hard-codes the carotenoid step
  count; CLI `has_proteins` lacks `dtype=str` (lowercase `true` auto-parsed) and has no CLI test for a bad
  value; task-7 report self-review text is stale; `phylum_summary` did not summarise `copy_flagged`
  (fixed in the final-review wave: `copy_*_total` columns).
- **Task 8:** `build_matrix.py` and `absence_check.py` each have their own `read_tsv` (consider a shared io
  helper); `gzip.open` handle not closed; `zstd -dc` error omits stderr; no zstandard-module branch;
  same-name `--locustags` files in different dirs collide (document); `select_absences` and the paralog
  union work per (query, subject), not per locus; sbatch has no `PC_EVALUATE` passthrough; the README
  merge awk was untested (replaced by `bin/merge_absence.py` in the final-review wave); each (genome, step) re-reads the genome.
- **Task 9:** README lacked the `comment="#"` hint for reading outputs (the chapter now does it); the
  `.PROVISIONAL` note is not in the `paired_table` file; `n_unpaired.values()` over-counts
  (`unknown_flag` is a subset); NA `has_proteins` is reachable only via the library; resampling is per
  step, not joint; carrying ASMID in `load_genome_table` is cleaner than the caller merge (done in the final-review wave).
- **Task 10:** the README may still say "input order" (grep found none); a soluble-step NA `tm` over-reports
  unevaluated by at most one; `_in_cluster` is O(n^2) per scaffold; `signal_peptide` is not dtype-validated
  in the library; the caller must build `signal_peptide`, `tm_helices`, `ref_identity`, `length`,
  `copy_filtered` and coordinates, and the real tmhmm / signalp / wolfpsort column names are unchecked.

## 6. Observations from Task 11 (no module code was changed)

- No bug found in module code during this task.
- The absence step writes per-chunk files only. The chapter expects a merged `carotenoid_absence_summary.tsv`.
  At Task 11 that merge had only an awk recipe, which could write a provisional merge under the canonical
  name. The final-review wave replaced it with `bin/merge_absence.py`. The name is a PROPOSED choice.
- Module outputs use `.PROVISIONAL` in the file name. The chapter loads only canonical names and refuses
  to show a `.PROVISIONAL` file as a result (it prints a notice instead). Since the final-review wave it
  also refuses a canonical name whose leading `#` lines contain `PROVISIONAL`.
- The chapter uses knitr + R because jupyter is not installed for `/usr/bin/python3.12`. If the lab
  later installs jupyter, the same logic ports to a Python cell. This is a convenience choice, not a design need.

## 7. What the DHN-melanin plan should reuse or change

Reuse:
- Package layout, the YAML definition loader with `verified` gate, provenance TSV with per-ref checks, and the
  PROVISIONAL tagging convention.
- Genome table (Task 2) as is; it now carries `asmid` and `folder`.
- Nullable-boolean flags, `clean_mask`, `unknown_counts`, and the copy-bucket invariant.
- Placement-required assignment rule and the `placements` TSV interface.
- Rules-first docs commits for any binning or threshold step. Absence check design, paired design with its stated limits.
- The chapter scaffold with `NOT YET RUN` placeholders and the render test.

Change:
- Decide the deposited-arm design (ID rewrite or `genome_of`) at the start. Do not hard-code the LOCUSTAG
  prefix format (it is 8 hex characters in all but 3 BFD genomes).
- Decide the contamination reference set and the haplotype pair-identity method before the matrix step.
- Pick a placement tool and test it before writing the assignment step.
- Hold the user review of the pathway YAML before Tasks 4-10, not in parallel with them. Fix the SLURM
  sizing with a real timing run (target 1-1.5 h per job, from the user's HPCC guidance in CLAUDE.md).
- Define a soluble flag in the YAML if a signal/TM check is wanted. Merge shared io helpers early.
- Name the final result files in the plan so the chapter and the CLIs agree.
