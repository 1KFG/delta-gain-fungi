# Design: pathway census for carotenoid and melanin biosynthesis across BFD

Written 2026-11-20. Status: APPROVED by the user 2026-11-20.
Template project: `analysis/enzymes_UPO/` (see `DESIGN.md`, `PLAN.md`).

Labels used in this document:
- **MEASURED**: checked against files in this repository.
- **PROPOSED**: my call. The user can reverse it.
- **UNVERIFIED**: from my memory of the literature. Must be checked against the cited paper before use.

---

## 1. Intent

### 1.1 Questions

Two new chapters in the DeltaGain-Fungi analyses. Each asks the same questions about one pathway:

1. **Completeness.** Which genomes carry all steps of the pathway? Which carry only some?
2. **Copy number.** How many copies of each step gene does each genome carry? How does that vary across phyla?
3. **Diversification.** Which steps show lineage-specific expansion, contraction, gain, loss or divergence of the active site?
4. **BFD effect.** Does BFD's uniform re-annotation of ~18k genomes without prior annotation reveal pathway differences that the NCBI-annotated genomes do not show?
5. **Prioritisation.** Which gene copies are the best candidates for later functional testing, and why?

Chapters:
- **Chapter A: carotenoid biosynthesis.**
- **Chapter B: melanin biosynthesis**, with two sub-pathways: DHN-melanin and L-DOPA melanin.

### 1.2 Why one shared tool

The deliverable for both chapters is the same shape: a table of pathway steps, and a copy-number matrix by genome and phylum. The census, tree-placement and prioritisation logic is identical. A shared tool, `analysis/pathway_census/`, takes a pathway definition file as input. Chapters A and B are reports on its output. A third pathway (for example pyomelanin or a siderophore) costs one new definition file.

## 2. Goals, non-goals, anti-goals

### Goals

- **G1.** A curated pathway definition for each of three pathways (carotenoid, DHN-melanin, L-DOPA melanin). Each lists steps, reference proteins with a citation, and expected domain architecture.
- **G2.** A per-genome, per-step copy-number matrix over all BFD genomes with proteins. Each row carries the flag `prior_annotation` (baseline_ok or none).
- **G3.** A phylum-level summary table: main synthesis genes × copy number across the kingdom. This is the headline table of each chapter.
- **G4.** Lineage-aware tests of expansion, contraction, gain and loss, with a species tree. The tests count related genomes as one lineage, not as independent observations.
- **G5.** A ranked candidate list of gene copies for functional testing, with the reason for each pick.
- **G6.** A direct test of the BFD effect (question 4), reported whether the answer is yes or no. Design decided by the user 2026-11-20: a **paired** test (§5.6). For baseline genomes, run the census on the deposited proteins and on the BFD proteins of the same assembly. Genomes with no prior annotation get descriptive counts only, with no inference about annotation effect.
- **G7.** Everything reproducible from the repository. MEASURED and UNVERIFIED numbers are labelled differently.

### Non-goals

- Tier-3 work: no prediction of substrate affinity, activity or kinetics. User decision 2026-11-20. The UPO calibration (r = +0.18 for the aliphatic pocket; clade explained most of the leave-one-out signal) is the reason. MEASURED in `analysis/enzymes_UPO/RESULT_T2_negative.md` and `STRUCTURE_VS_RESIDUES.md`.
- Structure prediction and docking as a ranking input. Structure may appear only as a descriptive figure.
- Pathway regulation (transcription factors, promoters). Out of scope for this version.
- A complete, publication-grade fungal phylogeny of all ~23.6k genomes. The tree built here is a working tree on a subsample (§5.7), not a deliverable in itself.
- Wet-lab work.

### Anti-goals (the work fails even if it runs)

- **Pfam counts presented as pathway counts.** Pfam families such as amine oxidase, PKS, multicopper oxidase and P450 are much broader than any pathway step. A step call needs HMM plus tree-based orthology.
- **"Missing step" read as defect.** DHN-melanin is expected in many Ascomycota and absent elsewhere. L-DOPA melanin is characterised in Basidiomycota such as Cryptococcus. Absence in a lineage is not an incomplete pathway. Completeness must be judged against the lineage's expected pathway. Whether this expectation holds across phyla is UNVERIFIED.
- **Annotation gap read as biology.** A step absent from a genome could be a missing prediction, a fragmented assembly or a true loss. The pipeline must check assembly quality (BUSCO completeness, N50 from `busco_genome` and `asm_stats` tables) and, for key absences, search the genome nucleotide sequence with tblastn before calling a loss.
- **Phylogeny read as biology.** Closely related genomes are not independent. Sampling in BFD is uneven (for example many Lecanoromycetes), so raw counts per phylum will mislead.
- **A confident ranking built on little data.** The prioritisation is triage into bins (§5.5), not a ranked score. No term is validated against activity data, and terms can conflict.
- **Completeness inflated by ubiquitous families.** Phytoene synthase (SQS_PSY) and polyprenyl synthase domains occur in nearly all fungal genomes (Fable review, VERIFIED by that reviewer). Completeness counts only pathway-specific clades, not the bare Pfam family.
- **Foreign genes read as fungal.** Bacterial or algal carotenoid genes can enter fungal assemblies, especially lichens. The pipeline checks contamination before counting a copy.
- **Duplicates from the assembly, not the genome.** Haplotype duplication inflates copy number. The pipeline checks for it.
- **A catalogue nobody can act on.** G5 must end in a list a testing partner can pick from.

## 3. Constraints

- Compute: UCR HPCC. Follow the existing project rules in `analysis/enzymes_UPO/DESIGN.md` §3: constraint list, `$SCRATCH` use, `/usr/bin/python3.12`, no `BASH_SOURCE` in SLURM scripts, zstd/gz compression, jobs sized toward 1-1.5 h.
- Data: BFD parquet tables and seqdb. See memory/`bfd-database-layout`. MEASURED gaps: no EggNOG, KEGG-KO, antiSMASH or Phobius table. Step calls therefore cannot rely on KO annotation.
- Available per-protein tables (MEASURED from the layout note): `pfam`, `cazy`, `merops`, `signalp` (positives only), `tmhmm`, `predgpi`, `wolfpsort`, `idp`, `gene_info` (coordinates, for neighbourhood analysis).
- Prior-annotation split: `prior_annotation_status.tsv`. MEASURED 2026-11-20 by recount of the `status` column: 5,703 `baseline_ok`, 18,145 `no_prior_annotation_present`, 2,359 `no_prior_annotation_not_downloaded`, 1,184 `annotated_NOT_downloaded`, 10 `annotated_folder_missing_protein`. (The earlier note's 5,511 is stale.) `samples.parquet` lists 24,363 genomes per the Fable review (not re-checked by me).
- Genomes not downloaded cannot be searched by tblastn. They are reported as "unconfirmed", not as losses.
- Whole-BFD scale: step families can contain very many proteins. Anything costlier than an alignment runs on dereplicated representatives.

## 4. Approaches considered

| Option | Description | Trade-off |
|---|---|---|
| A | Two separate projects copied from `enzymes_UPO` | Fast start. About 70% duplicated code. Fixes made twice. |
| **B (chosen by user)** | Shared `pathway_census` tool, pathway definition files, chapters as reports | More design up front. Cheap to add pathways. |
| C | Pfam-only census | Fast. Wrong for this question (see anti-goals). Use only as a first look. |

## 5. Design

### 5.1 Layout (PROPOSED)

```
analysis/pathway_census/
  README.md
  pathways/
    carotenoid.yaml
    melanin_dhn.yaml
    melanin_ldopa.yaml
  refs/            # reference FASTA per step, one provenance TSV
  hmms/            # per-step HMMs (plain text; hmmsearch cannot read compressed)
  bin/             # scripts (python3.12)
  sql/             # census queries against BFD parquet
  slurm/
  results/         # compressed
  chapters/
    A_carotenoid.qmd
    B_melanin.qmd
```

The chapters follow the existing `reports/*.qmd` style and publish through the existing Quarto site (`reports/_quarto.yml`). The exact hook-up is left to the builder.

### 5.2 Pathway definition file

Each step entry holds: step id, enzyme name, EC number where one exists, reference proteins (accession, species, PMID), expected domains, whether the step is expected to be single- or multi-copy, whether it is catalytic or auxiliary, and the lineages where it is expected. Every reference needs a PMID or accession checked by a person or by a cited database lookup before the definition is marked `verified: true`.

### 5.3 Initial step lists (UNVERIFIED, from memory)

**Carotenoid** (to be checked against the literature, including Neurospora, Mucorales, Phycomyces, Rhodotorula and Xanthophyllomyces work):
- GGPP synthase (crtE; al-3 in Neurospora).
- Phytoene synthase. Often fused with lycopene cyclase as a bifunctional crtYB/al-2 or carRA in several lineages. A fused gene must count toward both steps.
- Phytoene desaturase (crtI; al-1; carB in Mucorales).
- Lycopene cyclase (crtY or the fused domain).
- Carotenoid cleavage oxygenases, as separate steps: carX (retinal-forming) and carT (torulene cleavage) are different enzymes per the Fable review (PubMed-checked by the reviewer). Add cao-2 (Neurospora) and downstream steps such as the ylo-1 aldehyde dehydrogenase in neurosporaxanthin synthesis; the builder checks each one.
- Astaxanthin synthase (P450, crtS/asy) in Xanthophyllomyces.

**DHN-melanin** (UNVERIFIED, corrected after review): polyketide synthase (PKS1/ALB1/WA type), scytalone dehydratase (Scd1/Arp1; Rsy1 is the same enzyme per the Fable review, so it is not a separate step), THN reductase (Arp2/Buf1, short-chain dehydrogenase), multicopper oxidase (Abr1/Abr2). Ayg1 hydrolase belongs only to the heptaketide route, because some PKS1 enzymes (for example in Colletotrichum) make T4HN directly. Model two routes as variants, with Ayg1 optional. Pathway-specific transcription factors are auxiliary and not scored in v1. *Aspergillus terreus* lacks DHN-melanin while relatives have it, so "expected pathway per lineage" is not a fixed taxonomic rule. Completeness is therefore reported as observed step sets, not as a deficit against an assumed expectation.

**L-DOPA melanin** (UNVERIFIED): this is a one- or two-enzyme family census (tyrosinase PF00264 and laccase), not a pathway with a completeness score. It is reported as a short section inside Chapter B. Whether a tyrosine hydroxylation step is needed in fungi is to be checked.

The builder confirms each list, and may add or remove steps, before step 1 of the plan completes.

### 5.4 Pipeline

1. **Curate references.** Build `refs/` and the provenance table per §5.2. Record what could not be verified.
2. **Candidate retrieval.** Per step, search BFD proteins with the step HMM, using the reference set to build it. For large families, start from the Pfam hits in the `pfam` table to limit the search, then confirm with the step HMM. For small or rare families, search with the step HMM directly and skip the Pfam prefilter. Example: PF18916 (lycopene cyclase) has about 1,032 proteins BFD-wide per the Fable review, and a Pfam prefilter may miss the cyclase. Family sizes reported by the reviewer (VERIFIED by that reviewer, not re-checked by me): adh_short 2.6 M proteins, P450 2.0 M, ketoacyl-synt 374 k, Cu-oxidase about 220 k, tyrosinase 166 k. Use dereplicated representatives for these.
3. **Step assignment.** Place candidates on a per-step reference tree (as in `upo_place`) and assign each to a step clade. A candidate that lies outside every reference clade is `unassigned`, not dropped. Fused genes get multiple step labels.
4. **Domain and quality checks.** Keep expected-domain presence, catalytic-motif presence where the references define one, length relative to references, and flags for truncation. Count at domain level for fused genes. Add a contamination check (for example taxonomic placement of the best hits and of neighbouring genes on the scaffold) and a haplotype-duplication check.
5. **Presence/copy matrix.** Per genome and step: copy count, copy count after quality filtering, and genome quality flags.
6. **Absence check.** For genomes that lack a step that their close relatives have, run a tblastn of the step references against the genome. Report whether the gene is present but unannotated, fragmented, or truly absent.
7. **Lineage-aware tests.** Use the subsampled species tree (§5.7) to test gain, loss and copy-number change per step. Use a method that accounts for shared ancestry, for example phylogenetic comparative models or ancestral reconstruction with sensible sampling control. The method is left to the builder; it must state what it controls for. Per-genus or per-order downsampling is a required sensitivity check.
8. **BFD effect.** Paired test per §5.6. Run the census on both protein sets for each baseline genome and compare step presence and copy number. Report effect size and confidence interval. For no-prior-annotation genomes, report descriptive counts only.
9. **Prioritisation.** Assign copies to triage bins per §5.5. Output a table with a reason column.
10. **Chapters.** A and B reports, with the pathway-step table, the phylum copy-number table, tree figures, the absence-check summary and the candidate list.

### 5.5 Prioritisation (replaces tier 3): triage bins

Reviewer finding (Fable, opinion): this is defensible as triage into bins, not as a ranking. The terms can conflict and cannot be validated without activity data. The output is therefore bins (for example "tier 1: intact, in cluster, expressible, divergent from references") with the reasons shown, not a single score.

The bins use only terms measurable now:

- Domain and catalytic-motif integrity (hard filter, not a score term): truncated or motif-broken copies are excluded or flagged.
- Phylogenetic position: lineage-specific expansion, a paralog with an altered active-site residue relative to references, or the nearest characterised reference being far.
- Genomic context: other pathway genes within a set window on the same scaffold (coordinates from `gene_info`). The window is a parameter to fix with the builder.
- Expressibility: length, transmembrane helices (`tmhmm`), signal peptide (`signalp`), subcellular prediction (`wolfpsort`).
- Genome quality of the source genome.

Rules: each term must be defined before looking at results, and each is shown as its own column in the output. No fitted weights. The output is a list of bins, not an activity prediction. Bin definitions are left to the builder and fixed before results are viewed.

### 5.7 Species tree (user decision 2026-11-20)

MEASURED 2026-11-20: no species tree covering BFD exists in the repository or in the BFD run folder (searched `*.nwk`, `*.treefile` to depth 4; only transcription-factor gene trees were found). The user decided the tree will be built in parallel with the pathway work, on a subsample, because a tree of all genomes is too large.

- The census matrix (§5.4 step 5) is still computed for all genomes.
- The lineage-aware tests (step 7) and the tree-based BFD-effect comparison (step 8) run on the subsampled tree. Per-genome rows are collapsed to the tree's tips (copy number per tip, with the within-tip range kept).
- Tip unit: the user suggested major classes or genera. PROPOSED: one representative per genus, chosen by best assembly quality (BUSCO completeness, then N50), with the representative's `prior_annotation` status recorded and both statuses kept where a genus has both. A class-level backbone is the fallback if genus-level is too large. The builder confirms the size after counting genera in `species.parquet`.
- Inclusion rule (user decision 2026-11-20):
  - Prefer genomes with BUSCO completeness >= 80% as genus representatives.
  - A genus with no genome at >= 80% still gets a representative if its best genome is > 50%. Mark these tips `low_busco` and carry the flag into every test.
  - *Antarctolithica* is always included, whatever its BUSCO score. MEASURED 2026-11-20: `species.parquet` lists 8 genomes as "Ascomycota sp. Antarctolithica obscura" with an empty GENUS field. Grouping by GENUS would drop them. The builder must handle taxa with an empty GENUS (for example by parsing the SPECIES string) and check other unplaced taxa for the same problem. I have not checked the BUSCO values of these 8 genomes.
  - Genomes at <= 50% are excluded except for this named exception. The builder lists every genus lost to the 50% cut so the user can add exceptions.
- Tree method (user decision 2026-11-20): build the tree with `phyling` / `nf_phyling`, run on the one-per-genus list once that list is fixed. The tip list is therefore a hard input to the tree job. Invocation (user, 2026-11-20): `phyling` takes a folder of protein files, run as `nextflow run stajichlab/nf_phyling`. The user also requires the full marker information to be understood and documented before the tree is used: which marker set (HMM library and version), how many markers, how many markers each tip carries, and how missing markers are handled. UNVERIFIED: I have not checked the installed `phyling` version, its marker-set options, or its scale limits. The builder checks these and records them in `analysis/pathway_census/` before the tree job runs. Which marker set to use, and how the genus list is turned into phyling input, is left to the builder.
- Constraint: genera that are sparsely or unevenly sampled are the weak point. Tests must state the number of tips per phylum and flag phyla with too few tips for inference.
- Dependency: step 7 and step 8 cannot run until this tree exists. Steps 1–6 and 9 do not depend on it, so the two workstreams run in parallel.

### 5.6 What would count as a BFD finding

Why the design is paired (user decision 2026-11-20, on Fable's finding). The baseline and BFD-only genome sets differ before any biology: 13% of baseline genomes versus 21% of BFD-only genomes are below 80% BUSCO, only 835 of 7,568 species occur in both sets, and the baseline has 669 annotation providers (JGI is the largest). These figures are the reviewer's, VERIFIED by that reviewer with queries; I have not re-run them. A lineage-matched comparison would mix an annotation effect with a sampling effect.

Paired design: for each genome in `baseline_ok`, take the deposited proteins (`NCBI_ASM/<acc>/protein.faa` or provider equivalent) and the BFD proteins for the same assembly. Run the same census on both. The difference is an annotation effect on identical assemblies. Reviewer's count was 5,625 baseline genomes with both sets; the builder verifies how many have both. Report separately by annotation provider where numbers allow.

A BFD finding needs all of:
- a step-level difference in the paired comparison that is consistent across providers and survives genome-quality controls,
- and, for absences, a tblastn confirmation (or an "unconfirmed" label where the genome is not downloaded).

What the paired test does not show: whether BFD annotation of the 18,145 genomes without prior annotation reveals pathway genes that deposited annotation would have found or missed. For those genomes the report gives descriptive counts only. Anything else is an observation, not a finding. If no difference survives, that is reported as a negative result.

### 5.8 Order of work (user decision 2026-11-20)

1. **Chapter A, carotenoid, first, as a pilot.** It tests the method on one pathway.
2. **Chapter B, DHN-melanin, second.**
3. **L-DOPA** as a short family census inside Chapter B.

The shared tool is built during the carotenoid pilot. Lessons from the pilot go back into the design before DHN-melanin starts. The species-tree workstream (§5.7) runs in parallel.

Step lists are drafted by the builder and reviewed by the user before any census runs. Chapters start as draft `.qmd` files. Site publishing is a later step.

## 6. Risks

- Reference curation is the slowest and most error-prone step; wrong references corrupt every later table.
- Large families (PKS, multicopper oxidase, P450) make orthology assignment ambiguous. Expect many `unassigned`.
- BFD sampling is uneven across lineages.
- Gene-model quality differs between BFD and NCBI-annotated genomes in ways unrelated to biology. The funannotate and NCBI pipelines predict different gene models, which could produce an apparent BFD effect.
- Fused genes (carotenoid) break one-gene-one-step counting.

## 7. Left to the builder

- Final step lists and reference selection (§5.3).
- Tree-based orthology tooling, reusing `upo_place` where it fits.
- Comparative-method choice for §5.4 step 7.
- Neighbourhood window and score form (§5.5).
- Site integration for the chapters.
- Order of work and job sizing.

## 8. Reviews and open items

- Builder check 2026-11-20: done. Questions were settled with the user or by the edits above.
- Independent review by Fable: `2026-11-20-pathway-census-REVIEW_FABLE.md`. Its findings 1, 2, 4, 5, 6, 7 and 8 are applied here. Finding 3 (no species tree) is handled by §5.7. The review read an earlier version of this document.
- Open: the exact marker set and scale limits of `phyling` (§5.7); the count of baseline genomes that have both deposited and BFD proteins (§5.6); the BUSCO values of the 8 Antarctolithica genomes.
- Approved by the user 2026-11-20.
