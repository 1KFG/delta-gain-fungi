# Candidate effectors and secretion biology across Mucoromycota — analysis plan

**A design document.** Written 2026-11-16 to extend `analysis/mucor_secretome/REPORT.md`, which is
the completed first pass on the two MorenoLab *Mucor* strains (JEM1, JEM2) isolated from
*Marchantia polymorpha*.

**Revised 2026-11-18.** Two changes, both recorded where they bite rather than only here:
§4.5's prerequisite is **done** — the four missing predictors now exist for all four panel
proteomes (PR #15) — and **CotH has been re-scoped out of the effector work** into a separable
**Q4 / Phase F** on the family's evolution within Mucoromycota. No analysis phase has been run.

Written to the standard in `WRITING.md` (repo branch `upo-page`, commit `70f6dde` — it is **not**
on `main`; see §12.3). Standard of evidence taken from `analysis/enzymes_UPO/REPORT.md`.

::: {.callout-note}
## What this is, and where it stops being reliable

This is a plan, not a result. Every forward-looking number in it is a **cost or a sample size**,
never an outcome. Where a number describes something that exists now — how many Mucorales
genomes are in BFD, what a column contains — it was obtained by query on 2026-11-16 and the query
is named.

Three things it cannot do, stated before anything else:

1. **It cannot establish that any protein is an effector.** Nothing here is a bench measurement.
   The strongest claim available at the end of Phase D is "candidates for assay, ranked", which is
   what the first pass already says about its five AlphaFold 3 hits.
2. **The central comparative question — "is the MorenoLab pair typical of *Mucor*?" — may turn out
   to be unanswerable with BFD as it stands.** Counting secreted proteins across genomes compares
   gene-model predictors as much as it compares organisms. §2.5 and §7 Phase A specify the
   measurement that decides this, and §2.5 declares the kill switch in advance: if between-pipeline
   variance is of the same order as between-lineage variance, the composition question is reported
   as unanswerable rather than answered.
3. **The CotH *host-interaction* question is dead as posed, and the family has been re-scoped.**
   A Pfam census run while writing this plan finds PF08757 in 273 of 275 Mucorales genomes **and**
   in 25 of 26 Neocallimastigomycota genomes and 100 of 100 Mortierellales (§5.2a0), so
   "Mucorales-restricted" was never a premise this plan could build on. As of 2026-11-18 CotH is
   **not a focus of the effector/secretion work at all**; what remains is the evolutionary
   question — is PF08757 one family, and what is its history across Mucoromycota — carried as a
   separable **Q4 / Phase F** that depends on no other phase. The lifestyle test (§5.2e) and the
   plant-partner search (§5.2g) are deferred, not deleted. See the scope note at §5.2.

A fourth, specific to this dataset: **JEM1 has no conspecific genome anywhere in BFD**, because the
only public "*Mucor hiemalis*" assembly is mislabelled (first pass §1). Nothing in this plan can
separate strain-specific from species-specific for JEM1. That is a property of the world, not of
the method.
:::

---

## 0. What the first pass settles, and what it leaves open

### 0.1 Settled — do not redo these

| Settled | Evidence | Where |
|---|---|---|
| JEM1 is *M. hiemalis*, JEM2 is *M. circinelloides* | ITS vs UNITE v10 at 99.3 % / 98.6 %, with both BFD assemblies run as controls | REPORT.md §1 |
| BFD's `GCA_977110905.1_gzMucHiem1` is **not** *M. hiemalis* | its own ITS does not return *M. hiemalis* in the top 5 | REPORT.md §1 |
| Neither proteome carries detectable *Fusarium* | two independent tests; rDNA (positive control fires at 5 hits, cleaned assemblies 0) and protein-level diamond against 6 of 2,099 BFD *Fusarium* proteomes with JEM1 as the never-contaminated negative control (JEM2 shows no excess at any identity bin) | REPORT.md §2 |
| JEM1/JEM2 are annotated with **BFD's own invocations** for SignalP 6, TMHMM 2.0, PredGPI, TargetP 2.0, Pfam and dbCAN | invocations tabulated against `Fungi_BFD/nextflow/modules/BFD/*` | REPORT.md §3 |
| 1,344 secretion candidates across the two strains (863 secreted, 171 GPI, 266 membrane, 41 weak signal, 3 mitochondrial) | `bin/build_secretome.py` | REPORT.md §4 |
| Single-sequence structure prediction works on enzymes and fails on orphans | 35 of 37 chitin-active proteins assigned a fold; 52 of 346 `unknown_small` assigned | REPORT.md §6 |
| MSAs rescue about a sixth of the unfoldable set | 76 of 426 reach AF3 pTM ≥ 0.7; median pTM 0.30 | REPORT.md §9.5 |
| A structure-based host screen is feasible **and can be calibrated on its own controls** | positive 6Q76 pair 0.91, negative 0.13, panel background of 128 fungal × host pairs median 0.16 / p95 0.38; call line set at iPTM ≥ 0.5 | REPORT.md §9.1 |
| Domain-family compatibility inflates iPTM without any interaction | AtCERK1 LysM × *Marchantia* LysM scores **0.84**, above every fungal candidate in the panel | REPORT.md §9.3 |
| Eight silent pipeline failures, with their general shapes | `ERRORS_FOUND.md` | — |

The last two rows are the most load-bearing for this plan. §9.3 is the reason §6.4 below makes
domain-compatibility decoys a **designed** control rather than an accident, and `ERRORS_FOUND.md`
is the reason §7 Phase 0 exists at all.

### 0.2 Open, or never touched

- **Nothing comparative beyond a clustering background.** The first pass clustered JEM1/JEM2's
  candidates with 195,284 secreted candidates from 292 BFD Mucoromycota genomes to measure
  *reach*. It never compared genomes to each other on secretome **composition**, never normalised
  for annotation differences, and never asked whether JEM1/JEM2 are typical of *Mucor*.
- **No lineage-restricted surface families.** CotH is not mentioned anywhere in
  `analysis/mucor_secretome/`. Neither is any other Mucorales-restricted family.
- ~~**Four BFD tables unused**~~ — **resolved 2026-11-17 (PR #15).** `deeploc`, `wolfpsort`,
  `merops` and `idp_summary` now exist for JEM1, JEM2 and both *Marchantia* proteomes, run with
  BFD's own invocations and merge scripts. The first two are half the sorting evidence (§4);
  `idp_summary` bears directly on first-pass limitation 6 (low-confidence folds are ambiguous
  between prediction failure and real disorder); `merops` covers secreted peptidases, a candidate
  class the tiers do not reach. **See §4.5 for coverage rules before using them** — in two of the
  four, absence is a zero rather than a NULL.
- **No lifestyle axis** of any kind. BFD has no lifestyle column (verified: `species.parquet` has
  `PHYLUM/SUBPHYLUM/CLASS/ORDER/FAMILY/GENUS/SPECIES/STRAIN/NCBI_TAXONID/BUSCO_LINEAGE/SUBCLASS`
  and nothing else).
- **No expression data, no binding assay, no mutant.** Unchanged, and §10 names what each would
  settle.
- **The AF3 panels were built from the superseded `fair-esm` fold calls**, not the `transformers`
  ones that `REPORT.md` §6 now reports (`HANDOFF.md` §3). Provenance has to be stated or rebuilt
  before anything new is layered on them.

### 0.3 A ninth silent failure, found while writing this plan

`pfam.parquet` stores **versioned** accessions. Verified 2026-11-16: `WHERE pfam_acc IN
('PF08757','PF16731')` returns **0 rows**, while `WHERE pfam_id ILIKE '%CotH%'` returns 11,881
domain rows all carrying `pfam_acc = 'PF08757.18'`.

That matters here because the two sides of the first pass's comparison strip the version
differently:

| Source | Code | Accession as written |
|---|---|---|
| The 294 BFD Mucoromycota genomes | `sql/export.sql` — `select ... pfam_acc ... from read_parquet('pfam.parquet')`, no stripping | `PF01476.27`, `PF01522.27` |
| JEM1, JEM2, both *Marchantia* proteomes | `bin/annotate_to_tables.py:245` — `acc.split('.')[0]` | `PF01476`, `PF01522` |

`bin/build_secretome.py:134` assigns tier T1 on `set(pfam_accs) & set(CHITIN_PFAM)` where
`CHITIN_PFAM = {"PF01522", "PF01476"}` — an exact string intersection. One side is versioned and
the other is not, so **the Pfam branch of the tier rule can fire for the focal strains and never
for the background**, which is precisely the property `README.md` claims the design guarantees ("so that a candidate here means what it means across the 292 BFD
Mucoromycota genomes").

**The accession asymmetry is verified; its consequence is not yet.** Verified by query on
2026-11-16: `SELECT DISTINCT pfam_acc FROM read_parquet('pfam.parquet') WHERE pfam_acc LIKE
'PF01476%' OR pfam_acc LIKE 'PF01522%'` returns exactly two values, **`PF01476.27` and
`PF01522.27`** — both versioned, neither equal to the bare string the tier rule tests against. So
the intersection at `build_secretome.py:134` cannot match on the BFD side and can match on the JEM
side.

What remains to check, and is Phase 0's first task: **whether the T1 Pfam branch ever fired for a
BFD genome.** The test is two lines — count BFD-side T1 proteins whose only T1 evidence is a Pfam
accession rather than a CAZy family. If that is zero across 292 genomes while being non-zero for
JEM1/JEM2, the first pass's tier counts for the background need recomputing, and the stated
equivalence between focal and background tiers was not holding. It does **not** touch the
JEM1/JEM2 counts themselves, and it does not touch any structure or AF3 result, which are keyed on
protein IDs. The likely magnitude is small — CBM50 covers most of what PF01476 would have caught —
but "likely small" is not a measurement.

It has the shape `ERRORS_FOUND.md` names twice: **a check that can only ever return one answer**.

---

## 1. Three questions, kept apart

They are nested and they fail differently, so they get different evidence standards.

| # | Question | What would answer it | Dominant failure mode |
|---|---|---|---|
| **Q1 Composition** | How do Mucoromycota genomes differ in SSP / cysteine-rich / candidate-symbiosis secreted content, and is the MorenoLab pair typical of *Mucor*? | per-genome counts and fractions across a defined comparison set, with an annotation-error bar | **gene-model predictor differences and deposition sampling** (§2.4, §2.5) |
| **Q2 Lineage-specific surface factors** | Which protein families are restricted to Mucorales, and do they track lifestyle? | unbiased family discovery + copy-number counting against a false-discovery floor + a phylogenetically controlled lifestyle test | **assembly fragmentation inflating copy number; genus/lifestyle confounding** (§5.1) |
| **Q3 Interaction potential** | How does one get from "secreted, small, cysteine-rich" to a shortlist that would survive review? | a pre-declared ranking rule, structure, and a host co-fold screen with decoy controls | **domain-compatibility scores read as interactions** (first pass §9.3) |
| **Q4 CotH evolution in Mucoromycota** *(separable, see §5.2)* | Is PF08757 one family, and what is its duplication/loss history across Mucoromycota? | a gene tree with a kinase outgroup, reconciled against an in-house species tree | **a kinase-like profile over-recruiting, so that a "copy number" is not a family count** |

Q1 can be answered badly in a way that looks fine. Q2 can be answered in a way that is correct for
mammals and irrelevant for plants. Q3 can be answered in a way that is literally true and reads as
a finding. Each section below names the control that blocks its own failure mode.

**Q4 is deliberately kept out of Q2** (scope decision, 2026-11-18). The preliminary census in
§5.2(a0) killed the premise that made CotH a host-interaction story — it is not
Mucorales-restricted, and the anaerobic gut fungi carry an order of magnitude more per genome.
What it left behind is a genuine comparative-genomics question about the family's history in
Mucoromycota, which needs **none** of Q2's apparatus: no lifestyle curation, no *Marchantia* host
panel, no secretome pipeline. Treating it as a fourth question rather than a branch of the second
is what makes it cheap and independently publishable; see §5.2's scope note.

---

## 2. The comparison set, and why

### 2.1 What is actually in BFD

All counts below: `duckdb` 1.1.3 over
`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/species.parquet`, 2026-11-16.
Denominator for every percentage is the 23,680 rows of that table unless stated.

| BFD `PHYLUM` | genomes | species | note |
|---|---|---|---|
| Mucoromycota | **294** | 118 | Mucorales 275 / Umbelopsidales 14 / Endogonales 5 |
| Mortierellomycota | **100** | 55 | Mortierellales 100, 14 genera |
| Glomeromycota | **60** | 31 | Glomerales 36, Diversisporales 15, Archaeosporales 3, Paraglomerales 3, Entrophosporales 3 |
| Kickxellomycota | 190 | 169 | Kickxellales 171 / Harpellales 11 / Dimargaritales 7 / Ramicandelaberales 1 |
| Chytridiomycota | 89 | 62 | |
| Neocallimastigomycota | 26 | 9 | small, but it carries CotH at ~202 proteins per genome — §5.2a0 |
| Basidiobolomycota | 18 | 17 | |
| Blastocladiomycota | 14 | 9 | |
| Entomophthoromycota | 12 | 9 | |
| Zoopagomycota | 11 | 11 | |
| Ascomycota | 19,412 | 5,176 | background only |
| Basidiomycota | 3,275 | 1,666 | background only |

### 2.2 The circumscription trap: "Mucoromycota" means two different things here

**BFD's `PHYLUM` column uses a narrow circumscription.** Under Spatafora et al. 2016, Mucoromycota
comprises Mucoromycotina + Mortierellomycotina + Glomeromycotina. BFD files Mortierellomycota and
Glomeromycota as **separate phyla**, so `PHYLUM = 'Mucoromycota'` returns 294 genomes where the
modern phylum has **454** (294 + 100 + 60).

Anyone comparing this plan's numbers with the first pass will read two discrepancies. The first is
benign and explained above: the first pass used the narrow sense. The second is that the first
pass reports **292** Mucoromycota genomes in `REPORT.md` §5 and **294** in `README.md`; the reason
for the difference is **not recorded anywhere in the directory**, and it is a Phase 0 item to
determine it (most likely two genomes contribute no surviving candidate, but that is a guess and
must be checked, not asserted). **This plan uses "Mucoromycota *sensu lato*" (454 genomes, 204
species) and says so every time.** The narrow set is referred to as "Mucoromycotina *s.s.*" (294).

This is not pedantry: Glomeromycotina is the plant-symbiotic lineage of the phylum, and
*Marchantia paleacea* forms arbuscular mycorrhizal symbiosis while *M. polymorpha* has lost it
(first pass §7). Excluding Glomeromycota because of a column definition would remove the single
most relevant comparison from a plant-symbiosis study.

**And `SUBPHYLUM` cannot be used to fix it**: it is NULL for all 23,680 rows of both
`samples.parquet` and `species.parquet` (`analysis/delta_gain_funnel/REPORT.md`, "Upstream issue
found"; a fix is in flight). Group by `PHYLUM` + `ORDER` and reconstruct subphylum from the order
list, committed as an explicit mapping file.

### 2.3 Sampling bias, with numbers

| | |
|---|---|
| *Rhizopus* genomes in BFD | **110 of 275 Mucorales genomes (40.0 %)** — across only **8 species** |
| *Mucor* genomes | 50 of 275 (18.2 %), 22 species |
| The two genera together | 160 of 275 (58.2 %) |
| Largest BioProject, `PRJNA475137` | **57 of 454 Mucoromycota *s.l.* genomes (12.6 %)** |
| Four largest BioProjects | 119 of 454 (26.2 %) |
| GenBank-only (`GCA_`) assemblies | **436 of 454 (96.0 %)**; 18 are `GCF_` |
| `asm_stats.annotation_status` | `annotated` for **all 451** of 454 that join to `asm_stats` — the column carries no information here, and 3 genomes do not join at all |

Read together these say something specific. The **RefSeq-vs-GenBank axis that drove the UPO
novelty artefact is nearly degenerate in this clade** — at 96.0 % GCA there is almost nothing to
stratify on — so the deposition confounder here will not express itself as nr-novelty. It expresses
itself as **genus and project concentration**: four in ten Mucorales genomes are *Rhizopus*, eight
species deep, largely from clinical mucormycosis sequencing. Any statement of the form "Mucorales
have on average N copies of X" is, without correction, a statement about *Rhizopus* clinical
isolates.

### 2.4 The trap this project has already hit twice

Both instances are documented in this repository and both took the same shape: **a biological
reading of what is really a record of where genomes were deposited and how they were annotated.**

| Instance | The apparent signal | What it was | Where |
|---|---|---|---|
| Delta-gain funnel | Pfam families enriched in new genomes carry +14.6 to +14.9 points more novel clusters than an order-matched null | 1,033 Y1000+ genomes (5.8 % of 17,827 new-side genomes) deposited to GenBank unannotated. Counting them as known drops enriched-family novelty from 71.0 % to 17.5 % — **75.4 % of the novelty**, against 12.3 % for size-matched controls. The excess does not shrink, it **inverts** to 23.1 points *below* controls | `analysis/delta_gain_funnel/REPORT.md` |
| UPO census | 18 % of catalytic-core UPO representatives look sequence-novel | 2.0 % novelty for RefSeq-annotated genomes against 18.7 % for GenBank-only, z = 7.4. All 7,055 representatives have an `nr` hit | `analysis/enzymes_UPO/REPORT.md` §2 |

The funnel case is the better template because of **what made it an argument: the control**. Y1000+
adds ~883 distinct species at roughly one genome each, so it inflates cluster counts in *any*
family. The size-matched control families lost 12.3 % of their novelty to the same join. That
contrast — 75.4 % against 12.3 %, median family 86.6 % against 0.8 % — is the whole result. Without
it the join proves nothing.

**The operational rule this plan inherits: every count-based comparative claim is computed three
times.** Once as stated; once with one genome per species; once with the largest contributing
BioProject (and, separately, the largest contributing genus) removed. All three go in the table.
If the claim does not survive all three, it is reported as a deposition effect, in its own section,
with the magnitude — exactly as the funnel report did.

And a matched null, as the funnel had: the identical discovery pipeline run on a lineage where
"nothing special" is the expected answer (§5.1d).

### 2.5 The gene-model predictor confounder — the live one for Q1

This is the confounder that the two cases above do *not* cover, and it is the one that most
threatens Q1. Secretome size is a function of the annotation pipeline before it is a function of
the organism: predictors differ in how many short ORFs they call, and a short, single-exon,
signal-peptide-bearing ORF is exactly an SSP. A genome annotated by a pipeline that calls more
short genes has a larger SSP count for no biological reason.

Four measures, in increasing order of what they cost:

**(a) Report fractions, not counts, as the primary statistic.** SSPs as a fraction of the
proteome; cysteine-rich SSPs as a fraction of the secretome. A predictor that over-calls short
genes inflates numerator and denominator together. This does not remove the bias — short-gene
over-calling raises the SSP fraction specifically — but it removes the first-order effect, which
is proteome size.

**(b) Covariates, declared before fitting.** Regress the statistic on total protein count,
assembly N50 and contig count (`asm_stats`), and BUSCO completeness and duplication
(`busco_genome.parquet`), and report the Mucor pair's position as a **studentised residual with a
prediction interval**, not as a rank. "JEM1 is 1.4 SD above the Mucorales mean after adjusting for
proteome size and BUSCO completeness, 95 % PI [−0.9, +3.7]" is a statement; "JEM1 has more SSPs
than most Mucorales" is not.

**(c) The annotation-pipeline control — the decisive one.** Take **k = 8** Mucorales assemblies
spanning the genus range, re-annotate each with the same funannotate version and parameters that
produced JEM1 and JEM2, run the identical BFD predictor stack over the result, and compute the
same SSP fraction. The distribution of |ΔSSP fraction| between BFD's annotation and ours, **on the
same assembly**, is the measurement error of every count in Q1. Report it as an error bar on every
figure. This is the step whose absence is why the two earlier cases were only caught afterwards.

**(d) Kill switch, declared now.** Let *s* be the median |ΔSSP fraction| from (c) and *IQR* the
interquartile range of SSP fraction across Mucorales species in the primary set. **If *s* > 0.5 ×
*IQR*, Q1 is reported as unanswerable with BFD as it stands**, the measurement is published as the
result, and Phases B and C proceed on Q2/Q3 only. Declaring this before the measurement is what
makes a negative outcome a result rather than a disappointment — the pattern `enzymes_UPO` used
for the Boltz-2 run, whose 14 models were discarded when the pre-declared switch fired.

### 2.6 The comparison set itself

| Tier | Set | n | Why |
|---|---|---|---|
| **Focal** | JEM1, JEM2 | 2 | the subject |
| **Congeners** | *Mucor* in BFD | 50 genomes / 22 species | the only set that can say "typical of *Mucor*". JEM2 has 24 conspecifics here; **JEM1 has none** |
| **Order** | Mucorales | 275 / 106 / 36 genera | the family-restriction denominator |
| **Phylum *s.l.*** | + Mortierellomycota + Glomeromycota | 454 / 204 | adds the plant-symbiotic (Glomeromycotina) and soil-saprotroph (Mortierellomycotina) arms |
| **Near outgroup** | Kickxellomycota, Zoopagomycota, Entomophthoromycota, Basidiobolomycota | 231 / 206 | sister phyla — "restricted to Mucoromycota" needs them |
| **Far background** | all remaining BFD genomes | ~22,995 | the restriction test (§5.1b) is a presence/absence search, not a clustering, so the whole corpus is affordable |

**Primary unit is the species, not the genome.** One genome per species, chosen by a rule written
before the choice is made — highest BUSCO complete, ties broken by N50, ties broken by ASMID
ascending — giving 204 Mucoromycota *s.l.* species. Genome-level numbers appear only as a
secondary table with the *Rhizopus* share stated beside them.

### 2.7 Lifestyle labels: curated, cited, and confounded

There is no lifestyle column. Labels must be made, and the making is the weakest link in Q2, so it
is constrained:

- A committed TSV, `data/lifestyle.tsv`, one row per **species** in the comparison set, with
  columns `species, lifestyle, evidence, citation, curator, date`. One citation per row. No row
  without one.
- Classes: `plant_associated` (endophyte, rhizosphere, bryophyte-associated, plant pathogen),
  `animal_clinical`, `saprotroph` (soil, dung, litter), `mycoparasite`, `food_industrial`,
  `unknown`. **`unknown` is a class, not a deletion** — it is reported with its own n in every
  table, because dropping it silently would be a second sampling filter on top of deposition.
- Labels are assigned **blind to CotH copy number**, and the file is committed before the counts
  are joined to it. A `git log` on that file is the audit trail.
- **Lifestyle is confounded with genus and with BioProject.** *Rhizopus* is 40.0 % of Mucorales
  genomes and heavily clinical. A two-sample test of plant-associated against clinical would
  largely be a test of *Rhizopus* against everything else. §5.2e specifies the phylogenetic
  treatment this requires.

---

## 3. Operational definitions

Stated with thresholds, and with the **sensitivity to those thresholds treated as a quantity to
measure**. Defaults continue the first pass so the two are comparable; the grid is what makes the
defaults defensible.

### 3.1 Secreted

Unchanged from `bin/build_secretome.py` so the background stays valid:

| call | rule |
|---|---|
| signal peptide | SignalP 6 probability ≥ **0.7** (pool exported at ≥ 0.5) |
| `membrane` | a TMHMM helix starting > **5 aa** after the cleavage site — a helix inside the signal peptide *is* the signal peptide |
| `gpi_anchored` | a PredGPI `GPI-anchor` feature; kept, labelled |
| `mitochondrial` | TargetP mTP at probability ≥ **0.7**; excluded |
| `secreted` | signal peptide, none of the above |

Two mandatory derivations, both because of known BFD defects:

- **GPI status is derived, never read.** `is_gpi_anchored` in the DuckDB build is TRUE for every
  protein: PredGPI writes a `Chain` row per protein and a `GPI-anchor` row only for anchored ones
  (223,938,474 vs 3,299,703 rows). Use `read_parquet('predgpi.parquet') WHERE feature =
  'GPI-anchor'`, as `sql/export.sql` already does. A fix is in flight; re-check on rebuild.
- **5′ completeness from the sequence, never from the flag.** `gene_transcripts.parquet` predates
  Fungi_BFD `fb1d7cd`; ~10 % of models are wrongly flagged incomplete BFD-wide and about half in
  Mucorales and Glomeromycota (`analysis/enzymes_UPO/REPORT.md` §6). A model that genuinely lost
  its 5′ end lost its signal peptide and must not be in a secretome, so the check matters — it
  just has to come from `gene_CDS` + the assembly, as `bin/check_codons.py` does.

### 3.2 Small secreted protein (SSP)

**Mature length** — after removing residues 1..`peptide_end` from `signalp.signal_peptide` — not
full length. The first pass's T3 used full length ≤ 300 aa; switching to mature length moves the
boundary by the signal peptide (typically 18–25 aa) and is the more defensible definition. **The
change must be reported, with the count it moves**, because it makes the new T3 non-identical to
the old one.

Default **≤ 300 aa mature**. Grid: {150, 200, 250, 300, 400}.

### 3.3 Cysteine-rich

Cysteines counted in the **mature** sequence. Default ≥ **4** cysteines **and** ≥ **2 %** — both,
as in `build_secretome.py`, because either alone misbehaves at one end of the length range.
Grid: n_cys ∈ {2, 4, 6}, fraction ∈ {0.01, 0.02, 0.03, 0.04}.

### 3.4 Candidate

**"Candidate" is never a word on its own in any output.** It is always one of a small number of
named rules, each written here and each carrying its own count:

| Rule | Definition |
|---|---|
| `C-SSP` | secreted ∧ mature ≤ 300 aa |
| `C-CRP` | secreted ∧ ≥ 4 Cys ∧ ≥ 2 % Cys |
| `C-EFF` | `C-SSP` ∧ `C-CRP` ∧ no Pfam domain — the classical effector prior, the first pass's T3 |
| `C-CHI` | chitin/chitosan-active: CAZy CE4, GH18, GH20, GH46, GH75, CBM12/14/18/50, or Pfam LysM / polysaccharide deacetylase — the first pass's T1 |
| `C-PLA` | dbCAN substrate naming a plant polysaccharide, or cutinase / lipase / peroxygenase Pfam — the first pass's T2 |
| `C-RES` | secreted ∧ family restricted to Mucorales by §5.1's index |
| `C-PEP` | secreted peptidase: a `merops` hit at a declared identity and coverage, **new in this plan** |

A protein can satisfy several. Reporting is by rule, not by a precedence chain, with the overlap
matrix shown — the first-pass tiers are mutually exclusive by precedence, which hides how much
`C-CHI` and `C-EFF` overlap.

### 3.5 Threshold sensitivity is a measurement, not an assumption

For every headline comparative statistic, recompute it over the **full grid** of §3.2–3.3 — 5
lengths × 3 cysteine counts × 4 fractions × 3 SignalP cuts (0.5, 0.7, 0.9) = 180 cells — and emit a
**threshold-response table**. Cost is one pass over a table that already exists; this is cheap and
there is no excuse for not doing it.

The reporting rule: **a conclusion is stated as a conclusion only if it holds across the whole
grid.** Otherwise the plan states the sub-range over which it holds, as a range. "JEM1 is above the
*Mucor* median SSP fraction at every one of 180 threshold combinations" and "JEM1 is above the
median only when mature length ≤ 200 aa (36 of 180 cells)" are different claims, and the second is
still worth publishing — it localises the effect to short proteins.

---

## 4. Combining the four sorting predictors

### 4.1 They are not four independent votes, and the plan must not pretend they are

Available per protein: `signalp` (signal peptide + cleavage site), `targetP` (SP / mTP / cTP),
`deeploc` (DeepLoc 2.1; columns verified 2026-11-16 — `localizations`, `signals`,
`has_signal_peptide`, `has_mito_transit_peptide`, `membrane_types`, and 14 `p_*` probability
columns including `p_extracellular`), `wolfpsort` (`localization`, `score`). Plus `tmhmm` and `predgpi`, which are
topology modifiers, not sorting predictors.

SignalP 6, TargetP 2 and DeepLoc 2 are all deep models keying on the N-terminal signal and trained
on overlapping curated sets. **WoLF PSORT is the one genuinely different method** — amino-acid
composition and sorting-signal k-NN, not a neural sequence model. So an agreement count of 4/4 is
not four independent lines of evidence, and no output should describe it as one. How
non-independent they are is itself measurable and cheap: report the pairwise agreement rate
between each pair of predictors across the comparison set, with its denominator, instead of
asserting a correlation nobody has measured.

### 4.2 The rule

**The primary secretion call does not change.** It stays SignalP-based as in §3.1, because the
292-genome background and the 1,344 first-pass candidates were built that way and changing it
would discard that work.

The other three enter as a **0–3 agreement score**, reported alongside, never as a filter in the
primary analysis:

| predictor | counts as agreeing when |
|---|---|
| TargetP 2.0 | `prediction = 'SP'` and `probability ≥ 0.5` |
| DeepLoc 2.1 | `p_extracellular ≥ 0.5`, or `localizations` contains `Extracellular` |
| WoLF PSORT | top `localization` is `extr` |

Exact column semantics — whether `localizations` is a single label or a comma-joined list, what
`score` means for WoLF PSORT — are **Phase 0 schema probes**, not assumptions. The thresholds above
are provisional until those probes run.

### 4.3 Disagreement is a quantity to report, not noise to resolve

Emit the full 2×2×2×2 contingency table over the comparison set, denominator = every protein with
a SignalP call at p ≥ 0.5. Then the diagnostic that matters:

**Compare the disagreement rate in JEM1/JEM2 against the BFD Mucorales distribution.** If the
focal strains disagree more often than the background, that is an annotation or input-format
artefact in the locally-run predictors, not biology, and it invalidates the consensus before any
biology is read off it. This is the same logic as the first pass's backend-comparison step, which
is the only thing that caught the pLDDT units bug (`ERRORS_FOUND.md` §3).

### 4.4 Where a decision is forced

Only one place forces one: the bounded shortlist that goes to structure prediction (§6). There the
rule is **SignalP ≥ 0.7 AND agreement ≥ 2**, and the plan reports how many candidates that costs.
Proteins dropped by it go to a named `disagreement_set.tsv` with their four calls, published as a
section, not discarded — negative results get sections (`WRITING.md` rule 4).

### 4.5 What must be run before any of this — **DONE 2026-11-17**

This section recorded that DeepLoc, WoLF PSORT, MEROPS and the IDP predictor had never been run
on JEM1, JEM2 or either *Marchantia* proteome, while existing for all ~23,500 BFD genomes, so the
focal strains could not be placed on the consensus axis at all.

**All four have now been run** (PR #15; `slurm/deeploc_panels.sh` GPU, `slurm/annotate_panels_part2.sh`
CPU), each invoked as its `Fungi_BFD/nextflow/modules/BFD/*` module invokes it — DeepLoc 2.1 model
`Fast` on `cuda`, matching BFD's production profile rather than `test.config` — and merged by BFD's
**own** `merge_*.py` / `gather_AIUPred.py`. Column layouts were checked against
`idp_summary.parquet` and `wolfpsort.parquet` and match field-for-field.

Outputs in `out/panels/`. **Read the coverage rules before using them:**

| table | rows | coverage |
|---|---|---|
| `deeploc.parquet` | 61,950 | one row per protein, all four proteomes |
| `wolfpsort.csv` | 61,936 | all but 14 MPTAK proteins of <30 aa, which WoLF PSORT drops |
| `merops.csv` | 12,020 | 341 / 327 / 489 / 522 proteins with a protease hit |
| `idp_summary.csv` | 30,960 | one row per protein *having* an IDP region ≥30 aa |

**Absence is a zero, not a NULL** in `merops` and `idp_summary` — `gather_AIUPred.py` does
`if idp_residue_count == 0: continue`. §4.2's consensus rule must not read a missing
`idp_summary` row as a missing measurement. The 14 WoLF PSORT gaps *are* NULLs, all in MPTAK;
no protein in JEM1, JEM2 or MPALE is under 50 aa.

Sanity-checked against the lineage: the fraction of proteins carrying an IDP region ≥30 aa is
median 0.509 across BFD's 273 Mucorales genomes (p10 0.452, p90 0.603); JEM1 is 0.536 and JEM2
0.576. A broken AIUPred run would not land on the median.

Note the DeepLoc SIF is lab-built and local; **no sequence left the cluster** (§11).

---

## 5. Lineage-specific surface and host-association factors

### 5.1 Discovery first — CotH is the worked example, not the hypothesis

Starting from CotH and stopping there would answer a question nobody asked. The family-discovery
pass comes first and CotH is then one row in its output, which is also the only honest way to say
whether CotH is unusual.

**(a) Cluster.** MMseqs2 cascaded clustering over the **secreted** proteins of the 454
Mucoromycota *s.l.* genomes plus the 231 near-outgroup genomes, at 30 %, 50 % and 70 % identity
with 80 % coverage of the shorter sequence. The 50 %/80 % tier reproduces the first pass's
clustering so the two are comparable. Before clustering anything: check whether the corpus-wide
cluster sets built for the delta-gain funnel
(`/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/results/family_clusters_*`) can be
reused — they are c30/c50/c70/c90 over the whole corpus and would make step (b) nearly free. If the
layout has changed, fall back to the profile search in (b).

**(b) Test restriction against the whole corpus, not the comparison set.** For every cluster,
build an HMM from its alignment and `hmmsearch` it against all ~23,500 BFD proteomes. A
presence/absence search over the corpus is one pass; clustering 249 M proteins is not. Restriction
index per family:

```
restricted_to_Mucorales  =  (members in Mucorales) / (members anywhere in BFD)
breadth                  =  distinct Mucorales GENERA containing a member   (denominator 36)
depth                    =  distinct Mucorales SPECIES containing a member  (denominator 106)
secreted_fraction        =  members with a signal peptide / members
```

**(c) Call a family "Mucorales-restricted" only with breadth.** Require `restricted ≥ 0.95` **and**
`breadth ≥ 5` genera **and** `depth ≥ 10` species. A family in one genus is a genus-specific
expansion; a family in one BioProject's genomes is a deposition artefact. Rank by
`restricted × breadth × secreted_fraction` and annotate every cluster with its Pfam content — a
"new restricted family" that is a divergent member of a known Pfam clan is not new, and saying so
is cheaper than being corrected.

**(d) The false-discovery floor — the control that makes this an argument.** Run the identical
pipeline on **Mortierellales** (100 genomes, 55 species, 14 genera) and on a **random matched
Ascomycota set** (275 genomes drawn to match Mucorales' species-per-genome and genus counts), and
report how many "restricted families" the pipeline invents there. That number is the floor below
which Mucorales results mean nothing. This is the funnel's control structure (§2.4) transplanted.

**(e) The three-way recomputation** of §2.4 applies to every family count: as stated, one genome
per species, largest-BioProject-removed.

**Expected output:** a ranked table of candidate Mucorales-restricted families with their counts,
denominators and Pfam annotation, against a stated false-discovery floor.

**One caveat on CotH as a sanity check of this pass.** It is tempting to say "CotH should appear
in the output, and if it does not the pipeline is broken". That is only a valid check if CotH
proteins carry a classical signal peptide, and this pass runs over secreted proteins only —
surface-displayed does not imply SignalP-positive. So: report the fraction of PF08757 proteins
with a SignalP call at p ≥ 0.5 **before** using CotH's presence or absence to judge the pipeline.
If that fraction is low, CotH legitimately falls outside a secretome-based discovery pass, and the
restricted-family search has to be run a second time over **all** proteins, not just secreted ones,
with the cost that implies.

### 5.2 CotH in full

> **Scope, set 2026-11-18 — read before running any of this.**
>
> **CotH is not a focus of the effector/secretion work.** It entered this plan as the worked
> example of a Mucorales-restricted host-interaction factor, and the census in (a0) below
> established that it is not Mucorales-restricted: Neocallimastigales carries ~202 copies per
> genome against Mucorales' ~19, and Mortierellales carries it in 100 of 100 genomes. The
> premise that justified its place in Q2 is gone.
>
> What survives, and is **wanted as its own line of work**, is the evolutionary question (Q4):
> *is PF08757 one family, and what is its duplication and loss history across Mucoromycota?*
> The census makes that more interesting than the original framing did — a family sitting at
> ~19 copies per genome in one subphylum and ~202 in an unrelated early-diverging lineage is
> either one ancient family with two independent expansions or two things a kinase-like profile
> cannot tell apart. Either answer is worth having.
>
> **Run now (Q4, the evolution track):** (a) family definition, (b) defensible copy counting,
> (c) the gene tree, (d) the restriction test with denominators. These need an HMM, a species
> tree and IQ-TREE — no lifestyle curation, no host panel, no secretome.
>
> **Deferred (the host-interaction track):** (e) the lifestyle PGLS and (g) the *Marchantia*
> partner search. Both are downstream of a plant-association hypothesis that the census does
> not support, and (e) is gated on §2.7's lifestyle curation — days of a person's time. Do not
> spend that time on CotH. (f) is retained in a reduced form because its kill switch applies to
> any copy-number claim, including an evolutionary one.
>
> **Keep the honest framing below regardless of track.** Nothing here can establish that CotH
> does anything in a plant, and the evolutionary result must not be written as though it hints
> that way.

CotH ("spore coat protein H") is a kinase-like family; in Mucorales, CotH3 was identified as an
invasin binding host GRP78/BiP on endothelial cells, and CotH7 as binding integrin β1 on alveolar
epithelium. **All of that is mammalian-host work.** It is the starting literature, not a
transferable premise, and §5.2g is the part of this section that treats it as such.

Two premises travel with the family — that it is Mucorales-restricted, and that its copy number
varies meaningfully across Mucorales. (a0) below puts a first number on both, and **the first one
does not hold up.**

#### (a0) A preliminary census, run while writing this plan — and it does not say what the premise says

The Pfam CotH domain is `PF08757.18`. Counting **distinct proteins** carrying it at
`domain_i_evalue < 1e-5`, joined to `species.parquet` by LOCUSTAG (duckdb, 2026-11-16):

| Lineage | proteins | genomes with ≥ 1 | denominator (genomes in BFD) | species |
|---|---|---|---|---|
| Mucorales | **5,178** | 273 | 275 | 106 |
| **Neocallimastigomycota** (Neocallimastigales) | **5,062** | **25** | **26** | 9 |
| Mortierellales | 1,033 | 100 | 100 | 55 |
| Umbelopsidales | 77 | 14 | 14 | 8 |
| Basidiobolales | 45 | 17 | 18 | 16 |
| Helotiales (Ascomycota) | 46 | 42 | — | 27 |
| ~10 further Ascomycota orders | 3–14 each | 1–14 each | — | — |

**The "Mucorales-restricted" premise does not survive a first look.** Mucorales carries 5,178
PF08757 proteins across 273 genomes (~19 per genome); the anaerobic gut fungi carry 5,062 across
**25** genomes (~202 per genome), and Mortierellales carries it in 100 of 100. Scattered
Ascomycota hits reach 42 Helotiales genomes.

**Three reasons not to treat that table as a result, which is also the argument for §5.2a–d:**

1. It counts **proteins, not genes** — no transcript collapse, no domain-coverage requirement, no
   split-model detection. §5.2b exists because of exactly this.
2. It is **one profile with one threshold.** CotH is kinase-like, and kinase-like profiles
   over-recruit. The Neocallimastigomycota and Ascomycota hits may be homologs, may be a
   repeat-expanded gut-fungal family that happens to score, or may be protein kinases. **The
   tree in §5.2c is what decides, and nothing short of it does.**
3. It is **not corrected for anything** in §2 — not deposition, not annotation, not per-species
   aggregation.

What it does establish is a design constraint: **the plan must not be written as a test of whether
a known-restricted family tracks lifestyle.** The restriction itself is in question, §5.2d is
therefore a real test rather than a formality, and the headline finding of this workstream may end
up being about the family's distribution rather than about plant association.

Per-genome means within Mucorales, same query, same caveats — a provisional picture of the
variation the copy-number work has to explain, not a measurement:

| Genus | genomes | proteins | per genome |
|---|---|---|---|
| *Backusella* | 14 | 422 | 30.1 |
| *Rhizopus* | 109 | 2,395 | 22.0 |
| *Thamnidium* | 3 | 61 | 20.3 |
| *Mucor* | 50 | 924 | 18.5 |
| *Lichtheimia* | 17 | 305 | 17.9 |
| *Cunninghamella* | 10 | 121 | 12.1 |
| *Rhizomucor* | 5 | 35 | 7.0 |

A 4.3-fold spread between *Backusella* and *Rhizomucor* is the kind of number that invites a
story. It is also exactly the kind of number that §5.2b4 (copy number vs assembly fragmentation)
and §2.5c (the annotation error bar) exist to deflate first.

#### (a) Identify the family — do not trust one Pfam scan

Three independent detections, reported as three numbers and reconciled in the text if they differ:

1. **Pfam**, `pfam_acc LIKE 'PF08757%'` — note the `LIKE`; an exact match on `'PF08757'` returns
   zero rows (§0.3). Record the gathering-threshold behaviour: BFD runs `hmmsearch --cut_ga`, so
   divergent members below the GA cut are invisible to this route by construction.
2. **A curated profile.** Align the characterized *Rhizopus delemar* CotH proteins and their
   *Mucor* homologs, `hmmbuild`, and search the corpus at an e-value threshold declared in advance.
   This is the route that can find members the GA cut misses, and the one that must be used for
   the restriction test, because "absent outside Mucorales" measured with a GA cut is partly a
   statement about the cut.
3. **Profile-profile / structure.** For the members found by 1 and 2, a Foldseek search of their
   predicted structures against the corpus's folded set, to catch members too divergent in sequence.
   Optional; run only if 1 and 2 disagree by more than 10 %.

Reconcile the three in the text (`WRITING.md` rule 6). Two counts that differ and are not explained
look like an error.

#### (b) Count copies defensibly

Copy number is where this analysis is most likely to produce a confident artefact. Five rules:

1. **Unit is the gene, not the protein.** Collapse alternative transcripts through
   `gene_proteins.gene_id`.
2. **Require domain coverage.** A hit must cover ≥ 70 % of the HMM (`hmm_from`/`hmm_to` against
   `pfam_len`) to count as a copy. Partial hits are counted separately and reported separately.
3. **Detect split models.** Two models on the same contig within 10 kb, each covering complementary
   parts of the HMM and neither covering ≥ 70 %, are a probable single gene broken by an annotation
   or assembly error. **Count both ways and publish both columns** — `copies_strict` and
   `copies_merged`. Where they differ, the difference is the uncertainty.
4. **Correlate with assembly quality.** Regress copy number on contig count, N50 and BUSCO
   duplication. CotH genes are reported to sit in repeat-rich regions, where a fragmented assembly
   inflates counts and a collapsed one deflates them. **If copy number correlates with contig count
   across the comparison set, the lifestyle test in (e) is not interpretable and must say so.**
5. **Per species, not per genome.** Median across that species' genomes, with the range. Otherwise
   *Rhizopus*, at 110 genomes and 8 species, sets the Mucorales average by itself.

The annotation-pipeline control of §2.5c applies unchanged: re-annotate the same k = 8 assemblies
and see how much the CotH count moves. **That spread is the error bar on every copy number
reported.**

#### (c) A tree, because counts alone cannot tell expansion from loss — **this is now the centre of the work**

With (e) and (g) deferred, (c) stops being a supporting step for a lifestyle test and becomes the
deliverable. It answers Q4 directly, and it is the only thing that can.

Align all retained members (MAFFT L-INS-i or the HMM alignment), trim, and build a maximum
likelihood tree with IQ-TREE and model selection, with SH-aLRT and UFBoot support. Non-Mucorales
members exist — (a0) finds them in Neocallimastigomycota, Mortierellales, Basidiobolales and
several Ascomycota orders — so the tree includes them as ingroup, with fungal protein kinases as
the outgroup. Reconcile against a species tree built **in house** from BUSCO single-copy orthologs
of the comparison set, rather than taking a published tree.

**The first question the tree settles is whether there is one family at all.** CotH is
kinase-like and kinase-like profiles over-recruit; if the Neocallimastigomycota and Helotiales
hits fall outside the Mucorales clade on the kinase side of the outgroup, then "5,062 proteins in
25 genomes" is a profile artefact and not a copy number. **Resolve this before quoting any
cross-lineage comparison, including the one in (a0) and the one in `HANDOFF.md`** — both are
currently stated as raw Pfam hits, which is all they are.

**Then the history across Mucoromycota.** Reconcile the gene tree against the species tree and
report gains and losses per node, with the three-way recomputation of §2.4 applied to every count:

| Pattern | Reading |
|---|---|
| One clade spanning Mucorales + Mortierellales + Umbelopsidales, copy number rising toward Mucorales | an ancestral Mucoromycota family, expanded in the order — the headline result, and the one (a0) points at |
| Mucorales and Neocallimastigales members reciprocally monophyletic | two independent expansions of an old family in two unrelated early-diverging lineages; neither is "the CotH lineage" |
| Mucorales members nested *inside* a Neocallimastigales radiation, or vice versa | horizontal transfer is on the table and must be tested as such (topology test against a constrained tree), not asserted from a surprising placement |
| No consistent clade; members interdigitate with the kinase outgroup | the family is not a family at this profile; redefine from the tree per (a2) and recompute everything downstream |

A count says *Rhizopus* has more. A tree says whether that is a *Rhizopus*-specific duplication
clade or loss in everything else — and those are different events with different causes. The
species tree is reusable: §5.2e needs it if the lifestyle test is ever revived, and §5.1's
restricted-family pass needs one too, so it is not an extra cost.

**Denominator discipline carries over.** Per-genome copy number is meaningless without the
assembly check of (b4) — a fragmented assembly splits one gene into two models, and
Neocallimastigomycota assemblies are not uniformly good. Report copies against contig count
before reporting copies against lineage.

#### (d) The restriction test, with a denominator

Search all ~23,500 BFD proteomes with the curated profile from (a2), and repeat with the Pfam GA
cut, and report both. Output: proteins, genes, genomes and species carrying a hit, **by phylum,
each with its denominator**, per species rather than per genome.

Given (a0), the first row below is already unlikely and the test is a real one:

| Outcome | Reading |
|---|---|
| Hits confined to Mucorales (275 genomes) | the published restriction claim replicated at kingdom scale with a denominator. **(a0) makes this improbable** — it would require the Neocallimastigomycota, Mortierellales and Helotiales hits to fail the gene-level and coverage filters, which should be checked rather than hoped for |
| Hits across Mucoromycota *s.l.* and the near outgroups, in a single clade in (c)'s tree | the family is older than Mucorales; "Mucorales-restricted" is wrong as stated, and any lineage signal is in **copy number**, not presence. This is the outcome (a0) points at |
| Hits across Mucoromycota *s.l.*, Neocallimastigomycota **and** scattered Dikarya, **not** monophyletic in (c) | the profile is over-recruiting — probably onto kinase-like folds — and the Pfam count is not a family count. Then the family must be redefined from the tree and everything downstream recomputed on the redefined set |
| Neocallimastigomycota hits monophyletic and separate from the Mucorales clade | two independent expansions of an old family in two unrelated anaerobic/early-diverging lineages — an interesting result in its own right, and one that weakens any Mucorales-specific functional reading |

CotH is kinase-like, so the third row is a live risk and **the tree in (c) is the only thing that
separates "a distant homolog" from "a protein kinase"**. Concretely: include a set of
unambiguous fungal protein kinases as an outgroup in (c)'s alignment and check which side of the
tree the Neocallimastigomycota and Helotiales hits fall on.

#### (e) The lifestyle comparison, specified honestly — **DEFERRED (2026-11-18)**

> Not being run. It tests a plant-association hypothesis the (a0) census does not support, and it
> is gated on §2.7's lifestyle curation, which is days of a person's time that should not be spent
> on this family. **The specification is kept intact** because if the (c) tree shows a clean
> Mucorales-specific expansion clade, the question becomes askable again — and the power
> calculation below is the thing that should be run *first* if it ever is. Nothing in Q4 depends
> on this subsection.

**Not** a two-sample test of plant-associated against clinical. Lifestyle is confounded with genus
(§2.7) and the phylogeny is the confounder.

- **Power first, before the test.** Using the observed variance in `copies_merged` across Mucorales
  species, compute the minimum detectable difference at 80 % power for the actual n of each
  lifestyle class. **Report the MDE whether or not the test is run.** If the MDE exceeds the
  effect size the literature reports between *Rhizopus* and other Mucorales, the study cannot
  answer the question and the plan says so instead of reporting a p-value.
- **Primary test:** phylogenetic generalised least squares on the in-house species tree, copy
  number ~ lifestyle, with proteome size and contig count as covariates. PGLS, not a t-test,
  because related species share lifestyles and share copy numbers for reasons of descent.
- **Secondary:** sister-pair contrasts — congeneric species pairs differing in lifestyle. Report
  the number of usable pairs; if it is under 5, report it as anecdote and label it so.
- **Mandatory negative control:** run the identical PGLS on **three protein families chosen for
  having no plausible lifestyle association** (e.g. a housekeeping Pfam of comparable copy-number
  variance). If lifestyle "predicts" their copy number too, the test is detecting phylogenetic
  structure the model has not absorbed and the CotH result means nothing.

#### (f) What each result would mean — written before the result

| Result | What it supports | What it does **not** support |
|---|---|---|
| CotH expanded in plant-associated Mucorales, surviving §2.5c's error bar, §5.2b4's assembly check and (e)'s negative control | a hypothesis that the family has a surface role in plant association, worth a knockout | any claim about binding, about which plant molecule, or about symbiosis vs parasitism |
| No difference by lifestyle | consistent with CotH being a general Mucorales cell-surface family; the mammalian literature does not generalise to plants; the search for plant-relevant surface factors should go to §5.1's other families | that CotH is irrelevant — absence of a copy-number signal is not absence of a role |
| Expanded in clinical lineages only | consistent with the published animal-host story, and a positive argument that CotH is **not** where to look for plant factors | anything about the direction of causation; clinical lineages are also the most sequenced |
| Copy number correlates with contig count | **the measurement is dominated by assembly quality**; no lifestyle conclusion is available | — |

**Pre-declared kill switch:** if the §2.5c annotation control shows between-pipeline copy-number
variance ≥ 50 % of between-lifestyle variance, the lifestyle comparison is reported as
uninterpretable and the measurement is published in its place.

#### (g) The binding partner: does a plant homolog exist at all? — **DEFERRED (2026-11-18)**

> Not being run. This asks whether *Marchantia* has homologs of CotH's **mammalian** receptors,
> which is downstream of a plant-interaction hypothesis that is no longer the reason this family
> is being studied. Kept because it is cheap (two HMM searches against proteomes already
> annotated, §4.5) and because the integrin result is a clean reportable negative whenever
> anyone wants it. **Not** a prerequisite for Q4.

This is the question the plan is most at risk of answering lazily, so it is split into the part
genomics can settle and the part it cannot.

**What a search can settle (and the expected shapes of the answers — these are predictions, not
results):**

| Partner | Search | What is expected, and why it is not the answer |
|---|---|---|
| GRP78 / BiP (HSPA5) — CotH3's reported receptor | HMM/orthology search of both *Marchantia* proteomes, denominator = 20,412 (*M. polymorpha* Tak-1) and 18,566 (*M. paleacea*) proteins | BiP is a conserved ER-lumenal HSP70; plants have characterized BiP genes, so hits are expected. **Finding them settles nothing.** The mammalian mechanism depends on GRP78 being *displayed on the cell surface* under hyperglycaemia and ketoacidosis. A gene's presence is not a surface pool |
| Integrin β1 — CotH7's reported receptor | same | **Plants have no integrin genes.** If the search returns nothing, that is a clean, reportable negative with a denominator, and it removes CotH7's mechanism from the plant setting entirely |

**What a search cannot settle, and must be stated as such:** whether *Marchantia* has an
apoplastic or plasma-membrane-exposed BiP pool. That is a localisation question. The decisive
experiments are named in §10 and **are not being run**.

**The honest framing for the write-up**, stated now so it cannot drift later: CotH is the best
characterized Mucorales cell-surface adhesion family and its copy number varies across the order.
Everything known about its *binding* is from mammalian cell culture. This plan can establish what
the family's distribution and copy number are, with denominators and error bars, and whether
plant homologs of its reported partners exist. **It cannot establish that CotH does anything in a
plant, and no output of it should be worded as though it can.**

---

## 6. From "secreted, small, cysteine-rich" to a defensible shortlist

The first pass showed the structural route is feasible and gave it a calibration. This section
makes the ranking explicit and adds the controls the first pass found by accident.

### 6.1 A pre-declared ranking rule

Written before the data, each component with the weight of evidence it actually carries:

| Component | Source | What it is evidence of |
|---|---|---|
| Secretion consensus ≥ 2/3 (§4.4) | signalp + 3 | that the protein is exported — necessary, not sufficient |
| Confident fold assignment to a known effector scaffold | ESMFold/AF3 + Foldseek, TM ≥ 0.5 and pLDDT ≥ 70 | resemblance; **not function**. The first pass's cerato-platanin call sits at pLDDT 72.0, two points over the cut, and is threshold-fragile |
| **Disorder, separated from prediction failure** | `idp_summary.IDP_fraction` + AF3 MSA depth — **new** | a low-confidence model with high IDP fraction and deep MSA is predicted-disordered; low IDP fraction and shallow MSA is prediction failure. This closes first-pass limitation 6 for the ~70 % of 426 that the AF3 re-run left ambiguous, cheaply, with a table lookup |
| Lineage restriction, three-way recomputed (§2.4) | §5.1 index | lineage specificity, after sampling control |
| Secreted peptidase | `merops` — **new** | a candidate class the tiers miss |
| Host co-fold score | AF3, §6.4 | that two chains can be modelled together. **Not binding** |

Rank is a documented weighted sum with the weights in the code, and the **unweighted component
table is published beside it** so a reader can re-rank.

### 6.2 Bound the shortlist before folding

Fold the top **N = 300** by rule, not everything. Justification from measured throughput, not
guesswork: 640 proteins attempted / 559 folded in **13–20 minutes on one Blackwell**
(`HANDOFF.md` §4). Folding is not the bottleneck; **MSA generation is** (AF3 arrays ran 40 min to
5.5 h per task). N is bounded by MSA cost, so it is set before the MSA step is submitted.

### 6.3 Disorder gets its own answer, not a shrug

For every protein below the fold-confidence cut, emit `IDP_fraction`, AF3 MSA depth and the
pLDDT/pTM, and classify into `predicted_disordered`, `prediction_failure` and `ambiguous` by a
rule declared before looking. Report the three counts with their denominator (426 in the first
pass's unknown set). **Small, cysteine-rich, secreted and unstructured-until-bound is a real
effector phenotype**, so the disordered class is a result, not a residue.

### 6.4 The controls — the part the first pass learned the hard way

| Control | Why | Status |
|---|---|---|
| Positive: AVR-PikD × Pikp-1 HMA (6Q76) | a known effector–target pair; scored **0.91** | exists, reuse |
| Negative: AVR-PikD × an unrelated *Marchantia* protein | scored **0.13**, just below the panel background median of 0.16 | exists, reuse |
| Panel background | 128 fungal × host pairs, median 0.16, p95 0.38, p99 0.67 — the empirical null | exists, extend |
| **Domain-compatibility decoy** | AtCERK1 LysM × *Marchantia* LysM scores **0.84**, above every fungal candidate. This fell out of cross-pairing by accident and is the single most informative row in the first pass. **Make it designed**: for every candidate scored against a host receptor, score a non-partner protein of the same domain family against the same receptor | **new, mandatory** |
| **Length- and composition-matched fungal decoys** | 20 secreted proteins from a non-plant-associated Mucorales genome, matched on mature length and cysteine fraction, against the same host panel. A candidate must beat this distribution, not just the background | **new, mandatory** |
| Shuffled-sequence control | the same candidates with residues shuffled preserving composition | **new**, cheap, run on the top 10 only |

**The reporting rule stays the first pass's**: nothing in the panel reached the positive control,
so candidates are reported as ranked candidates and never as interactions. Any candidate whose
score is matched by its domain-compatibility decoy is **discounted in the table**, as the two LysM
candidates already were.

---

## 7. Phases, SLURM shapes, cost and dependencies

All compute on the UCR HPCC via SLURM. CPU: `-p epyc,intel --constraint="ryzen|broadwell|cascade|sapphire"`.
GPU: `-p short_gpu`, naming a card type — **never a bare `--gres=gpu:1`**, because the k80 (sm_37)
and p100 (sm_60) nodes are below torch 2.8's sm_70 floor (`HANDOFF.md` §1). Large temporaries to
`$SCRATCH`; **reusable intermediates do not go there** — the first pass lost its DuckDB export that
way and now caches under `analysis/mucor_secretome/cache/` (gitignored). Nothing is written to
`/bigdata/stajichlab/shared/...`.

Costs are **estimates**, anchored where possible on measured runs from the first pass.

### Phase 0 — probes and defect checks (prerequisite for everything)

| Step | Shape | Cost |
|---|---|---|
| 0.1 Schema probe: `deeploc`, `wolfpsort`, `merops`, `idp_summary`, `busco_genome` — exact column semantics, value vocabularies | 1 task, 4 cores, 32 GB | minutes |
| 0.2 Re-verify the three known defects on the current build: `is_gpi_anchored` all-TRUE, `SUBPHYLUM` all-NULL, codon flags. **Fixes are in flight; if the tables were rebuilt, the workarounds must be removed, not left in** | same job | minutes |
| 0.3 Close §0.3: the accessions are confirmed versioned (`PF01476.27`, `PF01522.27`); what remains is whether the T1 Pfam branch ever fired for a BFD genome, and audit every other exact-string match on a Pfam accession in `bin/` for the same defect | same job | minutes |
| 0.4 Reuse check on the funnel's corpus cluster sets | `ls` + one read | minutes |
| 0.5 Resolve the first pass's 292-vs-294 Mucoromycota discrepancy (§2.2) and record the reason | 1 task | minutes |

A filtered scan of the 12 GB `pfam.parquet` with 16 threads completed in ~70 s in testing, so
budget 10 minutes wall for the whole phase. **Dependency: everything.**

### Phase A — comparison set and the annotation error bar

| Step | Shape | Cost |
|---|---|---|
| A.1 Build the 454-genome Mucoromycota *s.l.* set + 231 near-outgroups; one-genome-per-species selection by the declared rule; commit `data/comparison_set.tsv` | 1 CPU task | minutes |
| A.2 Run DeepLoc 2.1, WoLF PSORT, MEROPS, AIUPred on JEM1, JEM2, both *Marchantia* proteomes with BFD's invocations (§4.5) | DeepLoc: 1 GPU task on `short_gpu` (name the card) or CPU array; others CPU array of 4 | hours |
| A.3 **The annotation-pipeline control**: re-annotate k = 8 Mucorales assemblies with JEM1/JEM2's funannotate version, run the BFD predictor stack, compute ΔSSP fraction | array of 8, 16 cores, 64 GB, `$SCRATCH` for intermediates | **the expensive step — days of wall clock; funannotate predict is the cost** |
| A.4 Evaluate the §2.5d kill switch and record the outcome | 1 task | minutes |

**Dependency: Phase 0.** A.3 gates the whole of Q1 and should be submitted first because it is the
long pole.

### Phase B — composition (Q1)

| Step | Shape | Cost |
|---|---|---|
| B.1 Secretome + definitions over the comparison set; GPI derived from `predgpi`, completeness from sequence | 1 task, 16 cores, 64 GB | ~1 h |
| B.2 Four-predictor consensus and the disagreement contingency table (§4.3) | same | minutes |
| B.3 The 180-cell threshold-response surface (§3.5) | same | minutes |
| B.4 Covariate models and the Mucor pair's studentised residual with its PI (§2.5b) | 1 task | minutes |
| B.5 Three-way recomputation of every headline statistic (§2.4) | same | minutes |

**Dependency: A.1, A.2, A.4.** If A.4's kill switch fired, B.4/B.5 still run but are reported as
bounded by the measured annotation error, with the error bar on every figure.

### Phase C — restricted families (Q2)

| Step | Shape | Cost |
|---|---|---|
| C.1 MMseqs2 cascaded clustering of the secreted sets at 30/50/70 % (or reuse per 0.4) | 1 task, 32 cores, 128 GB, `$SCRATCH` | hours |
| C.2 Per-cluster HMMs + corpus-wide `hmmsearch` restriction index | array over cluster shards, 16 cores each | **hours to a day; the array width is the lever** |
| C.3 The false-discovery floor on Mortierellales and a matched Ascomycota set (§5.1d) | same shape, 2 more runs | as C.2 |
| ~~C.7 Lifestyle curation (§2.7)~~ | — | **deferred with §5.2e — do not spend the person-days on CotH** |
| ~~C.8 MDE, PGLS, negative-control families (§5.2e)~~ | — | **deferred** |

**Dependency: Phase 0, A.1, A.3 (for the copy-number error bar).**

Lifestyle curation was the only step in this plan measured in *days of a person's time*. It
gated C.8 and nothing else, and both are now deferred, so Phase C is pure compute.

### Phase F — CotH evolution in Mucoromycota (Q4) — **separable, run independently**

Split out of Phase C on 2026-11-18. **Depends on Phase 0 only.** It does not need the comparison
set, the secretome, the consensus axis, the lifestyle table or the host panel, so it can run
before, after or alongside Phases A–D without touching them.

| Step | Shape | Cost |
|---|---|---|
| F.1 CotH: three detections, reconciled — Pfam GA, curated profile, structure-informed (§5.2a) | 1 task | ~1 h |
| F.2 Copy counting with the five rules, incl. split-model detection and the contig-count check (§5.2b) | 1 task | ~1 h |
| F.3 In-house species tree from BUSCO single-copy orthologs of Mucoromycota + Neocallimastigomycota + outgroups | 1 task, 32 cores | hours |
| F.4 CotH gene tree (IQ-TREE, model selection, SH-aLRT + UFBoot) **with a fungal protein-kinase outgroup** — settles whether PF08757 is one family (§5.2c) | 1 task, 32 cores | hours |
| F.5 Reconcile gene tree against species tree; gains/losses per node, §2.4's three-way recomputation on every count | 1 task | ~1 h |
| F.6 Topology test if and only if F.5 shows a nested placement suggesting transfer (§5.2c) | 1 task | ~1 h |
| F.7 Restriction test with denominators, both thresholds reported (§5.2d) | 1 task | ~1 h |

**F.4 gates everything after it.** If the gene tree shows PF08757 is not one family, F.5–F.7 are
recomputed on the redefined set and the (a0) census numbers — including the ones now quoted in
`HANDOFF.md` — are withdrawn rather than corrected. Budget for that outcome; a kinase-like
profile makes it a live possibility, not a formality.

### Phase D — shortlist and interaction (Q3)

| Step | Shape | Cost |
|---|---|---|
| D.1 Rank by §6.1; emit the component table | 1 task | minutes |
| D.2 Disorder triage (§6.3) | 1 task | minutes |
| D.3 Structure for the top N = 300 lacking a model | 1 GPU task, named card — **~15 min per 600 proteins measured** | under an hour |
| D.4 AF3 MSAs for the host panel + candidates + **all decoys** | CPU array; **40 min – 5.5 h per task measured** | the bottleneck |
| D.5 AF3 co-folding: candidates × host panel, plus domain-compatibility decoys, matched fungal decoys, shuffled controls | GPU array | hours |
| D.6 Report against the existing calibration (§6.4) | 1 task | minutes |

**Dependency: B.2 (consensus), C.2 (restriction component of the rank).** Decide and record
whether the AF3 panels are rebuilt from the `transformers` fold calls or keep `fair-esm`
provenance (`HANDOFF.md` §3) — **before** D.4, not after.

### Phase E — write-up

`REPORT.md` beside the code, a thin Quarto wrapper with `{{< include >}}`, figures from committed
tested code as PNG **and** PDF, per `WRITING.md`. The index's "Not published here" updated.

### Non-negotiable engineering, from `ERRORS_FOUND.md`

- Every derived directory is cleared by whatever derives it (`jobs/`, `pairs/`, `pred/`).
- Every resume check is tested against a real output path before being relied on — #7 silently
  disabled resumption for the life of those scripts because AF3 does not lowercase its job names.
- Controls are named by content hash, never by position in a fetch — #5 broke only the controls
  because only the controls were named positionally.
- Any two implementations of the same thing are run side by side and the difference printed; that
  step is nearly free and is the only thing that caught #3.
- A completeness gate: no report is written unless every requested job has a result.
- Every new piece of code gets a fixture test (the directory is at 100 tests).

---

## 8. Falsification controls, by claim

Each row: the claim this plan could make, the control that could kill it, and what outcome does.

| Claim | Control | Outcome that kills it |
|---|---|---|
| "JEM1/JEM2 are typical / atypical of *Mucor* in SSP content" | annotation-pipeline control (§2.5c) on k = 8 assemblies | median between-pipeline \|ΔSSP fraction\| > 0.5 × the across-species IQR |
| same | three-way recomputation (§2.4) | the claim reverses when restricted to one genome per species, or when the largest BioProject is dropped |
| same | threshold grid (§3.5) | holds in a minority of the 180 cells — then it is reported as a sub-range, not a conclusion |
| "Family F is Mucorales-restricted" | identical pipeline on Mortierellales and matched Ascomycota (§5.1d) | the floor run invents families at a comparable rate |
| same | corpus-wide profile search, not GA-cut Pfam (§5.2a2) | members appear outside Mucorales once the GA cut is relaxed |
| "CotH is Mucorales-restricted" (the literature premise) | corpus-wide census with denominators (§5.2a0, §5.2d) | **already failing at the Pfam level** — 25 of 26 Neocallimastigomycota and 100 of 100 Mortierellales genomes carry PF08757. It is rescued only if gene-level and coverage filters remove those hits, or if the tree shows they are not the same family |
| same | protein-kinase outgroup in the gene tree (§5.2c, §5.2d) | the outside-Mucorales hits fall on the kinase side — then the Pfam count was over-recruitment, not homology, and the premise survives for a different reason than it was asserted |
| "CotH copy number differs by lifestyle" | copy number vs contig count and BUSCO duplication (§5.2b4) | copy number correlates with assembly fragmentation |
| same | PGLS on three no-association control families (§5.2e) | lifestyle "predicts" those too |
| same | MDE computed before the test | MDE exceeds the literature effect size → underpowered, reported as such |
| "A plant homolog of CotH's partner exists" | denominator-carrying search of both proteomes | integrin search returns nothing (expected; a clean negative). A BiP hit is **not** a positive for the mechanism |
| "Protein P is a candidate host-interacting protein" | domain-compatibility decoy (§6.4) | the decoy matches or beats P — as LysM×LysM at 0.84 beat every fungal candidate |
| same | length/composition-matched fungal decoys | P sits inside the decoy distribution |
| same | the existing positive control at 0.91 | P is far below it — which was true of all five first-pass hits, and is why they are candidates |
| "The low-confidence folds are genuinely disordered" | `IDP_fraction` + MSA depth (§6.3) | low IDP fraction with shallow MSA → prediction failure, not disorder |
| Any consensus-secretion claim | focal-vs-background disagreement rate (§4.3) | JEM1/JEM2 disagree more than the BFD distribution → artefact in the local predictor runs |

---

## 9. Limitations

1. **Nothing here is measured at the bench.** Every output is a prediction or a count.
2. **Q1 may be unanswerable with BFD as deposited**, and the plan is built to find that out rather
   than to avoid it (§2.5d).
3. **JEM1 has no conspecific genome.** Strain-specific and species-specific cannot be separated
   for it, at any point, by any method in this plan.
4. **The isolation host is unrecorded** in the MorenoLab project, so both liverworts are treated
   symmetrically. Knowing it would let them be weighted and would roughly halve the host panel.
5. **Lifestyle labels are hand-curated from literature** and will be `unknown` for many species.
   The `unknown` class is reported, not dropped, but it is still missing data and it is almost
   certainly not missing at random — poorly studied species are poorly labelled.
6. **The CotH lifestyle test is likely underpowered**, with plant-associated Mucorales species in
   the single digits against 36 genera. §5.2e requires this to be quantified before the test.
7. **The mammalian CotH literature may simply not transfer**, and the family may not even be
   Mucorales-restricted (§5.2a0). A negative result on both counts is a likely outcome and is
   planned for as a reportable result (§5.2d, §5.2f).
8. **The §5.2a0 census is one Pfam scan.** It counts proteins rather than genes, applies no
   domain-coverage filter, no deduplication and no phylogeny, and is uncorrected for everything in
   §2. It is strong enough to change the plan's framing and far too weak to be quoted as a result.
9. **Copy number is an annotation statistic.** Split models, collapsed repeats and predictor
   differences all move it, and §5.2b measures rather than assumes that.
10. **`merops`, `deeploc`, `wolfpsort` and `idp_summary` are being used for the first time here**
    and their column semantics are Phase 0 probes. Thresholds in §4.2 are provisional until then.
11. **Three BFD defects are being worked around, and the fixes are in flight.** If the tables are
    rebuilt mid-analysis the workarounds become wrong in the other direction. Phase 0.2 re-checks
    them; whoever runs this must re-check again before publishing.
12. **The first pass's AF3 panels were built on superseded fold calls** (`fair-esm`, not
    `transformers`). Either provenance is acceptable; leaving it unstated is not.
13. **Glomeromycota is a weak comparator** despite being the most relevant lineage: 60 genomes,
    31 species, 58 of 60 GenBank-only, and obligate biotrophs whose gene content is shaped by
    unculturability as much as by symbiosis.

---

## 10. The decisive experiments, which are **not** being run

Named with what each would settle, per `WRITING.md` rule 5. None is in this plan.

| Experiment | What it would settle | Why it is not here |
|---|---|---|
| **RNA-seq of JEM1/JEM2 in contact with *Marchantia* thalli**, against axenic culture, 3 biological replicates per condition | Which candidates are *expressed in planta*. This is the single largest gap: every candidate in the first pass and in this plan is a sequence prediction with no expression evidence. Induction in contact would promote a handful of candidates above all the structural ranking in §6 | Requires the co-culture system and the isolation host, which is unrecorded |
| **Apoplastic / plasma-membrane proteomics of *Marchantia***, with an intracellular-contamination marker | Whether a surface-exposed BiP pool exists in a plant at all — the question on which CotH3's mechanism turns, and which no genome search can answer (§5.2g) | Not a genomics experiment |
| **Yeast two-hybrid or pull-down of the top AF3 candidates against their host partners** | Whether `JEM1MUCOR_000290` (the chitinase-domain protein that reached three of eight host LysM receptors at iPTM 0.73 / 0.67 / 0.54) binds anything. The iPTM ranking cannot distinguish binding from co-foldability, and the LysM×LysM control at 0.84 shows how large that gap is | Bench work |
| **CotH knockout or knockdown in *M. circinelloides***, scored on *Marchantia* colonisation | Whether CotH does anything in a plant. *M. circinelloides* is transformable, so this is the feasible decisive test — it is simply not a computational one | Bench work |
| **Chitin-oligomer binding assays (CO4 vs CO7)** on the chitin-active candidates and the Mp/Mpa LysM ectodomains | Whether the deacetylase/chitinase masking hypothesis in first-pass §8 holds. The AF3 ligand runs (16 of 32 at pTM ≥ 0.8, 4EBZ positive controls at 0.86/0.77) are consistent with it and establish nothing about affinity | Bench work |
| **Sequencing 5–10 additional *M. hiemalis* strains** | Would make JEM1's ~105 genome-specific secreted proteins interpretable, by separating strain-specific from species-specific. Currently impossible for lack of a conspecific reference | Sequencing cost; it is the cheapest thing on this list and the one most worth proposing |
| **An ML tree of the CotH family including bacterial CotH** | The family's origin and whether the Mucorales set is monophyletic. Partly addressed by §5.2c within fungi, but BFD holds no bacteria | Out of BFD's scope; would need external sequence retrieval. **Reconsider under Q4** — now that the family's history *is* the question rather than a side-check, the bacterial outgroup is the difference between "expanded within Mucoromycota" and "acquired by Mucoromycota". One UniProt/NCBI fetch of bacterial CotH, which is retrieval of public sequence and carries none of §11's constraints. **Decide before F.4, since it changes the alignment** |

---

## 11. Data handling and hard constraints

- **MorenoLab sequences never leave the cluster.** JEM1 and JEM2 are a collaborator's unpublished
  data. No hosted model endpoint, no web API, no online BLAST, no public structure server, no
  third-party folding service, at any step. ESMFold and AlphaFold 3 run on on-prem GPUs against
  local databases; DeepLoc runs from the lab-built apptainer SIF; MSAs are built locally. The
  first pass already holds to this (its ESMFold run is single-sequence specifically so that
  "nothing leaves the cluster") and nothing in this plan relaxes it. The same handling is applied
  to the public *Marchantia* proteomes, because a per-sequence exception is a rule that will
  eventually be got wrong.
- All compute via SLURM on the UCR HPCC; nothing runs on a laptop. CPU `-p epyc,intel
  --constraint="ryzen|broadwell|cascade|sapphire"`; GPU `-p short_gpu` with the card type named.
- Large temporaries to `$SCRATCH`; reusable intermediates to `analysis/mucor_secretome/cache/`.
- **Nothing is written into `/bigdata/stajichlab/shared/...`.** BFD tables and the MorenoLab
  genomes are read-only inputs. Defects found there are reported upstream, as the first pass did
  for `is_gpi_anchored` and the codon flags — and §0.3 adds a third to report.
- Code in the repo, documented, with a fixture test per new piece.

---

## 12. Deliverables, and what they are not

### 12.1 Tables

| File | Contents | Denominator it carries |
|---|---|---|
| `data/comparison_set.tsv` | the 454 + 231 genomes, with the one-per-species flag and the selection rule's inputs | 23,680 BFD rows |
| `data/lifestyle.tsv` | one row per species, citation per row, `unknown` included | 204 Mucoromycota *s.l.* species |
| `out/secretome_composition.tsv` | per genome and per species: proteome size, secretome size and fraction, `C-SSP`, `C-CRP`, `C-EFF`, `C-CHI`, `C-PLA`, `C-PEP`, agreement scores | proteome size, per row |
| `out/threshold_surface.tsv` | every headline statistic × 180 threshold cells | 180 |
| `out/annotation_control.tsv` | ΔSSP fraction per re-annotated assembly, and the kill-switch verdict | k = 8 |
| `out/restricted_families.tsv` | ranked candidate Mucorales-restricted families, with the false-discovery floor from the control runs | 36 genera / 106 species / ~23,500 genomes |
| `out/coth_copies.tsv` | `copies_strict`, `copies_merged`, assembly covariates, per genome and per species | 275 Mucorales genomes / 106 species |
| `out/shortlist.tsv` | ranked candidates with every component column unweighted | N = 300 folded, from the full candidate set |
| `out/decoy_controls.tsv` | domain-compatibility, matched-fungal and shuffled controls beside the candidates they qualify | per candidate |
| `out/disagreement_set.tsv` | candidates dropped by the §4.4 consensus rule, with their four calls | — |

### 12.2 What none of them is

**No binding is measured. No expression is measured. No activity is predicted.** The shortlist is a
list of proteins worth assaying, in an order — the wording the first pass settled on, and the
wording this plan keeps. A reviewer who reads a ranking as a finding has been misled by the
document, not by the method, and §12.1's denominators exist to make that harder.

### 12.3 One note for whoever runs this

`WRITING.md`, which this plan follows, is **not on `main`**. It exists on branch `upo-page` at
commits `92cc8e8` and `70f6dde`, together with a `CLAUDE.md` that is also absent from `main`.
Anyone looking for the repo's writing standard at the root of a `main` checkout will not find it.
Worth merging, or at least worth knowing before the next report is written against a standard that
is not where it is said to be.
