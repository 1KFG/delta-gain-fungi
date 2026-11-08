# Independent review — enzymes_UPO PLAN.md

Date: 2026-11-07
Scope: `PLAN.md`, `CENSUS.md`, `REFERENCES.md`, `sql/{census,hits,novel_lineages,sp}.sql`,
`data/{early_pf01328_hits.parquet,early_pf01328_hits.fa.zst,motifs.tsv}`
Nothing was edited; this file is the only thing written.

## Overall assessment

The census work is solid and the plan is unusually honest about what it has and has not verified,
but the plan is **structurally flawed in its framing, not merely incomplete**. Two things break it:

1. The UPO-vs-CPO call is built on the EGD/EHD motif dichotomy, which the literature has already
   shown is *not* a UPO/CPO discriminator, so the plan's central classification step cannot work
   as written.
2. The project has no engagement with Zamocky et al. 2022 (Antioxidants 11:1011, PMC9138132),
   which already publishes a Mucoromycota-basal heme-thiolate peroxidase clade, a distinct
   Glomeromycota short (~230 aa) subfamily, and Chytridiomycota HTPs — prior art on exactly this
   question.

The search, QC and stratification machinery is good. The scientific positioning and the calling
logic need replacing before anything is launched.

---

## 1. The EGD/EHD conflict is already resolved, and the resolution kills the `Type` step

Refs: PLAN §5 row "Type"; §9.1, §9.2; CENSUS §"Crude motif scan"; REFERENCES §3.3

The plan treats EHD-vs-EGD as an unresolved label conflict to be resolved against structures
before any hit can be called a UPO, and §9.2 worries that because ~95% of early-diverging hits
carry EHD, "the final set could be small." Both premises are wrong.

Faiza et al. 2019a (BMC Evol Biol, PMC6417270) states it explicitly: EGD is Glu122-Gly123-Asp124
of AaeUPO; in CPO the Gly is replaced by His, giving EHD; the summary patterns are
`-PCP-EGD-R--E` (AaeUPO) and `-PCP-EHD-E` (**MroUPO** and CPO) — "MroUPO was therefore placed
among the CPOs for analysis." MroUPO is a crystallised, well-characterised short-type
peroxygenase. That paper's own data therefore show EHD does not mean "not a UPO": the "CPO-type"
label is *defined by* the motif, and MroUPO is the counterexample proving motif != function.

Zamocky 2022 goes further: the distal motif is I/V/L-**E-H-D**-X-S, and His105 there is
**essential for Mn2+ binding** in CfumagHTP, while AaeAPO has Gly at that position and binds
Mg2+. EHD vs EGD is a divalent-metal-site difference, not a functional class boundary.

**Fix.** Delete the `Type` row as a UPO gate. Replace it with the actual discriminator, which
Faiza 2019a names and the plan is missing: the distal **acid-base pair**. AaeUPO uses
Arg189-Glu196; LfuCPO uses His105-Glu183. The UPO-vs-CPO call is the residue at the
Arg189-equivalent position (**Arg -> peroxygenase-like; His -> CPO-like**) plus the Glu. Report
EHD/EGD as a metal-site annotation, not a class call.

This also closes §10's "exact residue number of the acid-base Glu in AaeUPO — not verified": it
is **Glu196, paired with Arg189** (Faiza 2019a). Reading 2YOR is no longer needed to code §5.

## 2. No engagement with the direct prior art — a reviewer attacks this first

Refs: PLAN §0, §1; REFERENCES throughout

REFERENCES.md cites Kinner 2021 and the two Faiza 2019 papers, and flags a "Zamocky et al. 2025"
claim as unverified. The relevant paper is Zamocky et al. **2022**, "Evolution of Heme
Peroxygenases: Ancient Roots and Later Evolved Branches" (Antioxidants 11:1011, PMC9138132), and
it is absent from REFERENCES.md. It:

- places **Mucoromycota HTPs basal to the whole superfamily** ("deep ancient evolutionary origins
  in Mucoromycota and Chytridiomycota"), with ~10 representatives of an ancestral mucoromycetous
  clade;
- describes a **distinct Glomeromycota subfamily of short sequences averaging ~230 aa** — which
  CENSUS §4 independently reproduces (Glomeromycota median 225/233 aa). A published clade has
  been re-measured without citation;
- places two Chytridiomycota (*Spizellomyces*) HTPs **inside the oomycete family**;
- verifies the claim flagged as unverified: Mucoromycota, Glomeromycota and Chytridiomycota are
  "still under-represented in sequence databases," and no heme-thiolate peroxidase is known from
  Glomeromycota at the protein level despite its ancient roots;
- uses Family I (short) / Family II (long) **plus** a mixed fungal family and an oomycetous
  family — a four-way scheme Decision 1 does not accommodate;
- is built on **174 full-length sequences from RedoxiBase + GenBank/UniProt**.

**Consequence for §0.** The novelty claim as written ("which UPO genes exist in early-diverging
fungal genomes") is partly already answered. The defensible contribution is *scale and annotation
provenance*: 807 PF01328 proteins and 309 from BFD_only genomes, against ~15-20 early-diverging
sequences in Zamocky's 174-sequence set. Reframe §0 as "a ~40x expansion of the early-diverging
HTP sample, and the first with per-genome gene-model provenance," and state Zamocky's clades as
the baseline the tree must recover — a positive control currently absent for the target phyla
(§3.4's control is Dikarya-only).

Also: Faiza 2019a already includes *Absidia repens* and *Choanephora cucurbitarum* (both
Mucoromycota) in its CPO selection analysis. Check whether those two species' sequences are among
the 807 and how they were classified. A disagreement with a published call is either a bug or a
result; worth knowing which before publishing.

## 3. The UPObase fallbacks all miss the live alternative

Refs: PLAN §3.1

UPObase is dead (REFERENCES §5 measured this properly). The three listed fallbacks are email the
authors, Wayback, or rebuild from UniProt. **RedoxiBase** (peroxibase.toulouse.inra.fr, successor
to PeroxiBase; Redox Biol 2019, PMC6597783) is live and is the curated source Zamocky used for
this superfamily, and Zamocky's **Table S1 (174 curated full-length sequences)** is a downloadable
supplement.

Together they give a curated reference set, a classification scheme, and the prior-art baseline in
one step. This single action closes §9.1 (classification), §9.2 (is CPO-like a deliverable), §9.4
(UPObase replacement) and §10's unverified Zamocky claim. **Highest-value change in the plan.**

## 4. Tier 1 as defined is nearly empty

Refs: PLAN §8.4

Tier 1 = QC-clean + Cys + Glu + UPO call + **BFD_only** + **has SP**. Against CENSUS §4:

| Phylum | BFD_only hits | SP rate | Tier 1 ceiling |
| --- | --: | --: | --: |
| Mucoromycota | 177 | 1.1% | ~2 |
| Glomeromycota | 32 | 0.0% | 0 |
| Chytridiomycota | 4 | 50% | 2 |

Tier 1 across the three core phyla is **~4 proteins before any QC or motif filter**. The only
well-populated Tier 1 is Basidiobolomycota BFD_only (51.9% x 77 ~= 40), an "extended" phylum §9.6
has not decided to include. The ranking criterion and the target set are incompatible, and the
headline deliverable collapses.

This is expected, not bad luck. Kinner's Family I (short-type) is "mainly without predictable
signal peptides," and Family II (long-type, SP-bearing) is **restricted to Ascomycota and
Basidiomycota**. The plan mines a region of the tree where, by construction, long-type SP-bearing
UPOs do not exist. Gating Tier 1 on SP gates on a Dikarya trait.

**Fix.** Decouple the axes. Novelty tier (BFD_only / NCBI-but-not-annotated-as-UPO / already
annotated) x expression-feasibility score, where feasibility is scored on **mature-domain
integrity** (full HMM coverage, Cys + Glu + Arg present, no internal stop, start/stop intact),
free-vs-paired cysteine count, predicted N-glycosylation site count, absence of TM helices, and
length class — *not* native SP presence. Signal-peptide swapping is routine in UPO heterologous
expression (Pullmann & Weissenborn 2021, cited in REFERENCES §6 only as an accession source, not
as a method), so native SP should be a reported attribute rather than a filter.

Related gap: §8.4 ranks for heterologous expression but **never names a host**. Yeast vs *E. coli*
changes every criterion (disulfides, glycosylation). Decision 5 commits to expression within a
year; fix the host before defining the ranking.

## 5. The absence result is overstated — wrong denominator

Refs: PLAN §1 line 67, §8.6; CENSUS §2 vs §3

CENSUS §2 gives total genomes; CENSUS §3 gives *genomes with Pfam rows*. They differ:

| Phylum | Genomes (§2) | With Pfam rows (§3) | Gap |
| --- | --: | --: | --: |
| Kickxellomycota | 190 | 177 | **13** |
| Entomophthoromycota | 12 | 9 | **3** |
| Chytridiomycota | 89 | 83 | 6 |
| Neocallimastigomycota | 26 | 25 | 1 |
| Microsporidia | 163 | 133 | **30** |

PLAN line 67 states "Zero PF01328 proteins ... in Mortierellomycota (100 genomes),
Kickxellomycota (190) ... Entomophthoromycota (12)" — but the zero was measured over the smaller
searched set. The Entomophthoromycota absence rests on 9 of 12 genomes; Kickxellomycota on 177 of
190. Fix the sentence and add a `searched / not-searched` column to §8.6, otherwise the absence
note is attackable in one line. Mortierellomycota 100/100 and Blastocladiomycota 14/14 are clean —
those absences hold.

## 6. jackhmmer is the wrong instrument for pass 2

Refs: PLAN §4.2, §7

Three single seeds, 3 iterations, inclusion E 1e-3, against 10.17 M proteins, in a family where
PF01328 alone returns **97,520 proteins BFD-wide** and UPObase's Pog-g bin holds 1,509 members
with no motif beyond the heme Cys. jackhmmer will converge onto the Pog cloud and recruit
indiscriminately; "convergence" will be convergence onto noise, and the output will not be
interpretable without a decoy-calibrated threshold.

**Fix.**

- Build the Family I and Family II HMMs from curated multi-sequence alignments
  (RedoxiBase / Zamocky Table S1), not single seeds, and set the hmmsearch threshold by ROC
  against the 13 characterized references as positives and a sampled Pog/CPO decoy set as
  negatives.
- Restrict jackhmmer to the five zero-hit phyla, which is what §2's "Negative / sensitivity check"
  row actually says — §4.2's "all 813 proteomes" contradicts it.
- §4's dismissal of HHblits ("adds little here if jackhmmer converges") is a non-sequitur;
  jackhmmer convergence says nothing about HMM-HMM reach. Run profile-profile on the pass-2-only
  hits and on early-diverging proteins that fail PF01328 GA.

**Missing search route — the FunBGCeX trick, cited nowhere.** FunBGCeX's transferable method is:
take proteins with *no* domain hit, cluster all-vs-all BLASTp at E 0.01, align with Clustal Omega,
hmmbuild into new profiles (572 domainless proteins -> 59 new HMMs). §3.2/§6 builds custom HMMs
only from clades that form in the tree of hits already in hand — i.e. it only profiles what was
already found. Add the domainless arm: cluster early-diverging proteins that Foldseek or HMM-HMM
flags but that fail PF01328 GA, and build profiles from those clusters. This is the only route in
the plan to members that sequence HMMs genuinely miss.

**Foldseek should not be "optional" (§4.2).** Structure is the only orthogonal signal available.
Make it first-class on the shortlist and on all pass-2-only hits (ESMFold is affordable at that
scale), with 2YOR / 5FUJ / 1CPO as queries and a stated TM-score threshold.

## 7. Structure prediction is budgeted but never specified as a step

Refs: PLAN §7 row; absent from §5, §6, §8

§7 sizes "Structure prediction | shortlist only, 50-200 proteins | hours on one GPU." It appears
in no methods section and no deliverable, and has no acceptance criterion. Given self-generated
gene models, structure is the strongest artifact control available.

**Fix.** Specify it in §6: ESMFold/AF2 on the shortlist, then (a) Foldseek TM-score vs
2YOR/5FUJ/1CPO, (b) proximal Cys loop conformation and Cys-Sgamma accessibility consistent with a
thiolate ligand, (c) distal pocket residue set including the Arg/His position. Report pass/fail
counts per phylum, as §6's last bullet already promises for the sequence-level checks.

## 8. The novelty metric measures the right thing but is implemented poorly

Refs: PLAN §5 row "Distance"

Credit where due: "Max % identity to any of the **13 references**" is distance to nearest
*characterized* homolog, which is the informative axis. Three problems:

- **DIAMOND is wrong at this divergence.** Local identity without coverage is manipulable (a
  30-residue local hit at 60% identity outranks a full-length 28% match). Use global alignment
  (needle/parasail), require >=80% coverage of the reference, and report identity *and* coverage.
- **% identity will not discriminate.** All 13 references are Dikarya; every early-diverging hit
  will sit in a flat 20-30% band. Make tree **patristic distance** to the nearest characterized
  enzyme the primary continuous measure, with identity as a secondary descriptor.
- **Report both axes separately.** Distance to nearest characterized enzyme (biochemical-novelty
  axis) *and* distance to nearest RedoxiBase/Zamocky entry (sequence-novelty axis). With only the
  first, "biochemically unexplored" cannot be distinguished from "already sequenced, never
  assayed."

## 9. RNA-seq support is the load-bearing unknown and is unmeasured

Refs: PLAN §6, §10

§6 notes `species_with_real_rnaseq.txt` lists 52 species BFD-wide and that it is "not checked how
many are early-diverging." For a project whose entire vulnerability is self-generated gene models
in exactly the phyla where gene-model confidence is lowest, this is the first number to have. 52
BFD-wide across 23,680 genomes strongly suggests the early-diverging count is near zero — in which
case §6's RNA-seq bullet is not a control at all and must be replaced.

**Fix.** Measure it now (one query). If ~0: (a) pull public SRA RNA-seq for shortlist species and
align, or (b) **use miniprot as an independent gene-caller check, not just as rescue**. §4.3
already runs miniprot against the assemblies; compare its exon boundaries against the BFD model at
every candidate locus. Concordance between an independent homology-based projection and the
ab initio/evidence model is cheap, orthogonal evidence already being paid for in compute. Promote
it from "rescue" to a QC row in §5.

## 10. The Glomeromycota truncation worry is already answered — favourably

Refs: PLAN §6 bullet 1 vs CENSUS §4

§6 asks whether short Glomeromycota hits are biology or truncation. CENSUS §4 shows Glomeromycota
hits at 90.2% / 93.8% HMM coverage >=80% of a 186-aa model, with median length 225-233 aa. A
complete domain in a short protein is not a truncated long protein. Combined with Zamocky's
published Glomeromycota subfamily averaging ~230 aa, this is biology, and can be stated now. Keep
the cross-strain locus comparison as confirmation, but drop it from the critical path.

## 11. Chytridiomycota: the contamination check is testing a published hypothesis

Refs: PLAN §6 bullet 3, §9.3

Zamocky places *Spizellomyces* chytrid HTPs **inside the oomycete family**. The Chytridiomycota
set is 184 NCBI proteins with 56.5% SP — the highest SP fraction in the data and unlike the other
early-diverging phyla. Whether that is HGT, contamination, or genuine shared ancestry with
stramenopile peroxygenases is a first-order question, and §6's test (neighbour genes' best hits in
the same phylum, GC, codon usage) is the right test.

**Fix.** Reframe it as a hypothesis with a stated prediction rather than a QC footnote, and
include oomycete HTPs in the §8.3 tree outgroup, which currently lists only "CPO/Pog outgroup
sequences." This also answers §9.3: Chytridiomycota's value is this question, not novel sequence
count.

## 12. Reproducibility and internal consistency

- **The census is not reproducible from the repo.** `sql/census.sql` and `sql/hits.sql` both read
  `${SCRATCHPAD}/pf01328.parquet` and `${SCRATCHPAD}/pfam_genomes.parquet`, and the query that
  *produced* those from `pfam.parquet` is not in `sql/`. PLAN line 252 acknowledges the path
  problem but not the missing step. Add the extraction query — two lines, and the first thing a
  reviewer asks to see.
- **Three different denominators.** `novel_lineages.sql` defines "has BFD gene models" from
  `gene_info.parquet` (23,451); `hits.sql` uses `gene_proteins.parquet`; CENSUS totals use
  `samples.parquet` (23,680 / 827 / 813). Pick one definition and use it everywhere, or the
  per-phylum fractions in §6's artifact-rate table will not be comparable.
- **CENSUS §4 covers only 777 of the 807 hits.** Missing rows: Basidiobolomycota NCBI (9),
  Zoopagomycota NCBI (12), Monoblepharomycota (4), Aphelidiomycota (3),
  Calcarisporiellomycota (2). The 807 total and the per-phylum/annotation breakdown were verified
  directly from `data/early_pf01328_hits.parquet` — the parquet is complete and matches CENSUS §3
  exactly, so this is a reporting omission, not a data problem. But it omits SP rates for the
  Basidiobolomycota NCBI set, one of the groups the expression triage depends on.
- **`sp.sql` line 3** joins `busco_genome` and `asm_stats` on ASMID with no dedup. If BUSCO has
  multiple rows per ASMID (different lineage datasets), the medians and the `busco_lt80` counts in
  CENSUS §2 are wrong. Worth one `count(*) vs count(distinct ASMID)` check.
- **`sp.sql` `new_species_vs_ncbi`**: the `SPECIES not in (select SPECIES from j where ann='NCBI')`
  subquery is global across all phyla, not per-phylum. Harmless here, but the column answers a
  slightly different question than CENSUS §2's header claims.
- **The Decisions section contradicts §9.** PLAN line 8 says the Decisions "replace the matching
  open questions in section 9," but §9.1, §9.4 and §9.5 still ask questions Decisions 1, 4 and 5
  settle. Delete them, or the document disagrees with itself on settled points.
- **REFERENCES §3.1 is wrong.** It concludes "'113 putative UPOs from >800 genomes' is probably a
  mis-citation ... 113 is the number of *seed* sequences." It is not a mis-citation. Faiza 2019a
  (PMC6417270) reports 113 putative UPO sequences in 35 species from >800 genomes (812 Ensembl
  species/strains); UPObase (Faiza 2019b) then used those same 113 as seeds to reach 1,948. Both
  numbers are correct and consistent. Correspondingly, PLAN §3.4's instruction "**Do not use '113
  from >800 genomes' as the target**" should be softened: 113/812 ~= 0.14 UPO per genome *at
  UPO-strict stringency on a Dikarya-dominated set* is a legitimate order-of-magnitude precedent
  and the right sanity check on yield. The 807/827 figure is a family-level net (PF01328), not a
  UPO call — expect a one-to-two-order-of-magnitude drop after strict calling, and say so in §8 so
  the yield is not read as over-promising.

---

## What is genuinely good

- **§5 line 165, "Never pooled across phyla."** Exactly right, and most papers in this space get
  it wrong.
- **The novelty claim in §0 is stated narrowly and correctly** — BFD_only vs "new as UPO call" are
  separated, which is the distinction that makes a delta-gain application defensible.
- **§10 "Claims not verified" is the best part of the document.** Every item checked that was
  flagged unverified was genuinely unverified, and nothing checked that was asserted as measured
  was wrong. The parquet reproduces CENSUS exactly; the 827/826+1 genome arithmetic reconciles via
  the `annotated_folder_missing_protein` genome noted in CENSUS §1.
- **§5's catalytic checks are alignment-based, not substring-based**, with the crude substring scan
  explicitly labelled "not a UPO call" (CENSUS §4). This is the single most common failure mode in
  motif mining, and it was avoided.
- **§6 is a serious artifact-control section** — intron-position concordance, cross-strain locus
  comparison, GC/codon contamination checks, and an explicit per-phylum flagged-fraction report.
  The per-phylum artifact-rate bullet is what makes the whole thing auditable.
- **§3.4's requirement that all 13 references pass the pipeline** is the right smoke test.
- **REFERENCES §6 is careful** — catching that MroUPO's UniProt entry is a 67-aa fragment, and
  excluding the P450 BM3 / SyoA false positives from RCSB text search, is exactly the diligence
  this reference set needs.

---

## Prioritized changes

1. **Get Zamocky et al. 2022 (PMC9138132) and its Table S1 (174 curated sequences), plus
   RedoxiBase (peroxibase.toulouse.inra.fr).** One action closes §9.1, §9.2, §9.4, the
   dead-UPObase problem, and the unverified-Zamocky item, and supplies the prior-art baseline the
   tree must recover. Do this before anything else.
2. **Replace the `Type` (EGD/EHD) call with the distal acid-base pair call:** Arg vs His at the
   AaeUPO Arg189 position, plus Glu196. Report EHD/EGD as a metal-binding-site annotation. Delete
   §9.2 as a blocker — it is answered.
3. **Reframe §0's novelty claim** against Zamocky 2022: the contribution is scale (807 vs ~15-20)
   and gene-model provenance (309 BFD_only), not first discovery. Add recovery of his
   Mucoromycota-basal and Glomeromycota-short clades as a positive control for the target phyla.
4. **Measure how many of the 52 RNA-seq species are early-diverging, today.** If ~0, promote
   miniprot from rescue to an independent gene-model concordance check in §5, and plan SRA pulls
   for the shortlist.
5. **Rebuild the tiering (§8.4):** novelty axis x expression-feasibility axis, with feasibility
   scored on mature-domain integrity, cysteine pattern, glycosylation sites and TM absence — not
   native SP. Fix the expression host first. As written, Tier 1 for the three core phyla is ~4
   proteins.
6. **Fix pass 2 (§4.2):** curated-alignment Family I/II HMMs with decoy-calibrated thresholds;
   jackhmmer restricted to the five zero-hit phyla; HHblits on pass-2-only and PF01328-failing
   hits; add the FunBGCeX domainless-clustering arm; promote Foldseek from optional.
7. **Correct the absence claims (§1 line 67, §8.6)** to the searched denominators, with a
   `not-searched` column.
8. **Correct REFERENCES §3.1 and soften §3.4's "do not use 113"** — the citation is sound, and
   113/812 is a legitimate yield sanity check.
9. **Specify structure prediction as a §6 step with acceptance criteria**, not just a §7 budget
   line.
10. **Reproducibility cleanup:** add the `pf01328.parquet` extraction query to `sql/`; unify the
    genome denominator; complete CENSUS §4 to all 807 hits; check `busco_genome` for duplicate
    ASMID rows; delete the §9 questions the Decisions section already settles.

---

## Sources consulted

- Faiza et al. 2019a, BMC Evol Biol — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6417270/
- Zamocky et al. 2022, Antioxidants 11:1011 — https://pmc.ncbi.nlm.nih.gov/articles/PMC9138132/
- RedoxiBase — https://peroxibase.toulouse.inra.fr
- RedoxiBase paper, Redox Biol 2019 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6597783/
- UPObase, Database 2019 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6902001/
- Kinner et al. 2021, Front Bioeng Biotechnol — https://pmc.ncbi.nlm.nih.gov/articles/PMC8293615
- MweUPO structural paper — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11925332/
- Zamocky et al. 2025 (TU Delft record) — https://repository.tudelft.nl/record/uuid:4e78a3a7-94cd-4c19-afec-36c38368ae28
- FunBGCeX, Nat Commun 2024 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11109268/
