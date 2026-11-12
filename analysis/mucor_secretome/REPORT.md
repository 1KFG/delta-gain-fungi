# Candidate effectors of two *Mucor* strains isolated from *Marchantia*

Status: in progress. Last updated 2026-11-12.
Code: `analysis/mucor_secretome/`. All analysis ran on the UCR HPCC.

## The question

Two *Mucor* strains, JEM1 and JEM2, were isolated from *Marchantia* (MorenoLab,
`shared/projects/MorenoLab/Mucor_genomes`). The practical questions are:

1. **Which of their secreted proteins are candidate effectors?**
2. **Could any of them interact with host proteins?**

Species identity was a secondary question that turned out to matter, because it determines
which genomes count as "relatives" when measuring how lineage-restricted a candidate is.

## Summary of findings

| # | Finding | Confidence |
|---|---|---|
| 1 | JEM1 is *Mucor hiemalis*, JEM2 is *M. circinelloides* — `samples.csv` was right | High (ITS, controls passed) |
| 2 | BFD's only "*M. hiemalis*" genome (`GCA_977110905.1`) is **not** *M. hiemalis* | High (its own ITS) |
| 3 | JEM2's raw assembly is ~2,233 contigs of *Fusarium*; the cleaned assembly used for annotation is free of it | High (positive control fired) |
| 4 | 1,344 secreted candidates across the two strains; 99 are chitin/chitosan-active | High |
| 5 | Chitin-active candidates fold confidently to deacetylase, LysM and cerato-platanin scaffolds | High (TM ≥ 0.93, pLDDT ≥ 71) |
| 6 | A protein with **no Pfam domain** is a cerato-platanin by structure | Moderate (TM 0.84, pLDDT 71) |
| 7 | JEM1 has ~8× more strain-specific secreted proteins than JEM2 after controlling for reference sampling | Moderate — see caveat |
| 8 | Both *Marchantia* species carry comparable LysM receptor complements | High |

---

## 1. Species identity

`samples.csv` listed JEM1 as *M. hiemalis* and JEM2 as *M. circinelloides*. An earlier fastANI
run appeared to contradict this: JEM2 matched *M. circinelloides* at 99.8 %, but **JEM1 had no
relative above 80 %**, its best hit being BFD's *M. hiemalis* at 79.5 %.

ITS against UNITE v10 resolved it, with both BFD assemblies included as controls — which is
what made the result interpretable.

| Assembly | Top ITS hit | Identity | Verdict |
|---|---|---|---|
| **JEM1** (cleaned) | *Mucor hiemalis* | **99.3 %** | *M. hiemalis* |
| **JEM2** (cleaned) | *Mucor circinelloides* | **98.6 %** | *M. circinelloides* |
| BFD *M. circinelloides* (control) | *Mucor circinelloides* | 98.8 % | control passes |
| BFD *M. hiemalis* (control) | **"Mucorales sp."** | 98.5 % | ***M. hiemalis* not in its top 5** |

**The mislabelled genome is BFD's, not MorenoLab's.** JEM1's low ANI reflects a bad reference,
not an unusual query. No correction to `samples.csv` is needed; the public label on
`GCA_977110905.1_gzMucHiem1` is the one worth questioning.

Two bugs had to be fixed before this ran at all:

* UNITE's QIIME release keeps lineages in `sh_taxonomy_*.txt`, **not** in the FASTA headers, so
  grepping the FASTA for "Mucoromycota" matched nothing. `bin/unite_tax.py` joins the two.
* ITSx found **no ITS in any of the four assemblies**, controls included: the rDNA array is
  collapsed or broken in all of them. The MorenoLab ITSx result that an earlier version of this
  analysis would have fallen back on is not ITS either — JEM1's "full" sequence is SSU, reported
  with a negative ITS length. The method was changed to blast every Mucoromycota UNITE reference
  against each assembly, which needs no intact ITS.

## 2. Contamination

The ITS run found a **100 %-identity *Fusarium graminearum* ITS in JEM2's spades assembly**.
FCS-GX had already flagged 2,233 of its 13,555 contigs as *Fusarium venenatum* and marked them
EXCLUDE; the cleaned assembly has 2,584 contigs. Since the proteome feeding effector calling is
predicted from the cleaned assembly, a surviving contaminant would be scored as a *Mucor*
effector. Verified directly rather than assumed:

| Assembly | Ascomycota ITS hits ≥ 97 % over ≥ 300 bp | Contigs hit |
|---|---|---|
| JEM1 clean / masked | **0** | 0 |
| JEM2 clean / masked | **0** | 0 |
| JEM2 spades *(positive control)* | 5 | 1 |

The positive control firing while the cleaned assemblies stay at zero is what makes this a
result rather than an absence of evidence.

## 3. Annotation

The two strains had only `funannotate predict` output — no functional annotation — so the BFD
predictor set was run over them and over both *Marchantia* proteomes, **matching BFD's own
invocations** (`Fungi_BFD/nextflow/modules/BFD/*`) so that a candidate here means what it means
across the 292 BFD Mucoromycota genomes.

| Step | Invocation | Matches BFD module |
|---|---|---|
| SignalP 6 | `signalp6 -org euk --mode fast -format none -bs 100`, GPU | `RUN_SIGNALP` |
| TMHMM 2.0 | `tmhmm --short --noplot` | `RUN_TMHMM` |
| PredGPI | `predgpi.py -m gff3` | `RUN_PREDGPI` |
| TargetP 2.0 | `targetp -format short -org non-pl` | `RUN_TARGETP` |
| Pfam | `hmmsearch --cut_ga --noali` | `RUN_PFAM` |
| dbCAN | `dbcanlight search -m cazyme / -m sub`, `conclude` | `RUN_CAZY` |

| Proteome | Proteins | Signal peptides |
|---|---|---|
| JEM1 (*M. hiemalis*) | 11,487 | 668 |
| JEM2 (*M. circinelloides*) | 11,485 | 676 |
| *M. polymorpha* Tak-1 | 20,412 | 1,911 |
| *M. paleacea* | 18,566 | — |

Three gaps closed in the process, none cosmetic:

* **PredGPI and TargetP were not being run at all**, so two of the four secretion categories
  (`gpi_anchored`, `mitochondrial`) were unreachable for these genomes.
* **CAZymes came from a plain `hmmsearch` against `dbCAN.hmm`**, which yields family names and
  no substrates. Tier T2 is *defined* on dbCAN substrate strings, so it would have been
  systematically empty.
* **SignalP's table is built from the GFF3**, as `merge_signalp.py` does, not from
  `prediction_results.txt`.

## 4. Candidate effectors

Secretion call and tiers are `bin/build_secretome.py`'s, reused unchanged so the JEM strains and
the BFD genomes are scored identically.

| Category | JEM1 + JEM2 |
|---|---|
| secreted | 863 |
| gpi_anchored | 171 |
| membrane | 266 |
| weak_signal | 41 |
| mitochondrial | 3 |
| **total candidates** | **1,344** |

| Tier | Definition | n |
|---|---|---|
| T1 chitin/chitosan | CE4, GH18/20/46/75, LysM/CBM | **99** |
| T2 plant-facing | dbCAN plant substrate, cutinase, lipase, peroxygenase | 55 |
| T3 effector-like | no Pfam, ≤ 300 aa, ≥ 4 Cys, ≥ 2 % Cys | 155 |
| T4 other secreted | — | 1,035 |

**The fold set** is deliberately *not* T3. A protein enters it if it is secreted or
GPI-anchored and either carries no Pfam domain and is ≤ 300 aa (`unknown_small`), or is tier T1
(`chitin_active`). T3's cysteine thresholds are a reasonable prior for classical effectors but
would discard exactly the cysteine-poor unknowns a structure is most informative about:
filtering on a sequence feature is self-defeating when the premise is that sequence tells us
nothing.

## 5. Lineage restriction

Candidates were clustered with the 195,284 secreted candidates of all 292 BFD Mucoromycota
genomes (MMseqs2, 50 % identity, 80 % coverage) and each family's reach measured.

| Reach | JEM1 | JEM2 | JEM2, conspecifics hidden |
|---|---|---|---|
| **genome_specific** | **105** | **4** | **13** |
| species_restricted | 3 | 5 | 0 |
| genus_restricted | 22 | 19 | 18 |
| mucorales_wide | 103 | 173 | 170 |
| mucoromycota_wide | 53 | 47 | 47 |

105 against 4 looks like a large biological difference, and the first thing to ask is whether
it is an artefact of which genomes happen to have been sequenced: **JEM2 has 24 conspecific
*M. circinelloides* genomes in BFD, JEM1 has none**, because BFD's single "*M. hiemalis*" is the
mislabelled one. Hiding JEM2's conspecifics (third column) raises its `genome_specific` count
only from 4 to 13 — so sampling explains part of the gap but not most of it.

**Caveat that remains.** JEM1 is the only *M. hiemalis* genome available, so "strain-specific"
and "species-specific" cannot be distinguished for it. The honest statement is that JEM1 has
~105 secreted proteins with no homologue in any sequenced Mucoromycota genome, and that this
is ~8× JEM2's matched figure.

## 6. Predicted folds

ESMFold (single-sequence, so no MSA and nothing leaves the cluster) on 640 candidates; 559
modelled, 81 skipped as > 400 aa. Foldseek against PDB100, CATH50 and the AlphaFold proteome
set, fold called at TM-score ≥ 0.5 with pLDDT ≥ 70.

| Fold set | assigned | low_confidence | no_match |
|---|---|---|---|
| **chitin_active** | **35 (95 %)** | 2 | 0 |
| **unknown_small** | 52 (15 %) | 266 | 28 |

| pLDDT | n |
|---|---|
| ≥ 90 | 11 |
| 70–90 | 121 |
| 50–70 | 156 |
| < 50 | 271 |

**The split is the main methodological result.** ESMFold folds the chitin-active enzymes almost
perfectly and fails on most orphan small proteins — which is what single-sequence prediction is
expected to do, and why the unfoldable set is being re-run with AlphaFold 3 and real MSAs.

### Chitin-axis hits

| Protein | TM | pLDDT | Fold |
|---|---|---|---|
| JEM2MUCOR_005224 | 0.96 | 96.2 | PdaC peptidoglycan deacetylase |
| JEM1MUCOR_009286 | 0.97 | 94.6 | Chitin deacetylase, *Colletotrichum lindemuthianum* |
| JEM2MUCOR_010829 | 0.97 | 93.0 | Chitin deacetylase AngCDA, *Aspergillus niger* |
| JEM1MUCOR_002372 | 0.96 | 93.0 | Chitin deacetylase AngCDA |
| JEM1MUCOR_000290 | 0.92 | 83.7 | Chitinase domain-containing protein 1 |
| JEM1MUCOR_006585 | 0.94 | 76.2 | **LysM domain** |
| JEM1MUCOR_004460 | 0.92 | 73.5 | **LysM domain** |
| JEM2MUCOR_009931 | 0.93 | 72.5 | **LysM domain** |
| JEM2MUCOR_004445 | 0.94 | 72.5 | **LysM domain** |
| JEM1MUCOR_005242 | 0.93 | 72.0 | **LysM domain** |
| **JEM1MUCOR_001522** | 0.84 | 71.1 | **cerato-platanin** |

`JEM1MUCOR_001522` is the clearest demonstration of why structure was worth predicting: it has
**no Pfam domain**, was classed `unknown_small`, and is a cerato-platanin — a well-described
fungal chitin-binding effector/elicitor family — by fold alone.

## 7. Host target panel

Built identically for both liverworts, because the isolation host is not recorded anywhere in
the MorenoLab project. The two differ in the respect that matters: *M. paleacea* forms
arbuscular mycorrhizal symbiosis and *M. polymorpha* has lost it.

| Class | *M. polymorpha* | *M. paleacea* |
|---|---|---|
| **LysM receptor** | **5** | **4** |
| chitinase | 14 | 18 |
| chitin_binding | 7 | 3 |
| LRR receptor | 200 | 161 |
| cysteine protease | 11 | 23 |
| protease inhibitor | 2 | 4 |
| subtilase | 10 | 16 |
| PR protein | 50 | 22 |
| lectin | 11 | 4 |
| expansin | 38 | 31 |

Receptor kinases enter as their **ectodomain**, not the full receptor: the interaction surface
is outside the membrane, and a cytoplasmic kinase domain would dominate any structure
prediction.

**Both species carry a comparable LysM receptor complement**, so the unrecorded host species
does not block the chitin screen. The other differences should not be read biologically yet —
*M. polymorpha* has 20,412 annotated proteins against *M. paleacea*'s 18,566, and the
annotations come from different groups, so counts need normalising first. The panel also makes
no orthology claim between the two species; that is a phylogenetic question, not a Pfam one.

## 8. Why the chitin axis

Mucorales cell walls are chitin **and chitosan**. In *M. paleacea*, the LysM pair
MpaLYR/MpaCERK1 discriminates symbiosis from immunity by chitin-oligomer length — CO7 triggers
immunity, CO4 activates symbiosis genes ([Cell 2025](https://ouci.dntb.gov.ua/works/4VP1pDNb)).
In *M. polymorpha*, chitin and peptidoglycan are the patterns that trigger a ROS burst while
other MAMPs do not ([Curr Biol 2023](https://www.pubmed.gov/37619565/)).

That makes a mechanistic screen rather than a fishing expedition:

| Fungal side | Predicted action | Host side |
|---|---|---|
| CE4 chitin deacetylase | chitin → chitosan, poorly bound by LysM | MpLYK1 / MpLYR ectodomains |
| LysM / CBM50 protein | sequesters CO7 before the receptor sees it | same, competitively |
| cerato-platanin | binds chitin oligomers | same |
| GH18 / GH20 | trims CO7 → CO4: immunity → symbiosis signal | same |

For LysM biology the decisive interaction is protein–**ligand**, not protein–protein, so the
screen models effector + (GlcNAc)ₙ and receptor ectodomain + (GlcNAc)ₙ at n = 4 and 7 and
compares binding, rather than relying on protein–protein iPTM alone.

## 9. In progress

* **AlphaFold 3, chitin panel** — the folds above against both species' LysM ectodomains and
  chitinases, plus CO4/CO7 ligands. MSAs computed once per chain and reused across pairs.
* **AlphaFold 3, unknowns** — the ~300 candidates ESMFold left below pLDDT 70. Most are
  `mucorales_wide`, so homologues exist and an MSA-based predictor has something to work with.

**Controls are built in.** iPTM has no universal threshold, and the AlphaFold-Multimer figures
(> 0.8 confident, < 0.6 likely wrong) were calibrated on large complexes, not small
cysteine-rich effectors, where false negatives are common. The run includes 4EBZ (AtCERK1 LysM
ectodomain with chito-oligosaccharide) and 6Q76 (AVR-PikD with the Pikp-1 HMA domain) as
positives and a cross-pairing as negative. `bin/af3_report.py` prints the controls first and
applies no threshold of its own: a cut is only worth using if the positives clear it and the
negatives do not.

## 10. Limitations

1. **No expression data.** Every candidate is a prediction from sequence and structure; none is
   known to be expressed in planta.
2. **Structure prediction is prioritisation, not evidence of interaction.** A confident fold
   says what a protein resembles. A confident iPTM says two chains can be modelled together.
   Neither is a binding assay.
3. **JEM1's lineage restriction is confounded with species sampling** (§5).
4. **The *Marchantia* panel counts are not normalised** (§7), and assert no orthology.
5. **The host species of the isolates is unrecorded**, so both liverworts are treated
   symmetrically; knowing it would let them be weighted.
6. **ESMFold's low-pLDDT models are ambiguous.** Low confidence may mean prediction failure or
   genuine disorder — small, cysteine-rich, secreted and unstructured-until-bound is a real
   effector phenotype. The pTM values near 0.1 suggest both are present. The AF3 re-run is what
   separates them.

## 11. Upstream issues found

| Where | Issue | Status |
|---|---|---|
| `Fungi_BFD/nextflow/bin/build_BFD_duckDB.sh:316` | `SELECT DISTINCT protein_id FROM predgpi` makes `is_gpi_anchored` TRUE for **every** protein — PredGPI writes a `Chain` row per protein and a `GPI-anchor` row only for anchored ones (224 M vs 3.3 M rows) | **Reported, not fixed** — one-clause fix; `BFD_dashboard`'s `J4_function_counts.sql` already filters correctly |
| `tables/gene_transcripts.parquet` | start/stop codon flags wrong (case-sensitive matching on soft-masked genomes); fixed upstream in `fb1d7cd`, but the parquet predates it | Worked around — codons read from sequence |
| `GCA_977110905.1_gzMucHiem1` | labelled *M. hiemalis* in BFD; its ITS says otherwise | Reported here |
| `MorenoLab/Mucor_genomes/asm/AAFTF/` | JEM2's AAFTF output is in a stray `x/` subdirectory, leaving four broken `ITS/JEM2.*` symlinks | Worked around; MorenoLab's files untouched |

## 12. Notes for anyone re-running this

* **ESMFold must run on an a100.** The `esm/2.0.0` module ships torch 1.13.1+cu117, built for
  sm_37–sm_86. On an ada6000 (sm_89), h100 (sm_90) or blackwell6000 (sm_120) it does **not**
  refuse — it warns and returns garbage. One run produced 487 nonsense errors and 72
  "successful" models that were equally worthless. `bin/esmfold_batch.py` now checks the device
  capability against `torch.cuda.get_arch_list()` and exits.
* **exfab cannot run ESMFold or AF3** — its node has only ada6000.
* The `signalp`, `esm`, `predgpi`, `targetp` and `dbcanlight` modules all activate conda
  environments via the shell function, so `eval "$(conda shell.bash hook)"` must precede them.
* **Reusable intermediates do not belong in `$SCRATCH`.** The DuckDB export, candidate FASTA and
  MMseqs families were written there once and lost when the job ended; they live in
  `analysis/mucor_secretome/cache/` (gitignored) now and rebuild in under three minutes.
* Slurm's `START_TIME` for pending jobs here is a pessimistic placeholder, not a forecast.

## 13. Reproducing

```bash
sbatch analysis/mucor_secretome/slurm/secretome.sh        # Mucoromycota secretome + cache
sbatch analysis/mucor_secretome/slurm/signalp_gpu.sh      # SignalP for the non-BFD panels
sbatch analysis/mucor_secretome/slurm/annotate_panels.sh  # the other predictors
sbatch analysis/mucor_secretome/slurm/its_identify.sh     # species identification
sbatch analysis/mucor_secretome/slurm/contam_check.sh     # contamination verification
sbatch analysis/mucor_secretome/slurm/fold_effectors.sh   # ESMFold + Foldseek + host panel
sbatch analysis/mucor_secretome/slurm/conservation.sh     # lineage restriction
sbatch analysis/mucor_secretome/slurm/af3_prepare.sh      # then af3_msa.sh, af3_fold.sh
```

Tests: `python3 -m unittest discover -s analysis/mucor_secretome/test` (68 tests).
