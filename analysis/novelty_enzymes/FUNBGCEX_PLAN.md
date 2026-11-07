# Plan sketch: a FunBGCeX-style survey of the BFD genomes

Draft 2026-11-07, for review with Fable. This is a sketch. Nothing here has run yet.
Each fact is marked **[measured]** (checked on HPCC today), **[paper]** (from Tang & Matsuda 2024 or the FunBGCeX README), or **[assumption]** (not checked).

## 1. Goal

Repeat the FunBGCeX approach (Tang & Matsuda, *Nat. Commun.* 2024, [PMC11109268](https://pmc.ncbi.nlm.nih.gov/articles/PMC11109268/)) on the BFD genome set. Then extend it in the one direction that BFD makes possible: genomes that had no NCBI annotation.

The claim we want to test, stated so that it can fail:

> BFD annotation of ~18k previously unannotated genomes finds biosynthetic gene clusters (BGCs) built around domainless enzymes that a 1,990-genome NCBI-annotated survey could not see. Divergent BFD proteins with no Pfam domain can be turned into new HMMs that find further BGCs.

Three questions follow:

- Q1. Does the published FunBGCeX result reproduce on BFD gene models? (method check)
- Q2. How many BGCs, and how many novel ones, come from genomes that had no prior annotation? (the delta-gain claim)
- Q3. Do new HMMs built from domainless BFD proteins find BGCs that the published HMM library misses? (the method extension)

## 2. What FunBGCeX does [paper]

| Step | Detail |
|---|---|
| Input | A directory of GenBank files. CDS features must carry `/translation`. |
| Target detection | `hmmscan` against a target HMM. E-value and conditional E-value cutoff 1e-5. Published targets: Pyr4.hmm and PydY.hmm (domainless terpene cyclase families). |
| Core-enzyme library | Pfam-A plus manual, SMART and TIGRFAM profiles. PKS/NRPS HMMs built from FunBGCs sequences. Plus 59 HMMs built from domainless FunBGCs proteins. |
| Neighborhood | Up to 20 flanking genes per target hit. |
| Clustering | DIAMOND hits to FunBGCs at ≥25 % identity and ≥50 % coverage. Two genes cluster if hits share a parent BGC and lie within 15,000 bp (distance ignored at ≥50 % identity to the same protein). Stop after three uncorrelated genes. Biosynthetic genes within 2,500 bp join. Proteins <200 aa dropped unless biosynthetic. Self-resistance: a ≥50 % identity homolog elsewhere in the genome. |
| Known-BGC database | FunBGCs: ~700 curated fungal BGCs, 5,070 proteins. 572 proteins had no Pfam domain. |
| New HMMs | Domainless proteins grouped by BLASTp (E 0.01), aligned with Clustal Omega 1.2.3, `hmmbuild`. 59 profiles. |
| Similarity score | Mean % identity × (homologous genes / genes in the known BGC). Homologous = ≥45 % identity, ≥50 % coverage. |
| Published scale | 1,990 NCBI-annotated genomes. Pyr4 mode: 1,050 BGCs; 690 with no second core gene. antiSMASH 6.1.1 found 325 Pyr4-containing BGCs in the same genomes. Taxonomy of the 1,050: 1,010 Ascomycota, 25 Basidiomycota, 12 Chytridiomycota, 3 Mucoromycota. |
| Software | `pip install funbgcex`, Python 3.9–3.13, MIT license. Needs `diamond` (README tested 2.1.13), `hmmscan` and `hmmfetch` (HMMER 3.4) on `$PATH`. Modes: default (all BGCs), `-m target -f FASTA`, `-m pfam -p NAME`, `-m sre`, `--no_core`. Code: [github.com/ydmatsd/funbgcex](https://github.com/ydmatsd/funbgcex). Data and custom HMMs: Zenodo 10.5281/zenodo.8126803. |

Not yet checked: whether the CLI accepts a custom HMM as target (the README lists a protein FASTA for `-m target`, not an HMM), and which thresholds `funbgcex --help` exposes. This decides whether Q3 can use the stock tool or needs a fork. **[assumption]**

## 3. What we have [measured]

| Resource | State |
|---|---|
| Genome list | `species.parquet` / `samples.parquet`: 23,680 genomes. `gene_info` covers 23,452. |
| GenBank input | `Fungi_BFD_runs/genome_annotation/<out>/predict_results/<out>.gbk.gz`: 23,450 files. Spot check of one genome: 12,794 `/translation` qualifiers. So FunBGCeX can read BFD output with no conversion. |
| Pfam | `tables/pfam.parquet`, 368M domain rows (hmmscan domtblout). Gives the "domainless" set by anti-join with `gene_proteins`. |
| Gene coordinates | `gene_info`, `gene_intergenic_distances`. Enough to recompute neighborhoods outside FunBGCeX if needed. |
| SignalP, TMHMM, BUSCO, asm_stats | In `tables/`. Usable as QC covariates. |
| Prior-annotation status | `prior_annotation_status.tsv` (27,401 NCBI assemblies): 5,511 `baseline_ok`, 18,145 with no prior annotation. Join to BFD ASMID not yet built. |
| NCBI annotation of baseline genomes | `1KFG/2026/NCBI_fungi/source/NCBI_ASM/<acc>/`: `genomic.fna.gz`, `genomic.gff.gz`, `protein.faa.gz`. No `.gbff`. So the NCBI side needs GFF+FASTA → GenBank conversion or a `.gbff` download. |
| antiSMASH | BFD has `modules/BFD/ANTISMASH` and `MERGE_ANTISMASH`, but no outputs under `results/function/` and no table. Treat as not run. |
| DeltaGain novelty | Linclust tiers and no-hit sets from the funnel (see `HANDOFF.md`). The full funnel has not run yet. |

Phylum counts in `samples.parquet` (not all annotated): Ascomycota 19,412; Basidiomycota 3,275; Mucoromycota 294; Kickxellomycota 190; Microsporidia 163; Mortierellomycota 100; Chytridiomycota 89; Glomeromycota 60; others <30 each.

## 4. Work plan

### Phase 0 — Set up and reproduce (small)

1. Install FunBGCeX in a pixi or uv environment (`[workspace]` table if pixi). Pin DIAMOND and HMMER to the README versions. Record `funbgcex --help` in this folder.
2. Download Zenodo 8126803 (FunBGCs, custom HMMs, Pyr4.hmm, PydY.hmm) into `analysis/novelty_enzymes/ref/`. Record checksums.
3. Reproduction check. Pick ~20 genomes from the paper's Supplementary Data 4 that are also in BFD and have one of the published onoceroid / Pyr4 BGCs. Run FunBGCeX on (a) the NCBI annotation and (b) the BFD annotation of the same assembly. Report per-BGC recovery and boundary agreement.
   - Exit criterion: BFD models recover the published BGCs at a rate we report as a number. No pass/fail threshold set yet — **open decision**.
4. Timing pilot: 50 genomes, mixed sizes, in one SLURM array on `-p epyc,intel`. Measure wall time and memory per genome. The paper says "within a single day" for 1,990 genomes but gives no hardware, so we cannot extrapolate from it.

### Phase 1 — Full run on BFD (Q1, Q2)

5. Batch genomes so each job runs ~1–1.5 h, sized from the Phase 0 timing (per the HPCC sizing rule). Copy each `.gbk.gz` to `$SCRATCH`, decompress there, run, copy outputs back to `/bigdata`, compress with zstd.
6. Run three modes over all 23,450 genomes:
   - `-m target` with Pyr4 and PydY (direct repeat of the paper, at ~12× the genome count).
   - Default mode (all BGCs). Gives the global BGC set for the novelty analysis.
   - `-m sre` (self-resistance). Cheap if the default run already holds the hits — check whether it can reuse them.
7. Parse `allBGCs.csv` per genome into one Parquet table: genome LOCUSTAG, BGC id, coordinates, core enzymes, target hit, nearest FunBGC and similarity score, gene list.
8. Join to `prior_annotation_status.tsv` and to taxonomy. Report every count split three ways:
   - prior NCBI annotation (`baseline_ok`) vs BFD-only;
   - phylum;
   - FunBGCs similarity bin (known / related / novel; bin edges to be chosen and stated before looking at the result).

### Phase 2 — New HMMs from domainless BFD proteins (Q3)

9. Define the domainless pool: BFD proteins with no Pfam hit at the BFD cutoff, ≥200 aa, inside a Phase 1 BGC or within 20 genes of a core enzyme. Count it before clustering. The genome-wide domainless set is far too large to use without the BGC restriction. **[assumption — measure it]**
10. Cluster the pool (MMseqs2, e.g. 30 % identity / 80 % coverage; parameters are an open choice). Keep clusters seen in ≥N genomes and ≥2 genera, so one bad gene model cannot make a family.
11. Align (MAFFT or Clustal Omega to match the paper) and `hmmbuild`. Remove profiles that hit an existing Pfam or FunBGCeX HMM.
12. Use each new HMM as a target. If the CLI cannot take an HMM as target, either pass representative sequences via `-m target -f` or fork the target step. Decide after step 1.
13. Compare: BGCs found with new HMMs that the stock library did not find. That number is the Q3 result.

### Phase 3 — Gene-model artifact controls

These decide whether a reviewer believes the BFD-only clusters.

14. On the ~5,511 baseline genomes, compare NCBI and BFD annotations: same BGCs? same core-gene boundaries? Report the disagreement rate. This is the measured error of our gene calls in BGC context.
15. For each novel BGC: BUSCO of the genome, contig length vs BGC span (BGC at a contig edge = incomplete), fraction of genes with partial flags (`gene_transcripts`), intron-count outliers.
16. Conservation check: novel BGCs found in ≥2 independent genomes (different species) are stronger than singletons. Report both groups.
17. Structure check on a shortlist: predicted fold of the target and core enzymes (ESMFold or AlphaFold), TM-score to the family.
18. RNA-seq support where SRA data exist for the same species. Scope this only after the shortlist exists.

### Phase 4 — Prioritization and deliverables

19. Shortlist by the known/unknown framing: novel cluster + domainless core + conserved in ≥2 genomes + passes Phase 3 + host close to a tractable expression host.
20. Deliverables:
    - Parquet table of all BGCs, per genome, with provenance (FunBGCeX version, HMM set, DB checksums).
    - Counts table: BGCs per genome and novel BGCs, split by prior annotation and phylum.
    - New HMM set with build alignments, deposited (Zenodo).
    - Reproduction report (Phase 0) and artifact-rate report (Phase 3).

## 5. Comparators

- **antiSMASH.** The paper's main contrast is FunBGCeX vs antiSMASH (325 vs 1,050 Pyr4 BGCs). BFD has an antiSMASH module but no run. Options: run antiSMASH on a stratified subset only (cost unknown), or skip and compare against the published numbers. **Open decision.**
- **Published FunBGCeX run.** Supplementary Data 4 lists the 1,990 genomes. Overlap with BFD gives a direct like-for-like comparison.

## 6. Risks

| Risk | Effect | Mitigation |
|---|---|---|
| BFD gene models split or merge biosynthetic genes (large PKS/NRPS are common victims) | False novel BGCs, wrong boundaries | Phase 3 step 14 measures it on baseline genomes |
| Fragmented assemblies | BGCs cut at contig edges | Flag edge BGCs; report by N50 bin |
| Phylum imbalance (82 % Ascomycota) | Pooled counts hide early-diverging results | Always report per phylum |
| Domainless pool dominated by mis-predictions | New HMMs model noise | Require multi-genus support; structure check |
| FunBGCs database is ~700 BGCs, Ascomycota-heavy | "Novel" partly means "far from a small reference" | Also compare to MIBiG fungal entries; state the reference size with each novelty claim |
| Tool version drift (paper DIAMOND 2.1.7; README tests 2.1.13; DeltaGain uses 2.2.6) | DB incompatibility | Pin the README versions in a separate env |

## 7. Open questions for review

1. Scope of Phase 1: all three modes on all genomes, or Pyr4/PydY first and the global mode after?
2. Reproduction threshold in Phase 0: what recovery rate counts as "BFD models are fit for this"?
3. antiSMASH: run on a subset, or rely on the published comparison?
4. Domainless clustering parameters and minimum support (genomes, genera) for a new HMM.
5. Should Q2 use the DeltaGain funnel novelty (protein-level) as well as the FunBGCs similarity score (cluster-level)? They measure different things; both may be worth reporting side by side.
6. Is there an expression partner for validation? The paper validated in *Aspergillus oryzae* NSAR1. Without one, the output is a computational catalogue plus HMMs.

## 8. Not verified yet

- FunBGCeX runtime and memory per genome on HPCC.
- Whether `-m target` accepts an HMM.
- That all 23,450 `.gbk.gz` files carry translations (one genome checked).
- The ASMID ↔ NCBI accession join and how many BFD genomes it covers.
- Size of the domainless protein pool.
- The paper's Supplementary Data 4 overlap with BFD.
