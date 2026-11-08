# Project sketch: UPOs from early-diverging fungi (Mucoromycota, Glomeromycota, Chytridiomycota)

Status: plan sketch, 2026-11-07. Nothing launched. Census numbers are measured (CENSUS.md). Reference checks are in
REFERENCES.md. Every compute figure in section 7 is an **estimate**.

## Decisions (user, 2026-11-07)

These replace the matching open questions in section 9.

1. **Classification:** use the Kinner et al. 2021 Family I (short) / Family II (long) scheme
   ([PMC8293615](https://pmc.ncbi.nlm.nih.gov/articles/PMC8293615/)) for now. The EHD/EGD label conflict stays
   open and must be settled against structures.
2. **Scope of classes:** inclusive first. Report every class (UPO Family I, Family II, CPO-like, Pog, and
   "PF01328 but fails motif checks"). Decide later which classes are deliverables.
3. **Target lineages:** keep the early-diverging core, and add the newly sampled Dikarya lineages, for example
   lichens and Antarctic rock fungi. See "Newly sampled Dikarya lineages" below.
4. **Reference method:** UPObase is offline. Build our own standalone UPO identification tool that replicates
   the UPObase paper's method. Work is in `upo_id/` (in progress).
5. **Expression and assay:** expected to be arranged within about one year, once candidates exist. So the plan
   must end in a ranked candidate list fit for heterologous expression, not only a catalogue.

### Decisions on the upo_id tool (user, 2026-11-07, second round)

6. **UPO Family I vs CPO-like:** treat both as Family I for now. The split rests on 8 characterized enzymes
   against 16 uncharacterized CPOs, which is too few.
7. **Classification method:** replace HMM score competition with phylogenetic placement on a reference tree.
8. **UPObase data:** do not ask the authors. Build our own best reference set and method.
9. **CglUPO and CciUPO accessions:** still unverified. Look up later (future task).
10. **Residue probing:** use Allozymes' case study
    (https://www.allozymes.com/case-studies/engineering-specificity-into-unspecific-peroxygenases-upos),
    PMC12117455 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12117455/) and their references to define the
    key residue positions to score (see `RESIDUES.md`, in progress).

First result of the tool (2026-11-07, `upo_id/results/early807_*`): of 807 early-diverging PF01328 proteins,
9 are called UPO (3 Family I, 6 low-confidence Family II). All 9 are in NCBI-annotated genomes. None is in a
BFD-only genome. The newly sampled Dikarya lineages are not tested yet.

### Newly sampled Dikarya lineages (measured 2026-11-07, `sql/novel_lineages.sql`)

Counts use genomes with BFD gene models (23,451). Status from `prior_annotation_status.tsv`: 17,864
no prior annotation, 5,553 NCBI-annotated, 24 no match, 10 annotated without a protein file.

- Lineages with **no NCBI-annotated genome at all**:
  - Classes: 6 (4 Ascomycota, 2 Basidiomycota), 14 genomes in total. Examples: Sporopachydermiomycetes (4),
    Peribolosporomycetes (4), Alloascoideomycetes (3), Laboulbeniomycetes (1).
  - Orders: Ascomycota 19 orders / 54 genomes; Basidiomycota 3 / 9.
  - Families: Ascomycota 71 families / 170 genomes; Basidiomycota 36 / 111.
- Lichen classes (NCBI-annotated / BFD only): Lecanoromycetes 137 / 340; Arthoniomycetes 5 / 8;
  Lichinomycetes 4 / 6; Candelariomycetes 4 / 0; Coniocybomycetes 2 / 0.
- Antarctic rock fungi: *Antarctomyces pellizariae* is BFD only. *Cryomyces* has 4 NCBI-annotated genomes and
  1 BFD only. *Friedmanniomyces* genomes are all NCBI-annotated.

So, measured by genome count, Lecanoromycetes is the largest new Dikarya sample: 340 genomes have only BFD
annotation. No PF01328 counts exist yet for these lineages. Not checked: whether lichen genome assemblies
include photobiont or bacterial contigs, which would affect gene calls.

## 0. Question and claim

Question: which heme-thiolate peroxygenase (UPO) genes exist in early-diverging fungal genomes, and how many come
from genomes that only BFD has annotated?

Deliverable: a curated, phylum-stratified set of UPO candidates from Mucoromycota, Glomeromycota and
Chytridiomycota (plus related early-diverging phyla), with catalytic-motif status, gene-model QC, a phylogeny
against Dikarya UPOs, a candidate shortlist, and reusable HMMs.

Novelty claim, stated narrowly: candidates from genomes with no NCBI annotation (BFD_only) are new protein
sequences. Candidates from NCBI-annotated genomes are "new as UPO calls" at most, and only if the NCBI protein
was not already classified. The plan reports these two groups separately.

## 1. What the data already show (measured, CENSUS.md)

- 827 early-diverging genomes (Microsporidia excluded); 813 have proteins; 10.17 M proteins.
- BFD Pfam (`hmmsearch --cut_ga`) already finds **807 PF01328 proteins** in early-diverging genomes;
  **309** of them are in BFD_only genomes.
- Core three phyla (genomes NCBI / BFD_only; PF01328 proteins NCBI / BFD_only):
  - Mucoromycota: 145 / 148 genomes; 193 / 177 proteins; median 1 per genome; 49 BFD_only species absent from the
    NCBI-annotated set.
  - Glomeromycota: 43 / 17 genomes; 92 / 32 proteins; 6 new species.
  - Chytridiomycota: 81 / 8 genomes; 184 / 4 proteins; 3 new species.
- **Chytridiomycota is mostly NCBI-annotated already** (81 of 89). The BFD_only novelty in this phylum is small
  (4 proteins). For Chytridiomycota the value is a curated UPO call set, not new sequences.
- Zero PF01328 proteins at the gathering threshold in Mortierellomycota (100 genomes), Kickxellomycota (190),
  Blastocladiomycota (14), Neocallimastigomycota (26), Entomophthoromycota (12), Microsporidia (163).
- Hit shape differs by phylum: Mucoromycota hits are 250–350 aa with almost no signal peptide (1–3 %);
  Glomeromycota hits are mostly <250 aa with no SP; Chytridiomycota, Basidiobolomycota and Zoopagomycota hits have
  SP in ~50–56 %.
- A crude substring scan finds `PCP` + `EHD` in most hits and `EGD` in almost none (10 of 807). Interpretation is
  blocked by the label conflict in REFERENCES.md section 3.

What this means: the search itself is not the hard part. The hard part is the **UPO-vs-CPO-vs-Pog call**, because
PF01328 covers all of them (97,520 proteins BFD-wide, against 123 sequences in UPObase's UPO subfamilies I–V).

## 2. Target set

| Tier | Phyla | Reason |
|---|---|---|
| Core (user request) | Mucoromycota, Glomeromycota, Chytridiomycota | Named in the request. Have PF01328 hits. |
| Extended | Zoopagomycota, Basidiobolomycota, Monoblepharomycota, Aphelidiomycota, Calcarisporiellomycota | Have PF01328 hits (30, 86, 4, 3, 2 proteins). Basidiobolomycota is mostly BFD_only (16 of 18 genomes, 77 of 86 proteins). |
| Negative / sensitivity check | Mortierellomycota, Kickxellomycota, Blastocladiomycota, Neocallimastigomycota, Entomophthoromycota | 0 PF01328 hits at GA. Run the sensitive second pass on them. A true absence is a result. A hit found only by the second pass is a divergent-member candidate. |
| Excluded | Microsporidia | 0 hits; reduced obligate-intracellular genomes; not expected to carry UPOs (assumption, no data beyond the 0 hits). |
| Positive control | Ascomycota, Basidiomycota | Characterized UPOs come from here. See section 3.4. |

Taxonomy: BFD `samples.PHYLUM` already splits Mortierellomycota, Glomeromycota, Kickxellomycota etc. from
Mucoromycota (NCBI current usage). The plan uses `samples.PHYLUM` as is. Decision for the user: report under the
older "zygomycete" umbrella as well, or not.

Data paths (BFD layout, `Fungi_BFD_wt_csvlist/docs/BFD_Introduction.md`):
- Proteins: `Fungi_BFD_runs/seqdb/proteins/<sha1(LOCUSTAG)[0:2]>/<LOCUSTAG>.proteins.fa.gz` (bgzip + faidx);
  CDS in `seqdb/cds/`. `seqdb/genomes.tsv` maps LOCUSTAG to folder. Protein IDs are `<LOCUSTAG>_NNNNNN-T1`.
- Tables: `Fungi_BFD_runs/tables/` — `samples`, `pfam`, `signalp.signal_peptide` (positives only), `tmhmm`,
  `gene_transcripts` (`is_partial`, `has_start_codon`, `has_stop_codon`), `gene_introns`, `gene_exons`,
  `busco_genome`, `asm_stats`, `swissprot`.
- Prior annotation: `DeltaGain_Fungi/prior_annotation_status.tsv`, join `samples.ASMID = folder` (24 of 23,680
  unmatched; none early-diverging).
- Note: the BFD intro doc calls the Pfam table "hmmscan domtblout"; the PFAM module runs `hmmsearch --cut_ga`.
  Same output format; the GA threshold is what matters here.

## 3. Reference set

### 3.1 UPObase
- Offline from HPCC on 2026-11-07 (DNS failure). Not obtainable now.
- Fallbacks, in order: (a) email the UPObase authors (Faiza, Wang; Database 2019) for the FASTA;
  (b) Wayback Machine copy of the per-entry FASTA pages (not checked, API rate-limited);
  (c) rebuild an equivalent set: phmmer/jackhmmer from the characterized enzymes against UniProt fungi, then apply
  the UPObase motif strings (REFERENCES.md section 2). Option (c) is reproducible and does not depend on a dead site.

### 3.2 Profiles
- PF01328 (Pfam, GA cutoff) — broad family net. Already run in BFD.
- PROSITE PS51405 (heme haloperoxidase profile) — the profile Kinner Fig. 1 used. Run with `pfscan`/`ps_scan`
  as a second family net. Not yet run.
- Custom HMMs (built in this project): long-type UPO, short-type UPO, CPO, and one per early-diverging clade that
  forms in the tree (section 6). Build from MAFFT L-INS-i alignments with `hmmbuild`.

### 3.3 Characterized enzymes
13 with verified UniProt or PDB records (REFERENCES.md section 6): AaeUPO B9W4V6, CraUPO B9W4V8, LfuCPO P04963,
MroUPO (PDB 5FUJ/7ZBP; UniProt A0A023H437 is a 67-aa fragment), MweUPO A0ABF7PQ10, TteUPO G2R8M6, DcaUPO
A0AAN0NP18, HspUPO A0A1Y2TH07, PabUPO-II A0A4V1Q376, CviUPO (PDB 7ZCL), UPO1 P0DYC0, DchUPO (PDB 9TM7), DspUPO
(PDB 9WA6). CglUPO, CciUPO, MthUPO, HinUPO accessions are **not verified**; get them from the expression papers.

### 3.4 Positive control
- Run the full search and calling pipeline on (a) the 13 references, which must all pass; (b) a fixed Dikarya
  subset, for example one genome per order with BUSCO C >= 90 % (size not decided; estimate 300–600 genomes).
- Compare per-genome counts with UPObase-era numbers only for the same species. **Do not use "113 from >800
  genomes" as the target.** In the UPObase paper, 113 is the seed count (REFERENCES.md section 3.1).
- For NCBI-annotated early-diverging genomes, also search the NCBI protein set (DeltaGain baseline) and compare
  with the BFD calls. This measures what BFD re-annotation adds or loses for this family.

## 4. Search

1. **Pass 1 — family net (done).** Take PF01328 GA hits from `pfam.parquet` (807 early-diverging proteins).
   Add PS51405 hits. Union.
2. **Pass 2 — sensitive.** Targets: all 813 early-diverging proteomes, including the zero-hit phyla.
   - `jackhmmer` (3 iterations, inclusion E 1e-3) seeded with one long-type, one short-type and one CPO reference,
     against the early-diverging proteome set only.
   - `hmmsearch` with the custom HMMs at E <= 1e-3, no GA.
   - Optional: Foldseek against AlphaFold/ESMFold models of pass-2-only hits, with AaeUPO 2YOR, MroUPO 5FUJ and
     LfuCPO 1CPO as structural queries. Use only for the few hundred hits pass 2 adds, not the whole proteome.
3. **Genome-level rescue.** `miniprot` or `tblastn` of reference proteins against the genome assemblies of the
   target phyla. This finds UPO loci with no gene model or with a broken one. Report rescued loci separately.

Decision for the user: HHblits/HMM-HMM search adds little here if jackhmmer converges; skip unless pass 2 finds
very few divergent hits. (Opinion, no data.)

## 5. Calling and scoring each hit

Order matters: QC first, then motif call, then ranking.

| Step | Rule | Data source |
|---|---|---|
| Gene model QC | Flag `is_partial`, missing start/stop, length < 180 aa, PF01328 HMM coverage < 80 %, internal stop, fused models (> 600 aa or second unrelated Pfam domain) | `gene_transcripts`, `gene_proteins.length`, `pfam` |
| Genome QC | Flag genomes with BUSCO C < 80 % and contig N50 < 20 kb | `busco_genome`, `asm_stats` |
| Heme ligand | Cys in the PCP motif present at the aligned position of AaeUPO/MroUPO/LfuCPO proximal Cys | MAFFT alignment to references (not substring search) |
| Acid-base | Glu at the aligned position of the catalytic Glu (AaeUPO numbering; confirm the residue from the 2YOR / 5FUJ structure papers before coding) | alignment |
| Type | EGD-type vs EHD-type motif region, with the label conflict resolved on structures first (REFERENCES.md §3.3) | alignment + references |
| Short / long | Length, C-terminal Cys pair (long type), and tree placement | alignment, tree |
| Class call | UPO-long, UPO-short, CPO-like, Pog/unclassified. A hit without the proximal Cys is "not UPO" | combined |
| Secretion | SignalP 6 positive, TMHMM helices, GPI | `signalp`, `tmhmm`, `predgpi` |
| Distance | Max % identity to any of the 13 references (DIAMOND or global alignment); nearest reference name | computed |
| Novelty | BFD_only vs NCBI genome; for NCBI genomes, whether the NCBI protein exists and how it was annotated | `prior_annotation_status.tsv`, NCBI protein set |

All counts reported per phylum and per annotation group. Never pooled across phyla.

## 6. Gene-model artifact controls

- **Glomeromycota.** Genomes are large (median 109–140 Mb) and fragmented (median N50 16–34 kb). Hits are short
  (median 225–233 aa). Check whether short length is biology or truncation: compare to the HMM coverage, start/stop
  flags, and to the same locus in the other strain of the same species where one exists (25 + 9 species, 60 genomes,
  so many species have more than one genome).
- **Mucoromycota.** Many species have whole-genome duplications (assumption from literature, not checked here).
  Report paralog pairs; do not count them as independent discoveries.
- **Contamination.** A UPO in a Glomeromycota or Chytridiomycota assembly can come from a co-cultured or
  contaminating ascomycete. For each candidate: neighbor genes on the same contig must have best hits in the same
  phylum; flag contigs carrying only the candidate. Check GC and codon usage against the genome (`codon_freq`).
- **RNA-seq support.** `Fungi_BFD_runs/species_with_real_rnaseq.txt` lists 52 species BFD-wide. Not checked how many
  are early-diverging. Use RNA-seq evidence where it exists; mark the rest "no expression evidence".
- **Intron structure.** Compare intron positions of candidates with references and with each other (`gene_introns`).
  Shared intron positions support real gene models.
- **Per-phylum artifact rate.** Report the fraction of hits flagged at each QC step, per phylum.

## 7. Compute sizing (all estimates)

| Step | Size | Estimate |
|---|---|---|
| Pass 1 extraction (done) | 20 GB parquet scan | measured: ~1 min per duckdb scan; 3 min for two extractions |
| Sequence extraction from seqdb | 807 proteins, 398 genomes | measured: < 1 min with `samtools faidx` |
| jackhmmer, 3 seeds x 3 iter vs 10.17 M proteins | ~4 GB FASTA | estimate: 1–3 h on 16 CPUs; one SLURM job |
| hmmsearch, ~6 custom HMMs vs 10.17 M proteins | same | estimate: < 1 h on 16 CPUs |
| PS51405 scan | 10.17 M proteins | estimate: 1–3 h; can be limited to pass-1/pass-2 hits instead |
| miniprot rescue, 13 refs vs ~827 genomes (~40 GB of sequence, estimate) | – | estimate: 2–4 h as ~10 batched jobs of ~80 genomes |
| Structure prediction | shortlist only, 50–200 proteins | estimate: hours on one GPU, or use AFDB/ESM Atlas where present |
| Dikarya positive control | 300–600 genomes | estimate: < 2 h hmmsearch; reuse `pfam.parquet` for PF01328 |
| Trees | 1–3k sequences, IQ-TREE | estimate: hours, one job |

HPCC rules: `-p epyc,intel --constraint="ryzen|broadwell|cascade|sapphire"` (no abu_dhabi; no AVX2 there);
work in `$SCRATCH`, copy results to `/bigdata`; compress tabular output with `zstd`; no `BASH_SOURCE` in SLURM
scripts. Batch to ~1–1.5 h per job; this whole project fits in fewer than 20 jobs (estimate).

## 8. Deliverables

1. Census table (done, CENSUS.md) — genomes, species and PF01328 hits per phylum, split NCBI vs BFD_only.
2. Candidate table (one row per protein): phylum, genome, annotation group, length, motif status, class call,
   SP/TM, QC flags, nearest reference and identity, tier.
3. Phylogeny: early-diverging hits + references + a Dikarya sample + CPO/Pog outgroup sequences. Show where the
   early-diverging clades sit relative to long-type and short-type UPOs.
4. Tiers:
   - Tier 1: QC-clean, Cys + Glu present, UPO class call, BFD_only genome, has SP (expression-friendly).
   - Tier 2: as Tier 1 but NCBI genome, or no SP.
   - Tier 3: motif-complete but gene model flagged, or Pog/unclassified.
   - Rejected: no proximal Cys, or artifact.
5. HMMs: custom profiles per class and per early-diverging clade, with the seed alignments, deposited
   (Zenodo, estimate of what is useful; decision for the user).
6. A short note on absence: phyla with zero hits after the sensitive pass.

## 9. Open questions for the user

1. **Classification scheme.** Use Kinner Family I/II (long/short) or UPObase UPO I–V vs Pog? They give very
   different numbers. The EHD/EGD label conflict needs a decision based on structures.
2. **Is a CPO-like or Pog hit a deliverable?** Most early-diverging hits carry EHD. If only "true UPO" counts, the
   final set could be small. Decide before running.
3. **Chytridiomycota scope.** 81 of 89 genomes already have NCBI annotation, and BFD_only adds 4 proteins. Is a
   curated call set for this phylum enough, or should the novelty claim rest on Mucoromycota, Glomeromycota and
   Basidiobolomycota?
4. **UPObase replacement.** Ask the authors for the data, or rebuild from UniProt with their motif rules?
5. **Expression capacity.** Is there an expression/assay partner? This decides whether Tier 1 needs SP and
   long-type architecture, or only sequence quality.
6. Include the extended phyla (Basidiobolomycota has 77 BFD_only hits) in the headline set?
7. Size of the Dikarya positive-control subset.

## 10. Claims not verified

- "113 putative UPOs from >800 genomes in 35 species" (Kinner 2021 citing Faiza 2019a). The UPObase paper uses 113
  as its seed count. Faiza 2019a not read.
- Zámocký et al. 2025 statement that Mucoromycota, Glomeromycota and Chytridiomycota are under-investigated
  (repository unreachable). Kinner 2021 does not say this.
- UniProt/GenBank accessions for CglUPO, CciUPO, MthUPO, MfeUPO, MhiUPO, HinUPO, CveUPO.
- The exact residue number of the acid-base Glu in AaeUPO. InterPro confirms a Glu acid-base and Cys ligand in the
  family, but I did not read the structure paper.
- Whether UPObase still exists at another address.
- That `signalp.signal_peptide.parquet` holds positives only is from the BFD intro doc; I did not test it.
- BACKGROUND.md statements on expression yields (mg/L, ~10 % purification recovery) — not checked; not used here.
- Whole-genome duplication in Mucoromycota as a paralog source — assumption from literature, not measured.
- How many of the 52 RNA-seq species are early-diverging — not measured.

## Files

- `CENSUS.md` — measured counts.
- `REFERENCES.md` — source checks, reference accessions, conflicts.
- `sql/` — duckdb queries used (`${SCRATCHPAD}` was a session scratch dir; replace with a local path).
- `data/early_pf01328_hits.parquet` — 807 early-diverging PF01328 proteins with phylum, annotation group, length,
  SP flag, HMM coordinates.
- `data/early_pf01328_hits.fa.zst` — their sequences (from seqdb).
- `data/motifs.tsv` — crude PCP/EHD/EGD substring flags (first pass only).
