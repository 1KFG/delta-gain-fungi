# antiSMASH lite for BFD: what we measured and what we understand

Status: working note, 2026-11-22. Companion to `README.md` (the run recipe) and to `docs/BFD_ANTISMASH_DESIGN.md` in `Fungi_BFD` (the BFD side's design note, written by another session). Nothing here is a biological finding.

Labels: **MEASURED** = I ran or counted it, the method is named. **EXTRAPOLATION** = scaled from the genomes tested. **PROPOSED** = my view. **UNVERIFIED** = not checked.

## 1. Scope and limits

- 35 BFD genomes were run (Eurotiomycetes and Sordariomycetes: *Aspergillus*, *Penicillium*, *Nigrospora*, *Ophiocordyceps*, *Fusarium*). Yeasts, basidiomycetes and early-diverging fungi were **not** tested.
- Tool: antiSMASH 8.0.4 in the BFD apptainer image, `--taxon fungi --genefinding-gff3`, 8 CPUs, `short` partition.
- The test set was chosen by BUSCO and contig count and by the owner's species list. It is not a random sample of BFD.
- A `known` cluster call means similarity to a MIBiG entry. It does not show that the compound is made.

## 2. What the lite run keeps and drops

| Item | Old BFD module | Lite (default now) |
|---|---|---|
| Cluster detection (antiSMASH rules) | on | on |
| `--cb-knownclusters` (MIBiG) | off | on |
| `--cb-general` | on | on in most runs here; off in the tests of section 5 and proposed off for BFD |
| `--fullhmmer`, `--clusterhmmer` (Pfam annotation) | on | off |
| `--pfam2go` (GO terms) | on | off |
| Raw antiSMASH directory (JSON, GenBank, HTML, zip) | kept as JSON | not kept |
| Stored output | JSON gz | five small files per genome |

The five stored files per genome (21 to 46 KB gz, three genomes measured): `regions`, `cds` (genes inside regions), `known` (top-3 MIBiG hits per region), `known_pairs` (our gene to MIBiG gene) and `run.json` (antiSMASH and MIBiG versions, mode).

What lite loses: Pfam domain hits per gene. What it keeps: smCOG annotation, in the `cds.functions` text.

## 3. Cost (MEASURED, 35 genomes)

| Quantity | Value |
|---|---|
| Wall time per genome, lite with `--cb-general` | min 500 s, median 965 s, mean 1,343 s, max 4,408 s (8 CPUs) |
| Wall time, full flag set, Af293 | 6,535 s (same 44 regions as lite, 786 s) |
| Peak memory (sacct MaxRSS of the batch step) | 8.7 to 12.4 GB, median 11.2 GB |
| Raw antiSMASH directory per genome | 177 to 479 MB, median 264 MB |
| Raw directory with `--no-zip-output --no-summary-gbk --no-region-gbks` (A1163) | 102 MB versus 206 MB |
| Time with those three flags (A1163) | 951 s versus 881 s: no saving, within the node-to-node variation |
| Stored size, five files | 21 to 46 KB per genome |

EXTRAPOLATION to 24,000 genomes at the mean: 8,950 8-CPU hours (71,600 CPU hours). Raw output would be about 6 TB if kept; it never leaves `$SCRATCH`.

Run times vary by node: the same genome ran 3,450 and 3,896 s (*A. niger* CICC 2151) and 4,247 and 4,408 s (*N. sphaerica* CBS 166.26). A single timing has a spread of roughly 10 to 15%.

## 4. Where the time goes (MEASURED, per-module timings in the antiSMASH JSON)

*N. sphaerica* CBS 166.26 (1,129 contigs, the slowest genome), control run, wall 4,308 s:

| Module | Seconds | Share of module time |
|---|---|---|
| `clusterblast` | 3,660 | 94.8% |
| `genefunctions` | 64 | 1.7% |
| `hmm_detection` | 57 | 1.5% |
| all others | about 78 | about 2% |

Module time sums to 3,859 s. The other 449 s is start-up and output writing.

The JSON does not separate the general ClusterBlast (a large database of antiSMASH-predicted clusters) from the known-cluster ClusterBlast (MIBiG). I **suspected** that the general one dominates and tested it (section 5).

## 5. Dropping `--cb-general` (MEASURED on 3 genomes)

Same flags otherwise (`--cb-knownclusters` only, no `--cb-general`), compared with the earlier lite run of the same genome:

| Genome | Wall time with `--cb-general` | Without | Speed-up | Regions | `known` top-3 rows | `known_pairs` top 3 |
|---|---|---|---|---|---|---|
| *A. fumigatus* A1163 | 881 s (951 s on a repeat) | 233 s | 3.8x to 4.1x | 43 = 43, identical | 51 = 51, identical | identical |
| *A. niger* CICC 2151 | 3,896 s (3,450 s on a repeat) | 562 s | 6.9x (6.1x vs the repeat) | 97 = 97, identical | 51 = 51, identical | identical |
| *N. sphaerica* CBS 166.26 | 4,308 s (4,408 s on a repeat) | 807 s | 5.3x | 72 = 72, identical | 37 = 37, identical | identical |

"Identical" means the same coordinates, classes, contig-edge flags and gene counts for regions, and the same accession, similarity and gene hits for `known`. Module time in ClusterBlast without `--cb-general`: A1163 25 s of 163 s; CICC 2151 128 s of 420 s; CBS 166.26 156 s of 345 s.

So on these three genomes the general ClusterBlast accounts for most of the wall time and changes none of the stored tables. The `blast_top` table (general ClusterBlast hits) is the only output that needs it, and the BFD minimal set does not store it.

Limits: three genomes, all Pezizomycotina. Not tested: yeasts, basidiomycetes, Mucoromycota. The speed-up grows with the time of the control run, so the fragmented genomes gain most. The mean time without `--cb-general` across genomes is not measured.

## 6. Dropping short contigs does not help (MEASURED)

*N. sphaerica* CBS 166.26, wall time (the same genome, four runs):

| Run | Wall time | Regions |
|---|---|---|
| control | 4,308 s | 72 |
| contigs under 5 kb removed from FASTA and GFF3 (957 of 1,129 contigs kept) | 4,742 s | 72, identical |
| contigs under 10 kb removed (783 kept) | 4,061 s | 71 |
| antiSMASH `--minlength 10000` | 3,708 s | 71 |

- The differences are within the node-to-node spread and not in one direction.
- The 10 kb runs lost one region: a 5.9 kb terpene region on a contig edge.
- `--minlength` exists in 8.0.4 and defaults to 1000 bp (MEASURED, `--help-showall`).
- Across 33 genomes with 2,140 regions: no region sits on a contig under 5 kb; 3 sit under 10 kb and 6 under 20 kb. All 6 are on a contig edge and have a biosynthetic gene; 1 has a MIBiG hit of 75 or higher.
- Wall time tracks contig count (Spearman 0.80) and not genome size (0.11), but removing the short contigs did not remove the time: ClusterBlast scales with the number of regions and genes, not with contig count alone.

PROPOSED: do not filter contigs. `filter_contigs.py` stays in the folder for reuse.

## 7. Keeping more MIBiG hits per region (MEASURED, 12 genomes, 310 regions with a hit)

antiSMASH ranks hits by core-gene hits and score, not by the similarity number, so the best-similarity hit is not always rank 1.

| Question | Result |
|---|---|
| Best-similarity hit outside ranks 1 to 3 | 13 of 310 (4.2%) |
| Best-similarity hit outside ranks 1 to 5 | 5 of 310 (1.6%) |
| A rank 4 or 5 hit at 75 or higher that beats ranks 1 to 3 | 4 of 310 (1.3%); 3 of the 4 are on one-protein MIBiG entries |
| Regions that change from `weak` to `known` with top 5 | 2 of 310 (0.65%) |
| `known` gz TSV per genome, top 3 / top 5 | 1,862 / 2,247 bytes |
| `known_pairs` gz TSV per genome, top 3 / top 5 | 7,445 / 8,809 bytes |
| Parquet, 12 genomes, `known` top 3 / top 5 | 12 KB / 14 KB |
| Parquet, 12 genomes, `known_pairs` top 3 / top 5 | 71 KB / 82 KB |

`known_pairs` grows by 16%, as `known` does, so the mixed option (`known` top 5, `known_pairs` top 3) saves little. EXTRAPOLATION to 24,000 genomes (2,000 times the 12): about 165 MB at top 3 and 191 MB at top 5 for these two tables. PROPOSED: the gain is small; the decision on 3 or 5 is the owner's.

## 8. How to read the MIBiG columns (MEASURED, 35 genomes)

- Use the highest similarity among the kept hits, not rank 1 (the first table I sent used rank 1 and was wrong for 36 regions of 1,351).
- A similarity of 100 can rest on one gene: 188 of 430 `known` regions (at 32 genomes) have a best hit on a one-protein MIBiG entry. Filter on `n_entry_proteins` and look at `entry_fraction_hit`.
- All 430 best hits were marked `questionable` in the MIBiG quality field (the MIBiG value, not mine).
- Locus tags in `cds` match the BFD `gene_info.gene_id` for all 40,435 rows tested (0 unmatched).

## 9. Storage in BFD (MEASURED on the 32-genome test; EXTRAPOLATION to 24,000)

Three zstd parquet tables keyed by `LOCUSTAG` (same key as `cazy` and `pfam_hmmscan`), SQL in `build_bgc_parquet.sql`: `bgc_cluster` (one row per region), `bgc_cluster_gene` (one row per gene in a region), `bgc_mibig_hit` (one row per kept hit). No sequences.

| Table | 32 genomes | Rows | EXTRAPOLATION, 24,000 genomes |
|---|---|---|---|
| `bgc_cluster` | 52 KB | 2,066 | about 38 MB |
| `bgc_cluster_gene` | 390 KB | 40,435 | about 290 MB |
| `bgc_mibig_hit` | 25 KB | not counted at 32 genomes (610 at top 3 for 12 of them) | about 19 MB |

Queries on the real tables took 1 to 3 ms. A synthetic copy repeated 750 times (24,000 pseudo-genomes, 1.55 million clusters) in memory took 9 to 62 ms for a filter, a join and a group. The copies repeat the same 32 genomes, so this tests speed only.

## 10. Recommendation (PROPOSED, for the owner)

1. Run lite without `--cb-general` for BFD. On three genomes the time fell 3.8x to 6.9x with identical regions and MIBiG tables. I have not measured the mean without it across genomes, so no 24,000-genome figure is given. Test one yeast, one basidiomycete and one Mucoromycete first.
2. Do not add a contig-length filter.
3. Add `--no-zip-output --no-summary-gbk --no-region-gbks` for the scratch saving (about half of the raw size), knowing that it does not save time.
4. Keep `known` and `known_pairs` at top 3, or top 5 if the owner wants the extra 1.3% of regions.
5. Keep the key as `LOCUSTAG` and the three-table parquet design.
6. Test one non-Pezizomycotina genome (a yeast, a basidiomycete, a Mucoromycete) before a full run.

## 12. Non-Pezizomycotina genomes and the genetic code (MEASURED, 4 genomes, top 5 hits kept)

One genome each from three groups outside the Pezizomycotina, picked as the genome nearest the group median size among those with BUSCO 95 or higher and N50 of 500 kb or higher (the *D. hansenii* BUSCO is about 84 for every genome of that species against `fungi_odb12`). Each ran without and with `--cb-general`.

| Genome (LOCUSTAG) | BUSCO | Mb | Contigs | Wall time without / with | Peak memory without / with | Regions | `known` rows (top 5) | `known_pairs` rows |
|---|---|---|---|---|---|---|---|---|
| *Blastobotrys adeninivorans* TMCC 70007 (F634A7E4), yeast, code 1 | 95.0 | 12.1 | 16 | 82 s / 144 s | 1.6 / 7.5 GB | 4 | 0 | 0 |
| *Debaryomyces hansenii* CBS767 (F8C91A38), yeast, code 12 | 83.8 | 12.2 | 8 | 86 s / 247 s | 2.3 / 10.6 GB | 5 | 0 | 0 |
| *Flammulina filiformis* HNY6 (F67CDA18), Basidiomycota | 98.5 | 36.3 | 24 | 291 s / 787 s | 4.6 / 10.8 GB | 38 | 4 | 7 |
| *Mucor hiemalis* (F00CDB58), Mucoromycota | 98.0 | 35.8 | 26 | 147 s / 437 s | 3.7 / 9.3 GB | 11 | 1 | 1 |

- For all four genomes the regions, `known` and `known_pairs` tables are identical with and without `--cb-general`. The speed-up is 1.8x, 2.9x, 2.7x and 3.0x. Peak memory without `--cb-general` is 1.6 to 4.6 GB, against 7.5 to 10.8 GB with it.
- Yeasts give few regions (4 and 5) and no MIBiG hit at all. That is a result, not a failure: `--taxon fungi` works and the regions are found. MIBiG is mostly filamentous-fungus clusters, so yeast clusters have little to match.
- **Genetic code.** BFD records `TRANSL_TABLE` 12 (CUG codes Ser) for 1,164 genomes and table 26 for 1. antiSMASH 8.0.4 has no option for it. With the BFD GFF3, antiSMASH sets `transl_table` to 1 and translates CUG as Leu. For *D. hansenii* CBS767, 3,923 of 5,951 CDS (65.9%) differ from the BFD protein, all with the same length, and all 11,859 differing residues are Ser in BFD and Leu in antiSMASH. For the other three genomes the translations match, except 0 of 6,090 (*B. adeninivorans*), 8 of 14,976 (*F. filiformis*) and 7 of 10,434 (*M. hiemalis*); the cause of those 15 (CDS translated much shorter than the BFD protein) was not investigated.
- **Fix, MEASURED on one genome:** adding `;transl_table=12` to the attributes of every CDS line of the GFF3 before antiSMASH gives 5,951 of 5,951 identical translations and the same 5 regions, in 92 s against 86 s. The recipe is in `run_antismash_one.sbatch` (`TRANSL_TABLE`). Effect of the wrong code on regions and MIBiG hits was not measured separately: for *D. hansenii* the regions are the same with either code, and there are no MIBiG rows in either case.
- PROPOSED: the module should read `TRANSL_TABLE` from the sample sheet and add the attribute when it is not 1.

## 11. Mistakes made on the way (so they are not repeated)

- The first EBC table used the rank-1 hit only. antiSMASH rank is not similarity order. Corrected in message 011.
- A cancel command with comma-separated job names did not cancel; 9 jobs ran to the end with the old extractor and were overwritten. Check `squeue` after `scancel`.
- The first module test failed because my test script's `projectDir` had no `bin/`. It was a test error, not a module error.
- I committed module changes into the shared `Fungi_BFD` checkout before another session asked me not to touch it. They are listed in the message to that session.
