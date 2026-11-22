# antiSMASH output: the subset an MS-feature analysis needs

Status: working note, 2026-11-21. Nothing here is a biological finding. Labels: MEASURED = I measured it on the Af293 result described below; PROPOSED = my choice; UNVERIFIED = not checked.

## Why a subset

The one existing BFD-genome antiSMASH result is Af293 (antiSMASH 8.0.4, 42 regions, in `/bigdata/stajichlab/shared/projects/1KFG/common_annotate/annotate/Aspergillus_fumigatus_Af293/antismash_local/`). It is **228 MB** (MEASURED, `du`).

| Part | Size (MB) | Needed for MS-feature analyses? |
|---|---|---|
| `Aspergillus_fumigatus_Af293.json` | 68 | Source for the subset; then drop |
| `Aspergillus_fumigatus_Af293.gbk` (annotated copy of the whole genome) | 61 | No (BFD already holds the genome and annotation) |
| `Aspergillus_fumigatus_Af293.zip` | 38.5 | No (a copy of the rest) |
| `clusterblast/` (42 text files) | 21 | No |
| `regions.js` + `index.html` + `js` + `css` + `images` | 17 | No (HTML viewer) |
| 42 per-region `*.regionNNN.gbk` | about 22 | No (regions are in the tables below) |
| Subset tables (below) | 0.1 | Yes |

Inside the 68 MB JSON (MEASURED): genome sequence 29 MB; `PFAM_domain` features 8.5 MB; `CDS` features 7.4 MB (9,161 CDS, of which 819 lie inside regions); `modules.clusterblast` 18 MB; `modules.full_hmmer` 5.4 MB. The regions, NRPS/PKS domains, modules and predictions hold well under 1 MB.

## What the MS analyses need

From the EBC idea document (`IDEA_compound_to_gene_cluster_linking.md` in `EBC-KT`) and the pathway census, the analyses are: a species-by-cluster matrix, class-level linking (predicted cluster class against CANOPUS or NPClassifier class), and known cluster-to-compound pairs. They need:

| Need | Table or file | Source field |
|---|---|---|
| Cluster class and size, contig-edge flag (to drop fragments) | `regions.tsv.gz` | region products, protocluster category, `contig_edge`, core location |
| Mapping to BFD genes and to MIBiG loci | `regions.tsv.gz`, `cds.tsv.gz` | coordinates, locus tag, protein id, gene kind, function text |
| Cluster families without BiG-SCAPE (protein-level clustering of core genes) | `biosyn_proteins.faa.gz` | translations of genes with kind `biosynthetic` inside regions |
| Peptide and polyketide composition, a route to predicted mass | `domains.tsv.gz`, `modules.tsv.gz`, `adomain.tsv.gz` | NRPS/PKS domains and modules; A-domain substrate predictions (Stachelhaus top match with scores, physicochemical and cluster class) |
| Similarity to characterised clusters | `blast_top.tsv.gz` (general ClusterBlast, top 3) | **MIBiG similarity is not in the existing run.** It needs a new run with the known-cluster comparison on (see below) |

The existing Af293 run has only the general ClusterBlast, whose best hits are other fungal genome scaffolds (for example region 1: *Aspergillus tubingensis* scaffold, similarity 25), not MIBiG entries. It is of little use for compound linking.

## The extractor

`extract_antismash_subset.py INPUT.json OUTDIR GENOME_ID` writes seven files per genome: `regions`, `cds` (in-region only), `domains`, `modules`, `adomain`, `blast_top` (gzip TSV) and `biosyn_proteins.faa.gz`.

MEASURED on Af293: 42 regions, 819 CDS, 250 domains, 60 modules, 44 A-domain predictions, 65 biosynthetic proteins, **104 KB in total** (about 2,200 times smaller than 228 MB). Not yet written: a parser for the known-cluster (MIBiG) results, because the existing JSON has none to read. Write it from the first result that has them.

Extrapolation (PROPOSED, UNVERIFIED because one genome was measured and region count varies by genome): 100 genomes about 10 MB for the subset against about 23 GB raw; 425 genomes about 44 MB against about 97 GB.

## How to run (PROPOSED; flags UNVERIFIED because the command was not on the path on the login node after `module load antismash/8.0.4`)

1. Run in `$SCRATCH`, one genome per task. Input: the BFD GenBank file `genome_annotation/<name>/predict_results/<name>.gbk.gz` (decompress into `$SCRATCH`). Use fungal mode, no gene finding (BFD already has gene models), and the known-cluster comparison on. Check the exact flags with `antismash --help` first.
2. Skip what the subset does not use: the zip file and the HTML if antiSMASH offers a switch, and the general ClusterBlast.
3. Extract the subset into `/bigdata`. Delete the raw directory from `$SCRATCH` in the same job.
4. Time one run before any array. The Af293 JSON reports per-record times that sum to 7,905 s (UNVERIFIED how many CPUs). Size the array so each task takes about 1 to 1.5 hours.

## Which genomes (PROPOSED options, MEASURED counts from the BFD species table and the EBC corrected species table `ebc-021`)

The EBC side has metabolomics for 494 species. Exact name match gives 188 species with BFD genomes (4,802 genomes in all).

| Option | Species | Genomes |
|---|---|---|
| A. One genome per matched species | 188 | 188 |
| B. Up to 3 per matched species | 188 | 425 |
| C. Attributable single-species groups only, up to 3 each | 28 | 80 |
| D. As C, without the yeasts and other non-filamentous species listed by name | 20 | 56 (at most 3 each) |
| E. The 11 attributable species with two or more accessions | 10 with genomes | 29 (at most 3 each) |

My pick (PROPOSED): start with option D or E (56 or 29 genomes, a few MB of subset), because the EBC gate shows that only about 10 secondary-metabolism species have usable metabolomics. Choose genomes by assembly quality (BUSCO, contig count, from `busco_genome` and `asm_stats`), and for *A. fumigatus* include Af293 so that the new run can be compared with the existing one.

## Open

- The antiSMASH command and flags. The timing run. The known-cluster parser.
- Whether to keep raw JSON for a small reference set (for example 3 genomes, 200 MB) so that a field missing from the subset can be added later without a rerun.
