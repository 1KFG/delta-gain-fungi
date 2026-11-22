# antiSMASH output: the subset an MS-feature analysis needs

Status: working note, updated 2026-11-21. Nothing here is a biological finding. Labels: MEASURED = I measured it; PROPOSED = my choice; UNVERIFIED = not checked.

## Note for other agents (read this first)

- **What was done (MEASURED):** antiSMASH 8.0.4 was run in the BFD apptainer image on 24 BFD genomes (Af293 twice, so 24 runs on 23 distinct genomes at the time of writing; the last jobs were still running). Results are in `analysis/bgc_subset/results/<id>.{regions,cds,known,...}.tsv.gz`. The job script is `run_antismash_one.sbatch`.
- **Lite mode is enough to enumerate clusters (MEASURED on one genome):** on Af293 the lite run (`--cb-general --cb-knownclusters`) found the same 44 regions, with the same coordinates and classes, as the full run (`--fullhmmer --clusterhmmer --cb-general --cb-knownclusters --pfam2go`). I did not test other genomes or other flag sets.
- **Time (MEASURED, 8 CPUs):** full 6,535 s on Af293. Lite on 32 genomes: 500 to 4,408 s, median 1,028 s (about 17 min). The slowest are the most fragmented assemblies (about 70 min at over 1,000 contigs).
- **Extra tables (2026-11-22):** `known_pairs` (our gene to MIBiG gene: identity, coverage), `n_entry_proteins`/`synteny_score`/`core_bonus` in `known`, and `run.json` (versions, mode). `build_ebc_table.py` uses the highest similarity of the top-3 hits (rank 1 is not always the highest), and adds entry fraction and MIBiG quality. 44% of `known` regions in 32 genomes have a best hit on a 1-protein MIBiG entry.
- **Size (MEASURED):** raw antiSMASH output 242 MB for Af293. The subset tables are about 2 MB per genome. The three-table minimal set used by BFD is smaller (UNVERIFIED size).
- **What the flags do (background, not tested here):** cluster detection does not depend on `--fullhmmer` or `--clusterhmmer`; they add Pfam annotation. The match to known clusters is `--cb-knownclusters` (MIBiG). `--pfam2go` adds GO terms. No analysis here uses GO.
- **BFD module (committed, not run through the workflow):** `Fungi_BFD/nextflow/modules/BFD/ANTISMASH/main.nf`, commit 75dc089 on branch `unify-funannotate-sif`. Default `antismash_mode = lite`; `full` restores the old flags and keeps the raw JSON. Switched on with `--run_antismash true` (default false). Outputs: `results/function/antismash/<bucket>/<LOCUSTAG>.antismash.{regions,cds,known}.tsv.gz`, hash-bucketed by locustag like cazy and pfam. Mode, bucket width and merge wiring are untested end to end.
- **Bug found in existing code:** `bin/merge_antismash.py` read `records[].regions`, a key antiSMASH 8 JSON does not have, so it would have produced an empty table. It now reads the lite tables.
- **Not decided:** the database representation of a cluster (location, annotation, genes). The three minimal tables are a proposal. Whether to run all of BFD (about 23,680 genomes, about 800 s each, about 5,300 hours at 8 CPUs, UNVERIFIED for genomes larger than the ones tested).
- **Limits of the results:** a known-cluster hit means similarity to a MIBiG entry. It does not show that the compound is made.

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

## How to run (MEASURED: this is what ran)

`sbatch -J <name> run_antismash_one.sbatch <genome_folder_name> <short_id>` (short partition, 8 CPUs, 24 GB, 2 h). The script works in `$SCRATCH`, runs `antismash --taxon fungi --cb-general --cb-knownclusters` with `--genefinding-gff3` on the BFD gff3 and genome fasta, reduces the JSON with `extract_antismash_subset.py`, and deletes the scratch directory. `KEEP_RAW=1` also keeps the gzipped JSON.

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

- The database representation of a cluster, and the BFD merge step on the new tables.
- A run of the BFD module inside the workflow (the module test used an older version).
- Whether to keep raw JSON for a small reference set. The Af293 raw files were deleted on 2026-11-21.
