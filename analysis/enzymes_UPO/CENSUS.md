# UPO target census: early-diverging fungi in BFD

Measured 2026-11-07 with duckdb 1.1.3 on `Fungi_BFD_runs/tables/` (snapshot dated 2026-11-05/06) and
`DeltaGain_Fungi/prior_annotation_status.tsv`. Queries are in `sql/`. Hit data are in `data/`.

## 1. Join between BFD and the prior-annotation table

- Key: `samples.ASMID` = `prior_annotation_status.folder` (`<accession>_<asm_name>`).
  A join on the accession prefix (`GC[AF]_nnnnnnnnn.n`) gives the same count.
- 23,656 of 23,680 BFD genomes match. The 24 unmatched genomes are 22 Ascomycota and 2 Basidiomycota.
  **The join loses no early-diverging genome.**
- Label used here: `NCBI` = status `baseline_ok` (NCBI annotation with a protein file exists).
  `BFD_only` = status `no_prior_annotation_*` (only BFD has gene models).
- One Mucoromycota genome has status `annotated_folder_missing_protein`. It is left out of both columns.

## 2. Genomes, species and assembly quality per phylum

`species` = distinct `samples.SPECIES` strings. `new species` = BFD_only species whose name string does not occur
in the NCBI-annotated set. String matching can over-count (for example "sp." names). BUSCO = `busco_genome.complete_pct`.

| Phylum | Ann. | Genomes | Species | New species vs NCBI | BUSCO C median | BUSCO C <80 % | Median size (Mb) | Median N50 (kb) |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| Mucoromycota | NCBI | 145 | 69 | – | 96.4 | 6 | 35.6 | 100.9 |
| Mucoromycota | BFD_only | 148 | 66 | 49 | 97.5 | 7 | 35.1 | 73.0 |
| Glomeromycota | NCBI | 43 | 25 | – | 91.1 | 9 | 140.3 | 33.7 |
| Glomeromycota | BFD_only | 17 | 9 | 6 | 90.6 | 4 | 108.9 | 15.7 |
| Chytridiomycota | NCBI | 81 | 59 | – | 87.3 | 22 | 27.1 | 101.8 |
| Chytridiomycota | BFD_only | 8 | 7 | 3 | 79.6 | 3 | 16.1 | 38.9 |
| Mortierellomycota | NCBI | 85 | 51 | – | 93.4 | 4 | 38.3 | 35.3 |
| Mortierellomycota | BFD_only | 15 | 10 | 4 | 95.7 | 1 | 39.3 | 2196.7 |
| Zoopagomycota | NCBI | 6 | 6 | – | 74.1 | 4 | 18.4 | 77.4 |
| Zoopagomycota | BFD_only | 5 | 5 | 5 | 73.2 | 4 | 16.1 | 23.2 |
| Basidiobolomycota | NCBI | 2 | 2 | – | 93.9 | 0 | 73.1 | 58.2 |
| Basidiobolomycota | BFD_only | 16 | 16 | 15 | 94.8 | 1 | 87.4 | 19.2 |
| Kickxellomycota | NCBI | 187 | 167 | – | 70.1 | 123 | 14.2 | 9.6 |
| Kickxellomycota | BFD_only | 3 | 3 | 2 | 79.5 | 2 | 27.1 | 1874.1 |
| Entomophthoromycota | NCBI | 4 | 4 | – | 75.5 | 3 | 535.5 | 121.0 |
| Entomophthoromycota | BFD_only | 8 | 8 | 5 | 74.4 | 6 | 56.9 | 31.8 |
| Blastocladiomycota | NCBI | 11 | 9 | – | 76.6 | 7 | 26.4 | 217.8 |
| Blastocladiomycota | BFD_only | 3 | 2 | 0 | 57.7 | 3 | 33.1 | 11.3 |
| Neocallimastigomycota | NCBI | 6 | 6 | – | 75.2 | 6 | 129.4 | 573.3 |
| Neocallimastigomycota | BFD_only | 20 | 4 | 3 | 48.1 | 20 | 64.8 | 5.3 |
| Monoblepharomycota | NCBI | 3 | 3 | – | 81.2 | 1 | 33.0 | 347.3 |

Small phyla (genomes, NCBI / BFD_only): Aphelidiomycota 3 / 1, Calcarisporiellomycota 1 / 0, Olpidiomycota 1 / 0,
Rozellomycota 2 / 0, Sanchytriomycota 0 / 2. Microsporidia: 86 / 77 (excluded, see PLAN.md).

Totals: 827 early-diverging genomes excluding Microsporidia (samples.parquet). 813 of them have rows in
`gene_proteins.parquet`, with 10,166,518 proteins in total.

## 3. Existing Pfam PF01328 (Peroxidase_2) hits in BFD

Source: `pfam.parquet`, which BFD built with `hmmsearch --cut_ga` (Pfam gathering threshold;
`Fungi_BFD_runs/nextflow/modules/BFD/PFAM/main.nf` line 55). All PF01328 rows have domain i-Evalue <= 1.2e-5.
These are BFD gene models in both columns. For `NCBI` genomes, the NCBI protein set was **not** searched here.

Whole BFD: 98,159 PF01328 domain rows in 97,520 proteins (all phyla).

| Phylum | Ann. | Genomes with Pfam rows | PF01328 proteins | Genomes with >=1 hit | Median hits/genome | Max |
|---|---|--:|--:|--:|--:|--:|
| Mucoromycota | NCBI | 145 | 193 | 127 | 1 | 3 |
| Mucoromycota | BFD_only | 146 | 177 | 125 | 1 | 2 |
| Glomeromycota | NCBI | 42 | 92 | 40 | 2 | 8 |
| Glomeromycota | BFD_only | 17 | 32 | 17 | 2 | 3 |
| Chytridiomycota | NCBI | 78 | 184 | 53 | 2 | 13 |
| Chytridiomycota | BFD_only | 5 | 4 | 3 | 1 | 2 |
| Zoopagomycota | NCBI | 6 | 12 | 5 | 2 | 4 |
| Zoopagomycota | BFD_only | 5 | 18 | 5 | 1 | 11 |
| Basidiobolomycota | NCBI | 2 | 9 | 2 | 4.5 | 5 |
| Basidiobolomycota | BFD_only | 16 | 77 | 16 | 5 | 7 |
| Monoblepharomycota | NCBI | 3 | 4 | 2 | 1 | 3 |
| Aphelidiomycota | NCBI / BFD_only | 1 / 1 | 2 / 1 | 1 / 1 | – | – |
| Calcarisporiellomycota | NCBI | 1 | 2 | 1 | 2 | 2 |
| Mortierellomycota | both | 100 | **0** | 0 | 0 | 0 |
| Kickxellomycota | both | 177 | **0** | 0 | 0 | 0 |
| Blastocladiomycota | both | 14 | **0** | 0 | 0 | 0 |
| Neocallimastigomycota | both | 25 | **0** | 0 | 0 | 0 |
| Entomophthoromycota | both | 9 | **0** | 0 | 0 | 0 |
| Microsporidia | both | 133 | **0** | 0 | 0 | 0 |
| Ascomycota | NCBI / BFD_only | 4,019 / 15,187 | 20,239 / 57,971 | 3,065 / 8,980 | 4 / 2 | 55 / 42 |
| Basidiomycota | NCBI / BFD_only | 881 / 2,359 | 6,119 / 12,188 | 735 / 1,936 | 4 / 4 | 144 / 91 |

Early-diverging PF01328 proteins in total: 807 (BFD_only genomes: 309).

## 4. Hit properties (early-diverging, 807 proteins)

`has SP` = protein has a row in `signalp.signal_peptide.parquet` (SignalP 6, positives only).
`HMM cov >=80 %` = (hmm_to - hmm_from + 1) >= 0.8 x 186 (PF01328 model length).

| Phylum | Ann. | n | Median length (aa) | <250 aa | 250-350 aa | >350 aa | % with SP | % HMM cov >=80 % |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| Mucoromycota | NCBI | 193 | 302 | 8 | 184 | 1 | 2.6 | 98.4 |
| Mucoromycota | BFD_only | 177 | 302 | 10 | 159 | 8 | 1.1 | 98.9 |
| Glomeromycota | NCBI | 92 | 225 | 86 | 4 | 2 | 3.3 | 90.2 |
| Glomeromycota | BFD_only | 32 | 233 | 31 | 1 | 0 | 0.0 | 93.8 |
| Chytridiomycota | NCBI | 184 | 263 | 72 | 65 | 47 | 56.5 | 89.7 |
| Chytridiomycota | BFD_only | 4 | 275 | 0 | 3 | 1 | 50.0 | 100 |
| Basidiobolomycota | BFD_only | 77 | 264 | 18 | 47 | 12 | 51.9 | 92.2 |
| Zoopagomycota | BFD_only | 18 | 328 | 2 | 8 | 8 | 55.6 | 88.9 |

### Crude motif scan (first-pass only, not a UPO call)

Method: exact substring `PCP`; then `EHD` or `EGD` anywhere after the first `PCP`. No alignment. This over- and
under-counts. It is a feasibility check only. Sequences: `data/early_pf01328_hits.fa.zst`; per-protein flags:
`data/motifs.tsv`.

| Phylum | Ann. | n | PCP | EHD after PCP | EGD after PCP |
|---|---|--:|--:|--:|--:|
| Mucoromycota | NCBI | 193 | 191 | 190 | 2 |
| Mucoromycota | BFD_only | 177 | 175 | 167 | 0 |
| Glomeromycota | NCBI | 92 | 78 | 54 | 0 |
| Glomeromycota | BFD_only | 32 | 29 | 17 | 0 |
| Chytridiomycota | NCBI | 184 | 142 | 121 | 7 |
| Chytridiomycota | BFD_only | 4 | 4 | 4 | 0 |
| Basidiobolomycota | BFD_only | 77 | 60 | 60 | 1 |
| Zoopagomycota | BFD_only | 18 | 15 | 15 | 0 |

Observation: almost all early-diverging hits carry the EHD pattern, not EGD. The sources disagree on what EHD means
(see REFERENCES.md section 3). Do not interpret this until an alignment-based motif check is done.
