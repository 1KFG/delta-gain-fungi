# Three observations from the UPO shortlist build

2026-11-11. The three observations recorded in `CANDIDATES.md`, checked in more detail. Code is in
`bfd_findings/bin/` (tests in `bfd_findings/test/`, 8 passing) and the drivers are
`bfd_findings/slurm/report.sh` (job 31994904) and `bfd_findings/slurm/codons.sh` (job 31995187).
Tables are in `bfd_findings/out/`, and `out/summary.md` holds the generated tables these numbers
come from.

| # | Observation as first recorded | What the closer look shows |
|---|---|---|
| 1 | Mucorales UPO gene models in BFD are mostly incomplete (292 of 315) | **The gene models are fine; BFD's start/stop flags are wrong.** 99.3 % of starts and 99.8 % of stops flagged missing are present in the genome. The cause was found and fixed in Fungi_BFD today (fb1d7cd), but `tables/` predates the fix. |
| 2 | The glomeromycete clade attracts divergent outsiders | Confirmed. All 121 Glomeromycota members sit on short branches; 585 of the 663 outsiders carry `long_pendant`. |
| 3 | About 2,000 ascomycete proteins place in Basidiomycetous 1st | 2,265 proteins in 2,190 genomes, almost all single-copy, in 47.5 % of Dothideomycetes and all 67 Orbiliomycetes genomes, and placed like the clade's basidiomycete members. This looks like an old ascomycete lineage, not contamination. Unresolved whether it is vertical or one early transfer. |

---

## 1. BFD start/stop-codon flags (for BFD)

### What was seen

The shortlist's eligibility rule 2 required a complete gene model by BFD's flags
(`gene_transcripts.parquet`: `is_partial`, `has_start_codon`, `has_stop_codon`). Mucoromycota and
Glomeromycota failed it wholesale:

| phylum | genomes | UPO models | UPO flagged incomplete % | genome-wide flagged incomplete % |
|---|---|---|---|---|
| Ascomycota | 12,061 | 78,382 | 9.4 | 12.7 |
| Basidiomycota | 2,677 | 18,327 | 7.9 | 14.8 |
| Mucoromycota | 252 | 370 | 83.0 | 48.9 |
| Glomeromycota | 57 | 124 | 75.0 | 60.8 |
| Chytridiomycota | 56 | 188 | 22.9 | 27.7 |

Across all 23,452 BFD genomes the median genome has 14.5 % of its models flagged incomplete, and
936 genomes are at 40 % or more (`out/bfd_genome_gene_model_qc.tsv`). Within Mucoromycota,
Mucorales genomes have a median of 54 % and Umbelopsidales 16 %.

Three things argued against real truncation. Flagged-incomplete Mucoromycota UPOs are as long as
complete ones (median 302 vs 309 aa) and have the same number of exons (median 2). Only 19 of 308
sit within 1 kb of a contig end; most are 10 kb or more from one. And every CDS length is a
multiple of 3.

### Check against the genome sequence

`bin/check_codons.py` takes each model's CDS segments from `gene_CDS.parquet` and splices the CDS
from the input assembly (`input_clean_genomes/<ASMID>.fa.gz`, upper-cased). It then reads the first
codon, the last codon and the codon just downstream of the CDS. Five sets were checked: 5,816
models in 1,755 genomes, with random samples drawn with a fixed seed. The table is
`out/codon_check.tsv`.

| set | models | flagged no start | ...ATG actually present | flagged no stop | ...stop actually present | flagged present, actually absent |
|---|---|---|---|---|---|---|
| all Mucoromycota UPOs | 370 | 210 | 210 | 261 | 260 | 0 |
| random Mucoromycota models | 2,966 | 708 | 701 | 931 | 929 | 0 |
| random Glomeromycota models | 998 | 350 | 348 | 433 | 432 | 0 |
| random Ascomycota UPOs | 500 | 31 | 31 | 15 | 15 | 0 |
| random Ascomycota models | 982 | 70 | 69 | 72 | 72 | 0 |
| **total** | **5,816** | **1,369** | **1,359 (99.3 %)** | **1,712** | **1,708 (99.8 %)** | **0** |

- Of the 2,440 models with any FALSE flag, 2,426 have both an ATG start and an in-frame stop. The
  true incomplete rate in these samples is therefore about 0.2 % (14 of 5,816), not 10-60 %.
- All 870 UPO translations match BFD's protein exactly, and none has an internal stop.
- The error never runs the other way: no model flagged complete lacks its codon.
- It is strand-independent: 668 plus-strand and 701 minus-strand models carry false flags.
- TAA stops are flagged missing more often (39.8 %) than TAG (18.8 %) or TGA (16.1 %), which is
  what case-sensitive matching on soft-masked, AT-rich sequence would produce.

### Cause (already fixed in code, not yet in the tables)

`Fungi_BFD/nextflow/bin/build_genestats_table.py` compared the first and last codons to `ATG` and
the stop codons case-sensitively. Against a soft-masked genome, a lower-case `atg` or `taa` was
recorded as missing and the model was marked `is_partial`. Fungi_BFD commit **fb1d7cd** (2026-11-11
12:00, "read the genome upper-case so start/stop checks work on soft-masked genomes") fixes this
with `Fasta(dna, sequence_always_upper=True)`. On *Scedosporium aurantiacum* it reduced false start
flags from 593 to 2.

**`Fungi_BFD_runs/tables/gene_transcripts.parquet` is dated 2026-11-08 and still carries the old
flags.** The numbers above are an independent, cross-phylum measure of what the rebuild should
correct. The per-genome rates above therefore reflect how much soft-masked sequence falls on gene
ends, mostly in Mucorales and Glomeromycota, not annotation quality.

### Suggested actions for BFD

1. Re-run the gene-stats step with fb1d7cd and rebuild `gene_transcripts.parquet`, plus anything
   derived from `is_partial`, `has_start_codon`, `has_stop_codon` or the CDS md5. The commit notes
   that `CDS_sequence` and `md5checksum` change for masked CDS.
2. Verify the rebuild by re-running `bfd_findings/slurm/codons.sh` against it. The expected result
   is that flagged-missing codons drop to about 0.2 % of models and the "actually present" columns
   go to about 0. `out/bfd_genome_gene_model_qc.tsv` is the per-genome snapshot from before the
   rebuild to compare against.
3. Any consumer that filtered on these flags before the rebuild silently dropped about 10 % of
   models BFD-wide and about half in Mucorales and Glomeromycota.

### Consequence for this project

- The shortlist's gene-model rule removed proteins that were in fact complete: most of the
  Mucorales and glomeromycete pools and about 10 % elsewhere (12 % of Clade IV).
- The 37 picks themselves are unaffected, since every model the flags call complete really is
  complete. Pools B4 and B5 were drawn from a much smaller eligible set than they should have
  been, though.
- Re-run the shortlist after the tables are rebuilt, or switch rule 2 to the sequence check now.
- The catalytic-residue statistics never used the flags, so they stand.

## 2. Glomeromycete clade members

Table: `out/glomeromycete_clade_members.tsv` (784 proteins).

| phylum / class | members | `long_pendant` | short pendant | median pendant |
|---|---|---|---|---|
| Glomeromycota | 121 | 0 | 121 | 0.04 |
| Ascomycota: Dothideomycetes | 360 | 341 | 19 | 1.17 |
| Ascomycota: Sordariomycetes | 86 | 71 | 15 | 1.19 |
| Ascomycota: Eurotiomycetes | 66 | 65 | 1 | 1.22 |
| Basidiomycota | 62 | 61 | 1 | 1.36 |
| Basidiobolomycota | 29 | 12 | 17 | 1.00 |
| Ascomycota: Orbiliomycetes | 22 | 3 | 19 | 0.97 |
| Ascomycota: Lecanoromycetes | 20 | 20 | 0 | 1.17 |
| Ascomycota: Leotiomycetes | 17 | 12 | 5 | 1.13 |
| Ascomycota: Geoglossomycetes | 1 | 0 | 1 | 0.85 |

- The clade is clean at its core. Every Glomeromycota protein attaches close to the references (median
  pendant 0.04), and nearly every outsider attaches far out, at a median pendant of about 1.2, past the
  leave-one-out threshold. The `long_pendant` flag separates the two cleanly.
- 78 outsiders fall below the flag. The notable ones are 19 of 22 Orbiliomycetes and 17 of 29
  Basidiobolomycota; their median pendant of about 1.0 sits just under the threshold.
- These are the ones that could be real non-glomeromycete members, or the same attraction just under the cut.
  Shortlist pick B4-2 (*Basidiobolus*, pendant 0.96) is one of them.
- For counting: report this clade as 121 Glomeromycota proteins plus 78 short-pendant others,
  not 784.
- The `gene_model` column in this table uses BFD's flags and carries the error from section 1.

## 3. Ascomycete proteins placed in Basidiomycetous 1st

Dump: `out/asco_in_basidiomycetous_1st.tsv` (2,265 rows, with taxonomy, genome, copies per genome,
placement statistics, catalytic residues and the old caller's call) and
`out/asco_in_basidiomycetous_1st.faa` (sequences).

### Where they are

2,265 of the clade's 6,378 BFD members are Ascomycota, in 2,190 genomes. 2,123 genomes carry
one copy, 60 carry two, and 7 carry three or four.

| class | proteins | genomes carrying | BFD genomes in class | carrier % |
|---|---|---|---|---|
| Dothideomycetes | 1,364 | 1,298 | 2,731 | 47.5 |
| Eurotiomycetes | 680 | 675 | 3,220 | 21.0 |
| Leotiomycetes | 77 | 75 | 435 | 17.2 |
| Orbiliomycetes | 69 | 67 | 67 | **100** |
| Sordariomycetes | 65 | 65 | 6,016 | 1.1 |
| Arthoniomycetes | 3 | 3 | 13 | 23.1 |
| Lecanoromycetes | 2 | 2 | 479 | 0.4 |
| other / unplaced | 5 | 5 | | |

- **Orders:** Eurotiales 670, Pleosporales 446, Mycosphaerellales 401, Botryosphaeriales 310,
  Venturiales 111, Cladosporiales 82, Orbiliales 69, Hypocreales 56, Helotiales 40 and
  Thelebolales 35.
- **Eurotiales genera:** *Aspergillus* 614, *Penicillium* 34, *Talaromyces* 16. This bears on
  decision D1, Eurotiales sourcing: these proteins are a different UPO lineage from the Clade IV and
  short-UPO material.
- **Dothideomycetes genera:** *Macrophomina* 310, *Parastagonospora* 185, *Hortaea* 127, *Venturia*
  105, *Cladosporium* 80, *Zymoseptoria* 70. The counts follow how heavily each genus is sampled in
  BFD.

### How they place, compared with the clade's basidiomycete members

Quartiles are Q1 / median / Q3.

| | members | pendant length | clade LWR | residue at 180 | `long_pendant` |
|---|---|---|---|---|---|
| Basidiomycota | 4,097 | 0.171 / 0.279 / 0.377 | 0.995 / 1.000 / 1.000 | Asn 87.4 % | 14 |
| Ascomycota | 2,265 | 0.170 / 0.229 / 0.401 | 0.694 / 1.000 / 1.000 | Asn 87.4 % | 0 |

- **Branch lengths and residue 180:** the ascomycete proteins attach as close to the references as
  the basidiomycete members do. They carry the clade's Asn180 at the same frequency, and their length
  is tight (275 / 280 / 290 aa).
- **Lower confidence:** the clade LWR is lower, with 1,217 of 2,265 at or above 0.95 and none below
  0.5. Placement is less sure which clade they belong to, though it never doubts strongly.
- **A second catalytic substitution:** 785 (35 %) have Gln instead of Glu at position 109, the E of
  the EHD motif, beside the Asn180 substitution. 1,174 have Cys39, Glu109 and Asp111 intact.
- **All are "unassigned" in our Clades I-V,** as expected, since Basidiomycetous 1st is outside
  them. The old nearest-neighbour caller had put 727 in Clade V and 196 in Clade III, so they
  account for much of the extrapolated Clade V calls that placement removed (Figure 3).
- **Gene models:** by BFD's flags 2,186 are complete. Given section 1, nearly all of the other 79
  probably are too.

### Reading

- **Against contamination:** the proteins are spread across thousands of independently sequenced
  genomes, almost always as a single copy and with a consistent residue profile. They are present in
  every Orbiliomycetes genome, and Orbiliomycetes is an early-diverging Pezizomycotina class.
- **Against a recent transfer:** the same evidence argues against one recent horizontal transfer
  into a single ascomycete lineage.
- **What it fits:** an old ascomycete UPO lineage, either inherited vertically and largely lost in
  Sordariomycetes and Lecanoromycetes, or acquired by one early transfer and then inherited.
- **Why placement cannot decide:** the 187-leaf reference holds the basidiomycete diversity of this
  clade and the single ascomycete Zamocky identified. Placement can say these proteins belong with
  that clade. It cannot say whether they nest inside the basidiomycete diversity (transfer) or form
  a sister group to it (old lineage).

**Test that would settle it:** build an ML tree of a dereplicated sample of this clade, about 150
ascomycete and 150 basidiomycete members at 70 % identity, with Basidiomycetous 2nd and
Ascomycetous 1st as outgroups. Then ask:

- whether the ascomycete proteins are monophyletic;
- where they attach;
- whether their internal topology follows the species tree (Orbiliomycetes first).

Not done here.

## Files

| file | contents |
|---|---|
| `out/summary.md` | generated tables (do not edit) |
| `out/gene_model_qc_by_phylum.tsv` | UPO vs genome-wide flagged completeness by phylum (pre-rebuild flags) |
| `out/bfd_genome_gene_model_qc.tsv` | per-genome flag counts for all 23,452 BFD genomes (pre-rebuild snapshot) |
| `out/mucoromycota_upo_gene_models.tsv` | each Mucoromycota/Mucorales-clade UPO model: coordinates, flags, contig length, distance to contig end |
| `out/mucoromycota_genome_qc.tsv` | per Mucoromycota genome: UPO and genome-wide flagged rates, N50 |
| `out/codon_check.tsv` | 5,816 models: BFD flags vs the codons in the assembly |
| `out/glomeromycete_clade_members.tsv` | all 784 glomeromycete-clade proteins |
| `out/asco_in_basidiomycetous_1st.{tsv,faa}` | the 2,265 ascomycete proteins in Basidiomycetous 1st |
