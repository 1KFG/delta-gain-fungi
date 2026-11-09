# Zamocky et al. 2022 reference set: step 0 (source inventory)

Source: Zámocký M. et al. 2022, "Evolution of Heme Peroxygenases: Ancient Roots and Later Evolved
Branches", *Antioxidants* 11:1011, doi:10.3390/antiox11051011, PMC9138132 (CC BY 4.0).

`antioxidants-11-01011-s001.zip` is the unmodified supplement, from the PMC Open Access mirror
(`https://pmc-oa-opendata.s3.amazonaws.com/PMC9138132.1/antioxidants-11-01011-s001.zip`, fetched
2026-11-08). The PMC web page and MDPI both block scripted downloads; the AWS mirror is the
supported programmatic route. It contains:

- `Table S1 overview on peroxidase sequences abbreviations.xlsx`: 174 rows (abbreviation, ID,
  organism, phylum, length, PDB code, original annotation). **No clade or family column.**
- `Figure S1_20 HTPs aligned.pdf`: an alignment of **20** sequences, not the 174. There is no
  174-sequence alignment and no Newick tree in the supplement.

## `tableS1_resolved.tsv`

Sequences are fetched by `../../refset/bin/fetch_refset.py` (see `../../refset/slurm/fetch_refset.sh`) into
`fetched/`. Where a row has two possible sources, the one whose length matches Table S1 is kept.

One row per Table S1 entry, with the route to a sequence (`route` column):

| route | n | how step 1 fetches it |
|---|---|---|
| `uniprot_direct` | 73 | UniProt REST by accession. **15 are no longer in UniProtKB** (deleted as "not part of a reference proteome", or demerged); the script takes those from UniParc. (An earlier note here said all 73 were live; that check was wrong.) |
| `uniprot_xref` | 29 | PeroxiBase ID resolved through UniProt's PeroxiBase cross-reference (DB-0072) |
| `wayback_view_perox` | 59 | PeroxiBase ID with no UniProt link. From the Wayback Machine: the archived FASTA endpoint `tools/get_fasta/<id>/PEP` first (many entry pages do not inline the sequence), then the entry page `display_perox/view_perox/<id>`, on both the `inra.fr` and `inrae.fr` domains; capture closest to 2022-04-07, Zamocky's access date |
| `ncbi_efetch` | 8 | NCBI protein efetch (all 8 resolve) |
| `inferred_unique_PF01328_len_match` | 4 | Given ID unusable: `L2GBA` truncated; PeroxiBase 4233 not archived; PeroxiBase 12620 archived without sequence; PeroxiBase 6766 archived at 305 aa vs 386 in Table S1. Each resolves to the one UniProt PF01328 entry in that organism at the stated length (`L2GBA1`, `A0ACM9MD70`, `A0A0D0A0B1` by elimination against SlutHTP1/2, `B0DV60`). **Inferred, not given**; the `note` column says how |
| `pdb` | 1 | `5FUJ_A` is MroUPO, already in `clade_backbone/` |

Unresolved after two fetch runs (jobs 31880062, 31880552):

- **AteHTP2** (PeroxiBase 4192, *A. terreus*, 265 aa): archived without sequence, no UniProt entry
  at that length. **Missing: 173 of 174 recoverable.**
- **AorHTP2**: only a 205-aa UniProt gene model exists (Table S1: 259). Kept and flagged; step 3 QC
  decides.
- **CfumagHTP** (P04963, CPO): UniProt gives the 373-aa precursor; the mature chain (21-319, 299 aa)
  is what Table S1 lists. Kept as the precursor; the PF01328 alignment ignores the propeptide.

PeroxiBase IDs 4207 and 4222 (*Phytophthora infestans*) each map to two UniProt entries; keep one
after checking that the sequences are identical.

Overlap with the current 19-enzyme backbone, by accession: HspUPO, TteUPO, MroUPO (5FUJ_A), and
AaeAPO = PaDa-I parent (B9W4V6). Do the full deduplication by sequence in step 2, not by ID.

## Things step 2 must handle

- **No clade labels are published in machine-readable form.** Figure 1 (condensed tree, 11 named
  clades) and Figure 2 (circular tree, leaf labels unreadable at the published resolution) are
  the only clade information. **We rebuild the tree ourselves** (decided 2026-11-08): align the
  174 plus our 19 references, infer an ML tree with IQ-TREE (ModelFinder + UFBoot; their reported
  model was LG+G+I on 239 positions at 75 % site coverage), then name clades by matching their
  membership and phylum composition to Figure 1. The rebuilt tree is a check on Zamocky's topology
  as well as a source of labels; where the two disagree, report it rather than force a match.
- **Naming collision.** Zamocky's "Dikaryal Family I" is the *short* family and "Family II" the
  *long* one (the older Hofrichter usage). In our backbone, Clade I is *Long*. Never map by numeral.
- Taxon strings are inconsistent (`Ascomycota ` with a trailing space, `Amobeozoa`, three spellings
  of Mucoromycota); normalise them from NCBI taxonomy, not by hand.
- All 174 are curated, full-length, genome-predicted sequences, not characterized enzymes. Record
  `evidence = curated` for every row except those matching a characterized enzyme.
