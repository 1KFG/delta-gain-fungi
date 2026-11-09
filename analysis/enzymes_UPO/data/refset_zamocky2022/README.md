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

One row per Table S1 entry, with the route to a sequence (`route` column):

| route | n | how step 1 fetches it |
|---|---|---|
| `uniprot_direct` | 73 | UniProt REST by accession (all 73 live in release 2026_04) |
| `uniprot_xref` | 29 | PeroxiBase ID resolved through UniProt's PeroxiBase cross-reference (DB-0072) |
| `wayback_view_perox` | 61 | PeroxiBase ID with no UniProt link; the archived PeroxiBase entry page (Wayback Machine, `peroxibase.toulouse.inra.fr/display_perox/view_perox/<id>`, 2024 captures) holds the full protein sequence, class and orthogroup. Checked on ID 10810: 371 aa, matches Table S1 |
| `ncbi_efetch` | 8 | NCBI protein efetch (all 8 resolve) |
| `inferred_unique_PF01328_len_match` | 2 | ID unusable (`L2GBA` truncated; PeroxiBase 4233 not archived). Each has exactly one UniProt PF01328 entry in that organism at the stated length (`L2GBA1`, `A0ACM9MD70`). **Inferred, not given**; flag them |
| `pdb` | 1 | `5FUJ_A` is MroUPO, already in `clade_backbone/` |

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
