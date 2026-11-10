# Zamocky et al. 2022 reference set

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

## Steps 2-4 results: dedup, QC, in-house tree, clade names (job 31895597, 2026-11-08)

Driver: `../../refset/slurm/build_tree.sh`. Outputs in `tree/`. Code: `build_refset.py` (combine,
dedup, QC), `tree_outline.py` (outline / backbone ladder), `name_clades.py` (naming), all with
offline tests in `../../refset/test/`.

**Dedup (`tree/dedup.tsv`).** By sequence (identical, or one contained in the other). 5 of the 173
drop: HspHTP1 = HspUPO, TterHTP1 = TteUPO, AaeAPO = PaDa-I parent (backbone copy kept), and
BfucHTP1 = BcinHTP, ThetHTP1 = MtherHTP (within the set). MrotHTP1 (5FUJ chain A) is *not*
identical to backbone MroUPO; the two are sisters at UFBoot 100. **Tree input: 187 = 19 backbone +
168 Zamocky.**

**QC (`tree/qc.tsv`).** All 187 hit PF01328 at E <= 1e-5; none excluded. 66 flagged (flags do not
exclude): 65 `catalytic_residues_differ`, 2 length vs Table S1 (AorHTP2, CfumagHTP, above).

**Alignment and tree.** MAFFT L-INS-i, trimAl `-gt 0.75` (Zamocky's 75 % site coverage) -> 224
columns (theirs: 239). IQ-TREE 3.0.1, 1000 UFBoot + 1000 SH-aLRT, seed 20261108:

| run | model | lnL |
|---|---|---|
| `mfp` (the tree we use) | Q.PFAM+R7, ModelFinder/BIC | -53147.1 |
| `lgig` (check) | LG+I+G4, Zamocky's model | -53357.2 |

**Clade membership is identical under both models for all 187 leaves.**

### Zamocky Figure 1 clades, recovered

Rooted on the ten Mucorales sequences. Each clade = MRCA of anchors from the paper's *text*
(`fig1_anchors.tsv`), grown to the largest UFBoot >= 95 ancestor with no other clade's anchors.
The length check is independent of the anchors: the paper states a mean length per clade.

| Fig 1 clade | n | UFBoot (mfp / lgig) | mean aa (ours / paper) | composition |
|---|---|---|---|---|
| Mucoromycetous basal (root) | 10 | 100 / 100 | 301 / 301 | Mucorales 10 |
| Basidiomycetous 1st | 34 | 100 / 100 | 315 / 316 | Basidio 33, Asco 1 (PnoHTP07; the paper's "one short ascomycetous HGT") |
| Glomeromycetous ("Mucoromycetous small") | 6 | 94 / 93 | 230 / 230 | Gigaspora, Glomus, Rhizophagus |
| Ascomycetous 1st | 8 | 100 / 100 | 256 / 265 | Asco 8 (anchors from tree) |
| Basidiomycetous 2nd | 15 | 100 / 100 | 279 / 279 | Basidio 15 (anchors from tree) |
| Stramenopilous 1st | 4 | 100 / 100 | 260 / 276* | Saprolegnia, Aphanomyces |
| Stramenopilous 2nd | 12 | 99 / 100 | 258 / 276* | Phytophthora |
| Chytridiomycetous | 2 | 100 / 100 | 378 / - | Spizellomyces |
| Amoebozoan | 2 | 100 / 100 | 316 / - | Protostelium fungivorum |
| Basidiomycetous 3rd = **Family II, long** | 31 | 97 / 100 | 403 / 402 | Basidio 13, Asco 13, **our Clade I** (all 5) |
| Ascomycetous 2nd = **Family I, short** | 44 | 65 / 85 | 265 / 265 | Asco 31, Basidio 1, **our Clades III, IV, V** |

\* 276 is the paper's mean for the whole oomycete family, not per clade.

Where it differs from Zamocky: **the deep branching order is unresolved** (backbone UFBoot 8-60 in
both runs) and is not the same in the two runs, so the order of clades along the backbone is not
claimed. The terminal clades are robust. In our tree Ascomycetous 1st branches first after the
Mucorales; in the paper's text the Basidiomycetous 1st and glomeromycetous clades come first.

Not in any Figure 1 clade: 19 leaves, all Asco or backbone. 7 are our Clade II group (CviUPO,
PanUPO + 5 Asco), at the base of Family I (in a 60-leaf node with it, UFBoot 93 mfp) but outside the
anchors' MRCA; whether Zamocky's Family I includes them is not stated in the text. The other 12
(AfumHTP3, AniHTP5, AorHTP4, CgHTP4, FagHTP1, GzHTP2, MgrHTP1, MoHTP1, MscoHTP, NcHTP1, PsubHTP,
SAcoHTP2) fall on unsupported backbone branches.

**Table S1 error:** the three *Rhizophagus* sequences are listed as Mucoromycotina; *Rhizophagus*
is Glomeromycotina, and they sit in the glomeromycete clade. Their tree labels still say `Mucoro`.

### Our Clades I-V in the expanded tree

| clade | n | UFBoot (mfp / lgig) | Zamocky members added |
|---|---|---|---|
| I (long) | 31 | 97 / 100 | 26: all of Family II, including 13 ascomycete long HTPs |
| II | 7 | 98 / 98 | AniHTP4a, AniHTP4b, OpsiHTP, PnoHTP04, PnoHTP05 |
| III | 3 | 100 / 100 | MrotHTP1 |
| **IV** | **20** | **96 / 97** | **19**: AniHTP1-3, AorHTP2, AteHTP1, BcinHTP, BnarHTP1, BtulHTP1, CchloHTP1, CgHTP1, GzHTP1, PficHTP1, PnoHTP01, PnoHTP03, SborHTP1, ScepHTP1, SmusHTP, SsclHTP1, StriHTP1 |
| V | 14 | 86 / 85 | AnHTP2, AorHTP1, CgHTP2, CgHTP3, MtherHTP |

- **Clade IV is a real, supported clade**, not a one-member artefact: MphUPO plus 19 curated
  ascomycete sequences (Botrytis, Sclerotinia, Aspergillus, Fusarium, ...). This settles the
  topology question behind `+provisional_cladeIV`. It does not settle activity: MphUPO is still the
  only characterized member.
- Clades III, IV and V all fall inside Zamocky's Family I; Clade I is exactly Family II.
- Clade V itself is only moderately supported (UFBoot 85-86).
- 7 Family I sequences fall in none of our clades (Colletotrichum ChigHTP1/2, CfioHTP1,
  ColgloeoHTP; CfumagHTP = CPO; HjHTP1; PnoHTP02): candidate additional short sublineages.

### What this says about the current 19-enzyme caller (`upo_sites`)

Run on the Zamocky sequences, the caller agrees with the tree for **48/56** sequences that fall
inside Clades I-V. The 8 disagreements are 6 Clade IV members called V (identity 45-52 %), OpsiHTP
(II called III) and PnoHTP09 (I, below the 30 % floor).

**But 112 of the 168 sit outside all five clades, and the caller gave 81 of them a clade anyway**
(all but one in the `extrapolated` tier): the whole Basidiomycetous 1st and 2nd, oomycete,
Mucorales, glomeromycete, chytrid and amoebozoan clades are called II, III or V by nearest
neighbour. In BFD the `extrapolated` tier therefore contains lineages that are not in any of our
clades, which nearest-neighbour cannot detect and a placement on this tree can.

### Catalytic position 180 (HspUPO numbering) by clade

| clade | residue at 180 |
|---|---|
| Basidiomycetous 1st | **N x29, S x5; no Glu in 34/34** (the paper noted Glu->Asn in AbHTP1 only) |
| Family II (Clade I) | E x18 (incl. all 5 backbone), **A x10, S x2, T x1: all 13 ascomycete long HTPs** |
| every other clade | E, with rare Q |

Two consequences, both hypotheses until checked against a structure-based alignment of long UPOs:
(1) the low Glu180 rate in BFD Clade I calls (16.8 %, `G1_BFD_SUMMARY.md` §6, read there as
non-UPOs in the low-identity bin) may instead be genuine ascomycete long HTPs, which carry A/S/T
here in all 13 curated examples; (2) a whole early basidiomycete clade lacks the acid-base Glu.
