# UPO reference set and source checks

Checked 2026-11-07. "Verified" means I retrieved the record from the named source on that date.

## 1. Primary papers

- Kinner A, Rosenthal K, Lütz S (2021). Identification and expression of new unspecific peroxygenases — recent
  advances, challenges and opportunities. *Front. Bioeng. Biotechnol.* 9:705630.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8293615/
- Faiza M, Lan D, Huang S, Wang Y (2019). UPObase: an online database of unspecific peroxygenases.
  *Database* 2019:baz122. https://pmc.ncbi.nlm.nih.gov/articles/PMC6902001/

## 2. What the papers say (verified by reading)

Kinner et al. 2021:
- UPObase holds "more than 1,900 putative UPO protein sequences" (cites Faiza 2019b).
- Family I (short-type): >1,100 members, "occurring in most fungal phyla", mean ~30 kDa, mainly without
  predictable signal peptides. Exceptions with SP: CglUPO, MroUPO, MweUPO, LfuCPO.
- Family II (long-type): >900 members, mean 44.4 kDa, only Ascomycota and Basidiomycota, most with a C-terminal
  disulfide bridge and a predicted signal peptide. AaeUPO is the model.
- Motifs: "-EHD-S-E-" (Family I) and "-EGD-S-R-E-" (Family II), "along with the characteristic PCP motif in both
  families".
- Genome-mining sentence: "genome mining of more than 800 fungal genomes and identified 113 putative sequences in
  35 fungal species" (cites Faiza 2019a).
- 17 UPOs purified and characterized at the time (9 ascomycete, 8 basidiomycete).
- Putative UPOs in its Figure 1 come from matching PROSITE PS51405 (heme haloperoxidase profile) to UniProtKB.
- The paper does **not** mention Pfam PF01328, Mucoromycota, Glomeromycota, Chytridiomycota or early-diverging fungi.
- The paper does not name a glutamate acid-base residue or expand "PCP".

Faiza et al. 2019 (UPObase):
- 1,948 sequences (representatives of 95 % cd-hit clusters) from **812 genomes from Ensembl Fungi**.
- Seeds: **113 previously identified UPO sequences** (including AaeUPO, MroUPO, LfuCPO). Search: phmmer 3.1b2,
  E 10, inclusion E 0.01. Then cd-hit 90 %, MCL (I = 1.4), cd-hit 95 %, motif filters, ML phylogeny (MEGA7, WAG+G+F).
- Classes: UPO subfamilies I–V = 123 sequences (Ia 54, Ib 16, II 3, IIIa 29, IIIb 5, IV 6, V 10).
  Peroxidase-peroxygenase (Pog) superfamily = 1,825 (Pog-a..g; Pog-g = 1,509, no motif beyond the heme Cys).
- Motif labels in this paper: "PCP-EGD-R-E" = UPO-type; "PCP-EHD-E" = CPO-type. Per-subfamily motif strings are
  listed in the paper (for example subfamily II: `PCP NHG RGN S[IL]G VPPLPG IDG HGRF EGD SMTRXD RX{6}E TX{6}R`).
- No taxonomic breakdown. No early-diverging phyla named. No accession numbers.
- Access stated: upobase.bioinformaticsreview.com, with per-entry FASTA download.

## 3. Conflicts and corrections to BACKGROUND.md

1. **"113 putative UPOs from >800 genomes" is probably a mis-citation.** In the UPObase paper, 113 is the number of
   *seed* sequences; the >800-genome search produced 1,948. Kinner 2021 attributes "113 in 35 species from >800
   genomes" to Faiza 2019a, which I did not read. BACKGROUND.md repeats the Kinner sentence. Do not use "113 from
   >800 genomes" as a positive-control target until Faiza 2019a is checked.
2. **Family I vs II is not the same split as UPObase I–V vs Pog.** Kinner's 1,100 / 900 split (from Kellner 2014,
   Hofrichter 2015) and UPObase's 123 UPO / 1,825 Pog split are different classification schemes. The plan must
   say which one it uses.
3. **EHD vs EGD labels disagree.** Kinner: EHD = Family I (short UPOs). UPObase: EHD = CPO-type, EGD = UPO-type.
   The crude scan in CENSUS.md finds EHD in nearly all early-diverging hits. Under Kinner this reads "short-type UPO";
   under UPObase it reads "CPO-like". This must be resolved against structures (5FUJ MroUPO, 1CPO LfuCPO, 2YOR AaeUPO)
   before any hit is called a UPO.
4. **Zámocký et al. 2025 claim not verified.** BACKGROUND.md cites it for "Mucoromycota, Glomeromycota and
   Chytridiomycota remain insufficiently investigated". The TU Delft repository refused the connection
   (ECONNREFUSED, 2026-11-07). Kinner 2021 does not make this claim.

## 4. Pfam / InterPro / PROSITE (verified via InterPro API)

- **PF01328** "Peroxidase, family 2" (Peroxidase_2), type family. Integrated into **IPR000028 Chloroperoxidase**,
  together with **PROSITE PS51405** (heme haloperoxidase family profile).
- InterPro text: proximal heme ligand is Cys (as in P450); distal side polar; no distal Arg; the catalytic acid-base
  is a **glutamate**, not a histidine.
- PF01328 is broad: BFD has 97,520 PF01328 proteins in total (CENSUS.md). It captures CPO, UPO and Pog members.
  A PF01328 hit is a family call, not a UPO call.
- URL: https://www.ebi.ac.uk/interpro/entry/pfam/PF01328/ ; https://www.ebi.ac.uk/interpro/entry/InterPro/IPR000028/

## 5. UPObase availability (measured)

- `upobase.bioinformaticsreview.com`: DNS does not resolve from the HPCC login node (curl error 6, 2026-11-07,
  both http and https). Other hosts (ebi.ac.uk, rest.uniprot.org, rcsb.org) resolved from the same node.
- Wayback Machine availability API returned HTTP 429 (rate limit). Not checked further.
- Status: **UPObase is offline or moved. Not obtainable now.** Fallback options are in PLAN.md section 3.

## 6. Characterized reference enzymes with verified accessions

Source: UniProt REST and RCSB PDB REST, 2026-11-07. "PE" = UniProt protein existence.

| Short name | Organism | UniProt | Length | PDB (examples) | Note |
|---|---|---|--:|---|---|
| AaeUPO (AaP) | *Cyclocybe (Agrocybe) aegerita* | B9W4V6 (APO1_CYCAE, Swiss-Prot) | 371 | 2YOR, 2YP1, 5OXU | Long type; model enzyme |
| CraUPO | *Coprinellus radians* | B9W4V8 (APO1_COPRA, Swiss-Prot) | 261 | – | Characterized from native producer (Kinner Table 1) |
| LfuCPO | *Leptoxyphium fumago* | P04963 (PRXC_LEPFU, Swiss-Prot) | 373 | 1CPO, 2CPO | Chloroperoxidase; outgroup / CPO reference |
| MroUPO | *Marasmius rotula* | A0A023H437 (PE 1) | **67** | 5FUJ, 5FUK | UniProt entry is a 67-aa fragment. Full sequence must come from PDB 5FUJ or 7ZBP (239 aa, no UniProt link) |
| MweUPO | *Marasmius wettsteinii* | A0ABF7PQ10 (PE 1) | 240 | 8RNJ–8RNR | Short type |
| TteUPO | *Thermothielavioides terrestris* NRRL 8126 | G2R8M6 (PE 1) | 274 | 9IZ8 | |
| DcaUPO | *Daldinia caldariorum* | A0AAN0NP18 (PE 1) | 245 | 8IAG | |
| HspUPO | *Hypoxylon* sp. EC38 | A0A1Y2TH07 | 261 (PDB construct) | 7O1R, 7O2G | |
| PabUPO-II | *Candolleomyces (Psathyrella) aberdarensis* | A0A4V1Q376 | 335 (PDB construct) | 9HE6, 9IFP | Variant 9RFY maps to A0A4Q2DF39 |
| CviUPO | *Collariella virescens* | none linked | 226 (PDB construct) | 7ZCL | Sequence from PDB only |
| UPO1 | *Calonarius (Cortinarius) odorifer* | P0DYC0 (Swiss-Prot) | 381 | – | Biosynthetic UPO (phlegmacins cluster) |
| DchUPO | *Daldinia childiae* | none linked | 279 (PDB construct) | 9TM7 | |
| DspUPO | *Daldinia* sp. EC12 | none linked | 227 (PDB construct) | 9WA6 | |

Not found by my searches (accession **unverified**): CglUPO (*Chaetomium globosum*), CciUPO (*Coprinopsis cinerea*),
MthUPO / MfeUPO / MhiUPO (*Myceliophthora*), HinUPO (*Humicola insolens*), CveUPO (*Coprinus verticillatus*).
Their sequences are in the cited expression papers (for example Kiebist et al. 2017 for CglUPO; Püllmann et al. 2021
for the *Myceliophthora*/*Humicola* set). These papers were not opened in this session.

Excluded from the list: RCSB "unspecific peroxygenase" text hits that are P450 BM3 (P14779) or bacterial SyoA
(A0AAX7FLX8). They are peroxygenase-mode P450s, not heme-thiolate UPOs.

## 7. Sources

- https://pmc.ncbi.nlm.nih.gov/articles/PMC8293615/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6902001/
- https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01328
- https://www.ebi.ac.uk/interpro/api/entry/interpro/IPR000028
- https://rest.uniprot.org/uniprotkb/B9W4V6 (and the other accessions above)
- https://www.rcsb.org/structure/5FUJ , /7ZBP , /7ZCL , /7O1R , /9HE6 , /9TM7 , /9WA6 , /8IAG , /9IZ8
