# UPO expression candidates from under-characterized lineages

**Rules written 2026-11-11, before the selection was run.** Code: `candidates/bin/shortlist.py`,
driver `candidates/slurm/shortlist.sh`, tests `candidates/test/`. Results are appended below the
rules, not mixed into them.

## What "novel" means here

**Lineage novelty, not sequence novelty.** Every catalytic-core representative in BFD has a public
relative in nr, and apparent sequence novelty tracks whether the source genome's proteins were
deposited (G1 Finding 5, 5b). What *is* new is that whole supported clades of UPO-like proteins have
no characterized enzyme. A candidate is valuable because measuring it would be the first (or second)
measurement for its lineage.

## Rules

**Eligibility, all pools:**

1. Best placement not flagged `long_pendant` (farther from every reference than any correct
   leave-one-out call).
2. Complete gene model in BFD `gene_transcripts`: `is_partial` false, start and stop codons present.
3. Length within the 5th-95th percentile of the pool's otherwise-eligible proteins.
4. Pools A and B: full catalytic core (Cys39, Glu109, Asp111, Glu180; HspUPO numbering) and
   placement LWR >= 0.95 for the pool's clade.
5. Pool C: Cys39, Glu109, Asp111 and LWR >= 0.95. Glu180 is not required, because its absence is
   the lineage trait in question (G1 Finding 2). *(Tightened to "Glu180 absent"; see Deviations.)*

**Pools and quotas (37 in total; first written as 38, an arithmetic slip):**

| pool | lineage | quota | why |
|---|---|---|---|
| A | Clade IV | 10 | Most intact catalytic machinery (96 % Glu180); one characterized member (MphUPO), measured once |
| B1 | Ascomycetous 1st | 5 | No characterized member |
| B2 | Basidiomycetous 2nd | 5 | No characterized member |
| B3 | Family I (short), outside Clades I-V | 5 | No characterized member |
| B4 | Glomeromycete clade | 3 | No characterized member; small |
| B5 | Mucorales basal clade | 3 | No characterized member; small |
| C1 | Ascomycete long UPOs (Clade I) | 3 | Glu180 replaced by Ala/Ser/Thr; **contingent** on the structure check |
| C2 | Basidiomycetous 1st | 3 | Glu180 replaced by Asn/Ser; **contingent** on the structure check |

**Selection inside a pool:**

1. One protein per 70 % identity cluster (the G1 catalytic-core dereplication); proteins outside
   that set count as their own cluster.
2. At most one protein per genus.
3. Clusters taken largest first, round-robin across lineages (`LINEAGE_GROUP`), so no class fills a
   pool.
4. Within a cluster, prefer a protein from a RefSeq (`GCF_`) genome, then the highest LWR.

**Carried as information, not used to select:** best nr identity of the protein's 90 % cluster
representative, the T1 pocket match count (clade-V rule, not evidence outside Clade V; see
`MODEL_CARD.md`), and whether the protein is already in the skylark 96-protein panel
(`ASSAY_PANEL.md`, built to test the activity model rather than to cover lineages).

**What a result would mean.** Activity in a pool A or B candidate is a first measurement for a
lineage. Inactivity of a single candidate is weak evidence about its lineage (one gene model, one
expression attempt, one substrate). Pool C candidates should not be expressed until the Glu180
structural question is answered; they are listed so the structure check has concrete targets.

## Deviations from the rules, recorded when made

1. **Pool C must have Glu180 absent, not merely "not required"** (changed 2026-11-11 after the first
   run, job 31992843). As written, the rule let the selection prefer the 70 % catalytic-core
   clusters, which by construction carry Glu180, so all six pool C picks had Glu where their lineage
   almost never does (ascomycete long UPOs: 1.5 % Glu180). That defeats pool C's purpose, which is
   to give the structure check targets carrying the lineage trait. Rule 5 now reads: Cys39, Glu109,
   Asp111 present and **Glu180 absent**. Pools A and B are unaffected.

## Results (job 31993040, 2026-11-11)

Outputs: `candidates/out/shortlist.tsv` (one row per candidate, all columns), `shortlist.faa`
(sequences), `selection_summary.tsv` (exclusions per pool). **37 candidates in 8 pools; every pool met its
quota.**

### How many proteins each rule removed

| pool | in pool | long pendant | incomplete gene model | catalytic residues | LWR < 0.95 | length | eligible | picked |
|---|---|---|---|---|---|---|---|---|
| A | 8,437 | 0 | 1,043 | 909 | 6 | 465 | 6,014 | 10 |
| B1 | 6,228 | 6 | 548 | 2,476 | 48 | 207 | 2,943 | 5 |
| B2 | 4,394 | 38 | 303 | 1,811 | 37 | 176 | 2,029 | 5 |
| B3 | 6,291 | 199 | 321 | 2,268 | 5 | 349 | 3,149 | 5 |
| B4 | 784 | 585 | 108 | 54 | 4 | 4 | 29 | 3 |
| B5 | 315 | 0 | 292 | 1 | 0 | 4 | 18 | 3 |
| C1 | 33,783 | 342 | 3,670 | 2,357 | 7 | 2,430 | 24,977 | 3 |
| C2 | 6,378 | 14 | 395 | 1,673 | 116 | 406 | 3,774 | 3 |

### The candidates

| pool | # | protein | species | lineage | LWR | pendant | RefSeq genome | 70 % cluster | nr best id (90 % rep) | res. 180 | in skylark panel |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 1 | `F08DF36E_003893-T1` | *Fusarium proliferatum* | Sordariomycetes | 1.000 | 0.0783 | yes | 747 | 99.6 | E | no |
| A | 2 | `F391BB34_001299-T1` | *Neofusicoccum parvum* | Dothideomycetes | 1.000 | 0.1855 | yes | 727 | 99.2 | E | no |
| A | 3 | `F036E946_011951-T1` | *Aspergillus caelatus* | Eurotiomycetes | 1.000 | 0.0288 | yes | 483 | 90.2 | E | no |
| A | 4 | `F5362F83_010563-T1` | *Botrytis byssoidea* | Leotiomycetes | 1.000 | 0.0101 | yes | 138 | 100.0 | E | no |
| A | 5 | `F0C18633_012542-T1` | *Colletotrichum fioriniae* | Sordariomycetes | 1.000 | 0.0951 | yes | 350 | 100.0 | E | no |
| A | 6 | `F54975D0_011613-T1` | *Alternaria novae-zelandiae* | Dothideomycetes | 1.000 | 0.0134 | yes | 516 | 100.0 | E | no |
| A | 7 | `F1668EC8_010241-T1` | *Penicillium chrysogenum* | Eurotiomycetes | 1.000 | 0.1580 | yes | 136 | 100.0 | E | no |
| A | 8 | `F050C896_004813-T1` | *Pseudogymnoascus sp. VKM F-4519 FW-2642* | Leotiomycetes | 1.000 | 0.2152 | no | 32 | 100.0 | E | no |
| A | 9 | `FA681262_007606-T1` | *Cryphonectria parasitica* | Sordariomycetes | 1.000 | 0.2430 | yes | 126 | 87.5 | E | no |
| A | 10 | `F4A4F7CC_002961-T1` | *Dothidotthia symphoricarpi* | Dothideomycetes | 0.993 | 0.1786 | yes | 205 | 88.6 | E | no |
| B1 | 1 | `F0349401_004057-T1` | *Fusarium graminearum* | Sordariomycetes | 1.000 | 0.0116 | yes | 1295 | 99.3 | E | no |
| B1 | 2 | `F036E946_008369-T1` | *Aspergillus caelatus* | Eurotiomycetes | 1.000 | 0.0141 | yes | 558 | 100.0 | E | no |
| B1 | 3 | `F00BCD1C_008909-T1` | *Macroventuria anomochaeta* | Dothideomycetes | 1.000 | 0.2784 | yes | 94 | 100.0 | E | no |
| B1 | 4 | `F01E4952_000991-T1` | *Heterodermia albicans* | Lecanoromycetes | 1.000 | 0.2274 | no | 39 | 81.9 | E | no |
| B1 | 5 | `F515804B_006943-T1` | *Pezicula neosporulosa* | Leotiomycetes | 1.000 | 0.4479 | no | 3 | 59.1 | E | no |
| B2 | 1 | `F01F205C_007600-T1` | *Chalciporus chontae* | Agaricomycetes | 1.000 | 0.2099 | no | 289 | 73.2 | E | yes |
| B2 | 2 | `F020E2BE_002001-T1` | *Leccinum discolor* | Agaricomycetes | 1.000 | 0.2287 | no | 198 | 78.4 | E | no |
| B2 | 3 | `F0A381ED_010419-T1` | *Boletus coccyginus* | Agaricomycetes | 1.000 | 0.1742 | no | 121 | 76.5 | E | no |
| B2 | 4 | `F6EE9E7F_005487-T1` | *Psilocybe cubensis* | Agaricomycetes | 1.000 | 0.2377 | yes | 97 | 99.6 | E | no |
| B2 | 5 | `F082A955_007425-T1` | *Austroboletus rostrupii* | Agaricomycetes | 1.000 | 0.1780 | no | 76 | 79.7 | E | no |
| B3 | 1 | `F391BB34_011897-T1` | *Neofusicoccum parvum* | Dothideomycetes | 1.000 | 0.2541 | yes | 583 | 100.0 | E | no |
| B3 | 2 | `F34253C3_001499-T1` | *Colletotrichum siamense* | Sordariomycetes | 1.000 | 0.0501 | yes | 134 | 98.2 | E | no |
| B3 | 3 | `F15A4E68_003859-T1` | *Mycena rosella* | Agaricomycetes | 1.000 | 0.3215 | no | 74 | 100.0 | E | no |
| B3 | 4 | `F089B7A0_007745-T1` | *Hymenoscyphus koreanus* | Leotiomycetes | 1.000 | 0.4508 | no | 6 | 98.5 | E | yes |
| B3 | 5 | `F6EBB999_009246-T1` | *Leucodecton sp. Tripp_6156_NY-2796959* | Lecanoromycetes | 1.000 | 0.3171 | no | 4 | 73.5 | E | no |
| B4 | 1 | `FE9449CC_001707-T1` | *Rhizophagus prolifer* | Glomeromycota | 1.000 | 0.0384 | no | 32 | 100.0 | E | no |
| B4 | 2 | `F5BA5F3F_006438-T1` | *Basidiobolus sp. ARSEF 9900* | Basidiobolomycota | 1.000 | 0.9586 | no | 13 | 91.7 | E | no |
| B4 | 3 | `F2186B22_020194-T1` | *Dentiscutata heterogama* | Glomeromycota | 1.000 | 0.1087 | no | 14 | 100.0 | E | no |
| B5 | 1 | `F1707D45_012195-T1` | *Mucor velutinosus* | Mucoromycota | 1.000 | 0.0371 | yes | 41 | 92.6 | E | no |
| B5 | 2 | `F7B78768_009987-T1` | *Actinomucor elegans* | Mucoromycota | 1.000 | 0.1765 | no | 3 | 65.5 | E | no |
| B5 | 3 | `F07B100A_006870-T1` | *Absidia glauca* | Mucoromycota | 1.000 | 0.0359 | no | 2 | 100.0 | E | no |
| C1 | 1 | `F003EF48_000714-T1` | *Protomyces lactucae-debilis* | Taphrinomycetes | 1.000 | 0.6527 | yes | 1 | - | A | no |
| C1 | 2 | `F009D8E7_001268-T1` | *Neohortaea acidophila* | Dothideomycetes | 1.000 | 0.5821 | yes | 1 | - | A | no |
| C1 | 3 | `F00FD2C2_007881-T1` | *Blastomyces dermatitidis* | Eurotiomycetes | 1.000 | 0.3111 | yes | 1 | - | A | no |
| C2 | 1 | `F009D8E7_002898-T1` | *Neohortaea acidophila* | Dothideomycetes | 1.000 | 0.1774 | yes | 1 | - | N | no |
| C2 | 2 | `F0439E50_006510-T1` | *Suillus bovinus* | Agaricomycetes | 1.000 | 0.0948 | yes | 1 | - | N | no |
| C2 | 3 | `F0943129_016332-T1` | *Hyaloscypha bicolor* | Leotiomycetes | 1.000 | 0.1942 | yes | 1 | - | N | no |

### What the run showed, beyond the list

- **Mucorales UPO gene models in BFD are mostly incomplete.** 292 of 315 proteins in the Mucorales
  basal clade fail the gene-model check; 159 lack both start and stop codons. Clade IV, by contrast,
  is 88 % complete. Pool B5's 3 picks come from 18 eligible proteins. This is a BFD annotation
  problem for Mucorales and bears on HANDOFF open item 4 (gene-model QC): catalytic-residue
  statistics for Mucorales should not be read as biology until it is resolved.
- **The glomeromycete clade attracts divergent outsiders.** Of its 784 BFD members only 121 are
  Glomeromycota (360 are Dothideomycetes), and 585 carry `long_pendant`. Most are probably
  divergent proteins drawn to a small, isolated reference clade, not glomeromycete UPOs; read this
  clade's placement calls only together with the flag. The eligibility rules excluded them (29
  eligible, all short-pendant). Pick B4-2 is a *Basidiobolus* (Basidiobolomycota) with pendant 0.96,
  just under the flag; treat it as the weakest B pick.
- **Basidiomycetous 1st holds about 2,000 ascomycete proteins** (Dothideomycetes 1,364,
  Eurotiomycetes 680 of 6,378), and unlike the glomeromycete case they are *not* long-pendant (14 of
  6,378 flagged). Zamocky noted one ascomycete transfer into this clade; this is two orders of
  magnitude more. It may be transfer, or an ascomycete lineage the reference set does not resolve.
  Not interpreted here; it needs its own look before any claim.
- **Pool C, Basidiomycetous 1st:** two of three picks are ascomycetes (*Neohortaea*, *Hyaloscypha*)
  carrying Asn at 180, consistent with the horizontal transfer Zamocky noted into this clade.
  Pool C proteins fall outside the catalytic-core 70 % clusters (they lack Glu180), so cluster size
  does not rank them; after the RefSeq preference, the order is by protein ID, i.e. arbitrary.
- **Sequence novelty is low, as expected:** most picks have a 90 % representative at 74-100 %
  identity to an nr protein, and pool A picks come almost entirely from RefSeq genomes whose
  proteins are in nr. These are first measurements for their *lineage*, not unseen sequences.
- **Overlap with the skylark 96-protein panel: 2** (B2-1 *Chalciporus chontae*, B3-4 *Hymenoscyphus
  koreanus*). The two lists answer different questions and can be assayed together.

### Suggested first expression set (12)

Pool A ranks 1-4 (four classes of Clade IV), the top pick of each B pool (B1-B5), and B1-2, B2-2,
B3-2 for a second member in the three largest uncharacterized lineages. Pool C waits for the
Glu180 structure check.
