# Fungal pigments: carotenoid, melanin and mycosporine pathways in the fungal genome census

Status: pre-census chapter. Written 2026-11-21. Code and tables: `analysis/pathway_census/`.

::: {.callout-note}
## What this chapter is, what it is not, and where it stops being reliable

This chapter describes the **curated knowledge bases** for four fungal pigment and pigment-like pathways: carotenoid, DHN melanin, L-DOPA melanin and mycosporine. It describes what each table holds and how strong the evidence is for each step.

**It reports no gene counts.** The census has not been run on the BFD genomes for any of the four pathways. It waits for the owner's review of the curated step sets (section 6). Every statement of the form "fungi have X genes" is therefore missing on purpose.

**Nothing here was tested at the bench.** The curators read PubMed abstracts and database records. They did not read full texts, except where a table says so.

**Two evidence readings exist.** Under the strict quote rule used in the tables, no carotenoid route and no complete DHN route is above tier T3. Under the alternative context reading, two links (carG, Q9P885; SDH1, P56221) stay at T1a and the route tiers rise. The owner has not chosen between them.

Labels: MEASURED = counted from the committed tables on 2026-11-21. REPORTED = taken from a message of the EBC agent, not re-measured. UNVERIFIED = not checked. General biology that is not in the tables is marked *background* and was not verified here.
:::

## 1. Why these four pathways

The project began as two chapters for the fungal genome census, carotenoid and melanin biosynthesis, modelled on the earlier peroxygenase (UPO) chapter. Structure and affinity prediction was dropped. Mycosporines and mycosporine-like amino acids (MAAs) were added later because the metabolite corpus (held by the EBC agent) shows MAA library names in marine samples.

The working question was whether a genome that carries a pathway shows its product more often in the metabolite corpus. That comparison failed its pre-set detectability gate and was not run (section 5). The curated tables remain useful on their own as a reference for later census work.

## 2. How the tables are built

Each pathway has five flat tables in `kb/<pathway>/`: `compound`, `reaction`, `reaction_participant`, `step`, `step_reaction`. A step is a position in the pathway, defined by a protein family (Pfam, TIGRFAM or KO). A `step_reaction` row links a step to a reaction through one named protein of one named organism, with an evidence tier.

| Tier | Meaning |
|---|---|
| T1a | The fungal enzyme was assayed, or its product was identified after heterologous expression. The quoted abstract sentence must name that protein and state an assay or product |
| T1b | Gene deletion or complementation with analysis of intermediates in the fungus |
| T1c | Deletion with a colour phenotype only |
| T2 | Reaction shown for the same enzyme family in another organism |
| T3 | Inferred from family or orthology, or a T1 claim whose quote fails the rule above |
| T4 | Hypothesis with no direct evidence |

A route tier that holds only across several organisms is not reported as a tier for any single organism.

An independent script (`bin/verify_kb.py`) re-fetches each accession, gene name, PubMed identifier, quoted sentence, sequence hash, ChEBI or PubChem identifier and formula, and checks them against the files. MEASURED result: 1,677 checks, 1,647 pass, 6 fail, 24 unchecked (`kb/VERIFICATION_REPORT.md`). A pass means the public record agrees with the file. It does not mean the biology is right.

## 3. The four pathways

### 3.1 Carotenoids

**Steps in the table (MEASURED):** a core of four required steps (GGPP synthase, phytoene synthase, lycopene cyclase, phytoene desaturase) and eleven optional steps for branches: retinal (carX), neurosporaxanthin (carT, ylo1), astaxanthin (P450 synthase, CrtR reductase), trisporoids (CarS, AcaA, Tsp3, Tsp2, Tdh) and a torularhodin step with no enzyme identified.

**Tables:** 21 compounds, 23 reactions, 49 links.

| Tier | T1a | T1b | T1c | T3 | T4 |
|---|---|---|---|---|---|
| Links | 15 | 9 | 3 | 21 | 1 |

**What holds up:** individual enzyme links exist at T1a and T1b for several fungal species (MEASURED, T1a or T1b rows: *Neurospora crassa*, *Fusarium fujikuroi*, *Mucor circinelloides*, *Phaffia rhodozyma*, *Phycomyces blakesleeanus* and *Ustilago maydis*). **What does not:** under the strict rule every route that ends in a named product is T3, because the link carG (Q9P885) fails the rule that the quoted sentence names the protein, and it decides the carotenoid route tiers (curation ledger; I did not re-derive this per route).

**Evidence review changed tiers.** Phytoene desaturase tiers in the first review list differed from the abstracts by organism: four-step T1a for three organisms, five-step T1a for *Neurospora* AL-1 only, T1b for *Fusarium* carB. The owner has not yet decided these changes.

**Torularhodin** (the red yeast carotenoid) has no identified enzyme in the table. A census cannot place it.

### 3.2 DHN melanin

**Steps (MEASURED):** three alternative polyketide synthase steps (THN, hexaketide, heptaketide), a hydrolase (Ayg1), two naphthol reductases (4HNR and 3HNR), scytalone dehydratase and two polymerisation steps (an Abr2-type laccase and an Abr1-type multicopper oxidase).

**Tables:** 10 compounds, 11 reactions, 32 links.

| Tier | T1a | T1b | T1c | T3 |
|---|---|---|---|---|
| Links | 11 | 5 | 6 | 10 |

**Review found a swap.** In the pathway review the two DHN reductases were swapped. Q9HFV6 is the 4HNR and Q12634 the 3HNR. The curator caught it by fetching the UniProt records. The accession of every reference is now fetched, and a recalled accession counts only as a lead.

**Open items:** Q9HFV6 has no UniProt gene name, so the gene-name rule cannot be met; the vermelone dehydration link on P56221 (SDH1) is T3 under the strict reading; two routes have no product compound; the *Colletotrichum* and *Glomerella* organism key needs an owner decision.

### 3.3 L-DOPA and other melanins

**Steps (MEASURED):** tyrosinase, an exogenous-substrate laccase, the tyrosine-degradation enzymes HPPD and HmgA (the pyomelanin branch), an Asp-melanin pair (aspulvinone E synthase and a short tyrosinase-like protein), and a "dopachrome onward" step with no fungal enzyme identified.

**Tables:** 11 compounds, 13 reactions, 29 links.

| Tier | T1a | T1b | T1c | T3 | T4 |
|---|---|---|---|---|---|
| Links | 2 | 3 | 2 | 20 | 2 |

This is the weakest of the carotenoid and melanin tables: 20 of 29 links are T3. Two T1a links remain after the strict rule.

**Reference sequence error.** Q2UP46 was carried as the assayed melB protein. The assayed sequence is A0A1S9DK56 (615 of 616 positions agree with PDB 3W6W; Q2UP46 agrees at 585). The reference was moved. A0A1S9DK56 has an ORF name only, so it is tied to the paper by structure, not by gene name. This tie is an owner decision.

### 3.4 Mycosporines and MAAs

**Steps (MEASURED):** DDGS (a sugar-phosphate cyclase of the 3-dehydroquinate synthase fold), an O-methyltransferase, an ATP-grasp ligase, and two steps with no fungal enzyme identified (a glucosyltransferase for the glutaminol glucoside and glutamicol formation).

**Tables:** 9 compounds, 8 reactions, 10 links.

| Tier | T1a | T1c | T2 | T3 | T4 |
|---|---|---|---|---|---|
| Links | 2 | 1 | 2 | 3 | 2 |

**All fungal references come from one species.** Literature curation found fungal MAA biosynthetic genes in *Phaffia rhodozyma* only (PMID 36983003: a three-gene cluster, DDGS, OMT and an ATP-grasp ligase, GenBank OQ547787). BFD holds 12 *Phaffia* genomes (MEASURED). The two T2 links use the cyanobacterium *Nostoc punctiforme* as the reference organism.

**Six of eight references are `verified: false`.** The three *Phaffia* proteins (WEV88679.1, WEV88680.1, WEV88681.1) are GenBank records with no UniProtKB entry, so the gene-name rule cannot be met. The records cite the paper by DOI. Accepting that tie is an owner decision. The *Nostoc* references cite a genome paper, not a paper on these steps, so they were tied by gene name in another abstract.

**No kingdom-wide fungal gene list is claimed.**

### 3.5 Candidate search for mycosporine genes in other fungi

An earlier project folder (`/bigdata/stajichlab/shared/projects/Mycosporine`, an HGT candidate search) holds a 5,000-genome protein database. A DIAMOND search of six *Phaffia* and older query proteins against it was started on 2026-11-21. It **hit the 2-hour time limit and was cancelled**. `diamond_hits.tsv` is empty. The tree-placement job that depended on it never started. So no candidate list exists. Earlier, and from the owner's note: no horizontal gene transfer has been identified yet; the earlier work was about finding candidates first.

The search needs a smaller query set or a longer time limit (PROPOSED).

## 4. What the tables can be used for now

- As a reference of curated, quote-checked enzyme links with the tier attached. A T3 link is an annotation. It is not an assay.
- As the definition of steps for a census: protein family per step, optional steps, alternative routes, and decoy proteins (121 decoy references pass the existence and sequence checks).
- As a list of known gaps. Torularhodin, dopachrome onward, and the two mycosporine steps without an enzyme cannot be searched for. The census can report only the steps that have a defined family.

## 5. What happened when the corpus was compared with the genomes

Counts in this section are REPORTED by the EBC agent (corpus version `quant_v5`) unless marked otherwise. The gate was fixed before counts arrived: at least 6 capable accessions, at least 3 of them showing the compound, and a positive-control cohort at 3 accessions. The thresholds were not lowered afterward.

| Target | Fungal-reference accessions showing it | Gate |
|---|---|---|
| β-carotene, zeaxanthin, canthaxanthin | 0 of 84 each | Fail |
| Astaxanthin, retinal, torulene | 1, 2 and 0 of 84 | Fail |
| Melanin library names | 5 of 7 names zero; L-DOPA 0 fungal groups | Fail |
| Melanin formulas | 10 to 35% of groups in every cohort | Background, not evidence |
| MAAs in fungi | 0 of 82 | Fail |
| MAAs in marine samples | 36 of 87 marine accessions against 216 of 5,379 others; odds ratio 16.9 | Positive control only |

A weaker result than the gate is "not testable", not "no effect". The marine MAA signal shows that MAAs are present in marine samples. It does not show who made them; cyanobacteria, algae, fungi and animals that take MAAs up all occur there (UNVERIFIED).

Reasons a null cannot be read as absence (none resolved): the compound is not extracted or ionised by the corpus methods; the cultures do not express the pathway; melanin is polymeric, so the monomer is absent; the library lacks the names; the fungal reference set is small and dominated by yeasts that do not make these compounds. Full detail is in [the pathway and corpus report](pathway_corpus.html).

## 6. What is needed before a census

These are the owner's decisions. The full lists are in `refs/<pathway>/UNVERIFIED.md` and the review packet.

1. **Strict or context reading of the quote rule.** The strict reading is in the tables. Under the context reading, carG and SDH1 stay T1a and the carotenoid and DHN route tiers rise.
2. **Entries with no UniProt gene name** (Q9HFV6, Q2TL94, and ORF-name-only entries such as A0A1S9DK56): write a rule or drop them.
3. **GenBank-only *Phaffia* references:** accept a record that cites the paper as the tie, or not.
4. **Open per-pathway items** listed in section 3.

After these, the census can run: count genes per step per genome, report route completeness per organism, and keep tiers attached to every claim. I propose a start on the *Phaffia* genomes (a pathway with a known fungal reference) and on *Aspergillus* melanin (DHN) genes, since both have verified T1 references. This is a proposal and has not been tested.

## 7. Not done here

- No census on real data for any pathway.
- No bench work, no chemical identification, no full-text reading of the cited papers.
- No corpus comparison for fungal carotenoids, melanin or MAAs (gate failed).
- No fungal MAA gene list beyond *Phaffia rhodozyma*.
- No homology-based candidate search for mycosporine genes in other fungi (the job was cancelled).
- No antiSMASH cluster call for these pathways. The antiSMASH work for the EBC pilot (see the BGC subset note) does not target pigment clusters.
