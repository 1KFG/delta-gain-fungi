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

**Pools and quotas (38 in total):**

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

## Results

*(appended after the run)*
