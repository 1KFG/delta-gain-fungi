# G1 findings: fungal UPO diversity across BFD

**Measured 2026-11-08.** Numbers come from `G1_BFD_SUMMARY.md`, auto-generated from
`data/bfd_full/out` by `upo_sites/bin/summarize_bfd.py`. Per-protein table:
`data/bfd_full/upo_sites_bfd_all.tsv.zst`.

## The census

**97,520 PF01328 proteins in 15,137 genomes.** BFD holds 23,451 genomes with gene models, so
**8,314 fungal genomes carry no PF01328 at all.** Maximum copy number in one genome is **144**
(*Sphaerobolus stellatus*, Agaricomycetes), followed by seven *Mycena* species at 67-106.

## Finding 1: only about a third of fungal PF01328 is UPO-like

PF01328 is a broad heme-thiolate peroxidase family, not a UPO-specific one. Using the catalytic
residues that the 2025 panel shows are present in all 19 characterized UPOs — Cys39, Glu109,
Asp111, Glu180 — **34,561 of 97,520 proteins (35.4 %) carry the full set.**

The discriminating residue is **Glu180**, which the paper reports as universal:

| Set | Glu180 |
|---|---|
| Calls inside the validated identity range | **93.7 %** |
| Extrapolated calls | 39.3 % |
| Unassigned | 21.8 % |

**93.7 % in the validated set reproduces the published claim independently.** The low rate
elsewhere is therefore informative about those proteins, not about the method.

**Consequence: the catalytic core set of 34,561 is the defensible UPO population.** Every
downstream analysis — novelty, regulation, candidate selection — uses that set, not the raw
97,520.

## Finding 2: the Clade I bin at low identity is not trustworthy

Clade I calls carry Glu180 at only **16.8 %**, yet the same alignment gives all five characterized
Clade I enzymes a Glu at position 180. This is not an alignment artifact for long UPOs.

The most likely reading: at ~34 % median identity, the long characterized enzymes are simply the
nearest available match for divergent PF01328 proteins that are **not UPOs**, so the Clade I bin
absorbs them. **Unfiltered Clade I counts (33,424) are an upper bound.** After the catalytic
filter, Clade I falls to **5,101**.

## Finding 3 — the headline: the best-supported fungal UPO clade is represented by one uncharacterized enzyme

**Clade IV is unlike every other clade in this dataset.**

| Clade | n | median identity to nearest reference | share >= 52.7 % |
|---|---|---|---|
| **IV** | **6,922** | **66.5** | **71.8 %** |
| II | 3,113 | 39.1 | 12.9 % |
| V | 24,847 | 41.3 | 10.6 % |
| III | 10,154 | 33.5 | 10.6 % |
| I | 33,424 | 34.2 | 11.2 % |
| unassigned | 19,060 | 27.3 | 0 % |

Clade IV sits at **66.5 % median identity** where every other clade sits at 33-41 %. It is the only
clade where the majority of BFD assignments fall **inside** the range in which the caller was
actually validated. It is also catalytically intact: **6,002 of 6,922 pass the Glu180 filter
(86.7 %)**, the highest of any clade.

It spans Sordariomycetes (2,531), Dothideomycetes (2,464), Eurotiomycetes (1,485) and
Leotiomycetes (364).

**And its entire reference is MphUPO** — one protein, from *Macrophomina phaseolina*, listed in
Table S1 of the 2025 study as "**Putative UPO**", the single uncharacterized member of the panel.
Its measured peroxygenase activity on 4-propylguaiacol is 0-5 %, from n = 1.

So the most confidently delimited, most catalytically intact, and second-largest UPO group in
fungi is anchored on an enzyme nobody has characterized. **This is the clearest answer the project
has to the original goal of finding measurable novelty in UPO subtype lineages, and it is a
concrete experimental recommendation: characterize a Clade IV enzyme.**

### Caveats on Finding 3, stated plainly

- **Every Clade IV call is flagged `provisional_cladeIV` by the tool**, because leave-one-out
  cannot validate a clade with one member — removing MphUPO leaves no Clade IV reference at all.
  High identity to MphUPO is strong evidence of a coherent group; it is not evidence that the
  group's boundary is correctly drawn against clades not in the reference set.
- The tightness could partly reflect sampling: MphUPO is a Dothideomycete, and BFD samples
  Dothideomycetes and Sordariomycetes densely. Median 66.5 % identity is high regardless, but a
  denser reference set could split Clade IV.
- MphUPO's own low T1 (0-5 %) is a single measurement on one enzyme, and T1 is 4-propylguaiacol
  peroxygenation, not a general activity statement. **It does not license "Clade IV is inactive".**

## Finding 4: clade confidence across BFD

**Corrected 2026-11-08.** The first version of this section mixed two different definitions of
"validated range" — one table matched the confidence string exactly (which silently excluded every
Clade IV call, because those carry a `+provisional_cladeIV` suffix) and another used a prefix
match (which included them). The summariser now uses one definition everywhere. Corrected numbers:

| Set | n | % of 97,520 |
|---|---|---|
| Clade calls inside the validated identity range | **12,818** | 13.1 % |
| of which Clade IV (one characterized member, cannot be validated) | **4,968** | 38.8 % of the validated set |
| Validated **and** anchored on a clade with more than one reference | **7,850** | **8.1 %** |

**The honest state of fungal UPO subtype annotation**: with 19 characterized enzymes as the only
anchor, roughly 8 % of 97,520 fungal proteins have a clade call that is both inside the validated
identity range and anchored on a clade with more than one reference enzyme. A further 5 % sit at
high identity to MphUPO alone.

Improving this needs either more characterized enzymes or phylogenetic placement with proper
uncertainty, which is `PLAN.md` section 3.1 item 2.

## What this changes for the *Pichia* expression plan

Decision D1 asked whether Eurotiales sourcing was firm. It was previously moot. It is now
**actionable in a new way**: Clade IV contains **1,485 Eurotiomycetes proteins**, most at high
identity to MphUPO and with intact catalytic machinery. If a Eurotiales source is still wanted,
Clade IV Eurotiomycetes members are concrete, well-supported candidates — and expressing one would
simultaneously characterize the clade that most needs it.

This is a candidate-selection rationale based on **annotation confidence and catalytic integrity**,
not on a predicted activity score. No such score exists; see `RESULT_T2_negative.md`.

## Still open in G1

- **Novelty against public sequence space** (step 4). Dereplication of the 34,561 catalytic-core
  proteins is running; the DIAMOND search against `nr` follows. The deposition confound stands:
  absence from `nr` can mean "not deposited", so a nucleotide/WGS comparison is required before any
  novelty claim.
- **Contamination checks** for lichen and Glomeromycota lineages, untested, required before any
  lineage-level claim is reported.

---

## Finding 5: novelty against public protein space — and why the number is mostly about deposition

**Measured 2026-11-08.** DIAMOND blastp of the 7,055 catalytic-core representatives (90 %
dereplication) against NCBI `nr` (523 GB, Jan 2026 build). Job 31865316, 42 minutes.
Raw hits: `data/bfd_full/novelty/reps090_vs_nr.tsv.zst`.

**Every one of the 7,055 representatives has a hit in `nr`.** None is unprecedented.

| Best-hit identity | Representatives |
|---|---|
| = 100 % | **2,026** |
| >= 95 % | 2,781 |
| >= 90 % | 3,249 |
| >= 70 % | 5,788 |
| < 70 % ("novel" by this cut) | **1,267 (18.0 %)** |
| < 50 % | 160 |

Median best-hit identity is **87.3 %**.

### The confound is visible in the data, not just a caveat

**For 4,413 of 7,055 representatives (62.6 %), the best `nr` hit is to the same genus as the query
genome.** Those are most plausibly the BFD assembly itself, or a close relative, already deposited.
2,026 representatives match something in `nr` at exactly 100 % identity.

So this analysis measures **deposition status at least as much as biological novelty**, exactly as
flagged in `PLAN.md` before it was run. The lineage breakdown makes the point sharply:

| Lineage | novel (< 70 %) | n |
|---|---|---|
| Tremellomycetes | **54.0 %** | 87 |
| Lecanoromycetes (lichens) | **33.7 %** | 264 |
| Leotiomycetes | 27.6 % | 293 |
| Sordariomycetes | 16.8 % | 1,373 |
| Agaricomycetes | 16.6 % | 3,135 |
| Dothideomycetes | 15.1 % | 1,021 |
| **Eurotiomycetes** | **10.8 %** | 343 |
| Chytridiomycota | 6.2 % | 96 |

**Eurotiomycetes is the least novel lineage and Lecanoromycetes among the most.** That ordering
tracks how heavily each group has been sequenced and deposited — *Aspergillus* and *Penicillium*
are among the most deposited fungi in `nr`, lichens among the least. It is not evidence that
lichen UPOs are biologically more unusual.

**A nucleotide/WGS comparison is still required before any novelty claim is made.** It is not done.
Until it is, the defensible statement is: *18 % of dereplicated fungal UPOs have no close protein
relative in `nr`, and that figure is substantially driven by which genomes have been deposited.*

### Finding 5a: novelty does not track clade

| Clade | novel (< 70 %) | n |
|---|---|---|
| II | 19.8 % | 585 |
| V | 18.8 % | 2,466 |
| III | 17.2 % | 1,079 |
| I | 16.8 % | 1,851 |
| IV | 14.4 % | 789 |
| unassigned | 26.7 % | 285 |

**Range 14-20 % across all five clades — essentially flat.** No clade is meaningfully more novel
than another. This is a negative result and it is worth stating: the subtype a UPO belongs to does
not predict whether it has a close public relative.

It also tempers Finding 3. Clade IV is the *least* novel clade by this measure (14.4 %), which fits
its high internal identity: its members have close relatives, those relatives are simply other
uncharacterized Clade IV proteins. **Clade IV is under-characterized, not under-sequenced.** That
strengthens rather than weakens the case for characterizing one — the sequences are already public
and available; what is missing is biochemistry.
