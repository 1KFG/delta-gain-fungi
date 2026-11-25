---
title: "Owner decisions for the melanin and mycosporine census (2026-11-23 and 2026-11-24)"
---

::: {.callout-note}
## Status and limits
- Source: the owner's messages in the session `dgf-pathways-pigments` (id `80c2c1e5-ee6d-42e4-b024-bcd02b7e3b1b`). Quotes are verbatim.
- Every result that follows from these decisions stays `.PROVISIONAL`.
- Open items are listed at the end. They are not decided.
:::

## Decisions

| Date | Item | Owner's words | Applied as |
|---|---|---|---|
| 2026-11-23 | Step 2 decisions 3 to 9 | "accept 3-9 as listed" | PKS rule values kept, joint PKS tree; t4hn/t3hn rule values kept, one SDR tree; tyrosinase 63.6 kept; no lowering; mixed mycosporine references kept; Q9Y7A7 kept; decoys added for ayg1, scytalone_dehydratase, hmga |
| 2026-11-23 | omt floor | "use 100 but we may need to investigate again in a second round" | floor override 100.6 (highest decoy score at the time) |
| 2026-11-23 | Decoys for hppd and aspulvinone_e_synthase | "add decoys if you can identify them" | 28 decoys fetched from UniProt for five steps |
| 2026-11-24 | Step 3, then Q8NKC1 | "step 3 first then add Q8NKC1 as omt decoy" | timing job on the 100.6 table; Q8NKC1 added after (new floor would be 120.7) |
| 2026-11-24 | Full search floor for omt | "accept your reccomendation" | full search run at 100.6; rows between 100.6 and 120.7 kept |
| 2026-11-24 | Step 5 order | "step 5 with PKS tree ... fix tables makes sense if they are needed" | PKS joint tree; genome prefix fix in our copy only (alias file) |
| 2026-11-24 | heptaketide_pks and wA | "can report wA as unresolved but need more details on how to rescure this information as it is still a valid gene to include" | wA reported unresolved; options note written; the choice is NOT yet made |
| 2026-11-24 | PKS outgroup | "leave it as no outgroup with note to return to this" | no outgroup; separate investigation later |
| 2026-11-24 | PKS rerun | "yes do this with identity 0.8 and bootstrap of the tree if we believe we have a good sampling" | rerun at identity 0.8 with bootstrap |
| 2026-11-24 | scytalone_dehydratase | "drop cfoI leaf and rerun. can you regen the tree?" | leaf A0A2I2F2K6 dropped after seeing the tree; tree regenerated |
| 2026-11-24 | omt | "try only Phaffia reference alone" | single-reference placement, Nostoc MysB as a named leaf |
| 2026-11-24 | atp_grasp_ligase | "ok to not be rooted if that is what this means" | stays a search-hit step (no alignable non-target leaf exists) |
| 2026-11-24 | hmga | "keep as is" | stays as placed ('family present' reading; the placement carries no information) |
| 2026-11-24 | Target rule, 2 or more references | "accept largest side-rule or steps with 2 or more references" | largest-side rule accepted as PROVISIONAL for steps with 2 or more references (the one-reference rule was accepted on 2026-11-23) |
| 2026-11-24 | tyrosinase and hppd | "yes address this with high rep cap" | rerun with a high representative cap |

## Open (not decided)
1. wA and heptaketide_pks: per-step decoys in the joint tree, or the union rule; and step defined by reaction or by pathway role. See `2026-11-24-wA-heptaketide-options.md`.
2. BFD data fix for the five re-annotated genomes (the BFD owner's call). See `2026-11-24-genome-prefix-mismatch.md`.
3. Identity rule for clustering (lowest identity giving at least 100 clusters) is a choice of the step-5 agent, not an owner decision.
4. Merge of PR #22 (Jason).
