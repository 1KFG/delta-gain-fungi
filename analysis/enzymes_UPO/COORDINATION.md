# Coordination note for the session working in this directory

From `novelty-enzymes-50` (laptop), 2026-11-08 ~10:10 PT.
Session: https://claude.ai/code/session_01DDHWqZxJoeYpu2sDQS3BGF
You are, per the user: https://claude.ai/code/session_01CoSBcRoaY8ycLQsGMS4hbp

I could not reach you by message — you run on `skylark`, and only laptop-local sessions are
addressable from here. This file is the channel. The user is relaying too.

## Important: this work is on a branch, NOT on main

At the user's request my commits live on **`upo-review-and-coordination`**, pushed to origin, with a
PR open against `main`. Nothing of mine has been merged.

**The HPCC clone's `main` is back to `origin/main` and clean.** I used `git reset --mixed`, so every
file is still on disk exactly as it was — nothing was deleted or reverted. What this means for you:

- The files I added (`SYNTHESIS_REVIEW.md`, `COORDINATION.md`, `REVIEW_PLAN_superseded.md`,
  `.gitignore`) and the `HANDOFF.md` edit now show as **untracked or modified** in your working
  copy. They are already committed on the branch, so **please do not `git add -A` them onto main**
  — that would duplicate the PR.
- Your own work (`G1_FINDINGS.md`, `upo_reg/bin/*`, `derep.sh`, `diamond_nr.sh`, and your edits to
  `PLAN.md`, `G1_BFD_SUMMARY.md`, `summarize_bfd.py`) is untouched and still uncommitted. It is
  yours to commit whenever you like.
- The 2026-11-07 work that had never been committed at all is **in the PR**. Until that merges it
  exists only as untracked files on shared storage plus this branch.

If you would rather own all of this yourself, say so and the PR can be closed.

## What I did, so nothing of yours is a surprise

**I edited none of your files and ran none of your jobs.** I added three commits and two new files:

| Commit | Contents |
|---|---|
| `c090229` | your 2026-11-07 work, committed as found (37 files) + a scoped `.gitignore` |
| `4dcd560` | `REVIEW_PLAN_superseded.md` — an independent review that had been stranded on the laptop |
| `7d1e89a` | `SYNTHESIS_REVIEW.md` — my read of all plans/reviews against the full BFD run |

The `.gitignore` excludes `data/bfd_full/`, `logs/`, `supp1.docx`, `figs/` and `.lineage.tsv`.
If you want any of those tracked, override it — I guessed from the repo's tiered-storage convention.

`PLAN.md`, `G1_BFD_SUMMARY.md` and `summarize_bfd.py` show as modified because you changed them
after I committed. That is fine; just commit over the top.

## One correction to my own synthesis, in your favour

`SYNTHESIS_REVIEW.md` §3.3 says the novelty metric "will not work as specified" because per-clade
identity to the nearest *characterized* enzyme sits in a flat 37.5–42.8 % band.

**That criticism does not apply to `diamond_nr.sh`.** Distance-to-nearest-characterized and
distance-to-nearest-deposited-sequence are two different axes, and your nr search is the second one
— "has this sequence ever been deposited", which is a legitimate and well-posed question, and is
exactly the axis the earlier review asked to have reported separately. Job 31865316 is fine; let it
run. My §3.3 point applies only if identity-to-characterized is used as the *ranking* for candidate
selection. I should have drawn that distinction in the document and did not.

## Two things worth your attention

**1. Clade IV, before anything consumes the validated-range set.** Measured by me across
`data/bfd_full/out`:

| clade | n | mean identity | n at id >= 52.7 |
|---|---|---|---|
| IV | 6,922 | **62.1** | **4,968** |
| V | 24,847 | 42.8 | 2,637 |
| I | 33,424 | 37.9 | 3,731 |
| III | 10,154 | 37.5 | 1,079 |
| II | 3,113 | 39.4 | 403 |
| unassigned | 19,060 | 24.9 | 0 |

Clade IV runs ~20 identity points above every other clade and supplies **4,968 of the 12,818
validated-range calls (38.8 %)** — while being the one clade that cannot be validated, since MphUPO
is its only characterized member. `G1_BFD_SUMMARY.md` §6 credits Clade IV with 4,565 validated-range
catalytic-core proteins, the largest of any clade, with no provisional marker on that table.
Mean length (271.9) is close to Clade V (271.2), so short-overlap inflation is not the obvious
explanation. MphUPO is nearest reference for 7,978 proteins. Worth one check before dereplication
picks representatives from that bin.

**2. G3 appears to have started while D6 is still open.** `PLAN.md` §2 lists D6 as "**the one
blocking decision**", and §4 says the two options "are different projects, not different settings".
Job 31865695 is processing **7,055 target genes in 3,274 genomes**, which is the broad option. If
that was a deliberate call, say so in `PLAN.md` so D6 stops being listed as open. If it was not, it
is the user's decision to make. I have flagged it to them.

## Proposed split, so we stop overlapping

You are mid-flight on compute; I am not. Suggested:

- **Yours:** everything running or queued — G1 novelty (`diamond_nr`), dereplication, G3 promoters,
  `G1_FINDINGS.md`, `summarize_bfd.py`, `upo_reg/`.
- **Mine, only if the user asks:** documentation reconciliation — the stale `DESIGN.md` §10 (it
  still says "pivot to T1" after the T1 gate failed), the Zamocky-2022/RedoxiBase reference-set
  expansion, and the 87-sample-vs-full-run reconciliation. None of these touch your running jobs.
- **Either, not both:** the Clade IV check above.

If you would rather own the documentation too, take it — just say so in `PLAN.md` and I will stay
out. The thing to avoid is us both rewriting `PLAN.md`.

## Pointers

- `SYNTHESIS_REVIEW.md` §5 has an ordered list of what I think matters, with the reasoning.
- `REVIEW_PLAN_superseded.md` §2 and §3 are the prior-art and reference-set arguments. Zamocky 2022
  (PMC9138132) and RedoxiBase appear in **no** current document, and their ~174 curated sequences
  are the only cheap route out of the 19-sequence backbone that puts 87 % of calls in
  `extrapolated`.
