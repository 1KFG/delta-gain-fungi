# Heme access channel width: can Karich's ring-size rule be transferred?

**Measured 2026-11-08.** Tool: `upo_sites/bin/heme_channel.py`. Structures in
`data/structures/ref/`. This returns to the original goal 2 question — which rings a candidate can
open — which the T1 pivot had drifted away from.

## The rule we are trying to apply

Karich et al. 2017 explain PAH substrate range mechanically rather than statistically: **AaeUPO has
a ~7 Å heme access channel** and converts PAHs up to about 6 Å minimum width, while **MroUPO has a
~11 Å channel** and handles bulkier ones. That is a ring-size rule. Applying it across BFD needs a
channel width per enzyme.

## The measurement

CAVER's download URL had been guessed wrongly and 404'd, so this was implemented directly: protein
heavy atoms plus the heme as hard van der Waals spheres, a regular grid, an exact Euclidean
distance transform for clearance, and the bottleneck radius found by binary search on connected
components — exact to grid resolution and much faster than a path search. Width is reported as
diameter, to match how Karich reports it.

One implementation trap worth recording: the iron sits under the porphyrin, so the nearest free
voxel to it is a sub-ångström crevice between heme atoms. Seeding there made every structure report
one grid unit. The seed must be the **widest free point within ~9 Å of the iron** — the distal
pocket the channel actually serves.

## Result: the bottleneck does not reproduce the published contrast

| Structure | measured bottleneck | Karich published |
|---|---|---|
| AaeUPO (2YOR) | **6.00 Å** | ~7 Å |
| MroUPO (5FUJ) | **6.40 Å** | **~11 Å** |
| MroUPO (7ZBP) | 5.91 Å | ~11 Å |
| HspUPO (7O1R) | 6.98 Å | — |
| CviUPO (7ZCL) | 4.99 Å | — |

AaeUPO lands within 1 Å. **MroUPO is off by 4.6 Å, and the AaeUPO-versus-MroUPO contrast is the
whole basis of the rule.**

## But the profile shows the ordering is right, only the magnitude is wrong

A single bottleneck number hides a wide chamber behind a narrow neck, so the corridor was profiled:
the widest clearance at each distance shell from the iron.

| Distance from Fe | AaeUPO | MroUPO |
|---|---|---|
| 1 Å | 6.5 Å | 6.7 Å |
| 3 Å | 7.0 Å | 7.7 Å |
| 5 Å | 8.6 Å | **10.0 Å** |
| 7 Å | 10.8 Å | **12.9 Å** |

**MroUPO is wider than AaeUPO at every shell**, by 0.2 Å at the iron widening to about 2 Å at 5–7 Å
out. So the *qualitative ordering* Karich relies on is reproduced; the *magnitude* (4 Å) is not.

Beyond roughly 10 Å the profile grows without limit in both enzymes — that is the corridor leaving
the protein into bulk solvent, not a chamber, and it carries no information.

Notably, AaeUPO reads 7.0 Å at 3 Å from the iron and MroUPO reads 10.0 Å at 5 Å, which is close to
the published 7 and 11. That is **suggestive that the published "width" is measured near the active
site rather than as a global bottleneck**, but the distances are not matched and this is a
post-hoc observation on two enzymes. It is a hypothesis for the CAVER comparison to settle, not a
finding.

## Independent consistency checks on the other three structures

- **CviUPO is narrowest (4.99 Å) and has zero naphthalene activity.** Independently, Figure S13 of
  the 2025 study describes CviUPO as having four bulky residues crowding the heme. Two unrelated
  lines agreeing.
- **HspUPO is widest (6.98 Å) with naphthalene activity 100.**
- **But MroUPO (5.91 Å, naphthalene 100) and PaDa-I (6.00 Å, naphthalene 3) have the same channel
  and opposite activity.** Any apparent correlation across these four (r = +0.83) is carried
  entirely by CviUPO and is not evidence.

## Status and what is blocked

**The rule is not transferred.** Producing ring-size categories from these numbers as though they
were Karich's would repeat the aliphatic-pocket error: using a feature whose calibration had
already failed.

**CAVER 3.0.2 is now downloaded and installed** under `tools/caver_dl/` — no registration was
needed, the earlier 404 was a wrong URL path — with a wrapper at `upo_sites/bin/run_caver.py` that
cleans the structure, seeds from the heme iron and parses the tunnel bottleneck. **It has not been
run**: executing the downloaded jar was denied by the sandbox as external code and needs the user's
permission.

The decisive test is one command. If CAVER returns ~7 Å and ~11 Å for these two structures, the
rule transfers and channel width can be measured across clade representatives with AlphaFold3 and
Boltz-2 (both installed with weights). If CAVER also returns ~6 Å and ~6 Å, then the published
figures are not a bottleneck measure, and the ring rule cannot be applied structurally by anyone —
which is worth knowing before further work is built on it.

**Either way, the n = 2 calibration stands.** Only two enzymes anchor the width-to-activity
relation, so even a successful transfer yields a mechanistic prediction with a two-point
calibration, not a validated classifier.
