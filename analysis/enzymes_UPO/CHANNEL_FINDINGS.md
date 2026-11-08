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

## CAVER run: the published numbers are not a channel bottleneck

**Measured 2026-11-08 with CAVER 3.0.2**, the tool the published figures came from, installed
system-wide at `/opt/linux/centos/8.x/x86_64/pkgs/caver/3.02/caver` (the user installed it; my
downloaded copy was blocked by the sandbox as external code). Wrapper:
`upo_sites/bin/run_caver.py`, seeding from the heme iron on chain A plus heme, waters removed.

| Structure | CAVER tunnels | bottleneck **radius** | tunnel length |
|---|---|---|---|
| AaeUPO (2YOR) | 3 | 1.06, 1.06, 0.95 Å | 5.3, 5.8, 28.4 Å |
| MroUPO (5FUJ) | 1 | 1.00 Å | 3.9 Å |

**CAVER gives bottleneck radii of about 1 Å — roughly 2 Å diameter — for both enzymes, and AaeUPO
is marginally the wider of the two.** That is the opposite ordering to the published 7 versus 11 Å,
and an order of magnitude smaller.

A probe sweep confirms it is not a parameter choice. At probe radii of 1.5 Å and above, **CAVER
finds no tunnel at all in either structure**:

```
probe 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 A  ->  no tunnel, both enzymes
```

If the channels were 7 Å and 11 Å across, probes of 3.5 Å and 5.5 Å radius would pass.

### Conclusion: the ring-size rule cannot be transferred structurally

Two independent measurements — CAVER and the hard-sphere implementation here — agree that **the
published 7 Å and 11 Å are not static-structure channel bottlenecks.** They are most plausibly
approximate descriptions of the channel mouth, or a measure taken on a different representation.
Whatever they are, they cannot be recomputed for an uncharacterized enzyme, which is what applying
the rule across BFD would require.

The two methods also disagree sharply with each other (CAVER 1.0 Å radius versus 3.0 Å here), so
the hard-sphere measure is the more permissive of the two and should not be presented as
equivalent to CAVER.

**This also makes the planned AlphaFold3 / Boltz-2 modelling pointless for this purpose.** If the
measurement cannot reproduce the published contrast on experimentally determined crystal
structures, predicted structures will not either. The structure work is therefore **not run** —
that is a saving, not an omission.

**The one surviving observation** is that MroUPO's only tunnel is 3.9 Å long against AaeUPO's
5.3–28.4 Å. MroUPO's heme is markedly more solvent-exposed. That is consistent with it accepting
bulkier substrates and may be what the published "wider channel" describes, but tunnel length is
not channel width and no rule should be built on a two-enzyme comparison.

## Status and what is blocked

**The rule is not transferred.** Producing ring-size categories from these numbers as though they
were Karich's would repeat the aliphatic-pocket error: using a feature whose calibration had
already failed.

**CAVER has now been run** (see the section above) and the question is settled: the published
figures are not a bottleneck measure, so the ring rule cannot be applied structurally. The
AlphaFold3 / Boltz-2 modelling that would have followed is not run, because it cannot succeed where
the crystal structures already fail.

**Either way, the n = 2 calibration stands.** Only two enzymes anchor the width-to-activity
relation, so even a successful transfer yields a mechanistic prediction with a two-point
calibration, not a validated classifier.
