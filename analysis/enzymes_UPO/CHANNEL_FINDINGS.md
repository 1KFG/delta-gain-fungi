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

## CAVER run — CORRECTED 2026-11-08: version 3.0.2 carried a known bug

**An earlier version of this section concluded that the published 7 Å / 11 Å "are not channel
bottlenecks" and that the ring-size rule could not be transferred structurally. That conclusion
was wrong. It rested on CAVER 3.0.2 output affected by a documented bug.**

CAVER's own changelog for 3.01 states:

> "The starting point optimization procedure was corrected. Previously, the starting vertex was
> sometimes identified in **too narrow cavity** due to error in implementation."

That is exactly the symptom observed: 3.0.2 returned ~1 Å bottleneck radii for every structure and
found no tunnel above a 1.5 Å probe. The starting vertex was being planted in a sub-ångström
crevice under the porphyrin.

### CAVER 3.0.3 results

| Structure | bottleneck radius | **channel diameter** | naphthalene activity |
|---|---|---|---|
| MroUPO (5FUJ) | 2.59 Å | **5.19 Å** | 100 |
| MroUPO (7ZBP) | 2.45 Å | 4.91 Å | 100 |
| HspUPO (7O1R) | 2.44 Å | **4.87 Å** | 100 |
| AaeUPO/PaDa-I (2YOR) | 2.12 Å | **4.23 Å** | 3 |
| CviUPO (7ZCL) | 1.85 Å | **3.71 Å** | 0 |

**The ordering now matches the published claim: MroUPO's channel is wider than AaeUPO's**, in both
MroUPO structures, which is the comparison the whole ring-size rule rests on and which 3.0.2
reversed.

Against measured naphthalene activity across the four enzymes with both a structure and an assay:
**Pearson r = +0.933**, and the ranking is monotonic — Mro > Hsp > Aae > Cvi against 100, 100, 3, 0.

### What this does and does not re-open

**Re-opened:** a *relative* channel-width ranking is now defensible. Wider channel associates with
ring activity in the right direction, with a mechanistic rationale and a 4-point correlation.

**Still not reproduced:** the absolute values. Measured 4.23 Å and 5.19 Å against published 7 Å and
11 Å; the measured ratio is 1.23 against a published 1.57. So the published figures remain a
different quantity, and **Karich's absolute cutoff — "PAHs up to about 6 Å minimum width" — still
cannot be applied to these numbers.** What can be done is rank enzymes by channel width and treat
the ranking, not the absolute value, as the signal.

**n = 4.** Four enzymes have both a solved structure and a naphthalene measurement. An r of +0.933
on four points is a lead, not a validated relationship, and two of the four (Mro, Hsp) share the
same activity value so the correlation rests on Cvi and Aae at the low end.

### The methodological lesson, recorded

Two independent measurements agreeing is not sufficient evidence when one of them may be wrong for
a reason neither test can see. My hard-sphere tool and CAVER 3.0.2 agreed on ~6 Å and ~2 Å
respectively — both failing to reproduce the published contrast — and I read that agreement as
confirmation. It was two different errors pointing the same way. **Checking the tool's changelog
before trusting its output would have caught this.**

## Status and what is blocked

**The rule is not transferred.** Producing ring-size categories from these numbers as though they
were Karich's would repeat the aliphatic-pocket error: using a feature whose calibration had
already failed.

**CAVER 3.0.3 has been run and the structural route is open again** for relative ranking, though
not for Karich's absolute cutoff. AlphaFold3 and Boltz-2 are installed with weights, so modelling
clade representatives and ranking them by CAVER 3.0.3 channel width is now the sensible next step —
with the ranking, not the absolute width, as the output.

**Either way, the n = 2 calibration stands.** Only two enzymes anchor the width-to-activity
relation, so even a successful transfer yields a mechanistic prediction with a two-point
calibration, not a validated classifier.

---

## The n = 4 lead does not survive the fifth structure

**Measured 2026-11-08, immediately after the correction above.** Rather than leave r = +0.933 on
four points as the standing result, I searched the PDB for more structures of enzymes that also
have naphthalene measurements. One exists: **9IZ8 = TteUPO** (*Thielavia terrestris*), one of the
19, measured naphthalene 0.

| Enzyme | channel diameter (CAVER 3.0.3) | naphthalene |
|---|---|---|
| **TteUPO (9IZ8)** | **5.40 Å — the widest measured** | **0** |
| MroUPO (5FUJ) | 5.19 Å | 100 |
| HspUPO (7O1R) | 4.87 Å | 100 |
| AaeUPO (2YOR) | 4.23 Å | 3 |
| CviUPO (7ZCL) | 3.71 Å | 0 |

| Set | r |
|---|---|
| All five structures | **+0.452** |
| Excluding TteUPO (the original four) | +0.933 |

**TteUPO has the widest channel of any enzyme measured and no naphthalene activity at all.** It is
a direct counterexample, and adding it halves the correlation.

### Is excluding TteUPO defensible?

Partly, and the argument was made **before** this measurement, not after it. `G1_FINDINGS.md`
Finding 5 and `CALIBRATION_MATRIX.md` §4.1 both record that **TteUPO and AbrUPO-II are at or near
zero on all five assays** (ABTS, NBD, 2,6-DMP, naphthalene, 4-NA), and that an enzyme flat across
every assay is most simply poorly expressed protein rather than one selectively unable to oxidise
naphthalene. On that basis TteUPO's zero carries no information about substrate specificity.

But the honest position is that **this is still a post-hoc exclusion that rescues the correlation**,
and it reduces an already tiny sample to four. The pre-registered reasoning helps; it does not make
n = 4 into evidence.

### Where this leaves the channel-width route

**Not established, and weaker than it looked an hour ago.** The correct summary is:

- CAVER 3.0.3 reproduces the published *ordering* for AaeUPO versus MroUPO, which 3.0.2 reversed.
- Among enzymes whose activity measurement is interpretable, wider channel tracks naphthalene
  activity (r = +0.93, n = 4).
- The single enzyme whose measurement is **not** interpretable is also the widest channel and a
  zero, so the relationship cannot be tested on it either way.
- **Five structures, of which one must be set aside on independent grounds, is not a basis for
  scoring 97,520 proteins.**

Also measured, for completeness, though neither enzyme is among the 19 and so neither has
naphthalene data: **MweUPO (8RNJ) 5.21 Å** and **PabUPO-II (9HE6) 4.63 Å**.

### Consequence for the modelling plan

Modelling clade representatives and ranking them by channel width would produce a ranking whose
only calibration is four enzymes, one of them contested. **I am not running it on that basis.**
What would change this is more enzymes with both a structure and a ring-substrate measurement —
the same reference-set limit that bounds every other part of this project.

### Tooling note, corrected

I reported that the Boltz-2 environment was broken with a numpy/numba conflict. **That was my
error.** The environment's own numpy is 1.26.4, compatible with its numba 0.61.0. The 2.5.2 that
triggered the ImportError came from the user site directory
(`~/.local/lib/python3.12/site-packages/numpy`), which Python prepends ahead of an environment's
own packages. `PYTHONNOUSERSITE=1` makes Boltz-2 run normally.

Worth knowing beyond this project: that user-site numpy will shadow **any** module-provided Python
environment on this account, and the resulting error message blames the environment rather than the
shadowing.

So Boltz-2 is available. The decision not to model clade representatives is about the evidence —
four calibration enzymes, one contested — not about tooling.

---

## Boltz-2 validation: predicted structures fail the channel measurement

**Measured 2026-11-08.** All 19 calibration enzymes modelled with Boltz-2 (heme as an explicit
ligand, MSAs built from 2,000 real BFD PF01328 homologs). The 5 with crystal structures were
modelled **deliberately**, as the test of whether predicted structures can be measured this way at
all. Inputs: `upo_sites/bin/build_boltz_inputs.py`; job: `upo_sites/slurm/boltz_19.sh`.

| Enzyme | predicted | experimental | error | Boltz confidence |
|---|---|---|---|---|
| HspUPO | 4.86 Å | 4.87 Å | **−0.01** | 0.926 |
| MroUPO | 5.07 Å | 5.19 Å | **−0.12** | 0.964 |
| AaeUPO | 4.05 Å | 4.23 Å | **−0.18** | 0.941 |
| CviUPO | 2.36 Å | 3.71 Å | **−1.35** | **0.968** |
| TteUPO | 2.49 Å | 5.40 Å | **−2.91** | 0.873 |

**Three of five land within 0.2 Å. Two are badly wrong.** Overall r = **+0.332**, mean absolute
error 0.91 Å, worst 2.91 Å.

### Why this kills the n = 19 plan rather than bruising it

**TteUPO is measured as the widest channel of any enzyme (5.40 Å) and predicted as the narrowest
(2.49 Å).** The prediction is not noisy, it is backwards — on the single enzyme that matters most,
since TteUPO is the counterexample that collapsed the width-activity correlation in the first
place.

**Model confidence does not flag the failures.** CviUPO has the **highest** confidence of all five
(0.968) and the second-worst error. All five score 0.87–0.97. So there is no confidence threshold
that separates the usable models from the unusable ones — which means the other 14 predictions
cannot be triaged, only discarded.

### Conclusion

**The other 14 predicted structures are not used.** Channel width can be measured from a crystal
structure and not from a predicted one, at least at the accuracy this question needs, and no
available score tells you which predictions to trust.

The channel-width route therefore stays at **n = 5 experimental structures** (n = 4 once TteUPO is
set aside on independent grounds), which is what it was before. **This is the same reference-set
ceiling that bounds every other part of this project**, reached by a third independent road: not
enough characterized enzymes, and here not enough solved structures either.

What would move it: more crystal structures of UPOs that also have ring-substrate measurements.
Nothing computational substitutes.
