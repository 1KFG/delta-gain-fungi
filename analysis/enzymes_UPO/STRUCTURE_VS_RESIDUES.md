# Structure versus residues for T1 activity

19 enzymes with a Boltz structure, pocket residues and a T1 value. Active (> 10 %): 10.

## Leave-one-out, out of sample

| model | LOO R2 | LOO AUC (all) | **LOO AUC, clade V removed** |
|---|---|---|---|
| null | -0.114 | 0.000 | **0.000** |
| clade | +0.434 | 0.922 | **0.222** |
| pocket | +0.496 | 0.956 | **0.556** |
| channel | +0.108 | 0.689 | **0.333** |
| profile | -0.168 | 0.678 | **0.444** |
| pocket+chan | +0.490 | 0.944 | **0.444** |

Non-clade-V enzymes: 10, of which active: 1 (GmaUPO-II).

## Permutation test on the structural models

| model | LOO R2 | p |
|---|---|---|
| channel | +0.108 | 0.0233 |
| profile | -0.168 | 0.2040 |
| pocket+chan | +0.490 | 0.0004 |

## Per-enzyme channel geometry against activity

| enzyme | clade | active | activity | bottleneck | d3 | d5 |
|---|---|---|---|---|---|---|
| AluUPO | V | 1 | 67.7 | 2.39 | 4.9 | 4.3 |
| AtuUPO | V | 1 | 62.3 | 1.95 | 5.1 | 5.1 |
| AbrUPO-I | V | 1 | 54.0 | 4.85 | 11.2 | 15.0 |
| GmaUPO-II | I | 1 | 49.3 | 5.58 | 6.2 | 7.2 |
| DspUPO-I | V | 1 | 46.3 | 4.66 | 7.9 | 10.2 |
| HspUPO | V | 1 | 44.0 | 6.28 | 8.4 | 10.3 |
| TteUPO | V | 1 | 33.3 | 1.95 | 4.6 | 4.3 |
| AbrUPO-II | V | 1 | 27.7 | 1.95 | 4.9 | 4.5 |
| AacUPO | V | 1 | 25.3 | 6.28 | 12.2 | 14.9 |
| RneUPO | V | 1 | 25.0 | 3.98 | 4.9 | 4.6 |
| MphUPO | IV | 0 | 2.7 | 4.92 | 8.5 | 10.3 |
| PaDa-I | I | 0 | 1.7 | 8.79 | 9.4 | 10.9 |
| MorUPO | III | 0 | 1.7 | 4.78 | 5.6 | 7.0 |
| CmaUPO-I | I | 0 | 1.0 | 6.45 | 7.9 | 8.7 |
| MroUPO | III | 0 | 1.0 | 5.65 | 6.0 | 7.2 |
| PanUPO | II | 0 | 0.0 | 5.70 | 6.5 | 8.0 |
| CabUPO-III | I | 0 | 0.0 | 3.91 | 7.0 | 5.4 |
| LspUPO-II | I | 0 | 0.0 | 7.15 | 11.2 | 12.1 |
| CviUPO | II | 0 | 0.0 | 5.11 | 11.1 | 14.2 |

Bottleneck, active enzymes: mean 3.99 A (range 1.95-6.28). Inactive: mean 5.83 A (range 3.91-8.79).

## Verdict

Out-of-clade AUC by model: null 0.000, clade 0.222, pocket 0.556, channel 0.333, profile 0.444, pocket+chan 0.444.

The out-of-clade question rests on 1 active enzyme outside clade V. No AUC computed on 1 positive case is evidence, whichever model wins it. A structural feature that beat the residue score here would be a HYPOTHESIS worth testing on new enzymes, not a validated out-of-clade model.


---

## Conclusion (written after the measurement; the design was fixed before it)

**Structure does not close the out-of-clade gap. This was the last untried computational route,
and it is now closed by measurement rather than by argument.**

1. **Channel geometry carries a real but small signal.** Bottleneck alone reaches LOO R2 +0.108,
   permutation p = 0.023. That is better than chance and worse than the residue score's +0.496
   on the same folds.
2. **It is worse exactly where it was needed.** Out-of-clade AUC: channel 0.333, profile 0.444,
   against the residue score's 0.556. The hope was that geometry would transfer across clades
   because a channel is a physical constraint. It did not.
3. **Combining hurts.** `pocket+chan` scores LOO R2 +0.490 against pocket's +0.496, and
   out-of-clade AUC 0.444 against 0.556. Adding the structural feature makes the model worse on
   both axes, so it is not adopted.

### An observation that runs against the obvious reading

Active enzymes have **narrower** bottlenecks (mean 3.99 A, range 1.95-6.28) than inactive ones
(mean 5.83 A, range 3.91-8.79). The intuitive direction — a wider channel admits more substrate —
is the opposite of what is measured here.

This is not presented as a finding. It is confounded with clade in the same way everything else
in this calibration set is: clade V holds 9 of the 10 actives, and clade V enzymes may simply
have narrower channels for reasons unrelated to catalysis. It is also measured on apo predicted
models with no substrate bound. It is recorded because it is the kind of result that gets
quietly dropped for disagreeing with expectation.

### A bug worth recording

The first run of this analysis was wrong. `heme_channel.py` prints each structure's corridor
profile BEFORE its summary line, and the extraction attached every profile to the previous
enzyme, shifting all profile features by one structure and dropping the last enzyme entirely.
The `profile` model's first numbers were therefore meaningless. It was caught by noticing that
TteUPO had been excluded for a missing profile when its structure plainly existed.

Any parser of this tool's output must attach a profile block to the summary line that FOLLOWS
it, not the one preceding it.

### What this means for the model

The structural route was the one the Boltz-2 run was built for, and the one remaining reason to
think the out-of-clade problem was solvable without new biochemistry. It is now measured and it
does not work. Combined with `EXTERNAL_VALIDATION.md`, every computational route is closed:

| route | outcome |
|---|---|
| second panel, same assay | does not exist |
| external labels from other papers | invalid: labels are substrate-specific (PaDa-I) |
| applicability distance as confidence | failed, p = 0.88 |
| **structural features from 19 Boltz models** | **failed: out-of-clade 0.333, worse than 0.556** |
| assay 6 existing enzymes on 4-PG | open, cheap, predictions recorded |

The out-of-clade question still rests on ONE active enzyme outside clade V. No AUC computed on
one positive case is evidence, including the 0.556 the residue score posts. The honest statement
is unchanged: the score detects clade-V-like pockets, and whether that generalizes is a question
only measurement on new enzymes can answer.
