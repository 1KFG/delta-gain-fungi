# Applicability domain of the T1 pocket score

19 calibration enzymes with both an active-site profile and a LOO prediction.
Distance positions (disjoint from the scored POCKET5): 80, 84, 108, 110, 172, 173, 175, 176, 179, 180

## Does distance from the calibration set predict the model's error?

Nearest-neighbour distance among the 19: median 0.20, range 0.00-0.50.

Distance is confounded with clade: 7 of 7 enzymes at distance <= 0.10 are clade V, against 1 of 8 at >= 0.30. Clade V holds 9 of the 10 actives, whose activities sit far from the mean. Absolute error is larger there for ANY model, including the null. The raw correlation below therefore measures the activity distribution, not the model.

Raw Spearman(distance, |LOO error|) = -0.384, one-sided permutation p = 0.9435 (10000 shuffles). **Confounded - do not read this as a result.**

### The unconfounded version: does the model's SKILL decay with distance?

Skill = |null-model error| - |pocket-model error|. Positive means the pocket model beat the
baseline on that enzyme. Both models face the same response, so the magnitude of the
response cancels.

Mean skill, distance <= 0.10: +5.6 (n = 7); distance >= 0.30: +13.2 (n = 8).

Spearman(distance, -skill) = -0.281, one-sided permutation p = 0.8789.

**Distance does NOT predict error on this calibration set.** This is a null result and is
reported as one. The domain flag below is therefore a statement about COVERAGE - whether a
protein's active site resembles anything ever measured - and NOT a calibrated confidence.
It must not be presented as a probability of being right.

| enzyme | clade | active | NN distance | \|LOO error\| |
|---|---|---|---|---|
| LspUPO-II | I | 0 | 0.50 | 1.7 |
| RneUPO | V | 1 | 0.50 | 5.8 |
| CabUPO-III | I | 0 | 0.40 | 1.7 |
| CviUPO | II | 0 | 0.30 | 12.0 |
| MorUPO | III | 0 | 0.30 | 10.2 |
| MroUPO | III | 0 | 0.30 | 0.5 |
| MphUPO | IV | 0 | 0.30 | 19.1 |
| PanUPO | II | 0 | 0.30 | 12.0 |
| GmaUPO-II | I | 1 | 0.20 | 41.3 |
| CmaUPO-I | I | 0 | 0.20 | 0.5 |
| AbrUPO-I | V | 1 | 0.20 | 15.4 |
| PaDa-I | I | 0 | 0.20 | 0.3 |
| TteUPO | V | 1 | 0.10 | 3.0 |
| AtuUPO | V | 1 | 0.10 | 24.6 |
| HspUPO | V | 1 | 0.10 | 6.9 |
| AluUPO | V | 1 | 0.10 | 30.6 |
| DspUPO-I | V | 1 | 0.10 | 4.1 |
| AbrUPO-II | V | 1 | 0.00 | 13.7 |
| AacUPO | V | 1 | 0.00 | 29.2 |

The largest nearest-neighbour distance inside the calibration set is 0.50. A BFD protein further than that from all 19 enzymes sits outside the range of active-site divergence the score was ever calibrated over.

## How much of BFD lies outside that range?

97520 proteins read from the site layer.

| domain | n | % |
|---|---|---|
| inside | 28845 | 29.6 |
| outside | 67854 | 69.6 |
| unknown (no profile) | 821 | 0.8 |

**3284 of the 8192 proteins scoring >= 3/5 (40.1 %) sit OUTSIDE the calibration range.** For those, a high pocket score is an extrapolation on both axes at once: the clade is usually not V, and the rest of the active site resembles nothing that has been measured.

Outside-domain proteins by clade call:

| clade | n |
|---|---|
| unassigned | 62353 |
| V | 2521 |
| II | 1601 |
| III | 1028 |
| IV | 246 |
| I | 105 |

### How to read the two columns this adds

`calib_distance` is the fraction of differing residues, over the 10 non-scored active-site
positions, to the nearest of the 19 characterized enzymes.

`calib_domain` is `inside` when that distance is within the largest distance seen WITHIN the
calibration set (0.50), `outside` otherwise.

`outside` does NOT mean a protein is inactive, and `inside` does NOT mean the score is
right about it. Distance did not predict model error on 19 enzymes. The column says only
whether the protein's active site resembles something that has been measured. Use it to
decide what to assay, not to decide what is active.

