# T1 pocket score across BFD

Consensus `ILLAI` at HspUPO p73, p77, p81, p183, p187. Cut at >= 3 matches.

## Score distribution

| matches | n | % |
|---|---|---|
| 0/5 | 21532 | 22.1 |
| 1/5 | 25249 | 25.9 |
| 2/5 | 42547 | 43.6 |
| 3/5 | 7486 | 7.7 |
| 4/5 | 549 | 0.6 |
| 5/5 | 157 | 0.2 |

**8192 proteins (8.4 %) score >= 3.**

## Where those hits sit, by clade call

| clade | n |
|---|---|
| V | 4182 |
| unassigned | 3924 |
| II | 39 |
| IV | 31 |
| I | 11 |
| III | 5 |

## By trust tier

| tier | n |
|---|---|
| `unassigned` | 3924 |
| `extrapolated` | 3014 |
| `reliable` | 1223 |
| `reliable+provisional_cladeIV` | 26 |
| `extrapolated+provisional_cladeIV` | 5 |

## Hits the placement caller DECLINED, by phylum

These are the proteins the score exists for: no clade call, but a measured pocket.

| phylum | n |
|---|---|
| Ascomycota | 2555 |
| Basidiomycota | 1337 |
| Chytridiomycota | 13 |
| Mucoromycota | 9 |
| Basidiobolomycota | 5 |
| Glomeromycota | 3 |
| Zoopagomycota | 2 |

## How the score was validated

`upo_sites/bin/score_t1.py`. Response is mean peroxygenase conversion of 4-propylguaiacol
over three pH values, Figure S8, n = 19. Four models, leave-one-out. The active consensus,
the clade means and the regression coefficients are all refit on the 18 training enzymes
every fold, so the held-out enzyme never informs its own prediction.

| model | LOO R² | LOO ρ | LOO AUC | MAE |
|---|---|---|---|---|
| null (training mean) | −0.114 | −0.984 | 0.000 | 22.4 |
| clade only | +0.434 | 0.521 | 0.922 | 13.5 |
| **pocket only** | **+0.496** | **0.628** | **0.956** | **12.3** |
| clade + pocket | +0.469 | 0.532 | 0.911 | 12.7 |

Permutation test, 2,000 shuffles of the response: pocket **p = 0.0005**.

In-sample R² is 0.589 against LOO 0.496, so the overfit gap is small. The pocket beats the
clade label out of sample, and adding the clade to the pocket makes it slightly worse.

### The dose-response on the calibration set

| matches | n | mean activity | active (>10 %) |
|---|---|---|---|
| 0/5 | 5 | 0.7 | 0/5 |
| 1/5 | 4 | 12.8 | 1/4 |
| 2/5 | 1 | 2.7 | 0/1 |
| 3/5 | 2 | 29.2 | 2/2 |
| 4/5 | 4 | 52.9 | 4/4 |
| 5/5 | 3 | 38.6 | 3/3 |

Monotone except the n = 1 cell. A cut at ≥ 3 gives 9/9 precision and 9/10 recall.

## The limit, which decides how the score may be used

Nine of the ten active calibration enzymes are clade V. Removing clade V:

| model | AUC (clade V removed) | AUC (all) |
|---|---|---|
| clade only | 0.222 | 0.922 |
| pocket only | **0.556** | 0.956 |

0.556 is chance, and it rests on ONE remaining active enzyme. The score is therefore a
clade-V-like-pocket detector. It is **not** demonstrated to find actives outside clade V, and
it misses the one that exists: GmaUPO-II, clade I, 49.3 % conversion, 1 of 5 matches.

Two consequences, both binding:

1. A high score is a clade V hypothesis reached by a second route. Where the placement caller
   already says clade V, the score adds little.
2. **A low score is not evidence of inactivity.** GmaUPO-II proves this directly.

## What the score is actually for

The placement caller declines 78.9 % of BFD. Pocket residues are measured for every protein
regardless. For those 76,913 proteins the score is the only ranking signal that exists.

Enrichment against the clade call, as an internal consistency check:

| set | ≥3/5 | total | rate |
|---|---|---|---|
| clade V | 4,182 | 8,272 | 50.6 % |
| everything else | 4,010 | 89,248 | 4.5 % |

11.3-fold enrichment. The score is not a restatement of the clade label: it splits clade V
almost in half, and it flags 4,010 proteins outside it.

## The shortlist

Catalytic core present AND ≥ 3/5 pocket matches: **5,115 proteins across 1,378 species.**

3,077 of the 8,192 scoring proteins lack the catalytic core and are excluded here. A pocket
without the four catalytic residues cannot do peroxygenase chemistry whatever the score says.

| clade call | n |
|---|---|
| V | 3,869 |
| unassigned | 1,171 |
| II | 39 |
| IV | 23 |
| I | 8 |
| III | 5 |

1,246 of these sit outside clade V, including 29 in early-diverging phyla
(Chytridiomycota 10, Mucoromycota 9, Basidiobolomycota 5, Glomeromycota 3, Zoopagomycota 2).

This is a ranked candidate list for assay. It is not a prediction of activity, and it says
nothing about which rings any of these enzymes cleave. One substrate, one assay, one
publication, 19 calibration enzymes, 9 of 10 actives in a single clade.
