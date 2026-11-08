# G1: fungal PF01328 diversity across BFD

Generated from `data/bfd_full/out`. **n = 97520 proteins.**

Clade calls below 52.7 % identity to the nearest characterized enzyme are **extrapolations beyond the validated range** of the caller and are reported separately throughout. See `upo_sites/README.md`.

## 1. Confidence tiers

| tier | n | % |
|---|---|---|
| `extrapolated_low_margin` | 33275 | 34.1 |
| `extrapolated` | 30413 | 31.2 |
| `below_identity_floor` | 19060 | 19.5 |
| `validated_range` | 7829 | 8.0 |
| `validated_range+provisional_cladeIV` | 4960 | 5.1 |
| `extrapolated_low_margin+provisional_cladeIV` | 1671 | 1.7 |
| `extrapolated+provisional_cladeIV` | 283 | 0.3 |
| `validated_low_margin` | 21 | 0.0 |
| `validated_low_margin+provisional_cladeIV` | 8 | 0.0 |

**Only 12818 of 97520 (13.14 %) clade calls fall inside the validated identity range.** Everything else is a hypothesis.

**Of those 12818, 4968 (38.8 %) are Clade IV**, which has one characterized member and therefore cannot be validated by leave-one-out at all. Excluding Clade IV leaves **7850** calls (8.05 % of all proteins) that are both inside the validated range and anchored on a clade with more than one reference.

## 2. Clade counts

**Clade IV carries a single characterized member (MphUPO), so its `within validated range` count means 'high identity to one uncharacterized protein'. It is marked provisional everywhere below and must not be pooled with the other clades' validated counts.**

| clade | all calls | within validated range | in-sample T1 (n=19) |
|---|---|---|---|
| I | 33424 | 3731 | 0-53, bimodal |
| V | 24847 | 2637 | 33-82 |
| unassigned | 19060 | 0 | n/a |
| III | 10154 | 1079 | 0-3 |
| IV **(provisional, 1 ref)** | 6922 | 4968 | 0-5 |
| II | 3113 | 403 | 0 |

## 3. Pocket score (0-5) by clade -- all calls

| clade | n | mean | 5/5 | 0/5 |
|---|---|---|---|---|
| I | 33424 | 3.29 | 498 | 106 |
| V | 24847 | 3.59 | 6708 | 8 |
| unassigned | 19060 | 2.82 | 2662 | 1417 |
| III | 10154 | 3.40 | 1799 | 4 |
| IV | 6922 | 3.79 | 186 | 0 |
| II | 3113 | 3.96 | 1188 | 1 |

## 4. Clade distribution by lineage

| lineage | n | I | V | unassigned | III | IV | II |
|---|---|---|---|---|---|---|---|
| Sordariomycetes | 32590 | 7783 | 12571 | 4149 | 5065 | 2531 | 491 |
| Dothideomycetes | 29855 | 15054 | 4929 | 5378 | 875 | 2464 | 1155 |
| Agaricomycetes | 15976 | 5526 | 3184 | 4140 | 2400 | 33 | 693 |
| Eurotiomycetes | 10403 | 1274 | 3003 | 2803 | 1293 | 1485 | 545 |
| Leotiomycetes | 2734 | 972 | 504 | 759 | 110 | 364 | 25 |
| Lecanoromycetes | 2176 | 1231 | 189 | 518 | 187 | 18 | 33 |
| Pucciniomycetes | 498 | 8 | 3 | 482 | 2 | 0 | 3 |
| Microbotryomycetes | 475 | 276 | 28 | 123 | 10 | 1 | 37 |
| Tremellomycetes | 459 | 434 | 0 | 17 | 7 | 1 | 0 |
| Exobasidiomycetes | 438 | 357 | 1 | 80 | 0 | 0 | 0 |
| Mucoromycota | 370 | 15 | 147 | 119 | 2 | 8 | 79 |
| Ustilaginomycetes | 357 | 218 | 3 | 136 | 0 | 0 | 0 |
| Chytridiomycota | 188 | 13 | 82 | 51 | 25 | 9 | 8 |
| Orbiliomycetes | 170 | 0 | 24 | 77 | 55 | 0 | 14 |
| Arthoniomycetes | 132 | 76 | 26 | 11 | 11 | 1 | 7 |
| Taphrinomycetes | 127 | 80 | 0 | 25 | 10 | 0 | 12 |
| Glomeromycota | 124 | 0 | 71 | 11 | 40 | 0 | 2 |
| Ascomycota | 89 | 33 | 20 | 23 | 3 | 5 | 5 |
| Basidiobolomycota | 86 | 0 | 26 | 16 | 44 | 0 | 0 |
| Wallemiomycetes | 69 | 3 | 0 | 65 | 1 | 0 | 0 |
| Pezizomycetes | 65 | 31 | 1 | 30 | 1 | 1 | 1 |
| Zoopagomycota | 30 | 1 | 7 | 19 | 2 | 0 | 1 |
| Dacrymycetes | 24 | 15 | 4 | 4 | 1 | 0 | 0 |
| Basidiomycota | 23 | 9 | 6 | 8 | 0 | 0 | 0 |
| Geoglossomycetes | 18 | 1 | 6 | 3 | 7 | 0 | 1 |
| Candelariomycetes | 5 | 2 | 0 | 2 | 1 | 0 | 0 |
| Coniocybomycetes | 5 | 4 | 1 | 0 | 0 | 0 | 0 |
| Monoblepharomycota | 4 | 1 | 0 | 3 | 0 | 0 | 0 |
| ? | 4 | 0 | 2 | 2 | 0 | 0 | 0 |
| Xylonomycetes | 3 | 3 | 0 | 0 | 0 | 0 | 0 |
| Saitoellomycetes | 3 | 0 | 3 | 0 | 0 | 0 | 0 |
| Aphelidiomycota | 3 | 0 | 1 | 1 | 1 | 0 | 0 |
| Lipomycetes | 2 | 0 | 1 | 1 | 0 | 0 | 0 |
| Lichinomycetes | 2 | 0 | 1 | 0 | 0 | 0 | 1 |
| Calcarisporiellomycota | 2 | 0 | 1 | 0 | 0 | 1 | 0 |
| Cystobasidiomycetes | 2 | 1 | 0 | 1 | 0 | 0 | 0 |
| Malasseziomycetes | 2 | 1 | 0 | 1 | 0 | 0 | 0 |
| Mixiomycetes | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| Classiculomycetes | 1 | 0 | 0 | 1 | 0 | 0 | 0 |
| Geminibasidiomycetes | 1 | 0 | 0 | 1 | 0 | 0 | 0 |
| Pichiomycetes | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| Laboulbeniomycetes | 1 | 0 | 0 | 0 | 1 | 0 | 0 |
| Xylobotryomycetes | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| Atractiellomycetes | 1 | 0 | 1 | 0 | 0 | 0 | 0 |

## 5. Catalytic machinery

- Full catalytic set intact (Cys39, Glu109, Asp111, Glu180): **34561/97520 (35.4 %)**
- Proximal Cys39 present: **89121/97520 (91.4 %)**
- Residue at position 110 (His in short UPOs, Gly in long): H=48866, G=16519, A=9533, S=8584, V=4230, -=2512

## 6. The catalytic-core filter -- which of these are plausibly UPOs at all

PF01328 is a broad heme-thiolate peroxidase family. The 2025 paper states that **Glu180 is present in all 19** characterized UPOs, so it is a usable test of whether a PF01328 protein has UPO-like catalytic machinery.

| confidence tier | n | Glu180 |
|---|---|---|
| validated | 12818 | 93.7 % |
| extrapolated | 65642 | 39.3 % |
| unassigned | 19060 | 21.8 % |

| clade | n | Glu180 |
|---|---|---|
| I | 33424 | 16.8 % |
| V | 24847 | 69.5 % |
| unassigned | 19060 | 21.8 % |
| III | 10154 | 55.0 % |
| IV | 6922 | 95.7 % |
| II | 3113 | 88.1 % |

**Validated-range calls carry Glu180 at 93.7 %, matching the published claim. Clade I calls carry it at 16.8 %.** The same alignment gives all five characterized Clade I enzymes a Glu at 180, so this is not an alignment artifact for long UPOs. The most likely reading is that the low-identity Clade I bin is absorbing PF01328 proteins that are **not UPOs**, because the long characterized enzymes are the nearest available match for divergent sequences. Treat unfiltered Clade I counts as an upper bound.

### Catalytic core set

Requiring Cys39 + Glu109 + Asp111 + Glu180: **34561/97520 (35.4 %)** proteins.

| clade | catalytic core | of which validated-range |
|---|---|---|
| V | 15596 | 1951 |
| IV **(provisional, 1 ref)** | 6002 | 4565 |
| I | 5101 | 3418 |
| III | 4540 | 909 |
| II | 2419 | 399 |
| unassigned | 903 | 0 |

**This is the defensible UPO set for downstream work**, and it is the one to carry into novelty and regulatory analysis, not the raw 97,520.


## 7. Copy number per genome

- genomes with at least one PF01328: **15137**
- maximum copies in one genome: **144**

| genome | copies | species |
|---|---|---|
| FC847EB9 | 144 | *Sphaerobolus stellatus* |
| FA6CE80F | 106 | *Mycena vulgaris* |
| F05B6213 | 105 | *Mycena floridula* |
| F476BA95 | 91 | *Mycena purpureofusca* |
| F96223D5 | 87 | *Mycena leptocephala* |
| FE13C428 | 78 | *Mycena rebaudengoi* |
| F39AF692 | 67 | *Mycena venus* |
| FFDFC078 | 66 | *Coprinopsis sp. MPI-PUGE-AT-0042* |
| F41267C2 | 61 | *Candolleomyces eurysporus* |
| FA63C8CD | 57 | *Candolleomyces candolleanus* |
