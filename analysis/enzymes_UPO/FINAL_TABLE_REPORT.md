# Final BFD UPO table: composition

**97520 PF01328 proteins.**

## Clade (phylogenetic placement)

| clade | n |
|---|---|
| unassigned | 76913 |
| V | 8272 |
| IV | 5473 |
| II | 2737 |
| I | 2565 |
| III | 1560 |

## Trust tier — the column to filter on

| tier | n | % |
|---|---|---|
| `unassigned` | 76913 | 78.9 |
| `extrapolated` | 9603 | 9.8 |
| `reliable` | 5531 | 5.7 |
| `reliable+provisional_cladeIV` | 4924 | 5.0 |
| `extrapolated+provisional_cladeIV` | 549 | 0.6 |

## Independent corroboration of the trust boundary

The table carries two clade calls made by different methods: `clade` from EPA
phylogenetic placement, and `nn_clade` from identity to the nearest reference.
They are not derived from each other. Their agreement is therefore an external
check on where the `reliable` boundary sits.

| tier | agree | disagree | agreement |
|---|---|---|---|
| `reliable` | 10,411 | 44 | 99.6 % |
| `extrapolated` | 7,924 | 2,228 | 78.1 % |

The two methods converge above the boundary and diverge below it. This is
consistent with the leave-one-out precision measured on the references, and it
was measured on BFD proteins rather than on references.

It does NOT show the calls are correct. Two methods can agree and both be
wrong, and both rest on the same 19 characterized enzymes. It shows only that
the boundary drawn from reference tests reproduces on real data.

## Placement declines most of the database

76,913 of 97,520 proteins (78.9 %) are `unassigned`: no placement reached the
0.5 likelihood weight threshold. That is the caller refusing, not failing.
PF01328 is wider than the UPO clades, so most of these proteins are not
expected to belong to any of them.

10,455 proteins (10.7 %) carry a `reliable` call. That is the usable set.

## Catalytic core within the reliable set

| clade | n reliable | with catalytic core | % |
|---|---|---|---|
| I | 2,194 | 1,971 | 89.8 |
| II | 244 | 242 | 99.2 |
| III | 1,004 | 853 | 85.0 |
| IV | 4,924 | 4,567 | 92.7 |
| V | 2,089 | 1,461 | 69.9 |

Read this as composition, not as activity. `catalytic_core` means the four
catalytic residues are present, which is necessary for peroxygenase chemistry
and not sufficient for it. No assay supports a rate or a substrate range for
any of these proteins. Clade IV additionally has ONE characterized member, so
every clade IV row stays flagged `+provisional_cladeIV`.
