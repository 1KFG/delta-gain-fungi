# `upo_sites`

Find fungal unspecific peroxygenases (PF01328), annotate their active-site residues in HspUPO
numbering, assign them to Clades I-V, and score the aliphatic pocket.

## What this tool claims, and what it does not

It reports **measured sequence features** and a **clade assignment with a stated confidence tier**.

**It does not predict activity, and it will not be made to.** The reason is measured, not cautious:

| Question | Result on the 19-enzyme calibration set |
|---|---|
| Does the aliphatic pocket predict aromatic ring hydroxylation? | **No.** Pearson r = +0.176 |
| Does it predict peroxygenase-vs-peroxidase propensity? | r = +0.748 raw, but **partial r = +0.129 given clade** |
| Does clade predict peroxygenase propensity? | r = +0.835 |

Clade carries the signal; the pocket is nearly collinear with it and adds almost nothing. With
n = 19 across 5 clades, "this residue matters" cannot be separated from "this clade behaves this
way". A per-protein activity score would therefore be clade laundered through residue arithmetic,
and would rank 97,520 proteins by a feature with no measured predictive power. See
`../CALIBRATION_MATRIX.md` sections 4, 5 and 8.

## Install / environment

```
module load hmmer/3.4          # hmmsearch, hmmalign
/usr/bin/python3.12            # required; default python3 on HPCC is 3.9
```
No other dependencies. `zstd` is used only if the output path ends in `.zst`.

## Usage

Search a proteome and annotate what it finds:

```
/usr/bin/python3.12 upo_sites/bin/upo_sites.py \
    --proteome my_proteome.faa.gz \
    --hmm   upo_id/models/PF01328.hmm \
    --backbone data/clade_backbone/raw.fa \
    --out   results.tsv.zst
```

Annotate sequences already known to be PF01328 (for example pulled from the BFD pfam table), which
skips hmmsearch:

```
    --pre-selected
```

Accepts `.fa`, `.gz`, `.zst`. Multiple `--proteome` arguments are allowed.

## Method

1. `hmmsearch` PF01328 at the Pfam gathering threshold (or `--no-ga -E`).
2. `hmmalign` the hits in chunks, each chunk aligned **together with the 19-enzyme calibration
   backbone**, so every chunk carries its own references.
3. Read residues at fixed HspUPO positions. The HspUPO-position to HMM-match-state map is computed
   **once from the profile**, so it is identical in every chunk and at any input size. This is what
   makes chunking safe and the cost linear.
4. Assign clade by nearest neighbour to the backbone over match-state columns, reporting identity,
   margin over the best member of any other clade, and a confidence tier.

Only HMM **match states** are comparable across sequences; insert columns are not homologous.

## Positions

Numbering is HspUPO (UniProt A0A1Y2TH07 / NCBI OTA57433.1, PDB 7O1R, 261 aa). Location labels are
the paper's own, from Figure S2 of DOI 10.1186/s13068-025-02675-w.

| Location | Positions |
|---|---|
| Aliphatic pocket (AP) | 73, 77, 81, 183, 187, 228, 232; 223 is AP/E |
| Entrance (E) | 80, 84, 172, 175 |
| Catalytic pocket (CP) | 108, 176, 179 |
| Acid / base | 180 (Glu, universal); 110 (His, short) or 173 (Arg, long) |
| Proximal cysteine | 39 |

**21 requested, 18 scoreable. Positions 223, 228 and 232 fall in insert columns of PF01328** and
have no cross-sequence meaning. The tool detects this from the profile and reports it at startup.
The pocket score is therefore over the five reliable positions **73, 77, 81, 183, 187**.

## Validation

**Against the published structure-based matrix.** At the five pocket positions, across all 19
calibration enzymes, `hmmalign` agrees with Figure S2 at **95/95 = 100 %**. Over all 18 scoreable
positions, agreement is 95.5 %.

**Clade caller, leave-one-out: 18/19 = 94.7 %.** Run it yourself:

```
/usr/bin/python3.12 upo_sites/bin/validate_clade_loo.py \
    upo_id/models/PF01328.hmm data/clade_backbone/raw.fa
```

The single failure is MphUPO, the only characterized Clade IV enzyme: leave-one-out deletes its
entire clade, so it cannot be called correctly by construction. On clades with two or more members
the result is 18/18.

### The limitation that matters most

**The validation does not cover the regime the tool runs in.** Correct leave-one-out calls span
**52.7-96.2 %** identity to the nearest reference. BFD proteins sit at roughly **25-35 %**. On an
87-protein BFD sample, **only 1 of 87 reached 52.7 %**. Nearest-neighbour assignment is validated
at high identity and extrapolated almost everywhere it is actually used.

Margin does not rescue this: correct calls had margins from 9.9 upward while the single wrong call
had a margin of 17.7, so **margin did not separate correct from incorrect** at n = 19. The
`--min-margin` default is a heuristic, not a calibrated threshold.

The tool handles this by refusing to hide it, through `clade_confidence`:

| Tier | Meaning |
|---|---|
| `validated_range` | identity >= 52.7, inside the leave-one-out validated range |
| `extrapolated` | identity 30-52.7, below the validated range; treat as a hypothesis |
| `below_identity_floor` | identity < 30; clade reported as `unassigned` |
| `..._low_margin` suffix | margin below `--min-margin` (heuristic) |
| `+provisional_cladeIV` suffix | Clade IV has one characterized member; never settled |

**The planned fix is phylogenetic placement** (`epa-ng` / `pplacer`) against a fixed reference tree,
which gives a likelihood-weighted assignment and a real uncertainty. See `../PLAN.md` section 3.1.

## Output columns

`sample`, `protein_id`, `hmm_evalue` (`NA` under `--pre-selected`), `length`, `clade`,
`clade_identity`, `clade_margin`, `clade_confidence`, `nearest_ref`, `pocket_score`, `pocket_max`,
`pocket_residues`, `catalytic_ok`, `catalytic_detail`, `sites_mapped`, `clade_T1_insample`, then one
column per scoreable position named `p<N>_<location>`.

`clade_T1_insample` is the **observed** peroxygenase range for that clade in the 19-enzyme set
(Clade V 33-82, Clade I 0-53 and bimodal, Clades II-IV 0-5). It is a property of the calibration
set, not a prediction about the protein in that row.

## Files

```
bin/upo_sites.py            the annotator
bin/validate_clade_loo.py   leave-one-out validation of the clade caller
bin/extract_pf01328.py      pull PF01328 sequences out of the BFD seqdb into chunks
slurm/run_bfd.sh            SLURM array driver for the full BFD run
```
