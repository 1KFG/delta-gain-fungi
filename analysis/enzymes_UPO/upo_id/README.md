# upo_id

Classify fungal heme-thiolate peroxidase (PF01328) proteins into UPO classes.

- Primary scheme: Kinner et al. 2021 Family I (short) and Family II (long).
- Secondary scheme: UPObase (Faiza et al. 2019a,b) subfamilies I–V and Pog-a..g, in separate columns.
- Every check has its own column. The tool reports all classes; it does not discard non-UPOs.

The method, its sources and every deviation from UPObase are in `METHOD.md`.

## Install (UCR HPCC)

Nothing to install. The tool uses the system Python 3.12 (it has Biopython 1.88 and pyarrow) and HMMER/MAFFT modules.

```bash
module load hmmer/3.4          # hmmsearch, hmmalign, phmmer, hmmbuild
module load mafft/7.505        # only for 'build'
PY=/usr/bin/python3.12         # default python3 on HPCC is 3.9; do not use it
```

`pixi.toml` (workspace format) lists the same dependencies for use off HPCC. It was **not tested**.

## Usage

Classify proteins (plain, `.gz` or `.zst` FASTA):

```bash
$PY upo_id/upo_id.py classify \
    -i proteins.fa.zst \
    -m models \
    -o out.tsv.zst \
    --parquet out.parquet \
    --sp-table /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/signalp.signal_peptide.parquet \
    --cpu 4
```

- Default output: one row per HTP candidate (PF01328 GA hit, or a family-HMM hit with E ≤ 1e-5).
  Add `--all` for one row per input protein.
- `--sp-table`: Parquet with a `protein_id` column, or a TSV/list whose first column is the protein ID.
  Without it, `signal_peptide` is `NA`.
- Temporary files go to `--tmpdir`, else `$SCRATCH`, else `$TMPDIR`.
- `--no-phmmer` skips the UPObase-style phmmer seed search (reported only, not used in the decision).

Rebuild the models (about 25 s):

```bash
$PY upo_id/upo_id.py build --refs ref/refs_all.fa --train ref/train_sets.tsv \
    --manifest ref/manifest.tsv --pfam-hmm ref/PF01328.hmm --out models --cpu 4
```

Refetch the reference sequences (writes a new manifest with checksums and dates):

```bash
$PY scripts/fetch_refs.py --sources ref/ref_sources.tsv --out-fasta ref/refs_all.fa --out-manifest ref/manifest.tsv
```

`ref/PF01328.hmm` came from `hmmfetch $PFAM_DB/Pfam-A.hmm Peroxidase_2` with `module load db-pfam` (Pfam 38.2,
PF01328.23).

Leave-one-out test (one SLURM job; measured 6 min on 8 CPUs):

```bash
sbatch tests/run_loo.sh /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/upo_id
```

## Output columns

Positions are AaeUPO **precursor** numbering in the code (`SITES`); METHOD.md gives mature numbering (precursor − 43).

| Column(s) | Meaning |
|---|---|
| `protein_id`, `length`, `mw_kda` | ID, length, rough mass (110 Da per residue) |
| `pf01328_*` | PF01328 GA hit, E-value, score, domain count, HMM coverage |
| `FamII_*`, `FamI_UPO_*`, `CPO_*` | full-sequence E-value and score for the three family HMMs |
| `best_family_hmm`, `family_margin` | best family HMM and bit-score lead over the second |
| `pcp_triad`, `prox_cys`, `prox_cys_qpos` | residues at AaeUPO P35-C36-P37; query position of the Cys |
| `sig_motif` | residues at AaeUPO S94-I95-G96 |
| `exd_motif`, `mg_glu`, `exd_mid`, `exd_asp`, `mg_ser` | metal site at AaeUPO E122-G123-D124 and S126 (EGD vs EHD) |
| `sxxrxd_R`, `sxxrxd_D` | AaeUPO R129, D131 |
| `cs_arg`, `acid_base_glu`, `acid_base_qpos` | AaeUPO R189 and E196 columns |
| `rxe_motif` | first R-X(6,7)-E in the query segment aligned to AaeUPO 137–219, with query position |
| `acid_base_ok`, `acid_base_basis` | acid-base Glu found, and how (`E196_column` or `RxE_window`) |
| `cterm_cys1`, `cterm_cys2`, `cterm_ss_pair` | AaeUPO C278 and C319 (C-terminal disulfide) |
| `core_cov` | fraction of AaeUPO 17–237 covered by the alignment |
| `length_class` | long if ≥ 330 aa |
| `signal_peptide` | yes / no / NA |
| `upobase_catalytic_label` | PCP-EGD-R-E, PCP-EHD-E or neither (from the alignment) |
| `upobase_subfamily`, `upobase_group`, `upobase_all_matches` | UPObase Table 1 class (Ia..V, Pog-a..g, none) |
| `phmmer_best_seed`, `phmmer_seed_evalue`, `phmmer_included` | best UPObase seed by phmmer; included if E ≤ 0.01 |
| `nearest_characterized`, `nearest_pct_identity`, `nearest_query_cov` | closest of 13 characterized enzymes |
| `kinner_family` | I or II (empty if motif checks fail) |
| `final_class` | `UPO_FamilyII`, `UPO_FamilyI`, `CPO_like_FamilyI`, `HTP_fails_motif`, `no_HTP_hit` |
| `confidence`, `reason` | high/medium/low and the text reason, including each discordant check |

## Files

| Path | Content |
|---|---|
| `upo_id/upo_id.py` | the tool (build + classify) |
| `scripts/fetch_refs.py` | reference download with back-off |
| `ref/ref_sources.tsv` | what to fetch and why |
| `ref/refs_all.fa`, `ref/manifest.tsv` | 92 reference sequences; manifest with URL, header, SHA-256, date |
| `ref/train_sets.tsv` | training group for each reference, with the basis |
| `models/` | anchor HMM, three family HMMs, PF01328 HMM, reference alignment, position map, seeds |
| `tests/refs_selfclass.tsv` | all 92 references classified with the full models |
| `tests/loo/loo_results.tsv` | leave-one-out results (70 references) |
| `tests/proteomes/` | four whole BFD proteomes (negative-control check) |
| `results/early807_*` | the 807 early-diverging PF01328 proteins: per-protein table and class counts |

## Known limitations

- `UPO_FamilyI` vs `CPO_like_FamilyI` rests on 8 characterized short UPOs against 16 uncharacterized NCBI "CPO"
  sequences. Kinner puts LfuCPO in Family I. Treat both labels as Family I.
- Family assignment uses HMM bit-score competition, not a tree. Scores from HMMs of different size are not strictly
  comparable. Check `family_margin` and `confidence`; prefer a phylogeny for any claim about a lineage.
- The R-X(6,7)-E window rule can accept a chance match.
- Only 35 of the 113 UPObase stage-1 seeds are public. The UPObase database itself is offline.
- The 330-aa and 20-bit cutoffs are calibrated on 79 references only.
- Speed (measured): 807 candidates in 2.5 min on the login node; whole proteomes of 10,000–28,000 proteins in
  under 1 min each. Not yet tested on all of BFD (97,520 PF01328 proteins); split large inputs across jobs.
