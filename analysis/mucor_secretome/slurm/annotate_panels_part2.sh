#!/bin/bash
# The three CPU predictors that annotate_panels.sh never ran: WoLF PSORT, MEROPS and
# AIUPred. DeepLoc is the fourth and needs a GPU, so it is in slurm/deeploc_panels.sh.
#   sbatch analysis/mucor_secretome/slurm/annotate_panels_part2.sh
#
# Why this exists
# ---------------
# annotate_panels.sh's header states the contract these panels are built on: every predictor
# is invoked as the matching Fungi_BFD/nextflow/modules/BFD/* process invokes it, "so these
# three proteomes are annotated by the same tools, versions and flags as every BFD genome and
# build_secretome.py can classify them by one set of rules".
#
# Six predictors honoured that. Four did not exist for these proteomes at all, while existing
# for all ~23,500 BFD genomes: deeploc, wolfpsort, merops, idp_summary. Any comparison using
# them was therefore impossible, and substituting homologs from related BFD genomes would have
# been worse than useless: candidate effectors are small, fast-evolving and often
# lineage-specific, so they are exactly the proteins with no confident homolog. The missingness
# would have correlated with the outcome.
#
# This is a SEPARATE script rather than an edit to annotate_panels.sh on purpose:
# annotate_to_tables.py APPENDS, so re-running that script would duplicate every existing row.
#
# Invocations, each copied from the named BFD module (nextflow/modules/BFD/<NAME>/main.nf):
#   WOLFPSORT  cat X | runWolfPsortSummary fungi          (module wolfpsort)
#   MEROPS     blastp -db $MEROPS_DB/merops_scan.lib -seg yes -soft_masking true
#              -max_target_seqs 10 -evalue 1e-10 -outfmt 6 -use_sw_tback,
#              run from params.blastp_sif (blast 2.16.0) with db-merops/124
#   IDP        aiupred.py -i X -o X.aiupred.txt           (module aiupred)
# Tables are then built by BFD's OWN merge scripts (nextflow/bin/merge_*.py,
# gather_AIUPred.py) rather than by a reimplementation here, so the parsing is identical by
# construction and the column layouts match tables/{wolfpsort,merops,idp_summary}.parquet.
#
# One deliberate departure, and the reason for it:
#   WoLF PSORT takes an organism type. BFD runs `fungi` because BFD is all fungi. The two
#   Marchantia proteomes are plants, so they are run with `plant`. Using `fungi` on a plant
#   proteome would be wrong biology for the sake of a false equivalence -- and the comparison
#   that matters is Mucor against BFD's fungi, which is unaffected. The per-genome organism
#   type is recorded in out/panels/wolfpsort_organism.tsv so nobody has to guess later.
#
# Reads shared folders read-only; writes analysis/mucor_secretome/out/panels/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J annot_part2
#SBATCH -o logs/annot_part2_%j.out
#SBATCH -e logs/annot_part2_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

A=analysis/mucor_secretome
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-16}
BFD_BIN=${BFD_BIN:-$HOME/projects/Fungi_BFD/nextflow/bin}
OUT=$A/out/panels; mkdir -p $OUT
FAA=$OUT/faa
W="${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}/annot_part2"
mkdir -p "$W"

source /etc/profile.d/modules.sh 2>/dev/null || true
export CONDA_NO_PLUGINS=true
# aiupred and wolfpsort are conda environments behind their modulefiles, so their `conda
# activate` needs the shell hook first. A login shell has it from ~/.bashrc, which is why
# these modules load fine interactively; a batch shell does not, and the modulefile dies with
#   conda: error: unrecognized arguments: --stack
# 20 s in. annotate_panels.sh carries the same preamble for the same reason.
module load miniconda3
eval "$(conda shell.bash hook)"
module load wolfpsort aiupred db-merops/124 apptainer
BLASTP_SIF=${BLASTP_SIF:-docker://quay.io/biocontainers/blast:2.16.0--h66d330f_5}

# Fail here, with a message naming what is missing, rather than part-way through the first
# proteome. Each of these comes from a module loaded above.
for c in runWolfPsortSummary aiupred.py apptainer; do
    command -v "$c" >/dev/null || { echo "[ERROR] $c not on PATH after module load" >&2; exit 1; }
done
[ -n "${MEROPS_DB:-}" ] || { echo "[ERROR] MEROPS_DB unset (db-merops module)" >&2; exit 1; }
[ -s "$MEROPS_DB/merops_scan.lib.phr" ] || [ -s "$MEROPS_DB/merops_scan.lib.pin" ] || \
    echo "[warn] $MEROPS_DB/merops_scan.lib blast index not obviously present; continuing"

LOCUSTAGS="${PANEL_LOCUSTAGS:-JEM1MUCOR JEM2MUCOR MPTAK MPALE}"

# The proteomes are staged by annotate_panels.sh; this script does not re-derive them, so a
# missing .faa means the first script has not run and the two would go out of step.
for LT in $LOCUSTAGS; do
  [ -s "$FAA/$LT.faa" ] || { echo "[ERROR] $FAA/$LT.faa missing - run annotate_panels.sh first" >&2; exit 1; }
done

echo "locustag	organism" > $OUT/wolfpsort_organism.tsv

for LT in $LOCUSTAGS; do
  F=$FAA/$LT.faa
  N=$(grep -c '^>' "$F")
  echo "=== $LT  $N proteins  $(date +%T)"

  case "$LT" in
    MPTAK|MPALE) ORG=plant ;;
    *)           ORG=fungi ;;
  esac
  echo "$LT	$ORG" >> $OUT/wolfpsort_organism.tsv

  # --- WoLF PSORT (module WOLFPSORT) -------------------------------------------------
  runWolfPsortSummary "$ORG" < "$F" > "$W/$LT.wolfpsort.results.txt"
  gzip -f "$W/$LT.wolfpsort.results.txt"

  # --- MEROPS (module MEROPS) --------------------------------------------------------
  export TMPDIR="$W"
  PROT_REAL_DIR=$(dirname "$(readlink -f "$F")")
  apptainer exec \
      --bind "$PWD:$PWD,$MEROPS_DB:$MEROPS_DB,$PROT_REAL_DIR:$PROT_REAL_DIR,$TMPDIR:$TMPDIR" \
      "$BLASTP_SIF" \
      blastp -query "$F" \
             -db "$MEROPS_DB/merops_scan.lib" \
             -out "$W/$LT.blasttab" \
             -num_threads "$CPU" \
             -seg yes -soft_masking true \
             -max_target_seqs 10 \
             -evalue 1e-10 \
             -outfmt 6 \
             -use_sw_tback
  gzip -f "$W/$LT.blasttab"

  # --- AIUPred (module IDP) ----------------------------------------------------------
  aiupred.py -i "$F" -o "$W/$LT.aiupred.txt"
  [ -f "$W/$LT.aiupred.txt" ] && gzip -f "$W/$LT.aiupred.txt"
done

# --- tables, built by BFD's own merge scripts --------------------------------------
echo "=== merging $(date +%T)"
$PY "$BFD_BIN/merge_wolfpsort.py" $W/*.wolfpsort.results.txt.gz -o "$OUT/wolfpsort.csv"
$PY "$BFD_BIN/merge_merops.py"    $W/*.blasttab.gz              -o "$OUT/merops.csv"
$PY "$BFD_BIN/gather_AIUPred.py"  $W/*.aiupred.txt.gz \
        --outfile "$OUT/idp.csv" --outfilesum "$OUT/idp_summary.csv"

for f in wolfpsort merops idp_summary; do
  [ -s "$OUT/$f.csv" ] || { echo "[ERROR] $OUT/$f.csv is empty" >&2; exit 1; }
  echo "$f.csv  $(( $(wc -l < "$OUT/$f.csv") - 1 )) rows"
done

echo "done $(date +%T)"
