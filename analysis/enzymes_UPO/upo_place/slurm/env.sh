# Sourced by the upo_place SLURM scripts. Run from the repo root (sbatch's working directory).
UPO=analysis/enzymes_UPO
TREEDIR=$UPO/data/refset_zamocky2022/tree
PREF=$UPO/data/refset_zamocky2022/placement_ref   # small, committed: ref.hmm, ref.raxml.*
PY=/usr/bin/python3.12
TOOL=$UPO/upo_place/bin/upo_place.py
NCOL=224                                          # columns of tree/aln.trim.afa
MODEL=${MODEL:-Q.PFAM+R7}                         # ModelFinder's choice for mfp
# Pinned bioconda builds through pixi (no module exists for epa-ng / gappa).
PX="pixi exec -c conda-forge -c bioconda -s epa-ng==0.3.8 -s raxml-ng==2.0.4 --"
module load hmmer/3.4
