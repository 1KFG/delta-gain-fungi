# Sourced by the upo_place SLURM scripts. Run from the repo root (sbatch's working directory).
UPO=analysis/enzymes_UPO
TREEDIR=$UPO/data/refset_zamocky2022/tree
PREF=$UPO/data/refset_zamocky2022/placement_ref   # small, committed: ref.hmm, ref.raxml.*
PY=/usr/bin/python3.12
TOOL=$UPO/upo_place/bin/upo_place.py
NCOL=224                                          # columns of tree/aln.trim.afa
# Placement model. ModelFinder chose Q.PFAM+R7 for the tree, but epa-ng 0.3.8 rejects Q.PFAM
# ("Invalid model name"). LG+G4 on the same (mfp) topology: clade membership was identical under
# LG+I+G4 (lgig run), so the clades being placed into do not depend on the matrix.
MODEL=${MODEL:-LG+G4}
# Pinned bioconda builds through pixi (no module exists for epa-ng / gappa).
PX="pixi exec -c conda-forge -c bioconda -s epa-ng==0.3.8 -s raxml-ng==2.0.4 --"
module load hmmer/3.4
