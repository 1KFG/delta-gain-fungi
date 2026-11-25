#!/bin/bash
# Joint SDR tree: make the run folder of one step (t4hn_reductase or t3hn_reductase) from the joint prepare folder
# results/step5/melanin_dhn/sdr_joint: share the tree files (symlinks), keep the domains of this step's candidates,
# and define the target clade of this step by bin/step5_multi_ref_clades.py (largest side with this step's references).
# usage: derive_sdr_step.sh <step> ; env: PC_ROOT (analysis/pathway_census). Runs on a compute node (1 CPU) or the login node.
set -euo pipefail
S=$1; R=${PC_ROOT:?}
J=$R/results/step5/melanin_dhn/sdr_joint; D=$R/results/step5/melanin_dhn/$S; T=$R/results/step5/other_inputs/tree_inputs/sdr_joint
mkdir -p $D
for f in ref.raxml.bestTree ref.raxml.bestModel ref_set.faa ref_set.trim.afa rep_ids.txt reps_summary.json reps_by_phylum.tsv; do ln -sf $J/$f $D/$f; done
cp -L --remove-destination $J/clusters.tsv.zst $D/clusters.tsv.zst   # zstd refuses a symlink as input
/usr/bin/python3.12 - $R/results/search/melanin_myco_full $S $J/domains.faa $D/domains.faa $J/domain_coords.tsv $D/domain_coords.tsv <<'PY'
import glob, subprocess, sys, csv
d, step, fa_in, fa_out, co_in, co_out = sys.argv[1:7]
ids = set()
for p in sorted(glob.glob(f"{d}/candidates.chunk*.PROVISIONAL.tsv.zst")):
    for l in subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout.decode().splitlines():
        f = l.split("\t")
        if not l.startswith("#") and f[-1] == step:
            ids.add(f[0])
keep, n = False, 0
with open(fa_in) as fi, open(fa_out, "w") as fo:
    for line in fi:
        if line.startswith(">"):
            keep = line[1:].split()[0] in ids
            n += keep
        if keep:
            fo.write(line)
with open(co_in) as fi, open(co_out, "w") as fo:
    for i, line in enumerate(fi):
        if i == 0 or line.split("\t")[0] in ids:
            fo.write(line)
print(f"{step}: {len(ids)} candidate proteins, {n} domains kept")
PY
export PC_PATHWAY=melanin_dhn PC_STEP=$S PC_OUT=$D PC_DECOYS_FA=$T/decoys.faa PC_OUTGROUP_FA=$T/outgroups.faa PC_OUTGROUP_ROLE=tree_outgroup
if [ "$S" = t4hn_reductase ]; then export PC_MEMBERS_FA=$T/members.faa; else export PC_MEMBERS_FA=$T/members_for_t3hn.faa; fi
/usr/bin/python3.12 $R/bin/step5_multi_ref_clades.py --tree $D/ref.raxml.bestTree --refs $R/refs/melanin_dhn/$S.faa --decoys $PC_DECOYS_FA \
  --outgroup $PC_OUTGROUP_FA --members-fa $PC_MEMBERS_FA --outgroup-role tree_outgroup --reps $D/rep_ids.txt --step $S --out-dir $D
