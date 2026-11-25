#!/bin/bash
# Clusters at 0.95 identity (mmseqs2 easy-cluster, -c 0.8 --cov-mode 0) for the two cells marked "not run" in other_steps_RUN.md:
# hppd and the SDR union (t4hn + t3hn). Run from analysis/pathway_census inside a job or with at most 2 threads. Needs module mmseqs2/17-b804f.
set -euo pipefail
W=${SCRATCH:?}/cc95; mkdir -p $W
for spec in hppd=results/step5/melanin_ldopa/hppd/domains.faa sdr_union=results/step5/melanin_dhn/sdr_joint/domains.faa; do
  n=${spec%%=*}; f=${spec#*=}
  for id in 0.9 0.95; do
    mmseqs easy-cluster $f $W/$n.$id $W/tmp --min-seq-id $id -c 0.8 --cov-mode 0 --threads 2 -v 1 > /dev/null
    echo -e "$n\t$id\t$(cut -f1 $W/${n}.${id}_cluster.tsv | sort -u | wc -l)"
    rm -rf $W/tmp
  done
done
