#!/bin/bash
# Can AlphaFold3 use a card other than the a100?
#   for g in blackwell6000 ada6000 h100; do
#     sbatch --gres=gpu:$g:1 --export=ALL,GPUTYPE=$g -J af3probe_$g scripts/af3_gpu_probe.sh
#   done
#
# af3_fold.sh asks for --gres=gpu:a100:1. That was inherited, never tested. AF3 is a singularity
# container with its own jax (0.9.1), entirely separate from the ESMFold env in
# scripts/setup_gpu_env.sh, so what we proved about torch 2.8 says nothing about it.
#
# jax fails more honestly than torch did - an unsupported card usually means no CUDA device or a
# PTX error rather than silent garbage - but "usually" is not a basis for a production run, so
# this checks a computed result against a known answer instead of just looking for a crash.
# Writes nothing outside $SCRATCH.
#SBATCH -p short_gpu
#SBATCH -N 1 -n 1 -c 4
#SBATCH --mem=32G
#SBATCH -t 0:20:00
#SBATCH -J af3_gpu_probe
#SBATCH -o logs/af3_gpu_probe_%j.out
#SBATCH -e logs/af3_gpu_probe_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs
module load alphafold/3.0.2

echo "requested GPU type: ${GPUTYPE:-unspecified}"
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader || true

# --nv passes the host driver into the container; without it jax sees no GPU at all and would
# report a false negative here.
singularity exec --nv "$ALPHAFOLD_SING" python3 - <<'PYEOF'
import sys
import jax, jaxlib
import jax.numpy as jnp

print("jax", jax.__version__, "jaxlib", jaxlib.__version__)
devs = jax.devices()
print("devices", devs)
if not any(d.platform == "gpu" for d in devs):
    sys.exit("FAIL: jax sees no GPU - AF3 would fall back to CPU or refuse")
d = devs[0]
print("device_kind", getattr(d, "device_kind", "?"))

# A result with a known answer, not just "it ran": ones(n,n) @ ones(n,n) sums to n**3 exactly.
n = 1024
x = jnp.ones((n, n), dtype=jnp.float32)
got = float((x @ x).sum())
want = float(n) ** 3
rel = abs(got - want) / want
print(f"matmul got {got:.0f} want {want:.0f} rel_err {rel:.2e}", "OK" if rel < 1e-5 else "FAIL")
if rel >= 1e-5:
    sys.exit("FAIL: jax matmul is wrong on this card")
print("ALL CHECKS PASSED")
PYEOF
echo "done $(date +%T)"
