#!/bin/bash
# Build the project's own GPU python env, so structure prediction is not limited to a100.
#   bash scripts/setup_gpu_env.sh            # create or update
#   ESMFOLD_ENV=/some/other/prefix bash ...  # override the location
#
# Why this exists
# ---------------
# The cluster's `esm/2.0.0` module is a self-contained **Python 3.7** env
# (/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env) carrying torch 1.13.1+cu117, built for
# sm_37..sm_86. That excludes every GPU on this cluster except the a100 (sm_80): ada6000 is
# sm_89, h100 sm_90, blackwell6000 sm_120. On those cards torch 1.13 does not refuse - it warns
# and returns garbage (job 32033882: 487 errors and 72 worthless "successes"). Python 3.7 is also
# why the module cannot simply be upgraded in place; modern torch needs 3.9+.
#
# Counting GPUs across short_gpu / preempt_gpu / gpu / exfab: 24 a100 usable, and 42
# blackwell6000 + 16 ada6000 + 2 h100 = 60 locked out. This env unlocks those.
#
# The other half of the fix is in the code, not here: `esm.pretrained.esmfold_v1()` (fair-esm)
# needs the `openfold` package and its compiled CUDA kernels, which is what pins the old module
# together. HuggingFace `transformers.EsmForProteinFolding` vendors the openfold structure module
# in pure torch, so it needs no openfold install and follows whatever torch it is given. DeltaGain's
# own structure pilot already ran `facebook/esmfold_v1` this way on an RTX 6000 Ada (exfab),
# ~8-14 min per 400 proteins, which is the existence proof that this path works here.
#
# Not bit-identical to fair-esm. Models produced before and after the switch should be regenerated
# or kept separately, never silently mixed.
#
# Writes only under /bigdata/stajichlab/jstajich (never the shared group folder).
set -euo pipefail

ENVDIR="${ESMFOLD_ENV:-/bigdata/stajichlab/jstajich/envs/esmfold-torch28}"
CACHE="${ESMFOLD_CACHE:-/bigdata/stajichlab/jstajich/cache/huggingface}"

# Pins. cu128 is what reaches sm_120 (blackwell6000); torch 2.8.0 is the version already proven
# present and working in the miniconda3 base env on this cluster. transformers 4.57.6 likewise.
PY_VER=3.11
TORCH_VER=2.8.0
TORCH_INDEX=https://download.pytorch.org/whl/cu128
TRANSFORMERS_VER=4.57.6

echo "env    : $ENVDIR"
echo "cache  : $CACHE"
echo "pins   : python $PY_VER, torch $TORCH_VER (cu128), transformers $TRANSFORMERS_VER"

export CONDA_NO_PLUGINS=true
source /etc/profile.d/modules.sh 2>/dev/null || true
module load miniconda3
eval "$(conda shell.bash hook)"

# The miniconda3 module's package cache lives under /opt and is read-only, so mamba cannot take
# its lockfile ("Could not open lockfile .../pkgs/cache/cache.lock") and the solve stalls.
# Point the package cache somewhere we own.
export CONDA_PKGS_DIRS="${ESMFOLD_PKGS:-/bigdata/stajichlab/jstajich/cache/conda_pkgs}"
mkdir -p "$CONDA_PKGS_DIRS"

mkdir -p "$(dirname "$ENVDIR")" "$CACHE"

if [ ! -x "$ENVDIR/bin/python" ]; then
  echo "=== creating env ==="
  mamba create -y -p "$ENVDIR" "python=$PY_VER" pip
else
  echo "=== env exists, updating packages ==="
fi

conda activate "$ENVDIR"
python -V

echo "=== torch $TORCH_VER (cu128) ==="
# --index-url, not --extra-index-url: the pytorch index must win for the torch wheels.
pip install --no-cache-dir --index-url "$TORCH_INDEX" "torch==$TORCH_VER"

echo "=== transformers and friends ==="
# ESMFold's vendored openfold utils need scipy; biopython is for reading/writing structures.
pip install --no-cache-dir \
    "transformers==$TRANSFORMERS_VER" \
    "accelerate" \
    "scipy" \
    "biopython"

echo "=== pinned versions actually installed ==="
pip list 2>/dev/null | grep -E '^(torch|transformers|accelerate|scipy|biopython) ' || true

echo "=== import check (CPU; the GPU arch check belongs in a job) ==="
python - <<'PYEOF'
import torch, transformers
print("torch", torch.__version__, "cuda", torch.version.cuda)
print("transformers", transformers.__version__)
from transformers import AutoTokenizer, EsmForProteinFolding
print("EsmForProteinFolding import OK (no openfold package needed)")
PYEOF

echo "=== prefetching facebook/esmfold_v1 weights into the cache ==="
# Compute nodes may have no outbound network, so the weights are fetched here, on the login node,
# and read from $HF_HOME in the job. ~2.6 GB; skipped silently if already present.
HF_HOME="$CACHE" python - <<'PYEOF'
import os
from huggingface_hub import snapshot_download
p = snapshot_download("facebook/esmfold_v1")
print("weights at", p)
PYEOF

cat <<EOF

=== done ===
Use it in a job with:

    export CONDA_NO_PLUGINS=true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    conda activate $ENVDIR
    export HF_HOME=$CACHE

and request any GPU, not only a100.
EOF
