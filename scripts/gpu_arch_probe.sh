#!/bin/bash
# Does the new env actually work on each GPU type this cluster has?
#   for g in a100 ada6000 h100 blackwell6000; do
#     sbatch --gres=gpu:$g:1 --export=ALL,GPUTYPE=$g scripts/gpu_arch_probe.sh
#   done
#
# The old esm/2.0.0 torch (1.13.1+cu117, sm_37..sm_86) does not refuse to run on a card it was
# not built for - it warns and returns garbage. So "it ran" is not evidence. This probe checks
# three things in order, and the third is the only one that matters:
#   1. torch.cuda.get_arch_list() contains this card's sm_XX (or a PTX path that JITs).
#   2. a matmul gives the right answer, checked against a CPU reference, not just "no crash".
#   3. ESMFold folds a real sequence and the pLDDT is plausible (a garbage fold scores ~0.2-0.3
#      and has no secondary structure; the test peptide should land well above that).
# Writes nothing outside $SCRATCH.
#SBATCH -p short_gpu
#SBATCH -N 1 -n 1 -c 8
#SBATCH --mem=48G
#SBATCH -t 0:30:00
#SBATCH -J gpu_probe
#SBATCH -o logs/gpu_probe_%j.out
#SBATCH -e logs/gpu_probe_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

ENVDIR="${ESMFOLD_ENV:-/bigdata/stajichlab/jstajich/envs/esmfold-torch28}"
CACHE="${ESMFOLD_CACHE:-/bigdata/stajichlab/jstajich/cache/huggingface}"

export CONDA_NO_PLUGINS=true
source /etc/profile.d/modules.sh 2>/dev/null || true
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate "$ENVDIR"
export HF_HOME="$CACHE"
export TOKENIZERS_PARALLELISM=false

echo "requested GPU type: ${GPUTYPE:-unspecified}"
nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader || true

python - <<'PYEOF'
import sys, time
import torch

print("python", sys.version.split()[0])
print("torch", torch.__version__, "cuda", torch.version.cuda)
if not torch.cuda.is_available():
    sys.exit("FAIL: no CUDA device visible")

cap = torch.cuda.get_device_capability()
sm = f"sm_{cap[0]}{cap[1]}"
archs = torch.cuda.get_arch_list()
print("device", torch.cuda.get_device_name(0), sm)
print("built_for", archs)
print("arch_match", "YES" if sm in archs else f"NO - relies on PTX JIT from {archs[-1] if archs else 'nothing'}")

# 2. correctness, not just absence of a crash
torch.manual_seed(0)
a = torch.randn(1024, 1024, dtype=torch.float32)
b = torch.randn(1024, 1024, dtype=torch.float32)
ref = a @ b
got = (a.cuda() @ b.cuda()).cpu()
err = (ref - got).abs().max().item()
print(f"matmul max_abs_err {err:.3e}", "OK" if err < 1e-2 else "FAIL")
if err >= 1e-2:
    sys.exit("FAIL: GPU matmul disagrees with CPU - this is the silent-garbage mode")

# 3. a real fold. Ubiquitin (1UBQ, 76 aa) folds confidently in any working ESMFold.
from transformers import AutoTokenizer, EsmForProteinFolding
UBQ = ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")

t0 = time.time()
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
model = model.cuda().eval()
model.esm = model.esm.half()          # the ESM-2 trunk is fp16 on GPU
torch.backends.cuda.matmul.allow_tf32 = True
print(f"model loaded in {time.time()-t0:.0f}s")

t0 = time.time()
enc = tok([UBQ], return_tensors="pt", add_special_tokens=False)
with torch.no_grad():
    out = model(enc["input_ids"].cuda())
plddt = out["plddt"][0, :, 1].mean().item()
dt = time.time() - t0
print(f"ubiquitin pLDDT {plddt:.3f} in {dt:.1f}s")
# A correct ESMFold gives ubiquitin ~0.9. Garbage output sits near 0.2-0.3.
print("fold_sane", "YES" if plddt > 0.7 else "NO - output looks like the silent-garbage mode")
if plddt <= 0.7:
    sys.exit("FAIL: implausible pLDDT for ubiquitin")
print("ALL CHECKS PASSED")
PYEOF
echo "done $(date +%T)"
