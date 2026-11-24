"""Run tests/gainloss_null_models_toy_check.R (model-aware simulation and fit). Skipped when Rscript cannot load ape and castor."""
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent / "gainloss_null_models_toy_check.R"


def _rscript():
    r = shutil.which("Rscript")
    if r is None:
        return None
    p = subprocess.run([r, "-e", "library(ape); library(castor)"], capture_output=True, text=True)
    return r if p.returncode == 0 else None


RSCRIPT = _rscript()
pytestmark = pytest.mark.skipif(RSCRIPT is None, reason="Rscript with ape and castor not available")


@pytest.fixture(scope="module")
def toy():
    out = subprocess.run([RSCRIPT, str(SCRIPT)], capture_output=True, text=True, check=True).stdout
    return dict(l.split("=", 1) for l in out.splitlines() if "=" in l)


def test_fit_params_by_aic_columns_and_choice(toy):
    assert int(toy["params_rows"]) == 3 and int(toy["params_cols_ok"]) == 1
    assert int(toy["model_matches_aic"]) == 1 and int(toy["n_present_ok"]) == 1 and int(toy["er_equal_rates"]) == 1


def test_model_aware_flags_and_degenerate(toy):
    assert int(toy["flags_len"]) == 2 * 120 - 2 and int(toy["flags_same_len"]) == 1 and int(toy["const_degenerate"]) == 1


def test_replicate_returns_all_pairs_and_types(toy):
    assert int(toy["rep_rows"]) == 12 and int(toy["rep_steps"]) == 3 and int(toy["matched_rows"]) == 12
