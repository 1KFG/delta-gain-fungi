"""Run tests/gainloss_null_toy_check.R and assert the simulation helpers. Skipped when Rscript cannot load ape and castor."""
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent / "gainloss_null_toy_check.R"


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


def test_simulation_draws_right_number_of_named_tips(toy):
    assert int(toy["ntip"]) == 48 and int(toy["sim_len"]) == 48 and int(toy["sim_n"]) == 5 and int(toy["sim_names_ok"]) == 1


def test_simulation_is_reproducible_and_respects_zero_rates(toy):
    assert int(toy["same_seed_identical"]) == 1
    assert int(toy["zero_rate_all_present"]) == 1 and int(toy["zero_rate_all_absent"]) == 1


def test_fixed_seed_simulation_regression(toy):
    assert int(toy["seed1_sum"]) == SEED1_SUM


def test_known_event_counts_and_internal_scope(toy):
    assert int(toy["U_gains"]) == 3 and int(toy["W_gains"]) == 2 and int(toy["U_degenerate"]) == 0
    assert int(toy["UW_cogain_all"]) == 2 and int(toy["UW_cogain_internal"]) == 1
    assert int(toy["const_degenerate"]) == 1 and int(toy["const_events"]) == 0


def test_empirical_p_ties_and_zero_exceedances(toy):
    assert float(toy["p_zero_exceed"]) == pytest.approx(0.2)
    assert float(toy["p_ties"]) == pytest.approx(0.8)
    assert float(toy["p_all_exceed"]) == pytest.approx(1.0)


def test_replicate_pipeline_returns_all_pair_types(toy):
    assert int(toy["rep_rows"]) == 4 and int(toy["rep_steps"]) == 2


def test_prevalence_matched_simulation_in_window(toy):
    assert int(toy["matched_in_window"]) == 1 and int(toy["matched_len"]) == 48


SEED1_SUM = 8  # measured once with castor 1.7.10, set.seed(1); a regression value
