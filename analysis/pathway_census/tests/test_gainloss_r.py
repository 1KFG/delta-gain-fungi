"""Run tests/gainloss_toy_check.R (toy tree with known gains, losses and shared gains) and assert the values.

Skipped when Rscript cannot load ape and castor. On HPCC run `module unload R; module load R/4.2.0` first.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent / "gainloss_toy_check.R"


def _rscript():
    r = shutil.which("Rscript")
    if r is None:
        return None
    p = subprocess.run([r, "-e", "library(ape); library(castor)"], capture_output=True, text=True)
    return r if p.returncode == 0 else None


RSCRIPT = _rscript()


@pytest.fixture(scope="module")
def toy():
    out = subprocess.run([RSCRIPT, str(SCRIPT)], capture_output=True, text=True, check=True).stdout
    return dict(l.split("=", 1) for l in out.splitlines() if "=" in l)


pytestmark = pytest.mark.skipif(RSCRIPT is None, reason="Rscript with ape and castor not available")


def test_marginals_match_castor(toy):
    assert float(toy["marginal_max_abs_diff_vs_castor"]) < 1e-3


@pytest.mark.parametrize("key,val", [("ntip", 48), ("A_gains", 2), ("A_losses", 0), ("AB_gains", 3), ("AB_losses", 0),
                                      ("LB_gains", 0), ("LB_losses", 2)])
def test_known_gain_and_loss_counts(toy, key, val):
    assert int(toy[key]) == val


@pytest.mark.parametrize("key,val", [("AB_pars_changes", 2), ("AB_pars_changes_castor", 2), ("AB_pars_gains_late", 2),
                                      ("AB_pars_gains_early", 2), ("AB_pars_losses", 0), ("LB_pars_changes", 1),
                                      ("LB_pars_losses_late", 1), ("LB_pars_gains", 0), ("T2_pars_changes", 2),
                                      ("T2_sum_late", 2), ("T2_sum_early", 2)])
def test_parsimony_counts(toy, key, val):
    assert int(toy[key]) == val


def test_cogain_counts_and_permutation(toy):
    assert int(toy["U_V_cogain"]) == 3 and int(toy["U_W_cogain"]) == 2
    assert int(toy["perm_observed"]) == 3 and int(toy["perm_expected_below_observed"]) == 1
    assert int(toy["perm_p_below_0.05"]) == 1
    assert int(toy["perm_p_na_when_no_events"]) == 1


def test_prune_drops_na_tips(toy):
    assert int(toy["pruned_ntip"]) == 46
