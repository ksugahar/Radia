"""Golden test: AC levitation equilibrium of a Cu sphere above a coil.

Puts examples/levitation/sphere/coil_levitation_equilibrium.py
under CI.  The script composes Radia's open-boundary coil field with the
verified induced-dipole levitation force and finds the stable equilibrium
height.  Its own hard asserts are the primary lock (a non-zero exit means
one failed); this test re-runs it and double-checks the headline numbers
from the committed JSON.

Needs radia (the package under test); skipped cleanly if it is absent.
"""
import json
import os
import subprocess
import sys

import pytest

pytest.importorskip("radia")

_HERE = os.path.dirname(os.path.abspath(__file__))
_TEAM28 = os.path.join(_HERE, "..", "examples", "levitation", "sphere")
_SCRIPT = os.path.join(_TEAM28, "coil_levitation_equilibrium.py")
_JSON = os.path.join(_TEAM28, "coil_levitation_equilibrium_results.json")


@pytest.fixture(scope="module")
def results():
    # the script's internal asserts are the real verification; exit 0 = all pass
    proc = subprocess.run([sys.executable, _SCRIPT], cwd=_TEAM28,
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, (
        f"coil_levitation_equilibrium.py failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    with open(_JSON) as f:
        return json.load(f)


def test_coil_field_matches_analytic_loop(results):
    assert results["field_check"]["max_rel_error_percent"] < 0.1


def test_sphere_responds_diamagnetically(results):
    assert results["frequency"]["a_over_delta"] > 2.0
    assert results["sphere_response"]["Re_G"] < 0.0


def test_stable_equilibrium_exists(results):
    eq = results["equilibrium"]
    # a real equilibrium in the expected vicinity, lift balances weight
    assert 25.0 < eq["z_star_mm"] < 45.0
    assert abs(eq["residual_percent"]) < 2.0
    assert eq["peak_over_weight"] > 1.0


def test_equilibrium_is_vertically_stable(results):
    st = results["stability"]
    assert st["stable"] is True
    assert st["dFdz_N_per_m"] < 0.0
