"""Golden test: AC levitation equilibrium of a Cu sphere above a coil.

Puts docs/maglev/demos/sphere/coil_maglev_equilibrium.py
under the validation lane.  The script composes Radia's open-boundary coil field with the
verified induced-dipole levitation force and finds the stable equilibrium
height.  Its own hard asserts are the primary lock (a non-zero exit means
one failed); this test re-runs it and double-checks the headline numbers
from the committed JSON.

Needs radia (the package under test); skipped cleanly if it is absent.
"""
import json
import os
import shutil
import subprocess
import sys

import pytest

pytest.importorskip("radia")
pytestmark = pytest.mark.slow

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
_DEMOS = os.path.join(_ROOT, "docs", "maglev", "demos")
_SPHERE_DEMO = os.path.join(_DEMOS, "sphere")
_SCRIPT = os.path.join(_SPHERE_DEMO, "coil_maglev_equilibrium.py")
_JSON_NAME = "coil_maglev_equilibrium_results.json"


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    # the script's internal asserts are the real verification; exit 0 = all pass
    run_dir = tmp_path_factory.mktemp("coil_maglev_equilibrium")
    shutil.copytree(
        _SPHERE_DEMO,
        run_dir,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.json", "*.png"),
    )
    script = os.path.join(run_dir, os.path.basename(_SCRIPT))
    result_json = os.path.join(run_dir, _JSON_NAME)
    env = os.environ.copy()
    env["PYTHONPATH"] = _DEMOS + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, script],
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=300,
        env=env,
    )
    assert proc.returncode == 0, (
        f"coil_maglev_equilibrium.py failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    with open(result_json) as f:
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
