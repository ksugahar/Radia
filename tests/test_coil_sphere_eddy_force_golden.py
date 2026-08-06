"""Golden test: full eddy-current Lorentz force on a sphere-in-coil vs the
induced-dipole force -- the point-dipole approximation error at a/L ~ 0.5.

Locks docs/maglev/demos/sphere/coil_sphere_eddy_force.py:
  - the FULL axisymmetric eddy-current force and the induced-dipole force
    are the same order at a/L ~ 0.5 (ratio near 1, but not identical), and
  - shrinking the sphere (a/L 0.5 -> 0.1) drives the ratio TOWARD 1 -- the
    signature that the gap is genuinely the point-dipole approximation.

Runs four real axisymmetric NGSolve solves (~1-2 min); skipped cleanly if
ngsolve / netgen are unavailable.
"""
import json
import os
import shutil
import subprocess
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPHERE_DEMO = os.path.join(_HERE, "..", "docs", "maglev", "demos", "sphere")
_SCRIPT = os.path.join(_SPHERE_DEMO, "coil_sphere_eddy_force.py")
_JSON = os.path.join(_SPHERE_DEMO, "coil_sphere_eddy_force_results.json")


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    if not os.path.isfile(_SCRIPT):
        pytest.skip(f"coil sphere demo script not found: {_SCRIPT}")
    run_dir = tmp_path_factory.mktemp("coil_sphere_eddy_force")
    shutil.copytree(
        _SPHERE_DEMO, run_dir, dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.json", "*.png"))
    script = os.path.join(run_dir, os.path.basename(_SCRIPT))
    result_json = os.path.join(run_dir, os.path.basename(_JSON))
    proc = subprocess.run([sys.executable, script], cwd=run_dir,
                          capture_output=True, text=True, timeout=600,
                          check=False)
    assert proc.returncode == 0, (
        f"coil_sphere_eddy_force.py failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    with open(result_json) as f:
        return json.load(f)


def test_full_and_dipole_same_order_at_aL_half(results):
    r5 = results["case_5mm"]
    assert abs(r5["a_over_L"] - 0.5) < 0.1
    # the full eddy force and the dipole force agree to within a few % at
    # a/L ~ 0.5 (close, but the nonzero gap is the dipole-approximation error)
    assert 0.95 < r5["ratio_full_over_dipole"] < 1.10
    assert r5["F_full_N"] > 0 and r5["F_dipole_N"] > 0   # both lift


def test_dipole_becomes_exact_as_sphere_shrinks(results):
    r5, r1 = results["case_5mm"], results["case_1mm"]
    dev5 = abs(r5["ratio_full_over_dipole"] - 1.0)
    dev1 = abs(r1["ratio_full_over_dipole"] - 1.0)
    # a/L 0.5 -> 0.1 must move the ratio toward 1 (the dipole-approx signature)
    assert r1["a_over_L"] < r5["a_over_L"]
    assert dev1 < dev5
    assert dev1 < 0.01      # 1 mm sphere: dipole essentially exact


def test_coil_field_matches_analytic_loop(results):
    assert results["case_5mm"]["coil_field_max_err_percent"] < 5.0
