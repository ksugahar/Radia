"""Lock the radia-ih thermal "workpiece-only" contract (2026-05-31).

The Thermal step targets the WORKPIECE SOLID only.  ``calc_heat.py``
(and its axisym twin) MUST REJECT a multi-material (coil+workpiece)
mesh with a clear error rather than silently heating the coil region as
if it were steel (keiko 2026-05-31: a ``coil_wrkpiace_vol`` mesh was
passed to the thermal step).

These are SUBPROCESS tests: the pytest process never imports ngsolve /
PySide6 -- both the mesh generation and the calc run happen in child
processes.  That keeps them safe on the LAB box where importing PySide6
under pytest crashes (0xc0000139, MKL DLL shadow).  ``validation_test/panels`` is
CI-ignored, so this is a LAB/local regression lock (same tier as
``test_heat_chain_golden.py``).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURE_DIR = os.path.join(HERE, "fixtures")

WP_VOL = os.path.join(FIXTURE_DIR, "heat_workpiece_cylinder_R25_H25.vol")
GEN_WP = os.path.join(FIXTURE_DIR, "generate_heat_cylinder.py")
MULTI_VOL = os.path.join(FIXTURE_DIR, "heat_multimaterial_coil_wp.vol")
GEN_MULTI = os.path.join(FIXTURE_DIR, "generate_heat_multimaterial.py")

CALC_HEAT = os.path.join(ROOT, "src", "radia", "panels", "calc_heat.py")


def _ensure(gen_script, out_vol):
    """Return ``out_vol``, generating it via the OCC-only helper
    subprocess when missing.  Skips the test if generation fails (no
    ngsolve / netgen in this env)."""
    if os.path.isfile(out_vol):
        return out_vol
    proc = subprocess.run([sys.executable, gen_script],
                          capture_output=True, text=True, timeout=180)
    if proc.returncode != 0 or not os.path.isfile(out_vol):
        pytest.skip(
            f"could not generate {os.path.basename(out_vol)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return out_vol


def _run_heat(wp_vol, tmp_path):
    """Run calc_heat.py with a uniform q_surf (no .sol needed) and
    return the parsed JSON result.  calc_main always exits 0 and prints
    the result dict (errors included) as the last stdout line."""
    cmd = [sys.executable, CALC_HEAT,
           "--wp-vol", wp_vol,
           "--q-uniform", "1.0e6",
           "--dt", "1.0", "--t-end", "1.0",
           "--output", str(tmp_path / "heat.json"),
           "--msh-output", str(tmp_path / "heat_T.msh")]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    # calc_main exits non-zero on error (radia >= 4.92.0, commit
    # 5b88b67f, "No Fallbacks -- Fail Fast" policy).  Either a
    # successful run (returncode 0) or a clean JSON error
    # (returncode 1 + last stdout line is JSON with "error" key) is
    # acceptable here; a Python traceback on stderr (no JSON in
    # stdout) is not.
    assert proc.returncode in (0, 1), (
        f"calc_heat.py exited with unexpected code {proc.returncode} "
        f"(expected 0 success or 1 clean JSON error, NOT a crash):\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    # The result JSON is the last non-empty stdout line.
    last = [ln for ln in proc.stdout.splitlines() if ln.strip()][-1]
    return json.loads(last)


@pytest.mark.slow
def test_multi_material_mesh_rejected(tmp_path):
    """A coil+workpiece mesh must be rejected with a workpiece-only
    error -- NOT silently solved (which would heat the coil)."""
    multi = _ensure(GEN_MULTI, MULTI_VOL)
    result = _run_heat(multi, tmp_path)
    assert "error" in result, (
        f"thermal accepted a multi-material (coil+wp) mesh; expected a "
        f"workpiece-only rejection.  result keys={sorted(result)}")
    msg = result["error"].lower()
    assert "workpiece" in msg and "material" in msg, (
        f"rejection message should explain the workpiece-only contract; "
        f"got: {result['error']}")


@pytest.mark.slow
def test_single_material_mesh_accepted(tmp_path):
    """A clean single-material workpiece mesh must still solve."""
    wp = _ensure(GEN_WP, WP_VOL)
    result = _run_heat(wp, tmp_path)
    assert "error" not in result, (
        f"thermal rejected a clean single-material workpiece mesh: "
        f"{result.get('error')}")
    assert result.get("ne", 0) > 0, (
        f"expected a solved volume mesh (ne>0); result keys="
        f"{sorted(result)}")
