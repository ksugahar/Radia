"""Golden (a): REGCOIL vacuum parity + the deliverable REGCOIL does NOT produce.

`validation_test/stream_function/demo_regcoil_parity_deliverable.py` designs a PRODUCIBLE
vacuum target (uniform vertical B.n) on a torus winding surface and, from the
SAME run:
  - reaches B.n residual ~machine precision  (REGCOIL/NESCOIL forward-map parity)
  - emits a STEP CAD of the coil contours    (REGCOIL does NOT)
  - emits a PEEC L,R of a modular-coil turn  (REGCOIL does NOT)

This locks "design AT PARITY, deliverable BEYOND" -- the design-to-manufacture
vs design-only claim -- as a MEASUREMENT (per Repository-First, not a paper
claim).  Measured on LAB: B.n resid 4.9e-9, STEP 954 kB, L 3.09 uH.

Runs the demo in a SUBPROCESS (NGSolve/OCC heavy import); SKIPS on the LAB
MKL/Qt DLL shadow that breaks NGSolve inside the pytest process.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEMO = os.path.join(REPO, "validation_test", "stream_function",
                    "demo_regcoil_parity_deliverable.py")


def _run(work_dir):
    cmd = [sys.executable, DEMO, "--work-dir", work_dir,
           "--eval-max", "120", "--wind-maxh", "0.06",
           "--plasma-maxh", "0.05", "--n-levels", "12"]
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1200)
    if r.returncode != 0:
        low = (r.stderr or "").lower()
        if ("mkl" in low or "dll" in low or "libiomp" in low
                or (not r.stderr.strip() and abs(r.returncode) > 1000000)):
            pytest.skip("NGSolve/MKL subprocess env issue (LAB pytest); run "
                        "demo_regcoil_parity_deliverable.py directly to verify")
        raise AssertionError(
            f"demo failed (rc={r.returncode}):\nSTDERR:\n{r.stderr[-1500:]}")
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert lines, f"no stdout:\nSTDERR:\n{r.stderr[-1500:]}"
    return json.loads(lines[-1])


def test_regcoil_vacuum_parity_plus_deliverable(tmp_path):
    r = _run(str(tmp_path))
    d = r["design"]
    dl = r["deliverable"]

    # (1) DESIGN -- REGCOIL/NESCOIL forward-map PARITY on a producible target:
    # the winding-surface current potential reproduces uniform-vertical B.n to
    # machine precision (REGCOIL also hits this).  Measured ~4.9e-9; lock a
    # generous machine-precision-level band.
    assert d["bn_residual_rel"] < 1.0e-6, \
        f"vacuum B.n parity not at machine-precision level: {d['bn_residual_rel']}"
    assert d["n_contours"] >= 4, f"too few coils: {d['n_contours']}"

    # (2) DELIVERABLE -- what REGCOIL/NESCOIL/FOCUS do NOT emit, from the SAME run:
    #     a valid STEP CAD of the coil contours ...
    assert dl["step_ok"] is True, "STEP CAD was not produced"
    assert dl["step_bytes"] > 1000, f"STEP file implausibly small: {dl['step_bytes']}"
    assert dl["step_file"].endswith(".step")
    assert os.path.exists(dl["step_file"]), "STEP file missing on disk"
    #     ... and a PEEC inductance of a modular-coil turn (finite, positive L).
    assert math.isfinite(dl["peec_L_H"]) and dl["peec_L_H"] > 0.0, \
        f"PEEC L not a real positive inductance: {dl['peec_L_H']}"
    assert math.isfinite(dl["peec_R_ohm"]) and dl["peec_R_ohm"] >= 0.0
    assert dl["peec_n_nodes"] >= 3, "PEEC contour not a closed coil"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
