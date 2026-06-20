"""Golden (c): manufacturing END-TO-END -- the full design-to-manufacture pipeline
that NESCOIL/REGCOIL/FOCUS do NOT provide.

ONE calc_streamfunction.py --method manufacture run takes a target all the way:

  target B  ->  current potential psi  ->  iso-contours  ->  SINGLE-STROKE wire
  (field-aware chain)  ->  SHEET-METAL distort  ->  STEP CAD  ->  PEEC L,R

The design codes stop at psi.  This golden locks that radia emits, from ONE run,
a wound + distorted + CAD'd + inductance-characterized SINGLE conductor on the
canonical cylinder-former + DSV (an MRI-gradient-like single-conductor coil --
where the single-stroke wire IS the manufacturing deliverable):
  - a design (continuous) homogeneity ceiling,
  - a finite delivered single-stroke WIRE homogeneity,
  - distort that does not worsen the wire (it is a minimiser),
  - a valid STEP CAD file,
  - a positive PEEC inductance L.

Item (c) of "surpass NESCOIL/REGCOIL": the complete design-to-manufacture chain
in a single run.  Runs in a SUBPROCESS; SKIPS on the LAB MKL/Qt DLL shadow.
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
CALC = os.path.join(REPO, "src", "radia", "panels", "calc_streamfunction.py")
FIXTURE = os.path.join(HERE, "fixtures", "make_streamfunction_vol.py")


@pytest.fixture(scope="module")
def sample_vols(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("sf_e2e_vols"))
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([sys.executable, FIXTURE, d],
                       capture_output=True, text=True, env=env, timeout=600)
    if r.returncode != 0:
        pytest.skip(f"sample mesh gen failed (NGSolve/OCC unavailable?):\n"
                    f"{r.stderr[-600:]}")
    coil = os.path.join(d, "coil_cyl_surf.vol")
    evalv = os.path.join(d, "eval_dsv.vol")
    if not (os.path.exists(coil) and os.path.exists(evalv)):
        pytest.skip("sample .vol files not produced")
    return coil, evalv


def test_manufacture_end_to_end_deliverable(sample_vols, tmp_path):
    coil, evalv = sample_vols
    step = str(tmp_path / "coil_e2e.step")
    cmd = [sys.executable, CALC, "--coil-vol", coil, "--eval-vol", evalv,
           "--order", "1", "--target-cf", "x", "--eval-max", "40",
           "--confine", "abe", "--method", "manufacture", "--nlevels", "16",
           "--chain", "field_aware", "--distort", "--step-output", step,
           "--peec", "--peec-freq", "1e5"]
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1200)
    if r.returncode != 0:
        low = (r.stderr or "").lower()
        if ("mkl" in low or "dll" in low or "libiomp" in low
                or (not r.stderr.strip() and abs(r.returncode) > 1000000)):
            pytest.skip("NGSolve/MKL subprocess env issue (LAB pytest); run "
                        "calc_streamfunction.py directly to verify")
        raise AssertionError(
            f"manufacture e2e failed (rc={r.returncode}):\nSTDERR:\n{r.stderr[-1500:]}")
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert lines, f"no stdout:\nSTDERR:\n{r.stderr[-1500:]}"
    d = json.loads(lines[-1])
    assert "error" not in d, f"manufacture error: {d.get('error')}"
    assert d["method"] == "manufacture"

    # --- the full chain, end to end ---
    # design (continuous) ceiling
    assert 0.0 < d["homogeneity_rms"] < 0.1, \
        f"design homogeneity out of band: {d['homogeneity_rms']}"
    # delivered single-stroke WIRE (finite, positive) + distort recorded
    wa = d["wire_homogeneity_rms"]
    wb = d.get("wire_homogeneity_rms_before")
    assert math.isfinite(wa) and wa > 0.0, f"no delivered wire: {wa}"
    assert wb is not None, "distort did not run (no _before recorded)"
    # distort is a minimiser -> it never worsens the wire it started from
    assert wa <= wb * 1.001, f"distort worsened the wire: {wb} -> {wa}"
    # STEP CAD emitted from the SAME run (REGCOIL does not)
    assert d.get("step") == step, "STEP path not echoed"
    assert os.path.exists(step) and os.path.getsize(step) > 1000, \
        "STEP CAD not written / implausibly small"
    # PEEC inductance from the SAME run (REGCOIL does not)
    assert d["peec"]["L_H"] > 0.0 and d["peec"]["R_ohm"] >= 0.0, \
        f"PEEC L,R not physical: {d['peec']}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
