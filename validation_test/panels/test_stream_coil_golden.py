"""Golden test for calc_stream_coil.py (stream-function Gz gradient coil).

Stage-2 validation lock for the cylindrical z-gradient (Gz) gradient-coil
designer.  Runs the CLI with default args via subprocess, parses the JSON on
stdout, and asserts the key numbers lie inside a hard band:

- fitted_dBdz in [0.9, 1.1]: the discrete wire coil reproduces the target
  unit gradient (Gz = 1.0 T/m) to within +/-10%.
- gradient_nonlinearity < 0.05: on-axis Bz stays linear over the DSV.
- k_aca >= 1: ACA+ found at least one mode (solve is non-trivial).
- len(wire_z) == len(wire_I) > 0: the contour produced paired wire rings.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PANELS = os.path.join(REPO, "src", "radia", "panels")


def _run(extra_args):
    cmd = [
        sys.executable,
        os.path.join(PANELS, "calc_stream_coil.py"),
    ] + extra_args
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert r.returncode == 0, (
        f"calc_stream_coil.py failed (rc={r.returncode}):\n"
        f"STDERR:\n{r.stderr[-1500:]}"
    )
    # JSON is the last non-empty stdout line.
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert lines, f"no stdout from calc_stream_coil.py:\nSTDERR:\n{r.stderr[-1500:]}"
    return json.loads(lines[-1])


def test_stream_coil_default_golden():
    """Default-args design: gradient fidelity + linearity inside hard band."""
    result = _run([])

    assert "error" not in result, f"calc returned error: {result.get('error')}"

    fitted = result["fitted_dBdz"]
    nonlin = result["gradient_nonlinearity"]
    k_aca = result["k_aca"]
    wire_z = result["wire_z"]
    wire_I = result["wire_I"]

    assert 0.9 <= fitted <= 1.1, \
        f"fitted_dBdz out of band [0.9, 1.1]: {fitted}"
    assert nonlin < 0.05, \
        f"gradient_nonlinearity not < 0.05: {nonlin}"
    assert k_aca >= 1, f"k_aca must be >= 1: {k_aca}"
    assert len(wire_z) == len(wire_I) > 0, \
        f"wire_z/wire_I mismatch or empty: {len(wire_z)} vs {len(wire_I)}"


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
