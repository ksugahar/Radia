"""Golden (b): the IRON differentiator -- free-space design MISSES, material-aware HITS.

NESCOIL/REGCOIL/FOCUS design coils with the FREE-SPACE Biot-Savart kernel; they
cannot see magnetic material.  examples/kelvin_transformation/DtN_spectrum/
act8_03_general_iron_design.py SCENARIO B shows why that is a real limit: with
NON-CONCENTRIC iron (no closed-form Green function) the Kelvin-FEM material-aware
transfer M is the correct design kernel.

`_design_once` designs a coil two ways for the SAME producible target and then
re-evaluates each in an INDEPENDENT full Kelvin-FEM forward solve WITH the iron:
  - invert the material-aware M (iron-aware)  -> err_iron, HITS (~machine precision)
  - invert the free-space M (= REGCOIL)       -> err_free, MISSES by tens of %
The miss is a PHYSICAL shield effect (stable under mesh refinement), not an
artifact.  This golden locks the (b) claim -- "REGCOIL misses with iron, radia's
material-aware Kelvin-DtN matches" -- as a MEASUREMENT (Repository-First).

Runs _design_once in a SUBPROCESS (NGSolve Kelvin-FEM heavy import); SKIPS on the
LAB MKL/Qt DLL shadow that breaks NGSolve inside the pytest process.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DTN = os.path.join(REPO, "examples", "kelvin_transformation", "DtN_spectrum")

# import act8_03 + run ONE design at the coarsest resolution scenario_B uses
# (maxh 0.30 -- fast, and where the free-space miss is already large).
_RUN = (
    "import sys, json\n"
    "sys.path.insert(0, r'%s')\n"
    "import act8_03_general_iron_design as A\n"
    "ei, ef, nc, nt = A._design_once(('blob', (0.0, 0.0, 0.68), 0.16), 50.0, 2, 0.30)\n"
    "print(json.dumps({'err_iron': ei, 'err_free': ef, 'n_coil': nc, 'n_target': nt}))\n"
    % DTN.replace("\\", "/")
)


def test_iron_differentiator_free_space_misses():
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([sys.executable, "-c", _RUN], capture_output=True,
                       text=True, env=env, timeout=1200)
    if r.returncode != 0:
        low = (r.stderr or "").lower()
        if ("mkl" in low or "dll" in low or "libiomp" in low
                or (not r.stderr.strip() and abs(r.returncode) > 1000000)):
            pytest.skip("NGSolve/MKL subprocess env issue (LAB pytest); run "
                        "act8_03_general_iron_design.py directly to verify")
        raise AssertionError(
            f"act8_03 _design_once failed (rc={r.returncode}):\n"
            f"STDERR:\n{r.stderr[-1500:]}")
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert lines, f"no stdout:\nSTDERR:\n{r.stderr[-1500:]}"
    d = json.loads(lines[-1])

    # material-aware design HITS (inverts the exact iron operator).
    assert d["err_iron"] < 1.0e-2, \
        f"material-aware design did not hit the target: err_iron={d['err_iron']}"
    # free-space (REGCOIL) design MISSES by a lot when the iron is actually there
    # (scenario_B measured 31-39%; lock a conservative >20%).
    assert d["err_free"] > 0.20, \
        f"free-space design did not miss in iron (err_free={d['err_free']}); " \
        f"the iron reaction is not engaging?"
    # the differentiator: the free-space miss is orders above the material-aware
    # hit (measured ratio ~300x).
    assert d["err_free"] / max(d["err_iron"], 1.0e-12) > 10.0, \
        f"free-space miss not >> material-aware hit: " \
        f"{d['err_free']} vs {d['err_iron']}"
    assert d["n_coil"] >= 4 and d["n_target"] >= 8


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
