"""Golden lock for the radia-electromagnet "Clebsch hodograph" mode.

Runs ``calc_clebsch_hodograph.py`` as a subprocess (no panel) for the two
verified self-meshed geometries and asserts the BIDIRECTIONAL CONSISTENCY
(the a-posteriori field-quality metric -- B reconstructed from the flux psi
vs from the scalar Phi) stays inside a hard band, plus the standard
housekeeping keys and the visualization .msh.  Portable: needs only
NGSolve + radia.

Forward mode only.  The 3-D inverse pole-design (geometry-as-unknown
hodograph PDE) is a separate research line and is NOT covered here.

Verified 2026-06-12 (maxh=0.10, order=3):
  cylinder  consistency ~7e-4  field_error ~9e-4   (path B, 2-D dipole)
  sphere    consistency ~7e-4  field_error ~3e-3   (path A, axisym, axis-
            anchored psi; field_error = interior Hz vs 3/(mu_r+2) H0 with
            the far boundary truncated at r/a = 6, NOT exact-BC fed).
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

CALC = (Path(__file__).resolve().parents[2]
        / "src" / "radia" / "panels" / "calc_clebsch_hodograph.py")


def _run(tmp_path, geometry):
    out_json = tmp_path / f"{geometry}.json"
    out_msh = tmp_path / f"{geometry}.msh"
    cmd = [sys.executable, str(CALC),
           "--geometry", geometry,
           "--maxh", "0.10", "--fes-order", "3",
           "--output", str(out_json),
           "--msh-output", str(out_msh)]
    subprocess.run(cmd, check=True, timeout=600)
    return json.loads(out_json.read_text()), out_msh


@pytest.mark.parametrize("geometry,cons_max,ferr_max", [
    ("cylinder", 2e-3, 2e-3),
    ("sphere",   2e-3, 1e-2),
])
def test_clebsch_consistency_band(tmp_path, geometry, cons_max, ferr_max):
    got, out_msh = _run(tmp_path, geometry)
    assert "error" not in got, f"calc errored: {got.get('error')}"

    # Headline: hodograph self-consistency (B from psi vs B from Phi),
    # evaluated on the SAME FEM field -- locks that the reduced-Omega solve
    # and the conjugate psi / Phi reconstruction stay mutually consistent.
    assert 0.0 < got["consistency"] < cons_max, (
        f"{geometry} consistency {got['consistency']:.2e} "
        f"outside (0, {cons_max:.0e})")

    # FEM field accuracy vs the analytical solution.
    assert got["field_error"] < ferr_max, (
        f"{geometry} field_error {got['field_error']:.2e} >= {ferr_max:.0e}")

    # Standard housekeeping keys (Result Output Policy 2026-05-29).
    for k in ("ne", "ndof", "t_total_s", "consistency", "gmsh_file"):
        assert k in got, f"missing key {k!r}: {sorted(got)}"

    # Visualization .msh really got written and is non-empty.
    assert out_msh.exists() and out_msh.stat().st_size > 0, (
        f"{geometry}: .msh not written")
