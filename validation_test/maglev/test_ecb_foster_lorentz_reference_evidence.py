"""ECB plate force: replay the reference evidence, and hold the shipped kernel to it.

The reference lane (``ecb_foster_lorentz_reference.py``) solves the scalar
model directly and reconstructs the eddy current as (1/mu) curl(v z).  The
shipped ``compute_lorentz_force_via_foster`` reconstructs it as
-omega sigma Im(v), which gives no lift for a centred magnet and a horizontal
force that mirror symmetry forbids.  The last test runs the kernel on a small
mesh and is marked strict-xfail: it starts passing -- and therefore failing the
suite -- the moment the kernel is corrected, so the record cannot go stale.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SUMMARY = HERE / "ecb_foster_lorentz_reference_summary.json"


def _payload():
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def test_reference_summary_is_complete_and_physical():
    payload = _payload()
    assert payload["schema"] == "radia.maglev.ecb-foster-lorentz-reference.v1"
    for key in ("radia_version", "ngsolve_version", "python_version", "host"):
        assert payload["runtime"][key]
    assert payload["reference_passed"] is True
    assert all(payload["checks"].values()), payload["checks"]


def test_reference_lift_values_are_locked():
    """Centred lift on the finest mesh at 50 / 500 / 5000 Hz (conductor side)."""
    payload = _payload()
    fine = payload["problem"]["meshes"][-1]
    centred = {c["frequency_hz"]: c["reference_force_N"][2] for c in payload["cases"]
               if c["mesh"] == fine and c["x_pm_m"] == 0.0}
    assert centred[50.0] == pytest.approx(-2.5697e-2, rel=5e-3)
    assert centred[500.0] == pytest.approx(-1.0156, rel=5e-3)
    assert centred[5000.0] == pytest.approx(-1.9462, rel=1e-2)
    assert max(-v for v in centred.values()) < payload["problem"]["image_lift_bound_N"]


def test_shipped_kernel_defect_is_recorded():
    record = _payload()["shipped_kernel"]
    assert record["violates_mirror_symmetry"] is True
    assert record["produces_no_lift_when_centred"] is True
    assert record["exceeds_image_bound"] is True


@pytest.mark.xfail(
    strict=True,
    reason="compute_lorentz_force_via_foster reconstructs J_y = -omega sigma Im(v) "
           "instead of (1/mu) curl(v z): no lift when centred, forbidden horizontal "
           "force. Remove this marker when the kernel is corrected.",
)
def test_shipped_kernel_satisfies_the_centred_symmetry_checks():
    pytest.importorskip("ngsolve")
    from ngsolve import TaskManager

    import ecb_foster_lorentz_reference as lane
    from radia.maglev.ecb.lorentz import compute_lorentz_force_via_foster
    from radia.maglev.mixed_galerkin.alpha import _dirichlet_eigenmodes

    with TaskManager():
        mesh = lane.plate_mesh(20, 8, 2)
        lam, vecs, _mass, free, _fes, _volume = _dirichlet_eigenmodes(mesh, 60, "outer")
        horizontal, vertical = compute_lorentz_force_via_foster(
            mesh, lam, vecs, free, lane.SIGMA, lane.MU0,
            2j * math.pi * 500.0, lane.M_PM, lane.Z_PM, 0.0)
    assert abs(horizontal) < 1e-6 * max(abs(vertical), 1e-30)
    assert vertical < 0.0
