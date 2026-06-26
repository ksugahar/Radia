"""Golden regression test for BEM-A coil inductance on a VOLUME .vol.

Locks in the v4.40.0 fix for the "K matrix singular" pathology when
calc_inductance.py --coil-solver bem-a is given a Cubit-exported
volume mesh (with internal tetrahedra in addition to the boundary
triangles).

Before 2026-05-12 (Keiko mdx report): _build_bema_coil_mesh returned
the loaded volume mesh as-is, and compute_inductance_source_sink fed
it to HDivSurface(order=0) + LaplaceSL.  The saddle-point K matrix
gained extra null modes that the D[:-1, :] deflation could not
remove, and scipy.linalg.solve raised "A singular matrix detected:
slice(s) [0] are singular."

Fix: _build_bema_coil_mesh extracts a surface-only mesh via
_extract_surface_mesh_filtered before returning, restoring the
pure-surface assumption that compute_inductance_source_sink was
designed for.

Reference values (gapped torus R_major=30 mm, r_minor=3 mm, gap=5 deg):
  * PEEC golden hard band: [80, 92] nH (from
    validation_test/panels/golden/peec_inductance_torus_50kHz_Cu.json)
  * BEM-A on surface mesh (v4.25.0 reference): ~86.7 nH
  * BEM-A on volume .vol (this fix): ~87.3 nH

The two BEM-A values agree to ~0.7% — the same geometry, slightly
different mesh density (surface OCC mesh vs volume Cubit/Netgen mesh).
Both fall well inside the PEEC band.
"""
from __future__ import annotations

import os
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURE = os.path.join(HERE, "fixtures", "gapped_torus_2port.vol")
MAKE_FIXTURE = os.path.join(HERE, "fixtures", "make_gapped_torus_2port.py")

PEEC_BAND_NH = (80.0, 92.0)


def _have_ngsbem():
    try:
        import ngsolve.bem  # noqa: F401
        return True
    except Exception:
        return False


@pytest.fixture(scope="module")
def gapped_torus_vol():
    if not os.path.isfile(FIXTURE):
        # Generate on the fly via the fixture maker; the .vol is local
        # (not in git) since the Netgen mesh is ABI-tied to the runtime.
        if not os.path.isfile(MAKE_FIXTURE):
            pytest.skip(f"fixture maker missing: {MAKE_FIXTURE}")
        import subprocess
        rc = subprocess.run([sys.executable, MAKE_FIXTURE],
                             capture_output=True, text=True)
        if rc.returncode != 0 or not os.path.isfile(FIXTURE):
            pytest.skip(
                f"fixture make failed (rc={rc.returncode}): "
                f"{rc.stderr.splitlines()[-3:] if rc.stderr else ''}")
    return FIXTURE


@pytest.mark.slow
def test_bem_a_inductance_from_volume_vol(gapped_torus_vol):
    """End-to-end BEM-A inductance on the gapped-torus volume .vol.

    Drives the panel-side entry point _solve_coil_bem_a so that the
    volume → surface extraction in _build_bema_coil_mesh is exercised
    (NOT compute_inductance_source_sink directly).  This is the same
    code path the IH panel takes from the Run button.
    """
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")

    # Layer-4 entry-point shim: build a minimal argparse Namespace
    # matching what calc_inductance.run_inductance() would pass.
    from argparse import Namespace
    sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
    from calc_inductance import _solve_coil_bem_a

    args = Namespace(
        coil_vol=gapped_torus_vol,
        coil_source_name="source",
        coil_sink_name="sink",
        coil_rwg_quad_degree=5,
        coil_rwg_singular_nq=6,
        coil_bem_solver="auto",
        coil_aca_eps=1e-10,
        coil_gmres_tol=1e-10,
        coil_saddle_solver="auto",   # added in commit e4a30197 (BEM-A
                                     # saddle solver memory-efficient
                                     # path).  "auto" picks lu / minres
                                     # by n_tris (lu < 10000 < minres).
        coil_sigma=5.8e7,
        frequency=150_000.0,
        current=1.0,
        n_threads=0,
    )
    coil_data = _solve_coil_bem_a(args)

    # _solve_coil_bem_a raises on error; if we got here the solve
    # completed and "L_coil" (in Henrys) must be present.
    assert "L_coil" in coil_data, f"missing L_coil: {list(coil_data)}"
    L_nH = float(coil_data["L_coil"]) * 1e9
    lo, hi = PEEC_BAND_NH
    assert lo <= L_nH <= hi, (
        f"BEM-A on volume .vol: L = {L_nH:.3f} nH outside PEEC band "
        f"[{lo}, {hi}].  Regression: did _build_bema_coil_mesh stop "
        f"calling _extract_surface_mesh_filtered?")
    # Residual at machine precision after dense LU.
    res = float(coil_data["bem_a_residual"])
    assert res < 1e-10, f"div(J) residual too large: {res:.3e}"
