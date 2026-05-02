"""Order-convergence regression for BEM-A EFIE on gapped torus.

Checks that the Weggler stabilized EFIE (HDivSurface(p) x SurfaceL2(p-1)
product space) gives a stable L_self at fes_order = 0, 1.

Order=2 is currently NOT exercised here -- ngsolve.bem.LaplaceSL on
HDivSurface(2) at this mesh size requires several minutes of assembly
(O(N^2) column extraction) and the test would dominate CI time.  The
intree HACApK-based assembler (Phase D.4 / future) will let us re-add
order=2 to this regression cheaply.

Run::

    pytest tests/bem/test_coil_bem_a_order_convergence.py -v -s
"""
from __future__ import annotations

import math
import os
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BEM_REF_DIR = os.path.join(
    ROOT, "examples", "induction_heating", "bem_reference")
if BEM_REF_DIR not in sys.path:
    sys.path.insert(0, BEM_REF_DIR)


def _have_netgen_occ():
    try:
        from netgen.occ import Glue  # noqa: F401
        return True
    except Exception:
        return False


def _have_ngsbem():
    try:
        from ngsolve.bem import LaplaceSL  # noqa: F401
        return True
    except Exception:
        return False


def _make_coarse_gapped_torus():
    """A COARSE gapped-torus surface mesh so order>=1 finishes fast.

    maxh = 2 * r_minor (vs production maxh = r_minor), giving roughly
    1/4 the BND triangles and ~1/16 the assembly time for LaplaceSL.
    """
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh

    R_maj, r_min, gap_deg = 0.030, 0.003, 5.0
    wp = WorkPlane(Axes(p=Pnt(R_maj, 0, 0), n=Dir(0, 1, 0),
                         h=Dir(0, 0, 1)))
    profile = wp.Circle(r_min).Face()
    torus = profile.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)),
                              360.0 - gap_deg)
    for f in torus.faces:
        f.name = "conductor"
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"
    surf = Glue(torus.faces)
    geo = OCCGeometry(surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=2 * r_min, curvaturesafety=2.0,
                              segmentsperedge=2))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)
    return mesh


@pytest.fixture(scope="module")
def coarse_mesh():
    if not _have_netgen_occ():
        pytest.skip("netgen.occ not available")
    return _make_coarse_gapped_torus()


@pytest.mark.skipif(not _have_ngsbem(), reason="ngsolve.bem not available")
@pytest.mark.slow
@pytest.mark.parametrize("fes_order", [0, 1])
def test_order_in_PEEC_band(coarse_mesh, fes_order):
    """L stays inside the PEEC hard band [80, 92] nH at fes_order = 0, 1.

    The coarse mesh sacrifices a couple percent of accuracy vs the
    finer fixture, so we widen the band slightly to [78, 95] to allow
    for that without making the test brittle.
    """
    from bem_inductance import compute_inductance_source_sink

    res = compute_inductance_source_sink(
        coarse_mesh, source_label="source", sink_label="sink",
        fes_order=fes_order, solver="lu")
    assert "error" not in res, (
        f"order={fes_order}: solver error: {res.get('error')}")
    L_nH = res["L"] * 1e9
    band_lo, band_hi = 78.0, 95.0
    assert band_lo <= L_nH <= band_hi, (
        f"order={fes_order}: L = {L_nH:.3f} nH outside relaxed "
        f"hard band [{band_lo}, {band_hi}] (PEEC golden 85.10 nH).  "
        f"n_J={res['n_J']}, residual={res['residual']:.3e}")
    assert res["residual"] < 1e-9, (
        f"order={fes_order}: div(J) residual {res['residual']:.3e}")
