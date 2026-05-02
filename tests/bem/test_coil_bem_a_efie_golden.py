"""Golden regression tests for BEM-A EFIE coil self-inductance.

Locks in the saddle-point EFIE on HDivSurface x SurfaceL2 (Weggler
stabilized formulation, equivalent to Loop-Star / Lucy decomposition
at the variational level) on the gapped-torus geometry.

Cross-checks against the PEEC golden
``tests/panels/golden/peec_inductance_torus_50kHz_Cu.json``:
    R_maj=30 mm, r_min=3 mm, gap=5 deg
    expected: L = 85.10 nH  (PEEC golden hard band [80, 92])

The reference solver
``examples/induction_heating/bem_reference/bem_inductance.py``
::``compute_inductance_source_sink`` uses
``ngsolve.bem.LaplaceSL`` on ``HDivSurface`` -- slow at this mesh size
(~104s assembly at N_J=5064) but correct.  This is the validation
oracle until intree HACApK assembly lands.

Run::

    pytest tests/bem/test_coil_bem_a_efie_golden.py -v -s
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

# Make the legacy bem_inductance.py importable.
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


def _make_gapped_torus_surface_mesh(R_major=0.030, r_minor=0.003,
                                     gap_deg=5.0, maxh=0.003):
    """Build the BEM-A reference fixture (surface-only OCC torus with
    source / sink / conductor BND labels)."""
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh

    wp = WorkPlane(Axes(p=Pnt(R_major, 0, 0), n=Dir(0, 1, 0),
                         h=Dir(0, 0, 1)))
    circle = wp.Circle(r_minor).Face()
    sweep_deg = 360.0 - gap_deg
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), sweep_deg)
    for f in torus.faces:
        f.name = "conductor"
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"
    surf = Glue(torus.faces)
    geo = OCCGeometry(surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=maxh, curvaturesafety=2.0,
                              segmentsperedge=2))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)
    return mesh


@pytest.fixture(scope="module")
def gapped_torus_mesh():
    """Cached gapped-torus surface mesh for all tests in this module."""
    if not _have_netgen_occ():
        pytest.skip("netgen.occ not available")
    return _make_gapped_torus_surface_mesh(maxh=0.003)


@pytest.fixture(scope="module")
def efie_result(gapped_torus_mesh):
    """Run compute_inductance_source_sink ONCE per module so both tests
    in this module share the slow ngsolve.bem assembly (~110s)."""
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")
    from bem_inductance import compute_inductance_source_sink

    res = compute_inductance_source_sink(
        gapped_torus_mesh, source_label="source", sink_label="sink",
        fes_order=0, solver="lu")
    if "error" in res:
        pytest.fail(f"compute_inductance_source_sink: {res['error']}")
    return res


@pytest.mark.slow
def test_efie_gapped_torus_PEC_in_PEEC_band(efie_result):
    """BEM-A EFIE saddle-point at HDivSurface order=0 (RWG) gives L
    inside the PEEC golden hard band [80, 92] nH.

    PEEC golden 85.10 nH; analytical (Maxwell thin-wire, no gap) 89.80 nH.
    BEM-A is expected to land between these (typically ~86-87 nH at
    this mesh density), since Maxwell underestimates the gap effect
    and PEEC underestimates due to perimeter discretisation.
    """
    L_nH = efie_result["L"] * 1e9
    band_lo, band_hi = 80.0, 92.0
    assert band_lo <= L_nH <= band_hi, (
        f"L = {L_nH:.3f} nH outside PEEC golden hard band "
        f"[{band_lo}, {band_hi}]")

    # Current conservation residual on the divergence: ||D J - g||_inf
    # should be at machine precision after dense LU solve.
    assert efie_result["residual"] < 1e-10, (
        f"div(J) residual too large: {efie_result['residual']:.3e}")


@pytest.mark.slow
def test_efie_gapped_torus_port_areas_match(efie_result):
    """Source and sink cap areas should equal pi * r_minor^2 within
    mesh discretisation tolerance.  Catches regression in the OCC face
    classifier (smallest two faces -> source / sink)."""
    A_expected = math.pi * (0.003) ** 2   # = 28.27 mm^2
    A_src = efie_result["A_source"]
    A_snk = efie_result["A_sink"]
    rel_err_src = abs(A_src - A_expected) / A_expected
    rel_err_snk = abs(A_snk - A_expected) / A_expected
    assert rel_err_src < 0.10, (
        f"source area off by {rel_err_src*100:.1f}% "
        f"(got {A_src*1e6:.3f} mm^2, expected {A_expected*1e6:.3f} mm^2)")
    assert rel_err_snk < 0.10, (
        f"sink area off by {rel_err_snk*100:.1f}%")
    # Source and sink should be very close in area (geometry symmetric).
    rel_diff = abs(A_src - A_snk) / max(A_src, A_snk)
    assert rel_diff < 0.05, (
        f"source/sink area mismatch {rel_diff*100:.2f}% -- check OCC "
        f"face classifier")
