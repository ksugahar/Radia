"""Phase C.2 golden test: intree RWG saddle-point EFIE matches
ngsolve.bem.LaplaceSL on HDivSurface(0) end-to-end.

Asserts:
  - L_self (intree) and L_self (ngsbem oracle) agree to <= 1e-4
    relative on the gapped torus surface mesh (R=30, r=3, gap=5°,
    maxh = 5*r → ~700 BND tris).
  - intree div(J) residual is at machine precision (saddle solve).
  - intree L_self is in the PEEC golden hard band [80, 92] nH.

If this test ever fails, suspect the SS Duffy sign-handling
regression (non-CCW perm_b flips RWG basis sign on the shared edge --
see ``project_bem_a_intree_decision_shimizu.md`` BUG record).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
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


def _build_gapped_torus_mesh(R_major=0.030, r_minor=0.003, gap_deg=5.0,
                              maxh_factor=5.0):
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh

    wp = WorkPlane(Axes(p=Pnt(R_major, 0, 0), n=Dir(0, 1, 0),
                         h=Dir(0, 0, 1)))
    profile = wp.Circle(r_minor).Face()
    sweep_deg = 360.0 - gap_deg
    torus = profile.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), sweep_deg)
    for f in torus.faces:
        f.name = "conductor"
    fs = sorted(torus.faces, key=lambda f: f.mass)
    fs[0].name = "source"
    fs[1].name = "sink"
    surf = Glue(torus.faces)
    geo = OCCGeometry(surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=r_minor * maxh_factor,
                              curvaturesafety=1.0, segmentsperedge=1))
    mesh = Mesh(ngmesh)
    return mesh


def _arrays_from_mesh(mesh):
    from ngsolve import BND
    n_v = mesh.nv
    verts = np.empty((n_v, 3), dtype=np.float64)
    for i in range(n_v):
        verts[i] = [mesh.vertices[i].point[k] for k in range(3)]
    bnd_names = list(mesh.GetBoundaries())
    src_idx = {i for i, n in enumerate(bnd_names) if n == "source"}
    snk_idx = {i for i, n in enumerate(bnd_names) if n == "sink"}
    tris, src_mask, snk_mask = [], [], []
    for el in mesh.Elements(BND):
        tris.append([int(v.nr) for v in el.vertices])
        src_mask.append(el.index in src_idx)
        snk_mask.append(el.index in snk_idx)
    return (verts, np.array(tris, dtype=np.int64),
            np.array(src_mask, dtype=bool),
            np.array(snk_mask, dtype=bool))


@pytest.fixture(scope="module")
def gapped_torus_arrays():
    if not _have_netgen_occ():
        pytest.skip("netgen.occ not available")
    mesh = _build_gapped_torus_mesh()
    verts, tris, src_mask, snk_mask = _arrays_from_mesh(mesh)
    return mesh, verts, tris, src_mask, snk_mask


@pytest.fixture(scope="module")
def L_intree(gapped_torus_arrays):
    """Intree saddle-point EFIE, returns dict result."""
    from radia.bem.efie_rwg import solve_inductance_source_sink_intree

    _, verts, tris, src_mask, snk_mask = gapped_torus_arrays
    return solve_inductance_source_sink_intree(
        verts, tris, src_mask, snk_mask,
        regular_quad_degree=5, singular_n_q=8)


@pytest.fixture(scope="module")
def L_ngsbem(gapped_torus_arrays):
    """ngsolve.bem oracle, returns dict result."""
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")
    from bem_inductance import compute_inductance_source_sink

    mesh = gapped_torus_arrays[0]
    return compute_inductance_source_sink(
        mesh, source_label="source", sink_label="sink",
        fes_order=0, solver="lu")


@pytest.mark.slow
def test_intree_L_matches_ngsbem(L_intree, L_ngsbem):
    """Intree L_self matches ngsolve.bem oracle to <= 1e-4 relative."""
    L_in = L_intree["L"]
    L_ng = L_ngsbem["L"]
    rel = abs(L_in - L_ng) / abs(L_ng)
    assert rel < 1e-4, (
        f"intree {L_in*1e9:.4f} nH vs ngsbem {L_ng*1e9:.4f} nH: "
        f"rel diff {rel:.3e} > 1e-4")


@pytest.mark.slow
def test_intree_div_residual_machine_precision(L_intree):
    """Saddle-point LU should drive ||D J - g||_inf to ~1e-15."""
    assert L_intree["residual"] < 1e-12, (
        f"div(J) residual {L_intree['residual']:.3e} too large")


@pytest.mark.slow
def test_intree_L_in_PEEC_band(L_intree):
    """Intree L_self lies in PEEC golden hard band [80, 92] nH."""
    L_nH = L_intree["L"] * 1e9
    assert 80.0 <= L_nH <= 92.0, (
        f"intree L = {L_nH:.3f} nH outside PEEC golden band [80, 92]")
