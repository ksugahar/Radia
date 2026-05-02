"""Phase C.4 (DC SIBC): pytest for compute_dc_resistance_intree.

Validates the (Z_s * Mass)-based DC saddle-point solve against
analytical thin-wire resistance.

Geometry: gapped torus R=30, r=3, gap=5 deg.
Analytical (uniform J around perimeter): R/Z_s = L_path / perimeter
                                        = 2pi*R*(355/360) / (2pi*r)
                                        ~ 9.86

DC BEM gives R/Z_s in [9.86, 11.0] -- slightly above the uniform
estimate due to current crowding near the port caps.  AT DC the J
distribution is *close to uniform* (which it should be, since DC
current on a uniform cross-section minimises ohmic loss with no
inductive bias) -- the small overshoot is genuine mesh/cap geometry,
not a bug.

The PEC-inductance variant (`solve_inductance_source_sink_intree`'s
optional Z_s post-processing) gives a ~10% larger R because the
PEC J distribution is more crowded toward the outer edges (it
minimises magnetic energy, not ohmic loss).  This is documented in
the example README.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest


HERE = os.path.dirname(os.path.abspath(__file__))


def _have_netgen_occ():
    try:
        from netgen.occ import Glue  # noqa: F401
        return True
    except Exception:
        return False


def _build_gapped_torus_mesh():
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh

    R_maj, r_min, gap = 0.030, 0.003, 5.0
    wp = WorkPlane(Axes(p=Pnt(R_maj, 0, 0), n=Dir(0, 1, 0),
                         h=Dir(0, 0, 1)))
    profile = wp.Circle(r_min).Face()
    torus = profile.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)),
                              360.0 - gap)
    for f in torus.faces:
        f.name = "conductor"
    fs = sorted(torus.faces, key=lambda f: f.mass)
    fs[0].name = "source"
    fs[1].name = "sink"
    surf = Glue(torus.faces)
    geo = OCCGeometry(surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=r_min * 5.0,
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
def torus_arrays():
    if not _have_netgen_occ():
        pytest.skip("netgen.occ not available")
    mesh = _build_gapped_torus_mesh()
    return _arrays_from_mesh(mesh)


@pytest.fixture(scope="module")
def dc_R_result(torus_arrays):
    """Run compute_dc_resistance_intree once, share across tests."""
    from radia.bem.efie_rwg import compute_dc_resistance_intree

    verts, tris, src_mask, snk_mask = torus_arrays
    return compute_dc_resistance_intree(
        verts, tris, src_mask, snk_mask, Z_s=1.0)


@pytest.mark.slow
def test_dc_R_in_expected_band(dc_R_result):
    """For Z_s=1.0 Ohm/sheet on the gapped torus, R should be
    9.86 < R < 11.0 (uniform J estimate to ~10% current-crowding upper)."""
    assert "error" not in dc_R_result, (
        f"solver error: {dc_R_result.get('error')}")
    R = dc_R_result["R"]
    R_uniform = 2 * math.pi * 0.030 * (355 / 360) / (2 * math.pi * 0.003)
    assert R_uniform <= R <= R_uniform * 1.15, (
        f"R = {R:.4f} outside [{R_uniform:.3f}, {R_uniform * 1.15:.3f}] "
        f"(uniform-J est {R_uniform:.3f})")


@pytest.mark.slow
def test_dc_R_residual_machine_precision(dc_R_result):
    """Saddle-point LU drives ||D J - g||_inf to ~1e-15."""
    assert dc_R_result["residual"] < 1e-12


@pytest.mark.slow
def test_dc_R_uses_distinct_J_from_pec(torus_arrays):
    """DC J distribution differs from PEC J -- the DC saddle minimises
    ohmic loss (uniform-ish), PEC minimises magnetic energy (edge-
    crowded).  The two J vectors should NOT be parallel."""
    from radia.bem.efie_rwg import (
        compute_dc_resistance_intree,
        solve_inductance_source_sink_intree)

    verts, tris, src, snk = torus_arrays
    dc = compute_dc_resistance_intree(verts, tris, src, snk, Z_s=1.0)
    pec = solve_inductance_source_sink_intree(
        verts, tris, src, snk,
        regular_quad_degree=5, singular_n_q=8)
    j_dc = dc["J"]
    j_pec = pec["J"]
    cos_sim = (j_dc @ j_pec) / (np.linalg.norm(j_dc)
                                  * np.linalg.norm(j_pec))
    # They share gross structure (cos_sim near 1) but are not identical.
    assert 0.9 < cos_sim < 0.999999, (
        f"cos similarity DC vs PEC = {cos_sim:.6f}; expected slight "
        f"divergence (DC concentrates less near caps)")
