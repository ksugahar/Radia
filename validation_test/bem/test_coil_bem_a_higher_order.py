"""BEM-A impedance-EFIE beyond RT0, and on a volume mesh without extraction.

Until 2026-09-21 ``hacapk_cocr`` was RT0-only, because its H-matrix cluster
tree took one point per DOF from the edge midpoints under the assumption
``DOF i == edge i``.  The beak-fin comparison showed that RT0 does not
p-converge the current crowding at a rounded tip -- one order step moved R by
+1.69% where curving the geometry moved it by +0.067% -- so the compressed
solver had to follow the basis.  ``_dof_cluster_coords`` attributes DOFs
through the boundary elements (shared by two -> edge midpoint, owned by one
-> face centroid) and works for any order.

The same change lets the solver take a volume mesh directly: ``Compress``
removes the unused interior-edge DOFs whose null modes made the dense saddle
singular in 2026-05 (the fix then was to extract a flat surface, which throws
the curvature away; this is the cure).

These tests lock:

1. ``_dof_cluster_coords`` reproduces ``_edge_midpoint_coords`` exactly at
   order 0, so the historical RT0 H-matrix is unchanged.
2. ``hacapk_cocr`` and ``cocr`` agree at ``fes_order=1`` the way they are
   documented to agree at order 0.
3. The volume mesh, compressed inside the solver, gives the RT0 result of the
   extracted surface: same DOF count, same R and L.
"""
from __future__ import annotations

import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURE = os.path.join(HERE, "fixtures", "gapped_torus_2port.vol")
MAKE_FIXTURE = os.path.join(HERE, "fixtures", "make_gapped_torus_2port.py")

SIGMA_CU = 5.8e7
MU_0 = 4e-7 * math.pi


def _have_ngsbem():
    try:
        import ngsolve.bem  # noqa: F401
        return True
    except Exception:
        return False


def _zs(sigma, omega):
    delta = math.sqrt(2.0 / (omega * MU_0 * sigma))
    return (1.0 + 1.0j) / (sigma * delta)


@pytest.fixture(scope="module")
def gapped_torus_vol():
    if not os.path.isfile(FIXTURE):
        if not os.path.isfile(MAKE_FIXTURE):
            pytest.skip(f"fixture maker missing: {MAKE_FIXTURE}")
        import subprocess
        rc = subprocess.run([sys.executable, MAKE_FIXTURE],
                            capture_output=True, text=True)
        if rc.returncode != 0 or not os.path.isfile(FIXTURE):
            pytest.skip(f"fixture make failed (rc={rc.returncode})")
    return FIXTURE


def _surface_and_volume(volpath):
    from ngsolve import Mesh
    sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
    from surface_mesh_extract import _extract_surface_mesh_filtered
    volume = Mesh(volpath)
    surface = _extract_surface_mesh_filtered(volume, keep_label="")
    return surface, volume


def test_dof_cluster_coords_reproduce_edge_midpoints_at_rt0(gapped_torus_vol):
    import numpy as np
    from ngsolve import HDivSurface
    from radia.bem.coil_inductance_ngsolve import (
        _dof_cluster_coords, _edge_midpoint_coords)

    surface, _ = _surface_and_volume(gapped_torus_vol)
    fes = HDivSurface(surface, order=0)
    new = _dof_cluster_coords(surface, fes)
    old = _edge_midpoint_coords(surface, fes.ndof)
    assert new.shape == old.shape == (fes.ndof, 3)
    assert np.array_equal(new, old)


def test_dof_cluster_coords_cover_every_dof_at_higher_order(gapped_torus_vol):
    import numpy as np
    from ngsolve import HDivSurface
    from radia.bem.coil_inductance_ngsolve import _dof_cluster_coords

    surface, _ = _surface_and_volume(gapped_torus_vol)
    for order in (1, 2):
        fes = HDivSurface(surface, order=order)
        coords = _dof_cluster_coords(surface, fes)
        assert coords.shape == (fes.ndof, 3)
        assert np.isfinite(coords).all()


def test_uncompressed_volume_space_is_refused_for_the_cluster_tree(
        gapped_torus_vol):
    """An uncompressed HDivSurface on a volume mesh owns interior-edge DOFs
    that no boundary element reaches; the tree must refuse, not guess."""
    from ngsolve import HDivSurface
    from radia.bem.coil_inductance_ngsolve import _dof_cluster_coords

    _, volume = _surface_and_volume(gapped_torus_vol)
    if volume.ne == 0:
        pytest.skip("fixture is not a volume mesh")
    with pytest.raises(ValueError, match="Compress"):
        _dof_cluster_coords(volume, HDivSurface(volume, order=0))


@pytest.mark.slow
def test_hacapk_cocr_matches_cocr_at_order_one(gapped_torus_vol):
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")
    from ngsolve import TaskManager
    from radia.bem.coil_inductance_ngsolve import compute_inductance_source_sink

    surface, _ = _surface_and_volume(gapped_torus_vol)
    omega = 2 * math.pi * 150e3
    zs = _zs(SIGMA_CU, omega)
    with TaskManager():
        dense = compute_inductance_source_sink(
            surface, "source", "sink", fes_order=1, solver="cocr",
            omega=omega, Z_s_complex=zs)
        compressed = compute_inductance_source_sink(
            surface, "source", "sink", fes_order=1, solver="hacapk_cocr",
            omega=omega, Z_s_complex=zs)
    assert dense["n_J"] == compressed["n_J"]
    # The documented order-0 agreement: H-matrix matvec accuracy ~3e-7.
    assert abs(compressed["R"] / dense["R"] - 1.0) < 1e-4
    assert abs(compressed["L"] / dense["L"] - 1.0) < 1e-4
    # max|D J - g| is set by the div-free projector (sparse LU of D D^T
    # plus float accumulation over 13632 DOFs), not by the H-matrix, and
    # measured 2.0e-9 here; 1e-8 leaves room for mesh variation of the
    # auto-generated fixture without admitting a real constraint defect.
    assert float(compressed["residual"]) < 1e-8


@pytest.mark.slow
def test_volume_mesh_compresses_to_the_surface_result(gapped_torus_vol):
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")
    from ngsolve import TaskManager
    from radia.bem.coil_inductance_ngsolve import compute_inductance_source_sink

    surface, volume = _surface_and_volume(gapped_torus_vol)
    if volume.ne == 0:
        pytest.skip("fixture is not a volume mesh")
    omega = 2 * math.pi * 150e3
    zs = _zs(SIGMA_CU, omega)
    with TaskManager():
        on_surface = compute_inductance_source_sink(
            surface, "source", "sink", omega=omega, Z_s_complex=zs,
            solver="cocr")
        on_volume = compute_inductance_source_sink(
            volume, "source", "sink", omega=omega, Z_s_complex=zs,
            solver="cocr")
    # Compress strips the interior-edge DOFs: the boundary count is what the
    # extracted surface had.
    assert on_volume["n_J"] == on_surface["n_J"]
    assert on_volume["n_f"] == on_surface["n_f"]
    assert abs(on_volume["R"] / on_surface["R"] - 1.0) < 1e-8
    assert abs(on_volume["L"] / on_surface["L"] - 1.0) < 1e-8
