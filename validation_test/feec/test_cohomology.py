"""Golden: radia.cohomology -- pure-Python (gmsh-free) cohomology of a multiply-connected mesh.

Locks the topology engine that replaced the Gmsh computeHomology dependency:
  * b1 (first Betti number) = 0 / 1 / 2 on a solid cylinder / genus-1 washer / two-hole plate;
  * the realised HCurl(order=0) cut basis is genuinely CURL-FREE and has UNIT CIRCULATION around the hole
    (and ~0 around a contractible loop) -- i.e. a valid T-Omega "cut" / one-ampere-turn loop field;
  * for two holes the period matrix is NON-SINGULAR (the two cuts separate the two loops);
  * the gmsh-free CohomologyCutSolver.setup_from_mesh path runs and returns b1 (no gmsh import).

NGSolve + Netgen required (importorskip).  b1 is an integer topological invariant (mesh-robust); the
circulation / determinant checks use bands because meshing is non-deterministic.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, Cylinder, Pnt, Dir, OCCGeometry  # noqa: E402

from radia.cohomology import betti_numbers, cohomology_basis, circulation  # noqa: E402

_AX = Dir(0, 0, 1)


def _cyl(z0, r, h):
    return Cylinder(Pnt(0, 0, z0), _AX, r=r, h=h)


def _mesh(shape, h=0.02):
    return ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=h))


def _solid():
    return _mesh(_cyl(-0.015, 0.05, 0.03))


def _washer():
    return _mesh(_cyl(-0.015, 0.05, 0.03) - _cyl(-0.025, 0.02, 0.05))


def _two_hole():
    box = Box(Pnt(-0.08, -0.04, -0.015), Pnt(0.08, 0.04, 0.015))
    return _mesh(box - _cyl(-0.025, 0.015, 0.05).Move((-0.04, 0, 0))
                 - _cyl(-0.025, 0.015, 0.05).Move((0.04, 0, 0)))


def test_betti1_solid_washer_two_hole():
    """b1 = 0 / 1 / 2 for solid / washer (genus 1) / two-hole plate -- exact topological invariant."""
    with ng.TaskManager():
        assert betti_numbers(_solid()) == (1, 0)
        assert betti_numbers(_washer()) == (1, 1)
        assert betti_numbers(_two_hole()) == (1, 2)


def test_washer_cut_basis_curlfree_unit_circulation():
    """The washer cut function h_0 is curl-free, has UNIT circulation around the hole, ~0 on a contractible
    loop -- a genuine cohomology generator (one ampere-turn loop field)."""
    mesh = _washer()
    with ng.TaskManager():
        basis, b1, fes, _ctx, _loops = cohomology_basis(mesh)
        assert b1 == 1
        h = basis[0]
        curl_rel = (np.sqrt(ng.Integrate(ng.InnerProduct(ng.curl(h), ng.curl(h)), mesh))
                    / np.sqrt(ng.Integrate(ng.InnerProduct(h, h), mesh)))
        c_hole = circulation(h, mesh, 0.0, 0.0, 0.035)
        c_contractible = circulation(h, mesh, 0.035, 0.0, 0.006)
    assert curl_rel < 1e-6, f"cut basis not curl-free: ||curl h||/||h|| = {curl_rel:.2e}"
    assert abs(abs(c_hole) - 1.0) < 0.05, f"circulation around the hole not unit: {c_hole:.4f}"
    assert abs(c_contractible) < 1e-2, f"circulation on a contractible loop not ~0: {c_contractible:.2e}"


def test_two_hole_period_matrix_nonsingular():
    """Two holes -> the 2x2 period matrix (circulation of each cut around each hole) is NON-SINGULAR:
    the cohomology basis spans H^1 and distinguishes the two loops."""
    mesh = _two_hole()
    with ng.TaskManager():
        basis, b1, _fes, _ctx, _loops = cohomology_basis(mesh)
        assert b1 == 2
        P = np.array([[circulation(basis[k], mesh, cx, 0.0, 0.026) for cx in (-0.04, 0.04)]
                      for k in range(2)])
    assert abs(np.linalg.det(P)) > 0.1, f"period matrix singular (holes not separated):\n{P}"


def test_cutsolver_setup_from_mesh_is_gmsh_free():
    """CohomologyCutSolver.setup_from_mesh (the gmsh-free path) computes the cut basis from an NGSolve mesh
    and reports b1 -- no gmsh import, no .msh transfer."""
    from radia.cohomology_cut import CohomologyCutSolver

    solver = CohomologyCutSolver()
    with ng.TaskManager():
        n = solver.setup_from_mesh(_washer())
    assert n == 1, f"expected 1 coil loop on the washer, got {n}"
    assert len(solver.get_cohomology_basis()) == 1
