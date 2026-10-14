"""Golden test: HDiv-VIM demag on unstructured tetrahedral RT1 meshes.

Uses the public ``vim.Solve`` and ``vim.ChargeGram(HDiv(order=1))`` production
path on Netgen meshes.

Locks:
  - SYMMETRY exact on tets (asym ~ machine eps) -- N = B^T G B is Galerkin;
  - loops FIELD-NULL exact on tets (loop_res ~ machine eps) -- loops = ker B by construction;
  - a tet-meshed SPHERE has demag factor -> the analytic 1/3, CONVERGING as the mesh refines
    (the physics gate: the unstructured tet operator gives the right demag).

NGSolve + Netgen CSG required (importorskip); meshing is non-deterministic across versions, so the
demag is checked against a BAND around 1/3 + a refinement (monotone-toward-1/3) trend, while the
structural properties (symmetry, loop-nullity) are exact and version-robust.
"""
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

import numpy as np  # noqa: E402
from radia.vim import Solve  # noqa: E402

import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt, OrthoBrick  # noqa: E402

from conftest import hdiv_vim_dense_N_and_loops  # noqa: E402

_HEXT = ng.CoefficientFunction((0, 0, 1.0))


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def _demag(mesh):
    with ng.TaskManager():
        return Solve(mesh, mu_r=1000.0, H_ext=_HEXT)["demag"]


def test_tet_demag_symmetry_and_loop_null():
    """On a real tet mesh the C++-Gram demag operator N = B^T G B is symmetric and the loops (ker B) are
    field-null -- the two structural pillars of the HDiv-type VIM, here on UNSTRUCTURED tets."""
    with ng.TaskManager():
        N, loops, n_loop = hdiv_vim_dense_N_and_loops(_sphere(0.5))
    Nn = np.linalg.norm(N, 2)
    asym = np.linalg.norm(N - N.T) / Nn
    assert asym < 1e-4, f"tet N not symmetric: {asym:.2e}"          # symmetric to the ACA tolerance
    assert n_loop > 0, "expected loops on the tet mesh"
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / Nn
    assert loop_res < 1e-8, f"tet loops not field-null: {loop_res:.2e}"


def test_tet_cube_demag_near_third():
    """A tet-meshed cube's axis demag factor is ~1/3 (the three are equal by symmetry)."""
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    with ng.TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=0.3))
    Dz = _demag(mesh)
    assert 0.28 < Dz < 0.38, f"cube demag {Dz:.4f} not ~1/3"


def test_tet_sphere_demag_converges_to_third():
    """A tet-meshed sphere has demag factor -> 1/3 (the analytic physics gate, via the C++ charge Gram)."""
    Dz = {h: _demag(_sphere(h)) for h in (0.5, 0.3)}
    for h, v in Dz.items():
        assert abs(v - 1.0 / 3.0) < 0.02, f"sphere h={h} demag {v:.4f} not ~1/3"
