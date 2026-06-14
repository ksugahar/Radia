"""Golden test: HDiv-type VIM demag on UNSTRUCTURED TET meshes (production step #1a -- the critical
path from synthetic structured hex grids to REAL geometry).  Uses NGSolve HDiv(order=0) on Netgen
tet meshes (element-agnostic extraction: B = -div [volume charge] + Trace.n [surface charge];
M_mass = HDiv mass) with the tet/triangle Gram self-energy.

Locks:
  - SYMMETRY exact on tets (asym ~ machine eps) -- N = B^T G B is Galerkin;
  - loops FIELD-NULL exact on tets (loop_res ~ machine eps) -- loops = ker B by construction;
  - a tet-meshed SPHERE has demag factor -> the analytic 1/3, CONVERGING as the mesh refines
    (the physics gate: the unstructured tet operator gives the right demag).

NGSolve + Netgen CSG required (importorskip); meshing is non-deterministic across versions, so the
demag is checked against a BAND around 1/3 + a refinement (monotone-toward-1/3) trend, while the
structural properties (symmetry, loop-nullity) are exact and version-robust.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _core as tet  # noqa: E402

import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt, OrthoBrick  # noqa: E402


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_tet_demag_symmetry_and_loop_null():
    """On a real tet mesh the demag operator is symmetric and the loops are field-null -- the two
    structural pillars of the HDiv-type VIM, here on UNSTRUCTURED tets (not the structured hex grid)."""
    d = tet.build_demag(_sphere(0.5), nsub=4)
    import numpy as np
    Nn = np.linalg.norm(d["N"], 2)
    asym = np.linalg.norm(d["N"] - d["N"].T) / Nn
    assert asym < 1e-12, f"tet N not symmetric: {asym:.2e}"
    assert d["n_loop"] > 0, "expected loops on the tet mesh"
    loop_res = max(np.linalg.norm(d["N"] @ d["loops"][k]) for k in range(d["n_loop"])) / Nn
    assert loop_res < 1e-10, f"tet loops not field-null: {loop_res:.2e}"


def test_tet_cube_demag_near_third():
    """A tet-meshed cube's axis demag factor is ~1/3 (the three are equal by symmetry)."""
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    d = tet.build_demag(ng.Mesh(geo.GenerateMesh(maxh=0.3)), nsub=4)
    Dz = tet.demag_factor(d)
    assert 0.28 < Dz < 0.35, f"cube demag {Dz:.4f} not ~1/3"


def test_tet_sphere_demag_converges_to_third():
    """A tet-meshed sphere has demag factor -> 1/3, converging (increasing toward 1/3) as h refines.
    This is the analytic physics gate for the unstructured tet operator."""
    Dz = {h: tet.demag_factor(tet.build_demag(_sphere(h), nsub=4)) for h in (0.5, 0.3)}
    for h, v in Dz.items():
        assert 0.29 < v < 0.3333 + 1e-3, f"sphere h={h} demag {v:.4f} not approaching 1/3 from below"
    assert Dz[0.3] > Dz[0.5] - 1e-3, f"demag not converging toward 1/3 with refinement: {Dz}"
