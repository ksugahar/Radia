"""Golden test: the NEAR/FAR split of the analytic charge-Gram H-matrix (the build speedup).

The analytic entry is EXPENSIVE (PhiTet/TriPotential per outer point); the near/far split
(_ChargeGramHMatrix(..., near_factor=2.0)) uses it only for NEAR pairs and the cheap centroid-monopole
for FAR pairs.  Physically: the analytic entry matters for the NEAR (non-uniform-M, div M != 0)
interaction; the far field is monopole to O((size/r)^2).  The dense Python Gram path was removed, so the
reference here is the C++ ALL-ANALYTIC H-matrix (near_factor=1e30): the near/far split demag must match
the full-analytic C++ demag on BOTH a uniform body (sphere) AND a non-uniform body (cube).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import _core as tet  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402


def _sphere(h=0.5):
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


def _cube(h=0.5):
    g = CSGeometry(); g.Add(OrthoBrick(Pnt(-.5, -.5, -.5), Pnt(.5, .5, .5)))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


def _demag(d, near_factor):
    H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                               n_el=int(d["n_el"]), eps=1e-6, leaf=32, eta=2.0, near_factor=near_factor)
    B = d["B_csr"]; m = d["m_unit"]
    Nm = B.T @ np.asarray(H.matvec((B @ m).tolist()), float)
    return float((m @ Nm) / (m @ (d["M_mass"] @ m)))


def test_nearfar_split_sphere():
    d = tet.build_demag(_sphere())
    Dz_full = _demag(d, near_factor=1e30)               # C++ FULL-analytic reference
    Dz_split = _demag(d, near_factor=2.0)
    assert 0.30 < Dz_split < 0.36, f"split sphere demag {Dz_split:.4f} not ~1/3"
    assert abs(Dz_split - Dz_full) < 3e-2 * abs(Dz_full), \
        f"near/far split demag {Dz_split:.4f} drifts from full analytic {Dz_full:.4f}"


def test_nearfar_split_cube_nonuniform():
    """The non-uniform body: the split must keep the analytic accuracy on the cube (div M != 0)."""
    d = tet.build_demag(_cube())
    Dz_full = _demag(d, near_factor=1e30)
    Dz_split = _demag(d, near_factor=2.0)
    assert 0.0 < Dz_split < 1.0
    assert abs(Dz_split - Dz_full) < 3e-2 * abs(Dz_full), \
        f"near/far split cube demag {Dz_split:.4f} drifts from full analytic {Dz_full:.4f}"
