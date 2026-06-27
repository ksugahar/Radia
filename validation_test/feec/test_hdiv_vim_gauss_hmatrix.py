"""Gauss-point H-matrix ChargeGram path.

This locks the first production slice of

    G_charge ~= P^T K_point P + G_near_correction

where K_point is a HACApK H-matrix over quadrature points with the cheap
Laplace 1/r kernel, and the sparse near correction restores analytic
self/near charge Gram entries.  The path is opt-in while the larger
high-order/curved validation is built out.
"""

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import hdiv_demag_solve  # noqa: E402
from radia.vim._core import build_demag  # noqa: E402


def _sphere(h=0.7):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


def test_gauss_point_hmatrix_linear_sphere_matches_analytic_chargegram():
    mesh = _sphere()
    Hext = ng.CoefficientFunction((0, 0, 1000.0))
    with ng.TaskManager():
        ref = hdiv_demag_solve(mesh, 100.0, Hext, gram_eps=1e-6,
                               near_factor=2.0, gram_backend="analytic")
        got = hdiv_demag_solve(mesh, 100.0, Hext, gram_eps=1e-6,
                               near_factor=2.0, gram_backend="gauss", gauss_near_factor=2.0)
    assert got["gram_backend"] == "gauss"
    assert got["linear_solver"] == "mass-riesz-cg"   # default 'auto' = mass-riesz CG (gauss op's own C++ CG)
    assert got["hmat_stats"]["n_dof"] > got["n_charge"]   # point H-matrix, not charge-entry H-matrix
    rel_m = abs(got["M_avg"][2] - ref["M_avg"][2]) / abs(ref["M_avg"][2])
    assert rel_m < 5e-3
    assert abs(got["demag"] - ref["demag"]) < 2e-3


def test_analytic_entry_oracle_build_false_matches_built_entries():
    """The geometry-only (build=False) analytic charge-Gram is an exact ENTRY ORACLE: .entry() agrees
    bit-for-bit with the fully built analytic H-matrix.  This is the contract the Gauss near correction
    relies on -- it samples exact analytic near entries WITHOUT paying the full H-matrix build."""
    mesh = _sphere()
    with ng.TaskManager():
        d = build_demag(mesh)
    kw = dict(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
              n_el=int(d["n_el"]), near_factor=1e30)
    built = _rp._ChargeGramHMatrix(**kw, build=True)
    oracle = _rp._ChargeGramHMatrix(**kw, build=False)
    n = int(d["n_charge"])
    probe = [(0, 0), (0, 1), (1, 0), (n - 1, n - 1), (0, n - 1), (n // 2, n // 3)]
    for a, b in probe:
        assert oracle.entry(a, b) == built.entry(a, b)

