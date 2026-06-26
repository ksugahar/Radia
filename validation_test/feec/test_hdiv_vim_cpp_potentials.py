"""Golden test: the C++ analytic charge-Gram potentials (rad_hdiv::TriPotential / PhiTet, exposed
as the _hdiv_tri_potential / _hdiv_phi_tet probes) match the analytic reference to ~machine precision.

These are the building blocks: the accurate Wilton triangle 1/r surface potential and the
divergence-theorem tet Newtonian volume potential that the C++ ChargeGram H-matrix entry function uses
for near pairs.  The surface kernel is checked vs the kept radia.vim.tri_potential; the volume kernel vs
the divergence-theorem construction over tri_potential (the Python dense Gram path was removed, so the
volume reference is built inline here from the kept triangle potential).

Pure 1/r integrals (NO 1/4pi) -- the 1/(4pi) and the measure weights live in the Gram entry, not
in these geometric kernels.
"""
import numpy as np
import pytest

import radia._radia_pybind as _rp
from radia.vim._core import tri_potential


def _phi_tet(V, P):
    """Exact Newtonian potential INT_tet 1/|P-r'| dV' of a uniform tet via the divergence theorem
    (nabla'^2 R = 2/R -> (1/2) sum_{4 faces} d_face INT_face 1/R dA'), the inner face integral = the
    kept exact Wilton triangle potential tri_potential.  P: (3,) -> float."""
    V = np.asarray(V, float)
    R = np.atleast_2d(np.asarray(P, float))
    cen = V.mean(0)
    tot = np.zeros(len(R))
    for f in ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)):
        Fv = V[list(f)]
        n = np.cross(Fv[1] - Fv[0], Fv[2] - Fv[0])
        n = n / np.linalg.norm(n)
        if np.dot(Fv.mean(0) - cen, n) < 0:
            n = -n
        d = (Fv[0] - R) @ n
        tot += d * tri_potential(Fv, R)
    return float(0.5 * tot[0])

# Triangles: a flat in-plane (z=0) triangle and a tilted one (out-of-plane vertex).
_TRIS = [
    np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.9, 0.0]]),
    np.array([[0.0, 0.0, 0.0], [0.7, 0.1, 0.0], [0.1, 0.8, 0.05]]),
]
# Obs points: above the plane, ON the plane (z=0, the Wilton |z|->0 edge term), far, and below.
_TRI_OBS = [
    [0.3, 0.3, 0.5], [0.5, 0.4, 0.0], [2.0, 1.0, 0.3], [0.3, 0.3, -0.7],
]

# A non-degenerate tetrahedron and obs points outside (above/lateral/below) and inside.
_TET = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.9, 0.0], [0.3, 0.2, 1.1]])
_TET_OBS = [
    [0.3, 0.3, 1.6], [1.6, 0.3, 0.4], [0.3, 0.3, -0.5], [0.5, 0.5, 0.5],
]


def test_cpp_tripotential_matches_python():
    """C++ TriPotential == Python tri_potential to machine precision (incl. the on-plane z=0 case)."""
    maxrel = 0.0
    for V in _TRIS:
        for r in _TRI_OBS:
            c = _rp._hdiv_tri_potential(V.flatten().tolist(), [float(x) for x in r])
            p = float(tri_potential(V, np.asarray(r, float)))
            assert p > 0.0  # a 1/r integral over a real triangle is strictly positive
            maxrel = max(maxrel, abs(c - p) / abs(p))
    assert maxrel < 1e-12, f"C++ TriPotential drifted from Python: max rel = {maxrel:.3e}"


def test_cpp_phitet_matches_python():
    """C++ PhiTet == the divergence-theorem reference to machine precision (obs outside AND inside)."""
    maxrel = 0.0
    for P in _TET_OBS:
        c = _rp._hdiv_phi_tet(_TET.flatten().tolist(), [float(x) for x in P])
        p = _phi_tet(_TET, np.asarray(P, float))
        assert p > 0.0
        maxrel = max(maxrel, abs(c - p) / abs(p))
    assert maxrel < 1e-12, f"C++ PhiTet drifted from the reference: max rel = {maxrel:.3e}"


def test_cpp_tripotential_translation_invariant():
    """The 1/r triangle potential is invariant under a rigid translation of (triangle, obs)."""
    V, r = _TRIS[0], np.asarray(_TRI_OBS[0], float)
    shift = np.array([3.1, -2.4, 1.7])
    a = _rp._hdiv_tri_potential(V.flatten().tolist(), [float(x) for x in r])
    b = _rp._hdiv_tri_potential((V + shift).flatten().tolist(), [float(x) for x in (r + shift)])
    assert abs(a - b) / abs(a) < 1e-12
