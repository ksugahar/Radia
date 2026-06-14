"""Golden test: the C++ analytic charge-Gram potentials (rad_hdiv::TriPotential / PhiTet, exposed
as the _hdiv_tri_potential / _hdiv_phi_tet probes) match the dense Python reference
(radia.vim._core.tri_potential / phi_tet) to ~machine precision.

These are the M2 building blocks: the accurate Wilton triangle 1/r surface potential and the
divergence-theorem tet Newtonian volume potential that the C++ scalable ChargeGram H-matrix entry
function will use for near pairs (replacing the centroid-monopole approximation).  Because the C++
is a line-by-line port of the validated Python Gram, this test pins that they stay bit-identical:
any future edit to either side that drifts the formula is caught here, not at the (much coarser)
H-matrix-vs-dense level.

Pure 1/r integrals (NO 1/4pi) -- the 1/(4pi) and the measure weights live in the Gram entry, not
in these geometric kernels.
"""
import numpy as np
import pytest

import radia._radia_pybind as _rp
from radia.vim._core import tri_potential, phi_tet

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
    """C++ PhiTet == Python phi_tet to machine precision (obs outside AND inside the tet)."""
    maxrel = 0.0
    for P in _TET_OBS:
        c = _rp._hdiv_phi_tet(_TET.flatten().tolist(), [float(x) for x in P])
        p = float(phi_tet(_TET, np.asarray(P, float)))
        assert p > 0.0
        maxrel = max(maxrel, abs(c - p) / abs(p))
    assert maxrel < 1e-12, f"C++ PhiTet drifted from Python: max rel = {maxrel:.3e}"


def test_cpp_tripotential_translation_invariant():
    """The 1/r triangle potential is invariant under a rigid translation of (triangle, obs)."""
    V, r = _TRIS[0], np.asarray(_TRI_OBS[0], float)
    shift = np.array([3.1, -2.4, 1.7])
    a = _rp._hdiv_tri_potential(V.flatten().tolist(), [float(x) for x in r])
    b = _rp._hdiv_tri_potential((V + shift).flatten().tolist(), [float(x) for x in (r + shift)])
    assert abs(a - b) / abs(a) < 1e-12
