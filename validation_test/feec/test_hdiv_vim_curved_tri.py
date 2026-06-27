"""Golden: the C++ CURVED-panel surface-charge potential (rad_hdiv::CurvedTriPotential, exposed as
radia._radia_pybind._hdiv_curved_tri_potential) -- the first piece of curved-geometry support for the
HDiv-VIM charge Gram.

The flat charge Gram uses corner vertices; on a CURVED (mesh.Curve(p)) boundary the magnetic surface charge
sigma = M.n lives on the curved face, and the flat-facet approximation loses O(h^2).  The curved-panel Duffy
keeps the SAME reference-element quadrature but evaluates the curved isoparametric map X(xi) and the curved
area element J(xi)=|dX/dxi x dX/deta| at each reference Duffy point, with the singularity origin xi0 = the
reference point whose IMAGE is closest to the field point p.

This locks the C++ curved-tri potential  INT_curvedtri xi^e0 eta^e1 /|p-X(xi)| dA_curved  against an
INDEPENDENT converged plain-Gauss reference on a genuinely curved P2 triangle (bulged mid-edge nodes), at
~1e-3, for a near field point (the regime where the flat facet is wrong).
"""
import numpy as np
import pytest

rp = pytest.importorskip("radia._radia_pybind")


def _g01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def _N(xi, eta):
    L1, L2, L3 = 1 - xi - eta, xi, eta
    return np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1])


def _dN(xi, eta):
    L1 = 1 - xi - eta
    dxi = np.array([1-4*L1, 4*xi-1, 0.0, 4*(L1-xi), 4*eta, -4*eta])
    det = np.array([1-4*L1, 0.0, 4*eta-1, -4*xi, 4*xi, 4*(L1-eta)])
    return dxi, det


def _X(nodes, xi, eta):
    return _N(xi, eta) @ nodes


def _J(nodes, xi, eta):
    dxi, det = _dN(xi, eta)
    return np.linalg.norm(np.cross(dxi @ nodes, det @ nodes))


def _brute_curved(nodes, e, p, n):
    """plain collapsed-square Gauss on the reference triangle {xi+eta<=1}: xi=g[i], eta=g[j]*(1-g[i]),
    Jac (1-g[i]); curved map + curved area element.  Converges for p strictly OFF the panel."""
    g, w = _g01(n)
    s = 0.0
    for i in range(n):
        for j in range(n):
            xi = g[i]; eta = g[j] * (1 - g[i])
            X = _X(nodes, xi, eta); J = _J(nodes, xi, eta)
            r = np.linalg.norm(p - X)
            s += w[i] * w[j] * (1 - g[i]) * J * xi ** e[0] * eta ** e[1] / r
    return s


# a genuinely CURVED P2 triangle (mid-edge nodes pushed off the flat z=0 facet)
_NODES = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 1.0, 0.0],
    [0.5, -0.05, 0.30], [0.62, 0.52, 0.25], [0.08, 0.5, 0.20],
])


def test_cpp_curved_tri_potential_matches_brute():
    """C++ CurvedTriPotential == converged curved plain-Gauss reference to ~1e-3 (near field point)."""
    gx, gw = _g01(10)
    gl, gwl = gx.tolist(), gw.tolist()
    p = np.array([0.35, 0.33, 0.30])                          # near the bulged panel (flat-facet would be wrong)
    for e in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1)]:
        cpp = rp._hdiv_curved_tri_potential(_NODES.ravel().tolist(), e[0], e[1], p.tolist(), gl, gwl)
        ref = _brute_curved(_NODES, e, p, 80)
        rel = abs(cpp - ref) / (abs(ref) + 1e-30)
        assert rel < 1e-3, f"curved tri mono{e}: C++ {cpp:.6e} vs brute {ref:.6e} rel {rel:.2e}"


def test_cpp_curved_tri_far_field_high_accuracy():
    """Far from the panel the curved Duffy is high-accuracy (smooth integrand): ~1e-5 vs brute."""
    gx, gw = _g01(10)
    p = np.array([0.35, 0.33, 0.7])
    cpp = rp._hdiv_curved_tri_potential(_NODES.ravel().tolist(), 0, 0, p.tolist(), gx.tolist(), gw.tolist())
    ref = _brute_curved(_NODES, (0, 0), p, 80)
    assert abs(cpp - ref) / abs(ref) < 1e-4
