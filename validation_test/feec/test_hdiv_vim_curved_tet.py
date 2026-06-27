"""Golden: the C++ CURVED-panel VOLUME-charge potential (rad_hdiv::CurvedTetPotential, exposed as
radia._radia_pybind._hdiv_curved_tet_potential) -- the volume-charge half of curved-geometry support for the
HDiv-VIM charge Gram (companion to test_hdiv_vim_curved_tri).

On a curved (mesh.Curve(p)) mesh the volume charge rho = -div M lives on the curved tetrahedron; the
curved-panel Duffy keeps the reference-element quadrature but evaluates the curved P2 map X(xi) and the curved
VOLUME element Jv = |det dX/dxi| at each reference Duffy point, with the singularity origin xi0 = the
reference point whose image is closest to p.  The reference tet is always +oriented and Jv=|det| -> no
host-orientation sign issue (unlike the flat physical tet Duffy, which needed sign(host vol)).

Locks the C++ curved-tet potential against an independent converged plain-Gauss reference on a genuinely
curved P2 tet (bulged mid-edge nodes) at ~1e-3, for a near field point.
"""
import numpy as np
import pytest

rp = pytest.importorskip("radia._radia_pybind")


def _g01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def _N(xi, eta, zeta):
    L1, L2, L3, L4 = 1 - xi - eta - zeta, xi, eta, zeta
    return np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), L4*(2*L4-1),
                     4*L1*L2, 4*L2*L3, 4*L3*L1, 4*L1*L4, 4*L2*L4, 4*L3*L4])


def _dN(xi, eta, zeta):
    L1, L2, L3, L4 = 1 - xi - eta - zeta, xi, eta, zeta
    dL = [np.array([-1., -1, -1]), np.array([1., 0, 0]), np.array([0., 1, 0]), np.array([0., 0, 1])]
    d = np.zeros((10, 3))
    for k in range(3):
        d[0, k] = (4*L1-1)*dL[0][k]; d[1, k] = (4*L2-1)*dL[1][k]; d[2, k] = (4*L3-1)*dL[2][k]; d[3, k] = (4*L4-1)*dL[3][k]
        d[4, k] = 4*(dL[0][k]*L2+L1*dL[1][k]); d[5, k] = 4*(dL[1][k]*L3+L2*dL[2][k]); d[6, k] = 4*(dL[2][k]*L1+L3*dL[0][k])
        d[7, k] = 4*(dL[0][k]*L4+L1*dL[3][k]); d[8, k] = 4*(dL[1][k]*L4+L2*dL[3][k]); d[9, k] = 4*(dL[2][k]*L4+L3*dL[3][k])
    return d


def _X(nodes, xi, eta, zeta):
    return _N(xi, eta, zeta) @ nodes


def _Jv(nodes, xi, eta, zeta):
    Jac = (_dN(xi, eta, zeta).T @ nodes).T            # Jac[k][c] = dX_k/dxi_c
    return abs(np.linalg.det(Jac))


def _brute_tet(nodes, e, p, n):
    """plain Gauss on the reference tet {xi+eta+zeta<=1} (collapsed cube), curved map + curved volume element.
    Converges for p strictly OUTSIDE the curved tet."""
    g, w = _g01(n)
    s = 0.0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                xi = g[i]; eta = g[j]*(1-xi); zeta = g[k]*(1-xi-eta)
                jac = (1-xi)**2 * (1-g[j])
                X = _X(nodes, xi, eta, zeta); Jv = _Jv(nodes, xi, eta, zeta)
                r = np.linalg.norm(p - X)
                s += w[i]*w[j]*w[k]*jac*Jv*xi**e[0]*eta**e[1]*zeta**e[2]/r
    return s


_NODES = np.array([
    [0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1],
    [0.5, -0.05, 0.04], [0.55, 0.55, -0.04], [-0.04, 0.5, 0.05],
    [0.03, -0.03, 0.5], [0.52, 0.02, 0.52], [0.04, 0.5, 0.55],
])


def test_cpp_curved_tet_potential_matches_brute():
    """C++ CurvedTetPotential == converged curved plain-Gauss reference to ~1e-3 (near field point)."""
    gx, gw = _g01(10)
    gl, gwl = gx.tolist(), gw.tolist()
    p = np.array([0.25, 0.22, 0.22])                          # just outside the curved tet
    for e in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (1, 1, 0)]:
        cpp = rp._hdiv_curved_tet_potential(_NODES.ravel().tolist(), e[0], e[1], e[2], p.tolist(), gl, gwl)
        ref = _brute_tet(_NODES, e, p, 44)
        rel = abs(cpp - ref) / (abs(ref) + 1e-30)
        assert rel < 1e-3, f"curved tet mono{e}: C++ {cpp:.6e} vs brute {ref:.6e} rel {rel:.2e}"
