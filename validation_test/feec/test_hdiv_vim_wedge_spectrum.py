"""Golden: the WEDGE (PRISM) RT1 HDiv-VIM demag operator N = B^T G B is a VALID demag operator -- its
generalized spectrum eig(M_mass^{-1} N) lies in [0, 1] (a true magnetostatic demag factor is bounded).

This locks the C++ wedge-mode charge Gram (2026-07-04, memory hdiv-tet-hex-coupling-pyramid-gated): the
prism div-image is L2(prism,order=1) = tri-P1 (x) z-P1 = {1,x,y,z,xz,yz} (6/prism, a subset of the hex's 8
Q1 monomials); boundary faces are MIXED tri (SurfaceL2 P1) + quad (SurfaceL2 Q1); geometry is the 18-node
tri-P2 (x) z-P2 lattice.  The wedge mode mirrors the hex mode's both-domains-graded Duffy singular
quadrature on a 3-sub-tet / mixed 1-2 sub-tri decomposition.  numpy de-risk eig: 0.989/0.997 @ n=2/3; the
C++ reproduces it (0.992/0.998).  A prism-meshed cube's demag factor is the exact 1/3.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
import scipy.sparse as sp        # noqa: E402
import scipy.linalg as sla       # noqa: E402
import ngsolve as ng             # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh   # noqa: E402
from radia.vim import ChargeGram        # noqa: E402

_L = 0.02
_MP = lambda x, y, z: (_L * (x - 0.5), _L * (y - 0.5), _L * (z - 0.5))   # noqa: E731


def _prism_cube(n):
    try:
        mesh = MakeStructured3DMesh(prism=True, nx=n, ny=n, nz=n, mapping=_MP)
    except TypeError:
        pytest.skip("this NGSolve MakeStructured3DMesh has no prism= kwarg")
    if {len(el.vertices) for el in mesh.Elements(ng.VOL)} != {6}:
        pytest.skip("prism mesh did not produce pure wedges")
    return mesh


def _dense_N(B, G):
    B = sp.csr_matrix(B)
    nd = B.shape[1]
    N = np.zeros((nd, nd))
    for k in range(nd):
        e = np.zeros(nd); e[k] = 1.0
        N[:, k] = B.T @ np.array(G.matvec((B @ e).tolist()))
    return 0.5 * (N + N.T)


@pytest.mark.parametrize("n", [2, 3])
def test_wedge_demag_spectrum_in_unit_interval(n):
    """eig(M_mass^{-1} N) in [0, 1] for the wedge-mode charge-Gram demag operator on a prism cube."""
    mesh = _prism_cube(n)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = ChargeGram(fes)
        N = _dense_N(B, G)
        ev = sla.eigh(N, sp.csr_matrix(M_mass).toarray(), eigvals_only=True)
    assert ev.min() > -1e-6, f"n={n}: wedge demag operator not PSD (min eig {ev.min():.3e})"
    assert ev.max() < 1.0 + 2e-2, \
        f"n={n}: UNPHYSICAL wedge demag eig {ev.max():.4f} > 1 (the wedge charge-Gram regressed)"


def test_wedge_cube_demag_factor_is_one_third():
    """A prism-meshed cube's uniform-M demag factor D = Rayleigh quotient of N is the exact cube value 1/3."""
    mesh = _prism_cube(3)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = ChargeGram(fes)
        N = _dense_N(B, G)
        Md = sp.csr_matrix(M_mass).toarray()
        gfu = ng.GridFunction(fes); gfu.Set(ng.CoefficientFunction((0, 0, 1)))
        mu = gfu.vec.FV().NumPy().copy()
        D = float((mu @ (N @ mu)) / (mu @ (Md @ mu)))
    assert abs(D - 1.0 / 3.0) < 5e-3, f"wedge cube demag {D:.4f} off 1/3"
