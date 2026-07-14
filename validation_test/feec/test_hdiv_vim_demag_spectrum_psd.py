"""Golden: the high-order HDiv-VIM demag operator N = B^T G B is a VALID demag operator -- its generalized
spectrum eig(M_mass^{-1} N) lies in [0, 1] (a true magnetostatic demag factor is bounded: INT_all |H_demag|^2
<= INT_body |M|^2).

This locks the PER-ELEMENT change-of-basis fix (radia.vim._vim._charge_basis): the old code reused a single
element-0 L2->monomial change-of-basis S for ALL elements, but NGSolve orients the (high-order) shape
functions + the element map by each element's GLOBAL vertex order, so a single S scrambles the monomial
charge of differently-oriented elements -> a SPURIOUS NET CHARGE (broken neutrality) -> the demag operator
got UNPHYSICAL eigenvalues > 1 (p=1 ~1.2, p=2 ~8) and the high-order MATERIAL solve over-magnetized ~2x/~4x.
Per-element S restores neutrality and the [0,1] spectrum.  Invisible to the uniform-M demag-FACTOR test
(constant charge -> the monopole error cancels), so this dedicated spectrum gate is needed.
See memory hdiv-highorder-material-solve-wrong round 4.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import scipy.sparse as sp        # noqa: E402
import scipy.linalg as sla       # noqa: E402
import ngsolve as ng             # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt   # noqa: E402
from radia.vim import ChargeGram    # noqa: E402


def _dense_N(B, G):
    B = sp.csr_matrix(B)
    nd = B.shape[1]
    N = np.zeros((nd, nd))
    for k in range(nd):
        e = np.zeros(nd); e[k] = 1.0
        N[:, k] = B.T @ np.array(G.matvec((B @ e).tolist()))
    return 0.5 * (N + N.T)


@pytest.mark.parametrize("p", [1])
def test_highorder_demag_spectrum_in_unit_interval(p):
    """eig(M_mass^{-1} N) in [0, 1] for the order-p charge-Gram demag operator (was max ~1.2 (p1) with the
    single-element-0 change-of-basis bug).  The public production path exercises p=1."""
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.8))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=p)
        B, G, M_mass = ChargeGram(fes, ho_far_factor=float("inf"))   # exact Gram
        N = _dense_N(B, G)
        ev = sla.eigh(N, sp.csr_matrix(M_mass).toarray(), eigvals_only=True)
    # a true demag spectrum is in [0,1]; allow a small PSD/quadrature margin
    assert ev.min() > -1e-3, f"p={p}: demag operator not PSD (min eig {ev.min():.3e})"
    assert ev.max() < 1.0 + 5e-2, \
        f"p={p}: UNPHYSICAL demag eig {ev.max():.4f} > 1 -- the per-element change-of-basis regressed " \
        f"(spurious net charge / broken neutrality)"


def test_rt1_material_solve_is_consistent_fixed_point():
    """The RT1 uniform-field cube material solve sits at the consistent linear fixed point
    M_avg ~ chi*H0/(1+chi*D) (D = the demag factor = Rayleigh quotient of N for uniform M) -- NOT the old
    ~1.8x blow-up.  The reference is the analytic fixed point.
    (The strong M4 guard is the eig-in-[0,1] spectrum test above; this is the end-to-end material-solve check.)"""
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.7))
    H0, mu_r = 1000.0, 100.0
    chi = mu_r - 1.0
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = ChargeGram(fes, ho_far_factor=float("inf"))
        N = _dense_N(B, G)
        Md = sp.csr_matrix(M_mass).toarray()
        gfH = ng.GridFunction(fes); gfH.Set(ng.CoefficientFunction((0, 0, H0)))
        h = gfH.vec.FV().NumPy().copy()
        m = np.linalg.solve((1.0 / chi) * Md + N, Md @ h)
        gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
        vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
        mavg = ng.Integrate(gfM[2], mesh) / vol
        gfu = ng.GridFunction(fes); gfu.Set(ng.CoefficientFunction((0, 0, 1)))
        mu = gfu.vec.FV().NumPy().copy()
        D = float((mu @ (N @ mu)) / (mu @ (Md @ mu)))      # demag factor (uniform-M Rayleigh quotient)
    m_fixed = chi * H0 / (1.0 + chi * D)
    rel = abs(mavg - m_fixed) / abs(m_fixed)
    assert rel < 0.25, \
        f"RT1 M_avg {mavg:.0f} vs consistent fixed point {m_fixed:.0f} (D={D:.3f}, rel {rel:.2f}) -- M4 regressed"
