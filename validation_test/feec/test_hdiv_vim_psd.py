"""Golden: the HDiv-VIM charge Gram N = B^T G B is POSITIVE-SEMIDEFINITE at every order.

N is a demag SELF-ENERGY (1/2 m^T N m, the magnetic charge interaction energy), so it MUST be PSD -- a
negative eigenvalue is unphysical and silently corrupts any ENERGY / EIGENVALUE use of N (e.g. the nonlinear
energy-minimisation solve slides down the spurious negative mode -> garbage high-order magnetisation).

This guards a regression the demag-FACTOR golden could not catch: the factor only probes a couple of specific
M, which happened to land on N's positive spectrum, so an under-integrated (non-PSD) Gram passed it.  Here we
check the FULL spectrum.  The near/self quadrature floor in ChargeGram is quad >= 3*p exactly so that N
stays PSD (quad=4,5 give min eig ~ -3.2e-4 at p=2 -- NON-PSD; quad>=6 PSD).

NGSolve + Netgen required (importorskip).  N is formed densely on a small mesh and its symmetric-part spectrum
is checked; PSD is a quadrature-accuracy property (mesh-robust), so a coarse sphere suffices.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Sphere, OCCGeometry, Pnt  # noqa: E402

from radia.vim import ChargeGram  # noqa: E402


def _N_dense(fes):
    """N = B^T G B as a dense symmetric matrix (default near/self quad = 3*p -> PSD)."""
    B, G, _M = ChargeGram(fes)                       # default intorder -> quad = max(3p, 4)
    n = fes.ndof
    N = np.zeros((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        N[:, j] = B.T @ np.asarray(G.matvec((B @ e).tolist()), float)
        e[j] = 0.0
    return 0.5 * (N + N.T)                                   # symmetric part (the physical energy operator)


@pytest.mark.parametrize("order", [1, 2])
def test_charge_gram_is_psd(order):
    """The production pure-TET RT1/RT2 charge Gram is positive semidefinite."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=1.0))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        w = np.linalg.eigvalsh(_N_dense(fes))
    assert w.min() >= -1e-9 * w.max(), (
        f"charge Gram NOT PSD at order={order}: min eig {w.min():.3e} vs max {w.max():.3e} "
        f"(near/self quad floor too low -- needs quad >= 3*p)")
