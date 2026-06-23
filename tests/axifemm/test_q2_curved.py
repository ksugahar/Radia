# -*- coding: utf-8 -*-
r"""Regression test for AxiHenrotteFE_Q2_Curved -- the 9-node curved biquadratic
Henrotte quad element (Phase B6 C++ port of axifemm_quad_q2_curved.py).

Selected via H1Henrotte(..., order=2, curvedquad=True): order-2 QUADS sample 9
curved node positions from the (Curve(2)-d) ElementTransformation and use the
isoparametric element, so general / curved-boundary quads work where the
axis-aligned closed-form Q2 (curvedquad=False, the default) is invalid.

Gates (eddy-decay eigenvalue K v = lambda M v, V-DOF stiffness + sigma-mass):
  1. STRAIGHT axis-aligned quads: curvedquad=True must reproduce the validated
     axis-aligned closed form to quadrature precision (8x8 Gauss is exact for the
     polynomial integrand).
  2. GENERAL (annular-sector, skewed) quads: curvedquad=True converges; the
     axis-aligned form is built for rectangles and is wrong / NaN there.

Requires the rebuilt axifem extension (Build.ps1 -AxiFemmOnly).  The repo build
copies axifem.pyd to src/radia, which this test loads explicitly.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from ngsolve import BilinearForm, CoefficientFunction as CF, TaskManager
from ngsolve.meshes import MakeStructured2DMesh
from radia.axifem import (H1Henrotte, AxiHenrotteStiffnessBFI,
                          AxiHenrotteSigmaMassBFI)


def _lam_min(mesh, curvedquad):
    fes = H1Henrotte(mesh, order=2, dirichlet=".*", curvedquad=curvedquad)
    aK = BilinearForm(fes, symmetric=True)
    aK += AxiHenrotteStiffnessBFI(CF(1.0))
    aM = BilinearForm(fes, symmetric=True)
    aM += AxiHenrotteSigmaMassBFI(CF(1.0))
    with TaskManager():
        aK.Assemble()
        aM.Assemble()
    n = fes.ndof
    free = np.array([i for i in range(n) if fes.FreeDofs()[i]], dtype=int)

    def dense(mat):
        rr, cc, vv = mat.COO()
        A = sp.csr_matrix((np.asarray(vv, float),
                           (np.asarray(rr, int), np.asarray(cc, int))),
                          shape=(n, n)).toarray()[np.ix_(free, free)]
        return 0.5 * (A + A.T)

    w = eigh(dense(aK.mat), dense(aM.mat), eigvals_only=True, subset_by_index=[0, 0])
    return float(w[0])


def _rect(nx, ny):
    return MakeStructured2DMesh(quads=True, nx=nx, ny=ny,
                                mapping=lambda x, y: (0.1 + 0.9 * x, -0.5 + y))


def _annulus(nrho, nphi, R1=1.0, R2=2.0, dphi=math.radians(30.0)):
    def mp(x, y):
        rho = R1 + (R2 - R1) * x
        phi = -dphi + 2 * dphi * y
        return (rho * math.cos(phi), rho * math.sin(phi))
    return MakeStructured2DMesh(quads=True, nx=nrho, ny=nphi, mapping=mp)


def test_curved_matches_axis_aligned_on_straight_quads():
    for nx in (4, 6, 8):
        m = _rect(nx, nx)
        la = _lam_min(m, False)   # axis-aligned closed form
        lc = _lam_min(m, True)    # curved isoparametric
        rel = abs(lc - la) / abs(la)
        print(f"  straight nx={nx}: axis={la:.8f} curved={lc:.8f} rel={rel:.2e}")
        assert rel < 1e-9, f"curved Q2 != axis-aligned on straight quad (rel {rel:.2e})"
    print("[OK] Q2_Curved reproduces the axis-aligned closed form on straight quads")


def test_curved_converges_on_skewed_annular_quads():
    lams = []
    for (nr, npi) in [(3, 6), (4, 8), (6, 12), (8, 16)]:
        lc = _lam_min(_annulus(nr, npi), True)
        assert np.isfinite(lc), f"curved Q2 NaN on annular {nr}x{npi}"
        lams.append(lc)
        print(f"  annular {nr}x{npi}: curved lam_min={lc:.6f}")
    # monotone-ish convergence: successive changes shrink and stay bounded
    d1 = abs(lams[1] - lams[0]); d3 = abs(lams[3] - lams[2])
    assert d3 < d1, f"curved Q2 not converging on skewed quads ({d1:.3e} -> {d3:.3e})"
    print("[OK] Q2_Curved converges on general (skewed annular) quads")


if __name__ == "__main__":
    test_curved_matches_axis_aligned_on_straight_quads()
    test_curved_converges_on_skewed_annular_quads()
    print("\nAll AxiHenrotteFE_Q2_Curved tests passed.")
