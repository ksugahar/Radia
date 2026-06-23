# -*- coding: utf-8 -*-
r"""Regression test: P2 triangle (order=2) is full-rank on an axis-touching mesh.

The AxiHenrotteFESpace reserves a per-VOL-element face-center DOF (nv+ne+ei.Nr())
for Q2 quads, but P2 TRIANGLES never use theirs.  Before the fix those dead slots
stayed in the free-dof set, leaving the custom-BFI K and M rank-deficient (a zero
eigenvalue) so the eddy generalized eigenvalue solve  K v = lambda M v  failed on
any triangle mesh -- even with correct axis|outer Dirichlet.  The FESpace now marks
the trig-center slots UNUSED (FinalizeUpdate), restoring full rank.

Gates:
  1. Every P2-triangle face-center slot is UNUSED and not free.
  2. The custom-BFI K, M assemble full-rank so the eddy generalized eigenvalue
     K v = lambda M v solves (a zero eigenvalue would make M not positive-definite
     and eigh raise -- the pre-fix failure), and P1 and P2 agree on the slowest
     mode (P2 is physically consistent, not a spurious null space).

Requires the rebuilt axifem extension (Build.ps1 -AxiFemOnly); loads src/radia.
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
from ngsolve import BilinearForm, CoefficientFunction as CF, Mesh
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, X
from radia.axifem import (H1Henrotte, AxiHenrotteStiffnessBFI,
                          AxiHenrotteSigmaMassBFI)

A = 1.0


def _sphere_meridian_mesh(h):
    """Half-disk (r>=0) of radius A: the meridian section of a sphere, axis at r=0."""
    f = WorkPlane().Circle(0, 0, A).Face()
    half = MoveTo(0, -A).Rectangle(A, 2 * A).Face()
    g = f * half
    g.faces.name = "cond"
    g.edges.Min(X).name = "axis"
    g.edges.Max(X).name = "outer"
    return Mesh(OCCGeometry(g, dim=2).GenerateMesh(maxh=h))


def _assemble(fes):
    aK = BilinearForm(fes, symmetric=True); aK += AxiHenrotteStiffnessBFI(CF(1.0)); aK.Assemble()
    aM = BilinearForm(fes, symmetric=True); aM += AxiHenrotteSigmaMassBFI(CF(1.0)); aM.Assemble()
    n = fes.ndof
    free = np.array([i for i in range(n) if fes.FreeDofs()[i]], dtype=int)

    def dense(mat):
        rr, cc, vv = mat.COO()
        Amat = sp.csr_matrix((np.asarray(vv, float),
                              (np.asarray(rr, int), np.asarray(cc, int))),
                             shape=(n, n)).toarray()[np.ix_(free, free)]
        return 0.5 * (Amat + Amat.T)

    return dense(aK.mat), dense(aM.mat)


def test_p2_trig_center_dofs_unused():
    mesh = _sphere_meridian_mesh(0.20)
    fes = H1Henrotte(mesh, order=2, dirichlet="axis|outer")
    nv, ne = mesh.nv, mesh.nedge
    from ngsolve import COUPLING_TYPE
    n_centers = 0
    for i in range(mesh.ne):                      # one center slot per VOL element
        dof = nv + ne + i
        assert fes.CouplingType(dof) == COUPLING_TYPE.UNUSED_DOF, \
            f"P2 trig-center DOF {dof} should be UNUSED, got {fes.CouplingType(dof)}"
        assert not fes.FreeDofs()[dof], f"trig-center DOF {dof} must not be free"
        n_centers += 1
    assert n_centers == mesh.ne
    print(f"[OK] all {n_centers} P2 trig-center slots are UNUSED and not free")


def test_p2_eddy_eigenvalue_full_rank():
    # The core fix: P2 K/M are full rank, so the generalized eig solves (pre-fix
    # this raised LinAlgError -- M had a zero eigenvalue from the dead DOFs).
    # P1 and P2 must agree on the slowest mode (P2 is consistent, not spurious).
    lams = {}
    for order in (1, 2):
        K, M = _assemble(H1Henrotte(_sphere_meridian_mesh(0.10), order=order,
                                    dirichlet="axis|outer"))
        w = eigh(K, M, eigvals_only=True, subset_by_index=[0, 0])   # raises if M not PD
        lam = float(w[0])
        assert np.isfinite(lam) and lam > 0, f"order={order}: bad lam_min {lam}"
        lams[order] = lam
        print(f"  order={order}: lambda_min = {lam:.4f}")
    rel = abs(lams[2] - lams[1]) / lams[1]
    assert rel < 0.15, \
        f"P1 ({lams[1]:.3f}) and P2 ({lams[2]:.3f}) disagree by {rel:.1%} -- P2 may be spurious"
    print(f"[OK] P2 triangle is full-rank on axis; P1 vs P2 agree to {rel:.1%}")


if __name__ == "__main__":
    test_p2_trig_center_dofs_unused()
    test_p2_eddy_eigenvalue_full_rank()
    print("\nAll P2 axis-eddy regression tests passed.")
