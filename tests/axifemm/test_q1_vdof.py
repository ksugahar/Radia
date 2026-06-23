# -*- coding: utf-8 -*-
r"""Regression test: the Q1 axis-aligned quad stiffness is a TRUE V-DOF operator.

solve_axi_magnetostatic dispatches order=1 to the V-DOF custom-BFI path, where the
quad stiffness must be  K_ij = 2pi/mu r_i r_j INT grad(psi_i).grad(psi_j)/r
(Q1StiffnessVDof), NOT the old A=psi energy form (Q1StiffnessCoef, which cannot
represent a uniform axial B_z at order 1).  A correct V-DOF stiffness has the
uniform field A_phi = B0 r/2 in its kernel, so imposing it on the boundary of a
structured axis-aligned quad mesh and solving K V = 0 in the interior must
reproduce B0 r/2 to solver precision.

AxiHenrotteFE_Q1_AxisAligned supports only axis-aligned rectangles, so this uses a
structured quad mesh (general/curved quads at order 1 are not supported).

Requires the rebuilt axifem extension (Build.ps1 -AxiFemmOnly); loads src/radia.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from ngsolve import BilinearForm, GridFunction, CoefficientFunction
from ngsolve.meshes import MakeStructured2DMesh
from radia.axifem import H1Henrotte, AxiHenrotteStiffnessBFI

MU0, B0 = 4e-7 * math.pi, 1.3
R, Z0, Z1 = 1.0, -0.5, 0.5


def _uniform_field_residual(nx, ny):
    mesh = MakeStructured2DMesh(quads=True, nx=nx, ny=ny,
                                mapping=lambda x, y: (R * x, Z0 + (Z1 - Z0) * y))
    fes = H1Henrotte(mesh, order=1)
    aK = BilinearForm(fes, symmetric=True)
    aK += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    aK.Assemble()
    n = fes.ndof
    r = np.array([mesh[v].point[0] for v in mesh.vertices])
    z = np.array([mesh[v].point[1] for v in mesh.vertices])
    Vex = B0 * r / 2.0                       # exact uniform-field potential A_phi

    # dense K via mat-vec (avoids COO symmetric-storage ambiguity)
    e = GridFunction(fes); res = GridFunction(fes)
    K = np.zeros((n, n))
    for i in range(n):
        e.vec.FV().NumPy()[:] = 0.0; e.vec.FV().NumPy()[i] = 1.0
        res.vec.data = aK.mat * e.vec
        K[:, i] = res.vec.FV().NumPy()

    onb = (r < 1e-9) | (np.abs(r - R) < 1e-9) | (np.abs(z - Z0) < 1e-9) | (np.abs(z - Z1) < 1e-9)
    free = np.where(~onb)[0]; bnd = np.where(onb)[0]
    Vf = spla.spsolve(sp.csr_matrix(K[np.ix_(free, free)]), -K[np.ix_(free, bnd)] @ Vex[bnd])
    return np.linalg.norm(Vf - Vex[free]) / np.linalg.norm(Vex[free])


def test_q1_vdof_reproduces_uniform_field():
    for (nx, ny) in [(4, 4), (8, 8), (12, 12)]:
        rel = _uniform_field_residual(nx, ny)
        print(f"  Q1 {nx}x{ny}: ||V_int - B0 r/2||/||.|| = {rel:.2e}")
        assert rel < 1e-8, \
            f"Q1 {nx}x{ny} does not reproduce a uniform B_z (rel {rel:.2e}) -- " \
            f"the stiffness is not a V-DOF operator"
    print("[OK] Q1 axis-aligned quad stiffness is a true V-DOF operator")


if __name__ == "__main__":
    test_q1_vdof_reproduces_uniform_field()
    print("\nQ1 V-DOF test passed.")
