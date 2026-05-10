"""Phase 2-D numerical correctness test:
   Solve the Cu disk axisymmetric eddy-current eigenvalue problem on
   H1Henrotte and check that the leading time constant matches the Python
   reference implementation (axifemm/axifemm_quad.py Q1 = 223.06 us).

   Formulation (option (a), BilinearForm-side 1/r weighting):
     unknown u = psi = 2*pi*r*A_phi      (FEMM stream function)
     stiffness:  K_ij = (1/(2*pi*mu)) * int( grad(shape_i).grad(shape_j) / r dx )
     sigma mass: M_ij = (sigma/(2*pi))  * int( shape_i shape_j / r dx )
     eigenvalue: K v = lambda M v,  tau = 1/lambda
     Dirichlet:  u = 0 on axis (r=0) and on far-field box.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ngsolve import (
    Mesh, BilinearForm, InnerProduct, grad, dx, x, CoefficientFunction,
)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia.radia_axifemm import H1Henrotte


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
MU0 = 4 * np.pi * 1e-7

# Python axifemm Q1 reference (memory: project_axifemm_q1_validation)
PYTHON_Q1_REF = [223.06, 87.60, 48.86, 31.49]  # us


def build_disk_mesh(R_air, Z_air, maxh_disk, maxh_air):
    air = MoveTo(0, -Z_air).Rectangle(R_air, 2 * Z_air).Face()
    disk = MoveTo(0, -T_DISK / 2).Rectangle(R_DISK, T_DISK).Face()
    disk.faces.name = "conductor"
    air.faces.name = "air"
    disk.maxh = maxh_disk
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "right"
    air.edges.Min(Y).name = "bot"
    air.edges.Max(Y).name = "top"
    shape = Glue([air, disk])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_air))


def to_scipy(mat):
    rows, cols, vals = mat.COO()
    rows = np.asarray(rows, dtype=np.int64)
    cols = np.asarray(cols, dtype=np.int64)
    vals = np.asarray(vals, dtype=np.float64)
    n = max(int(rows.max()), int(cols.max())) + 1
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def solve_disk(R_air, Z_air, maxh_disk, maxh_air, n_eigs=4):
    mesh = build_disk_mesh(R_air, Z_air, maxh_disk, maxh_air)

    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
    print(f"  mesh: nv={mesh.nv} ne={mesh.ne}   ndof={fes.ndof}   "
          f"freedofs={sum(1 for f in fes.FreeDofs() if f)}")
    u, v = fes.TnT()

    # CoefficientFunction: sigma in conductor, 0 elsewhere.
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    # x is the radial coordinate r (NGSolve 2D x-axis).
    a = BilinearForm(fes, symmetric=True)
    a += (1.0 / (2 * np.pi * MU0 * x)) * InnerProduct(grad(u), grad(v)) * dx
    a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += (sigma_cf / (2 * np.pi * x)) * u * v * dx
    m.Assemble()

    # Reduce to free DOFs and convert to scipy.
    K = to_scipy(a.mat)
    M = to_scipy(m.mat)

    free_mask = np.array([fes.FreeDofs()[i] for i in range(fes.ndof)], dtype=bool)
    free = np.where(free_mask)[0]
    K = K[free[:, None], free[None, :]]
    M = M[free[:, None], free[None, :]]
    K = (K + K.T) / 2.0
    M = (M + M.T) / 2.0

    eigs, _ = spla.eigsh(K, k=n_eigs, M=M, sigma=0.0, which="LM",
                         tol=1e-10, maxiter=3000)
    taus_us = np.sort(1.0 / eigs)[::-1] * 1e6
    return taus_us, mesh.nv, mesh.ne


def main():
    print("=== Phase 2-D: H1Henrotte axisymmetric Cu disk eigenvalue ===")
    print(f"    Python axifemm Q1 reference tau_1..tau_4 = {PYTHON_Q1_REF} us")

    cases = [
        (50e-3,  50e-3,  0.5e-3, 5e-3,  "coarse"),
        (100e-3, 100e-3, 0.3e-3, 5e-3,  "medium"),
        (200e-3, 200e-3, 0.2e-3, 5e-3,  "fine"),
    ]

    print(f"\n{'case':<8} {'tau_1 us':>10} {'tau_2':>8} {'tau_3':>8} {'tau_4':>8} "
          f"{'r1/Q1':>8}")
    for R_air, Z_air, h_d, h_a, label in cases:
        try:
            taus, _, _ = solve_disk(R_air, Z_air, h_d, h_a)
            r1 = taus[0] / PYTHON_Q1_REF[0]
            print(f"{label:<8} {taus[0]:>10.3f} {taus[1]:>8.2f} {taus[2]:>8.2f} "
                  f"{taus[3]:>8.2f} {r1:>8.4f}")
        except Exception as e:
            print(f"{label:<8} ERROR: {e}")
            import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()
