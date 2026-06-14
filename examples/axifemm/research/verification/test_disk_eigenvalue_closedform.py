"""Phase 2-E disk validation:
   Cu disk eddy-current eigenvalue problem solved with the new Custom
   Integrators (closed-form Mathematica element matrices, no Gauss
   quadrature). Target: tau_1 within 1% of Python axifemm Q1 = 223.06 us.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager,
)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
MU0 = 4 * np.pi * 1e-7

# Reference: axifemm/disk_hiruma_quad.py Q1 result (memory: project_axifemm_q1_validation)
PYTHON_Q1_TAU = [223.06, 87.60, 48.86, 31.49]  # us


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


def to_scipy(mat, n):
    rows, cols, vals = mat.COO()
    return sp.csr_matrix((np.asarray(vals, dtype=np.float64),
                           (np.asarray(rows, dtype=np.int64),
                            np.asarray(cols, dtype=np.int64))),
                          shape=(n, n))


def solve_disk(R_air, Z_air, maxh_disk, maxh_air, n_eigs=4):
    mesh = build_disk_mesh(R_air, Z_air, maxh_disk, maxh_air)

    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne}  ndof={fes.ndof}  free={n_free}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m.Assemble()

    K = to_scipy(a.mat, fes.ndof)
    M = to_scipy(m.mat, fes.ndof)

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]])
    K = K[free[:, None], free[None, :]]
    M = M[free[:, None], free[None, :]]
    K = (K + K.T) / 2.0
    M = (M + M.T) / 2.0

    eigs, _ = spla.eigsh(K, k=n_eigs, M=M, sigma=0.0, which="LM",
                         tol=1e-10, maxiter=3000)
    taus = np.sort(1.0 / eigs)[::-1] * 1e6
    return taus, mesh.ne


def main():
    print("=== Phase 2-E: H1Henrotte + closed-form Integrators on Cu disk ===")
    print(f"    Python axifemm Q1 reference tau_1..tau_4 = {PYTHON_Q1_TAU} us\n")

    cases = [
        (50e-3,  50e-3,  0.5e-3, 5e-3,  "coarse"),
        (100e-3, 100e-3, 0.3e-3, 5e-3,  "medium"),
        (200e-3, 200e-3, 0.2e-3, 5e-3,  "fine"),
        (500e-3, 500e-3, 0.15e-3, 10e-3, "very fine"),
    ]

    print(f"{'case':<10} {'ne':>6} {'tau_1':>10} {'tau_2':>8} {'tau_3':>8} "
          f"{'tau_4':>8} {'r1/Q1':>8}")
    for R_air, Z_air, h_d, h_a, label in cases:
        try:
            taus, ne = solve_disk(R_air, Z_air, h_d, h_a)
            r1 = taus[0] / PYTHON_Q1_TAU[0]
            print(f"{label:<10} {ne:>6} {taus[0]:>10.3f} {taus[1]:>8.2f} "
                  f"{taus[2]:>8.2f} {taus[3]:>8.2f} {r1:>8.4f}")
        except Exception as e:
            print(f"{label:<10} ERROR: {e}")
            import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()
