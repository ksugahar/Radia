"""Axisymmetric Cu disk eddy-current eigen reference (AXIFEM Q1 / Q2 quads).

Solid copper disk (radius R_DISK, thickness T_DISK, conductivity SIGMA_CU)
in an air box, discretised on a structured axis-aligned quad mesh so that
the closed-form AXIFEM quad elements are dispatched:

    * ``solve_disk``    -- H1Henrotte(order=1), Q1 quads
    * ``solve_disk_q2`` -- H1Henrotte(order=2), Q2 quads (9 DOFs / quad)

Both assemble the stiffness K (1/mu) and the sigma-mass M (conductor only)
and return the slowest eddy-current decay time constants
``tau_n = 1 / lambda_n`` of the generalized eigenproblem ``K u = lambda M u``
(scipy ``eigsh``, shift-invert at sigma=0), in microseconds, plus mesh size
information. ``BEM_TAU`` is the stored axisymmetric integral-equation (BEM)
modal spectrum of the same disk, used as the comparison reference.

The mesh builder reads the module globals ``R_DISK``/``T_DISK`` at call time
and the solvers read ``SIGMA_CU``/``MU0`` at call time, so callers may
temporarily override them on the module object for other disk variants.

Requires the native ``radia.axifem`` extension plus ngsolve and scipy.
"""

import os
from math import pi

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from netgen.meshing import (
    Mesh as NgMesh, EdgeDescriptor, Element1D, Element2D, FaceDescriptor,
    MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager,
    ngsglobals,
)
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
B0 = 1.0  # imposed axial flux density [T] (excitation A_phi = B0 r / 2)

# Stored axisymmetric BEM modal spectrum of the disk [us].
BEM_TAU = [224.3070587702379, 88.42395618426703, 49.53986759314118,
           32.090204181064834, 23.339598470067966, 22.588437092711622]
# Pure-Python Q1 prototype leading time constants [us].
PYTHON_Q1_TAU = [223.06, 87.60, 48.86, 31.49]

# Location of the refined BEM reference result (regenerable, not shipped).
BEM_TAU_REF_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "bem_disk",
    "bem_disk_axisym_v3_refined.json",
))


def make_structured_disk_quad_mesh(NR_disk=40, Nz_disk=8, NR_air=15, Nz_air=15,
                                    R_air=200e-3, Z_air=200e-3):
    """Structured axis-aligned quad mesh of the disk + air box, as NGSolve Mesh.

    Layout in (r, z):
        - Disk: r in [0, R_disk], z in [-T/2, T/2], NR_disk x Nz_disk uniform cells.
        - Air right of disk: r in [R_disk, R_air], geometric grading.
        - Air above disk: z in [T/2, Z_air], geometric grading.
        - Air below disk: z in [-Z_air, -T/2], geometric grading.

    Boundary names: "axis" (r=0), "right" (r=R_air), "top" (z=Z_air), "bot" (z=-Z_air).
    Material names: "conductor" (in disk), "air" (elsewhere).
    """
    r_disk = np.linspace(0, R_DISK, NR_disk + 1)
    r_air  = np.geomspace(R_DISK, R_air, NR_air + 1)[1:]
    r_grid = np.concatenate([r_disk, r_air])

    z_disk = np.linspace(-T_DISK / 2, T_DISK / 2, Nz_disk + 1)
    z_above = np.geomspace(T_DISK/2 + 1e-7, Z_air, Nz_air + 1)[1:]
    z_below = -np.geomspace(T_DISK/2 + 1e-7, Z_air, Nz_air + 1)[1:][::-1]
    z_grid = np.concatenate([z_below, z_disk, z_above])

    NR = len(r_grid) - 1
    NZ = len(z_grid) - 1

    ngmesh = NgMesh()
    ngmesh.dim = 2

    # Materials: index 1 = air, 2 = conductor.
    ngmesh.SetMaterial(1, "air")
    ngmesh.SetMaterial(2, "conductor")

    # FaceDescriptors for the 4 boundary edges.
    fd_axis  = ngmesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    fd_right = ngmesh.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))
    fd_top   = ngmesh.Add(FaceDescriptor(surfnr=3, domin=0, bc=3))
    fd_bot   = ngmesh.Add(FaceDescriptor(surfnr=4, domin=0, bc=4))
    ngmesh.SetBCName(0, "axis")
    ngmesh.SetBCName(1, "right")
    ngmesh.SetBCName(2, "top")
    ngmesh.SetBCName(3, "bot")
    for boundary, name in enumerate(("axis", "right", "top", "bot"), start=1):
        edge = EdgeDescriptor()
        edge.edgenr = boundary
        edge.surfnr = (boundary, -1)
        edge.domin = 1
        edge.domout = 0
        edge.name = name
        ngmesh.Add(edge)

    # Add nodes.
    pids = np.empty((NZ + 1, NR + 1), dtype=object)
    for j in range(NZ + 1):
        for i in range(NR + 1):
            pids[j, i] = ngmesh.Add(MeshPoint(Pnt(r_grid[i], z_grid[j], 0)))

    # Add quads with material index. Vertex order MUST be CCW in (r, z):
    # (i,j) -> (i+1,j) -> (i+1,j+1) -> (i,j+1).
    for j in range(NZ):
        for i in range(NR):
            r_c = 0.5 * (r_grid[i] + r_grid[i + 1])
            z_c = 0.5 * (z_grid[j] + z_grid[j + 1])
            in_disk = (r_c < R_DISK) and (abs(z_c) < T_DISK / 2)
            mat_index = 2 if in_disk else 1
            ngmesh.Add(Element2D(mat_index, [
                pids[j,     i    ],
                pids[j,     i + 1],
                pids[j + 1, i + 1],
                pids[j + 1, i    ],
            ]))

    # Add boundary edges. Each Element1D needs a face index pointing at a
    # FaceDescriptor (whose bc index gives the BC name).
    # Axis: r=0 -> column i=0
    for j in range(NZ):
        ngmesh.Add(Element1D([pids[j, 0], pids[j + 1, 0]], index=1))
    # Right: r=R_air -> column i=NR
    for j in range(NZ):
        ngmesh.Add(Element1D([pids[j, NR], pids[j + 1, NR]], index=2))
    # Top: z=Z_air -> row j=NZ
    for i in range(NR):
        ngmesh.Add(Element1D([pids[NZ, i], pids[NZ, i + 1]], index=3))
    # Bot: z=-Z_air -> row j=0
    for i in range(NR):
        ngmesh.Add(Element1D([pids[0, i], pids[0, i + 1]], index=4))

    return Mesh(ngmesh)


def to_scipy_csr(mat, n):
    """NGSolve sparse matrix -> symmetrised scipy CSR (n x n)."""
    rs, cs, vs = mat.COO()
    K = sp.csr_matrix(
        (np.asarray(vs, dtype=np.float64),
         (np.asarray(rs, dtype=np.int64), np.asarray(cs, dtype=np.int64))),
        shape=(n, n))
    return (K + K.T) * 0.5


def _eigen_time_constants(fes, mesh, N_stages, fmt):
    """Assemble K, M on ``fes`` and return eigsh time constants [us] + free count."""
    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m.Assemble()

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]

    eigs, _ = spla.eigsh(K_red, k=min(N_stages, len(free) // 2),
                         M=M_red, sigma=0.0, which="LM",
                         tol=1e-10, maxiter=3000)
    eigsh_taus = sorted((1.0/e) * 1e6 for e in eigs)[::-1]
    print(f"  scipy eigsh tau_n[:{len(eigsh_taus)}] = "
          f"{[format(t, fmt) for t in eigsh_taus]} us")
    return eigsh_taus, int(len(free))


def solve_disk(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air, N_stages=6,
               label=""):
    """Q1 quad disk: slowest ``N_stages`` eddy time constants [us] + mesh info."""
    print(f"\n=== Q1 quad: NR_d={NR_disk} Nz_d={Nz_disk} NR_a={NR_air} "
          f"Nz_a={Nz_air} R_air={R_air*1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                            R_air, Z_air)
    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne} ndof={fes.ndof}  free={n_free}  "
          f"materials={mesh.GetMaterials()}")
    eigsh_taus, free = _eigen_time_constants(fes, mesh, N_stages, ".3f")
    return {"eigsh_tau_us": eigsh_taus,
            "mesh": {"ne": mesh.ne, "ndof": fes.ndof, "free": free}}


def solve_disk_q2(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air, N_stages=6,
                  label=""):
    """Q2 quad disk: slowest ``N_stages`` eddy time constants [us] + mesh info."""
    print(f"\n=== Q2 quad: NR_d={NR_disk} Nz_d={Nz_disk} NR_a={NR_air} "
          f"Nz_a={Nz_air} R_air={R_air*1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                          R_air, Z_air)
    fes = H1Henrotte(mesh, order=2, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne} ndof={fes.ndof}  free={n_free}  "
          f"materials={mesh.GetMaterials()}")
    eigsh_taus, free = _eigen_time_constants(fes, mesh, N_stages, ".4f")
    return {"eigsh_tau_us": eigsh_taus,
            "mesh": {"ne": mesh.ne, "ndof": fes.ndof, "free": free}}


def main():
    cases_q1 = [
        (20,  4,  10, 10, 100e-3, 100e-3, "coarse"),
        (40,  8,  15, 15, 200e-3, 200e-3, "medium"),
        (80, 16,  20, 20, 500e-3, 500e-3, "fine"),
        (160, 32, 25, 25, 500e-3, 500e-3, "very fine"),
    ]
    cases_q2 = [
        (10, 4,  8,  8, 100e-3, 100e-3, "coarse"),
        (20, 8, 12, 12, 200e-3, 200e-3, "medium"),
        (40, 16, 15, 15, 500e-3, 500e-3, "fine"),
    ]
    rows = []
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases_q1:
        rows.append(("Q1 " + lbl, solve_disk(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
                                              N_stages=6, label=lbl)))
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases_q2:
        rows.append(("Q2 " + lbl, solve_disk_q2(NR_d, Nz_d, NR_a, Nz_a, R_a,
                                                 Z_a, N_stages=6, label=lbl)))

    print("\n" + "="*70)
    print("Summary: first eddy time constant vs stored BEM modal reference")
    print("="*70)
    print(f"{'case':<14} {'ne':>8} {'free':>6}  {'tau_1 us':>10}  {'/BEM':>7}")
    for lbl, r in rows:
        t1 = r["eigsh_tau_us"][0]
        print(f"{lbl:<14} {r['mesh']['ne']:>8} {r['mesh']['free']:>6}  "
              f"{t1:>10.4f}  {t1 / BEM_TAU[0]:>7.4f}")


if __name__ == "__main__":
    main()
