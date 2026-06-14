"""Direct K, M comparison: Python reference vs NGSolve+C++ assembly.

Build the same 3-triangle off-axis mesh, assemble K and M via:
  (a) Python reference: per-triangle AxifemmP2Triangle + manual global assembly.
  (b) NGSolve: H1Henrotte + AxiHenrotteStiffnessBFI/SigmaMassBFI.
Permute DOFs to match (NGSolve numbers vertex DOFs first, then edges by edge
index; Python references nodes_p2 ordering which is "vertex first, then mid-
edges in dictionary insertion order"). Compare K, M entry by entry.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.sparse as sp
import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
)
from math import pi

from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from axifemm_p2_triangle import AxifemmP2Triangle, MU0

SIGMA = 5.8e7


def make_mesh():
    """3-vertex off-axis 1-triangle test mesh, conductor everywhere."""
    m = NgMesh()
    m.dim = 2
    m.SetMaterial(1, "conductor")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    m.SetBCName(0, "outer")
    p0 = m.Add(MeshPoint(Pnt(1e-3, 0.0,  0)))
    p1 = m.Add(MeshPoint(Pnt(5e-3, 0.0,  0)))
    p2 = m.Add(MeshPoint(Pnt(3e-3, 2e-3, 0)))
    m.Add(Element2D(1, [p0, p1, p2]))
    m.Add(Element1D([p0, p1], index=1))
    m.Add(Element1D([p1, p2], index=1))
    m.Add(Element1D([p2, p0], index=1))
    coords = [(1e-3, 0.0), (5e-3, 0.0), (3e-3, 2e-3)]
    return Mesh(m), coords


def cpp_assembly(mesh):
    ngsglobals.msg_level = 0
    fes = H1Henrotte(mesh, order=2, dirichlet=None)
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m_bf.Assemble()
    def to_dense(bfmat, n):
        rs, cs, vs = bfmat.COO()
        return sp.coo_matrix((np.array(vs), (np.array(rs), np.array(cs))),
                             shape=(n, n)).toarray()
    return to_dense(a.mat, fes.ndof), to_dense(m_bf.mat, fes.ndof), fes


def py_assembly(coords):
    rn = np.array([
        coords[0][0], coords[1][0], coords[2][0],
        0.5*(coords[0][0]+coords[1][0]),
        0.5*(coords[1][0]+coords[2][0]),
        0.5*(coords[2][0]+coords[0][0]),
    ])
    zn = np.array([
        coords[0][1], coords[1][1], coords[2][1],
        0.5*(coords[0][1]+coords[1][1]),
        0.5*(coords[1][1]+coords[2][1]),
        0.5*(coords[2][1]+coords[0][1]),
    ])
    tri = AxifemmP2Triangle(rn, zn)
    Ke = tri.stiffness(mu_r=MU0)
    Me = tri.sigma_mass(sigma=SIGMA)
    return Ke, Me


def main():
    mesh, coords = make_mesh()
    print(f"mesh.ne={mesh.ne}  mesh.nv={mesh.nv}")
    Kc, Mc, fes = cpp_assembly(mesh)
    print(f"C++ ndof={fes.ndof}")
    Kp, Mp = py_assembly(coords)
    print(f"Python local element: K shape {Kp.shape}")

    # Get C++ ordering: vertex DOFs 0..nv-1, then edge DOFs nv..nv+ne-1.
    # For a single triangle, C++ assembly gives a 6x6 nonzero block + maybe
    # zero rows/cols.
    print(f"\nC++ K (nonzero block):")
    nz_rows = np.where(np.any(np.abs(Kc) > 1e-15, axis=1))[0]
    print(f"  nonzero rows: {nz_rows}")
    Kc_nz = Kc[np.ix_(nz_rows, nz_rows)]
    Mc_nz = Mc[np.ix_(nz_rows, nz_rows)]
    print(f"  C++ K block shape {Kc_nz.shape}")
    print(f"\nPython K (6x6):")
    print(Kp)
    print(f"\nC++ K (6x6 nz block):")
    print(Kc_nz)
    print(f"\nDiff K (max abs): {np.max(np.abs(Kp - Kc_nz)):.3e}")
    print(f"Diff K (relative): "
          f"{np.max(np.abs(Kp - Kc_nz)) / np.max(np.abs(Kp)):.3e}")

    print(f"\nPython M (6x6):")
    print(Mp)
    print(f"\nC++ M (6x6 nz block):")
    print(Mc_nz)
    print(f"\nDiff M (max abs): {np.max(np.abs(Mp - Mc_nz)):.3e}")

    # Permutation check: maybe the DOF ordering differs. Test all 6! permutations
    # is overkill; try just swapping vertex order vs edge order.
    print(f"\n=== Test permutations to find matching ordering ===")
    from itertools import permutations
    best_err = np.inf
    best_perm = None
    for perm in permutations(range(6)):
        P = np.array(perm)
        Kp_perm = Kp[np.ix_(P, P)]
        err = np.max(np.abs(Kp_perm - Kc_nz))
        if err < best_err:
            best_err = err
            best_perm = perm
    print(f"Best permutation: {best_perm}, max err = {best_err:.3e}")


if __name__ == "__main__":
    main()
