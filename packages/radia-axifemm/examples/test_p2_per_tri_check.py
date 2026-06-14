"""Print element K for tri_a (1,0)-(5,0)-(5,2) and compare directly between
Python (AxifemmP2Triangle) and what C++ NGSolve should produce.

If element matrices DIFFER, the C++ integrator has a bug that wasn't caught
by test_p2_assembly_diff.py (which used a different triangle).
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.sparse as sp
from math import pi

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from axifemm_p2_triangle import AxifemmP2Triangle

MU0 = 4 * pi * 1e-7
SIGMA = 5.8e7


def py_K(V):
    rn = np.array([V[0][0], V[1][0], V[2][0],
                   0.5*(V[0][0]+V[1][0]), 0.5*(V[1][0]+V[2][0]), 0.5*(V[2][0]+V[0][0])])
    zn = np.array([V[0][1], V[1][1], V[2][1],
                   0.5*(V[0][1]+V[1][1]), 0.5*(V[1][1]+V[2][1]), 0.5*(V[2][1]+V[0][1])])
    return AxifemmP2Triangle(rn, zn).stiffness(mu_r=MU0)


def cpp_K_single(V):
    """Build single-triangle mesh and extract the 6x6 nonzero block."""
    ngsglobals.msg_level = 0
    m = NgMesh()
    m.dim = 2
    m.SetMaterial(1, "conductor")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    m.SetBCName(0, "outer")
    pids = [m.Add(MeshPoint(Pnt(v[0], v[1], 0))) for v in V]
    m.Add(Element2D(1, pids))
    m.Add(Element1D([pids[0], pids[1]], index=1))
    m.Add(Element1D([pids[1], pids[2]], index=1))
    m.Add(Element1D([pids[2], pids[0]], index=1))
    mesh = Mesh(m)
    fes = H1Henrotte(mesh, order=2, dirichlet=None)
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()
    rs_, cs_, vs_ = a.mat.COO()
    K = sp.coo_matrix((np.array(vs_), (np.array(rs_), np.array(cs_))),
                     shape=(fes.ndof, fes.ndof)).toarray()
    # Order DOFs by coord matching local indices [v0, v1, v2, m01, m12, m20].
    target_coords = list(V) + [
        ((V[0][0] + V[1][0]) / 2, (V[0][1] + V[1][1]) / 2),
        ((V[1][0] + V[2][0]) / 2, (V[1][1] + V[2][1]) / 2),
        ((V[2][0] + V[0][0]) / 2, (V[2][1] + V[0][1]) / 2),
    ]
    dof_coords = [(v.point[0], v.point[1]) for v in mesh.vertices]
    for e in mesh.edges:
        ev = list(e.vertices)
        v0 = mesh.vertices[ev[0].nr]; v1 = mesh.vertices[ev[1].nr]
        dof_coords.append(((v0.point[0]+v1.point[0])/2, (v0.point[1]+v1.point[1])/2))
    perm = []
    for tc in target_coords:
        diffs = [(tc[0]-d[0])**2 + (tc[1]-d[1])**2 for d in dof_coords]
        perm.append(int(np.argmin(diffs)))
    return K[np.ix_(perm, perm)]


def main():
    tri_a = [(1e-3, 0), (5e-3, 0), (5e-3, 2e-3)]
    print(f"Tri_a vertices: {tri_a}")
    Kp = py_K(tri_a)
    Kc = cpp_K_single(tri_a)
    print(f"\nPython K_a (6x6):")
    print(Kp)
    print(f"\nC++ K_a (6x6, permuted to match local order):")
    print(Kc)
    print(f"\nDiff max abs: {np.max(np.abs(Kp - Kc)):.3e}")
    print(f"Diff rel: {np.max(np.abs(Kp - Kc))/np.max(np.abs(Kp)):.3e}")
    print(f"\nKp[4, 2] = {Kp[4, 2]:.4e}")
    print(f"Kc[4, 2] = {Kc[4, 2]:.4e}")


if __name__ == "__main__":
    main()
