"""Compare C++ vs Python element matrix on the WORST-conditioned sphere
triangle (cond ~ 1e14 near polar cap). If they diverge here, that confirms
C++ CalcInverse<Mat<6,6>> instability is the root cause.
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

# Worst tri from test_p2_vandermonde_cond.py output:
WORST = {
    "rn": [0.0006470476127563026, 0.0, 0.0, 0.0003235238063781513, 0.0, 0.0003235238063781513],
    "zn": [-0.0024148145657226705, -0.005, -0.0025, -0.0037074072828613353, -0.00375, -0.002457407282861335],
}


def main():
    ngsglobals.msg_level = 0
    rn = np.array(WORST["rn"])
    zn = np.array(WORST["zn"])
    print(f"Vertices (3): {list(zip(rn[:3], zn[:3]))}")
    print(f"Mid-edges (3): {list(zip(rn[3:], zn[3:]))}")
    V = np.zeros((6, 6))
    for i in range(6):
        s = rn[i] * rn[i]
        V[i] = [1.0, s, zn[i], s*s, s*zn[i], zn[i]**2]
    print(f"Vandermonde cond: {np.linalg.cond(V):.3e}")

    # Python reference
    tri = AxifemmP2Triangle(rn, zn)
    Kp = tri.stiffness(mu_r=MU0)
    Mp = tri.sigma_mass(sigma=SIGMA)
    print(f"\nPython K (6x6): max|K|={np.max(np.abs(Kp)):.3e}")
    print(Kp)

    # C++ via 1-triangle netgen mesh
    m = NgMesh()
    m.dim = 2
    m.SetMaterial(1, "conductor")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    m.SetBCName(0, "outer")
    p0 = m.Add(MeshPoint(Pnt(rn[0], zn[0], 0)))
    p1 = m.Add(MeshPoint(Pnt(rn[1], zn[1], 0)))
    p2 = m.Add(MeshPoint(Pnt(rn[2], zn[2], 0)))
    m.Add(Element2D(1, [p0, p1, p2]))
    m.Add(Element1D([p0, p1], index=1))
    m.Add(Element1D([p1, p2], index=1))
    m.Add(Element1D([p2, p0], index=1))
    mesh = Mesh(m)
    fes = H1Henrotte(mesh, order=2, dirichlet=None)
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()
    rs_, cs_, vs_ = a.mat.COO()
    Kc = sp.coo_matrix((np.array(vs_), (np.array(rs_), np.array(cs_))),
                       shape=(fes.ndof, fes.ndof)).toarray()
    # Pick 6 nonzero rows.
    nz = np.where(np.any(np.abs(Kc) > 1e-25, axis=1))[0]
    Kc_block = Kc[np.ix_(nz, nz)]
    print(f"\nC++ K (nonzero block {Kc_block.shape}): max|K|={np.max(np.abs(Kc_block)):.3e}")
    print(Kc_block)

    print(f"\nDiff max abs: {np.max(np.abs(Kp - Kc_block)):.3e}")
    print(f"Diff rel: {np.max(np.abs(Kp - Kc_block)) / max(np.max(np.abs(Kp)), 1e-30):.3e}")


if __name__ == "__main__":
    main()
