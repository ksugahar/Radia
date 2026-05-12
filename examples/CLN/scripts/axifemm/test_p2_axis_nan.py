"""Diagnostic: trigger NaN in C++ P2 Vandermonde for axis-touching tri.

Set up a single triangle similar to those near the sphere's pole (one vertex
on axis r=0, two slightly off-axis), instantiate AxiHenrotteFE_P2_Triangle
directly via the C++ binding, and check whether the K matrix it builds
(through an assembly on a 1-element mesh) contains NaN.

Also do the same via Python reference AxifemmP2Triangle for cross-check.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
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
import scipy.sparse as sp

MU0 = 4 * pi * 1e-7


def make_single_tri_mesh(r0, z0, r1, z1, r2, z2):
    """1-triangle 2D netgen mesh, conductor material, all-outer BC."""
    m = NgMesh()
    m.dim = 2
    m.SetMaterial(1, "conductor")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    m.SetBCName(0, "outer")
    p0 = m.Add(MeshPoint(Pnt(r0, z0, 0)))
    p1 = m.Add(MeshPoint(Pnt(r1, z1, 0)))
    p2 = m.Add(MeshPoint(Pnt(r2, z2, 0)))
    m.Add(Element2D(1, [p0, p1, p2]))
    m.Add(Element1D([p0, p1], index=1))
    m.Add(Element1D([p1, p2], index=1))
    m.Add(Element1D([p2, p0], index=1))
    return Mesh(m)


def cpp_K_M(r0, z0, r1, z1, r2, z2):
    ngsglobals.msg_level = 0
    mesh = make_single_tri_mesh(r0, z0, r1, z1, r2, z2)
    fes = H1Henrotte(mesh, order=2, dirichlet=None)
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()
    rows, cols, vals = a.mat.COO()
    K = sp.coo_matrix(
        (np.array(vals), (np.array(rows), np.array(cols))),
        shape=(fes.ndof, fes.ndof),
    ).toarray()

    sigma_cf = mesh.MaterialCF({"conductor": 5.8e7}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m_bf.Assemble()
    rows, cols, vals = m_bf.mat.COO()
    M = sp.coo_matrix(
        (np.array(vals), (np.array(rows), np.array(cols))),
        shape=(fes.ndof, fes.ndof),
    ).toarray()
    return K, M


def py_K_M(r0, z0, r1, z1, r2, z2):
    """Python reference K, M for the same triangle (with midpoints)."""
    from axifemm_p2_triangle import AxifemmP2Triangle
    rn = np.array([r0, r1, r2, 0.5*(r0+r1), 0.5*(r1+r2), 0.5*(r2+r0)])
    zn = np.array([z0, z1, z2, 0.5*(z0+z1), 0.5*(z1+z2), 0.5*(z2+z0)])
    tri = AxifemmP2Triangle(rn, zn)
    return tri.stiffness(mu_r=MU0), tri.sigma_mass(sigma=5.8e7)


def report(label, K, M):
    finite_K = np.isfinite(K).all()
    finite_M = np.isfinite(M).all()
    print(f"  {label}:")
    print(f"    K finite={finite_K}  max|K|={np.nanmax(np.abs(K)):.3e}")
    print(f"    M finite={finite_M}  max|M|={np.nanmax(np.abs(M)):.3e}")


def main():
    print("=" * 64)
    print("P2 axis-touching triangle: Python vs C++ comparison")
    print("=" * 64)

    # Sphere mesh polar-region triangles: TWO vertices on axis (r=0 at theta=0).
    # These come from the (j=0, i=0..) cells in make_sphere_axisym_trig_mesh.
    cases = [
        ("fully off-axis",       (1e-3, 0,    5e-3, 0,    3e-3, 2e-3)),
        ("one axis vertex (v0)", (0.0,  0.0,  5e-3, 0,    3e-3, 2e-3)),
        ("TWO axis vertices (polar cap of sphere)",
         (0.0, 0.0,    0.0, 1.25e-3,    1.25e-3 * np.sin(pi/12), 1.25e-3 * np.cos(pi/12))),
        ("TWO axis vertices, dr larger",
         (0.0, 0.0,    0.0, 1.25e-3,    1.25e-3 * np.sin(pi/6), 1.25e-3 * np.cos(pi/6))),
    ]

    for name, coords in cases:
        print(f"\n[{name}]")
        print(f"  coords = {coords}")
        try:
            Kp, Mp = py_K_M(*coords)
            report("python", Kp, Mp)
        except Exception as e:
            print(f"  python: raised {type(e).__name__}: {e}")
        try:
            Kc, Mc = cpp_K_M(*coords)
            report("cpp   ", Kc, Mc)
        except Exception as e:
            print(f"  cpp   : raised {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
