"""Minimal 2-triangle P2 reproducer.

Builds the smallest possible 2D netgen mesh that exercises the C++ P2 TRIG
dispatch path: 2 triangles forming a single quad (general, not axis-aligned).
Walks through each step with explicit prints to localize where the segfault
happens (mesh ctor / NGSolve Mesh() wrap / GetMaterials / H1Henrotte ctor /
fes.Update / Assemble).
"""
from __future__ import annotations
import sys, traceback
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
)


def step(msg):
    print(f"[STEP] {msg}", flush=True)


def build_2tri_mesh():
    """2 triangles sharing one edge, off-axis, single material 'conductor'.
    Vertices:
        0: (1e-3, 0)
        1: (5e-3, 0)
        2: (5e-3, 2e-3)
        3: (1e-3, 2e-3)
    Triangles: [0,1,2] and [0,2,3].
    """
    m = NgMesh()
    m.dim = 2
    step("SetMaterial")
    m.SetMaterial(1, "conductor")
    step("Add FaceDescriptor(domin=0)")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    m.SetBCName(0, "outer")

    step("Add 4 points")
    pts = [
        m.Add(MeshPoint(Pnt(1e-3, 0,    0))),
        m.Add(MeshPoint(Pnt(5e-3, 0,    0))),
        m.Add(MeshPoint(Pnt(5e-3, 2e-3, 0))),
        m.Add(MeshPoint(Pnt(1e-3, 2e-3, 0))),
    ]
    step("Add 2 triangles")
    m.Add(Element2D(1, [pts[0], pts[1], pts[2]]))
    m.Add(Element2D(1, [pts[0], pts[2], pts[3]]))
    step("Add boundary edges (4)")
    m.Add(Element1D([pts[0], pts[1]], index=1))
    m.Add(Element1D([pts[1], pts[2]], index=1))
    m.Add(Element1D([pts[2], pts[3]], index=1))
    m.Add(Element1D([pts[3], pts[0]], index=1))
    return m


def main():
    ngsglobals.msg_level = 0
    step("Build netgen mesh")
    ngmesh = build_2tri_mesh()
    step(f"netgen ne (2D)?  {ngmesh}")

    step("Wrap in NGSolve Mesh()")
    mesh = Mesh(ngmesh)
    step(f"NGSolve mesh ne={mesh.ne}  nv={mesh.nv}")

    step("Call mesh.GetMaterials()")
    mats = mesh.GetMaterials()
    print(f"  Materials: {mats}", flush=True)

    step("Call mesh.GetBoundaries()")
    bnds = mesh.GetBoundaries()
    print(f"  Boundaries: {bnds}", flush=True)

    step("Import H1Henrotte")
    from radia.radia_axifemm import (
        H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
    )

    step("Create H1Henrotte(mesh, order=2)")
    fes = H1Henrotte(mesh, order=2, dirichlet="outer")
    step(f"ndof={fes.ndof}")

    step("Make BFI K")
    from math import pi
    MU0 = 4 * pi * 1e-7
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    step("Assemble K")
    with TaskManager(): a.Assemble()
    print(f"  K nnz reported by NGSolve: {a.mat.nze}", flush=True)

    step("Make BFI sigma-M")
    sigma_cf = mesh.MaterialCF({"conductor": 5.8e7}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    step("Assemble M")
    with TaskManager(): m_bf.Assemble()
    print(f"  M nnz reported by NGSolve: {m_bf.mat.nze}", flush=True)

    step("DONE — no segfault")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
