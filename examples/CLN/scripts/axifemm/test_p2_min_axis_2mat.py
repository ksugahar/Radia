"""Minimal axis-touching + 2-material reproducer.

Sphere-like setup: a few triangles, with some vertices at r=0 (axis), and
two materials air/conductor. This is the next-smaller case after
test_p2_min_2tri (which passed) toward the failing test_p2_cpp_sphere.
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


def build_axis_2mat_mesh():
    """4 triangles:
        - 2 with axis-touching vertices (r=0)
        - 2 fully off-axis
    Materials: index 1=air, 2=conductor.
    """
    m = NgMesh()
    m.dim = 2
    step("SetMaterial 1=air, 2=conductor")
    m.SetMaterial(1, "air")
    m.SetMaterial(2, "conductor")
    step("Add 2 FaceDescriptors (working test_hiruma_disk_q1 pattern: domin=0)")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))  # idx 1 -> Element2D(1, ...)
    m.Add(FaceDescriptor(surfnr=2, domin=0, bc=1))  # idx 2 -> Element2D(2, ...)
    m.SetBCName(0, "outer")

    step("Add 6 points (some on axis)")
    # Axis vertices (r=0)
    p0 = m.Add(MeshPoint(Pnt(0.0,    -1e-3, 0)))  # axis bottom
    p1 = m.Add(MeshPoint(Pnt(0.0,     1e-3, 0)))  # axis middle
    p2 = m.Add(MeshPoint(Pnt(0.0,     3e-3, 0)))  # axis top
    # Off-axis vertices
    p3 = m.Add(MeshPoint(Pnt(2e-3,   -1e-3, 0)))
    p4 = m.Add(MeshPoint(Pnt(2e-3,    1e-3, 0)))
    p5 = m.Add(MeshPoint(Pnt(2e-3,    3e-3, 0)))

    step("Add 4 triangles (2 axis-touching, 2 off-axis)")
    # Conductor (inner): mat 2
    m.Add(Element2D(2, [p0, p3, p1]))       # axis-touching
    m.Add(Element2D(2, [p1, p3, p4]))       # off-axis (shared edge with axis-touching)
    # Air (outer): mat 1
    m.Add(Element2D(1, [p1, p4, p2]))       # axis-touching
    m.Add(Element2D(1, [p2, p4, p5]))       # off-axis

    step("Add boundary edges (only outer r-side)")
    m.Add(Element1D([p3, p4], index=1))
    m.Add(Element1D([p4, p5], index=1))
    return m


def main():
    ngsglobals.msg_level = 0
    step("Build netgen mesh")
    ngmesh = build_axis_2mat_mesh()
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

    step("Create H1Henrotte(mesh, order=2, dirichlet='outer')")
    fes = H1Henrotte(mesh, order=2, dirichlet="outer")
    step(f"ndof={fes.ndof}")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    step(f"free={n_free}")

    step("Make BFI K + assemble")
    from math import pi
    MU0 = 4 * pi * 1e-7
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()
    print(f"  K nnz: {a.mat.nze}", flush=True)

    step("Make BFI sigma-M + assemble")
    sigma_cf = mesh.MaterialCF({"conductor": 5.8e7}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m_bf.Assemble()
    print(f"  M nnz: {m_bf.mat.nze}", flush=True)

    step("DONE — no segfault")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
