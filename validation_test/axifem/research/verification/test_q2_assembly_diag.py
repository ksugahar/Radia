"""Diagnose Q2 multi-element assembly: build a tiny 2-quad mesh and dump the
9x9 element matrix from each element + the assembled global matrix vs a
manual reference.

Goal: find the source of the factor-~8 error in test_hiruma_disk_q2.py.
"""
import sys
from pathlib import Path
import numpy as np

REPO = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "tests" / "axifem" / "_reference_python").is_dir()
)
sys.path.insert(0, str(REPO / "tests" / "axifem" / "_reference_python"))

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager,
)
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


def build_two_quad_mesh():
    """Two stacked axis-aligned quads: r in [1mm, 2mm], z in [0, 1mm] and
    z in [1mm, 2mm]. Material: 'cu' for both."""
    ngmesh = NgMesh()
    ngmesh.dim = 2
    ngmesh.SetMaterial(1, "cu")
    fd = ngmesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    ngmesh.SetBCName(0, "bnd")
    pts = [(1e-3, 0), (2e-3, 0), (1e-3, 1e-3), (2e-3, 1e-3),
           (1e-3, 2e-3), (2e-3, 2e-3)]
    pids = [ngmesh.Add(MeshPoint(Pnt(r, z, 0))) for r, z in pts]
    # Bottom quad: vertices in CCW (ra,za)(rb,za)(rb,zb)(ra,zb)
    ngmesh.Add(Element2D(1, [pids[0], pids[1], pids[3], pids[2]]))
    # Top quad
    ngmesh.Add(Element2D(1, [pids[2], pids[3], pids[5], pids[4]]))
    return Mesh(ngmesh)


def main():
    mesh = build_two_quad_mesh()
    print(f"mesh: nv={mesh.nv} ne={mesh.ne} nedge={mesh.nedge}")

    print("\nElement / Vertex / Edge maps (raw NGSolve):")
    for el in mesh.Elements():
        print(f"  elem {el.nr} type={el.type}  "
              f"verts={[v.nr for v in el.vertices]}  "
              f"edges={[e.nr for e in el.edges]}")
        for v in el.vertices:
            p = mesh[v].point
            print(f"     v{v.nr} = ({p[0]*1e3:.3f}, {p[1]*1e3:.3f}) mm")
        for e in el.edges:
            ev = list(mesh[e].vertices)
            p0 = mesh[ev[0]].point
            p1 = mesh[ev[1]].point
            mid = ((p0[0]+p1[0])*0.5e3, (p0[1]+p1[1])*0.5e3)
            print(f"     e{e.nr} = v{ev[0].nr}-v{ev[1].nr} mid=({mid[0]:.3f}, {mid[1]:.3f}) mm")

    fes = H1Henrotte(mesh, order=2)
    print(f"\nfes.ndof = {fes.ndof} (expected nv+ne+nq = {mesh.nv}+{mesh.nedge}+{mesh.ne} = "
          f"{mesh.nv + mesh.nedge + mesh.ne})")

    print("\nGetDofNrs per element:")
    for ei in [(0, 'VOL'), (1, 'VOL')]:
        from ngsolve.comp import ElementId, VOL
        eid = ElementId(VOL, ei[0])
        dnums = fes.GetDofNrs(eid)
        print(f"  elem {ei[0]} dnums = {list(dnums)}")

    # Assemble K and dump
    mu_cf = CoefficientFunction(4*np.pi*1e-7)
    sig_cf = CoefficientFunction(5.8e7)

    a = BilinearForm(fes, symmetric=False)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()
    K = np.zeros((fes.ndof, fes.ndof))
    rs, cs, vs = a.mat.COO()
    for r, c, v in zip(rs, cs, vs):
        K[r, c] += v

    m = BilinearForm(fes, symmetric=False)
    m += AxiHenrotteSigmaMassBFI(sig_cf)
    with TaskManager(): m.Assemble()
    M = np.zeros((fes.ndof, fes.ndof))
    rs, cs, vs = m.mat.COO()
    for r, c, v in zip(rs, cs, vs):
        M[r, c] += v

    print(f"\nK trace = {np.trace(K):.6e}")
    print(f"K Frob norm = {np.linalg.norm(K):.6e}")
    print(f"K symm? max|K-K^T| = {np.max(np.abs(K - K.T)):.3e}")

    print(f"\nM trace = {np.trace(M):.6e}")
    print(f"M Frob norm = {np.linalg.norm(M):.6e}")

    # Compare with single-element K from Python prototype, summed appropriately
    from axifem_quad_q2 import (
        element_matrices_q2_numerical, element_sigma_mass_q2_numerical,
    )
    K1, _ = element_matrices_q2_numerical(1e-3, 2e-3, 0, 1e-3, 4*np.pi*1e-7, 4*np.pi*1e-7)
    K2, _ = element_matrices_q2_numerical(1e-3, 2e-3, 1e-3, 2e-3, 4*np.pi*1e-7, 4*np.pi*1e-7)
    M1 = element_sigma_mass_q2_numerical(1e-3, 2e-3, 0, 1e-3, 5.8e7)
    M2 = element_sigma_mass_q2_numerical(1e-3, 2e-3, 1e-3, 2e-3, 5.8e7)

    print(f"\nPython per-element K1 trace = {np.trace(K1):.6e}")
    print(f"Python per-element K2 trace = {np.trace(K2):.6e}")
    print(f"Sum: {np.trace(K1)+np.trace(K2):.6e}")
    print(f"NGSolve assembled K trace  = {np.trace(K):.6e}")
    # If duplicate vertex/edge DOFs are summed correctly, NGSolve trace should be
    # less than (or equal to) sum of per-element traces because shared DOFs only
    # add once on the diagonal... actually they ADD on the diagonal since shared
    # vertex/edge DOFs accumulate.

    print(f"\nPython per-element M1 trace = {np.trace(M1):.6e}")
    print(f"Python per-element M2 trace = {np.trace(M2):.6e}")
    print(f"Sum: {np.trace(M1)+np.trace(M2):.6e}")
    print(f"NGSolve assembled M trace  = {np.trace(M):.6e}")


if __name__ == "__main__":
    main()
