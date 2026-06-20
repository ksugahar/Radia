"""2D rectangular bar via 3D thin slab + HCurl + tree-cotree gauge.

Demonstrates that the tree-cotree machinery from 3D works identically in
quasi-2D — solving the same problem as rect_2D_kameari.py (H1 scalar)
to high precision.

Setup:
  - 3D mesh: a x b x dz (dz << a, b makes it quasi-2D)
  - HCurl elements (vector A)
  - Dirichlet on x, y faces (A × n = 0 → A_z = 0 there)
  - Free on z faces (A_z unconstrained, allowing z-independence)
  - Tree-cotree gauge masks gradient kernel
  - Tanimoto pattern: accumulated R-weighted potential

The dz length cancels out when extracting per-unit-length R, L:
  R_per_length = R_3D * dz
  L_per_length = L_3D * dz
"""
import json
from collections import deque
from math import pi
from pathlib import Path

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, BitArray,
    Integrate, TaskManager, ngsglobals,
)

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
DZ = 0.1e-3        # thin slab thickness in z (must be small for quasi-2D)
H_MESH = 0.10e-3
ORDER = 2
N_STAGES = 12


def build_slab():
    box = Box(Pnt(0, 0, 0), Pnt(ax, ay, DZ))
    box.mat("conductor")
    # Name only x and y faces so we Dirichlet only those
    for face in box.faces:
        n = face.center
        # Determine which face by checking normal-aligned coordinate extreme
        if abs(n.x - 0) < 1e-9:
            face.bc("xleft")
        elif abs(n.x - ax) < 1e-9:
            face.bc("xright")
        elif abs(n.y - 0) < 1e-9:
            face.bc("yleft")
        elif abs(n.y - ay) < 1e-9:
            face.bc("yright")
        elif abs(n.z - 0) < 1e-9:
            face.bc("zbot")
        elif abs(n.z - DZ) < 1e-9:
            face.bc("ztop")
    box.maxh = H_MESH
    return OCCGeometry(box)


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


def kameari_treecotree(mesh, n_stages):
    """HCurl 3D + tree-cotree + Tanimoto. Quasi-2D via thin slab."""
    fes = HCurl(mesh, order=ORDER,
                dirichlet="xleft|xright|yleft|yright",
                nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    tree_edges = build_spanning_tree(mesh)
    free = BitArray(fes.FreeDofs())
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and free[dofs[0]]:
            free[dofs[0]] = False
            masked += 1
    print(f"  Tree masked {masked}, active = {sum(free)}")

    u, v = fes.TnT()
    sigma_inv = 1.0 / sigma
    J = CoefficientFunction((0, 0, sigma))
    Apot = None

    R_list, L_list = [], []
    for k in range(n_stages):
        a = BilinearForm(fes)
        a += (1.0 / mu0) * curl(u) * curl(v) * dx
        f = LinearForm(fes)
        f += J * v * dx

        with TaskManager():
            a.Assemble()
            f.Assemble()

        gfA = GridFunction(fes, name=f"A_{k}")
        with TaskManager():
            inv = a.mat.Inverse(freedofs=free, inverse="sparsecholesky")
            gfA.vec.data = inv * f.vec

        R_inv = Integrate(J * J * sigma_inv * dx, mesh)
        if R_inv < 1e-30:
            print(f"  R_inv too small, stop")
            break
        Rk = 1.0 / R_inv

        if Apot is None:
            Apot = Rk * gfA
        else:
            Apot = Apot + Rk * gfA

        Lk = Rk * Integrate(J * Apot * dx, mesh)
        # Convert from 3D total to per-unit-length z
        Rk_per_len = Rk * DZ
        Lk_per_len = Lk * DZ

        R_list.append(Rk_per_len)
        L_list.append(Lk_per_len)
        tau_us = Lk_per_len / Rk_per_len * 1e6 if Rk_per_len > 0 and Lk_per_len > 0 else 0
        print(f"[{k}] R = {Rk_per_len:.10e} Ω/m,  L = {Lk_per_len*1e9:.10e} nH/m,  τ = {tau_us:.4f} µs")

        if Lk <= 0:
            print("  L<=0, stop")
            break

        J = J - sigma * Apot / Lk

    return R_list, L_list


def main():
    ngsglobals.msg_level = 0
    print(f"Quasi-2D rect bar {ax*1000}x{ay*1000} mm Cu (z-thin slab dz={DZ*1000} mm)")
    print(f"  h={H_MESH*1000} mm, order={ORDER}, sigma={sigma:.2e}")

    geo = build_slab()
    mesh = Mesh(geo.GenerateMesh(maxh=H_MESH))
    print(f"  ne={mesh.ne}, nv={mesh.nv}, nedge={mesh.nedge}")
    print(f"  boundaries: {mesh.GetBoundaries()}")

    R0_anal = 1 / (sigma * ax * ay)
    print(f"  R_DC analytic = 1/(σ·a·b) = {R0_anal:.10e} Ω/m")
    print()

    R, L = kameari_treecotree(mesh, N_STAGES)

    print()
    print(f"Summary ({len(R)} stages, HCurl + tree-cotree, quasi-2D)")
    print(f"  R_0 / R_DC = {R[0]/R0_anal:.10f}")

    out = {
        "method": "Quasi-2D (thin 3D slab) HCurl + tree-cotree + Tanimoto",
        "geometry_mm": [ax * 1000, ay * 1000, DZ * 1000],
        "h_mm": H_MESH * 1000,
        "order": ORDER,
        "ne": mesh.ne,
        "R_DC_analytic": R0_anal,
        "R_per_length": [float(r) for r in R],
        "L_per_length": [float(l) for l in L],
    }
    out_path = Path(__file__).parent / "rect_2D_kameari_treecotree_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
