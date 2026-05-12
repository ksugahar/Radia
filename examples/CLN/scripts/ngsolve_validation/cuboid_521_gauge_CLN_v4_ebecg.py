"""5x2x1 mm Cu cuboid CLN — v4: EBE+CG instead of direct sparse Cholesky.

v3 used direct sparse Cholesky. Memory-limited at h=0.07mm
(factorization needed ~13 GB).

v4 replaces:
  - HCurl curl-curl solve: Preconditioner(a, "local") + CG  (matrix-free)
  - H1 Helmholtz-Hodge: Preconditioner(a, "local") + CG     (matrix-free)
Memory dominated by sparse matrix storage (no factorization),
so h=0.05mm with ndof~4M should fit.

Trade-offs:
  - CG iteration count grows with mesh refinement; each stage uses
    fresh RHS so CG must be re-run each time
  - Helmholtz-Hodge auxiliary solve also needs CG (same preconditioner,
    different operator)

Verification path:
  - At h=0.20mm and h=0.10mm, EBE+CG should reproduce v3 direct sparse
    to within CG tolerance (set 1e-12)
  - Then push h=0.07mm and h=0.05mm to extend Stage 5+ precision

Reference: cuboid_521_mathematica_reference.json (Nagamine high-precision).
"""
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, H1, BilinearForm, LinearForm, GridFunction, Preconditioner,
    CoefficientFunction, curl, grad, dx, x, y, z, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)
from ngsolve.krylovspace import CGSolver
from collections import deque
from math import pi
import json
import time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_box = ax * ay * az
N_STAGES = 12
H_COND = 0.10e-3   # start at 0.10 to verify EBE+CG matches v3, then push
ORDER = 3
CG_TOL = 1e-12
CG_MAXITER = 100000


def build_cuboid():
    box = Box(Pnt(0, 0, 0), Pnt(ax, ay, az))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    return OCCGeometry(box)


def classify_vertices(mesh):
    boundary_v = set()
    for el in mesh.Elements(BND):
        for v in el.vertices:
            boundary_v.add(v.nr)
    return set(range(mesh.nv)) - boundary_v, boundary_v


def build_spanning_tree_interior(mesh):
    interior_v, _ = classify_vertices(mesh)
    nv = mesh.nv
    visited = [False] * nv
    for v in range(nv):
        if v not in interior_v:
            visited[v] = True
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        if v0 in interior_v and v1 in interior_v:
            adj[v0].append((v1, ed.nr))
            adj[v1].append((v0, ed.nr))
    for start in interior_v:
        if visited[start]:
            continue
        visited[start] = True
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for vn, edn in adj[v]:
                if not visited[vn]:
                    visited[vn] = True
                    tree_edges.append(edn)
                    queue.append(vn)
    return tree_edges


def kameari_v4(mesh, n_stages, case="A"):
    fes_HC = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    fes_H1 = H1(mesh, order=ORDER, dirichlet="conductor_surface")
    print(f"  HCurl ndof = {fes_HC.ndof}, H1 ndof = {fes_H1.ndof}")

    tree_edges = build_spanning_tree_interior(mesh)
    free_HC = BitArray(fes_HC.FreeDofs())
    masked_count = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes_HC.GetDofNrs(edge)
        if dofs and free_HC[dofs[0]]:
            free_HC[dofs[0]] = False
            masked_count += 1
    interior_v, boundary_v = classify_vertices(mesh)
    print(f"  Interior vertices = {len(interior_v)}, boundary = {len(boundary_v)}")
    print(f"  Tree-cotree masked {masked_count} interior-edge DoFs")
    print(f"  Active HCurl DoFs = {sum(free_HC)}")

    u, v = fes_HC.TnT()
    uu, vv = fes_H1.TnT()

    a_HC = BilinearForm(fes_HC)
    a_HC += (1.0/mu0) * curl(u) * curl(v) * dx
    pre_HC = Preconditioner(a_HC, type="local")
    a_H1 = BilinearForm(fes_H1)
    a_H1 += grad(uu) * grad(vv) * dx
    pre_H1 = Preconditioner(a_H1, type="local")

    print("  Assembling matrices + preconditioners...")
    t0 = time.time()
    with TaskManager():
        a_HC.Assemble()
        a_H1.Assemble()
    print(f"  Assembled in {time.time()-t0:.1f}s")

    # CGSolver respects FreeDofs of the preconditioner (built from FES).
    # We need to enforce tree-cotree mask separately. NGSolve "local" pre
    # doesn't see free_HC directly; we apply tree mask by zeroing those
    # entries in solution after solve (a posteriori projection).
    # Equivalent: project gfA.vec on free_HC = True subspace.
    inv_HC = CGSolver(a_HC.mat, pre_HC.mat, maxiter=CG_MAXITER, tol=CG_TOL,
                      printrates=False)
    inv_H1 = CGSolver(a_H1.mat, pre_H1.mat, maxiter=CG_MAXITER, tol=CG_TOL,
                      printrates=False)

    if case == "A":
        J = CoefficientFunction((0, 0, sigma_Cu))
    else:
        J = CoefficientFunction((-y, x, 0)) * (sigma_Cu / 2.0)

    sigma_inv = 1.0 / sigma_Cu
    R_list, L_list, cg_iters_HC, cg_iters_H1 = [], [], [], []
    Apot = None

    for n in range(n_stages):
        # Step 1: curl-curl A_raw via EBE+CG
        f_HC = LinearForm(fes_HC)
        f_HC += J * v * dx
        with TaskManager():
            f_HC.Assemble()
        gfA_raw = GridFunction(fes_HC, name=f"A_raw_{n}")
        t0 = time.time()
        with TaskManager():
            gfA_raw.vec.data = inv_HC * f_HC.vec
        # Apply tree-cotree mask: zero non-free DoFs (project to gauge subspace)
        for d in range(fes_HC.ndof):
            if not free_HC[d]:
                gfA_raw.vec[d] = 0.0
        cg_HC = inv_HC.iterations
        t_HC = time.time() - t0

        # Step 2: Helmholtz-Hodge correction via EBE+CG
        f_H1 = LinearForm(fes_H1)
        f_H1 += grad(vv) * gfA_raw * dx
        with TaskManager():
            f_H1.Assemble()
        gf_phi = GridFunction(fes_H1, name=f"phi_{n}")
        t0 = time.time()
        with TaskManager():
            f_H1.vec.data = -f_H1.vec
            gf_phi.vec.data = inv_H1 * f_H1.vec
        cg_H1 = inv_H1.iterations
        t_H1 = time.time() - t0

        A_corr = gfA_raw + grad(gf_phi)

        R_inv = Integrate(J * J * sigma_inv * dx, mesh)
        if R_inv < 1e-25:
            print(f"[{n}] R_inv too small, stop")
            break
        Rn = 1.0 / R_inv
        if Apot is None:
            Apot = Rn * A_corr
        else:
            Apot = Apot + Rn * A_corr
        Ln = Rn * Integrate(J * Apot * dx, mesh)
        R_list.append(float(Rn))
        L_list.append(float(Ln))
        cg_iters_HC.append(cg_HC)
        cg_iters_H1.append(cg_H1)
        tau_us = Ln/Rn * 1e6 if Rn > 0 and Ln > 0 else 0
        print(f"[{n}] R={Rn:.6e}, L={Ln*1e9:.6e} nH, "
              f"CG: HC {cg_HC} ({t_HC:.1f}s), H1 {cg_H1} ({t_H1:.1f}s)")
        if Ln <= 0:
            print(f"[{n}] L<=0, stop")
            break
        J = J - sigma_Cu * Apot / Ln

    return R_list, L_list, cg_iters_HC, cg_iters_H1


def main():
    ngsglobals.msg_level = 0
    print(f"5x2x1 mm Cu cuboid (v4: EBE+CG instead of direct sparse)")
    print(f"  h = {H_COND*1000} mm, order = {ORDER}, CG tol = {CG_TOL}")
    geo = build_cuboid()
    mesh = Mesh(geo.GenerateMesh(maxh=H_COND))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}, nedge = {mesh.nedge}")

    R0_A_anal = 1 / (sigma_Cu * V_box)
    print(f"\nMathematica R_0 reference:")
    print(f"  Case A: 1/(σV) = {R0_A_anal:.16e}")

    print()
    print("=" * 60)
    print(f"CASE A (uniform J_z), EBE+CG")
    print("=" * 60)
    R_A, L_A, cg_HC, cg_H1 = kameari_v4(mesh, N_STAGES, case="A")

    out = {
        "method": "EBE+CG (Preconditioner type=local) + interior tree-cotree + Helmholtz-Hodge + Tanimoto accumulated (v4)",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "h_mm": H_COND*1000, "order": ORDER, "ne": mesh.ne,
        "cg_tol": CG_TOL,
        "Mathematica_R0": {"caseA": R0_A_anal},
        "case_A": {"R": R_A, "L": L_A,
                   "cg_iters_HC": cg_HC, "cg_iters_H1": cg_H1},
    }
    out_path = (Path(__file__).parent /
                f"cuboid_521_gauge_v4_h{int(H_COND*1000*100):03d}_results.json")
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
