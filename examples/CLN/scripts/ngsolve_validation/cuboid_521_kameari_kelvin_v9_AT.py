"""Kameari + Kelvin v9: Tanimoto A-T formulation adapted for radiation BC.

Key insight from Tanimoto (CLN_AT.ipynb): Stage 0 source is a BOUNDARY
integral on conductor surface, NOT a volumetric J = sigma * A_ext.

Volumetric A_ext = (B_0/2)(-y, x, 0) is unbounded at infinity, breaking
Kelvin pullback. But the boundary version E_s = A_ext on conductor faces
is BOUNDED (|A_ext| finite on cuboid surface). So Stage 0 driving
function is well-defined regardless of Kelvin.

Tanimoto's structure:
  Stage 0: solve T (auxiliary HCurl in conductor) from boundary E_s
           J_0 = curl T  in conductor
  The subsequent historical CLN stage update has been retired.

For Kelvin compatibility:
  fesT = HCurl in conductor only, dirichlet on conductor surface (PEC)
  fesA = Periodic(HCurl on full domain, dirichlet=GND vertex/sphere)
  The historical CLN accumulator and inductance recipe have been retired.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v9 (A-T formulation + Kelvin)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, ds, x, y, z, Cross, specialcf,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import make_kelvin_nu_cf, NU_0
from collections import deque
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

R_K = 25e-3
OFFSET = (2.5 * R_K, 0, 0)
H_COND = 0.5e-3
H_AIR = 2.0e-3
ORDER = 2
N_STAGES = 6
BONUS_INT = 8
GAUGE_EPS = 1e-8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


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


def build_geo():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.maxh = H_AIR

    inner_air = sphere_inner - cuboid
    inner_air.name = "air"

    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid], offset=OFFSET, R_K=R_K, inner_maxh=H_AIR)
    return OCCGeometry(geo)


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari + Kelvin v9: Tanimoto A-T formulation ===", flush=True)
    print(f"  Boundary E_s source on conductor surface (not volumetric)",
          flush=True)
    print(f"  Stage 0: T from boundary E, J_0 = curl T", flush=True)
    print(f"  Stage n>0: A from body J_n on full Kelvin domain\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, mats = {mesh.GetMaterials()}, "
          f"bdry = {set(mesh.GetBoundaries())} ({time.time()-t0:.1f}s)\n",
          flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})

    # E_s = A_ext at conductor surface (boundary source for T solve)
    # For uniform B_z: A_ext = (B_0/2)(-y, x, 0), bounded on cuboid surface
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    bdry_set = set(mesh.GetBoundaries())
    cond_bdry = None
    for b in bdry_set:
        if b and "default" in b.lower():
            cond_bdry = b
            break
    print(f"  Conductor boundary detected: '{cond_bdry}'", flush=True)

    # ============================================================
    # Stage 0: solve T (HCurl in conductor with surface PEC)
    # to get J_0 = curl T as the Kameari starting current
    # ============================================================
    print("=== Stage 0: T solve with boundary E source ===", flush=True)
    # Tanimoto: T has NO Dirichlet (free space; boundary integral drives it)
    fesT = HCurl(mesh, order=ORDER,
                 definedon=mesh.Materials("conductor"),
                 nograds=True)
    print(f"  fesT (cond-only) ndof = {fesT.ndof}", flush=True)

    T, W = fesT.TnT()
    n_normal = specialcf.normal(mesh.dim)
    a_T = BilinearForm(fesT)
    a_T += (1.0/sigma_Cu) * curl(T) * curl(W) * dx("conductor",
                                                    bonus_intorder=BONUS_INT)
    # Gauge regularization for T (no Dirichlet → kernel of curl)
    a_T += GAUGE_EPS * (1.0/sigma_Cu) * T * W * dx("conductor",
                                                    bonus_intorder=BONUS_INT)
    f_T = LinearForm(fesT)
    # Tanimoto's source: -Cross(E_s, W.Trace()) * n * ds(conductor surface)
    f_T += -Cross(A_ext, W.Trace()) * n_normal * ds(cond_bdry,
                                                     bonus_intorder=BONUS_INT)

    print("  Assembling+factoring T system...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_T.Assemble()
        f_T.Assemble()
        try:
            inv_T = a_T.mat.Inverse(fesT.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_T = a_T.mat.Inverse(fesT.FreeDofs(), inverse="sparsecholesky")
    gfT = GridFunction(fesT)
    with TaskManager():
        gfT.vec.data = inv_T * f_T.vec
    print(f"  T solved ({time.time()-t0:.1f}s)", flush=True)

    # J_0 = curl T (in conductor)
    J_0_cf = curl(gfT)
    R_inv_0 = float(Integrate(J_0_cf * J_0_cf * sigma_inv_cf
                               * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    R_0 = 1.0 / R_inv_0 if R_inv_0 > 1e-30 else float('inf')
    print(f"  R_0 = {R_0:.4e} (vs Parseval analytical {R0_anal:.4e}, "
          f"ratio {R_0/R0_anal:.4f})", flush=True)

    # ============================================================
    # Stage 0+: Kelvin A solve with body J source from T
    # ============================================================
    print("\n=== Setting up global Kelvin A operator ===", flush=True)
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))
    fesA = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                          nograds=True))
    print(f"  fesA (full Kelvin) ndof = {fesA.ndof}", flush=True)
    tree_edges = build_spanning_tree(mesh)
    fd = fesA.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fesA.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges", flush=True)

    A_var, N_test = fesA.TnT()
    a_A = BilinearForm(fesA)
    a_A += nu_cf * curl(A_var) * curl(N_test) * dx(bonus_intorder=BONUS_INT)
    a_A += GAUGE_EPS * NU_0 * A_var * N_test * dx(bonus_intorder=BONUS_INT)

    print("  Assembling+factoring A system...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_A.Assemble()
        try:
            inv_A = a_A.mat.Inverse(fesA.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_A = a_A.mat.Inverse(fesA.FreeDofs(), inverse="sparsecholesky")
    print(f"  A operator factored ({time.time()-t0:.1f}s)\n", flush=True)

    # ============================================================
    # Kameari iteration
    # ============================================================
    diag = []
    # J_n stored as a CF chain; initial J_0 = curl(gfT) extended by zero
    J_cf = curl(gfT) * sigma_cf / sigma_Cu   # ensures J=0 outside conductor
    # Apot accumulator
    gfApot = GridFunction(fesA)
    gfApot.vec[:] = 0.0

    raise RuntimeError("The historical CLN field-reduction implementation was retired.")


if __name__ == "__main__":
    main()
