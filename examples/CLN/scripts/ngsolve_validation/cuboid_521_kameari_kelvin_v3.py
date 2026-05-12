"""Kameari + Kelvin v3: explicit tree-cotree gauge mask + nograds=True
(in addition to the production solver's gauge_eps + Periodic + GND vertex).

The single-vertex Dirichlet + 1e-8 nu_0 * u*v gauge regularization is too
weak for our cuboid + sigma A_ext source — Kelvin solve produced |A| ~ 10^13
(L_0 ~ 10^15 H, should be ~ 33 H).

Adding a BFS spanning-tree mask on first-DoF edges + nograds=True kills
the lowest-order Whitney gradient kernel + higher-order gradient kernel
explicitly, which is what was needed in my air-box Kameari version
(cuboid_521_vacuum_kameari_breakdown.py).
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
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
H_AIR = 4.0e-3
ORDER = 3
N_STAGES = 8
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


def kameari_kelvin_v3(mesh, n_stages):
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))
    # nograds + tree-cotree (was the only working config; restore)
    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                         nograds=True))
    print(f"  Periodic HCurl(nograds=True) ndof = {fes.ndof}")
    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} edges, active = {sum(fd)}")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0, "kelvin": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=8)  # was 4
    a += GAUGE_EPS * NU_0 * u * v * dx(bonus_intorder=8)

    print("  Assembling...")
    t0 = time.time()
    with TaskManager():
        a.Assemble()
    print(f"  Assembled in {time.time()-t0:.1f}s")
    print("  Factorizing...")
    t0 = time.time()
    try:
        with TaskManager():
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        print(f"  Pardiso: {time.time()-t0:.1f}s")
    except Exception:
        with TaskManager():
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        print(f"  SparseCholesky: {time.time()-t0:.1f}s")

    diag = []
    J_cf = sigma_cf * A_ext
    Apot = None

    for n in range(n_stages):
        print(f"\n[Stage {n}]")
        t0 = time.time()
        f = LinearForm(fes)
        f += J_cf * v * dx("conductor", bonus_intorder=8)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec
        # Diagnostic: ||A||_L2
        normA2 = float(Integrate(gfA * gfA * dx, mesh))

        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx("conductor", bonus_intorder=8), mesh))
        if abs(R_inv) < 1e-30:
            break
        Rn = 1.0 / R_inv

        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA

        Ln_int = float(Integrate(J_cf * Apot * dx("conductor", bonus_intorder=8), mesh))
        Ln = Rn * Ln_int

        tau_us = Ln/Rn * 1e6 if Rn > 0 and Ln > 0 else float('nan')
        signL = "+" if Ln > 0 else "-" if Ln < 0 else "0"
        print(f"  ||A||_L2^2 = {normA2:.4e}")
        print(f"  R = {Rn:.4e}, L = {Ln:.4e} H ({signL}, tau = {tau_us:.4f} us)")
        print(f"  ({time.time()-t0:.1f}s)")

        diag.append({"n": n, "R_n": Rn, "L_n": Ln,
                     "tau_us": tau_us if not (tau_us != tau_us) else None,
                     "energy_norm": R_inv, "normA2": normA2,
                     "sign_L": signL})

        if Ln <= 0 or Rn <= 0:
            print(f"  ** SIGN FLIP **")
            if n >= 2:
                break

        J_cf = J_cf - sigma_cf * Apot / Ln

    return diag


def main():
    ngsglobals.msg_level = 1
    print(f"=== Kameari + Kelvin v3 (nograds=True + tree-cotree) ===")
    print(f"  R_K = {R_K*1000} mm, offset = {[o*1000 for o in OFFSET]} mm")
    print(f"  h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, order = {ORDER}\n")

    geo = build_geo()
    print("Generating mesh (grading=0.5)...")
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, mats = {mesh.GetMaterials()}  ({time.time()-t0:.1f}s)\n")
    print(f"Reference R_0 = {R0_anal:.6e}")
    print(f"ELF Foster B_z tau_lead = 11.51 us, expected L_0 ~ 33 H\n")

    diag = kameari_kelvin_v3(mesh, N_STAGES)

    print("\n" + "=" * 78)
    print("KAMEARI + KELVIN v3 SUMMARY")
    print("=" * 78)
    print(f"{'n':>3} {'R_n':>14} {'L_n [H]':>14} {'tau [us]':>14} {'sign(L)':>8}")
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] else " N/A"
        print(f"{d['n']:>3} {d['R_n']:>14.4e} {d['L_n']:>14.4e} {tau_s:>14} {d['sign_L']:>8}")

    if diag and diag[0]['tau_us']:
        tau0 = diag[0]['tau_us']
        elf = 11.51
        print(f"\n=== Head-to-head ===")
        print(f"  Kameari + Kelvin v3  tau_0 = {tau0:.3f} us")
        print(f"  Kameari + air-box    tau_0 = 104.37 us")
        print(f"  Vector Fit + ELF     tau   = {elf:.3f} us")
        print(f"  Kelvin / ELF = {tau0/elf:.4f}")

    out = {
        "method": "Kameari + Kelvin v3 (nograds + tree-cotree)",
        "R_K_mm": R_K*1000, "OFFSET_mm": [o*1000 for o in OFFSET],
        "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000, "order": ORDER,
        "ne": mesh.ne, "stages": diag,
        "elf_foster_tau_lead_us": 11.51,
    }
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v3.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
