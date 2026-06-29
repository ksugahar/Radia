"""tanimoto_AT_kelvin.py — Tanimoto A-T with Kelvin two-sphere domain.

DIAGNOSTIC TEST 2026-05-09: Tanimoto+H-H on conductor-only mesh gives
225.77 us for cylinder (vs BEM/v5/axifem 4-way consensus 218.69 us,
+3% gap). Hypothesis: gap is caused by missing vacuum coupling.

This script ports v5's Kelvin two-sphere mesh setup (cylinder + inner
air sphere + Kelvin exterior sphere with periodic identification + GND
vertex) to the Tanimoto A-T alternating Poisson iteration.

  fesA = HCurl(full mesh, nograds=True, dirichlet_bbnd="GND") with tree-cotree
  fesT = HCurl on conductor only, nograds=True
  K_A = ∫(nu_cf) curl·curl on full mesh (with Kelvin transformation)
  K_T = ∫(1/σ) curl·curl on conductor only
  Magnetic energy integrates nu_cf * |B|² over the full Kelvin domain.
  Source: H-H projected A_s on conductor, Galerkin volume integral

Predicted outcome: stage 1 τ ≈ 218-219 us if vacuum coupling is the cause
of the 3% gap.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import deque
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")

from netgen.occ import Cylinder, Box, Sphere, Pnt, X, Y, Z, OCCGeometry
from ngsolve import (
    Mesh, HCurl, H1, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, grad, dx, x, y, z,
    Integrate, InnerProduct, TaskManager, ngsglobals,
)

from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import (
    make_kelvin_nu_cf,
    make_reduced_potential_background_cf,
    NU_0,
)


# === Physical parameters ====================================================
mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
B0 = 1.0

# Default cylinder
R_DISK = 10e-3
T_DISK = 2e-3

# Kelvin parameters (match v5)
R_K_DEFAULT = 50e-3
OFFSET_DEFAULT = (5 * R_K_DEFAULT, 0.0, 0.0)

# FE parameters
ORDER = 2
H_COND = 1.0e-3
H_AIR = 12.0e-3
N_STAGES = 4
BONUS_INT = 8


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


GEOMETRIES = {
    "cylinder":   {"kind": "cyl",  "R": 10e-3, "t": 2e-3},
    "A1_square":  {"kind": "rect", "Lx": 17.72e-3, "Ly": 17.72e-3,
                    "Lz": 2e-3},   # disk-like, V/t/area = cylinder R=10 t=2
    "cuboid_521": {"kind": "rect", "Lx": 5e-3, "Ly": 2e-3, "Lz": 1e-3},
    "sphere":     {"kind": "sphere", "R": 10e-3},
}


def build_geo(geom_name, R_K, OFFSET):
    p = GEOMETRIES[geom_name]
    if p["kind"] == "cyl":
        cond = Cylinder(Pnt(0, 0, -p["t"] / 2), Z, r=p["R"], h=p["t"])
    elif p["kind"] == "sphere":
        cond = Sphere(Pnt(0, 0, 0), p["R"])
    else:
        Lx, Ly, Lz = p["Lx"], p["Ly"], p["Lz"]
        cond = Box(Pnt(-Lx/2, -Ly/2, -Lz/2),
                   Pnt( Lx/2,  Ly/2,  Lz/2))
    cond.mat("conductor")
    cond.maxh = H_COND
    inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_sphere.maxh = H_AIR
    inner_air = inner_sphere - cond
    inner_air.name = "air"
    inner_air.maxh = H_AIR
    geo, _ = add_kelvin_exterior_domain(
        [inner_air, cond], offset=OFFSET, R_K=R_K, inner_maxh=H_AIR,
    )
    return OCCGeometry(geo)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", choices=list(GEOMETRIES.keys()),
                         default="cylinder")
    parser.add_argument("--stages", type=int, default=N_STAGES)
    parser.add_argument("--order", type=int, default=ORDER)
    parser.add_argument("--out-tag", type=str, default=None)
    args = parser.parse_args()

    n_stages = args.stages
    ORDER_use = args.order
    R_K = R_K_DEFAULT
    OFFSET = OFFSET_DEFAULT

    ngsglobals.msg_level = 0
    print("=" * 78, flush=True)
    print(" Tanimoto A-T + Kelvin two-sphere (vacuum coupling diagnostic)",
          flush=True)
    print("=" * 78, flush=True)
    print(f"  Cylinder R={R_DISK*1000} mm, t={T_DISK*1000} mm,"
          f" σ={sigma_Cu:.2e} S/m", flush=True)
    print(f"  Kelvin: R_K = {R_K*1000} mm, OFFSET = {OFFSET}", flush=True)
    print(f"  ORDER={ORDER_use} HCurl, H_COND={H_COND*1000}mm,"
          f" H_AIR={H_AIR*1000}mm, N_STAGES={n_stages}",
          flush=True)
    print()

    print(" Building geometry + mesh ...", flush=True)
    geo = build_geo(args.geom, R_K, OFFSET)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER_use)
    print(f"  GenerateMesh + Curve(p={ORDER_use}): {time.time()-t0:.1f} s",
          flush=True)
    print(f"  ne={mesh.ne}, nv={mesh.nv}", flush=True)
    print(f"  materials={mesh.GetMaterials()}", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0,
                                 "kelvin": 0.0})
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                               kelvin_mats=("kelvin",))

    A_s_cf_raw = make_reduced_potential_background_cf(
        mesh,
        F_inner_factory=lambda xc, yc, zc: CoefficientFunction(
            (-yc, xc, 0)) * (B0 / 2.0),
        R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3,
    )

    # Spaces -------------------------------------------------------------
    fesA = Periodic(HCurl(mesh, order=ORDER_use, dirichlet_bbnd="GND",
                           nograds=True))
    print(f"  fesA HCurl ndof = {fesA.ndof}", flush=True)

    fesT = HCurl(mesh, order=ORDER_use, nograds=True,
                  definedon=mesh.Materials("conductor"))
    print(f"  fesT HCurl(cond only) ndof = {fesT.ndof}"
          f" active = {sum(fesT.FreeDofs())}",
          flush=True)

    tree_edges = build_spanning_tree(mesh)
    fdA = fesA.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fesA.GetDofNrs(edge)
        if dofs and fdA[dofs[0]]:
            fdA[dofs[0]] = False
            masked += 1
    print(f"  fesA tree-cotree masked {masked}, active = {sum(fdA)}",
          flush=True)

    A_tr, A_te = fesA.TnT()
    T_tr, T_te = fesT.TnT()

    # σ-weighted H-H projection of A_s globally
    # User instruction: "σ weighted, div(σA)=0 を満たすべき"
    # For uniform σ in conductor and σ=0 in air/Kelvin, this naturally
    # restricts to conductor (σ-factor in air/Kelvin makes contribution 0)
    # and enforces div(σA_s')=0 in cond, σA_s'·n=0 on ∂cond.
    # Use fes_phi on full mesh with σ_cf weighting; air DOFs are pinned.
    print(" σ-weighted Helmholtz-Hodge projection of A_s ...",
          flush=True)
    fes_phi_proj = H1(mesh, order=ORDER_use + 1,
                       definedon=mesh.Materials("conductor"))
    phi_tr, phi_te = fes_phi_proj.TnT()
    a_phi_proj = BilinearForm(fes_phi_proj)
    a_phi_proj += sigma_cf * grad(phi_tr) * grad(phi_te) * dx(
        "conductor", bonus_intorder=BONUS_INT)
    f_phi_proj = LinearForm(fes_phi_proj)
    f_phi_proj += sigma_cf * A_s_cf_raw * grad(phi_te) * dx(
        "conductor", bonus_intorder=BONUS_INT)
    with TaskManager():
        a_phi_proj.Assemble()
        f_phi_proj.Assemble()
        fd_phi = fes_phi_proj.FreeDofs()
        pinned_phi = -1
        for ip in range(fes_phi_proj.ndof):
            if fd_phi[ip]:
                fd_phi[ip] = False
                pinned_phi = ip
                break
        if pinned_phi >= 0:
            f_phi_proj.vec[pinned_phi] = 0.0
        try:
            inv_phi_proj = a_phi_proj.mat.Inverse(fd_phi, inverse="pardiso")
        except Exception:
            inv_phi_proj = a_phi_proj.mat.Inverse(
                fd_phi, inverse="sparsecholesky")
        gf_phi_proj = GridFunction(fes_phi_proj)
        gf_phi_proj.vec.data = inv_phi_proj * f_phi_proj.vec
    A_s_cf = A_s_cf_raw - grad(gf_phi_proj)
    print(f"  pinned phi vertex {pinned_phi}", flush=True)

    # Setup σ-weighted H-H operator for T-iterates (clean higher-order
    # grad-zero modes in HCurl(ORDER>=2)). Reuses fes_phi_proj structure.
    # For each T_n: solve ∫σ ∇φ ∇q = ∫σ T_n · ∇q on cond, then T_n' = T_n - ∇φ
    def hodge_project_T(gfT_in):
        l_phi = LinearForm(fes_phi_proj)
        l_phi += sigma_cf * gfT_in * grad(phi_te) * dx(
            "conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            l_phi.Assemble()
        if pinned_phi >= 0:
            l_phi.vec[pinned_phi] = 0.0
        gf_phi_T = GridFunction(fes_phi_proj)
        with TaskManager():
            gf_phi_T.vec.data = inv_phi_proj * l_phi.vec
        return grad(gf_phi_T)  # CoefficientFunction

    # Stage 0: J_0 = σA_s' directly (skip K_T solve for T_0) ------------
    # User instruction (2026-05-10):
    #   "T_0 を skip して J_0 = σA_s' で initialize"
    # Physically: at t=0+ after step B applied, instantaneous eddy current
    # is J(0+) = σ × E_induced = σ × A_s (after H-H projection)
    # R_0 = 1/∫|J_0|²/σ = 1/∫σ|A_s|² (initial Joule loss per unit current²)
    print("\n Stage 0: J_0 = σ A_s' (direct, no T_0 solve)", flush=True)
    J_acc = sigma_cf * A_s_cf   # CoefficientFunction directly
    R_inv_0 = float(Integrate(
        sigma_cf * A_s_cf * A_s_cf * dx(
            "conductor", bonus_intorder=BONUS_INT), mesh))
    R_0 = 1.0 / R_inv_0
    print(f"  R_inv_0 = ∫σ|A_s'|² = {R_inv_0:.6e},  R_0 = {R_0:.6e}",
          flush=True)

    # Still need fes_T inversion for subsequent T equations (stage 1+)
    a_T = BilinearForm(fesT, check_unused=False)
    a_T += (1.0 / sigma_Cu) * curl(T_tr) * curl(T_te) * dx(
        "conductor", bonus_intorder=BONUS_INT)
    t0 = time.time()
    with TaskManager():
        a_T.Assemble()
        try:
            inv_T = a_T.mat.Inverse(fesT.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_T = a_T.mat.Inverse(fesT.FreeDofs(),
                                     inverse="sparsecholesky")
    print(f"  K_T assemble + factor (for stage 1+):"
          f" {time.time()-t0:.1f} s", flush=True)

    Rn = [R_0]
    Ln = []
    diag = []
    R_curr = R_0

    # Setup A bilinear form (with Kelvin nu_cf) -------------------------
    a_A = BilinearForm(fesA, check_unused=False)
    a_A += nu_cf * curl(A_tr) * curl(A_te) * dx(bonus_intorder=BONUS_INT)
    print(" Assembling K_A (Kelvin) ...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_A.Assemble()
        try:
            inv_A = a_A.mat.Inverse(fdA, inverse="pardiso")
            inverse_used_A = "pardiso"
        except Exception as e:
            print(f"  Pardiso failed ({e}), sparsecholesky", flush=True)
            inv_A = a_A.mat.Inverse(fdA, inverse="sparsecholesky")
            inverse_used_A = "sparsecholesky"
    print(f"  K_A assemble + factor ({inverse_used_A}):"
          f" {time.time()-t0:.1f} s", flush=True)

    raise RuntimeError("Historical CLN field-stage reduction is retired.")


if __name__ == "__main__":
    main()
