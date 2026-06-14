"""Kameari + Kelvin v6: project J source to tangential subspace at conductor surface.

Key insight: J_0 = σ A_ext has J·n ≠ 0 at conductor surface, which is unphysical
for eddy current. Closed PEC's Dirichlet BC automatically enforces this via
the HCurl space; Kelvin (no Dirichlet at conductor) doesn't.

v6 fix: project J source onto HCurl(conductor; tangential at surface) =
        closed-PEC HCurl space. Use this projected J as source for the
        global Kelvin curl-curl solve.

Expected: τ_0 should now approach the radiation-BC physics (~14 μs from
mpmath Cauer Stage 0 of Y_Foster).
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting...", flush=True)

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
H_AIR = 2.0e-3
ORDER = 2
N_STAGES = 8
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
    # Tag conductor surface specifically
    return OCCGeometry(geo)


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari + Kelvin v6: J-projection to tangential HCurl ===",
          flush=True)
    print(f"Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"ELF Foster B_z tau_lead = 11.51 us", flush=True)
    print(f"mpmath Cauer Stage 0 tau_0 (target) ≈ 14.4 us\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, mats = {mesh.GetMaterials()}, "
          f"bdry = {set(mesh.GetBoundaries())}  ({time.time()-t0:.1f}s)\n",
          flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})
    A_ext_cf = CoefficientFunction((-y, x, 0)) * 0.5

    # === Two FE spaces ===
    # fes_global = full HCurl on cond+air+kelvin, periodic, GND Dirichlet
    fes_global = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                                 nograds=True))
    print(f"  fes_global ndof = {fes_global.ndof}", flush=True)

    tree_edges = build_spanning_tree(mesh)
    fd = fes_global.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes_global.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges", flush=True)

    # fes_cond = closed-PEC HCurl on conductor with tangential Dirichlet
    # Conductor surface needs to be Dirichlet. The mesh's conductor surface
    # boundary should be auto-named via NetGen's conductor material faces.
    bdry_set = set(mesh.GetBoundaries())
    cond_bdry = None
    for b in bdry_set:
        if b and ("conductor" in b.lower() or "default" in b.lower()):
            cond_bdry = b
            break
    print(f"  Available boundaries: {bdry_set}", flush=True)
    print(f"  Selected conductor boundary for Dirichlet: '{cond_bdry}'",
          flush=True)

    # Use full mesh HCurl with definedon=conductor for "conductor space"
    # but with full-mesh Dirichlet on conductor surface (= tangential J=0)
    # NOTE: NGSolve's HCurl with definedon material gives conductor-only DOFs.
    fes_cond_full = HCurl(mesh, order=ORDER,
                          definedon=mesh.Materials("conductor"),
                          dirichlet=cond_bdry if cond_bdry else "",
                          nograds=True)
    print(f"  fes_cond (cond-only HCurl, tangential Dirichlet) ndof = "
          f"{fes_cond_full.ndof}", flush=True)

    # === Project J source = sigma * A_ext to fes_cond ===
    # Use L2 projection: find gfJ_proj in fes_cond minimizing ||gfJ_proj - sigma*A_ext||
    # This automatically gives J·n=0 at conductor surface (tangential Dirichlet)
    print("\n=== Projecting J_0 = sigma * A_ext to closed-PEC HCurl ===",
          flush=True)
    gfJ_proj = GridFunction(fes_cond_full)
    try:
        gfJ_proj.Set(sigma_cf * A_ext_cf,
                     definedon=mesh.Materials("conductor"))
    except Exception as e:
        print(f"  Set failed: {e}", flush=True)
        gfJ_proj.Set(sigma_cf * A_ext_cf)

    # Compute R_0 from projected J
    R_inv_proj = float(Integrate(gfJ_proj * gfJ_proj * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_proj = 1.0 / R_inv_proj if R_inv_proj > 1e-30 else float('inf')
    print(f"  R_inv from projected J = {R_inv_proj:.6e}", flush=True)
    print(f"  R_0_projected = {R0_proj:.6e}", flush=True)
    print(f"  R_0_analytical (Parseval) = {R0_anal:.6e}", flush=True)
    print(f"  ratio R_0_proj/R_0_anal = {R0_proj/R0_anal:.6f} "
          f"(if 1, projection is conservative; if <1, projection lost mass)\n",
          flush=True)

    # === Set up Kelvin operator on fes_global ===
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))
    u, v = fes_global.TnT()
    a = BilinearForm(fes_global)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    a += GAUGE_EPS * NU_0 * u * v * dx(bonus_intorder=BONUS_INT)
    print("  Assembling+factor (full Kelvin operator)...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes_global.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes_global.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    # === Kameari iteration ===
    diag = []
    # Initialize: gfJ in fes_cond_full, gfApot in fes_global
    gfJ = GridFunction(fes_cond_full)
    gfJ.vec.data = gfJ_proj.vec.data
    gfApot = GridFunction(fes_global)
    gfApot.vec[:] = 0.0

    for n in range(N_STAGES):
        t_stage = time.time()

        # Solve a(gfA, w) = (gfJ, w) on fes_global, where (gfJ, w) integrates
        # over conductor (gfJ defined on conductor only)
        f = LinearForm(fes_global)
        f += gfJ * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes_global)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(gfJ * gfJ * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm bad ({R_inv:.2e}), STOP", flush=True)
            break
        Rn = 1.0 / R_inv

        # Update Apot
        gfApot.vec.data += Rn * gfA.vec

        # L_n = R_n × <J, A_pot>_cond
        Ln_int = float(Integrate(gfJ * gfApot * dx("conductor",
                                                    bonus_intorder=BONUS_INT),
                                 mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        tau_us = Ln/Rn*1e6 if Rn > 0 and Ln > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.3f}" if tau_us else "N/A"
        print(f"  [{n}] R={Rn:.4e} (anal {R0_anal:.4e}, ratio "
              f"{Rn/R0_anal:.4f}), L={Ln:.4e}({signL}), tau={tau_str}us "
              f"({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": Rn, "L_n": Ln, "tau_us": tau_us,
                     "energy_norm": R_inv, "sign_L": signL})

        if Ln <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break

        # Update J: J_{n+1} = J_n - sigma * Apot / Ln, projected to fes_cond
        # Direct: gfJ.vec -= (sigma_Cu/Ln) * gfApot.vec
        # But gfJ is on fes_cond_full and gfApot is on fes_global — need projection.
        # Project sigma * gfApot/Ln onto fes_cond_full:
        gfApot_proj = GridFunction(fes_cond_full)
        try:
            gfApot_proj.Set((sigma_Cu / Ln) * gfApot,
                            definedon=mesh.Materials("conductor"))
        except Exception as e:
            print(f"      Apot projection failed: {e}", flush=True)
            break
        gfJ.vec.data -= gfApot_proj.vec

    print("\n" + "=" * 78, flush=True)
    print("KAMEARI + KELVIN v6 (J-projected to tangential) SUMMARY", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'R_n':>14} {'L_n [H]':>14} {'tau [us]':>12}",
          flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] else "N/A"
        print(f"{d['n']:>3} {d['R_n']:>14.4e} {d['L_n']:>14.4e} {tau_s:>12}",
              flush=True)

    if diag and diag[0]['tau_us']:
        print(f"\n  tau_0 = {diag[0]['tau_us']:.3f} us "
              f"(target ~14.4 us, ELF lead 11.51 us)", flush=True)
        ratio = diag[0]['tau_us'] / 14.4
        print(f"  ratio tau_0 / 14.4 us = {ratio:.4f}", flush=True)

    out = {"R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
           "order": ORDER, "ne": mesh.ne, "stages": diag,
           "R0_analytic": R0_anal,
           "target_tau_0_us": 14.4, "elf_tau_lead_us": 11.51}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v6.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
