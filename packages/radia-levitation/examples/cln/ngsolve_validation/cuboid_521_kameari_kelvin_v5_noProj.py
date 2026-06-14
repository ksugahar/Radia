"""Kameari + Kelvin v5: NO L2 projection of J. Use CoefficientFunction throughout.

v3/v4 had R_0 = 4.5e6 (vs analytical Parseval 2.85e6, 60% over). The closed-PEC
isolation test had R_0 = 2.85e6 exact. Difference: closed PEC has no material
interface, but Kelvin (cond/air/kelvin) has TWO interfaces (cond-air, air-kelvin
via periodic). gfJ.Set L2-projects sigma_cf * A_ext (which jumps at cond-air
boundary) onto continuous HCurl basis -> 60% energy norm loss.

v5 fix: keep J as a CF chain, never project. Apot remains a GridFunction for
the iteration update, but cross-products and energy norms use CF expressions
directly.

Hypothesis: this should recover R_0 = 2.85e6 (analytical Parseval) and
allow the iteration to find physically-correct tau values.
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
    return OCCGeometry(geo)


def kameari_kelvin_v5(mesh, n_stages):
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))
    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                         nograds=True))
    print(f"  HCurl ndof = {fes.ndof}", flush=True)

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0, "kelvin": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    a += GAUGE_EPS * NU_0 * u * v * dx(bonus_intorder=BONUS_INT)
    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s", flush=True)

    diag = []

    # J as CF chain (never projected)
    J_cf = sigma_cf * A_ext   # initial impressed source

    # Apot as GridFunction accumulator
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    # History of accumulated Apot/Ln chunks (for cross-product reconstruction)
    # J_n = J_0 - sigma * (sum_{m<n} Apot_m / L_m)
    # We store Apot_history (cumulative gfApot snapshots)
    # and L_history (the L_n values)
    # Then <J_n, J_m/sigma>_cond can be reconstructed by expanding the CF tree
    # But this gets complex. For diagnostic purpose, also project to gfJ for tracking only.

    for n in range(n_stages):
        t_stage = time.time()

        # Solve curl-curl A = J source (CF, no projection)
        f = LinearForm(fes)
        f += J_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        # R_n from analytical CF (NO projection error)
        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf
                                 * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm bad, STOP", flush=True)
            break
        Rn = 1.0 / R_inv

        # Update Apot: gfApot += Rn * gfA
        gfApot.vec.data += Rn * gfA.vec

        # L_n = R_n * <J_cf, gfApot>_cond  (gfApot is GridFunction; J_cf is CF)
        Ln_int = float(Integrate(J_cf * gfApot
                                 * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        tau_us = Ln/Rn*1e6 if Rn > 0 and Ln > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.3f}" if tau_us else "N/A"
        print(f"  [{n}] R={Rn:.4e} (anal {R0_anal:.4e}, ratio "
              f"{Rn/R0_anal:.4f}), L={Ln:.4e}({signL}), tau={tau_str}us "
              f"({elapsed:.1f}s)", flush=True)

        diag.append({"n": n, "R_n": Rn, "L_n": Ln, "tau_us": tau_us,
                     "energy_norm": R_inv, "sign_L": signL,
                     "stage_time_s": elapsed})

        if Ln <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break

        # Update J_cf: J_{n+1} = J_n - sigma * Apot_total / Ln
        # CF tree: subtract sigma_cf * gfApot / Ln from previous J
        J_cf = J_cf - sigma_cf * gfApot * (1.0 / Ln)

    return diag


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari+Kelvin v5: NO J projection (CF chain) ===", flush=True)
    print(f"R_K = {R_K*1000} mm, h_air = {H_AIR*1000} mm, h_cond = {H_COND*1000} mm",
          flush=True)
    print(f"order = {ORDER}, bonus_intorder = {BONUS_INT}", flush=True)
    print(f"Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"ELF Foster B_z tau_lead = 11.51 us\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, mats = {mesh.GetMaterials()} ({time.time()-t0:.1f}s)\n",
          flush=True)

    diag = kameari_kelvin_v5(mesh, N_STAGES)

    print("\n" + "=" * 78, flush=True)
    print("KAMEARI + KELVIN v5 (no projection) SUMMARY", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'R_n':>14} {'L_n [H]':>14} {'tau [us]':>12}", flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] else "N/A"
        print(f"{d['n']:>3} {d['R_n']:>14.4e} {d['L_n']:>14.4e} {tau_s:>12}",
              flush=True)

    if diag and diag[0]['R_n']:
        ratio = diag[0]['R_n'] / R0_anal
        print(f"\n  R_0 / Parseval-analytical = {ratio:.6f} "
              f"(should be 1.0 within machine precision)", flush=True)
    if diag and diag[0]['tau_us']:
        print(f"  tau_0 = {diag[0]['tau_us']:.3f} us "
              f"(ELF tau_lead = 11.51 us, ratio {diag[0]['tau_us']/11.51:.4f})",
              flush=True)

    out = {"R_K_mm": R_K*1000, "OFFSET_mm": [o*1000 for o in OFFSET],
           "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
           "order": ORDER, "ne": mesh.ne, "stages": diag,
           "R0_analytic": R0_anal, "elf_foster_tau_lead_us": 11.51}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v5.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
