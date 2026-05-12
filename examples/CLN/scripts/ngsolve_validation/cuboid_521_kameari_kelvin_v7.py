"""Kameari + Kelvin v7: gauge mass with Kelvin-aware metric ν_0² / nu_cf.

Hypothesis: my v3-v6 used gauge_eps × ν_0 × |A|² (constant ν_0) for gauge fix,
which is TOO WEAK near GND in the Kelvin material where physical mass scales
as ν_0² / nu_kelvin = ν_0 × (R_K/ρ')². This means GND is under-constrained,
and A can grow spuriously there, polluting the iteration's L_n.

v7 fix: use gauge_eps × ν_0² / nu_cf × |A|² for the gauge mass term, which
in the Kelvin material correctly scales as the physical mass via the
inverse Jacobian.

This should constrain A → 0 at GND properly and recover the radiation BC.

Expected: τ_0 should now approach mpmath Cauer Stage 0 ~ 14.4 μs.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v7...", flush=True)

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


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari + Kelvin v7: Kelvin-aware gauge mass ν_0²/nu_cf ===",
          flush=True)
    print(f"Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"mpmath Cauer Stage 0 (target): tau_0 ~ 14.4 us\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

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

    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))

    # Kelvin-aware gauge mass weight: ν_0² / nu_cf
    # = ν_0 in non-Kelvin (where nu_cf = ν_0)
    # = ν_0 * (R_K/ρ')² in Kelvin material (which is the physical-pullback weight)
    gauge_weight = NU_0 * NU_0 / nu_cf

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    a += GAUGE_EPS * gauge_weight * u * v * dx(bonus_intorder=BONUS_INT)

    print("  Assembling+factor (Kelvin-aware gauge mass)...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    # Kameari iteration with CF chain (no projection of J)
    diag = []
    J_cf = sigma_cf * A_ext
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    for n in range(N_STAGES):
        t_stage = time.time()

        f = LinearForm(fes)
        f += J_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm bad ({R_inv:.2e}) STOP", flush=True)
            break
        Rn = 1.0 / R_inv

        gfApot.vec.data += Rn * gfA.vec

        Ln_int = float(Integrate(J_cf * gfApot
                                 * dx("conductor", bonus_intorder=BONUS_INT),
                                 mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        tau_us = Ln/Rn*1e6 if Rn > 0 and Ln > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.3f}" if tau_us else "N/A"
        print(f"  [{n}] R={Rn:.4e}, L={Ln:.4e}({signL}), tau={tau_str}us "
              f"({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": Rn, "L_n": Ln, "tau_us": tau_us,
                     "energy_norm": R_inv, "sign_L": signL})

        if Ln <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break

        # Update J_cf as CF chain (no projection)
        J_cf = J_cf - sigma_cf * gfApot * (1.0 / Ln)

    print("\n" + "=" * 78, flush=True)
    print("KAMEARI + KELVIN v7 (Kelvin-aware gauge) SUMMARY", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'R_n':>14} {'L_n [H]':>14} {'tau [us]':>12}",
          flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] else "N/A"
        print(f"{d['n']:>3} {d['R_n']:>14.4e} {d['L_n']:>14.4e} {tau_s:>12}",
              flush=True)

    if diag and diag[0]['tau_us']:
        ratio = diag[0]['tau_us'] / 14.4
        print(f"\n  tau_0 = {diag[0]['tau_us']:.3f} us "
              f"(target ~14.4, ratio {ratio:.3f}x)", flush=True)

    out = {"R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
           "order": ORDER, "ne": mesh.ne, "stages": diag,
           "R0_analytic": R0_anal,
           "target_tau_0_us": 14.4, "elf_tau_lead_us": 11.51}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v7.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
