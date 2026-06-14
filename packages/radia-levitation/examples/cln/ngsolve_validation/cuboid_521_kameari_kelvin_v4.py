"""Kameari + Kelvin v4: 3 fixes to address GND-vertex gauge contamination.

v3 had τ_0 = 327-668 us (28-58× ELF). Hypotheses:
(F1) gauge_eps × NU_0 (constant) dominates near GND where nu_kelvin → 0.
     Fix: use gauge_eps × nu_cf (modulated) so gauge mass scales with curl-curl.
(F2) gauge_eps = 1e-8 too large; try 1e-12.
(F3) R_K too small (25 mm); try 50 mm to put more headroom.

Tests Stage 0 only — fast, decisive: τ_0 should approach radiation-BC physics.

Comparison target: ELF + VF Foster lead τ_lead = 11.51 us.
But for Kameari Cauer-I, τ_0 != Foster lead in general (different decomp).
A weaker but indicative target: τ_0 should be of ORDER ~10 us (not ~300+).
Or: m_z(0+) slope should match Parseval analytically (already verified for R_0).

After Stage 0 OK, also run N=8 to compare ladder.
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
N_STAGES = 8
ORDER = 2
BONUS_INT = 8

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


def build_geo(R_K, h_air, h_cond=0.5e-3):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = h_cond

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.maxh = h_air

    inner_air = sphere_inner - cuboid
    inner_air.name = "air"

    OFFSET = (2.5 * R_K, 0, 0)
    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid], offset=OFFSET, R_K=R_K, inner_maxh=h_air)
    return OCCGeometry(geo), OFFSET


def kameari_kelvin_v4(mesh, OFFSET, R_K, gauge_eps, modulated_gauge,
                      n_stages=N_STAGES):
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))
    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                         nograds=True))
    print(f"    HCurl ndof = {fes.ndof}", flush=True)

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"    tree-cotree masked {masked} edges", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    if modulated_gauge:
        a += gauge_eps * nu_cf * u * v * dx(bonus_intorder=BONUS_INT)
    else:
        a += gauge_eps * NU_0 * u * v * dx(bonus_intorder=BONUS_INT)

    print(f"    assembling+factor (gauge_eps={gauge_eps}, "
          f"modulated={modulated_gauge})...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"      {time.time()-t0:.1f}s", flush=True)

    diag = []
    gfJ = GridFunction(fes)
    gfJ.Set(sigma_cf * A_ext)
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0
    J_history = []

    for n in range(n_stages):
        f = LinearForm(fes)
        f += gfJ * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(gfJ * gfJ * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"    [{n}] energy_norm bad ({R_inv:.2e}) STOP", flush=True)
            break
        Rn = 1.0 / R_inv

        cross = []
        for gfJ_m in J_history:
            xnm = float(Integrate(gfJ * gfJ_m * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT), mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / R_inv

        gfApot.vec.data += Rn * gfA.vec
        Ln_int = float(Integrate(gfJ * gfApot
                                 * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        tau_us = Ln/Rn*1e6 if Rn > 0 and Ln > 0 else None
        tau_str = f"{tau_us:.3f}" if tau_us else "N/A"
        print(f"    [{n}] R={Rn:.3e}, L={Ln:.3e}({signL}), tau={tau_str}us, "
              f"drift={rel_drift:.2e}, En={R_inv:.3e}", flush=True)
        diag.append({"n": n, "R_n": Rn, "L_n": Ln, "tau_us": tau_us,
                     "rel_drift": rel_drift, "energy_norm": R_inv,
                     "sign_L": signL})

        gfJ_save = GridFunction(fes)
        gfJ_save.vec.data = gfJ.vec.data
        J_history.append(gfJ_save)

        if Ln <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break
        gfJ.vec.data -= (sigma_Cu / Ln) * gfApot.vec

    return diag


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari+Kelvin v4: 3 setup variations, Stage 0..{N_STAGES} ===",
          flush=True)
    print(f"Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"ELF Foster B_z tau_lead = 11.51 us", flush=True)
    print(f"order={ORDER}, bonus_intorder={BONUS_INT}\n", flush=True)

    configs = [
        # (name, R_K, h_air, gauge_eps, modulated_gauge)
        ("v3 baseline (control)",       25e-3, 2.0e-3, 1e-8,  False),
        ("F1 modulated_gauge=True",     25e-3, 2.0e-3, 1e-8,  True),
        ("F2 gauge_eps=1e-12",          25e-3, 2.0e-3, 1e-12, False),
        ("F3 R_K=50mm",                 50e-3, 4.0e-3, 1e-8,  False),
        ("F1+F2 modulated + 1e-12",     25e-3, 2.0e-3, 1e-12, True),
    ]

    all_results = []
    for name, R_K, h_air, gauge_eps, modulated in configs:
        print("=" * 78, flush=True)
        print(name, flush=True)
        print("=" * 78, flush=True)
        try:
            geo, OFFSET = build_geo(R_K, h_air)
            mesh = Mesh(geo.GenerateMesh(maxh=h_air, grading=0.5))
            mesh.Curve(ORDER + 1)
            print(f"    ne = {mesh.ne}, R_K = {R_K*1000} mm, "
                  f"h_air = {h_air*1000} mm", flush=True)
            diag = kameari_kelvin_v4(mesh, OFFSET, R_K, gauge_eps, modulated)
            all_results.append({"config": name, "R_K_mm": R_K*1000,
                                "h_air_mm": h_air*1000,
                                "gauge_eps": gauge_eps,
                                "modulated_gauge": modulated,
                                "ne": mesh.ne, "stages": diag})
        except Exception as e:
            print(f"    FAILED: {e}", flush=True)
        print("", flush=True)

    print("=" * 78, flush=True)
    print("SUMMARY: which fix recovers physical tau_0?", flush=True)
    print("=" * 78, flush=True)
    print(f"{'config':>40} {'tau_0 [us]':>14} {'L_0 [H]':>14}", flush=True)
    for r in all_results:
        if r['stages'] and r['stages'][0]['tau_us']:
            tau0 = r['stages'][0]['tau_us']
            L0 = r['stages'][0]['L_n']
        else:
            tau0 = float('nan'); L0 = float('nan')
        print(f"{r['config']:>40} {tau0:>14.3f} {L0:>14.3e}", flush=True)

    out = {"reference_R0_Parseval": R0_anal,
           "elf_foster_tau_lead_us": 11.51,
           "configurations": all_results}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v4.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
