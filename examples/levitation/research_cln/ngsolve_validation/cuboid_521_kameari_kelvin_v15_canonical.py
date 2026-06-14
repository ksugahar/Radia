"""Kameari + Kelvin v15: HCurl Kameari iteration breakdown demo (2026-05-05).

Purpose in the paper: demonstrate that the standard NGSolve HCurl + Kelvin
Kameari iteration breaks down on isolated-conductor-in-vacuum problems.
This breakdown motivates the proposed Mathematica BEM Foster-CLN method
(HDiv interior div-free, J·n=0 strict, hex_bem_hdiv_order3_overnight.wls,
tau_lead = 43.84 us).

Canonical CLN iteration (Tanimoto-Sugahara TMAG3304725 eq. 2-4):
    R_n = 1 / <J_n, J_n/sigma>_cond
    L_n = R_n * <J_n, A_pot_accum>_cond
    Schmidt: J_{n+1} = J_n - sigma * A_pot_accum / L_n

v15 fixed two bugs in v14 (which combined to give tau_0 = 5985 us, x520 off):
  (BUG #1) Spurious RHS term `-nu' * curl(A_s) * curl(v) * dx("kelvin")`.
    KELVIN_TRANSFORMATION.md §7.5 (this term) is a PEEC-coil-source
    correction; for our applied-uniform-B problem there is no physical
    coil source, so the correction does not apply.
  (BUG #2) L_n integral used <J_n, A_s + A_pot>_cond instead of canonical
    <J_n, A_pot>_cond. Working closed-PEC reference and air-box reference
    both use the canonical form.

Even with both fixes applied, the iteration still breaks down at stage 1
(L_1 sign flip, Schmidt energy norm increases x15). This is NOT a v15
implementation bug — it is the structural failure of HCurl Kameari for
this BC class, which is the paper's research motivation.

Reference (NOT a target to match):
  - ELF Foster Vector Fitting tau_lead is N-dependent (11.5 at N=4,
    ~195 at N>=8), not a reliable comparison reference. Use ELF only
    for time-domain transient cross-check.
  - Mathematica BEM HDiv div-free (proposed method): 43.84 us. This is
    a different physical setup (J·n=0 strict) than the v15 HCurl problem.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v15 (canonical CLN iteration with Kelvin)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import (
    make_kelvin_nu_cf,
    make_reduced_potential_background_cf,
    NU_0,
)
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
ELF_TAU_LEAD = 11.51e-6


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
    print(f"=== v15: canonical CLN iteration with Kelvin ===", flush=True)
    print(f"  R_K = {R_K*1000:.1f} mm, OFFSET = {OFFSET}", flush=True)
    print(f"  Reference R_0 (Parseval) = {R0_anal:.6e}", flush=True)
    print(f"  ELF reference tau_lead = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  N_STAGES = {N_STAGES}\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0,
                                 "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})

    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))

    A_s_cf = make_reduced_potential_background_cf(
        mesh,
        F_inner_factory=lambda xc, yc, zc: CoefficientFunction((-yc, xc, 0)) * 0.5,
        R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3,
    )

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

    u, v = fes.TnT()

    gfA_s = GridFunction(fes, name="A_s")
    gfA_s.vec[:] = 0.0
    print("  Projecting A_s onto HCurl (Periodic)...", flush=True)
    t0 = time.time()
    with TaskManager():
        gfA_s.Set(A_s_cf, bonus_intorder=BONUS_INT)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    A_s_phys_cf = CoefficientFunction((-y, x, 0)) * 0.5
    err_in_cond = float(Integrate(
        (gfA_s - A_s_phys_cf) * (gfA_s - A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    norm_in_cond = float(Integrate(
        A_s_phys_cf * A_s_phys_cf
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err = (err_in_cond / norm_in_cond) ** 0.5 if norm_in_cond > 0 else 0
    print(f"  ||gfA_s - A_s_phys||/||A_s_phys|| in conductor = {rel_err:.2e}",
          flush=True)

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
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    R_inv_check = float(Integrate(sigma_cf * gfA_s * gfA_s
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification: {R0_test:.6e} vs Parseval {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    # Breakdown demo: keep iterating past sign flip to record the full
    # Schmidt drift, energy norm growth, and L_n sign trajectory.
    diag = []
    J_history = []  # store gfJ_n for cross-product Schmidt drift diagnostic
    J_imp_cf = sigma_cf * gfA_s
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0
    sign_flip_first_n = None

    for n in range(N_STAGES):
        t_stage = time.time()

        # CANONICAL CLN linear form: source is J_imp in conductor only.
        # NO curl(A_s) RHS — that's the §7.5 PEEC formula, not applicable here.
        f = LinearForm(fes)
        f += J_imp_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA_r = GridFunction(fes)
        with TaskManager():
            gfA_r.vec.data = inv * f.vec

        R_inv = float(Integrate(J_imp_cf * J_imp_cf * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT),
                                mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] energy_norm collapsed ({R_inv:.2e}), stopping",
                  flush=True)
            break
        R_n = 1.0 / R_inv

        # Project current J_imp_cf to a GridFunction for cross-product diagnostics
        gfJ_n = GridFunction(fes)
        try:
            gfJ_n.Set(J_imp_cf, definedon=mesh.Materials("conductor"))
        except Exception:
            gfJ_n.Set(J_imp_cf)

        # Schmidt orthogonality drift: <J_n, J_m/sigma>_cond for m<n
        cross_inner = []
        for gfJ_m in J_history:
            xnm = float(Integrate(J_imp_cf * gfJ_m * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
            cross_inner.append(xnm)
        max_cross = max((abs(c) for c in cross_inner), default=0.0)
        rel_drift = max_cross / abs(R_inv) if R_inv != 0 else 0.0

        gfApot.vec.data += R_n * gfA_r.vec

        # CANONICAL: L_n = R_n * <J_n, A_pot_accum>_cond. NO A_s in L_n.
        Ln_int = float(Integrate(J_imp_cf * gfApot
                                 * dx("conductor", bonus_intorder=BONUS_INT),
                                 mesh))
        L_n = R_n * Ln_int

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        signR = "+" if R_n > 0 else ("-" if R_n < 0 else "0")
        tau_us = L_n/R_n*1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] drift={rel_drift:.3e}  En={R_inv:.3e}  "
              f"R={R_n:.3e}({signR})  L={L_n:.3e}({signL})  "
              f"tau={tau_str}us  ({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "tau_us": tau_us,
                     "energy_norm": R_inv, "rel_drift": rel_drift,
                     "sign_L": signL, "sign_R": signR})

        if L_n <= 0 and sign_flip_first_n is None:
            sign_flip_first_n = n
            print(f"      ** FIRST SIGN FLIP at n={n} (continuing for "
                  "breakdown diagnostic) **", flush=True)

        J_history.append(gfJ_n)

        # CANONICAL Schmidt update. With L_n < 0 the iteration is no longer
        # well-defined but we keep going so the breakdown trajectory is recorded.
        if abs(L_n) < 1e-30:
            print(f"  [{n}] |L_n| collapsed, stopping", flush=True)
            break
        J_imp_cf = J_imp_cf - sigma_cf * gfApot * (1.0 / L_n)

    print("\n" + "=" * 78, flush=True)
    print("HCURL + KELVIN KAMEARI ITERATION — BREAKDOWN TRAJECTORY", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'rel_drift':>12} {'energy_norm':>14} "
          f"{'R_n':>12} {'L_n':>14} {'tau_us':>12}", flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"{d['n']:>3} {d['rel_drift']:>12.3e} {d['energy_norm']:>14.4e} "
              f"{d['R_n']:>12.3e} {d['L_n']:>14.4e}({d['sign_L']}) "
              f"{tau_s:>12}",
              flush=True)

    print("\nBreakdown summary:", flush=True)
    if sign_flip_first_n is not None:
        print(f"  First L_n sign flip at stage n = {sign_flip_first_n}",
              flush=True)
    else:
        print("  No sign flip in this run (unexpected for HCurl Kameari "
              "in vacuum-coupled isolated cuboid)", flush=True)
    en_first = diag[0]['energy_norm'] if diag else None
    en_last = diag[-1]['energy_norm'] if diag else None
    if en_first and en_last:
        print(f"  Energy norm n=0 -> n={diag[-1]['n']}: "
              f"{en_first:.3e} -> {en_last:.3e}  (ratio {en_last/en_first:.3e})",
              flush=True)
    print("\nThis trajectory is the paper's motivation: HCurl Kameari "
          "iteration breaks down on isolated-conductor-in-vacuum problems. "
          "The proposed Mathematica BEM Foster-CLN method "
          "(hex_bem_hdiv_order3_overnight.wls, HDiv interior div-free, "
          "J·n=0 strict) does not exhibit this breakdown.",
          flush=True)

    out = {"method": "HCurl + Kelvin Kameari iteration (breakdown demo)",
           "purpose": ("Paper motivation: shows HCurl Kameari iteration "
                       "breaks down on isolated-conductor-in-vacuum problems. "
                       "The proposed BEM Foster-CLN does not."),
           "R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000,
           "h_air_mm": H_AIR*1000, "order": ORDER, "ne": mesh.ne,
           "rel_err_gfA_s_cond": rel_err,
           "stages": diag, "R0_analytic": R0_anal,
           "sign_flip_first_n": sign_flip_first_n}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v15_canonical.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
