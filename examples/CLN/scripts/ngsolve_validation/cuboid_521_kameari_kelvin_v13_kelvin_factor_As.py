"""Kameari + Kelvin v13: A_s に Kelvin_factor をかけて全領域で定義 + gfA_s 投影.

ご指導: 「A_s は外部領域でも、Kelvin_factor をかけて出すでいけるはず」

v11 失敗: A_s_cf を analytical CF にして、curl も analytical で作ったが、
         offset 点で 1/rho'^4 の特異性があり (nu - nu_0) curl(A_s) が積分破綻。
v12 失敗: A_s = 0 in Kelvin にしたら reduced-A の引き算項が論理破綻
         (stage 0 で sign flip)。

v13 方針:
  1. A_s は inner + Kelvin の全領域で定義 (make_kelvin_aware_A_s_cf):
     - inner: A_s = A_phys = (B_0/2)(-y, x, 0)
     - Kelvin: A_s_comp = (R/rho')^2 * Householder * A_phys(r_phys)  [pullback]
  2. それを HCurl GridFunction gfA_s に投影 (Set with bonus_intorder).
     -> 投影が特異性を mesh レベルで正則化 (有限値で打ち切られる)
  3. curl は NGSolve 内蔵の curl(gfA_s) を使用 (analytical を回避)。
     -> 離散 curl が双線形形式と整合、A_s の弱い特異性を mesh で吸収。

Reduced-A 弱形式 (A_s 全領域で定義):
    a(A_r, v) = (J_imp, v)_cond
              - integral_full of (nu - nu_0) curl(gfA_s) . curl(v) dV

  - inner: nu = nu_0 -> 第 2 項 = 0 ✓
  - Kelvin: nu = nu_0 (rho'/R)^2 -> 引き算が物理的に正しく効く
  - Kelvin offset 付近: gfA_s 投影で正則化済 -> 特異性消失

2 段試験 (N_STAGES = 2): L_0, L_1 共に正かつ tau 妥当なら成功。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v13 (A_s with Kelvin pullback, gfA_s projection)...",
      flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, sqrt,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import (
    make_kelvin_nu_cf,
    make_kelvin_aware_A_s_cf,
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
N_STAGES = 2
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
    print(f"=== Kameari + Kelvin v13: A_s with Kelvin pullback, 2-stage ===",
          flush=True)
    print(f"  A_s defined on full mesh (inner + Kelvin pullback)", flush=True)
    print(f"  Projected onto HCurl gfA_s, builtin curl(gfA_s)", flush=True)
    print(f"  Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"  Target tau_0 (mpmath) ~ 14.4 us\n", flush=True)

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

    # === A_s analytical CF: Kelvin-aware (pullback in Kelvin region) ===
    def A_phys_factory(xc, yc, zc):
        return CoefficientFunction((-yc, xc, 0)) * 0.5
    A_s_cf = make_kelvin_aware_A_s_cf(
        mesh, A_phys_factory, R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",))

    # === FE space ===
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

    # === gfA_s: project A_s_cf onto HCurl on FULL mesh ===
    # The Kelvin pullback A_s has a 1/rho'^3 singularity at offset center.
    # Set() integrates A_s * basis over each element via quadrature; the
    # singularity is "regularized" at mesh resolution (not sampled at offset).
    # The discrete curl(gfA_s) will be consistent with the bilinear form.
    gfA_s = GridFunction(fes, name="A_s")
    gfA_s.vec[:] = 0.0
    print("  Projecting A_s onto FULL mesh (inner + Kelvin pullback)...",
          flush=True)
    t0 = time.time()
    with TaskManager():
        gfA_s.Set(A_s_cf, bonus_intorder=BONUS_INT)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    # Diagnostic: in conductor, gfA_s should match A_s_phys exactly
    A_s_phys_cf = CoefficientFunction((-y, x, 0)) * 0.5
    err_in_cond = float(Integrate(
        (gfA_s - A_s_phys_cf) * (gfA_s - A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    norm_in_cond = float(Integrate(
        A_s_phys_cf * A_s_phys_cf
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err_cond = (err_in_cond / norm_in_cond) ** 0.5 if norm_in_cond > 0 else 0
    print(f"  ||gfA_s - A_s_phys||/||A_s_phys|| in conductor = "
          f"{rel_err_cond:.2e}", flush=True)

    curl_gfAs = curl(gfA_s)
    curl_target = CoefficientFunction((0, 0, 1))
    err_curl = float(Integrate(
        (curl_gfAs - curl_target) * (curl_gfAs - curl_target)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err_curl = (err_curl / V_cond) ** 0.5
    print(f"  ||curl(gfA_s) - z_hat|| / sqrt(V_cond) in conductor = "
          f"{rel_err_curl:.2e}", flush=True)

    # Diagnostic: ||gfA_s||^2 on Kelvin (singular pullback magnitude)
    normAs2_kelvin = float(Integrate(
        gfA_s * gfA_s * dx("kelvin", bonus_intorder=BONUS_INT), mesh))
    normCurlAs2_kelvin = float(Integrate(
        curl_gfAs * curl_gfAs * dx("kelvin", bonus_intorder=BONUS_INT), mesh))
    print(f"  ||gfA_s||^2 in Kelvin = {normAs2_kelvin:.3e}", flush=True)
    print(f"  ||curl(gfA_s)||^2 in Kelvin = {normCurlAs2_kelvin:.3e}\n",
          flush=True)

    # === Bilinear form ===
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

    # === Verify R_0 ===
    R_inv_check = float(Integrate(sigma_cf * gfA_s * gfA_s
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification: {R0_test:.6e} vs Parseval {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    # === Kameari iteration ===
    diag = []
    J_imp_cf = sigma_cf * gfA_s
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    for n in range(N_STAGES):
        t_stage = time.time()

        # Reduced-A linear form: 2nd term over FULL mesh
        # (nu - nu_0) is nonzero only in Kelvin -> auto-restricted there
        f = LinearForm(fes)
        f += J_imp_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        f += -(nu_cf - NU_0) * curl(gfA_s) * curl(v) \
                * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA_r = GridFunction(fes)
        with TaskManager():
            gfA_r.vec.data = inv * f.vec

        R_inv = float(Integrate(J_imp_cf * J_imp_cf * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT),
                                mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm bad ({R_inv:.2e}) STOP", flush=True)
            break
        R_n = 1.0 / R_inv

        gfApot.vec.data += R_n * gfA_r.vec

        Ln_int = float(Integrate(J_imp_cf * (gfA_s + gfApot)
                                 * dx("conductor", bonus_intorder=BONUS_INT),
                                 mesh))
        L_n = R_n * Ln_int

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        tau_us = L_n/R_n*1e6 if R_n > 0 and L_n > 0 else None
        normAr2 = float(Integrate(gfA_r * gfA_r
                                  * dx(bonus_intorder=BONUS_INT), mesh))
        normAr2_cond = float(Integrate(gfA_r * gfA_r
                                       * dx("conductor", bonus_intorder=BONUS_INT),
                                       mesh))
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] R={R_n:.4e}, L={L_n:.4e}({signL}), tau={tau_str}us",
              flush=True)
        print(f"      ||A_r||^2 full={normAr2:.3e}, in_cond={normAr2_cond:.3e} "
              f"({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "tau_us": tau_us,
                     "energy_norm": R_inv, "sign_L": signL,
                     "normAr2_full": normAr2, "normAr2_cond": normAr2_cond})

        if L_n <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break

        J_imp_cf = J_imp_cf - sigma_cf * (gfA_s + gfApot) * (1.0 / L_n)

    print("\n" + "=" * 78, flush=True)
    print("KAMEARI + KELVIN v13 (A_s Kelvin pullback projected) 2-stage TEST",
          flush=True)
    print("=" * 78, flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"  n={d['n']}: R={d['R_n']:.4e}, L={d['L_n']:.4e}({d['sign_L']}), "
              f"tau={tau_s} us", flush=True)

    if diag and diag[0]['tau_us']:
        print(f"\n  tau_0 = {diag[0]['tau_us']:.4f} us "
              f"(target ~14.4, ratio {diag[0]['tau_us']/14.4:.3f}x)",
              flush=True)
    if len(diag) >= 2:
        all_pos = all(d['L_n'] > 0 for d in diag)
        if all_pos:
            print("  [OK] Both stages positive L -- formulation passes 2-stage",
                  flush=True)
        else:
            print("  [NG] Sign flip -- formulation breakdown", flush=True)

    out = {"method": "Kameari + Kelvin v13 A_s pullback projected 2-stage",
           "R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000,
           "h_air_mm": H_AIR*1000, "order": ORDER, "ne": mesh.ne,
           "rel_err_gfA_s_cond": rel_err_cond,
           "rel_err_curl_gfA_s_cond": rel_err_curl,
           "normAs2_kelvin": normAs2_kelvin,
           "normCurlAs2_kelvin": normCurlAs2_kelvin,
           "stages": diag, "R0_analytic": R0_anal,
           "target_tau_0_us": 14.4, "elf_tau_lead_us": 11.51}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v13_kelvin_factor_As.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
