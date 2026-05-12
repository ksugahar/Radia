"""Kameari + Kelvin v14: 修正された定式 (2026-05-04 完成)

v11/v12/v13 全て破綻した原因と v14 での修正:

(1) (nu - nu_0) reduced-A 形は Kelvin と非互換
    → 直接形 -nu' curl(A_s) curl(v) on kext only を使用
    詳細: docs/kelvin/KELVIN_TRANSFORMATION.md §7.5

(2) A_phys = (B_0/2)(-y, x, 0) は無限遠で発散 → Convention A pullback
    で 1/rho'^3 特異性
    → Convention B (-(rho'/R)^2 × A_phys) を使用 (offset-local 評価、消滅)
    詳細: make_reduced_potential_background_cf in radia.kelvin_material

(3) Kelvin 領域での A_s 評価は **offset 中心の local 座標** を使う
    → make_reduced_potential_background_cf が自動で対応

期待値: ELF reference tau_lead = 11.51 μs (4-stage Cauer ladder)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v14 (CORRECTED: Convention B + direct form)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, sqrt,
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
N_STAGES = 4   # 4 stages for ELF cross-check (tau_lead from ladder)
BONUS_INT = 8
GAUGE_EPS = 1e-8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))
ELF_TAU_LEAD = 11.51e-6   # seconds


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
    print(f"=== Kameari + Kelvin v14: CORRECTED formulation ===", flush=True)
    print(f"  Convention B (offset-local A_s) + direct form RHS", flush=True)
    print(f"  Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
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

    # === A_s using Convention B (CORRECTED): -(rho'/R)^2 × A_phys at offset-local ===
    # A_phys = (B_0/2)(-y, x, 0) for uniform B_z = 1.
    # Convention B factor: -(rho'/R)^2 in Kelvin region (offset-local eval).
    # Vanishes at offset, no singularity. See docs/kelvin/KELVIN_TRANSFORMATION.md §7.4.
    A_s_cf = make_reduced_potential_background_cf(
        mesh,
        F_inner_factory=lambda xc, yc, zc: CoefficientFunction((-yc, xc, 0)) * 0.5,
        R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3,
    )

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

    # === Project A_s to HCurl GridFunction (also use Periodic to enforce
    #     tangential continuity at the kelvin_int <-> kelvin_ext interface) ===
    gfA_s = GridFunction(fes, name="A_s")
    gfA_s.vec[:] = 0.0
    print("  Projecting A_s onto HCurl (Periodic)...", flush=True)
    t0 = time.time()
    with TaskManager():
        gfA_s.Set(A_s_cf, bonus_intorder=BONUS_INT)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    # Diagnostic: gfA_s in conductor matches A_phys exactly (linear)
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

    # === Kameari iteration with CORRECTED reduced-A ===
    diag = []
    J_imp_cf = sigma_cf * gfA_s
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    for n in range(N_STAGES):
        t_stage = time.time()

        # CORRECTED reduced-A linear form (2026-05-04):
        #   a(A_r, v) = (J_imp, v)_cond - int_kext nu' curl(A_s) curl(v) dV
        # The OLD -(nu - nu_0) form is INVALID with Kelvin pullback.
        # See docs/kelvin/KELVIN_TRANSFORMATION.md §7.5.
        f = LinearForm(fes)
        f += J_imp_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        f += -nu_cf * curl(gfA_s) * curl(v) \
                * dx("kelvin", bonus_intorder=BONUS_INT)
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

        # L_n = R_n * <J_n, A_total>_cond, with A_total = gfA_s + gfApot
        Ln_int = float(Integrate(J_imp_cf * (gfA_s + gfApot)
                                 * dx("conductor", bonus_intorder=BONUS_INT),
                                 mesh))
        L_n = R_n * Ln_int

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        # tau_n = L_n / R_n (Tanimoto convention; see KELVIN_TRANSFORMATION.md §7.4)
        tau_us = L_n/R_n*1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] R={R_n:.4e}, L={L_n:.4e}({signL}), tau={tau_str}us "
              f"({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "tau_us": tau_us,
                     "energy_norm": R_inv, "sign_L": signL})

        if L_n <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break

        # Schmidt orthogonalization
        J_imp_cf = J_imp_cf - sigma_cf * (gfA_s + gfApot) * (1.0 / L_n)

    print("\n" + "=" * 78, flush=True)
    print("KAMEARI + KELVIN v14 (CORRECTED) RESULTS", flush=True)
    print("=" * 78, flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"  n={d['n']}: R={d['R_n']:.4e}, L={d['L_n']:.4e}({d['sign_L']}), "
              f"tau={tau_s} us", flush=True)

    if diag and diag[0]['tau_us']:
        ratio = diag[0]['tau_us'] / (ELF_TAU_LEAD * 1e6)
        print(f"\n  tau_0 = {diag[0]['tau_us']:.4f} us "
              f"(ELF tau_lead = {ELF_TAU_LEAD*1e6:.2f}, ratio {ratio:.3f}x)",
              flush=True)

    # Find dominant tau (max over stages, since tau_lead may not be tau_0)
    valid_taus = [d['tau_us'] for d in diag if d['tau_us'] is not None]
    if valid_taus:
        tau_max = max(valid_taus)
        ratio_max = tau_max / (ELF_TAU_LEAD * 1e6)
        print(f"  max tau_n = {tau_max:.4f} us "
              f"(vs ELF {ELF_TAU_LEAD*1e6:.2f}, ratio {ratio_max:.3f}x)",
              flush=True)

    if len(diag) >= 2:
        all_pos = all(d['L_n'] > 0 for d in diag)
        if all_pos:
            print("  [OK] All stages positive L -- formulation stable",
                  flush=True)
        else:
            print("  [NG] Sign flip detected", flush=True)

    out = {"method": "Kameari + Kelvin v14 (Convention B + direct form)",
           "R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000,
           "h_air_mm": H_AIR*1000, "order": ORDER, "ne": mesh.ne,
           "rel_err_gfA_s_cond": rel_err,
           "stages": diag, "R0_analytic": R0_anal,
           "elf_tau_lead_us": ELF_TAU_LEAD*1e6}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v14_corrected.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
