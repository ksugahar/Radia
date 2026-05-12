"""Kameari + Kelvin v11: Reduced-A formulation (proper energy subtraction).

Per ご指導:
  「一様磁場のエネルギーは ∞ だが、それを差し引いて、CLN を構築する」
  「2段試験をしたら定式化が破綻しているかどうかすぐわかる」

Reduced-A formulation (S:/Radia/01_GitHub/examples/kelvin_transformation/
A-formulation/docs/reduced_A_formulation.md):

  A = A_s + A_r
  - A_s: source (uniform B applied), KNOWN, energy ∞ (excluded from FE)
  - A_r: perturbation, decays at infinity, FE space target (Kelvin OK)

  Weak form for A_r:
    a(A_r, v) = (J_imp, v)_cond - ∫_Ω (ν - ν_0) curl(A_s) · curl(v) dV

  - First term: standard impressed J source on conductor
  - Second term: "subtracts the infinite uniform-B energy" via Kelvin
    material's ν modulation (nonzero only where ν ≠ ν_0, = Kelvin region)

Kameari iteration:
  Stage 0: J_0 = σ A_s in conductor → solve for A_r_0 (reduced).
  L_0 = R_0 × ⟨J_0, A_total_0⟩_cond where A_total = A_s + A_r_0 in conductor

For Kelvin pullback of A_s:
  A_s_phys = (B_0/2)(-y, x, 0) at physical point r_phys
  A_s_kelvin (computational) = (R_K/ρ')² × Householder × A_s_phys(R_K²/ρ'/...)
  Implemented via make_kelvin_aware_A_s_cf

2-stage test: if L_0, L_1 both positive and τ_n reasonable, formulation OK.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v11 (reduced-A formulation)...", flush=True)

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
N_STAGES = 2   # 2-stage test per ご指導
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
    print(f"=== Kameari + Kelvin v11: REDUCED-A formulation, 2-stage test ===",
          flush=True)
    print(f"  Per ご指導: 一様磁場の ∞ エネルギーを差し引く", flush=True)
    print(f"  Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"  Target tau_0 (mpmath) ~ 14.4 us\n", flush=True)

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

    # === A_s: Kelvin-aware source for uniform B_0 z_hat ===
    # In non-Kelvin material: A_s = (B_0/2)(-y, x, 0) directly
    # In Kelvin material: A_s = (R_K/ρ')² × Householder × A_s_phys(r_phys)
    def A_phys_factory(xc, yc, zc):
        return CoefficientFunction((-yc, xc, 0)) * 0.5

    A_s_cf = make_kelvin_aware_A_s_cf(
        mesh, A_phys_factory, R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",))

    # === curl(A_s) as analytical CF (since curl() not overloaded for CF) ===
    # In non-Kelvin material: curl A_s = curl((B_0/2)(-y, x, 0)) = B_0 ẑ = (0, 0, 1)
    # In Kelvin material: B_kelvin = -(R_K/ρ')⁴ × Householder × B_phys
    # For uniform B_phys = (0, 0, 1), Householder × ẑ = ẑ - 2(ẑ·n̂)n̂
    ox, oy, oz = OFFSET
    dxp = x - ox
    dyp = y - oy
    dzp = z - oz
    rho_p_sq = dxp**2 + dyp**2 + dzp**2 + 1e-24
    rho_p = sqrt(rho_p_sq)
    nx = dxp / rho_p
    ny = dyp / rho_p
    nz = dzp / rho_p
    # Householder of ẑ = (0,0,1):  H ẑ = ẑ - 2(ẑ·n̂)n̂ = (0-2 nz nx, 0-2 nz ny, 1-2 nz²)
    Hz_x = -2 * nz * nx
    Hz_y = -2 * nz * ny
    Hz_z = 1 - 2 * nz * nz
    # Kelvin B factor: -(R_K/ρ')⁴
    factor_kelvin = -(R_K / rho_p) ** 4
    curl_A_s_kelvin = CoefficientFunction((factor_kelvin * Hz_x,
                                            factor_kelvin * Hz_y,
                                            factor_kelvin * Hz_z))
    curl_A_s_inner = CoefficientFunction((0.0, 0.0, 1.0))  # uniform B_0 = 1

    def _switch(kelvin_comp, inner_comp):
        d = {}
        for m in mesh.GetMaterials():
            ml = m.lower()
            is_kelvin = "kelvin" in ml
            d[m] = kelvin_comp if is_kelvin else inner_comp
        return mesh.MaterialCF(d, default=inner_comp)

    curl_A_s_cf = CoefficientFunction((
        _switch(curl_A_s_kelvin[0], curl_A_s_inner[0]),
        _switch(curl_A_s_kelvin[1], curl_A_s_inner[1]),
        _switch(curl_A_s_kelvin[2], curl_A_s_inner[2]),
    ))

    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))

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

    # === Bilinear form (same as full-A): nu_kelvin curl(u) · curl(v) ===
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

    # === Verify R_0 (impressed J^2/sigma in conductor) ===
    # In conductor: A_s = (B_0/2)(-y, x, 0) (no Kelvin pullback in non-Kelvin)
    # J_0 = sigma * A_s
    R_inv_check = float(Integrate(sigma_cf * A_s_cf * A_s_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification: {R0_test:.6e} vs Parseval {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    # === Kameari iteration with REDUCED-A solves ===
    diag = []
    # J_imp_cf accumulates the iteration's J source (only modified in conductor)
    # Stage 0: J_imp = sigma * A_s_cf (in conductor, A_s_cf = A_s analytical)
    J_imp_cf = sigma_cf * A_s_cf
    # Apot accumulates A_r contributions (NOT A_s — A_s is fixed background)
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    for n in range(N_STAGES):
        t_stage = time.time()

        # === Reduced-A linear form ===
        # f(v) = (J_imp, v)_cond - ((nu - nu_0) curl(A_s), curl(v))_full
        f = LinearForm(fes)
        f += J_imp_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        # The Kelvin energy-subtraction term (only nonzero in Kelvin material)
        f += -(nu_cf - NU_0) * curl_A_s_cf * curl(v) * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        # Solve for A_r (perturbation; A_total = A_s + A_r in conductor)
        gfA_r = GridFunction(fes)
        with TaskManager():
            gfA_r.vec.data = inv * f.vec

        # R_n = 1 / <J_n, J_n/sigma>_cond  (impressed J energy norm)
        R_inv = float(Integrate(J_imp_cf * J_imp_cf * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm bad ({R_inv:.2e}) STOP", flush=True)
            break
        R_n = 1.0 / R_inv

        # Apot accumulates: gfApot += R_n * gfA_r
        # Note: only A_r accumulates (A_s is fixed background)
        gfApot.vec.data += R_n * gfA_r.vec

        # L_n = R_n × <J_n, A_total>_cond where A_total = A_s + Apot in conductor
        # In conductor: A_s_cf = A_s (no Kelvin pullback in non-Kelvin region)
        Ln_int = float(Integrate(J_imp_cf * (A_s_cf + gfApot)
                                 * dx("conductor", bonus_intorder=BONUS_INT),
                                 mesh))
        L_n = R_n * Ln_int

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        tau_us = L_n/R_n*1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] R={R_n:.4e}, L={L_n:.4e}({signL}), tau={tau_str}us "
              f"({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "tau_us": tau_us,
                     "energy_norm": R_inv, "sign_L": signL})

        if L_n <= 0:
            print(f"      ** SIGN FLIP at n={n} (formulation may be broken) **",
                  flush=True)
            break

        # Update J for next stage: J_{n+1} = J_n - sigma * (A_s + Apot) / L_n
        # Equivalent Schmidt orthogonalization on full A_total
        J_imp_cf = J_imp_cf - sigma_cf * (A_s_cf + gfApot) * (1.0 / L_n)

    print("\n" + "=" * 78, flush=True)
    print("KAMEARI + KELVIN v11 (REDUCED-A) 2-stage TEST", flush=True)
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
            print("  ✓ Both stages have positive L — formulation appears OK",
                  flush=True)
        else:
            print("  ✗ Sign flip detected — formulation breakdown", flush=True)

    out = {"method": "Kameari + Kelvin v11 reduced-A 2-stage",
           "R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000,
           "h_air_mm": H_AIR*1000, "order": ORDER, "ne": mesh.ne,
           "stages": diag, "R0_analytic": R0_anal,
           "target_tau_0_us": 14.4, "elf_tau_lead_us": 11.51}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v11_reduced_A.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
