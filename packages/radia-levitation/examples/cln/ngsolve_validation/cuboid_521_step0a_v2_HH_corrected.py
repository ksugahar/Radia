"""Step 0a v2: A-formulation Kameari Cauer ladder with Tanimoto-canonical
recipe (Helmholtz-Hodge gauge correction, NO GAUGE_EPS, NO tree-cotree).

Geometry: cuboid 5x2x1 mm Cu in vacuum sphere, Dirichlet on outer sphere.
Drive: J_0 = (sigma * B_0/2) (-y, x, 0) for uniform B_z = 1 T applied.

Reference recipe (Tanimoto canonical, 20240917_A_ICCG_最新版.ipynb):
  fes  = HCurl(mesh, order, nograds=True OR type1=True, dirichlet=...)
  gauge= H1(mesh, order, dirichlet=...)
  for nStage:
      R = 1/Integrate(J*J/sigma*dx, mesh)
      Solve a(A,v)=(1/mu)curl(A)curl(v),  f(v)=(J,v)              # NO penalty
      Solve aa(phi,psi)=grad(phi)grad(psi), ff(psi)=-grad(psi)·A  # Helmholtz-Hodge
      A_corr = A + grad(phi)
      Apot += R * A_corr            # CF accumulation
      L = Integrate(R*J*Apot*dx, mesh)
      J = J - sigma * Apot / L

Key fixes vs Step 0a v1 (which had sign-flip at stage 1):
- REMOVE penalty term (GAUGE_EPS * NU_0 * u*v) — Tanimoto: bad
- REMOVE tree-cotree mask (use Helmholtz-Hodge instead)
- ADD H1 gauge solve and grad(phi) correction
- Use Tanimoto-canonical L formula

Compare:
- Closed-PEC limit (R_outer -> 0): tau_lead = 25.46 us
- Vacuum (ELF): tau_lead = 11.51 us
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
print("Starting Step 0a v2 (Helmholtz-Hodge corrected)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, grad, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from math import pi
import json, time, argparse
from pathlib import Path

mu0 = 4 * pi * 1e-7
NU_0 = 1.0 / mu0
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

H_COND = 0.5e-3
H_AIR = 2.0e-3
ORDER = 2
N_STAGES = 8
BONUS_INT = 8

ELF_TAU_LEAD = 11.51e-6
TAU_PEC_LIMIT = mu0 * sigma_Cu * ax**2 * ay**2 / (pi**2 * (ax**2 + ay**2))
R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_geo(R_outer):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_outer = Sphere(Pnt(0, 0, 0), R_outer)
    for f in sphere_outer.faces:
        f.name = "outer"
    sphere_outer.maxh = H_AIR

    air = sphere_outer - cuboid
    air.name = "air"

    geo = Glue([cuboid, air])
    return OCCGeometry(geo)


def main(R_outer, order=ORDER, n_stages=N_STAGES):
    ngsglobals.msg_level = 0
    print(f"=== Step 0a v2: A-form + air-box + Helmholtz-Hodge, R_outer = {R_outer*1000:.1f} mm ===", flush=True)
    print(f"  Closed-PEC analytic tau_lead     = {TAU_PEC_LIMIT*1e6:.3f} us", flush=True)
    print(f"  ELF (vacuum) reference tau_lead  = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  Analytic R_0 (Case B)            = {R0_anal:.6e} Ohm", flush=True)
    print(f"  N_STAGES = {n_stages}, ORDER = {order}\n", flush=True)

    geo = build_geo(R_outer)
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}\n", flush=True)

    # Regularize air with small conductivity to keep operator well-conditioned
    SIGMA_AIR_RATIO = 1e-6   # eps * sigma_Cu in air (Tanimoto-style regularization)
    sigma_air = SIGMA_AIR_RATIO * sigma_Cu
    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": sigma_air})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 1.0/sigma_air})
    print(f"  sigma_air = {sigma_air:.2e} S/m  (regularization {SIGMA_AIR_RATIO:.0e}*sigma_Cu)\n", flush=True)

    # === FE spaces (Tanimoto canonical) ===
    fes = HCurl(mesh, order=order, dirichlet="outer", nograds=True)
    gauge = H1(mesh, order=order, dirichlet="outer")
    print(f"  HCurl ndof = {fes.ndof}", flush=True)
    print(f"  H1    ndof = {gauge.ndof}\n", flush=True)

    u, v = fes.TnT()
    uu, vv = gauge.TnT()

    # === Bilinear forms (NO penalty term) ===
    a_HC = BilinearForm(fes)
    a_HC += NU_0 * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)

    a_HH = BilinearForm(gauge)
    a_HH += grad(uu) * grad(vv) * dx(bonus_intorder=BONUS_INT)

    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_HC.Assemble()
        a_HH.Assemble()
        try:
            inv_HC = a_HC.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_HC = a_HC.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        try:
            inv_HH = a_HH.mat.Inverse(gauge.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_HH = a_HH.mat.Inverse(gauge.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    # === Initial drive: J_0 = sigma_cf * A_s, A_s = (B_0/2)(-y, x, 0) ===
    A_s_cf = CoefficientFunction((-y, x, 0)) * 0.5
    J = sigma_cf * A_s_cf

    R0_check = float(Integrate(J * J * sigma_inv_cf
                               * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    R0_init = 1.0 / R0_check if R0_check > 1e-30 else float('inf')
    print(f"  R_0 init: {R0_init:.6e} Ohm vs analytic {R0_anal:.6e}, "
          f"ratio {R0_init/R0_anal:.6f}\n", flush=True)

    diag = []
    Apot = None   # CF accumulator (Tanimoto canonical)

    for n in range(n_stages):
        t_stage = time.time()

        # Resistance: R_n = 1 / <J, J/sigma>_cond
        R_inv = float(Integrate(J * J * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] R_inv too small ({R_inv:.2e}), STOP", flush=True)
            break
        R_n = 1.0 / R_inv

        # === Solve A_n: curl-curl(A) = mu * J  ===
        f = LinearForm(fes)
        f += J * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv_HC * f.vec

        # === Helmholtz-Hodge correction: solve grad-grad(phi) = -div(A) ===
        # Weak: (grad(uu), grad(vv)) = -(grad(vv), gfA)  -> phi such that
        #       div(A + grad(phi)) = 0 weakly
        ff = LinearForm(gauge)
        ff += grad(vv) * gfA * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            ff.Assemble()
            ff.vec.data = -ff.vec
        gfu = GridFunction(gauge)
        with TaskManager():
            gfu.vec.data = inv_HH * ff.vec

        # === Corrected A and accumulate Apot (Tanimoto canonical) ===
        if Apot is None:
            Apot = R_n * (gfA + grad(gfu))
        else:
            Apot = Apot + R_n * (gfA + grad(gfu))

        # === L_n = Integrate(R*J*Apot*dx, conductor) ===
        # Note: J·A_pot vanishes outside conductor (sigma=0), but A_pot may be
        # nonzero in air. Restrict integration to conductor where J is supported.
        L_n = float(Integrate(R_n * J * Apot
                              * dx("conductor", bonus_intorder=BONUS_INT), mesh))

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        tau_us = L_n / R_n * 1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] R={R_n:.4e}, L={L_n:.4e}({signL}), "
              f"tau=L/R={tau_str}us  ({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "tau_us": tau_us,
                     "energy_norm_inv": R_inv, "sign_L": signL})

        if L_n <= 0:
            print(f"      ** SIGN FLIP at n={n} ** STOP", flush=True)
            break

        # Schmidt update: J = J - sigma * Apot / L
        J = J - sigma_cf * Apot / L_n

    print("\n" + "=" * 78, flush=True)
    print(f"STEP 0a v2 (Helmholtz-Hodge, R_outer={R_outer*1000:.1f} mm) RESULTS", flush=True)
    print("=" * 78, flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"  n={d['n']}: R={d['R_n']:.4e}, L={d['L_n']:.4e}({d['sign_L']}), "
              f"tau={tau_s} us", flush=True)

    valid_taus = [d['tau_us'] for d in diag if d['tau_us'] is not None]
    if valid_taus:
        tau_max = max(valid_taus)
        ratio_pec = tau_max / (TAU_PEC_LIMIT * 1e6)
        ratio_elf = tau_max / (ELF_TAU_LEAD * 1e6)
        print(f"\n  max tau_n = {tau_max:.4f} us", flush=True)
        print(f"    vs closed-PEC limit {TAU_PEC_LIMIT*1e6:.2f} us  ratio {ratio_pec:.3f}", flush=True)
        print(f"    vs ELF (vacuum)     {ELF_TAU_LEAD*1e6:.2f} us  ratio {ratio_elf:.3f}", flush=True)

    out = {"method": "A-form + air-box + Helmholtz-Hodge (Tanimoto canonical)",
           "R_outer_mm": R_outer*1000,
           "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
           "order": order, "ne": mesh.ne, "ndof": fes.ndof,
           "stages": diag, "R0_analytic": R0_anal,
           "tau_pec_limit_us": TAU_PEC_LIMIT*1e6,
           "tau_elf_us": ELF_TAU_LEAD*1e6}
    out_path = Path(__file__).parent / f"cuboid_521_step0a_v2_Router_{int(R_outer*1000):03d}mm.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--R_outer_mm", type=float, default=25.0)
    parser.add_argument("--order", type=int, default=ORDER)
    parser.add_argument("--n_stages", type=int, default=N_STAGES)
    args = parser.parse_args()
    main(R_outer=args.R_outer_mm * 1e-3, order=args.order, n_stages=args.n_stages)
