"""Step 0a v3: A-form Kameari with Tanimoto canonical L = Integrate(B^2/mu dx)
over FULL mesh (magnetic energy).

Difference from v2:
- v2 used L = Integrate(R*J*Apot*dx, conductor) — works for closed PEC,
  fails for air-box because of conductor-surface boundary term.
- v3 uses L = Integrate(B*B/mu*dx, full_mesh) — physical magnetic energy,
  always positive, valid for any geometry.

Reference: CLN_AT.ipynb (Tanimoto thesis, 修論).

In closed PEC: int_par_part B^2 = J*A so v2's formula = v3's formula.
In air-box: gauss-bound integration introduces conductor-surface boundary
term: int_full B^2/mu = int_cond J*A + boundary_term. The boundary term
matters and is what v3 captures.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
print("Starting Step 0a v3 (Tanimoto B^2/mu energy)...", flush=True)

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
    print(f"=== Step 0a v3: A-form + air-box + B^2 energy, R_outer = {R_outer*1000:.1f} mm ===", flush=True)
    print(f"  Closed-PEC analytic tau_lead     = {TAU_PEC_LIMIT*1e6:.3f} us", flush=True)
    print(f"  ELF (vacuum) reference tau_lead  = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  N_STAGES = {n_stages}, ORDER = {order}\n", flush=True)

    geo = build_geo(R_outer)
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})

    fes = HCurl(mesh, order=order, dirichlet="outer", nograds=True)
    gauge = H1(mesh, order=order, dirichlet="outer")
    print(f"  HCurl ndof = {fes.ndof}, H1 ndof = {gauge.ndof}\n", flush=True)

    u, v = fes.TnT()
    uu, vv = gauge.TnT()

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

    A_s_cf = CoefficientFunction((-y, x, 0)) * 0.5
    J = sigma_cf * A_s_cf

    R0_check = float(Integrate(J * J * sigma_inv_cf
                               * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    print(f"  R_0 init: {1.0/R0_check:.6e} Ohm vs analytic {R0_anal:.6e}\n", flush=True)

    diag = []
    Apot = None   # CF accumulator for A
    B_pot = None  # CF accumulator for B = curl(A)
    J_history = []  # GridFunctions for orthogonality diagnostic

    for n in range(n_stages):
        t_stage = time.time()

        R_inv = float(Integrate(J * J * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] R_inv too small ({R_inv:.2e}), STOP", flush=True)
            break
        R_n = 1.0 / R_inv

        # Solve A_n
        f = LinearForm(fes)
        f += J * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv_HC * f.vec

        # Helmholtz-Hodge correction
        ff = LinearForm(gauge)
        ff += grad(vv) * gfA * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            ff.Assemble()
            ff.vec.data = -ff.vec
        gfu = GridFunction(gauge)
        with TaskManager():
            gfu.vec.data = inv_HH * ff.vec

        # Accumulate Apot and B_pot (Tanimoto canonical, with R_n weighting)
        # Note: gauge correction phi affects A but not curl(A)=B
        if Apot is None:
            Apot = R_n * (gfA + grad(gfu))
            B_pot = R_n * curl(gfA)
        else:
            Apot = Apot + R_n * (gfA + grad(gfu))
            B_pot = B_pot + R_n * curl(gfA)

        # === L_n = magnetic energy = Integrate(B^2/mu, full_mesh) ===
        L_n = float(Integrate(B_pot * B_pot / mu0 * dx(bonus_intorder=BONUS_INT), mesh))

        # Cross-check: L via J*Apot integration over conductor
        L_n_alt = float(Integrate(R_n * J * Apot
                                  * dx("conductor", bonus_intorder=BONUS_INT), mesh))

        # Orthogonality drift
        gfJ_curr = GridFunction(fes)
        gfJ_curr.Set(J, definedon=mesh.Materials("conductor"), bonus_intorder=BONUS_INT)
        max_drift = 0.0
        for m, gfJ_m in enumerate(J_history):
            cross = float(Integrate(gfJ_curr * gfJ_m * sigma_inv_cf
                                    * dx("conductor", bonus_intorder=BONUS_INT), mesh))
            drift_rel = abs(cross) / R_inv if R_inv > 0 else 0
            max_drift = max(max_drift, drift_rel)
        J_history.append(gfJ_curr)

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        tau_us = L_n / R_n * 1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] R={R_n:.4e}, L_B2={L_n:.4e}({signL}), L_alt={L_n_alt:.4e}, "
              f"tau={tau_str}us, drift={max_drift:.2e}  ({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n_B2": L_n, "L_n_alt": L_n_alt,
                     "tau_us": tau_us, "drift": max_drift, "sign_L": signL})

        if L_n <= 0:
            print(f"      ** L_B2 <= 0 ** STOP", flush=True)
            break

        J = J - sigma_cf * Apot / L_n

    print("\n" + "=" * 78, flush=True)
    print(f"STEP 0a v3 (B^2 energy, R_outer={R_outer*1000:.1f} mm) RESULTS", flush=True)
    print("=" * 78, flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"  n={d['n']}: R={d['R_n']:.4e}, L_B2={d['L_n_B2']:.4e}({d['sign_L']}), "
              f"L_alt={d['L_n_alt']:.4e}, tau={tau_s} us, drift={d['drift']:.2e}",
              flush=True)

    valid_taus = [d['tau_us'] for d in diag if d['tau_us'] is not None]
    if valid_taus:
        tau_max = max(valid_taus)
        print(f"\n  max tau_n = {tau_max:.4f} us "
              f"(vs ELF {ELF_TAU_LEAD*1e6:.2f}, ratio {tau_max/(ELF_TAU_LEAD*1e6):.3f})",
              flush=True)

    out = {"method": "A-form + air-box + Helmholtz-Hodge + Tanimoto B^2/mu energy",
           "R_outer_mm": R_outer*1000,
           "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
           "order": order, "ne": mesh.ne, "ndof": fes.ndof,
           "stages": diag,
           "tau_pec_limit_us": TAU_PEC_LIMIT*1e6, "tau_elf_us": ELF_TAU_LEAD*1e6}
    out_path = Path(__file__).parent / f"cuboid_521_step0a_v3_Router_{int(R_outer*1000):03d}mm.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--R_outer_mm", type=float, default=25.0)
    parser.add_argument("--order", type=int, default=ORDER)
    parser.add_argument("--n_stages", type=int, default=N_STAGES)
    args = parser.parse_args()
    main(R_outer=args.R_outer_mm * 1e-3, order=args.order, n_stages=args.n_stages)
