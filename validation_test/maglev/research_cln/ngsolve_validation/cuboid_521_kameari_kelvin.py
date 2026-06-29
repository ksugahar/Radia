"""Kameari iteration on cuboid 5x2x1 vacuum problem,
with Sugahara two-sphere KELVIN transform as inner solver.

Replaces NGSolve air-box truncation (which gave tau_0 = 104 us, divergent
with AIR_SCALE) with the Kelvin conformal map (radiation BC exact).

Head-to-head with Vector Fitting + ELF (tau_lead = 11.51 us).

Uses the production Kelvin helpers in S:/Radia/01_GitHub/src/radia/.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, Periodic, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import make_kelvin_nu_cf, NU_0
from math import pi
import json
import time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

# Conductor
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Kelvin sphere parameters
R_K = 25e-3              # 25 mm sphere = 5x cuboid bounding sphere (~5x scale)
OFFSET = (4 * R_K, 0, 0) # Outer Kelvin sphere: well-separated from inner

H_COND = 0.5e-3          # 0.5 mm in conductor
H_AIR = 4.0e-3           # 4 mm in air (between cuboid and inner sphere boundary)
H_KELVIN = 8.0e-3        # 8 mm in Kelvin exterior
ORDER = 1                # production setting (Sugahara papers)
N_STAGES = 8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_geo():
    """Cuboid + inner air sphere + Kelvin outer sphere (Sugahara two-sphere)."""
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.name = "air"
    sphere_inner.maxh = H_AIR

    inner_air = sphere_inner - cuboid
    for f in inner_air.faces:
        # Re-tag the spherical face (after subtraction, naming is preserved)
        if f.name not in ("kelvin_int",):
            pass

    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid],
        offset=OFFSET,
        R_K=R_K,
        inner_maxh=H_AIR,
        outer_maxh_factor=H_KELVIN/H_AIR,
    )
    return OCCGeometry(geo), info


def kameari_kelvin(mesh, n_stages):
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))
    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND"))
    print(f"  Periodic HCurl ndof = {fes.ndof}")
    print(f"  Materials = {mesh.GetMaterials()}")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0, "kelvin": 0.0})

    A_ext = CoefficientFunction((-y, x, 0)) * 0.5  # B_0 = 1, z-direction

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
    a += 1e-8 * NU_0 * u * v * dx  # gauge regularization (production setting)

    print("  Assembling Kelvin curl-curl operator...")
    t0 = time.time()
    with TaskManager():
        a.Assemble()
    print(f"  Assembled in {time.time()-t0:.1f}s")

    print("  Factorizing (sparsecholesky fallback if pardiso missing)...")
    t0 = time.time()
    try:
        with TaskManager():
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        print(f"  Pardiso factor: {time.time()-t0:.1f}s")
    except Exception:
        with TaskManager():
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        print(f"  SparseCholesky factor: {time.time()-t0:.1f}s")

    diag = []
    J_cf = sigma_cf * A_ext
    Apot = None

    for n in range(n_stages):
        print(f"\n[Stage {n}]")
        t0 = time.time()
        f = LinearForm(fes)
        f += J_cf * v * dx("conductor")
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx("conductor"), mesh))
        if abs(R_inv) < 1e-30:
            print(f"  R_inv too small, stop")
            break
        Rn = 1.0 / R_inv

        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA

        Ln_int = float(Integrate(J_cf * Apot * dx("conductor"), mesh))
        Ln = Rn * Ln_int
        tau_us = Ln / Rn * 1e6 if Rn > 0 and Ln > 0 else float('nan')
        elapsed = time.time() - t0

        signL = "+" if Ln > 0 else "-" if Ln < 0 else "0"
        signR = "+" if Rn > 0 else "-" if Rn < 0 else "0"
        print(f"  R = {Rn:.4e} ({signR})")
        print(f"  L = {Ln:.4e} ({signL}, {Ln*1e9:.4e} nH)")
        print(f"  tau = {tau_us:.4f} us")
        print(f"  ({elapsed:.1f}s)")

        diag.append({
            "n": n, "R_n": Rn, "L_n": Ln,
            "tau_us": tau_us if not (tau_us != tau_us) else None,
            "energy_norm": R_inv,
            "sign_L": signL, "sign_R": signR,
        })

        if Ln <= 0 or Rn <= 0:
            print(f"  ** SIGN FLIP at stage {n}, breaking **")
            break

        J_cf = J_cf - sigma_cf * Apot / Ln

    return diag


def main():
    ngsglobals.msg_level = 1
    print(f"=== Kameari + Kelvin (production Sugahara two-sphere) ===")
    print(f"  Conductor: {ax*1000}x{ay*1000}x{az*1000} mm Cu")
    print(f"  Inner air sphere: R_K = {R_K*1000} mm")
    print(f"  Outer Kelvin offset: {OFFSET}")
    print(f"  h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, h_kelvin = {H_KELVIN*1000} mm")
    print(f"  order = {ORDER}, N_stages = {N_STAGES}\n")

    geo, info = build_geo()
    print("Generating mesh...")
    t0 = time.time()
    mesh = Mesh(geo.GenerateMesh(maxh=H_KELVIN))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}  ({time.time()-t0:.1f}s)\n")

    print(f"Reference R_0 (Parseval) = {R0_anal:.6e}")
    print(f"ELF Foster B_z tau_lead   = 11.51 us\n")

    diag = kameari_kelvin(mesh, N_STAGES)

    print("\n" + "=" * 78)
    print("KAMEARI + KELVIN SUMMARY")
    print("=" * 78)
    print(f"{'n':>3} {'R_n':>14} {'L_n [nH]':>14} {'tau [us]':>12} {'sign(L)':>8}")
    print("-" * 78)
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] else " N/A"
        print(f"{d['n']:>3} {d['R_n']:>14.4e} {d['L_n']*1e9:>14.4e} {tau_s:>12} "
              f"{d['sign_L']:>8}")

    if diag and diag[0]['tau_us']:
        tau0 = diag[0]['tau_us']
        elf_lead = 11.51
        print(f"\n=== Head-to-head ===")
        print(f"  Kameari + Kelvin     tau_0    = {tau0:.3f} us")
        print(f"  Kameari + air-box    tau_0    = 104.37 us  (5x box, my prior run)")
        print(f"  Vector Fit + ELF     tau_lead = {elf_lead:.3f} us")
        print(f"  Kameari Kelvin / ELF = {tau0/elf_lead:.4f}")
        print(f"  Kameari Kelvin / air-box = {tau0/104.37:.4f}")
        if abs(tau0 - elf_lead) / elf_lead < 0.10:
            print(f"\n  >>> SUCCESS: Kameari + Kelvin recovers ELF tau_lead within 10% <<<")
        elif tau0 < 50:
            print(f"\n  [Kameari + Kelvin gives smaller tau_0 than air-box; better than air-box]")
        else:
            print(f"\n  [Still off from ELF; mesh refinement needed?]")

    out = {
        "method": "Kameari + Kelvin (Sugahara two-sphere) on cuboid 5x2x1 vacuum",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "R_K_mm": R_K*1000,
        "OFFSET_mm": [o*1000 for o in OFFSET],
        "h_cond_mm": H_COND*1000,
        "h_air_mm": H_AIR*1000,
        "h_kelvin_mm": H_KELVIN*1000,
        "order": ORDER,
        "ne": mesh.ne, "nv": mesh.nv,
        "R0_analytic": R0_anal,
        "n_stages_attempted": N_STAGES,
        "n_stages_completed": len(diag),
        "stages": diag,
        "elf_foster_tau_lead_us": 11.51,
    }
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
