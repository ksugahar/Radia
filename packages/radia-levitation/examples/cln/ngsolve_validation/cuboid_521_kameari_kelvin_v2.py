"""Kameari + Kelvin v2: use solve_full_A_kelvin directly (production solver)
   with mesh.Curve(2), order=2, grading=0.5 — mirroring the canonical example
   Coil_3D_A_HCurl_with_Kelvin.py.

Key fixes vs v1:
  - mesh.Curve(2) for spherical face accuracy
  - order = 2 (1 was too low for curl-curl + conformal map)
  - grading = 0.5 in mesh generation
  - call solve_full_A_kelvin (= production code path) for each Kameari stage
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, CoefficientFunction, x, y, z, dx,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_solver import solve_full_A_kelvin, NU_0
from math import pi
import json
import time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

R_K = 25e-3                     # 25 mm
OFFSET = (2.5 * R_K, 0, 0)       # = (62.5 mm, 0, 0); offset/R_K = 2.5 like the example

H_COND = 0.5e-3
H_AIR = 4.0e-3
ORDER = 2
N_STAGES = 8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


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
        [inner_air, cuboid],
        offset=OFFSET,
        R_K=R_K,
        inner_maxh=H_AIR,
    )
    return OCCGeometry(geo)


def kameari_stage_solve(mesh, J_cf):
    """Solve curl-curl A = J on the Kelvin geometry, return GridFunction."""
    res = solve_full_A_kelvin(
        mesh,
        J_source_cf=J_cf,
        R_K=R_K,
        offset=OFFSET,
        source_material="conductor",   # cuboid material name
        order=ORDER,
    )
    return res["gfu"], res["nu_cf"]


def main():
    ngsglobals.msg_level = 1
    print(f"=== Kameari + Kelvin v2 (production solver path) ===")
    print(f"  Conductor: {ax*1000}x{ay*1000}x{az*1000} mm Cu")
    print(f"  Kelvin: R_K = {R_K*1000} mm, offset = {[o*1000 for o in OFFSET]} mm")
    print(f"  Mesh: h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, order = {ORDER}")
    print(f"  N_stages = {N_STAGES}\n")

    geo = build_geo()
    print("Generating mesh (grading=0.5)...")
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}, "
          f"mats = {mesh.GetMaterials()}  ({time.time()-t0:.1f}s)\n")

    print(f"Reference R_0 (Parseval) = {R0_anal:.6e}")
    print(f"ELF Foster B_z tau_lead   = 11.51 us\n")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0, "kelvin": 0.0})

    A_ext = CoefficientFunction((-y, x, 0)) * 0.5  # B_0 = 1, z-direction

    diag = []
    J_cf = sigma_cf * A_ext
    Apot = None

    for n in range(N_STAGES):
        print(f"\n[Stage {n}] solve_full_A_kelvin...")
        t0 = time.time()
        gfA, nu_cf = kameari_stage_solve(mesh, J_cf)
        t_solve = time.time() - t0

        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx("conductor"), mesh))
        if abs(R_inv) < 1e-30:
            print(f"  R_inv too small ({R_inv:.2e}), stop")
            break
        Rn = 1.0 / R_inv

        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA

        Ln_int = float(Integrate(J_cf * Apot * dx("conductor"), mesh))
        Ln = Rn * Ln_int

        tau_us = Ln/Rn * 1e6 if Rn > 0 and Ln > 0 else float('nan')

        signL = "+" if Ln > 0 else "-" if Ln < 0 else "0"
        signR = "+" if Rn > 0 else "-" if Rn < 0 else "0"
        print(f"  R = {Rn:.4e} ({signR})")
        print(f"  L = {Ln:.4e} ({signL}, {Ln*1e9:.4e} nH)")
        print(f"  tau = {tau_us:.4f} us")
        print(f"  ({t_solve:.1f}s)")

        diag.append({
            "n": n, "R_n": Rn, "L_n": Ln,
            "tau_us": tau_us if not (tau_us != tau_us) else None,
            "energy_norm": R_inv,
            "sign_L": signL, "sign_R": signR,
            "solve_time_s": t_solve,
        })

        if Ln <= 0 or Rn <= 0:
            print(f"  ** SIGN FLIP at stage {n} **")
            if n >= 2:
                break

        J_cf = J_cf - sigma_cf * Apot / Ln

    print("\n" + "=" * 78)
    print("KAMEARI + KELVIN v2 SUMMARY")
    print("=" * 78)
    print(f"{'n':>3} {'R_n':>14} {'L_n [H]':>14} {'tau [us]':>12} {'sign(L)':>8}")
    print("-" * 78)
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] else " N/A"
        print(f"{d['n']:>3} {d['R_n']:>14.4e} {d['L_n']:>14.4e} {tau_s:>12} "
              f"{d['sign_L']:>8}")

    if diag and diag[0]['tau_us']:
        tau0 = diag[0]['tau_us']
        elf_lead = 11.51
        print(f"\n=== Head-to-head ===")
        print(f"  Kameari + Kelvin v2  tau_0    = {tau0:.3f} us")
        print(f"  Kameari + air-box    tau_0    = 104.37 us  (5x box)")
        print(f"  Vector Fit + ELF     tau_lead = {elf_lead:.3f} us")
        print(f"  Kameari Kelvin / ELF = {tau0/elf_lead:.4f}")
        if abs(tau0 - elf_lead) / elf_lead < 0.10:
            print(f"\n  >>> Kameari + Kelvin recovers ELF tau_lead within 10% <<<")

    out = {
        "method": "Kameari + Kelvin v2 (production solver)",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "R_K_mm": R_K*1000,
        "OFFSET_mm": [o*1000 for o in OFFSET],
        "h_cond_mm": H_COND*1000,
        "h_air_mm": H_AIR*1000,
        "order": ORDER,
        "ne": mesh.ne, "nv": mesh.nv,
        "n_stages_attempted": N_STAGES,
        "n_stages_completed": len(diag),
        "stages": diag,
        "elf_foster_tau_lead_us": 11.51,
    }
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v2.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
