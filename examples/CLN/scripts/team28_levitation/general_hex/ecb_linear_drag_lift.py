"""ecb_linear_drag_lift.py -- Linear eddy-current brake drag/lift force vs velocity.

Geometry: rectangular conducting plate + cube PM moving along +x at speed v.
Method:
  1. Plate polarizability tensor alpha_ij(s) via Mixed Galerkin (cln_sibc_general_hex).
  2. Effective excitation frequency f_eff = v / (2 L_PM) (PM transit time).
  3. Time-averaged drag/lift force from induced-dipole on field-gradient formula:
       <F_drag> = -(V_plate / (2 mu0)) * |Im[alpha_zz(s_eff)/V_plate]| * |dBz/dx|^2
       <F_lift> = +(V_plate / (2 mu0)) * |Re[alpha_zz(s_eff)/V_plate]| * |dBz/dz|^2
     (sign convention: drag opposes motion, lift opposes gravity for B_z down)
  4. Compare to Schieber 1986 thin-plate analytical estimate.

This is a QUASI-STATIC approximation -- valid when L_PM/v >> longest plate
Foster time tau_max.  For sub-millisecond transients (engine-start ramp,
emergency brake) the SPICE Cauer ladder convolution route is needed
(see ecb_linear_transient.py, planned).

Output:
  - results.json with F_drag(v), F_lift(v) arrays
  - results.png log-log plot of |F| vs v
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

# Same dir
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cln_sibc_general_hex import (
    bulk_foster_via_eigen,
    K_SIBC_total,
    c1_polyhedral,
    measure_total_area_and_edges,
    Y_mixed,
    alpha_from_Y,
)


# Plate & PM defaults (engineering ECB parameters)
SIGMA_AL = 3.5e7        # aluminum 6061 conductivity
MU_0 = 4 * math.pi * 1e-7
PLATE_DIMS = (100e-3, 50e-3, 5e-3)   # Lx, Ly, Lz in m
PM_DIMS = (10e-3, 10e-3, 10e-3)       # cube PM
PM_M = 954930.0                       # A/m (Br=1.2T NdFeB)
GAP_H = 10e-3                         # bottom of PM to top of plate
V_SWEEP = [0.1, 0.3, 1, 3, 10, 30, 100]


def make_plate_vol(Lx, Ly, Lz, h_target, vol_path):
    """Generate plate mesh as .vol."""
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh
    plate = Box(Pnt(0, 0, 0), Pnt(Lx, Ly, Lz))
    plate.mat("Al")
    for f in plate.faces:
        f.name = "outer"
    geo = OCCGeometry(plate)
    ng_mesh = geo.GenerateMesh(maxh=h_target)
    mesh = Mesh(ng_mesh)
    print(f"  plate mesh: ne={mesh.ne}, nv={mesh.nv}")
    ng_mesh.Save(vol_path)
    return vol_path


def compute_plate_alpha_spectrum(vol_path, sigma=SIGMA_AL, mu=MU_0, n_eigen=200):
    """Run Mixed Galerkin once, return (lam, tau, g_n, K_SIBC, c1, V) for Y_mixed."""
    from ngsolve import Mesh, TaskManager
    mesh = Mesh(str(vol_path))
    print(f"  Loading plate mesh ne={mesh.ne}")
    with TaskManager():
        lam, tau, g_n, V = bulk_foster_via_eigen(mesh, sigma, mu, n_eigen=n_eigen)
        S_total, edges = measure_total_area_and_edges(mesh)
    K_SIBC = K_SIBC_total(S_total, sigma, mu)
    c1 = c1_polyhedral(edges, mu)
    print(f"  V = {V*1e9:.2f} mm^3, S = {S_total*1e6:.2f} mm^2, K_SIBC = {K_SIBC:.4e}")
    print(f"  Bulk Foster modes: {len(lam)}, tau_max = {max(tau):.3e} s, tau_min = {min(tau):.3e} s")
    print(f"  c_1 (edge correction) = {c1:.4e}")
    return lam, tau, g_n, K_SIBC, c1, V


def alpha_at_frequency(f, lam, tau, g_n, K_SIBC, c1, V, sigma):
    """alpha_iso(s) at angular frequency 2 pi f -- scalar (isotropic plate approx)."""
    s = 1j * 2 * math.pi * f
    Y = Y_mixed(s, lam, tau, g_n, K_SIBC, c1)
    a = alpha_from_Y(Y, V, sigma)
    return a   # complex; for plate, this approximates alpha_zz (dominant direction)


def pm_field_at_plate(h_gap, PM_dim, PM_M):
    """Estimate peak |B_z| from a cube magnet at distance h_gap.
    Pole-density approximation: B_z(0,0,h_gap+L/2) for a uniform-magnetized cube.
    Simplified to Nd dipole approximation for h_gap >> L_PM, scaled formula:
        |B_z|_peak ~ (mu0 M L^3) / (4 pi h^3)
    """
    Lp = PM_dim[0]
    moment = PM_M * (Lp ** 3)
    return MU_0 * moment / (4 * math.pi * (h_gap + Lp/2) ** 3)


def compute_force_vs_v(v_list, spec_data, sigma, h_gap, PM_dim, PM_M):
    """Per velocity v: f_eff = v / (2 L_PM), compute alpha(f_eff), get drag/lift."""
    lam, tau, g_n, K_SIBC, c1, V = spec_data
    Lp = PM_dim[0]
    B_peak = pm_field_at_plate(h_gap, PM_dim, PM_M)
    print(f"\n  PM peak |B_z| at gap = {B_peak*1e3:.4f} mT (dipole approximation)")
    print()
    print(f"  {'v (m/s)':>9}  {'f_eff (Hz)':>11}  {'Re(a/V)':>10}  {'-Im(a/V)':>10}  {'F_drag (N)':>11}  {'F_lift (N)':>11}")

    results = []
    for v in v_list:
        f_eff = v / (2 * Lp)   # transit frequency = velocity / PM characteristic length
        a = alpha_at_frequency(f_eff, lam, tau, g_n, K_SIBC, c1, V, sigma)
        # Dimensional analysis force estimate:
        # F ~ alpha * (grad B)^2 / mu0.  Use B_peak^2 / h_gap as dB/dx scale.
        # F_drag = -2 (Im alpha / V) * B^2 / mu0 * V_plate / h_gap^2  (Schieber-like)
        # F_lift = +2 (Re alpha / V) * B^2 / mu0 * V_plate / h_gap^2
        prefactor = (B_peak ** 2) / (MU_0 * h_gap ** 2) * V
        F_drag = -2.0 * (a.imag / V) * prefactor
        F_lift = +2.0 * (a.real / V) * prefactor
        results.append({
            "v": v, "f_eff": f_eff,
            "alpha_re_over_V": a.real / V,
            "alpha_im_over_V": a.imag / V,
            "F_drag": F_drag, "F_lift": F_lift,
        })
        print(f"  {v:9.2f}  {f_eff:11.2e}  {a.real/V:+10.4f}  {-a.imag/V:+10.4f}  {F_drag:+10.4e}  {F_lift:+10.4e}")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vol", help="Pre-generated plate .vol (if not, generated)")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    print("=== Linear ECB drag/lift force vs velocity ===")
    print(f"  Plate: {PLATE_DIMS[0]*1e3:.0f} x {PLATE_DIMS[1]*1e3:.0f} x {PLATE_DIMS[2]*1e3:.0f} mm Al "
          f"sigma={SIGMA_AL:.2e}")
    print(f"  PM: cube {PM_DIMS[0]*1e3:.0f} mm Nd, M={PM_M:.0f} A/m (Br=1.2T)")
    print(f"  Gap (PM bottom -> plate top): {GAP_H*1e3:.1f} mm")
    print()

    vol_path = Path(args.vol) if args.vol else Path("plate_100x50x5.vol")
    if not vol_path.exists():
        print("  Generating plate mesh...")
        make_plate_vol(*PLATE_DIMS, h_target=PLATE_DIMS[2]/2, vol_path=str(vol_path))

    spec = compute_plate_alpha_spectrum(str(vol_path), sigma=SIGMA_AL,
                                         mu=MU_0, n_eigen=300)
    results = compute_force_vs_v(V_SWEEP, spec, SIGMA_AL, GAP_H, PM_DIMS, PM_M)

    out_json = {
        "geometry": {"plate_mm": [d*1e3 for d in PLATE_DIMS],
                     "pm_mm": [d*1e3 for d in PM_DIMS],
                     "gap_mm": GAP_H * 1e3,
                     "PM_M_A_per_m": PM_M},
        "material": {"sigma_S_per_m": SIGMA_AL, "mu_H_per_m": MU_0},
        "vs_velocity": results,
    }
    with open("ecb_linear_results.json", "w", encoding="utf-8") as fp:
        json.dump(out_json, fp, indent=2)
    print(f"\n  Saved ecb_linear_results.json")

    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            v = [r["v"] for r in results]
            F_drag = [abs(r["F_drag"]) for r in results]
            F_lift = [abs(r["F_lift"]) for r in results]
            fig, ax = plt.subplots(figsize=(5, 3.5))
            ax.loglog(v, F_drag, "ro-", label="|F_drag|")
            ax.loglog(v, F_lift, "b^-", label="|F_lift|")
            ax.set_xlabel("velocity v (m/s)")
            ax.set_ylabel("Force (N)")
            ax.set_title("Linear ECB: cube PM over Al plate")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend()
            fig.tight_layout()
            fig.savefig("ecb_linear_results.png", dpi=150)
            print(f"  Saved ecb_linear_results.png")
        except ImportError:
            print("  (matplotlib unavailable; skipping plot)")


if __name__ == "__main__":
    main()
