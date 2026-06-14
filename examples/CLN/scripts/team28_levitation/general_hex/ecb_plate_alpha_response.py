"""ecb_plate_alpha_response.py -- conducting plate alpha(s) frequency response for ECB design.

Linear eddy-current brake: a PM moves at velocity v over a flat conducting plate.
The plate's polarizability tensor alpha(s) determines how it responds at the
"effective transit frequency" f_eff = v / (2 L_PM).  Re[alpha] -> lift component,
-Im[alpha] -> drag component.

This script focuses on the CORE building block:

  - plate Mixed-Galerkin alpha(s) over the design frequency band
  - identify the resonance / breakpoint frequency f_break (where drag peaks)
  - identify the levitation crossover f_lift (where alpha switches sign)
  - export the SPICE-extractable Foster spectrum + Mellin tail

For computing the actual brake FORCE F_drag(v), F_lift(v), the standard
recipes are:

  1. Image-method (analytical, thin-plate limit): Schieber 1986,
     F_drag = (3 mu0 m_PM^2 / 32 pi) g(xi) / h^4
     where g(xi) is a dimensionless transit-Reynolds function tied to
     -Im[alpha(s_eff)] / V_plate.
  2. Maxwell stress tensor integration on a closed surface around the PM
     (the rigorous answer), using B_ext + B_induced (from alpha(s)).
  3. Schur-Cauer time-domain (transient): inverse-Laplace alpha(s) into a
     SPICE Cauer ladder + Warburg block, drive with B_ext(t), integrate.

The alpha(s) computed here is the SINGLE input to all three routes;
this script just delivers it cleanly so the force calculation can be
plugged in.
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cln_sibc_general_hex import (
    bulk_foster_via_eigen,
    K_SIBC_total,
    c1_polyhedral,
    measure_total_area_and_edges,
    Y_mixed,
    alpha_from_Y,
)


SIGMA_AL = 3.5e7
MU_0 = 4 * math.pi * 1e-7
PLATE_DIMS = (100e-3, 50e-3, 5e-3)
PM_LEN = 10e-3                       # for effective-frequency mapping


def make_plate_vol(Lx, Ly, Lz, vol_path):
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh
    plate = Box(Pnt(0, 0, 0), Pnt(Lx, Ly, Lz))
    plate.mat("Al")
    for f in plate.faces:
        f.name = "outer"
    geo = OCCGeometry(plate)
    ng_mesh = geo.GenerateMesh(maxh=Lz / 2)
    mesh = Mesh(ng_mesh)
    ng_mesh.Save(vol_path)
    return mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vol")
    parser.add_argument("--n-eigen", type=int, default=400)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    print("=== ECB plate alpha(s) frequency response ===")
    print(f"  Plate: {PLATE_DIMS[0]*1e3:.0f} x {PLATE_DIMS[1]*1e3:.0f} x {PLATE_DIMS[2]*1e3:.0f} mm Al")
    print(f"  PM characteristic length L_PM = {PM_LEN*1e3:.1f} mm")
    print(f"  Effective frequency mapping: f_eff = v / (2 L_PM)")
    print()

    vol_path = Path(args.vol) if args.vol else Path("plate_100x50x5.vol")
    if not vol_path.exists():
        print("  Generating plate mesh...")
        from ngsolve import Mesh
        m = make_plate_vol(*PLATE_DIMS, vol_path=str(vol_path))
        print(f"  ne = {m.ne}, nv = {m.nv}")

    from ngsolve import Mesh, TaskManager
    mesh = Mesh(str(vol_path))
    print(f"  Mesh loaded: ne = {mesh.ne}")
    with TaskManager():
        lam, tau, g_n, V = bulk_foster_via_eigen(mesh, SIGMA_AL, MU_0,
                                                   n_eigen=args.n_eigen)
        S_total, edges = measure_total_area_and_edges(mesh)
    K_SIBC = K_SIBC_total(S_total, SIGMA_AL, MU_0)
    c1 = c1_polyhedral(edges, MU_0)

    print(f"  V = {V*1e9:.2f} mm^3, S = {S_total*1e6:.2f} mm^2")
    print(f"  K_SIBC = {K_SIBC:.4e}, c_1 = {c1:.4e}")
    print(f"  Foster modes: {len(lam)}, completeness = {np.sum(g_n)/(SIGMA_AL*V):.4f}")
    print(f"  tau range: [{min(tau)*1e6:.2f}, {max(tau)*1e6:.2f}] us  -->  ")
    print(f"  break freq range: [{1/(2*math.pi*max(tau)):.2f}, {1/(2*math.pi*min(tau)):.2f}] Hz")
    print()

    # Frequency sweep
    f_arr = np.logspace(-1, 8, 73)   # 0.1 Hz to 100 MHz
    alpha_re = np.zeros_like(f_arr)
    alpha_im = np.zeros_like(f_arr)
    for i, f in enumerate(f_arr):
        s = 1j * 2 * math.pi * f
        Y = Y_mixed(s, lam, tau, g_n, K_SIBC, c1)
        a = alpha_from_Y(Y, V, SIGMA_AL)
        alpha_re[i] = a.real / V
        alpha_im[i] = a.imag / V

    # Find drag peak (where -Im[alpha]/V is max) and lift crossover (Re[alpha]/V = 0.5)
    i_drag_peak = int(np.argmax(-alpha_im))
    f_drag_peak = f_arr[i_drag_peak]
    v_drag_peak = f_drag_peak * 2 * PM_LEN
    re_at_drag = alpha_re[i_drag_peak]
    im_at_drag = alpha_im[i_drag_peak]

    # Lift crossover: smallest f where Re[alpha/V] >= 0.5 (half perfect-conductor)
    i_lift = next((i for i, x in enumerate(alpha_re) if x >= 0.5), None)
    if i_lift is not None:
        f_lift = f_arr[i_lift]
        v_lift = f_lift * 2 * PM_LEN
    else:
        f_lift, v_lift = None, None

    print("=== ECB design break points ===")
    print(f"  Drag peak f = {f_drag_peak:.2f} Hz  =>  v_peak = {v_drag_peak:.2f} m/s")
    print(f"    at this f: Re(alpha/V) = {re_at_drag:+.3f}, -Im(alpha/V) = {-im_at_drag:+.3f}")
    if f_lift:
        print(f"  Lift crossover (Re/V = 0.5) f = {f_lift:.2f} Hz  =>  v_lift = {v_lift:.2f} m/s")
    else:
        print(f"  Lift crossover (Re/V = 0.5) NOT reached in scanned f-range")
    print()
    print("  v sweep (m/s):       0.1     1.0    10.0   100.0  1000.0")
    print("  f_eff (Hz):       {:5.2e}".format(0.1/(2*PM_LEN)),
          "{:5.2e}".format(1.0/(2*PM_LEN)),
          "{:5.2e}".format(10.0/(2*PM_LEN)),
          "{:5.2e}".format(100.0/(2*PM_LEN)),
          "{:5.2e}".format(1000.0/(2*PM_LEN)))
    sample_v = [0.1, 1.0, 10.0, 100.0, 1000.0]
    re_at_v, im_at_v = [], []
    for v in sample_v:
        f = v / (2 * PM_LEN)
        s = 1j * 2 * math.pi * f
        a = alpha_from_Y(Y_mixed(s, lam, tau, g_n, K_SIBC, c1), V, SIGMA_AL)
        re_at_v.append(a.real / V)
        im_at_v.append(a.imag / V)
    print("  Re(a/V):       {:+5.3f}  {:+5.3f}  {:+5.3f}  {:+5.3f}  {:+5.3f}".format(*re_at_v))
    print("  -Im(a/V):      {:+5.3f}  {:+5.3f}  {:+5.3f}  {:+5.3f}  {:+5.3f}".format(*[-x for x in im_at_v]))

    out = {
        "plate_dims_mm": [d*1e3 for d in PLATE_DIMS],
        "sigma_S_per_m": SIGMA_AL,
        "PM_characteristic_length_mm": PM_LEN * 1e3,
        "Foster_modes": int(len(lam)),
        "completeness": float(np.sum(g_n) / (SIGMA_AL * V)),
        "K_SIBC": K_SIBC,
        "c_1": c1,
        "f_drag_peak_Hz": f_drag_peak,
        "v_drag_peak_m_per_s": v_drag_peak,
        "f_lift_crossover_Hz": f_lift,
        "v_lift_crossover_m_per_s": v_lift,
        "freq_sweep_Hz": f_arr.tolist(),
        "alpha_re_over_V": alpha_re.tolist(),
        "alpha_im_over_V": alpha_im.tolist(),
        "tau_n": tau.tolist(),
        "g_n": g_n.tolist(),
    }
    with open("ecb_plate_alpha.json", "w", encoding="utf-8") as fp:
        json.dump(out, fp, indent=2)
    print(f"\n  Saved ecb_plate_alpha.json")

    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.semilogx(f_arr, alpha_re, "b-", label=r"Re[$\alpha$]/V (lift response)")
            ax.semilogx(f_arr, -alpha_im, "r-", label=r"$-$Im[$\alpha$]/V (drag response)")
            ax.axvline(f_drag_peak, color="r", ls="--", alpha=0.5, label=f"drag peak f={f_drag_peak:.1f} Hz")
            if f_lift:
                ax.axvline(f_lift, color="b", ls="--", alpha=0.5, label=f"lift crossover f={f_lift:.1f} Hz")
            ax.set_xlabel("frequency f (Hz)")
            ax.set_ylabel(r"$\alpha(s) / V_{plate}$")
            ax.set_title(f"ECB plate {PLATE_DIMS[0]*1e3:.0f}x{PLATE_DIMS[1]*1e3:.0f}x{PLATE_DIMS[2]*1e3:.0f} mm Al, alpha(s) response")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig("ecb_plate_alpha.png", dpi=150)
            print(f"  Saved ecb_plate_alpha.png")
        except ImportError:
            print("  matplotlib unavailable; skipping plot")


if __name__ == "__main__":
    main()
