r"""The SPECTROMETER end pack, NONLINEAR: the pole-tip corner is a saturable throat.

THE PHYSICS (why a spectrometer end needs the NONLINEAR design)
---------------------------------------------------------------
A large bending spectrometer dipole runs near the iron knee.  endpack_two_plane.py
found (LINEAR, the high-mu equipotential limit) that a hard-cut pole END concentrates
flux at the tip CORNER: the corner field is

    B_corner = kappa * B_gap ,   kappa = tip_enhancement ~ 1.11   (geometry-only)

In a saturating iron that kappa > 1 means the CORNER reaches the iron knee FIRST --
the corner saturates at a gap field

    B_gap_knee = B_K / kappa ~ 1.5 / 1.11 ~ 1.35 T

i.e. ~11% BELOW the bulk iron knee B_K = 1.5 T.  Above that the corner mu_r collapses,
the flux there can no longer follow the pole edge, and the EFFECTIVE FIELD BOUNDARY
(EFB ~ the effective magnetic length L_eff) DRIFTS with excitation.  For a spectrometer
that is fatal: the pole-edge angle's EDGE FOCUSING (the vertical focusing tan(beta)/rho)
depends on the EFB, so a drifting EFB means the optics CHANGE with the field setting.

THE DESIGN (the Rogowski chamfer is the corner-throat width knob)
----------------------------------------------------------------
The corner is exactly a Chaplygin saturable THROAT (clebsch_dipole_saturation_2d): the
flux concentration kappa is the throat's inverse cross-section, and it saturates first.
The Rogowski end chamfer (endpack_two_plane's s-y design) WIDENS that throat -- it lowers
kappa -- so the corner knee B_gap_knee = B_K/kappa RISES toward the bulk knee B_K.  The
SAME chamfer that zeroes the LINEAR corner over-field (the cosmetic field-quality lever)
REMOVES the PREMATURE corner saturation (the hard engineering lever) and keeps the EFB
-- and the edge focusing -- STABLE up to the bulk limit.  Linear and nonlinear levers
POINT THE SAME WAY; saturation gives the chamfer its hard justification.

WHAT IS COMPUTED
----------------
1. design(...) -- the nonlinear end-pack MAP at LINEAR cost: reuse the linear
   equipotential corner concentration kappa(chamfer) (endpack_two_plane), overlay the
   Froehlich iron BH (clebsch_dipole_saturation_3d) -> the corner knee B_gap_knee(chamfer)
   and the chamfer DEPTH needed to keep the corner unsaturated at a target operating
   field B_op.  Milliseconds of overlay on one ~8 s linear sweep -- the Chaplygin
   "nonlinear analysis done linearly" applied to the END corner.

This is a DESIGN-GRADE lumped model (the same lumped-magnetic-circuit class as
clebsch_dipole_saturation_2d, ~10% vs FEM), and its two ingredients are
independently VALIDATED elsewhere in the repo:
  - the corner concentration kappa is the LINEAR equipotential tip_enhancement
    (endpack_two_plane.py, golden-tested) -- geometry-only, drive-agnostic;
  - the Froehlich BH + the well-conditioned B-input A-formulation that backs the
    iron saturation are validated in clebsch_dipole_saturation_3d.py (the documented
    cure: the reduced-Omega mu(|H|) Picard STALLS at high mu, the A-formulation does
    not).
So the composition "B_corner = kappa*B_gap until B_K" is well-founded.  A fully
coil-driven 3-D corner-saturation FEM (the equipotential/MMF drive forces flux
across the gap -- a uniform applied field does NOT reproduce the corner concentration,
and the coil's Biot-Savart B_s projection is the expensive serial step) is the
documented expensive extension, not run here.

run:  python endpack_spectrometer_saturation.py            # the linear-cost design map
      python endpack_spectrometer_saturation.py --fig        # + figure
      python endpack_spectrometer_saturation.py --b-op 1.6   # a higher operating field
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import endpack_two_plane as ep                       # noqa: E402  linear corner concentration
import clebsch_dipole_saturation_3d as sat3          # noqa: E402  Froehlich BH + A-formulation
import accel_pole_ends_fem as fem                    # noqa: E402  H-frame mesh + coil

BK = sat3.BK                 # 1.5 T  iron half-mu (Froehlich) knee
MUR0 = sat3.MUR0             # 2000   unsaturated mu_r
G2 = ep.G2                   # 0.020  half-gap
L_BEAM = fem.L_BEAM


# ============================================================
# the nonlinear end-pack design MAP at linear cost
# ============================================================

def design(B_op=1.45, fast=True):
    """The nonlinear end-pack design map.  Reuse the LINEAR corner concentration
    kappa(chamfer) (endpack_two_plane's equipotential depth sweep) and overlay the
    Froehlich iron knee: the corner saturates at B_gap = B_K/kappa, BELOW the bulk
    knee B_K.  Find the chamfer depth that keeps the corner unsaturated at B_op."""
    d = ep.design_two_plane(fast=fast)
    sweep = d["reflect_3d"]["depth_sweep"]
    depth = np.array([s["chamfer_depth_m"] for s in sweep])
    kappa = np.array([s["tip_enhancement"] for s in sweep])   # corner concentration (geom-only)
    knee = BK / kappa                                          # corner gap-field knee per chamfer

    order = np.argsort(kappa)                                  # kappa decreases with depth
    k_sorted, d_sorted = kappa[order], depth[order]
    kappa_req = BK / B_op                                      # corner not the limiter at B_op
    in_rng = k_sorted.min() <= kappa_req <= k_sorted.max()
    depth_req = float(np.interp(kappa_req, k_sorted, d_sorted)) if in_rng else float("nan")
    in_cos = k_sorted.min() <= 1.0 <= k_sorted.max()
    depth_cosmetic = float(np.interp(1.0, k_sorted, d_sorted)) if in_cos else float("nan")

    flat_kappa = float(kappa[0])
    flat_knee = float(knee[0])
    bulk_margin = float((BK - flat_knee) / BK)                 # how much earlier the corner saturates

    rows = [{"chamfer_depth_m": float(dd), "kappa": float(kk),
             "corner_knee_Bgap_T": float(BK / kk),
             "corner_saturated_at_Bop": bool(kk > kappa_req)}
            for dd, kk in zip(depth, kappa)]

    return {
        "B_K_T": float(BK), "B_op_T": float(B_op),
        "flat_corner_kappa": flat_kappa,
        "flat_corner_knee_Bgap_T": flat_knee,           # ~1.35 T -- corner saturates first
        "corner_premature_fraction": bulk_margin,       # ~0.11 below the bulk knee
        "kappa_required_for_Bop": float(kappa_req),
        "depth_required_for_Bop_m": depth_req,          # chamfer to clear B_op
        "depth_linear_cosmetic_m": depth_cosmetic,      # chamfer that zeroes the over-field
        "flat_corner_saturates_below_Bop": bool(flat_knee < B_op),
        "sweep": rows,
        "n_linear_solves": int(len(sweep)),             # the only FEM-grade cost (the linear sweep)
    }


# ============================================================
# figure + main
# ============================================================

def _figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.0), dpi=140)

    # LEFT: the corner knee B_gap_knee = B_K/kappa vs chamfer depth
    sw = res["sweep"]
    dmm = [r["chamfer_depth_m"] * 1e3 for r in sw]
    knee = [r["corner_knee_Bgap_T"] for r in sw]
    ax[0].plot(dmm, knee, "o-", color="C3", label="corner knee $B_{gap}=B_K/\\kappa$")
    ax[0].axhline(res["B_K_T"], color="0.4", lw=1, ls="--",
                  label=f"bulk iron knee $B_K$={res['B_K_T']:.2f} T")
    ax[0].axhline(res["B_op_T"], color="C2", lw=1, ls=":",
                  label=f"operating $B_{{op}}$={res['B_op_T']:.2f} T")
    ax[0].set_xlabel("end chamfer depth [mm]")
    ax[0].set_ylabel("gap field at corner saturation [T]")
    ax[0].set_title(f"The corner saturates EARLY (flat: {res['flat_corner_knee_Bgap_T']:.2f} T,\n"
                    f"{res['corner_premature_fraction']*100:.0f}% below $B_K$); "
                    f"the chamfer raises it")
    ax[0].legend(fontsize=8)

    # MIDDLE: kappa (corner concentration) vs depth -- the throat-width knob
    kap = [r["kappa"] for r in sw]
    ax[1].plot(dmm, kap, "s-", color="C0", label="corner concentration $\\kappa$")
    ax[1].axhline(1.0, color="k", lw=0.6)
    ax[1].axhline(res["kappa_required_for_Bop"], color="C2", lw=1, ls=":",
                  label=f"$\\kappa$ needed for $B_{{op}}$ = {res['kappa_required_for_Bop']:.3f}")
    if res["depth_required_for_Bop_m"] == res["depth_required_for_Bop_m"]:
        ax[1].axvline(res["depth_required_for_Bop_m"] * 1e3, color="C2", lw=1.2, ls="--",
                      label=f"required depth {res['depth_required_for_Bop_m']*1e3:.1f} mm")
    ax[1].set_xlabel("end chamfer depth [mm]")
    ax[1].set_ylabel("corner concentration $\\kappa=B_{corner}/B_{gap}$")
    ax[1].set_title("The Rogowski chamfer is the corner-throat\n"
                    "width knob ($\\kappa\\downarrow$ raises the knee)")
    ax[1].legend(fontsize=8)

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", action="store_true")
    ap.add_argument("--b-op", type=float, default=1.45, help="operating gap field [T]")
    args = ap.parse_args()

    print("=" * 76)
    print("The SPECTROMETER end pack, NONLINEAR -- the pole-tip corner is a saturable throat")
    print("=" * 76)

    res = design(B_op=args.b_op)
    print(f"\niron: bulk knee B_K = {res['B_K_T']:.2f} T (Froehlich mu_r0={MUR0:.0f}); "
          f"operating B_op = {res['B_op_T']:.2f} T")
    print(f"\nFLAT end -- the corner CONCENTRATES flux (kappa = "
          f"{res['flat_corner_kappa']:.3f}) and SATURATES FIRST:")
    print(f"  corner knee  B_gap = B_K/kappa = {res['flat_corner_knee_Bgap_T']:.2f} T "
          f"-- {res['corner_premature_fraction']*100:.0f}% BELOW the bulk knee B_K")
    print(f"  -> at B_op = {res['B_op_T']:.2f} T the flat corner is "
          f"{'SATURATED' if res['flat_corner_saturates_below_Bop'] else 'OK'} "
          f"(the EFB / edge focusing drifts).")
    print(f"\nthe Rogowski chamfer is the corner-throat width knob (kappa -> knee = B_K/kappa):")
    print(f"  {'depth [mm]':<12}{'kappa':<9}{'corner knee [T]':<18}{'corner sat @ B_op?'}")
    for r in res["sweep"]:
        print(f"  {r['chamfer_depth_m']*1e3:<12.1f}{r['kappa']:<9.3f}"
              f"{r['corner_knee_Bgap_T']:<18.2f}{'YES' if r['corner_saturated_at_Bop'] else 'no'}")
    if res["depth_required_for_Bop_m"] == res["depth_required_for_Bop_m"]:
        print(f"  -> need chamfer depth >= {res['depth_required_for_Bop_m']*1e3:.1f} mm "
              f"(kappa <= {res['kappa_required_for_Bop']:.3f}) to clear B_op without "
              f"corner saturation;")
        print(f"     the LINEAR cosmetic optimum (zero over-field, kappa=1) is "
              f"{res['depth_linear_cosmetic_m']*1e3:.1f} mm -- the SAME lever, now with a "
              f"hard saturation justification.")
    print(f"\n  => the whole nonlinear end-pack map = {res['n_linear_solves']} linear "
          f"equipotential solves + a BH overlay (the Chaplygin 'nonlinear done linearly'\n"
          f"     lever applied to the END corner).  Design-grade (lumped-circuit class,\n"
          f"     ~10% vs FEM like clebsch_dipole_saturation_2d); its components are\n"
          f"     independently validated -- kappa by the linear equipotential\n"
          f"     (endpack_two_plane), the BH + A-formulation by clebsch_dipole_saturation_3d.")

    jpath = os.path.splitext(os.path.abspath(__file__))[0] + ".json"
    with open(jpath, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(res)


if __name__ == "__main__":
    main()
