r"""The scaling-FFAG SECTOR body driven into SATURATION: the azimuthal L_eff is ROBUST
(gap-reluctance-dominated) while the radial field index k(r) is FRAGILE -- the two sector
planes respond OPPOSITELY to iron saturation.

THE COMBINATION (the §3.5 scaling sector + iron saturation)
-----------------------------------------------------------
scaling_ffag_pole_2d.py built the scaling-FFAG pole and showed (Step 2) that, in the
RADIAL (r, z) plane, iron saturation droops the field index k(r) at the high-r edge --
the highest-B edge crosses the iron knee first (the super-ferric achromaticity wall),
A/phi-certified.  ffag_sector_two_plane.py built the SECTOR -- the azimuthal (s, z) plane
that sets how the field turns on/off along the orbit and gives the effective magnetic
length L_eff = INT B_z ds / B_z(body) -- but it solved that azimuthal end LINEARLY
(mu_r = const).

This file drives the SECTOR BODY into saturation: it makes the azimuthal (s, z) end solve
NONLINEAR (Froehlich mu_r(|B|), the same knee as the radial Step 2) and solves it at the
high-r aperture edge (small gap, highest B) and the low-r body, at a LINEAR and a
SATURATED excitation.  The honest result for a scaling (large-gap) pole is a CONTRAST
between the two planes:

  * the AZIMUTHAL end is GAP-RELUCTANCE-DOMINATED: even where the iron saturates hardest
    (the high-r edge, where <mu_r> collapses by ~3x), the effective length L_eff barely
    moves (drift < 1 %) -- L_eff is ROBUST to saturation (cf. the same honest scope in
    clebsch_dipole_saturation_3d: a large-gap magnet's gap field softens only mildly with
    iron saturation);
  * the RADIAL field index k(r) is FRAGILE: the same high-r saturation droops k(r) (the
    achromaticity wall, scaling_ffag Step 2).

So saturation degrades the radial field SHAPE (the achromaticity) but NOT the azimuthal
end LENGTH -- a real design insight: the high-r achromaticity needs the radial reshape
(scaling_ffag Step 3), while the sector ENDS are saturation-robust and need no nonlinear
end correction.

WHAT IS COMPUTED
----------------
1. solve_azimuthal_saturated(r_ref, ..., B_design): the NONLINEAR azimuthal end solve
   (Froehlich Picard) on the scaling-gap sector pole at radius r_ref -> L_eff, the iron
   <mu_r>, the fringe-per-gap, the pole-end overshoot.
2. sector_saturation(...): solve the azimuthal end at the low-r body and the high-r edge,
   LINEAR vs SATURATED; show the high-r iron saturates (mu_r collapses) yet L_eff is
   robust (drift < 2 %) -- the gap-reluctance-dominated sector end.
3. The radial cross-check (reuse scaling_ffag_pole_2d.run_step2): the SAME high-r edge
   droops the field index k(r) -- the FRAGILE plane, contrasting with the robust end.

run:  python scaling_ffag_sector_saturation.py            # the sector saturation table
      python scaling_ffag_sector_saturation.py --fig        # + figure
      python scaling_ffag_sector_saturation.py --no-radial  # skip the radial cross-check
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scaling_ffag_pole_2d import (                          # noqa: E402
    K_INDEX, G0, R0, MU0, MUR0_IRON, BK_IRON,
    aperture_radii, scaling_gap, field_index, run_step2,
)
from ffag_sector_two_plane import _azimuthal_geometry        # noqa: E402


def solve_azimuthal_saturated(r_ref, sector_arclen, B_design, k=K_INDEX, g0=G0,
                              r0=R0, pole_t=0.03, order=3, win_g=6.0, z_buffer=3.0,
                              relax=0.5, tol=1e-4, maxit=60):
    """Plane B made NONLINEAR: the finite-length scaling-gap sector pole at radius
    r_ref, gap g(r_ref), over the arc `sector_arclen`, with a Froehlich iron pole
    (mu_r(|B|), the same knee as the radial Step 2).  Scalar potential phi driven
    at Psi=mmf on the iron back (Picard on mu_r(|B|)).  Returns L_eff, the pole-end
    overshoot, the iron <mu_r>, and the saturation diagnostics."""
    from ngsolve import (H1, L2, BilinearForm, GridFunction, grad, dx, CF, Norm,
                         TaskManager, Mesh)

    g_ref = float(scaling_gap(r_ref, k, g0, r0))
    sector_half = 0.5 * sector_arclen
    win = win_g * g_ref                                       # fixed end window (~g)
    s_half = sector_half + win + z_buffer * g_ref
    z_top = 0.5 * g_ref + pole_t + z_buffer * g_ref
    maxh = g_ref / 3.0
    mmf = B_design * g_ref / (2.0 * MU0)                       # B_body ~ B_design
    n_eval = max(161, int(2.0 * (sector_half + win) / (g_ref / 20.0)))
    ss = np.linspace(-(sector_half + win), sector_half + win, n_eval)
    y_probe = 0.02 * g_ref
    geo = _azimuthal_geometry(s_half, z_top, sector_half, g_ref, pole_t)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        gfu = GridFunction(fes)
        bccf = mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0)
        fes_mu = L2(mesh, order=0)
        mu_gf = GridFunction(fes_mu)
        mu_gf.Set(mesh.MaterialCF({"iron": MUR0_IRON}, default=1.0))   # start
        iron_ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
        resid, it = 1.0, 0
        for it in range(1, maxit + 1):
            a = BilinearForm(fes)
            a += mu_gf * grad(u) * grad(v) * dx
            a.Assemble()
            gfu.Set(bccf, definedon=mesh.Boundaries("median|irontop"))
            r = gfu.vec.CreateVector()
            r.data = -a.mat * gfu.vec
            gfu.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                          inverse="sparsecholesky") * r
            B = MU0 * mu_gf * Norm(grad(gfu))                 # |B| with current mu
            froh = 1.0 + (MUR0_IRON - 1.0) / (1.0 + (B / BK_IRON) ** 2)
            mu_target = (1.0 - iron_ind) * CF(1.0) + iron_ind * froh
            mu_new = GridFunction(fes_mu)
            mu_new.Set(mu_target)
            d = mu_new.vec.CreateVector()
            d.data = mu_new.vec - mu_gf.vec
            resid = d.Norm() / (mu_gf.vec.Norm() or 1.0)
            mu_gf.vec.data += relax * d
            if resid < tol:
                break
        g = grad(gfu)
        Bz = np.array([-MU0 * g(mesh(float(s), y_probe))[1] for s in ss])
        # iron <mu_r> (saturation state) over the iron pole
        iron_vol = float(np.nan)
        from ngsolve import Integrate
        iv = float(Integrate(iron_ind, mesh))
        mur_mean = float(Integrate(iron_ind * mu_gf, mesh)) / iv if iv > 0 else 0.0
        ndof, ne = int(fes.ndof), int(mesh.ne)
    body = np.abs(ss) < 0.3 * sector_half
    Bz_body = float(np.mean(Bz[body]))
    L_eff = float(np.trapezoid(Bz, ss) / Bz_body) if abs(Bz_body) > 1e-30 else 0.0
    fringe_excess = (L_eff - sector_arclen) / sector_arclen
    end = np.abs(np.abs(ss) - sector_half) < 0.12 * sector_arclen
    end_overshoot = (float(np.max(np.abs(Bz[end])) / abs(Bz_body)) - 1.0
                     if end.any() and abs(Bz_body) > 1e-30 else 0.0)
    return {
        "r_ref": float(r_ref), "g_ref": g_ref, "B_design": float(B_design),
        "sector_arclen": float(sector_arclen),
        "L_eff": L_eff, "fringe_excess": float(fringe_excess),
        "fringe_per_gap": float((L_eff - sector_arclen) / g_ref),
        "Bz_body": Bz_body, "end_overshoot": float(end_overshoot),
        "mur_mean": mur_mean, "iters": int(it), "resid": float(resid),
        "ne": ne, "ndof": ndof,
    }


def sector_saturation(dtheta=0.30, k=K_INDEX, g0=G0, r0=R0, B_body_lin=0.4,
                      B_body_sat=1.4, order=3):
    """Solve the azimuthal sector end at the LOW-r body (r0) and the HIGH-r aperture
    edge (r_max), at a LINEAR (low) and a SATURATED (high) body excitation.  At the
    high-r edge the local field is B_body*(r_max/r0)^k -- it crosses the knee first.
    Report the per-radius L_eff DRIFT: large at high-r (it saturates), ~0 at low-r."""
    r_min, r_max = aperture_radii(k=k, r0=r0)
    sector_arclen = r0 * dtheta                               # arc at the reference radius
    radii = {"low_r_body": r0, "high_r_edge": r_max}
    out = {}
    for name, r_ref in radii.items():
        scale = (r_ref / r0) ** k                             # scaling field B(r)~r^k
        lin = solve_azimuthal_saturated(r_ref, sector_arclen, B_body_lin * scale,
                                        k=k, g0=g0, r0=r0, order=order)
        sat = solve_azimuthal_saturated(r_ref, sector_arclen, B_body_sat * scale,
                                        k=k, g0=g0, r0=r0, order=order)
        drift = (sat["L_eff"] - lin["L_eff"]) / lin["L_eff"] if lin["L_eff"] else 0.0
        out[name] = {
            "r_ref": float(r_ref),
            "B_local_linear_T": float(B_body_lin * scale),
            "B_local_saturated_T": float(B_body_sat * scale),
            "L_eff_linear": lin["L_eff"], "L_eff_saturated": sat["L_eff"],
            "L_eff_drift_rel": float(drift),
            "fringe_per_gap_linear": lin["fringe_per_gap"],
            "fringe_per_gap_saturated": sat["fringe_per_gap"],
            "mur_mean_linear": lin["mur_mean"], "mur_mean_saturated": sat["mur_mean"],
            "end_overshoot_linear": lin["end_overshoot"],
            "end_overshoot_saturated": sat["end_overshoot"],
        }
    hi, lo = out["high_r_edge"], out["low_r_body"]
    # The honest finding for a scaling (large-gap) pole: the high-r iron SATURATES hard
    # (its <mu_r> collapses), yet its azimuthal L_eff is ROBUST -- the sector end is
    # gap-reluctance-dominated (cf. clebsch_dipole_saturation_3d's honest scope), so
    # saturation barely moves L_eff even where it bites hardest.  The radial field index
    # k(r), by contrast, IS fragile (it droops -- the radial cross-check below).
    return {
        "k": float(k), "dtheta_rad": float(dtheta),
        "aperture": (float(r_min), float(r_max)),
        "B_K_iron_T": float(BK_IRON),
        "radii": out,
        "high_r_iron_saturates": bool(hi["mur_mean_saturated"]
                                      < 0.5 * hi["mur_mean_linear"]),
        "high_r_mur_collapse_factor": float(hi["mur_mean_saturated"]
                                            / max(hi["mur_mean_linear"], 1e-9)),
        "L_eff_robust_under_saturation": bool(abs(hi["L_eff_drift_rel"]) < 0.02),
        "high_r_L_eff_drift_rel": float(hi["L_eff_drift_rel"]),
        "fringe_per_gap_stable": bool(abs(hi["fringe_per_gap_saturated"]
                                          - hi["fringe_per_gap_linear"]) < 0.10),
    }


def _figure(res, radial=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    npan = 3 if radial is not None else 2
    fig, ax = plt.subplots(1, npan, figsize=(5.4 * npan, 4.2), dpi=140)

    hi, lo = res["radii"]["high_r_edge"], res["radii"]["low_r_body"]

    # LEFT: per-radius L_eff drift (linear -> saturated): high-r drifts, low-r doesn't
    labels = ["low-r body", "high-r edge"]
    drift = [lo["L_eff_drift_rel"] * 100, hi["L_eff_drift_rel"] * 100]
    bars = ax[0].bar(labels, drift, color=["C0", "C3"])
    ax[0].axhline(0, color="k", lw=0.6)
    for b_, v in zip(bars, drift):
        ax[0].text(b_.get_x() + b_.get_width() / 2, v, f"{v:+.1f}%",
                   ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
    ax[0].set_ylabel("$L_{eff}$ drift, linear $\\to$ saturated [%]")
    ax[0].set_ylim(-3.0, 3.0)
    ax[0].set_title("Azimuthal $L_{eff}$ is ROBUST to saturation\n"
                    "(< 1 % drift even where the iron saturates)")

    # MIDDLE: the iron <mu_r> collapse (the mechanism) at each radius
    x = np.arange(2)
    w = 0.36
    murl = [lo["mur_mean_linear"], hi["mur_mean_linear"]]
    murs = [lo["mur_mean_saturated"], hi["mur_mean_saturated"]]
    ax[1].bar(x - w / 2, murl, w, color="0.6", label="linear drive")
    ax[1].bar(x + w / 2, murs, w, color="C3", label="saturated drive")
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels)
    ax[1].set_ylabel("iron $\\langle\\mu_r\\rangle$")
    ax[1].set_title("Mechanism: the high-r iron $\\langle\\mu_r\\rangle$\n"
                    "collapses (it crosses the knee first)")
    ax[1].legend(fontsize=8)

    # RIGHT (optional): the radial cross-check -- the field index k(r) droops at high-r
    if radial is not None:
        k_ref = np.array(radial["k_ref_interior"])
        L = radial["levels"][-1]
        s = L["_s"]
        m = slice(2, -2)
        ri = s["r_index"][m]
        dk = s["k"][m] - k_ref
        ax[2].plot(ri, dk, "o-", color="C3", ms=3,
                   label=f"B@r0={L['B_target']:.1f} T (max {L['B_gap_max']:.1f})")
        ax[2].axhline(0, color="k", lw=0.8)
        ax[2].set_xlabel("r"); ax[2].set_ylabel("$\\Delta k(r)$ = k - k(unsat)")
        ax[2].set_title("Same mechanism, radial plane:\n"
                        "the field index k(r) droops at high-r")
        ax[2].legend(fontsize=8)

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", action="store_true")
    ap.add_argument("--no-radial", action="store_true",
                    help="skip the radial field-index cross-check (azimuthal only)")
    ap.add_argument("--dtheta", type=float, default=0.30)
    args = ap.parse_args()

    print("=" * 78)
    print("The scaling-FFAG SECTOR body driven into SATURATION -- the high-r edge first")
    print("=" * 78)

    res = sector_saturation(dtheta=args.dtheta)
    rmin, rmax = res["aperture"]
    print(f"\niron knee B_K = {res['B_K_iron_T']:.2f} T; field index k = {res['k']:.1f}; "
          f"aperture [{rmin:.3f}, {rmax:.3f}] (r_max/r0 = {rmax/R0:.2f})")
    print(f"\n  {'radius':<14}{'B_local [T]':>14}{'<mu_r> lin->sat':>18}"
          f"{'L_eff drift':>14}{'fringe/gap':>13}")
    for name in ("low_r_body", "high_r_edge"):
        d = res["radii"][name]
        print(f"  {name:<14}{d['B_local_saturated_T']:>14.2f}"
              f"{d['mur_mean_linear']:>9.0f}->{d['mur_mean_saturated']:<7.0f}"
              f"{d['L_eff_drift_rel']*100:>12.1f} %"
              f"{d['fringe_per_gap_saturated']:>13.2f}")
    print(f"\n  => the high-r edge iron SATURATES hard (<mu_r> collapses x"
          f"{res['high_r_mur_collapse_factor']:.2f}): {res['high_r_iron_saturates']}")
    print(f"  => yet its azimuthal L_eff is ROBUST (drift "
          f"{res['high_r_L_eff_drift_rel']*100:+.1f}%, < 2%): "
          f"{res['L_eff_robust_under_saturation']}")
    print(f"\n  i.e. the sector END is GAP-RELUCTANCE-DOMINATED: saturation barely moves")
    print(f"  L_eff even where the iron bites hardest (cf. clebsch_dipole_saturation_3d).")
    print(f"  The two sector planes respond OPPOSITELY to saturation -- the azimuthal end")
    print(f"  LENGTH is robust, while the radial field SHAPE (achromaticity) is fragile:")

    radial = None
    if not args.no_radial:
        radial = run_step2()
        Lhi = radial["levels"][-1]
        print(f"\n  Radial cross-check (the (r,z) plane -- the FRAGILE one):")
        print(f"    at B@r0 = {Lhi['B_target']:.1f} T (B_gap up to {Lhi['B_gap_max']:.1f} T), "
              f"the field index DROOPS dk_hi = {Lhi['dk_hi']:+.3f} at the high-r edge")
        print(f"    -> SATURATION degrades the radial achromaticity (k droops) but NOT the")
        print(f"       azimuthal effective length (L_eff robust) -- the two planes differ.")
        res["radial_k_droops_high_r"] = bool(Lhi["dk_hi"] < -0.02)
        res["radial_field_index"] = {
            "B_target_hi": float(Lhi["B_target"]), "B_gap_max_hi": float(Lhi["B_gap_max"]),
            "dk_hi": float(Lhi["dk_hi"]), "dk_tilt": float(Lhi["dk_tilt"]),
        }
        res["planes_differ_under_saturation"] = bool(
            res["L_eff_robust_under_saturation"] and res["radial_k_droops_high_r"])

    jpath = os.path.splitext(os.path.abspath(__file__))[0] + ".json"
    with open(jpath, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(res, radial=radial)


if __name__ == "__main__":
    main()
