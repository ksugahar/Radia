"""Rung 1: SynRM flux-channel iron sizing under a saturation cap.

The channel that carries d-axis flux around a flux barrier is, at Rung-1
abstraction, a 90-degree pure turn: walls are flux lines, terminals are MMF
equipotentials (gap entry idealised + q-axis symmetry plane), carried flux
Phi_d fixed, |B| <= B_cap everywhere, and the design question is the iron.

STRUCTURAL FACT exploited and verified here: a pure turn with flux-line walls
is EXACTLY solvable for ANY material law, because curl H = 0 with azimuthal H
forces H = C/r regardless of mu(B).  The whole cap-binding channel family is
therefore 1-D quadrature:

    H_cap = H(B_cap),   rho = r_out/r_in = H_cap / H(B_out)
    f(rho) = INT_1^rho B_of_H(H_cap/s) ds,     Phi  = r_in f(rho)
    body area = (Theta/2) r_in^2 (rho^2 - 1),  MMF  = r_in H_cap * angle

Consequences, stated honestly up front:
  (a) the hodograph design with BOTH walls at constant B has a theta-
      independent solution, i.e. it IS this annulus -- shape freedom adds
      nothing for the pure symmetric turn.  What the run buys is a machinery
      lock against an EXACT NONLINEAR reference (far stronger than the
      constant-mu sanity) plus the design chart engineers actually need;
  (b) the engineering payoff at this rung is aspect + sizing WITH the real
      saturating curve versus the linear picture (B ~ 1/r), quantified below;
  (c) genuinely non-annulus optima require distributed flux collection along
      the gap-side wall (the real SynRM channel) or asymmetric ends -- that
      is Rung 1.5, where no quadrature exists and the hodograph is the only
      linear-cost tool.

Golden bands asserted at the end of the run (2026-07-28 baseline, LAB):
  hodograph vs exact family : r_in, rho_fit, MMF rel err < 2e-3 per point;
                              body area rel err < 5e-3 (the residual 3.86e-3
                              is the angular quantization of the body mask --
                              the sampled body arc spans 1.5647 rad instead of
                              pi/2, exactly a 0.385% deficit -- not machinery)
  orientation               : J single-signed on every design
  FEM verify (rho=3, rho*)  : body inner wall vs the cap, mean < 1.0 %,
                              max < 2.0 %; MMF < 1.0 %; h/8 = h/16
  linear-designed annulus   : FEM peak matches the quadrature prediction

Measured at that baseline: machinery lock 2.6e-6..4.6e-5 on r_in / 2.8e-7..
1.5e-5 on rho / circularity ~1e-8; chart optimum rho* = 5.83 (5%-flat over
rho in [3.8, 9.5]); the linear-designed baseline carries the same flux with
peak 1.391 T (26.8% of the cap unused) and 2.67x the optimal iron
(sizing-with-the-real-curve x0.518, then aspect x0.722); the linear 1/r width
rule overestimates the needed channel width by +18% (rho=1.5) to +59% (rho=3).

Run:  python verify_synrm_channel_annulus_lock.py
Writes results_synrm_channel_annulus_lock.json next to this file (committed).
"""
import datetime
import json
import math
import os
import platform
import sys

import numpy as np
from scipy.optimize import brentq
from ngsolve import SetNumThreads, TaskManager

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_ipm_bridge_free_boundary as vb            # noqa: E402
from verify_ipm_bridge_free_boundary import (            # noqa: E402
    B_KNEE, DFLUX, MARGIN, THETA, TRIM, forward_solve, mu_s_of,
)

H_CAP = B_KNEE / mu_s_of(B_KNEE)
SPAN = THETA + 2.0 * MARGIN

# ---------------- material inversion (dense monotone table) ----------------
_BG = np.geomspace(1e-6, 6.0, 6000)
_HG = np.array([b / mu_s_of(b) for b in _BG])
if not np.all(np.diff(_HG) > 0):
    raise RuntimeError("H(B) is not monotone; the inversion table is invalid")


def B_of_H(H):
    return np.interp(H, _HG, _BG)


def f_of_rho(rho, n=4000):
    s = np.geomspace(1.0, float(rho), n)
    return float(np.trapezoid(B_of_H(H_CAP / s), s))


def family(rho):
    """Exact cap-binding annulus channel (any material law): H = C/r."""
    f = f_of_rho(rho)
    r_in = DFLUX / f
    return {
        "rho": float(rho),
        "B_out_T": float(B_of_H(H_CAP / rho)),
        "f": f,
        "r_in_mm": 1e3 * r_in,
        "width_mm": 1e3 * r_in * (rho - 1.0),
        "w_over_phi_per_T": (rho - 1.0) / f,
        "body_area_mm2": 1e6 * 0.5 * THETA * r_in ** 2 * (rho ** 2 - 1.0),
        "area_over_phi2": 0.5 * THETA * (rho ** 2 - 1.0) / f ** 2,
        "mmf_full_span_A": r_in * H_CAP * SPAN,
    }


# ---------------- helpers ----------------
def fit_circle(pts):
    A = np.c_[2 * pts[:, 0], 2 * pts[:, 1], np.ones(len(pts))]
    bb = (pts ** 2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, bb, rcond=None)
    cx, cy = sol[0], sol[1]
    r = math.sqrt(sol[2] + cx * cx + cy * cy)
    dev = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - r).max() / r
    return (cx, cy), r, dev


def hodo_point(rho, rep):
    """One hodograph design at constant walls (B_out(rho), B_cap)."""
    b_out = float(B_of_H(H_CAP / rho))
    old = vb.B_OUT0
    vb.B_OUT0 = b_out            # design(ramp=False) reads the module global
    try:
        inner, outer, inlet, outlet, ths, body = vb.design(
            rep, ramp=False, tag=f"rho{rho:.4g}")
    finally:
        vb.B_OUT0 = old
    d = rep["design"][f"rho{rho:.4g}"]
    _, ri, devi = fit_circle(inner)
    _, ro, devo = fit_circle(outer)
    ex = family(rho)
    out = {
        "rho": float(rho), "B_out_T": b_out,
        "r_in_fit_mm": 1e3 * ri, "rho_fit": ro / ri,
        "circ_dev": max(devi, devo),
        "r_in_rel_err": abs(1e3 * ri - ex["r_in_mm"]) / ex["r_in_mm"],
        "rho_rel_err": abs(ro / ri - rho) / rho,
        "area_rel_err": abs(d["body_iron_area_mm2"] - ex["body_area_mm2"])
        / ex["body_area_mm2"],
        "mmf_rel_err": abs(abs(d["mmf_design_A"]) - ex["mmf_full_span_A"])
        / ex["mmf_full_span_A"],
        "J_single_sign": d["J_single_sign"],
    }
    print(f"  [lock rho={rho:6.3f}] r_in {1e3*ri:9.5f} mm vs exact "
          f"{ex['r_in_mm']:9.5f} (rel {out['r_in_rel_err']:.2e})  "
          f"rho_fit rel {out['rho_rel_err']:.2e}  area rel "
          f"{out['area_rel_err']:.2e}  MMF rel {out['mmf_rel_err']:.2e}  "
          f"circ {out['circ_dev']:.2e}")
    return out, (inner, outer, inlet, outlet, ths, body)


def annulus_outline(r_in, rho, nsamp=175):
    """Directly-built annulus outline (for the linear-designed baseline)."""
    phis = np.linspace(-MARGIN, THETA + MARGIN, nsamp)
    ur = np.c_[np.cos(phis), np.sin(phis)]
    inner = r_in * ur
    outer = (rho * r_in) * ur
    rr = np.linspace(rho * r_in, r_in, 60)          # outer -> inner
    inlet = np.c_[rr * math.cos(phis[0]), rr * math.sin(phis[0])]
    outlet = np.c_[rr * math.cos(phis[-1]), rr * math.sin(phis[-1])]
    body = (phis >= -1e-12) & (phis <= THETA + 1e-12)
    return inner, outer, inlet, outlet, phis, body


def fem_verify(curves, r_in, rho, label, failures, peak_expect, mmf_expect,
               peak_band, hdiv=16.0):
    """Independent nonlinear FEM on an outline; check peak and MMF."""
    inner, outer, inlet, outlet, ths, body = curves
    w = r_in * (rho - 1.0)
    maxh = min(w / hdiv, r_in * SPAN / hdiv)
    res, bi, bo = forward_solve(inner, outer, inlet, outlet, maxh, label)
    body_c = body[TRIM:-TRIM]
    peak = float(bi[body_c].max())
    e_in = np.abs(bi[body_c] - peak_expect) / peak_expect
    mmf_rel = abs(res["mmf_fem_A"] - mmf_expect) / mmf_expect
    print(f"  [fem {label}] body inner |B| {bi[body_c].min():.3f}.."
          f"{peak:.3f} T (expect {peak_expect:.3f}; mean dev "
          f"{e_in.mean()*100:.3f}% max {e_in.max()*100:.3f}%)  "
          f"MMF {res['mmf_fem_A']:.4f} vs {mmf_expect:.4f} A "
          f"(rel {mmf_rel:.2e})  ne={res['n_elements']}")
    if e_in.mean() > peak_band[0]:
        failures.append(f"{label}: inner wall mean dev {e_in.mean():.4f} > "
                        f"{peak_band[0]}")
    if e_in.max() > peak_band[1]:
        failures.append(f"{label}: inner wall max dev {e_in.max():.4f} > "
                        f"{peak_band[1]}")
    if mmf_rel > 0.01:
        failures.append(f"{label}: MMF rel {mmf_rel:.4f} > 1.0%")
    return {"label": label, "n_elements": res["n_elements"],
            "picard_iterations": res["picard_iterations"],
            "inner_body_min_T": float(bi[body_c].min()),
            "inner_body_peak_T": peak,
            "inner_dev_mean": float(e_in.mean()),
            "inner_dev_max": float(e_in.max()),
            "mmf_fem_A": res["mmf_fem_A"], "mmf_expect_A": mmf_expect,
            "mmf_rel_err": mmf_rel}


def main():
    SetNumThreads(4)
    failures = []
    report = {"case": {
        "abstraction": "pure 90-degree turn, flux-line walls, equipotential "
                       "terminals, 20-degree lead-in/out",
        "B_cap_T": B_KNEE, "H_cap_A_per_m": H_CAP,
        "flux_Wb_per_m": DFLUX, "body_turn_deg": math.degrees(THETA),
        "material": "same representative curve as the promoted bridge driver "
                    "(NOT a datasheet fit; percentages are model-dependent)",
    }}

    # ---------- (1) exact chart ----------
    rhos = np.geomspace(1.05, 60.0, 400)
    aof = np.array([0.5 * THETA * (r ** 2 - 1.0) / f_of_rho(r) ** 2
                    for r in rhos])
    i0 = int(np.argmin(aof))
    rho_star = float(rhos[i0])
    fam_star = family(rho_star)
    flat = rhos[aof <= 1.05 * aof[i0]]
    print(f"(1) exact cap-binding family: area/Phi^2 minimum at rho* = "
          f"{rho_star:.3f} (value {aof[i0]:.4f}); within 5% of min over "
          f"rho in [{flat[0]:.2f}, {flat[-1]:.2f}]")
    chart_rows = [family(r) for r in
                  (1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, rho_star, 10.0, 20.0)]
    print("    rho   B_out   w/Phi[1/T]  area/Phi^2   r_in[mm]  w[mm]")
    for c in chart_rows:
        print(f"  {c['rho']:6.3f}  {c['B_out_T']:.3f}   {c['w_over_phi_per_T']:.4f}"
              f"      {c['area_over_phi2']:.4f}     {c['r_in_mm']:7.4f} "
              f"{c['width_mm']:7.4f}")
    report["chart"] = chart_rows
    report["optimum"] = {"rho_star": rho_star,
                         "area_over_phi2": float(aof[i0]),
                         "flat_5pct_rho_range": [float(flat[0]), float(flat[-1])]}

    # ---------- (2) linear-rule baseline ----------
    rho_lin = brentq(lambda r: r * r * math.log(r) - (r * r - 1.0), 1.3, 4.0)
    r_in_lin = DFLUX / (B_KNEE * math.log(rho_lin))
    rgrid = np.linspace(r_in_lin, rho_lin * r_in_lin, 4000)

    def phi_of_C(C):
        return float(np.trapezoid(B_of_H(C / rgrid), rgrid))

    C_lin = brentq(lambda C: phi_of_C(C) - DFLUX, 1e-8, 1e3)
    peak_lin = float(B_of_H(C_lin / r_in_lin))
    area_lin = 1e6 * 0.5 * THETA * r_in_lin ** 2 * (rho_lin ** 2 - 1.0)
    fam_capbound_lin = family(rho_lin)
    dec_sizing = area_lin / fam_capbound_lin["body_area_mm2"]
    dec_aspect = fam_capbound_lin["body_area_mm2"] / fam_star["body_area_mm2"]
    print(f"\n(2) linear-designed baseline: rho_lin = {rho_lin:.4f}, "
          f"r_in = {1e3*r_in_lin:.4f} mm, area = {area_lin:.4f} mm^2")
    print(f"    on that FIXED shape the real material carries Phi with peak "
          f"|B| = {peak_lin:.3f} T -> cap slack {B_KNEE - peak_lin:.3f} T "
          f"({100*(B_KNEE-peak_lin)/B_KNEE:.1f}% unused)")
    print(f"    decomposition at same Phi & cap: sizing-with-real-curve "
          f"x{1/dec_sizing:.3f} area, then aspect rho_lin->rho* "
          f"x{1/dec_aspect:.3f}; total linear-designed / optimal = "
          f"{area_lin/fam_star['body_area_mm2']:.3f}x")
    wcomp = []
    print("    width rule check   rho | w/Phi exact | w/Phi linear(1/r) | "
          "linear overestimates")
    for r in (1.5, 2.0, rho_lin, 3.0, 5.0):
        we = family(r)["w_over_phi_per_T"]
        wl = (r - 1.0) / (B_KNEE * math.log(r))
        wcomp.append({"rho": r, "w_over_phi_exact": we,
                      "w_over_phi_linear_rule": wl,
                      "linear_overestimate_pct": 100 * (wl / we - 1)})
        print(f"                    {r:5.3f} |    {we:.4f}   |      {wl:.4f}"
              f"       |   +{100*(wl/we-1):.1f} %")
    report["linear_baseline"] = {
        "rho_lin": rho_lin, "r_in_mm": 1e3 * r_in_lin,
        "area_mm2": area_lin, "peak_B_on_shape_T": peak_lin,
        "cap_slack_T": B_KNEE - peak_lin,
        "area_ratio_vs_capbound_same_rho": dec_sizing,
        "area_ratio_capbound_vs_optimal": dec_aspect,
        "area_ratio_vs_optimal": area_lin / fam_star["body_area_mm2"],
        "width_rule_rows": wcomp,
        "mmf_full_span_A": C_lin * SPAN,
    }

    with TaskManager():
        # ---------- (3) hodograph machinery lock vs the exact family ------
        print("\n(3) hodograph (const-B walls) vs the exact nonlinear annulus")
        rep = {}
        locks, curves_at = [], {}
        for rho in (1.5, rho_lin, 3.0, 5.0, rho_star, 10.0, 20.0):
            out, curves = hodo_point(rho, rep)
            locks.append(out)
            curves_at[round(rho, 4)] = curves
            if not out["J_single_sign"]:
                failures.append(f"rho={rho:.3f}: inverse map folds")
            for k, band in (("r_in_rel_err", 2e-3), ("rho_rel_err", 2e-3),
                            ("area_rel_err", 5e-3), ("mmf_rel_err", 2e-3)):
                if out[k] > band:
                    failures.append(f"rho={rho:.3f}: {k} {out[k]:.2e} > {band}")
        report["hodograph_lock"] = locks

        # ---------- (4) independent nonlinear FEM verifications -----------
        print("\n(4) independent nonlinear FEM on three outlines")
        fems = []
        fam3 = family(3.0)
        for hdiv, lab in ((8.0, "rho3_h8"), (16.0, "rho3_h16")):
            fems.append(fem_verify(curves_at[3.0], 1e-3 * fam3["r_in_mm"], 3.0,
                                   lab, failures, B_KNEE,
                                   fam3["mmf_full_span_A"], (0.010, 0.020),
                                   hdiv=hdiv))
        fems.append(fem_verify(curves_at[round(rho_star, 4)],
                               1e-3 * fam_star["r_in_mm"], rho_star,
                               f"rhostar{rho_star:.3g}_h16", failures, B_KNEE,
                               fam_star["mmf_full_span_A"], (0.010, 0.020)))
        # linear-designed annulus: expect its own (slack) peak, not the cap
        lin_curves = annulus_outline(r_in_lin, rho_lin)
        fems.append(fem_verify(lin_curves, r_in_lin, rho_lin, "linear_h16",
                               failures, peak_lin, C_lin * SPAN,
                               (0.015, 0.030)))
        report["fem_verify"] = fems

    report["verdict"] = {
        "pure_turn_optimum_is_annulus": True,
        "rho_star": rho_star,
        "iron_linear_designed_over_optimal": area_lin / fam_star["body_area_mm2"],
        "cap_slack_of_linear_design_T": B_KNEE - peak_lin,
        "note": "shape freedom adds nothing for the pure symmetric turn; "
                "the payoff here is designing aspect+size with the real "
                "curve; non-annulus optima need distributed collection "
                "(Rung 1.5)",
    }
    report["meta"] = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "hostname": platform.node(), "python_version": platform.python_version(),
        "purpose": "correctness validation only (no timing claims)",
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results_synrm_channel_annulus_lock.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"\nresults -> {out}")
    if failures:
        for f_ in failures:
            print("CHECK FAIL:", f_)
        raise SystemExit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
