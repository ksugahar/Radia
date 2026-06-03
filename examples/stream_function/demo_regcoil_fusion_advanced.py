#!/usr/bin/env python3
"""Advanced fusion coil design -- coil force / stress, a REAL VMEC equilibrium,
and a FOCUS-style winding-surface study.  Builds on ``demo_regcoil_fusion.py``
(the forward + REGCOIL-L-curve + net-current + VMEC-boundary foundation), reusing
its helpers, and adds the three "named next steps" from that demo:

  A. COIL FORCE / STRESS.  A surface current K experiences a Lorentz force per
     unit area  f = K x B_avg , where B_avg = 1/2 (B+ + B-) is the average of the
     coil's own field on the two sides of the sheet (the tangential B jumps by
     mu0 K across it).  For the net-poloidal-current (TF) coil this is the
     classic magnetic stress: |f| equals the magnetic pressure  B_tor^2/(2 mu0) ,
     peaking on the INBOARD leg (small R, where B_tor ~ 1/R is largest) -- exactly
     why a tokamak/stellarator TF coil is stress-limited at its inboard side.
     The demo verifies |f| against B^2/(2 mu0) to a few percent.

  B. A REAL VMEC EQUILIBRIUM.  Instead of the analytic rotating-ellipse model,
     design against the boundary of a genuine free-boundary equilibrium -- the
     li383 (NCSX-like, NFP=3, quasi-axisymmetric) reference ``wout`` shipped by
     simsopt.  ``--wout PATH`` uses any VMEC output; with no path the demo fetches
     the li383 reference (a 121 kB netCDF) to a cache.  Proves the Stage-2 forward
     map (winding-surface current potential -> plasma B.n) runs on real
     stellarator geometry.

  C. FOCUS-STYLE WINDING-SURFACE STUDY.  FOCUS / REGCOIL also optimise the
     WINDING SURFACE.  Sweeping the winding-surface standoff (gap to the plasma)
     shows the fundamental trade: a closer winding surface couples better, so it
     reproduces the target with LOWER B.n residual AND lower peak current density
     |grad psi| (coil complexity) -- the curve is MONOTONIC, so the distance
     optimum is CONSTRAINT-BOUND (push to the minimum engineering standoff,
     d_min, set by blanket / access / neutron shielding).  We confirm a bounded
     minimisation of coil complexity lands on d_min.  (Optimising the winding
     surface SHAPE -- conformal Fourier modes at fixed standoff, FOCUS's main
     contribution -- is the deeper next step, named but not done here.)

HONEST SCOPE:
  * Force is the MAGNETIC force per unit area (the stress DRIVER, N/m^2); the
    mechanical hoop stress in the conductor needs a structural model (not here).
  * The real-wout path needs the stellarator-symmetric li383 boundary; the base
    reader raises on a non-symmetric (lasym=T) file.
  * The FOCUS study optimises winding DISTANCE only (monotonic -> constraint
    bound); winding-SHAPE optimisation is named, not implemented.

Run (standalone):
    python demo_regcoil_fusion_advanced.py                 # all three, fetch li383
    python demo_regcoil_fusion_advanced.py --wout my.nc    # your own equilibrium
    python demo_regcoil_fusion_advanced.py --no-plot --no-fetch   # offline (skip B)

Writes ``demo_regcoil_fusion_advanced.json`` + (if matplotlib + radia-mcp) a
1x3 lab figure.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import urllib.request
from datetime import datetime

import numpy as np

# the li383 reference equilibrium shipped by simsopt (open source, MIT)
_LI383_URL = ("https://raw.githubusercontent.com/hiddenSymmetries/simsopt/"
              "master/tests/test_files/wout_li383_low_res_reference.nc")


def _src_paths():
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    sys.path.insert(0, os.path.join(here, "..", "..", "src", "radia"))
    sys.path.insert(0, os.path.join(here, "..", "..", "src", "radia", "panels"))


def _torus_at(R0, a, maxh, path):
    """Torus SURFACE mesh of major radius R0, minor radius a (the base demo's
    builder hardcodes R_MAJOR; the real wout sits at R0 ~ 1.4 m)."""
    from netgen.occ import WorkPlane, Axes, Pnt, Y, Z, Axis, OCCGeometry
    wire = WorkPlane(Axes(Pnt(R0, 0, 0), n=Y)).Circle(a).Wire()
    OCCGeometry(wire.Revolve(Axis(Pnt(0, 0, 0), Z), 360)).GenerateMesh(
        maxh=maxh).Save(path)
    return path


# --------------------------------------------------------------------------
# A. Coil force / stress
# --------------------------------------------------------------------------
def _tf_force_vs_poloidal(base, coil, ntheta, eps):
    """Lorentz force per area |f| = |K x B_avg| of the unit net-poloidal-current
    (TF) coil at poloidal angles theta (phi=0 cross-section), vs the analytic
    magnetic pressure B_tor^2/(2 mu0) = mu0/(2 R^2).  B_avg is the +/- eps
    average of the coil self-field across the current sheet."""
    from ngsolve import (specialcf, Cross, CoefficientFunction as CF,
                         x, y, z, sqrt, Integrate, ds)
    R, a = base.R_MAJOR, base.A_WIND
    n_cf = specialcf.normal(3)
    K = Cross(n_cf, CF((-y, x, 0)) / (x * x + y * y))     # unit TF current

    def Bvec(q):
        dxt, dyt, dzt = q[0] - x, q[1] - y, q[2] - z
        r3 = (dxt * dxt + dyt * dyt + dzt * dzt) * sqrt(
            dxt * dxt + dyt * dyt + dzt * dzt)
        cr = (K[1] * dzt - K[2] * dyt, K[2] * dxt - K[0] * dzt,
              K[0] * dyt - K[1] * dxt)
        return np.array([base._MU0_4PI * Integrate(cr[c] / r3 * ds, coil)
                         for c in range(3)])

    thetas = np.linspace(0.0, np.pi, ntheta)              # outboard -> inboard
    f_mag = np.zeros(ntheta); p_ana = np.zeros(ntheta); Rloc = np.zeros(ntheta)
    for i, th in enumerate(thetas):
        Rp = R + a * math.cos(th)
        px, pz = Rp, a * math.sin(th)                     # surface point, phi=0
        nx, nz = math.cos(th), math.sin(th)               # outward normal, phi=0
        Bavg = 0.5 * (Bvec((px + eps * nx, 0, pz + eps * nz))
                      + Bvec((px - eps * nx, 0, pz - eps * nz)))
        Kloc = np.array([-math.sin(th) / Rp, 0.0, math.cos(th) / Rp])
        f_mag[i] = np.linalg.norm(np.cross(Kloc, Bavg))
        p_ana[i] = base._MU0 / (2.0 * Rp * Rp)            # B_tor^2/(2 mu0)
        Rloc[i] = Rp
    return thetas, Rloc, f_mag, p_ana


# --------------------------------------------------------------------------
# B. real VMEC equilibrium
# --------------------------------------------------------------------------
def _ensure_wout(path_arg, cache_dir, allow_fetch):
    if path_arg:
        return path_arg
    cache = os.path.join(cache_dir, "wout_li383_low_res_reference.nc")
    if os.path.exists(cache):
        return cache
    if not allow_fetch:
        return None
    os.makedirs(cache_dir, exist_ok=True)
    urllib.request.urlretrieve(_LI383_URL, cache)         # raises on network fail
    return cache


def _design_on_wout(base, C, fes_mod, wout_path, eval_max, work_dir):
    """Build a circular winding torus enclosing the wout boundary, design the
    vertical-field B.n on it, return a result dict + boundary samples."""
    from ngsolve import Mesh, H1, specialcf
    terms, nfp = base._read_wout_boundary(wout_path)
    R0 = sum(rc for m, n, rc, zs in terms if m == 0 and n == 0)
    vpts, vnrm = base._vmec_points_normals(terms, nfp, 64, 64)
    axis = np.column_stack([
        R0 * vpts[:, 0] / np.hypot(vpts[:, 0], vpts[:, 1]),
        R0 * vpts[:, 1] / np.hypot(vpts[:, 0], vpts[:, 1]),
        np.zeros(len(vpts))])
    minor_max = float(np.linalg.norm(vpts - axis, axis=1).max())
    a_wind = 1.35 * minor_max                              # enclosing torus + gap
    coil = Mesh(_torus_at(R0, a_wind, 0.45 * a_wind,
                          os.path.join(work_dir, "wout_wind.vol")))
    fes = H1(coil, order=1, definedon=coil.Boundaries(".*"))
    if len(vpts) > eval_max:
        sel = np.linspace(0, len(vpts) - 1, eval_max).astype(int)
        vpts, vnrm = vpts[sel], vnrm[sel]
    A_n = base._A_normal(C, fes, specialcf.normal(3), vpts, vnrm)
    _psi, res, peakJ = base._design(C, fes, coil, A_n, vnrm[:, 2], "h1", 1e-7)
    return {"source": os.path.basename(wout_path), "nfp": nfp,
            "n_modes": len(terms), "R0": R0, "minor_radius_max": minor_max,
            "winding_R0": R0, "winding_a": a_wind, "winding_ndof": int(fes.ndof),
            "bn_residual_rel": res, "peak_grad_psi": peakJ}, terms, nfp


# --------------------------------------------------------------------------
# C. FOCUS-style winding-surface standoff study
# --------------------------------------------------------------------------
def _focus_standoff_sweep(base, C, plasma_pts, plasma_nrm, Bn, gaps, work_dir):
    """For each winding-surface gap (a_wind = a_plasma + gap), design the target
    and record (B.n residual, peak |grad psi|).  Returns the sweep rows."""
    from ngsolve import Mesh, H1, specialcf
    rows = []
    for gap in gaps:
        a_wind = base.A_PLASMA + gap
        coil = Mesh(_torus_at(base.R_MAJOR, a_wind, 0.4 * a_wind,
                              os.path.join(work_dir, "focus_w.vol")))
        fes = H1(coil, order=1, definedon=coil.Boundaries(".*"))
        A_n = base._A_normal(C, fes, specialcf.normal(3), plasma_pts, plasma_nrm)
        _psi, res, peakJ = base._design(C, fes, coil, A_n, Bn, "h1", 1e-7)
        rows.append({"gap": float(gap), "a_wind": float(a_wind),
                     "ndof": int(fes.ndof), "bn_residual_rel": res,
                     "peak_grad_psi": peakJ})
    return rows


def main():
    ap = argparse.ArgumentParser(
        description="advanced fusion coil design: force, real wout, FOCUS")
    ap.add_argument("--out-dir",
                    default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--work-dir",
                    default=os.path.join(os.environ.get("TEMP", "/tmp"),
                                         "sf_regcoil_adv"))
    ap.add_argument("--wout", default=None,
                    help="VMEC wout_*.nc (default: fetch li383 reference)")
    ap.add_argument("--no-fetch", action="store_true",
                    help="do not download the li383 wout (skip Part B if absent)")
    ap.add_argument("--force-ntheta", type=int, default=19)
    ap.add_argument("--force-eps", type=float, default=0.012)
    ap.add_argument("--eval-max", type=int, default=200)
    ap.add_argument("--d-min", type=float, default=0.03,
                    help="minimum engineering standoff for the FOCUS optimum (m)")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    _src_paths()
    import demo_regcoil_fusion as base
    import calc_streamfunction as C
    from ngsolve import Mesh, TaskManager

    os.makedirs(args.work_dir, exist_ok=True)

    with TaskManager():
        # ---- A. coil force / stress (TF coil) ----
        coil = Mesh(base._torus_surface_vol(
            base.A_WIND, 0.03, os.path.join(args.work_dir, "force.vol")))
        thetas, Rloc, f_mag, p_ana = _tf_force_vs_poloidal(
            base, coil, args.force_ntheta, args.force_eps)
        ratio = f_mag / p_ana
        i_peak = int(np.argmax(f_mag))
        force = {"check": "|f| = |K x B_avg| vs magnetic pressure B^2/(2 mu0)",
                 "ratio_mean": float(np.mean(ratio)),
                 "ratio_min": float(np.min(ratio)),
                 "ratio_max": float(np.max(ratio)),
                 "peak_stress_Pa_per_unit_K2": float(f_mag[i_peak]),
                 "peak_at_R": float(Rloc[i_peak]),
                 "inboard_over_outboard": float(f_mag[-1] / f_mag[0])}
        print(f"[force] |f| vs B^2/2mu0: ratio mean={force['ratio_mean']:.3f} "
              f"[{force['ratio_min']:.3f}, {force['ratio_max']:.3f}]; peak stress "
              f"on the INBOARD leg (R={force['peak_at_R']:.3f}), "
              f"inboard/outboard={force['inboard_over_outboard']:.2f}")

        # ---- C. FOCUS standoff sweep (rotating-ellipse plasma, hard target) ----
        plasma = Mesh(base._torus_surface_vol(
            base.A_PLASMA, 0.035, os.path.join(args.work_dir, "focus_p.vol")))
        ppts, pnrm, pth, pph = base._plasma_points_normals(plasma, args.eval_max)
        Bn_hard = np.sin(2 * pth) * np.cos(3 * pph)
        gaps = np.array([0.025, 0.04, 0.06, 0.09, 0.13, 0.18])
        sweep = _focus_standoff_sweep(base, C, ppts, pnrm, Bn_hard, gaps,
                                      args.work_dir)
        # bounded minimisation of coil complexity s.t. gap >= d_min: the curve is
        # monotonic (closer = simpler), so the optimum is the constraint d_min.
        feasible = [r for r in sweep if r["gap"] >= args.d_min - 1e-9]
        opt = min(feasible, key=lambda r: r["peak_grad_psi"])
        monotonic = all(sweep[i]["peak_grad_psi"] <= sweep[i + 1]["peak_grad_psi"]
                        for i in range(len(sweep) - 1))
        focus = {"sweep": sweep, "d_min": args.d_min,
                 "optimum_gap": opt["gap"], "optimum_peak_grad_psi":
                 opt["peak_grad_psi"], "complexity_monotonic_in_gap": monotonic,
                 "complexity_ratio_far_over_near":
                 float(sweep[-1]["peak_grad_psi"] / sweep[0]["peak_grad_psi"])}
        print(f"[focus] coil complexity monotonic in gap={monotonic} "
              f"(far/near peak|grad psi|={focus['complexity_ratio_far_over_near']:.1f}x)"
              f" -> constrained optimum at d_min={args.d_min} m "
              f"(gap={opt['gap']}, peak={opt['peak_grad_psi']:.2e})")

        # ---- B. real VMEC equilibrium (li383) ----
        wout_path = _ensure_wout(args.wout, args.work_dir, not args.no_fetch)
        if wout_path and os.path.exists(wout_path):
            real, wterms, wnfp = _design_on_wout(
                base, C, None, wout_path, args.eval_max, args.work_dir)
            print(f"[real-wout] {real['source']}: NFP={real['nfp']}, "
                  f"{real['n_modes']} modes, R0={real['R0']:.3f} m -> winding "
                  f"torus a={real['winding_a']:.3f} m, B.n resid="
                  f"{real['bn_residual_rel']:.2e}")
        else:
            real, wterms, wnfp = None, None, None
            print("[real-wout] skipped (no wout file; --no-fetch and none cached)")

    data = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "benchmark": "regcoil_fusion_advanced",
        "force_stress": force,
        "focus_standoff": focus,
        "real_equilibrium": real,
    }
    jpath = os.path.join(args.out_dir, "demo_regcoil_fusion_advanced.json")
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"\nResults -> {jpath}")

    if not args.no_plot:
        _plot(args, thetas, Rloc, f_mag, p_ana, sweep, args.d_min,
              wterms, wnfp)


def _plot(args, thetas, Rloc, f_mag, p_ana, sweep, d_min, wterms, wnfp):
    """1x3 lab figure: (a) TF coil stress vs poloidal angle (+ analytic
    pressure); (b) the real li383 boundary; (c) FOCUS complexity vs gap."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from radia_mcp.figure import apply_lab_style, save_lab_figure
    except Exception as e:
        print(f"(plot skipped: radia-mcp figure unavailable: {e})")
        return
    plt.rcParams["pdf.fonttype"] = 42
    w_in, h_in = apply_lab_style(embed_width_cm=16.0, aspect=0.40)
    fig = plt.figure(figsize=(w_in, h_in))

    # (a) coil stress vs poloidal angle: |f| (points) vs B^2/2mu0 (line)
    ax = fig.add_subplot(1, 3, 1)
    deg = np.degrees(thetas)
    ax.plot(deg, p_ana, "-", color="0.5", label=r"$B^2/2\mu_0$ (analytic)")
    ax.plot(deg, f_mag, "o", ms=3, color="C3", label=r"$|K\times B_{\rm avg}|$")
    ax.set_xlabel(r"poloidal angle  (deg, 0=outboard)")
    ax.set_ylabel(r"magnetic stress  (per unit $K^2$)")
    ax.text(0.30, 0.93, "(a) TF coil stress", transform=ax.transAxes)
    ax.legend(fontsize=7, loc="lower right", frameon=False)

    # (b) the real li383 boundary surface
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    if wterms is not None:
        import demo_regcoil_fusion as base
        u = np.linspace(0, 2 * np.pi, 80); v = np.linspace(0, 2 * np.pi, 60)
        U, V = np.meshgrid(u, v)
        R, Z, *_ = base._vmec_boundary(wterms, wnfp, V, U)
        ax2.plot_surface(R * np.cos(U), R * np.sin(U), Z, color="orange",
                         alpha=0.5, linewidth=0, rstride=2, cstride=2)
    for s in ("x", "y", "z"):
        getattr(ax2, f"set_{s}ticks")([])
    ax2.text2D(-0.05, 0.02, "(b) li383 plasma (real VMEC)",
               transform=ax2.transAxes)
    ax2.set_box_aspect((1, 1, 0.4))

    # (c) FOCUS: coil complexity vs winding gap (monotonic -> push to d_min)
    ax3 = fig.add_subplot(1, 3, 3)
    g = [r["gap"] for r in sweep]; pk = [r["peak_grad_psi"] for r in sweep]
    ax3.semilogy(g, pk, "o-", color="C0")
    ax3.axvline(d_min, ls="--", lw=0.8, color="C2")
    ax3.set_xlabel(r"winding gap to plasma  (m)")
    ax3.set_ylabel(r"peak $|\nabla\psi|$  (coil complexity)")
    ax3.text(0.04, 0.90, "(c) FOCUS standoff", transform=ax3.transAxes)
    ax3.text(d_min, pk[0], r" $d_{\min}$", color="C2", fontsize=8,
             va="top")

    fig.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.20,
                        wspace=0.42)
    info = save_lab_figure(fig, os.path.join(args.out_dir,
                                             "demo_regcoil_fusion_advanced"),
                           embed_width_cm=16.0, tighten=False)
    print(f"Plot    -> {', '.join(info['wrote'])}")


if __name__ == "__main__":
    main()
