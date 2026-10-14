"""Leaf-stacking + perturbation diagnostic for straight accelerator dipoles:
WHEN does "foliate the 3-D magnet into 2-D leaves + connect by a perturbation"
actually land?  Answer: it is governed by the aspect ratio L/gap, and the
inter-leaf coupling scales as ~ gap/L.

THE IDEA (the "foliate-and-perturb" program for quasi-2-D magnets).  A beamline
magnet varies slowly along the beam; only the ENDS are genuinely 3-D.  So:

    0th order : STACK the 2-D cross-section ("leaf") solution along the beam;
    1st order : CONNECT adjacent leaves by a beam-direction perturbation.

This pays off only if the inter-leaf coupling is small AND localised to the
ends.  This file MEASURES that on real finite-length C-frame dipoles solved by
the reduced-Omega + CoilBuilder forward engine (beam = y, gap = z, width = x;
geometry adapted from accel_pole_ends_fem.py but parametrised by the iron
length L so the ASPECT RATIO L/gap can be swept).

WHAT IS MEASURED (no separate 2-D solve needed -- for a STRAIGHT, constant-gap
magnet the BODY slice y=0 IS the 2-D infinite-long leaf, so the body field is
the 0th-order leaf-stack):

  * delta(y) = || B_perp(.,y) - B_perp(.,body) || / || B_perp(.,body) ||
        the 0th-ORDER LEAF-STACKING ERROR over a good-field aperture grid
        (~0 in the body, grows at the ends).
  * eps(y) = (g/2) |dBz/dy| / |Bz_body|
        the local PERTURBATION PARAMETER = (transverse scale)/(beam-variation
        scale).  The genuine smallness parameter of the leaf coupling (NOT an
        operator-norm ratio, which is trapped at 1 by grad_perp^2 = -d^2/dy^2
        in current-free air).
  * fringe_excess = (L_eff - L_iron)/L_iron
        the integrated 1st-order correction: how much the fringe (the inter-leaf
        coupling) adds to the 0th-order leaf-stack integral Bz_body * L_iron.

THE RESULT (aspect_sweep): the fringe / leaf-coupling DECAYS as ~ gap/L (a
log-log slope near -1).  So:
  - a COMPACT magnet (the lab test dipole, L/gap = 3) is firmly NON-perturbative
    (fringe ~ +69 %, the 0th-order stack is off by ~40 %, the 3-D-ness is NOT
    end-localised);
  - foliate-and-perturb only lands for LONG magnets (L/gap >> 1, few-% fringe).
The scaling law is the transferable answer: it tells you, for any L/gap, whether
the body-2-D + end-perturbation scheme is worth it for your magnet.

run:  python leaf_coupling_perturbation_3d.py              # single magnet (L=L_iron)
      python leaf_coupling_perturbation_3d.py --sweep       # aspect-ratio scaling
      python leaf_coupling_perturbation_3d.py --sweep --fig # + figure
"""

from _validation_output import validation_output
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from accel_pole_ends_fem import (GAP, POLE_W, Z_OUT, T_LEG, L_BEAM, AIR,  # noqa: E402
                                 COIL_W, COIL_H, COIL_R, COIL_STRAIGHT, NI, MU0)


def _air_half(L):
    """Air-box half-size that comfortably contains the iron + coil + fringe."""
    return max(AIR, 0.5 * L + 4 * GAP + 0.05)


def _build_mesh_L(L, maxh_air=0.05, maxh_iron=0.025):
    """H-frame (window) dipole of iron length L, x-symmetric, netgen.occ (no
    Cubit).  Adapted from accel_pole_ends_fem.build_mesh, parametrised by L so
    the aspect ratio L/gap can be swept."""
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    a_half = _air_half(L)
    hL = L / 2
    top = Box(Pnt(-POLE_W, -hL, GAP / 2), Pnt(POLE_W, hL, Z_OUT))
    bot = Box(Pnt(-POLE_W, -hL, -Z_OUT), Pnt(POLE_W, hL, -GAP / 2))
    leg_l = Box(Pnt(-POLE_W, -hL, -GAP / 2), Pnt(-POLE_W + T_LEG, hL, GAP / 2))
    leg_r = Box(Pnt(POLE_W - T_LEG, -hL, -GAP / 2), Pnt(POLE_W, hL, GAP / 2))
    iron = (top + bot + leg_l + leg_r)
    iron.mat("iron")
    iron.maxh = maxh_iron
    air = Box(Pnt(-a_half, -a_half, -a_half), Pnt(a_half, a_half, a_half)) - iron
    air.mat("air")
    for f in air.faces:
        c = f.center
        if max(abs(c.x), abs(c.y), abs(c.z)) > 0.9 * a_half:
            f.name = "outer"
    shape = Glue([air, iron])
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh_air))
    return mesh


def _build_coil_L(L):
    """Racetrack coil pair (CoilBuilder) -> Radia, with the straight section made
    longer than the iron (end-turns OUTSIDE the iron ends so the body sees a
    uniform source).  Adapted from accel_pole_ends_fem.build_coil."""
    import radia as rad
    from radia.coil_builder import CoilBuilder
    rad.UtiDelAll()
    straight = max(COIL_STRAIGHT, 1.6 * L)
    z_coil = GAP / 2 + COIL_H / 2
    upper = (CoilBuilder(current=NI)
             .set_start([COIL_R, -straight / 2, z_coil])
             .set_cross_section(width=COIL_W, height=COIL_H)
             .add_straight(straight)
             .add_arc(radius=COIL_R, arc_angle=180)
             .add_straight(straight)
             .add_arc(radius=COIL_R, arc_angle=180))
    lower = upper.mirror("xy")
    return rad.ObjCnt(upper.to_radia() + lower.to_radia())


def _solve_field(L, mu_r=1000.0, order=2, maxh_air=0.05, maxh_iron=0.025):
    """Reduced-Omega forward solve of the length-L finite dipole; returns
    (mesh, B, rad) with B = mu (H_s - grad Omega).  RadiaField source LinearForm
    + field readout are SERIAL (Radia is not thread-safe); only the stiffness
    solve is wrapped in TaskManager (the accel_pole_ends_fem.solve pattern)."""
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)
    import radia as rad
    mesh = _build_mesh_L(L, maxh_air, maxh_iron)
    coils = _build_coil_L(L)
    Hs = rad.RadiaField(coils, "h")
    mu = mesh.MaterialCF({"iron": mu_r * MU0}, default=MU0)
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    f = LinearForm(fes)
    f += mu * Hs * grad(v) * dx
    f.Assemble()                                     # serial (RadiaField source)
    with TaskManager():
        a = BilinearForm(fes)
        a += mu * grad(u) * grad(v) * dx
        a.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                     inverse="sparsecholesky") * f.vec
    B = mu * (Hs - grad(gfu))
    return mesh, B, rad


def run(L=L_BEAM, mu_r=1000.0, order=2, maxh_air=0.05, maxh_iron=0.025,
        r_ap=0.008, n_ap=5, n_y=81, delta_thresh=0.05):
    """Leaf-stacking diagnostic for the dipole of iron length L: delta(y), the
    perturbation parameter eps(y), and the integrated fringe excess; returns the
    per-magnet record (with the L/gap aspect ratio)."""
    mesh, B, rad = _solve_field(L, mu_r, order, maxh_air, maxh_iron)
    y_max = min(0.5 * L + 5 * GAP, 0.95 * _air_half(L))
    ys = np.linspace(-y_max, y_max, n_y)
    xs = np.linspace(-r_ap, r_ap, n_ap)
    zs = np.linspace(-r_ap, r_ap, n_ap)

    def bz_grid(yv):                                  # dipole field over the aperture leaf
        return np.array([float(B(mesh(float(xx), float(yv), float(zz)))[2])
                         for xx in xs for zz in zs])

    g_body = bz_grid(0.0)                             # body slice = the 2-D leaf
    n_body = float(np.linalg.norm(g_body))
    delta = np.array([np.linalg.norm(bz_grid(yv) - g_body) / n_body for yv in ys])
    bz_axis = np.array([float(B(mesh(0.0, float(yv), 0.0))[2]) for yv in ys])
    rad.UtiDelAll()

    body = np.abs(ys) < 0.30 * L
    bz_body = float(np.mean(bz_axis[body]))
    dbz = np.gradient(bz_axis, ys)
    eps = (GAP / 2.0) * np.abs(dbz) / abs(bz_body)

    L_iron = float(L)
    bbar1_true = float(np.trapezoid(bz_axis, ys))
    L_eff = bbar1_true / bz_body
    fringe_excess = (L_eff - L_iron) / L_iron
    iron = np.abs(ys) < 0.5 * L_iron
    return {
        "L_iron_m": L_iron, "aspect_L_over_gap": float(L / GAP),
        "mu_r": float(mu_r), "order": int(order), "ne": int(mesh.ne),
        "ys": ys.tolist(), "delta": delta.tolist(), "eps": eps.tolist(),
        "bz_body_T": bz_body,
        "delta_body_max": float(np.max(delta[body])),
        "eps_body_max": float(np.max(eps[body])),
        "delta_farend_max": float(np.max(delta[np.abs(ys) > 0.85 * y_max])),
        "end_localisation_frac": float(np.mean(delta[iron] > delta_thresh)),
        "L_eff_m": L_eff, "fringe_excess": float(fringe_excess),
        "bbar1_0th_leafstack_Tm": bz_body * L_iron, "bbar1_true_Tm": bbar1_true,
    }


def aspect_sweep(Ls=(0.080, 0.120, 0.200, 0.320), mu_r=1000.0, order=2,
                 maxh_air=0.055, maxh_iron=0.028):
    """Sweep the iron length -> aspect ratio L/gap, and establish the leaf-
    coupling SCALING.  The fringe excess (inter-leaf perturbation) decays as
    ~ (gap/L)^p; a log-log fit gives the exponent p (~1)."""
    rows = [run(L=L, mu_r=mu_r, order=order, maxh_air=maxh_air,
                maxh_iron=maxh_iron) for L in Ls]
    ar = np.array([r["aspect_L_over_gap"] for r in rows])
    fr = np.array([r["fringe_excess"] for r in rows])
    db = np.array([r["delta_body_max"] for r in rows])
    # log-log slope of fringe_excess vs L/gap (expect ~ -1: fringe ~ gap/L)
    slope = float(np.polyfit(np.log(ar), np.log(fr), 1)[0])
    # extrapolate the aspect ratio where the fringe drops to 10% (a typical
    # "good-enough perturbation" tolerance) from the power-law fit.
    c = float(np.exp(np.polyfit(np.log(ar), np.log(fr), 1)[1]))
    ar_10 = float((0.10 / c) ** (1.0 / slope)) if slope != 0 else float("nan")
    return {
        "aspect_L_over_gap": ar.tolist(),
        "fringe_excess": fr.tolist(),
        "delta_body_max": db.tolist(),
        "fringe_scaling_exponent": slope,         # ~ -1 => fringe ~ gap/L
        "aspect_for_10pct_fringe": ar_10,         # L/gap needed for a 10% fringe
        "monotone_decay": bool(np.all(np.diff(fr) < 0)),
        "rows": rows,
    }


def _figure_single(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ys = np.array(res["ys"]) / (0.5 * res["L_iron_m"])
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for a in ax:
        a.axvspan(-0.6, 0.6, color="0.85")
        a.axvline(1.0, color="r", ls="--", lw=1)
        a.axvline(-1.0, color="r", ls="--", lw=1)
        a.set_xlabel("y / (L_iron/2)")
    ax[0].plot(ys, res["delta"], "o-", ms=3)
    ax[0].axhline(0.05, color="k", ls=":", lw=1)
    ax[0].set_ylabel("leaf-stacking error delta(y)")
    ax[0].set_title(f"L/gap={res['aspect_L_over_gap']:.1f}: body "
                    f"{res['delta_body_max']:.1e}, end-loc "
                    f"{res['end_localisation_frac']*100:.0f}%")
    ax[1].plot(ys, res["eps"], "s-", ms=3, color="C1")
    ax[1].set_ylabel("perturbation eps(y)=(g/2)|dBz/dy|/Bz")
    ax[1].set_title(f"fringe excess {res['fringe_excess']*100:+.0f}%")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "leaf_coupling_perturbation_3d.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


def _figure_sweep(sw):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ar = np.array(sw["aspect_L_over_gap"])
    fr = np.array(sw["fringe_excess"])
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(ar, fr, "o-", ms=5, label="measured fringe excess")
    ax[0].loglog(ar, fr[0] * (ar / ar[0]) ** (-1.0), "k--", lw=1,
                 label="~ (gap/L)^1 (slope -1)")
    ax[0].axhline(0.10, color="r", ls=":", lw=1, label="10% (perturbation lands)")
    ax[0].set_xlabel("aspect ratio  L_iron / gap")
    ax[0].set_ylabel("fringe excess (L_eff-L_iron)/L_iron")
    ax[0].set_title(f"leaf coupling ~ gap/L  (fit slope {sw['fringe_scaling_exponent']:.2f};"
                    f" 10% at L/gap~{sw['aspect_for_10pct_fringe']:.0f})")
    ax[0].legend(fontsize=8)
    ax[1].semilogx(ar, np.array(sw["delta_body_max"]) * 100, "s-", ms=5)
    ax[1].set_xlabel("aspect ratio  L_iron / gap")
    ax[1].set_ylabel("body leaf-stacking error (%)")
    ax[1].set_title("body 2-D-exactness improves with length")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "leaf_coupling_perturbation_3d_sweep.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", action="store_true")
    ap.add_argument("--sweep", action="store_true",
                    help="aspect-ratio (L/gap) scaling of the leaf coupling")
    ap.add_argument("--order", type=int, default=2)
    args = ap.parse_args()

    print("=" * 74)
    print("Leaf-stacking + perturbation diagnostic -- straight C-frame dipole")
    print("=" * 74)
    res = run(order=args.order)
    print(f"aspect ratio L/gap          : {res['aspect_L_over_gap']:.2f}"
          f"   (iron {res['L_iron_m']*1e3:.0f} mm, gap {GAP*1e3:.0f} mm)")
    print(f"body flat-top field Bz      : {res['bz_body_T']*1e3:.2f} mT")
    print(f"leaf-stack error  body max  : {res['delta_body_max']:.2e}")
    print(f"                  far-end max: {res['delta_farend_max']:.2e}"
          f"   (the ends ARE 3-D)")
    print(f"perturbation eps  body max  : {res['eps_body_max']:.2e}")
    print(f"3-D-ness end-localisation   : {res['end_localisation_frac']*100:.0f}%"
          f" of the iron span > 5%")
    print(f"fringe excess (1st-order)   : {res['fringe_excess']*100:+.1f}%"
          f"   (0th-stack {res['bbar1_0th_leafstack_Tm']*1e3:.2f} vs true "
          f"{res['bbar1_true_Tm']*1e3:.2f} mT.m)")
    here = os.path.dirname(os.path.abspath(__file__))
    with validation_output("leaf_coupling_perturbation_3d.json").open("w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("saved leaf_coupling_perturbation_3d.json")
    if args.fig:
        _figure_single(res)

    if args.sweep:
        print("\n" + "=" * 74)
        print("ASPECT-RATIO SWEEP -- leaf coupling vs L/gap")
        print("=" * 74)
        sw = aspect_sweep(order=args.order)
        with validation_output("leaf_coupling_perturbation_3d_sweep.json").open("w", encoding="utf-8") as fh:
            json.dump(sw, fh, indent=2)
        print("saved leaf_coupling_perturbation_3d_sweep.json")
        print(f"{'L/gap':<8}{'fringe excess':<16}{'body delta':<12}")
        for a, fr, db in zip(sw["aspect_L_over_gap"], sw["fringe_excess"],
                             sw["delta_body_max"]):
            print(f"{a:<8.1f}{fr*100:<+16.1f}{db*100:<12.2f}")
        print("-" * 74)
        print(f"fringe scaling exponent     : {sw['fringe_scaling_exponent']:.2f}"
              f"   (~ -1 => leaf coupling ~ gap/L)")
        print(f"L/gap for a 10% fringe      : ~{sw['aspect_for_10pct_fringe']:.0f}"
              f"   (foliate-and-perturb lands above this)")
        print(f"VERDICT: compact magnets (L/gap~3) are NON-perturbative;"
              f" the scheme lands for long (L/gap>>1) dipoles.")
        if args.fig:
            _figure_sweep(sw)


if __name__ == "__main__":
    main()
