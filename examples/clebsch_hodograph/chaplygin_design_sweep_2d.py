r"""Maximally exploiting "nonlinear analysis done LINEARLY": a saturable-magnet
DESIGN SWEEP whose every nonlinear operating point is ONE linear hodograph solve.

The companion `chaplygin_hodograph_2d.py` established the engine: the 2-D hodograph
turns nonlinear current-free magnetostatics `div(nu(|B|) grad A) = 0` into a LINEAR
variable-coefficient problem (Molenbroek-Chaplygin), and for a slender saturable flux
guide that linear problem collapses to a single quadrature -- the "1-shot" drive

    DeltaPhi(Psi) = INT_0^L  nu(|B|(x)) |B|(x) dx ,     |B|(x) = Psi / w(x) ,

which reproduces the full nonlinear FEM Picard loop (to the slenderness error).

This file EXPLOITS that: a whole DESIGN SPACE -- a grid of (throat width x operating
field) -- is computed by ONE linear quadrature PER POINT, mesh-free, in milliseconds,
where the equivalent nonlinear FEM would need a Picard loop (many linear solves on a
curved mesh) at every point.  So the entire nonlinear design map is obtained at LINEAR
cost.  Concretely:

  * design_map(...)  -- the drive DeltaPhi(w_throat, B_throat) over a grid: a family of
    drive-vs-field curves, one per throat width, each BENDING UP as the throat
    saturates (the design content).  ~M quadratures, milliseconds.

  * transfer_curve(...) -- the engineering payoff: invert the map to the flux
    REGULATOR / field-limiter transfer Psi(drive).  As the throat saturates the curve
    CLAMPS (a given extra MMF pushes ever less extra flux), and a NARROWER throat
    clamps at LOWER flux -- the throat width is the design knob for the clamp level.

  * validate_vs_fem(...) -- run the full nonlinear FEM Picard at a few sample points
    (reusing chaplygin_hodograph_2d.solve_chaplygin): the 1-shot agrees to ~1-2 %
    (slenderness-limited), and each FEM point costs N (~10-15) Picard iterations -- the
    cost the linear sweep avoids at EVERY map point.

  * cost_summary(...) -- M map points at 1 quadrature each vs M x N_iter FEM linear
    solves: the design map is ~N_iter x (Picard) x (mesh-solve / quadrature) cheaper.

Honest scope: the 1-shot is the SLENDER-guide (fixed simple hodograph image) limit;
its error vs the true nonlinear field is the slenderness error (~1 % here, tightening
as the guide is made more slender -- see chaplygin_hodograph_2d.slenderness_trend).
The design map is exact to that error.  Genuinely 2-D fields (the field TURNS, the
hodograph image is a real region) need the linear hodograph PDE solve, not a quadrature
(`chaplygin_turning_guide_2d.py`); those too are LINEAR (one solve per case), so the
same "design space at linear cost" exploit applies one rung up.

run:  python chaplygin_design_sweep_2d.py            # design map + clamp (fast)
      python chaplygin_design_sweep_2d.py --fem      # + nonlinear FEM validation (slow)
"""
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chaplygin_hodograph_2d as ch                 # noqa: E402  the engine + FEM loop

MU0 = ch.MU0
_trap = np.trapezoid


# --------------------------------------------------------------------------- #
# the mesh-free 1-shot drive (the linear hodograph quadrature)
# --------------------------------------------------------------------------- #
def _w_of_x(xc, L, H, depth, Rc):
    """Analytic local guide width w(x) of the circular-notch throat (vectorised,
    mesh-free) -- the same profile chaplygin_hodograph_2d's _throat_geom returns."""
    yc = H / 2.0 + Rc - depth
    half = math.sqrt(max(0.0, Rc * Rc - (Rc - depth) ** 2))
    xx = np.asarray(xc) - L / 2.0
    return np.where(np.abs(xx) >= half, H,
                    2.0 * (yc - np.sqrt(np.maximum(0.0, Rc * Rc - xx * xx))))


def one_shot_drive(B_throat, depth, L=0.20, H=0.04, Rc=0.11, mur0=200.0, Bk=1.0,
                   n_quad=2001):
    """The drive (magnetomotive force) DeltaPhi to hold the throat at field B_throat,
    by ONE linear hodograph quadrature -- no mesh, no Picard iteration."""
    w_throat = H - 2.0 * depth
    Psi = B_throat * w_throat                         # flux carried at this throat field
    xs = np.linspace(0.0, L, n_quad)
    w = _w_of_x(xs, L, H, depth, Rc)
    Bx = Psi / w                                      # |B|(x) = Psi / w(x)
    nu = 1.0 / (MU0 * ch._mu_r(Bx, mur0, Bk))         # reluctivity nu(|B|)
    return float(_trap(nu * Bx, xs))


# --------------------------------------------------------------------------- #
# the design map: drive over (throat width x throat field), all 1-shot
# --------------------------------------------------------------------------- #
def design_map(depths=(0.008, 0.012, 0.016), B_throats=None, H=0.04, **kw):
    """A grid of one-shot drives DeltaPhi(w_throat, B_throat).  Returns the map, the
    saturation BEND per width, and the wall-clock cost (milliseconds for the whole
    map)."""
    if B_throats is None:
        B_throats = np.linspace(0.2, 2.4, 12)
    B_throats = np.asarray(B_throats, float)
    t0 = time.perf_counter()
    M = np.array([[one_shot_drive(float(B), float(d), H=H, **kw) for B in B_throats]
                  for d in depths])
    dt = time.perf_counter() - t0
    widths = [H - 2.0 * d for d in depths]
    # saturation bend: reluctance-per-field at the top vs bottom of the field range.
    bend = [(M[i, -1] / B_throats[-1]) / (M[i, 0] / B_throats[0]) for i in range(len(depths))]
    return {
        "depths": list(depths), "widths": widths, "B_throats": B_throats.tolist(),
        "drive": M.tolist(), "bend_per_width": [float(b) for b in bend],
        "n_points": int(M.size), "map_seconds": float(dt),
    }


def transfer_curve(depth, drives=None, H=0.04, B_max=3.0, n=160, **kw):
    """Invert the 1-shot to the flux REGULATOR transfer Psi(drive): for each drive
    (MMF) the flux that passes.  Saturation makes it CLAMP (concave).  Returns the
    curve + the knee (where the incremental permeance dPsi/d(drive) has fallen to half
    its small-signal value -- the clamp onset)."""
    w_throat = H - 2.0 * depth
    Bg = np.linspace(1e-3, B_max, n)                  # throat-field grid
    dr = np.array([one_shot_drive(float(B), float(depth), H=H, **kw) for B in Bg])
    Psi = Bg * w_throat                               # flux at each throat field
    perm = np.gradient(Psi, dr)                       # incremental permeance dPsi/d(drive)
    p0 = float(perm[0])                               # small-signal (unsaturated) permeance
    knee_idx = int(np.argmin(np.abs(perm - 0.5 * p0)))
    out = {
        "depth": float(depth), "w_throat": float(w_throat),
        "drive": dr.tolist(), "Psi": Psi.tolist(), "B_throat": Bg.tolist(),
        "perm0": p0, "knee_drive": float(dr[knee_idx]), "knee_Psi": float(Psi[knee_idx]),
        "knee_B": float(Bg[knee_idx]),
    }
    if drives is not None:                            # sample Psi at requested drives
        out["Psi_at"] = [float(np.interp(d, dr, Psi)) for d in drives]
        out["drives"] = list(drives)
    return out


# --------------------------------------------------------------------------- #
# validation against the full nonlinear FEM Picard loop (the cost the sweep avoids)
# --------------------------------------------------------------------------- #
def validate_vs_fem(samples=((0.012, (0.008, 0.024)),), mur0=200.0, Bk=1.0,
                    L=0.20, H=0.04, Rc=0.11, order=2, maxh=0.006):
    """Run the full 2-D nonlinear FEM Picard at sample (depth, Psi-list) points and
    compare to the 1-shot.  Returns the agreement + the Picard iteration counts (the
    per-point cost the linear sweep replaces with one quadrature)."""
    rows = []
    for depth, Psi_list in samples:
        r = ch.solve_chaplygin(Psi_list=Psi_list, mur0=mur0, Bk=Bk, L=L, H=H,
                               depth=depth, Rc=Rc, order=order, maxh=maxh)
        for Psi, Bt, nit, fem, one, rel in r["rows"]:
            rows.append({"depth": float(depth), "B_throat": float(Bt),
                         "fem_drive": float(fem), "one_shot_drive": float(one),
                         "rel_err": float(rel), "picard_iters": int(nit)})
    return {
        "rows": rows,
        "max_rel_err": float(max(r["rel_err"] for r in rows)),
        "mean_picard_iters": float(np.mean([r["picard_iters"] for r in rows])),
    }


def cost_summary(mp, val):
    """The linear-cost win: M map points at 1 quadrature each vs the equivalent
    nonlinear FEM (M points x mean Picard iterations linear solves)."""
    M = mp["n_points"]
    iters = val["mean_picard_iters"]
    return {
        "map_points": M, "map_ms": mp["map_seconds"] * 1e3,
        "mean_picard_iters": iters,
        "fem_equiv_linear_solves": int(round(M * iters)),
        "one_shot_solves": M,
        "iters_saved_per_point": iters,
    }


def run(with_fem=False):
    mp = design_map()
    curves = [transfer_curve(d) for d in mp["depths"]]
    out = {"design_map": mp, "transfer": curves}
    if with_fem:
        val = validate_vs_fem(samples=((0.008, (0.008, 0.024)),
                                       (0.016, (0.004, 0.012))))
        out["validation"] = val
        out["cost"] = cost_summary(mp, val)
    return out


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.6, 4.2), dpi=150)
    mp = out["design_map"]
    B = np.array(mp["B_throats"])
    colors = ["C0", "C1", "C2", "C3"]

    # Panel A: the design map -- drive vs throat field, one curve per throat width
    for i, (w, drv) in enumerate(zip(mp["widths"], mp["drive"])):
        axA.plot(B, np.array(drv), "o-", color=colors[i % 4], ms=3,
                 label=f"w={w*1e3:.0f} mm (bend x{mp['bend_per_width'][i]:.1f})")
    val = out.get("validation")
    if val is not None:
        bb = [r["B_throat"] for r in val["rows"]]
        ff = [r["fem_drive"] for r in val["rows"]]
        axA.plot(bb, ff, "kx", ms=9, mew=2, label="nonlinear FEM (validate)")
    axA.set_xlabel("throat field  $|B|_{throat}$  [T]"); axA.set_ylabel("drive  $\\Delta\\Phi$  [A]")
    axA.set_title("design map: drive vs field, per throat width\n(each point = ONE linear 1-shot)")
    axA.legend(fontsize=7, loc="upper left")

    # Panel B: the flux-regulator transfer curve Psi(drive) -- the clamp
    for i, c in enumerate(out["transfer"]):
        axB.plot(np.array(c["drive"]), np.array(c["Psi"]) * 1e3, "-", color=colors[i % 4],
                 label=f"w={c['w_throat']*1e3:.0f} mm")
        axB.plot(c["knee_drive"], c["knee_Psi"] * 1e3, "o", color=colors[i % 4], ms=6)
    axB.set_xlabel("drive  $\\Delta\\Phi$  [A]"); axB.set_ylabel("flux  $\\Psi$  [mWb/m]")
    axB.set_title("flux regulator: $\\Psi$(drive) CLAMPS at saturation\n(narrower throat clamps lower; o = knee)")
    axB.set_xlim(0, max(c["knee_drive"] for c in out["transfer"]) * 3.0)
    axB.legend(fontsize=8, loc="lower right")

    # Panel C: the linear-cost win
    cost = out.get("cost")
    if cost is not None:
        bars = axC.bar(["1-shot\n(this sweep)", "nonlinear FEM\n(equivalent)"],
                       [cost["one_shot_solves"], cost["fem_equiv_linear_solves"]],
                       color=["C0", "C3"])
        axC.set_yscale("log"); axC.set_ylabel("linear solves for the whole map")
        for b, v in zip(bars, [cost["one_shot_solves"], cost["fem_equiv_linear_solves"]]):
            axC.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=9)
        axC.set_title(f"{cost['map_points']} nonlinear points in {cost['map_ms']:.0f} ms\n"
                      f"(FEM: ~{cost['mean_picard_iters']:.0f} Picard iters EACH)")
    else:
        axC.text(0.5, 0.5, "cost (vs nonlinear FEM):\nrun with --fem", ha="center", va="center",
                 transform=axC.transAxes, fontsize=11)
        axC.set_xticks([]); axC.set_yticks([]); axC.set_title("the linear-cost win")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    with_fem = "--fem" in sys.argv
    print("Nonlinear-as-linear: a saturable-magnet DESIGN SWEEP at linear cost\n")
    out = run(with_fem=with_fem)
    mp = out["design_map"]
    print(f"  design map: {mp['n_points']} nonlinear operating points "
          f"({len(mp['depths'])} throat widths x {len(mp['B_throats'])} fields) "
          f"in {mp['map_seconds']*1e3:.1f} ms (mesh-free 1-shot each):")
    for w, drv, bend in zip(mp["widths"], mp["drive"], mp["bend_per_width"]):
        print(f"    w_throat = {w*1e3:4.0f} mm:  drive {drv[0]:7.1f} A (B={mp['B_throats'][0]:.1f} T) "
              f"-> {drv[-1]:7.1f} A (B={mp['B_throats'][-1]:.1f} T)  bend x{bend:.2f} (saturation)")
    print(f"  flux-regulator clamp (the field-limiter design knob):")
    for c in out["transfer"]:
        print(f"    w_throat = {c['w_throat']*1e3:4.0f} mm:  clamp knee at Psi = {c['knee_Psi']*1e3:.2f} "
              f"mWb/m  (B_throat = {c['knee_B']:.2f} T, drive = {c['knee_drive']:.0f} A)")
    print(f"    -> a narrower throat clamps the flux at a LOWER level (knee Psi ~ Bk*w_throat).")
    if with_fem:
        val, cost = out["validation"], out["cost"]
        print(f"  nonlinear FEM validation (Picard loop at {len(val['rows'])} sample points):")
        for r in val["rows"]:
            print(f"    w={ (0.04-2*r['depth'])*1e3:.0f} mm, B={r['B_throat']:.2f} T:  "
                  f"FEM {r['fem_drive']:7.1f} A vs 1-shot {r['one_shot_drive']:7.1f} A  "
                  f"(rel.err {r['rel_err']:.1e}, {r['picard_iters']} Picard iters)")
        print(f"    max rel.err = {val['max_rel_err']:.1e} (slenderness-limited)")
        print(f"  THE LINEAR-COST WIN: {cost['map_points']} nonlinear points in "
              f"{cost['map_ms']:.0f} ms = {cost['one_shot_solves']} quadratures; the equivalent")
        print(f"    nonlinear FEM = {cost['map_points']} x ~{cost['mean_picard_iters']:.0f} Picard "
              f"iters = ~{cost['fem_equiv_linear_solves']} curved-mesh linear solves.")
    else:
        print("  nonlinear FEM validation + cost: run with  --fem  (slow)")
    print("\n  => the hodograph made saturation a COEFFICIENT, so an ENTIRE nonlinear design")
    print("     space is one linear quadrature per point -- explored at linear cost, where the")
    print("     nonlinear FEM Picard would be prohibitive.")
    _plot(out)


if __name__ == "__main__":
    main()
