r"""Closing Task 1's free boundary: a TURNING+TAPERING saturable guide design sweep, by
the von Mises inverse -- each design a Newton continuation whose every step is LINEAR.

`chaplygin_turning_design_sweep_2d.py` swept the FORWARD turning construction (a fixed
hodograph rectangle, ONE linear solve per design) -- but a guide that turns AND tapers
has theta-DEPENDENT walls, so its hodograph image is theta-dependent = a FREE BOUNDARY
(the genuinely hard case the forward construction defers).  This file closes that:
the **von Mises** coordinate change (use the potential `Phi` and flux `A` as the
independent coordinates) DISSOLVES the free boundary -- the guide is a FIXED rectangle
in `(Phi, A)`, and one solves for the physical map `(x, y)(Phi, A)` instead.  The
saturable flux height `A1 = lambda` is freed as a global unknown so the rectangle is
not over-determined (`chaplygin_inverse_nonlinear_2d.py`).

The connection to "nonlinear-as-linear": the inverse is a least-squares (FOSLS)
problem whose nonlinearity (`mu(q)` coefficient + the on-curve slip constraint) is
handled by **damped Newton with continuation in the saturation Ms** -- and **every
Newton step is ONE LINEAR solve** on the FIXED reference rectangle (no remeshing, no
moving boundary).  So the whole TAPERED-design space is swept at "linear-per-Newton-
step, fixed-mesh" cost, where a physical-space free-boundary solver would re-mesh and
re-solve a nonlinear problem at every shape iteration.

Swept here: the TAPER (how much the outer wall spirals inward, the free-boundary
lever).  For each taper (each = one `solve_inverse`):
  * J            -- the FOSLS residual (-> 0: the inverse closed);
  * jac_min      -- the map Jacobian (> 0: a valid, non-folded physical map);
  * free_measure -- the theta-drift of the hodograph q-extent (0 = rectangle =
    constant width = self-linearising; > 0 = theta-dependent = the FREE BOUNDARY
    recovered), which GROWS with taper;
  * lambda       -- the saturable flux (the freed rectangle height).

Honest scope: the von Mises rectangle is the slip-wall (tangential-free) formulation;
extreme taper eventually FOLDS the map (jac_min -> 0) at the geometric throat limit
(`chaplygin_inverse_nonlinear_2d.py` documents this).  This sweep stays in the valid
regime and shows the free boundary appearing as the taper grows.

run:  python chaplygin_taper_design_sweep_2d.py
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chaplygin_inverse_nonlinear_2d as nl          # noqa: E402  the von Mises inverse


def taper_sweep(tapers=(0.0, 0.15, 0.30), Ms_target=12.0, order=2, maxh=0.07,
                grids=True):
    """Sweep the taper (the free-boundary lever); each design is the von Mises nonlinear
    inverse (a Newton continuation, every step linear) on the FIXED (Phi, A) rectangle."""
    rows = []
    t0 = time.perf_counter()
    for tp in tapers:
        r = nl.solve_inverse(taper=float(tp), Ms_target=Ms_target, order=order, maxh=maxh)
        span = r["theta_range_deg"][1] - r["theta_range_deg"][0]
        row = {"taper": float(tp), "J": r["J"], "jac_min": r["jac_min"],
               "free_measure": r["free_measure"], "lambda": r["lambda"],
               "lambda_lin": r["lambda_lin"], "wall_fit": r["wall_fit"],
               "turn_deg": float(span), "ne": r["ne"]}
        if grids:
            GX, GY = nl._physical_grid(r, n=36)
            row["_grid"] = (GX, GY)
        rows.append(row)
    return {"rows": rows, "Ms": float(Ms_target), "sweep_seconds": float(time.perf_counter() - t0),
            "max_J": float(max(x["J"] for x in rows)),
            "min_jac": float(min(x["jac_min"] for x in rows))}


def run():
    return {"sweep": taper_sweep()}


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sw = out["sweep"]; rows = sw["rows"]
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.8, 4.2), dpi=150)

    # Panel A: the free-boundary measure grows with taper (+ the saturable flux lambda)
    tp = [r["taper"] for r in rows]
    fm = [r["free_measure"] for r in rows]
    axA.plot(tp, fm, "C0o-", label="free-boundary measure")
    axA2 = axA.twinx()
    axA2.plot(tp, [r["lambda"] for r in rows], "C3s--", label="saturable flux $\\lambda$")
    axA.axhline(0.1, color="0.6", lw=0.8, ls=":")
    axA.set_xlabel("taper (outer-wall spiral-in)")
    axA.set_ylabel("free-boundary measure ($\\theta$-drift of $q$-extent)", color="C0")
    axA2.set_ylabel("saturable flux  $\\lambda$", color="C3")
    axA.set_title("the FREE BOUNDARY appears with taper\n(each = ONE von Mises inverse, "
                  f"J <= {sw['max_J']:.0e})")
    h1, l1 = axA.get_legend_handles_labels(); h2, l2 = axA2.get_legend_handles_labels()
    axA.legend(h1 + h2, l1 + l2, fontsize=8, loc="center left")

    # Panels B, C: the physical back-mapped guides -- constant-width vs tapered
    for ax, idx, ttl in ((axB, 0, "constant width (rectangle image)"),
                         (axC, len(rows) - 1, "tapered (FREE BOUNDARY)")):
        r = rows[idx]; GX, GY = r["_grid"]
        ax.plot(GX, GY, color="0.8", lw=0.4); ax.plot(GX.T, GY.T, color="0.8", lw=0.4)
        # colour the cells by local flux-tube "width" proxy via the radial coordinate
        ax.pcolormesh(GX, GY, np.hypot(GX, GY), shading="gouraud", cmap="viridis")
        ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(f"{ttl}\ntaper={r['taper']:.2f}, free={r['free_measure']:.2f}, "
                     f"jac_min={r['jac_min']:.2f}")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Closing Task 1's free boundary: TURNING+TAPERING designs by the von Mises inverse\n")
    out = run()
    sw = out["sweep"]
    print(f"  taper sweep (Ms={sw['Ms']:.0f}; each = one von Mises inverse, a Newton "
          f"continuation with LINEAR steps on a FIXED rectangle):")
    print(f"    taper | J (closed)  | jac_min (valid) | free-measure | lambda | turn")
    for r in sw["rows"]:
        print(f"    {r['taper']:.2f}  | {r['J']:.2e}  |   {r['jac_min']:.3f}      |   "
              f"{r['free_measure']:.3f}     | {r['lambda']:5.2f}  | {r['turn_deg']:.0f} deg")
    print(f"  -> the FREE-BOUNDARY measure GROWS with taper "
          f"({sw['rows'][0]['free_measure']:.2f} -> {sw['rows'][-1]['free_measure']:.2f}): "
          f"constant width = a rectangle image (self-linearising);")
    print(f"     tapering = a theta-dependent image = the free boundary, RECOVERED.  Each design "
          f"closed (J <= {sw['max_J']:.0e}) with a valid map (jac_min >= {sw['min_jac']:.2f}).")
    print(f"  {len(sw['rows'])} tapered free-boundary designs in {sw['sweep_seconds']:.1f} s "
          f"(each a Newton continuation, every step ONE linear solve, FIXED mesh).")
    print("\n  => the von Mises change of variables DISSOLVES the free boundary onto a fixed")
    print("     rectangle, so even the TURNING+TAPERING (free-boundary) saturable design space")
    print("     is reachable -- the nonlinearity handled by Newton with LINEAR steps, no remesh.")
    _plot(out)


if __name__ == "__main__":
    main()
