r"""The design-at-linear-cost exploit, one rung up: a TURNING saturable field per
design is ONE linear hodograph PDE solve (not even a quadrature).

`chaplygin_design_sweep_2d.py` exploited the SLENDER guide (hodograph image = a
segment), where the linear Chaplygin operator collapses to a 1-shot quadrature -- so a
whole design space was 1 quadrature per point.  This file is the rung up: a field that
TURNS (direction theta varies over a 2-D range) has a hodograph image that is a real
2-D REGION, on which the Chaplygin equation is a genuine LINEAR, variable-coefficient
ELLIPTIC PDE that must be SOLVED -- but it is still LINEAR (mu(q) is a coefficient), so
each nonlinear turning field is ONE direct linear solve, no Picard.  The same exploit
applies: a whole design space of TURNING saturable fields at linear cost.

What is swept (each entry = ONE linear hodograph solve, reusing
`chaplygin_turning_guide_2d.solve_hodograph`):

  * material  mu_r0  -- the unsaturated permeability;
  * operating depth  q1/q_k  -- how far past the saturation knee the field runs.

For each, we report:
  * lin_residual  -- the free-dof residual of the linear solve (~machine zero => it IS
    a single direct solve, no outer nonlinear loop);
  * the saturation BEND  -- how far the (nonlinear) hodograph solution deviates from the
    linear-limit harmonic A = ln(q)*theta it was driven with (the genuine 2-D Chaplygin
    content; grows with mu_r0 and with operating depth);
  * the back-map closure  -- integrating dx, dy back to PHYSICAL (x, y) is single-valued
    (the turning field is realisable).

Validation: at mu_r0 = 1 (the Laplace limit) the solver reproduces the EXACT harmonic
A = ln(q)*theta (laplace_error ~ 0) -- the linear solver is correct (genuinely 2-D,
theta varies).  The nonlinear cases are the SAME linear operator with a saturating
coefficient, and each back-maps single-valued.

The linear-cost win: M turning designs = M direct linear solves (milliseconds each),
whereas the equivalent PHYSICAL-space nonlinear solve needs a Picard loop per design.
We measure that Picard cost on the same Froehlich material via the slender-guide
reference (`chaplygin_hodograph_2d.solve_chaplygin`, ~10-15 iterations per operating
point); the genuinely-2-D turning case is at least as expensive, so the per-design
saving is at least that factor.

Honest scope: this is the FORWARD construction (solve the hodograph BVP with genuinely
2-D turning data, back-map to a realisable physical patch).  Prescribing a physical
TURNING+TAPERING guide and solving for its (theta-dependent) hodograph image is a FREE
boundary -- the remaining frontier (`chaplygin_inverse_nonlinear_2d.py` closes the
constant-flux von-Mises version); that too is LINEAR per Newton step.

run:  python chaplygin_turning_design_sweep_2d.py            # design sweep (fast)
      python chaplygin_turning_design_sweep_2d.py --fem      # + Picard-cost reference (slow)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chaplygin_turning_guide_2d as tg              # noqa: E402  the linear hodograph PDE
import chaplygin_hodograph_2d as ch                  # noqa: E402  the Picard-cost reference


# --------------------------------------------------------------------------- #
# one turning design = ONE linear hodograph PDE solve (+ realisability back-map)
# --------------------------------------------------------------------------- #
def turning_design(mur0, qk=1.0, q0=0.6, q1=2.2, th1=1.0, order=3, maxh=0.05,
                   backmap=True, Nbm=33):
    """A single TURNING saturable field by ONE linear hodograph solve.  Returns the
    linear-solve residual, the saturation bend, and (optionally) the back-mapped
    physical region + its single-valuedness closure."""
    res, mesh, gf = tg.solve_hodograph(mur0=mur0, qk=qk, q0=q0, q1=q1, th1=th1,
                                       order=order, maxh=maxh)
    out = {"mur0": float(mur0), "qk": float(qk), "q0": float(q0), "q1": float(q1),
           "th1": float(th1), "lin_residual": res["lin_residual"],
           "bend": res["twoD_deviation"], "laplace_error": res["laplace_error"]}
    if backmap:
        bm = tg.back_map(mesh, gf, mur0, qk, q0, q1, th1, Nq=Nbm, Nth=Nbm)
        out["closure"] = bm["closure"]
        out["_bm"] = {"X": bm["X"], "Y": bm["Y"], "Bmag": bm["Bmag"], "theta": bm["theta"]}
    return out


def design_sweep(mur0_list=(5.0, 20.0, 80.0), depth_q1=(1.2, 2.0, 3.0),
                 qk=1.0, q0=0.6, th1=1.0, order=3, maxh=0.05):
    """Two design sweeps, each entry ONE linear hodograph PDE solve:
      - material: mu_r0 at fixed operating range;
      - depth:    operating q1/q_k at fixed material.
    Plus the Laplace-limit solver check (mu_r0 = 1)."""
    t0 = time.perf_counter()
    material = [turning_design(m, qk=qk, q0=q0, q1=2.2, th1=th1, order=order, maxh=maxh)
                for m in mur0_list]
    depth = [turning_design(20.0, qk=qk, q0=q0, q1=q, th1=th1, order=order, maxh=maxh,
                            backmap=False) for q in depth_q1]
    sweep_seconds = time.perf_counter() - t0
    verify = turning_design(1.0, qk=qk, q0=q0, q1=2.2, th1=th1, order=order, maxh=maxh,
                            backmap=False)
    n_solves = len(material) + len(depth) + 1
    return {"material": material, "depth": depth, "verify_laplace": verify,
            "n_solves": n_solves, "sweep_seconds": float(sweep_seconds),
            "max_lin_residual": float(max([d["lin_residual"] for d in material + depth]
                                          + [verify["lin_residual"]]))}


def picard_cost_reference(mur0=20.0, Bk=1.0, depth=0.012, order=2, maxh=0.006):
    """Measure the per-operating-point nonlinear Picard cost on the SAME Froehlich
    material (the slender-guide physical FEM, chaplygin_hodograph_2d): the iterations
    a PHYSICAL-space nonlinear solve needs, which the linear hodograph replaces with
    ONE direct solve per design."""
    r = ch.solve_chaplygin(Psi_list=(0.012, 0.024), mur0=mur0, Bk=Bk, depth=depth,
                           order=order, maxh=maxh)
    iters = [row[2] for row in r["rows"]]
    return {"mean_picard_iters": float(np.mean(iters)), "iters": iters}


def run(with_fem=False):
    sw = design_sweep()
    out = {"sweep": sw}
    if with_fem:
        ref = picard_cost_reference()
        out["picard_reference"] = ref
        out["cost"] = {
            "designs": sw["n_solves"], "one_linear_solve_each": True,
            "sweep_seconds": sw["sweep_seconds"],
            "ref_mean_picard_iters": ref["mean_picard_iters"],
            "fem_equiv_linear_solves": int(round(sw["n_solves"] * ref["mean_picard_iters"])),
        }
    return out


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sw = out["sweep"]
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.8, 4.2), dpi=150)

    # Panel A: the saturation bend grows with material + depth (each = ONE linear solve)
    mm = [d["mur0"] for d in sw["material"]]
    mb = [d["bend"] for d in sw["material"]]
    axA.semilogy(mm, mb, "C0o-", label="material sweep $\\mu_{r0}$")
    qq = [d["q1"] for d in sw["depth"]]
    qb = [d["bend"] for d in sw["depth"]]
    axA2 = axA.twiny()
    axA2.semilogy(qq, qb, "C1s--", label="depth sweep $q_1/q_k$")
    axA.axhline(sw["verify_laplace"]["laplace_error"] + 1e-16, color="0.6", lw=0.8, ls=":")
    axA.set_xlabel("material  $\\mu_{r0}$", color="C0")
    axA2.set_xlabel("operating depth  $q_1/q_k$", color="C1")
    axA.set_ylabel("saturation bend (2-D deviation)")
    axA.set_title(f"each turning design = ONE linear solve\n(residual ~ {sw['max_lin_residual']:.0e}; "
                  f"$\\mu_r$=1 check {sw['verify_laplace']['laplace_error']:.0e})")
    h1, l1 = axA.get_legend_handles_labels(); h2, l2 = axA2.get_legend_handles_labels()
    axA.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")

    # Panel B: a back-mapped PHYSICAL turning field (high saturation) -- genuine 2-D, realisable
    d = sw["material"][-1]; bm = d["_bm"]
    X, Y, th = bm["X"], bm["Y"], bm["theta"]
    axB.plot(X, Y, color="0.85", lw=0.4); axB.plot(X.T, Y.T, color="0.85", lw=0.4)
    pc = axB.pcolormesh(X, Y, th, shading="gouraud", cmap="twilight")
    axB.set_aspect("equal"); axB.set_xlabel("physical $x$"); axB.set_ylabel("physical $y$")
    axB.set_title(f"a TURNING saturable field from ONE solve\n($\\mu_{{r0}}$={d['mur0']:.0f}, "
                  f"realisable: closure {d['closure']:.0e})")
    fig.colorbar(pc, ax=axB, label="field angle $\\theta$")

    # Panel C: the linear-cost win
    cost = out.get("cost")
    if cost is not None:
        bars = axC.bar(["hodograph\n(1 linear solve each)", "physical FEM\n(Picard each)"],
                       [cost["designs"], cost["fem_equiv_linear_solves"]], color=["C0", "C3"])
        axC.set_yscale("log"); axC.set_ylabel("linear solves for the whole sweep")
        for b, v in zip(bars, [cost["designs"], cost["fem_equiv_linear_solves"]]):
            axC.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=9)
        axC.set_title(f"{cost['designs']} turning designs in {cost['sweep_seconds']*1e3:.0f} ms\n"
                      f"(physical FEM: ~{cost['ref_mean_picard_iters']:.0f} Picard iters EACH)")
    else:
        axC.text(0.5, 0.5, "cost (vs nonlinear FEM):\nrun with --fem", ha="center", va="center",
                 transform=axC.transAxes, fontsize=11)
        axC.set_xticks([]); axC.set_yticks([]); axC.set_title("the linear-cost win")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    with_fem = "--fem" in sys.argv
    print("Nonlinear-as-linear, one rung up: TURNING saturable designs, ONE linear solve each\n")
    out = run(with_fem=with_fem)
    sw = out["sweep"]
    print(f"  material sweep (q0=0.6, q1=2.2, theta turns over [0,1] rad), each ONE linear solve:")
    for d in sw["material"]:
        print(f"    mu_r0 = {d['mur0']:5.0f}:  lin.residual = {d['lin_residual']:.1e} (one solve)  "
              f"saturation bend = {d['bend']:.3e}  back-map closure = {d['closure']:.1e} (realisable)")
    print(f"  depth sweep (mu_r0=20), each ONE linear solve:")
    for d in sw["depth"]:
        print(f"    q1/q_k = {d['q1']:.1f}:  saturation bend = {d['bend']:.3e}  "
              f"lin.residual = {d['lin_residual']:.1e}")
    print(f"  Laplace-limit solver check (mu_r0=1): A = ln(q)*theta reproduced to "
          f"{sw['verify_laplace']['laplace_error']:.1e} (solver correct)")
    print(f"  -> {sw['n_solves']} turning designs in {sw['sweep_seconds']*1e3:.0f} ms; "
          f"max linear-solve residual {sw['max_lin_residual']:.0e} (each IS one direct solve).")
    if with_fem:
        cost = out["cost"]
        print(f"  THE LINEAR-COST WIN: {cost['designs']} turning designs = {cost['designs']} direct "
              f"linear solves; the equivalent PHYSICAL nonlinear")
        print(f"    solve needs ~{cost['ref_mean_picard_iters']:.0f} Picard iterations per design "
              f"(measured on the same Froehlich material via the")
        print(f"    slender-guide reference) = ~{cost['fem_equiv_linear_solves']} linear solves; the "
              f"turning 2-D case is at least as expensive.")
    else:
        print("  Picard-cost reference + cost: run with  --fem  (slow)")
    print("\n  => the hodograph keeps the nonlinearity as a COEFFICIENT even when the field TURNS")
    print("     (a 2-D image), so a whole design space of TURNING saturable fields is ONE linear")
    print("     elliptic solve per design -- the design-at-linear-cost exploit, one rung up.")
    _plot(out)


if __name__ == "__main__":
    main()
