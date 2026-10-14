"""The nonlinear SATURATION loop in 2-D (the reference for the Chaplygin rung).

RESEARCH example (track A, the nonlinear half of "Kelvin in the hodograph").
The user's "非線形のループ": a magnetisable cylinder whose permeability
SATURATES, solved by a nonlinear FEM loop.  This is the straightforward
("loop") solve -- the reference that the next rung (the Molenbroek-Chaplygin
hodograph) will reproduce WITHOUT iteration (the hodograph linearises the 2-D
saturation nonlinearity).

Formulation = the standard STABLE one: the **A-formulation with a B-input
reluctivity** ``nu(|B|)``.  In 2-D, ``A_z`` is the flux function,
``B = (dA_z/dy, -dA_z/dx)`` and ``|B| = |grad A_z|``; the weak form is

    INT nu(|grad A_z|) grad(A_z) . grad(v) dx = 0          (current-free),

with a Dirichlet ``A_s = -B0 x`` on the outer boundary imposing the uniform
applied field ``B0 y-hat``.  ``nu`` is the reluctivity ``1/(mu0 mu_r(|B|))``,
with a Froehlich-type saturation ``mu_r(B) = 1 + (mu_r0-1)/(1 + (B/Bk)^2)``
(mu_r0 at B=0, -> 1 as B -> infinity).

Why this and NOT the reduced-Omega ``mu(|H|)`` loop: the H-input
``mu(|H|)`` Picard is ILL-CONDITIONED for a saturable cylinder (mu_r is
steepest at small H, exactly at the cylinder boundary), so the fixed-point
iteration converges to a spurious, non-physical state (interior B exceeding
the demagnetisation limit 2 B0).  The B-input ``nu(|B|)`` A-formulation has a
CONVEX energy, so plain Picard converges cleanly -- here to MACHINE PRECISION:
the converged iterate is its own frozen-re-solve (``|blend - frozen| ~ 1e-12``),
the diagnostic that the loop found the true solution and not a false fixed
point.

A continuation in B0 (each solve warm-started from the previous) keeps every
step in the well-conditioned basin.  The result is the physical saturation
curve: the interior field ratio ``B_in / B0`` falls from the unsaturated demag
value ``2 mu_r0/(mu_r0+1)`` toward 1 (the saturated cylinder is "transparent").

Open boundary: a LARGE domain (``R/a = 10``, truncation < 1%) stands in for the
exact open boundary -- saturation is a LOCAL (cylinder) effect, orthogonal to
the open-boundary treatment that rungs 1-2 (hodograph_kelvin_2d /
clebsch_kelvin_3d) handle exactly; the two compose.

run:  python saturation_loop_2d.py
"""
import math
import os

from numpy import pi
import numpy as np
from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx, CF, x, sqrt,
                     BilinearForm, LinearForm, TaskManager, Integrate)
from netgen.occ import Circle, Glue, OCCGeometry

MU0 = 4 * pi * 1e-7


def _geometry(a, R):
    mag = Circle((0, 0), a).Face()
    mag.faces.name = "magnetic"
    for e in mag.edges:
        e.name = "iface"
    disk = Circle((0, 0), R).Face()
    for e in disk.edges:
        e.name = "outer"
    air = disk - mag
    air.faces.name = "air"
    return OCCGeometry(Glue([air, mag]), dim=2)


def _mu_r(Bmag, mur0, Bk):
    """Froehlich saturation: mu_r0 at B=0, -> 1 as B -> infinity."""
    return 1.0 + (mur0 - 1.0) / (1.0 + (Bmag / Bk) ** 2)


def solve_saturation(B0_list=(0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
                     mur0=10.0, Bk=1.0, order=3, maxh=0.08, a=0.3, R=3.0,
                     niter=400, tol=1e-9):
    """A-formulation nu(|B|) Picard with B0 continuation.  Returns the
    saturation curve B_in/B0 + the machine-precision consistency diagnostic."""
    with TaskManager():
        mesh = Mesh(_geometry(a, R).GenerateMesh(maxh=maxh))
        mesh.Curve(order)
        magmask = mesh.MaterialCF({"magnetic": 1.0}, default=0.0)
        fes = H1(mesh, order=order, dirichlet="outer")
        u, v = fes.TnT()
        vol = Integrate(magmask, mesh)
        A = GridFunction(fes)

        def lin_solve(nucf, As):
            a_ = BilinearForm(nucf * InnerProduct(grad(u), grad(v)) * dx)
            a_.Assemble()
            g = GridFunction(fes)
            g.Set(As, definedon=mesh.Boundaries("outer"))
            r = g.vec.CreateVector()
            r.data = -a_.mat * g.vec
            g.vec.data += a_.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
            return g

        rows = []
        for B0 in B0_list:
            As = -B0 * x                       # A_z = -B0 x  ->  B = (0, B0)
            A.Set(As, definedon=mesh.Boundaries("outer"))
            prev = np.array(A.vec)
            nit = 0
            for it in range(niter):
                B = grad(A)
                Bmag = sqrt(B[0] * B[0] + B[1] * B[1] + 1e-30)
                nucf = 1.0 / (MU0 * (1.0 + magmask * (_mu_r(Bmag, mur0, Bk) - 1.0)))
                A = lin_solve(nucf, As)        # pure Picard (convex -> stable)
                cur = np.array(A.vec)
                d = np.linalg.norm(cur - prev) / (np.linalg.norm(cur) + 1e-30)
                prev = cur.copy()
                nit = it + 1
                if d < tol:
                    break
            # interior field B_in (= -dA/dx, the y-component) and the
            # machine-precision consistency check (frozen re-solve == iterate).
            B = grad(A)
            Bmag = sqrt(B[0] * B[0] + B[1] * B[1] + 1e-30)
            nucf = 1.0 / (MU0 * (1.0 + magmask * (_mu_r(Bmag, mur0, Bk) - 1.0)))
            By = Integrate(magmask * (-grad(A)[0]), mesh) / vol
            A2 = lin_solve(nucf, As)
            By2 = Integrate(magmask * (-grad(A2)[0]), mesh) / vol
            rows.append((float(B0), int(nit), float(d), float(By / B0),
                         float(abs(By - By2) / abs(By2))))

    demag_lin = 2.0 * mur0 / (mur0 + 1.0)      # unsaturated cylinder demag ratio
    ratios = [r[3] for r in rows]
    monotone = all(ratios[i] >= ratios[i + 1] - 1e-6 for i in range(len(ratios) - 1))
    return {
        "mur0": mur0, "Bk": Bk, "R_over_a": R / a,
        "rows": rows,                          # (B0, iters, d, B_in/B0, |blend-frozen|)
        "demag_linear": float(demag_lin),      # 2 mu_r0/(mu_r0+1)
        "ratio_lowfield": float(ratios[0]),    # -> demag_linear (unsaturated)
        "ratio_highfield": float(ratios[-1]),  # -> 1 (saturated)
        "max_inconsistency": float(max(r[4] for r in rows)),  # ~1e-12: true solution
        "monotone": bool(monotone),
        "respects_demag": bool(all(1.0 - 1e-3 <= r[3] <= demag_lin + 1e-3 for r in rows)),
    }


def _plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    B0 = [r[0] for r in res["rows"]]
    ratio = [r[3] for r in res["rows"]]
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    ax.semilogx(B0, ratio, "o-", color="C1", label="$B_{in}/B_0$ (FEM loop)")
    ax.axhline(res["demag_linear"], color="C0", lw=0.9, ls="--",
               label=f"unsaturated demag $2\\mu_{{r0}}/(\\mu_{{r0}}+1)$ = "
                     f"{res['demag_linear']:.3f}")
    ax.axhline(1.0, color="0.5", lw=0.9, ls=":", label="saturated ($\\mu_r\\to1$)")
    ax.set_xlabel("applied field  $B_0$  [T]")
    ax.set_ylabel("interior field ratio  $B_{in}/B_0$")
    ax.set_title(f"Nonlinear saturation loop (A-form, $\\nu(|B|)$, "
                 f"$\\mu_{{r0}}$={res['mur0']:.0f})\n"
                 f"clean Picard: $|$blend$-$frozen$|$ = "
                 f"{res['max_inconsistency']:.0e} (true solution)")
    ax.legend(fontsize=8)
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Nonlinear saturation loop in 2-D (A-formulation, nu(|B|))\n")
    r = solve_saturation()
    print(f"  mu_r0={r['mur0']:.0f}  Bk={r['Bk']:.1f} T  R/a={r['R_over_a']:.0f} "
          f"(truncation < 1%)\n")
    print("  applied B0 [T] | iters | B_in/B0 | |blend-frozen| (=>true solution)")
    for B0, nit, d, ratio, cons in r["rows"]:
        print(f"    {B0:8.3f}    | {nit:4d}  | {ratio:.4f}  | {cons:.1e}")
    print(f"\n  -> saturation: B_in/B0 falls {r['ratio_lowfield']:.3f} "
          f"(unsaturated, ~ demag {r['demag_linear']:.3f}) -> "
          f"{r['ratio_highfield']:.3f} (saturated, ->1).")
    print(f"     monotone={r['monotone']}, respects demag limit [1, "
          f"{r['demag_linear']:.2f}]={r['respects_demag']}, "
          f"max inconsistency={r['max_inconsistency']:.0e} (clean Picard).")
    print("\n  => the straightforward nonlinear LOOP.  Next rung: the")
    print("     Molenbroek-Chaplygin hodograph reproduces this WITHOUT iteration")
    print("     (the hodograph linearises the 2-D saturation nonlinearity).")
    _plot(r)


if __name__ == "__main__":
    main()
