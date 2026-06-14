"""Frontier 2, the inverse direction: DISSOLVING the turning-guide free boundary
with the von Mises (Phi, A) coordinate change (linear case solved + verified;
the nonlinear case is the genuine open wall, characterised honestly below).

RESEARCH example (track A).  chaplygin_free_boundary_2d.py showed that a
tapering turning guide has a theta-DEPENDENT hodograph image -- a FREE BOUNDARY
in the (q, theta) plane.  The key idea for the INVERSE problem (prescribe the
physical guide, find its image) is that the free boundary is an artefact of the
COORDINATE CHOICE: take the potential Phi and the flux function A as the
INDEPENDENT variables (the von Mises coordinates).  A flux guide is bounded by
two flux lines (A = 0, A = Psi = the walls) and two equipotentials (Phi = const
= the ports), so in (Phi, A) the domain is ALWAYS the fixed rectangle
[0, Phi1] x [0, A1] -- the free boundary is GONE.  One solves instead for the
physical map (x, y)(Phi, A) on that fixed rectangle.

The governing system (mu0 = 1; derived from the same first-order Chaplygin pair
as chaplygin_turning_guide_2d.py, here with Phi, A independent):

    x_A = -y_Phi / mu(q) ,   y_A =  x_Phi / mu(q) ,   q = 1 / |grad_Phi (x,y)| ,

solved in least squares: minimise INT (x_A + y_Phi/mu)^2 + (y_A - x_Phi/mu)^2.
For mu = 1 this is the Cauchy-Riemann system -- (x + i y) is analytic in
(Phi + i A) -- so the EXACT map for an annular bend is the conformal
``f = e^{i(Phi + iA)}``:  x = e^{-A} cos Phi, y = e^{-A} sin Phi (r = e^{-A},
angle = Phi).  The least-squares solver RECOVERS it to ~1e-8 (rel. error), with
the residual J -> 0 -- the free boundary has been dissolved into a fixed-domain
solve.

THE NONLINEAR WALL (honest -- the genuinely open part).  With mu = mu(q) the
A-spacing 1/(mu q) changes, so the consistent DISTRIBUTION of Phi, A along the
walls is mu-DEPENDENT.  Prescribing the full boundary map (x, y) from a fixed
GEOMETRIC parametrisation (uniform angle) then over-constrains the problem and
the nonlinear least-squares map FOLDS (Jacobian < 0).  The correct boundary
condition is a SLIP condition -- the boundary point lies ON the prescribed wall
curve with its tangential position free -- which for CURVED walls is a NONLINEAR
constraint (e.g. x^2 + y^2 = r_in^2 on an arc).  Combined with the mu(q) Picard
this is a nonlinear PDE with nonlinear boundary constraints: the genuine
free-boundary inverse, still open here.  (Straight-wall guides have LINEAR slip
constraints but are self-linearising = trivial image.)  So this file SOLVES the
inverse in the linear case and pins down EXACTLY what makes the nonlinear case
hard.

run:  python chaplygin_inverse_vonmises_2d.py
"""
import math
import os

import numpy as np
from ngsolve import (Mesh, VectorH1, GridFunction, grad, dx, CF, x, y, exp, cos,
                     sin, BilinearForm, TaskManager, Integrate, InnerProduct)
from netgen.occ import WorkPlane, OCCGeometry


def _rect(Phi1, A1):
    f = WorkPlane().MoveTo(0, 0).Rectangle(Phi1, A1).Face()
    f.faces.name = "dom"
    for e in f.edges:
        e.name = "bnd"
    return OCCGeometry(f, dim=2)


def solve_inverse(Phi1=1.2, A1=0.8, order=3, maxh=0.05, mu=1.0):
    """Least-squares von Mises inverse on the fixed (Phi, A) rectangle.  For
    mu=1 the Dirichlet data is the exact conformal annular-bend map; the solver
    must recover it (the free boundary dissolved)."""
    exact = CF((exp(-y) * cos(x), exp(-y) * sin(x)))    # Phi=x, A=y
    with TaskManager():
        mesh = Mesh(_rect(Phi1, A1).GenerateMesh(maxh=maxh))
        fes = VectorH1(mesh, order=order, dirichlet="bnd")
        U, V = fes.TnT()
        # grad(W)[i,j] = d W_i / d x_j ; x_0 = Phi, x_1 = A
        def R1(W): return grad(W)[0, 1] + grad(W)[1, 0] / mu   # x_A + y_Phi/mu
        def R2(W): return grad(W)[1, 1] - grad(W)[0, 0] / mu   # y_A - x_Phi/mu
        a = BilinearForm(fes)
        a += (R1(U) * R1(V) + R2(U) * R2(V)) * dx
        a.Assemble()
        gf = GridFunction(fes)
        gf.Set(exact, definedon=mesh.Boundaries("bnd"))
        r = gf.vec.CreateVector()
        r.data = -a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
        rel_err = math.sqrt(Integrate(InnerProduct(gf - exact, gf - exact), mesh)
                            / Integrate(InnerProduct(exact, exact), mesh))
        J = float(Integrate(R1(gf) ** 2 + R2(gf) ** 2, mesh))
        # the recovered hodograph image q = 1/|grad_Phi|, theta = arg(grad_Phi)
        qs, ths = [], []
        for px in np.linspace(0.1 * Phi1, 0.9 * Phi1, 24):
            for py in np.linspace(0.1 * A1, 0.9 * A1, 12):
                mp = mesh(px, py)
                xP = float(grad(gf)[0, 0](mp)); yP = float(grad(gf)[1, 0](mp))
                qs.append(1.0 / math.hypot(xP, yP))
                ths.append(math.atan2(yP, xP))
    qs = np.array(qs); ths = np.array(ths)
    return {
        "ne": int(mesh.ne), "rel_err": float(rel_err), "J": J,
        "q_range": (float(qs.min()), float(qs.max())),
        "theta_range_deg": (float(np.degrees(ths.min())),
                            float(np.degrees(ths.max()))),
        "gf": gf, "mesh": mesh, "Phi1": Phi1, "A1": A1,
    }


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    gf, mesh, Phi1, A1 = r["gf"], r["mesh"], r["Phi1"], r["A1"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.2), dpi=150)
    # the fixed (Phi, A) rectangle (the free boundary is GONE)
    gp, ga = np.linspace(0, Phi1, 13), np.linspace(0, A1, 9)
    for a_ in ga:
        ax.plot([0, Phi1], [a_, a_], color="0.8", lw=0.6)
    for p_ in gp:
        ax.plot([p_, p_], [0, A1], color="0.8", lw=0.6)
    ax.add_patch(plt.Rectangle((0, 0), Phi1, A1, fill=False, color="C0", lw=2))
    ax.set_xlabel(r"$\Phi$ (potential)"); ax.set_ylabel("$A$ (flux)")
    ax.set_title("von Mises domain: a FIXED rectangle\n(the free boundary is "
                 "dissolved)")
    ax.set_aspect("equal")
    # the recovered physical map (image of the rectangle grid)
    n = 40
    GP, GA = np.meshgrid(np.linspace(0, Phi1, n), np.linspace(0, A1, n))
    X = np.zeros_like(GP); Y = np.zeros_like(GP)
    for i in range(n):
        for j in range(n):
            v = gf(mesh(float(GP[i, j]), float(GA[i, j])))
            X[i, j] = float(v[0]); Y[i, j] = float(v[1])
    ax2.plot(X, Y, color="0.8", lw=0.5)
    ax2.plot(X.T, Y.T, color="0.8", lw=0.5)
    ax2.set_xlabel("physical $x$"); ax2.set_ylabel("physical $y$")
    ax2.set_aspect("equal")
    ax2.set_title(f"recovered physical map (annular bend)\n"
                  f"rel.err vs conformal = {r['rel_err']:.1e}, residual "
                  f"J = {r['J']:.0e}")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Frontier 2 inverse: dissolving the free boundary (von Mises)\n")
    r = solve_inverse()
    print(f"  linear (mu=1): ne={r['ne']}, recover the conformal annular-bend "
          f"map to rel.err {r['rel_err']:.2e}, residual J={r['J']:.2e}")
    print(f"    -> the field turns (theta in "
          f"[{r['theta_range_deg'][0]:.0f}, {r['theta_range_deg'][1]:.0f}] deg), "
          f"q in [{r['q_range'][0]:.2f}, {r['q_range'][1]:.2f}]")
    print("  => the free boundary is DISSOLVED: the inverse is a fixed-domain")
    print("     solve in (Phi, A).  Verified in the linear case.  The NONLINEAR")
    print("     free-boundary inverse additionally needs a SLIP boundary")
    print("     condition (point on the wall curve, tangential free) -- a")
    print("     nonlinear constraint for curved walls -- and is the open wall")
    print("     (full Dirichlet folds the nonlinear map; see the docstring).")
    _plot(r)


if __name__ == "__main__":
    main()
