"""Chaplygin, one rung past the 1-shot: a TURNING field solved by ONE LINEAR
PDE on the genuine 2-D hodograph plane (not a quadrature).

RESEARCH example (track A, the Chaplygin frontier).  Rung 1.5b
(chaplygin_hodograph_2d.py) did the SLENDER flux guide whose hodograph image is
a thin SEGMENT (theta ~ const) -- there the linear Chaplygin system collapses to
a 1-shot quadrature.  This file goes one step further: a field that TURNS, so the
direction theta varies over a genuine 2-D range and the hodograph image is a
2-D region.  On it the Chaplygin equation is a genuine LINEAR, variable-
coefficient ELLIPTIC PDE that must be SOLVED (one linear solve), not integrated.

THE LINEAR CHAPLYGIN PDE (derived; verified against the linear Laplace limit)
---------------------------------------------------------------------------
Magnetostatics H = grad(Phi), |H| = q, direction theta; B = mu(q) H; flux
function A with B = (A_y, -A_x).  Taking (q, theta) as the INDEPENDENT
coordinates and inverting

    dPhi = q cos th dx + q sin th dy ,   dA = -mu q sin th dx + mu q cos th dy ,

the integrability of dx, dy (x, y single-valued) gives the first-order pair

    Phi_q   = -q (mu q)' / (mu q)^2 * A_theta ,     A_q = (mu/q) Phi_theta ,

whose elimination is the **linear, q-variable-coefficient, self-adjoint**
elliptic PDE for the flux function A(q, theta):

    d/dq( P(q) A_q ) + Q(q) A_theta_theta = 0 ,
        P(q) = q / mu(q) ,     Q(q) = (mu(q) q)' / (mu(q)^2 q) .

mu(q) is a KNOWN COEFFICIENT (a function of the independent coordinate q) -- the
nonlinearity is gone.  Linear limit mu = const: P = q/mu, Q = 1/(mu q) ->
A_qq + (1/q) A_q + (1/q^2) A_theta_theta = 0 = LAPLACE in (ln q, theta).  So
``A = ln(q) * theta`` is an EXACT solution -- the solver's verification.

WHY THIS IS THE FRONTIER (honest)
---------------------------------
A turning flux GUIDE (iron walls = flux lines) on a FIXED hodograph rectangle is
exactly the CONSTANT-WIDTH bend, which is theta-independent (1-D, self-
linearising -- the field |H| ~ 1/r is forced by geometry, mu only reshapes B).
A GENUINELY turning+tapering guide has theta-DEPENDENT wall positions -> a
theta-dependent hodograph image = a FREE BOUNDARY (unknown image), the genuinely
hard case.  So the clean "one linear solve" demonstration is the FORWARD
construction: solve the linear Chaplygin BVP on the rectangle with genuinely 2-D
(theta-dependent) data, then MAP BACK (integrate dx, dy) to a physical patch --
a genuine 2-D nonlinear (saturating) TURNING field obtained from ONE linear
solve.  The free-boundary closure (prescribe the physical guide, solve for its
hodograph image) is the remaining frontier, noted but not claimed here.

run:  python chaplygin_turning_guide_2d.py
"""
import math
import os

from numpy import pi
import numpy as np
from ngsolve import (Mesh, H1, GridFunction, grad, dx, CF, x, y, sqrt, log,
                     BilinearForm, LinearForm, TaskManager, Integrate)
from netgen.occ import WorkPlane, OCCGeometry


def _mu_r(q, mur0, qk):
    """Froehlich saturation in q=|H|: mu_r0 at q=0, -> 1 as q -> infinity."""
    return 1.0 + (mur0 - 1.0) / (1.0 + (q / qk) ** 2)


def _PQ(q, mur0, qk):
    """The Chaplygin coefficients P(q)=q/mu, Q(q)=(mu q)'/(mu^2 q) as CFs of q
    (mu0 cancels -> use mu_r).  (mu_r q)' = mu_r + q mu_r'."""
    mur = _mu_r(q, mur0, qk)
    # d mu_r/dq for Froehlich
    murp = (mur0 - 1.0) * (-2.0 * q / qk ** 2) / (1.0 + (q / qk) ** 2) ** 2
    mq_prime = mur + q * murp
    P = q / mur
    Q = mq_prime / (mur * mur * q)
    return P, Q


def _rect(q0, q1, th1, maxh):
    wp = WorkPlane().MoveTo(q0, 0.0).LineTo(q1, 0.0).LineTo(q1, th1) \
        .LineTo(q0, th1).Close()
    face = wp.Face()
    face.faces.name = "hodo"
    for e in face.edges:
        c = e.center
        if abs(c[1]) < th1 * 1e-3:
            e.name = "th0"
        elif abs(c[1] - th1) < th1 * 1e-3:
            e.name = "th1"
        elif c[0] < 0.5 * (q0 + q1):
            e.name = "q0"
        else:
            e.name = "q1"
    return OCCGeometry(face, dim=2)


def solve_hodograph(mur0=1.0, qk=1.0, q0=0.5, q1=2.0, th1=1.0, order=3,
                    maxh=0.05, bc="lnq_theta"):
    """ONE linear solve of the self-adjoint Chaplygin PDE
    div(P(q) grad_weighted A) = 0 on the hodograph rectangle [q0,q1]x[0,th1],
    P=q/mu, Q=(mu q)'/(mu^2 q), with Dirichlet data on all four sides.

    bc='lnq_theta': Dirichlet A = ln(q)*theta everywhere on the boundary.  For
    mu_r=const this is the EXACT harmonic (Laplace) solution, so the interior
    FEM A must reproduce ln(q)*theta -- the solver verification.  For saturating
    mu_r it drives the SAME genuinely-2-D turning data through the nonlinear
    (variable-coefficient) operator: the interior then DEVIATES from ln(q)*theta
    by the saturation, and that deviation is the genuine 2-D Chaplygin content.
    """
    with TaskManager():
        mesh = Mesh(_rect(q0, q1, th1, maxh).GenerateMesh(maxh=maxh))
        P, Q = _PQ(x, mur0, qk)                          # x is the q-coordinate
        gbnd = log(x) * y                                 # A = ln(q)*theta (y=theta)

        fes = H1(mesh, order=order, dirichlet="q0|q1|th0|th1")
        u, v = fes.TnT()
        a = BilinearForm((P * grad(u)[0] * grad(v)[0]
                          + Q * grad(u)[1] * grad(v)[1]) * dx)
        a.Assemble()
        gf = GridFunction(fes, name="A")
        gf.Set(gbnd, definedon=mesh.Boundaries("q0|q1|th0|th1"))
        r = gf.vec.CreateVector()
        r.data = -a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                     inverse="sparsecholesky") * r

        Aexact = log(x) * y
        # error vs the exact Laplace harmonic (meaningful only for mu_r=const)
        num = Integrate((gf - Aexact) ** 2, mesh)
        den = Integrate(Aexact ** 2, mesh)
        laplace_error = float(math.sqrt(num / den))
        # genuine 2-D content: how far the (nonlinear) solution bends away from
        # the separable ln(q)*theta data it was driven with.
        twoD_deviation = laplace_error
        # the free-dof residual of the LINEAR solve (machine-zero => it IS a
        # single direct linear solve, no Picard / no outer nonlinear loop).
        resvec = gf.vec.CreateVector()
        resvec.data = a.mat * gf.vec
        fd = np.array(fes.FreeDofs(), dtype=bool)
        rv = np.array(resvec)
        lin_residual = float(np.linalg.norm(rv[fd])
                             / (np.linalg.norm(np.array(gf.vec)) + 1e-30))

        return {
            "mur0": mur0, "qk": qk, "order": order, "ne": int(mesh.ne),
            "q_range": (q0, q1), "theta_range": (0.0, th1),
            "laplace_error": laplace_error,          # ~0 when mu_r=const (solver OK)
            "twoD_deviation": twoD_deviation,        # saturation bend (mu_r>1)
            "lin_residual": lin_residual,            # ~machine 0 (one linear solve)
        }, mesh, gf


def verify_solver():
    """Milestone 1: the linear Chaplygin solver reproduces the EXACT Laplace
    harmonic A=ln(q)*theta for mu_r=const (genuinely 2-D, theta varies)."""
    res, _, _ = solve_hodograph(mur0=1.0, order=3, maxh=0.05)
    return res


def nonlinear_turn(mur0=20.0, qk=1.0, order=3, maxh=0.05):
    """The genuine case: the SAME turning data driven through the SATURATING
    variable-coefficient Chaplygin operator -- ONE linear solve (mu(q) is a
    coefficient), the solution bends away from the linear harmonic."""
    res, mesh, gf = solve_hodograph(mur0=mur0, qk=qk, order=order, maxh=maxh)
    return res, mesh, gf


def _mur_np(q, mur0, qk):
    return 1.0 + (mur0 - 1.0) / (1.0 + (q / qk) ** 2)


def _murp_np(q, mur0, qk):
    return (mur0 - 1.0) * (-2.0 * q / qk ** 2) / (1.0 + (q / qk) ** 2) ** 2


def back_map(mesh, gf, mur0, qk, q0, q1, th1, Nq=49, Nth=49):
    """Map the hodograph solution back to PHYSICAL (x, y) by integrating
    dx, dy, and VERIFY single-valuedness (the physical field is realisable).

    dx = (cos th / q) dPhi - (sin th /(mu q)) dA ,
    dy = (sin th / q) dPhi + (cos th /(mu q)) dA ,
    with Phi recovered from  Phi_q = -q (mu q)'/(mu q)^2 A_th ,  Phi_th = (mu/q) A_q
    (mu0 = 1: the physical region up to an overall scale).  Single-valuedness:
    integrate each of x, y by TWO paths (q-then-theta vs theta-then-q); the
    corner mismatch is the path-independence residual (-> 0 iff realisable)."""
    qs = np.linspace(q0, q1, Nq)
    ths = np.linspace(0.0, th1, Nth)
    gA = grad(gf)
    Aq = np.zeros((Nth, Nq))
    Ath = np.zeros((Nth, Nq))
    for i, th in enumerate(ths):
        for j, q in enumerate(qs):
            ga = gA(mesh(float(q), float(th)))
            Aq[i, j] = float(ga[0])
            Ath[i, j] = float(ga[1])

    mur = _mur_np(qs, mur0, qk)[None, :]                  # (1, Nq)
    murp = _murp_np(qs, mur0, qk)[None, :]
    qq = qs[None, :]
    m = mur * qq                                          # mu_r q = |B|
    mprime = mur + qq * murp                              # (mu_r q)'
    Phi_q = -qq * (mprime / m ** 2) * Ath                # Phi_q = -q (mu q)'/(mu q)^2 A_th
    Phi_th = (qq / mur) * Aq                              # Phi_th = (q/mu) A_q  (mu0=1)
    cth = np.cos(ths)[:, None]
    sth = np.sin(ths)[:, None]
    x_q = (cth / qq) * Phi_q - (sth / m) * Aq
    x_th = (cth / qq) * Phi_th - (sth / m) * Ath
    y_q = (sth / qq) * Phi_q + (cth / m) * Aq
    y_th = (sth / qq) * Phi_th + (cth / m) * Ath

    def integ(fq, fth):
        # path A: along theta=0 (use fq), then up each column (use fth)
        row = np.concatenate(([0.0], np.cumsum(
            0.5 * (fq[0, 1:] + fq[0, :-1]) * np.diff(qs))))
        XA = np.zeros((Nth, Nq))
        XA[0, :] = row
        for i in range(1, Nth):
            XA[i, :] = XA[i - 1, :] + 0.5 * (fth[i, :] + fth[i - 1, :]) \
                * (ths[i] - ths[i - 1])
        # path B: up theta=q0 column (use fth), then along each row (use fq)
        col = np.concatenate(([0.0], np.cumsum(
            0.5 * (fth[1:, 0] + fth[:-1, 0]) * np.diff(ths))))
        XB = np.zeros((Nth, Nq))
        XB[:, 0] = col
        for j in range(1, Nq):
            XB[:, j] = XB[:, j - 1] + 0.5 * (fq[:, j] + fq[:, j - 1]) \
                * (qs[j] - qs[j - 1])
        return XA, XB

    XA, XB = integ(x_q, x_th)
    YA, YB = integ(y_q, y_th)
    scale = max(np.ptp(XA), np.ptp(YA), 1e-30)
    closure = float(max(np.max(np.abs(XA - XB)), np.max(np.abs(YA - YB))) / scale)
    Bmag = m + 0.0 * cth                                  # |B| = mu_r q (broadcast)
    return {"qs": qs, "ths": ths, "X": XA, "Y": YA, "Bmag": np.broadcast_to(
        m, (Nth, Nq)).copy(), "theta": np.broadcast_to(ths[:, None],
        (Nth, Nq)).copy(), "closure": closure}


def _plot(lin, nl, bm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.2), dpi=150)
    # left: the hodograph rectangle (the 2-D image -- NOT a segment)
    qs, ths = bm["qs"], bm["ths"]
    QQ, TT = np.meshgrid(qs, ths)
    ax.plot(QQ, TT, color="0.85", lw=0.5)
    ax.plot(QQ.T, TT.T, color="0.85", lw=0.5)
    pc = ax.pcolormesh(QQ, TT, bm["Bmag"], shading="gouraud", cmap="viridis")
    ax.set_xlabel("hodograph coordinate  $q=|H|$")
    ax.set_ylabel(r"field direction  $\theta$  [rad]")
    ax.set_title("the 2-D hodograph image (a region, not a segment)\n"
                 f"linear solve residual {nl['lin_residual']:.0e}; "
                 f"Laplace-limit check {lin['laplace_error']:.0e}")
    fig.colorbar(pc, ax=ax, label="$|B|=\\mu_r q$")
    # right: the physical region the rectangle maps to (the TURNING field)
    X, Y = bm["X"], bm["Y"]
    ax2.plot(X, Y, color="0.85", lw=0.5)
    ax2.plot(X.T, Y.T, color="0.85", lw=0.5)
    pc2 = ax2.pcolormesh(X, Y, bm["theta"], shading="gouraud", cmap="twilight")
    ax2.set_aspect("equal")
    ax2.set_xlabel("physical  $x$")
    ax2.set_ylabel("physical  $y$")
    ax2.set_title("back-mapped PHYSICAL region: a turning field\n"
                  f"single-valued to closure {bm['closure']:.1e} (realisable)")
    fig.colorbar(pc2, ax=ax2, label=r"field angle $\theta$")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Chaplygin turning guide: ONE linear PDE solve on the 2-D hodograph\n")
    v = verify_solver()
    print(f"  [verify] mu_r=1 (Laplace limit): A=ln(q)*theta reproduced to "
          f"laplace_error = {v['laplace_error']:.2e}")
    print(f"           linear-solve residual = {v['lin_residual']:.2e} "
          f"(one linear solve, no Picard); theta turns over "
          f"[{v['theta_range'][0]:.1f}, {v['theta_range'][1]:.1f}] rad\n")

    r, mesh, gf = nonlinear_turn(mur0=20.0, qk=1.0)
    print(f"  [nonlinear] mu_r0={r['mur0']:.0f}: SAME turning data through the "
          f"saturating operator")
    print(f"           q in [{r['q_range'][0]:.1f}, {r['q_range'][1]:.1f}] "
          f"(spans saturation, q_k={r['qk']:.1f}); ONE linear solve "
          f"(residual {r['lin_residual']:.1e})")
    print(f"           saturation bend vs the linear harmonic = "
          f"{r['twoD_deviation']:.2e} (genuine 2-D Chaplygin content)")

    bm = back_map(mesh, gf, r["mur0"], r["qk"], *r["q_range"], r["theta_range"][1])
    print(f"           back-mapped to PHYSICAL space, single-valued to closure "
          f"{bm['closure']:.1e} (the field is realisable)\n")
    print("  => the nonlinearity is a COEFFICIENT mu(q): the turning field is a")
    print("     single LINEAR elliptic solve on the 2-D hodograph plane (NOT a")
    print("     quadrature), and it back-maps to a realisable physical TURNING")
    print("     field.  The free-boundary closure (prescribe the physical guide")
    print("     -> solve for its hodograph image) is the remaining frontier.")
    _plot(v, r, bm)


if __name__ == "__main__":
    main()
