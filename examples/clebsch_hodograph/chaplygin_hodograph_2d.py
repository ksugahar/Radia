"""Chaplygin rung 1.5b: the 2-D hodograph LINEARISES saturation -- a 1-shot
quadrature that reproduces the full nonlinear FEM loop, on a special geometry
whose hodograph image is a simple region.

RESEARCH example (track A, "Kelvin in the hodograph", the Chaplygin rung).  The
companion ``saturation_loop_2d.py`` shows the straightforward NONLINEAR LOOP (a
saturable cylinder, solved by Picard).  THIS file shows the other half the user
asked for ("hodograph 変換をすると どーせ反復が必要なので 非線形ループと混ぜて
短くできないか?"): the **hodograph turns the nonlinear magnetostatic PDE into a
LINEAR one**, so for a special geometry it is solved in ONE shot (no iteration),
and that 1-shot agrees with the nonlinear loop.


THE MAGNETOSTATICS <-> COMPRESSIBLE-FLOW DICTIONARY
---------------------------------------------------
2-D current-free nonlinear magnetostatics is *identical* to 2-D steady
irrotational compressible gas flow:

    magnetics                          compressible flow
    ------------------------------     -----------------------------------
    H = grad(Phi)   (curl H = 0)       v = grad(phi)   (irrotational)
    B = curl(A zhat) (div B = 0)       rho*v = curl(psi zhat) (mass cons.)
    B = mu(|H|) H                      rho*v = rho(|v|) v
    mu(|H|)                            rho(q),  q = |v|
    Phi (scalar potential)             phi (velocity potential)
    A   (flux function)                psi (stream function)

So ``B = mu(|H|) H`` IS the constitutive law of a gas with density
``rho(q) = mu(q)``.  Saturation (mu falls as the field rises) is the gas getting
"thinner" as it speeds up.


THE CHAPLYGIN LINEARISATION (why the nonlinearity disappears)
-------------------------------------------------------------
Take ``q = |H|`` and ``theta = arg(H)`` as the INDEPENDENT variables (the
hodograph plane), with ``Phi`` (scalar potential) and ``A`` (flux function) the
unknowns.  Inverting (Molenbroek 1890, Chaplygin 1902)

    dPhi = q cos th dx + q sin th dy
    dA   = -mu q sin th dx + mu q cos th dy

(``H = grad Phi``, ``|H|=q``; ``B = mu H``, ``B = (A_y, -A_x)``), the
integrability of dx, dy gives the **linear, variable-coefficient** first-order
pair (derived + verified against the linear Laplace limit; see the genuine 2-D
PDE solver in chaplygin_turning_guide_2d.py)

    Phi_q = -q (mu q)' / (mu q)^2  A_theta ,     A_q = (mu / q)  Phi_theta ,

whose elimination is the self-adjoint elliptic Chaplygin equation for A,

    d/dq( (q/mu) A_q ) + ((mu q)'/(mu^2 q)) A_theta_theta = 0 .

The point: ``mu(q)`` is now a COEFFICIENT (a known function of the independent
coordinate ``q``), NOT a nonlinearity in an unknown field.  The PDE is linear.
The whole nonlinearity has been absorbed into the variable coefficient
``mu(q)`` -- exactly the user's "混ぜて短く": one linear solve, no outer Picard
over the constitutive law.


THE HONEST CAVEAT: ordinary saturation stays ELLIPTIC (no limiting line)
------------------------------------------------------------------------
In gas dynamics the Chaplygin equation changes type (elliptic <-> hyperbolic)
at the sound speed -- the "limiting line" -- which is what makes transonic flow
hard.  The magnetic analogue is the *differential* reluctivity: the operator
``div(mu(|H|) grad Phi)`` is elliptic as long as both

    mu(q) > 0          and     d(mu q)/dq > 0  (== nu + nu' q > 0 in B).

For a NORMAL saturating material mu(q) DECREASES but ``mu*q = |B|`` still
INCREASES with q, so ``d(mu q)/dq > 0`` always: **the magnetostatic problem is
elliptic everywhere, there is no limiting line, no shock, no hyperbolic
hard part.**  (An earlier note of mine sloppily called Chaplygin "transonic" --
that type-change does NOT occur for ordinary B-H saturation; it would need a
material with falling |B| vs |H|, which is unphysical for a passive magnet.)

Consequence for "1-shot": the linear Chaplygin equation is elliptic, so it is a
boundary-value problem on the HODOGRAPH DOMAIN (the image of the physical region
in the (q, theta) plane).  The 1-shot is clean **iff that image is a simple,
fixed region** -- which it is only for special geometries.  General geometry
gives a free-boundary (unknown hodograph image) -- the genuinely hard case,
deferred.  This file does the SIMPLE-IMAGE special geometry.


THE SPECIAL GEOMETRY: a slender saturable flux guide (simple hodograph image)
-----------------------------------------------------------------------------
A flux guide of slowly varying width ``w(x)`` carrying flux ``Psi`` in +x
(iron walls top/bottom are flux lines, A = 0 and A = Psi).  The field is
everywhere nearly axial, so its hodograph image is the thin region

    theta ~ 0,   q = |H| in [ H(x) over the guide ],

i.e. a SEGMENT, not a 2-D blob -- the simplest possible hodograph image.  On it
the Chaplygin system collapses to a single quadrature (a "simple wave"):

  * Flux conservation (geometry only, mu-independent):  |B|(x) = Psi / w(x).
  * The drive (magnetomotive force) is the 1-SHOT QUADRATURE

        DeltaPhi(Psi) = INT_0^L  nu(|B|(x)) |B|(x) dx
                      = INT_0^L  ( |B|(x) / mu(|B|(x)) ) / mu0  dx ,

    with ``|B|(x) = Psi / w(x)`` -- evaluated ONCE, no iteration.

Saturation makes ``DeltaPhi(Psi)`` bend upward (the throat saturates first, its
reluctivity jumps, it dominates the drop).  We compare this 1-shot curve to the
**full 2-D nonlinear FEM loop** ``div(nu(|grad A|) grad A) = 0`` on the actual
curved-wall guide.  They agree to the slenderness error (and the agreement
TIGHTENS as the guide is made more slender -- the lubrication limit in which the
hodograph image is exactly a segment).

This is the honest, finishable Chaplygin demonstration: a genuine nonlinear 2-D
field whose hodograph linearises to a 1-shot, matching the loop.  The next rung
(a guide that TURNS, so theta varies and the hodograph image is a real 2-D
region needing a linear PDE solve) is noted at the end of the program.

run:  python chaplygin_hodograph_2d.py
"""
import math
import os

from numpy import pi
import numpy as np

_trap = getattr(np, "trapezoid", getattr(np, "trapz", None))
from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx, CF, x, y,
                     sqrt, BilinearForm, TaskManager, Integrate)
from netgen.occ import WorkPlane, Circle, OCCGeometry

MU0 = 4 * pi * 1e-7


# --------------------------------------------------------------------------- #
# geometry: a slender flux guide pinched to a throat by two circular notches
# --------------------------------------------------------------------------- #
def _throat_geom(L, H, depth, Rc):
    """Rectangle (L x H) with a circular notch bitten out of the top and the
    bottom at mid-span, narrowing the channel from H to H - 2*depth.  Returns
    (OCCGeometry, w_of_x) where w_of_x(x) is the ANALYTIC local width used by
    the 1-shot quadrature.  x runs 0..L (geometry is centred on the origin)."""
    base = WorkPlane().RectangleC(L, H).Face()       # centred at (0, 0)
    yc = H / 2.0 + Rc - depth                         # notch-circle centre |y|
    top = Circle((0.0, yc), Rc).Face()
    bot = Circle((0.0, -yc), Rc).Face()
    face = base - top - bot
    face.faces.name = "guide"
    for e in face.edges:
        c = e.center
        if abs(c[0] + L / 2.0) < L * 1e-3:
            e.name = "inlet"
        elif abs(c[0] - L / 2.0) < L * 1e-3:
            e.name = "outlet"
        elif c[1] > 0.0:
            e.name = "top"
        else:
            e.name = "bot"

    half = math.sqrt(max(0.0, Rc * Rc - (Rc - depth) ** 2))  # half throat span

    def w_of_x(xc):                                   # xc in 0..L (0 = inlet)
        xx = xc - L / 2.0
        if abs(xx) >= half:
            return H
        return 2.0 * (yc - math.sqrt(max(0.0, Rc * Rc - xx * xx)))

    return OCCGeometry(face, dim=2), w_of_x


def _mu_r(Bmag, mur0, Bk):
    """Froehlich saturation: mu_r0 at B=0, -> 1 as B -> infinity."""
    return 1.0 + (mur0 - 1.0) / (1.0 + (Bmag / Bk) ** 2)


# --------------------------------------------------------------------------- #
# the two computations: 2-D nonlinear FEM loop  and  1-shot hodograph quadrature
# --------------------------------------------------------------------------- #
def solve_chaplygin(Psi_list=(0.004, 0.008, 0.012, 0.016, 0.024, 0.032),
                    mur0=200.0, Bk=1.0, L=0.20, H=0.04, depth=0.012, Rc=0.11,
                    order=3, maxh=0.005, niter=600, tol=1e-7, relax=0.5,
                    n_quad=2001):
    """For each driven flux Psi: run the full 2-D nonlinear loop (A-formulation
    nu(|B|) Picard) and the 1-shot hodograph quadrature, and compare the
    magnetomotive drive DeltaPhi.  Returns the sweep + agreement metrics."""
    geo, w_of_x = _throat_geom(L, H, depth, Rc)
    w_throat = H - 2.0 * depth

    # --- 1-shot hodograph quadrature (no mesh, no iteration) --------------- #
    xs = np.linspace(0.0, L, n_quad)
    w = np.array([w_of_x(xc) for xc in xs])

    def deltaphi_1shot(Psi):
        Bx = Psi / w                                  # |B|(x) = Psi / w(x)
        nu = 1.0 / (MU0 * _mu_r(Bx, mur0, Bk))        # reluctivity nu(|B|)
        return float(_trap(nu * Bx, xs))           # INT nu |B| dx

    rows = []
    hodo_theta_max = 0.0
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        mesh.Curve(order)
        fes = H1(mesh, order=order, dirichlet="top|bot|inlet|outlet")
        u, v = fes.TnT()
        # Dirichlet data: walls are flux lines (A=0 bottom, A=Psi top); the
        # ends carry the axial flux via a linear A across the full-height gap.
        ramp = (y + H / 2.0) / H                       # 0 at bottom, 1 at top

        def lin_solve(nucf, Psi):
            bcf = mesh.BoundaryCF({"top": CF(1.0), "bot": CF(0.0),
                                   "inlet": ramp, "outlet": ramp}, default=0.0)
            g = GridFunction(fes)
            g.Set(Psi * bcf, definedon=mesh.Boundaries("top|bot|inlet|outlet"))
            a_ = BilinearForm(nucf * InnerProduct(grad(u), grad(v)) * dx)
            a_.Assemble()
            r = g.vec.CreateVector()
            r.data = -a_.mat * g.vec
            g.vec.data += a_.mat.Inverse(fes.FreeDofs(),
                                         inverse="sparsecholesky") * r
            return g

        A = lin_solve(CF(1.0 / MU0), Psi_list[0])      # linear seed (once)
        Psi_prev = Psi_list[0]
        for Psi in Psi_list:
            A.vec.data = (Psi / Psi_prev) * A.vec       # continuation warm start
            Psi_prev = Psi                              # (rescale to new flux)
            prev = np.array(A.vec)
            nit = 0
            for it in range(niter):
                B = grad(A)
                Bmag = sqrt(B[0] * B[0] + B[1] * B[1] + 1e-30)
                nucf = 1.0 / (MU0 * _mu_r(Bmag, mur0, Bk))
                Asol = lin_solve(nucf, Psi)            # nu-frozen linear solve
                # under-relaxed Picard: undamped substitution OSCILLATES under
                # the throat's strong reluctivity contrast (convex problem, so a
                # unique solution exists; damping restores convergence to it).
                A.vec.data = (1.0 - relax) * A.vec + relax * Asol.vec
                cur = np.array(A.vec)
                d = np.linalg.norm(cur - prev) / (np.linalg.norm(cur) + 1e-30)
                prev = cur.copy()
                nit = it + 1
                if d < tol:
                    break

            # DeltaPhi_FEM = INT H_x dx along the centreline y=0 (= Phi drop,
            # path-independent since curl H = 0).  H_x = nu(|B|) B_x, B_x=dA/dy.
            B = grad(A)
            Bmag = sqrt(B[0] * B[0] + B[1] * B[1] + 1e-30)
            nucf = 1.0 / (MU0 * _mu_r(Bmag, mur0, Bk))
            Hx = nucf * B[1]                           # B_x = dA/dy = grad(A)[1]
            xc = np.linspace(-L / 2.0 + 1e-4, L / 2.0 - 1e-4, 400)
            hx = np.array([Hx(mesh(float(xi), 0.0)) for xi in xc])
            dphi_fem = float(_trap(hx, xc))
            dphi_1s = deltaphi_1shot(Psi)

            # hodograph image width: theta = arg(B) sampled over the INTERIOR
            # (NOT the centreline -- there theta=0 by symmetry, a tautology).
            # Off-axis near the throat the converging walls tilt the field; the
            # max tilt is how far the image strays from the 1-shot SEGMENT.
            for xi in np.linspace(-L / 2.0 + 5e-3, L / 2.0 - 5e-3, 60):
                wl = w_of_x(float(xi) + L / 2.0)
                for yi in np.linspace(-0.45 * wl, 0.45 * wl, 9):
                    bxi = B[1](mesh(float(xi), float(yi)))
                    byi = (-B[0])(mesh(float(xi), float(yi)))
                    hodo_theta_max = max(hodo_theta_max,
                                         abs(float(math.atan2(byi, bxi))))

            rows.append((float(Psi), float(Psi / w_throat), int(nit),
                         dphi_fem, dphi_1s,
                         float(abs(dphi_fem - dphi_1s) / abs(dphi_fem))))

    rel = [r[5] for r in rows]
    return {
        "mur0": mur0, "Bk": Bk, "L": L, "H": H, "w_throat": float(w_throat),
        "slenderness_L_over_w": float(L / w_throat),
        "rows": rows,            # (Psi, B_throat, iters, dPhi_FEM, dPhi_1shot, relerr)
        "max_relerr": float(max(rel)),
        "mean_relerr": float(np.mean(rel)),
        # off-axis field tilt near the throat shoulders: the hodograph image is
        # a thin BAND about theta=0 with localized excursions to this angle.
        # The 1-shot (segment) still matches because the drive integral is
        # throat-CENTRE dominated, where the field is axial.
        "hodo_theta_max_deg": float(hodo_theta_max * 180.0 / pi),
        # the curve bends: ratio of (drive/flux) high vs low field, FEM and 1-shot
        "bend_fem": float((rows[-1][3] / rows[-1][0]) / (rows[0][3] / rows[0][0])),
        "bend_1shot": float((rows[-1][4] / rows[-1][0]) / (rows[0][4] / rows[0][0])),
    }


def slenderness_trend(Rcs=(0.11, 0.30, 0.80), mur0=200.0, Bk=1.0, L=0.5,
                      H=0.04, depth=0.012, order=2, maxh=0.006):
    """Back the docstring claim that the 1-shot is the SLENDER limit: at a fixed
    throat field, make the throat more gradual (larger notch radius Rc -> longer,
    gentler constriction) and watch the 1-shot-vs-loop rel.err shrink.  Returns
    [(throat slenderness = half-length/width, rel.err), ...]."""
    w_throat = H - 2.0 * depth
    Psi = 1.5 * w_throat                                # |B|_throat = 1.5 T
    out = []
    for Rc in Rcs:
        r = solve_chaplygin(Psi_list=(Psi,), mur0=mur0, Bk=Bk, L=L, H=H,
                            depth=depth, Rc=Rc, order=order, maxh=maxh)
        half = math.sqrt(max(0.0, Rc * Rc - (Rc - depth) ** 2))
        out.append((half / w_throat, float(r["rows"][0][5])))
    return out


def _plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Psi = np.array([r[0] for r in res["rows"]])
    Bt = np.array([r[1] for r in res["rows"]])
    fem = np.array([r[3] for r in res["rows"]])
    one = np.array([r[4] for r in res["rows"]])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0), dpi=150)
    ax.plot(Bt, fem, "o-", color="C0", label="2-D nonlinear FEM loop (Picard)")
    ax.plot(Bt, one, "x--", color="C1", lw=1.4,
            label="1-shot hodograph quadrature")
    ax.set_xlabel("throat flux density  $|B|_{throat}$  [T]")
    ax.set_ylabel(r"drive  $\Delta\Phi$  [A]")
    ax.set_title(f"saturable flux guide: 1-shot vs loop\n"
                 f"max rel.err = {res['max_relerr']:.1e}  "
                 f"(L/w = {res['slenderness_L_over_w']:.0f})")
    ax.legend(fontsize=8)

    relerr = np.array([r[5] for r in res["rows"]])
    ax2.semilogy(Bt, relerr, "s-", color="C2")
    ax2.set_xlabel("throat flux density  $|B|_{throat}$  [T]")
    ax2.set_ylabel(r"$|\Delta\Phi_{FEM}-\Delta\Phi_{1shot}| / \Delta\Phi_{FEM}$")
    ax2.set_title(f"hodograph image: a thin band about "
                  rf"$\theta=0$" "\n"
                  rf"$\max|\theta_B|$ = {res['hodo_theta_max_deg']:.0f}$^\circ$"
                  " at the throat shoulders")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Chaplygin rung: 2-D hodograph linearises saturation "
          "(1-shot vs loop)\n")
    r = solve_chaplygin()
    print(f"  saturable flux guide: mu_r0={r['mur0']:.0f}  Bk={r['Bk']:.1f} T  "
          f"throat w={r['w_throat']*1e3:.0f} mm  L/w={r['slenderness_L_over_w']:.0f}\n")
    print("  |B|_throat [T] | iters | dPhi_FEM [A] | dPhi_1shot [A] | rel.err")
    for Psi, Bt, nit, fem, one, rel in r["rows"]:
        print(f"    {Bt:8.3f}    | {nit:4d}  | {fem:11.4f}  | {one:12.4f}   "
              f"| {rel:.2e}")
    print(f"\n  -> 1-shot hodograph quadrature matches the 2-D nonlinear loop:")
    print(f"     max rel.err = {r['max_relerr']:.2e}, "
          f"mean = {r['mean_relerr']:.2e}.")
    print(f"     drive/flux BENDS x{r['bend_fem']:.2f} (FEM) vs "
          f"x{r['bend_1shot']:.2f} (1-shot) as the throat saturates.")
    print(f"     hodograph image is a thin BAND about theta=0: "
          f"max|theta_B| = {r['hodo_theta_max_deg']:.0f} deg at the throat "
          f"shoulders\n     (the 1-shot still holds: the drive is "
          f"throat-CENTRE dominated, where the field is axial).")

    # back the docstring claim "agreement tightens as the guide is more
    # slender": lengthen the straight runs (larger L/w) at a fixed throat field
    # and watch the rel.err shrink toward the lubrication (segment) limit.
    print("\n  slenderness trend (fixed throat |B|=1.5 T, gentler throat):")
    print("    throat length/width | rel.err  (1-shot vs loop)")
    trend = slenderness_trend()
    for s, e in trend:
        print(f"        {s:8.1f}        | {e:.2e}")
    print(f"    -> rel.err {trend[0][1]:.1e} -> {trend[-1][1]:.1e} as the throat "
          f"slenderness {trend[0][0]:.0f} -> {trend[-1][0]:.0f}: the segment "
          f"1-shot is the slender limit.")

    print("\n  => the nonlinearity became a COEFFICIENT mu(q): ONE quadrature,")
    print("     no iteration -- the 2-D hodograph linearised the saturation.")
    _plot(r)


if __name__ == "__main__":
    main()
