r"""Self-contained compute for hodograph_feasibility_2d.ipynb (docs/clebsch_hodograph).

2D LINEAR hodograph feasibility of a BENDING MAGNET.  Imported by the notebook and by the
golden tests/feec/test_hodograph_feasibility_2d.py.

Setup (2D, current-free air, constant mu -> Laplace).  Demand the mid-plane field
    B_y(x,0) = g(x) = [tanh((x+x0)/d) - tanh((x-x0)/d)] / (2 tanh(x0/d))   (g(0)=1),
a flat-top of half-width ~x0 with an edge of width d.  The scalar potential is HARMONIC in
the gap, so g must extend upward as its unique smooth (analytic) continuation.  That
continuation has a CLOSED FORM (Omega(z) = (i d/2)[logcosh((z+x0)/d) - logcosh((z-x0)/d)] /
tanh(x0/d), phi = Re Omega):

    phi(x,y) = (d / (2 tanh(x0/d))) [ atan2(sinh a-  sin b, cosh a-  cos b)
                                     - atan2(sinh a+  sin b, cosh a+  cos b) ],
    a+=(x+x0)/d, a-=(x-x0)/d, b=y/d,   with  -d phi/dy|_{y=0} = g(x),  phi(x,0)=0.

tanh's nearest singularity is at height  y_sing = pi d / 2, so the field is realizable by an
iron pole at gap h ONLY IF  y_sing > h, i.e.

    d > d* = (2/pi) h ~ 0.64 h        -- the field edge is no sharper than ~0.64 x gap.

The iron pole is the equipotential {phi = phi0}; below d* the continuation is singular inside
the gap and no single smooth pole can produce the field.  NOTE the (2/pi) factor is specific
to the tanh edge; what is UNIVERSAL is the gap scaling (edge sharpness ~ gap).  This bound is
a HARMONIC-CONTINUATION / analyticity fact (Cauchy-Kovalevskaya flavour), NOT a hodograph
result and NOT a "limit line": ordinary passive saturation stays ELLIPTIC (mu falls but
|B|=mu*|H| still rises, d(mu*q)/dq>0), so there is no type change / shock for this problem.

The genuinely hodograph-specific content is the PARTIAL von-Mises transform (the Sugahara-lab
hodograph: keep one coordinate, transform ONE potential).  Here keep x and replace y by the
scalar-potential coordinate  s = -phi  (monotone-increasing in y since B_y = -phi_y > 0 across
the gap: s=0 on the mid-plane, s=s0 at the pole).  The unknown iron-pole SHAPE -- a free
boundary in physical space -- becomes the FIXED top edge {s=s0} of the rectangle [-xr,xr]x[0,s0]
in the (x,s) chart; the physical map y(x,s) solves the von-Mises PDE

    y_s^2 y_xx - 2 y_x y_s y_xs + (1 + y_x^2) y_ss = 0            (mu = const, linear)

on that fixed rectangle, with a strictly positive Jacobian y_s>0 (single-valued, no fold).
This is the "design in field space" transparency: the pole is a coordinate line, and saturation
would enter as a coefficient mu(q) on the SAME fixed chart (no new free boundary).  vonmises_chart()
verifies the PDE residual (-> 0 under refinement) and the fold-free Jacobian on the closed form.
"""
import numpy as np

X0, H, D = 2.0, 1.0, 0.9        # good-field half-width, gap, a feasible edge width
YT = 1.15                       # rectangle height for the FEM check (< y_sing = pi D/2)


def _norm(d):
    return np.tanh(X0 / d)


def g(x, d=D):
    x = np.asarray(x, dtype=float)
    return (np.tanh((x + X0) / d) - np.tanh((x - X0) / d)) / (2.0 * _norm(d))


def phi_np(x, y, d=D):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    ap, am, b = (x + X0) / d, (x - X0) / d, y / d
    raw = (np.arctan2(np.sinh(am) * np.sin(b), np.cosh(am) * np.cos(b))
           - np.arctan2(np.sinh(ap) * np.sin(b), np.cosh(ap) * np.cos(b)))
    return 0.5 * d * raw / _norm(d)


def y_sing(d):
    return np.pi * d / 2.0            # nearest singularity height of the continued field


def d_star(h=H):
    return 2.0 * h / np.pi           # feasibility boundary: d > d* is realizable


def feasibility_table(ds, h=H):
    return [{"d": float(d), "y_sing": float(y_sing(d)), "feasible": bool(y_sing(d) > h)}
            for d in ds]


def analytic_pole(xg, phi0, d=D, ytop=None):
    ytop = YT if ytop is None else ytop
    yy = np.linspace(0, ytop, 2500)
    yp = np.full(len(xg), np.nan)
    for j, x in enumerate(xg):
        s = phi_np(x, yy, d) - phi0
        k = np.where(np.sign(s[:-1]) != np.sign(s[1:]))[0]
        if len(k):
            i = k[0]; t = -s[i] / (s[i + 1] - s[i]); yp[j] = yy[i] + t * (yy[1] - yy[0])
    return yp


def fem_verify(d=D, order=4, maxh=0.06, w=3.5):
    """Independent ngsolve linear-FEM (manufactured-solution) check on a rectangle BELOW the
    singularity: impose the analytic phi on the whole boundary, solve Laplace, and confirm
    (a) interior B_y(x,0)=g(x), (b) FEM phi = analytic phi, (c) FEM equipotential = read-off
    pole.  Returns metrics + arrays for plotting."""
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry
    from ngsolve import x as X, y as Y, sinh, cosh, sin, cos, atan2

    phi0 = float(phi_np(0.0, H, d))
    ap, am, b = (X + X0) / d, (X - X0) / d, Y / d
    phi_cf = (0.5 * d / float(_norm(d))) * (atan2(sinh(am) * sin(b), cosh(am) * cos(b))
                                            - atan2(sinh(ap) * sin(b), cosh(ap) * cos(b)))
    face = WorkPlane().MoveTo(-w, 0).Rectangle(2 * w, YT).Face()
    for e in face.edges:
        e.name = "bnd"
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=maxh))
        fes = ng.H1(mesh, order=order, dirichlet="bnd")
        u, v = fes.TnT()
        a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
        gfu = ng.GridFunction(fes)
        gfu.Set(phi_cf, definedon=mesh.Boundaries("bnd"))
        r = gfu.vec.CreateVector(); r.data = -a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
        grad = ng.grad(gfu)
        xs = np.linspace(-2.8, 2.8, 141)
        By = np.array([-grad(mesh(float(x), 0.02))[1] for x in xs])
        phi_l2 = float(ng.sqrt(ng.Integrate((gfu - phi_cf) ** 2, mesh)
                               / ng.Integrate(phi_cf ** 2, mesh)))
        xp = np.linspace(-2.6, 2.6, 61)
        yy = np.linspace(0.02, YT - 0.02, 300)
        yp_fem = np.full(len(xp), np.nan)
        for j, xx in enumerate(xp):
            col = np.array([gfu(mesh(float(xx), float(t))) for t in yy]) - phi0
            k = np.where(np.sign(col[:-1]) != np.sign(col[1:]))[0]
            if len(k):
                i = k[0]; t = -col[i] / (col[i + 1] - col[i]); yp_fem[j] = yy[i] + t * (yy[1] - yy[0])
    By_t = g(xs, d)
    rel = float(np.max(np.abs(By - By_t)) / np.max(By_t))
    yp_an = analytic_pole(xp, phi0, d)
    m = np.isfinite(yp_fem) & np.isfinite(yp_an)
    pole_err = float(np.max(np.abs(yp_fem[m] - yp_an[m]))) if m.any() else float("nan")
    return {
        "d": float(d), "gap_h": float(H), "phi0": phi0,
        "d_star": float(d_star(H)), "feasible": bool(y_sing(d) > H),
        "By_rel_err": rel, "phi_L2_err": phi_l2, "pole_match_err": pole_err,
        "pole_gap_at_x0": float(yp_an[np.argmin(np.abs(xp))]),
        "_xs": xs, "_By_fem": By, "_By_target": By_t,
        "_xp": xp, "_yp_fem": yp_fem, "_yp_an": yp_an,
    }


def vonmises_chart(d=D, xr=1.5, nX=161, nS=121, ytop=1.25, ny=80001):
    """Partial von-Mises hodograph (keep x, transform ONE potential): the (x, s) chart with
    s = -phi (monotone-increasing scalar potential; s=0 on the mid-plane, s=s0 at the pole).

    Invert phi(x,.) = -s on a fine y-grid to get the physical height y(x,s).  The unknown
    iron-pole shape -- a FREE BOUNDARY in physical (x,y) space -- becomes the FIXED top edge
    {s=s0} of the rectangle [-xr,xr] x [0,s0].  Verifies, on the closed-form map:
      - pde_resid : relative residual of the von-Mises PDE
                    y_s^2 y_xx - 2 y_x y_s y_xs + (1+y_x^2) y_ss = 0   (-> 0 under refinement)
      - jac_min   : min of y_s  (>0 => single-valued diffeomorphism, no fold)
      - pole_match: top edge y(x,s0) vs the equipotential read-off {phi = -s0}
    """
    s0 = -float(phi_np(0.0, H, d))
    xs = np.linspace(-xr, xr, nX)
    Ss = np.linspace(0.0, s0, nS)
    Y = np.zeros((nX, nS))
    for i, xx in enumerate(xs):
        yy = np.linspace(0.0, ytop, ny)
        sv = -phi_np(np.full_like(yy, xx), yy, d)          # s = -phi, increasing in y
        Y[i, :] = np.interp(Ss, sv, yy)
    dx, dS = xs[1] - xs[0], Ss[1] - Ss[0]
    yX = np.gradient(Y, dx, axis=0); yS = np.gradient(Y, dS, axis=1)
    yXX = np.gradient(yX, dx, axis=0); ySS = np.gradient(yS, dS, axis=1)
    yXS = np.gradient(yX, dS, axis=1)
    r = yS ** 2 * yXX - 2.0 * yX * yS * yXS + (1.0 + yX ** 2) * ySS
    sc = yS ** 2 * np.abs(yXX) + 2.0 * np.abs(yX * yS * yXS) + (1.0 + yX ** 2) * np.abs(ySS)
    rel = np.abs(r) / np.maximum(sc, 1e-30)
    core = rel[2:-2, 2:-2]                                  # drop the 2-cell FD border
    poleh = Y[:, -1]
    yp = analytic_pole(xs, -s0, d, ytop=ytop)
    m = np.isfinite(yp)
    return {
        "d": float(d), "s0": float(s0), "xr": float(xr), "nX": int(nX), "nS": int(nS),
        "pde_resid": float(core.max()), "pde_resid_median": float(np.median(core)),
        "jac_min": float(yS.min()), "single_valued": bool(yS.min() > 0.0),
        "pole_match": float(np.max(np.abs(poleh[m] - yp[m]))) if m.any() else float("nan"),
        "_xs": xs, "_Ss": Ss, "_Y": Y, "_poleh": poleh,
    }


def figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, RED, GREEN, AMBER = "#1f6feb", "#d1495b", "#2e8b57", "#e0851e"
    fig, ax = plt.subplots(1, 3, figsize=(16.2, 4.7), dpi=140)

    # (1) feasibility boundary
    dd = np.linspace(0.15, 1.6, 300); ds = d_star(H)
    ax[0].plot(dd, y_sing(dd), color=BLUE, lw=2.6, label=r"blow-up height $\pi d/2$")
    ax[0].axhline(H, color="#41506a", ls="--", lw=1.4, label="gap h")
    ax[0].axvspan(0.15, ds, color=RED, alpha=0.10)
    ax[0].plot([ds], [H], marker="*", ms=18, color=AMBER, mec="#1a2230", zorder=5)
    ax[0].text(0.32, 2.05, "INFEASIBLE", color=RED, fontsize=12, fontweight="bold")
    ax[0].text(1.05, 0.4, "feasible", color=GREEN, fontsize=12, fontweight="bold")
    ax[0].annotate(fr"$d^*=(2/\pi)h\approx{ds:.2f}$", (ds, H), (ds + 0.14, 1.6), fontsize=11.5,
                   arrowprops=dict(arrowstyle="->", lw=1.1))
    ax[0].set_xlabel("demanded edge width d"); ax[0].set_ylabel("height in the gap")
    ax[0].set_xlim(0.15, 1.6); ax[0].set_ylim(0, 2.5); ax[0].legend(fontsize=9.5, loc="lower right")
    ax[0].set_title("Feasibility: edge no sharper than (2/pi) x gap")

    # (2) read-off pole vs FEM equipotential (they overlap)
    ax[1].plot(res["_xp"], res["_yp_an"], color=BLUE, lw=3.0, label="pole read off (hodograph)")
    ax[1].plot(res["_xp"], res["_yp_fem"], color=AMBER, lw=1.4, ls="--",
               label="FEM equipotential")
    ax[1].plot([-2.8, 2.8], [0, 0], color="#1a2230", lw=1.6)
    ax[1].axhline(H, color="#9aa6b6", ls=":", lw=1.0)
    ax[1].set_xlabel("x"); ax[1].set_ylabel("y (gap)")
    ax[1].set_xlim(-2.8, 2.8); ax[1].set_ylim(0, 1.25)
    ax[1].set_title(f"Inverse design: pole = equipotential\n(FEM matches to {res['pole_match_err']:.0e})")
    ax[1].legend(fontsize=9.5, loc="upper center")

    # (3) FEM B_y(x,0) vs the demanded g(x)
    ax[2].plot(res["_xs"], res["_By_target"], color=BLUE, lw=3.0, label="demanded g(x)")
    ax[2].plot(res["_xs"], res["_By_fem"], color=RED, lw=1.3, ls="--", label="FEM B_y(x,0)")
    ax[2].set_xlabel("x"); ax[2].set_ylabel(r"$B_y(x,0)/B_0$")
    ax[2].set_xlim(-4, 4); ax[2].set_ylim(0, 1.15)
    ax[2].set_title(f"FEM verification\n(reproduces the demand to {res['By_rel_err']*100:.2f}% of B0)")
    ax[2].legend(fontsize=9.5, loc="upper right")

    fig.tight_layout(); fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    import os
    return os.path.abspath(path)


def figure_vonmises(vm, path):
    """The partial von-Mises hodograph: the curved physical gap (free-boundary iron pole)
    re-charted onto a FIXED rectangle where the pole is the top coordinate line s=s0."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, RED, GREEN, AMBER = "#1f6feb", "#d1495b", "#2e8b57", "#e0851e"
    xs, Ss, Y, poleh, s0 = vm["_xs"], vm["_Ss"], vm["_Y"], vm["_poleh"], vm["s0"]
    lv = [int(round(f * (len(Ss) - 1))) for f in (0.25, 0.5, 0.75)]   # interior s-levels
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.9), dpi=140)

    # (left) physical gap: mid-plane, curved pole (free boundary), curved s-level lines
    ax[0].fill_between(xs, 0, poleh, color=BLUE, alpha=0.06)
    for k in lv:
        ax[0].plot(xs, Y[:, k], color="#9aa6b6", lw=1.0)
    ax[0].plot(xs, np.zeros_like(xs), color="#1a2230", lw=1.8, label="mid-plane  s=0")
    ax[0].plot(xs, poleh, color=RED, lw=3.0, label="iron pole  s=s0 (free boundary)")
    ax[0].set_xlabel("x"); ax[0].set_ylabel("y (physical gap)")
    ax[0].set_xlim(xs.min(), xs.max()); ax[0].set_ylim(0, poleh.max() * 1.12)
    ax[0].legend(fontsize=9.5, loc="lower center")
    ax[0].set_title("Physical space: the pole shape is an unknown curve")

    # (right) the (x, s) chart: fixed rectangle, pole = straight top edge
    ax[1].fill_between([xs.min(), xs.max()], 0, s0, color=GREEN, alpha=0.06)
    for k in lv:
        ax[1].plot([xs.min(), xs.max()], [Ss[k], Ss[k]], color="#9aa6b6", lw=1.0)
    ax[1].plot([xs.min(), xs.max()], [0, 0], color="#1a2230", lw=1.8, label="mid-plane  s=0")
    ax[1].plot([xs.min(), xs.max()], [s0, s0], color=RED, lw=3.0, label="iron pole  s=s0 (fixed edge)")
    for xv in (xs.min(), xs.max()):
        ax[1].plot([xv, xv], [0, s0], color="#41506a", lw=1.2, ls="--")
    ax[1].set_xlabel("x"); ax[1].set_ylabel(r"scalar potential  $s=-\phi$")
    ax[1].set_xlim(xs.min() * 1.05, xs.max() * 1.05); ax[1].set_ylim(-0.03 * s0, 1.12 * s0)
    ax[1].legend(fontsize=9.5, loc="lower center")
    ax[1].set_title(f"(x,s) chart: free boundary -> fixed edge\n"
                    f"(PDE resid {vm['pde_resid']:.1e}, Jacobian y_s>0 min {vm['jac_min']:.2f})")

    fig.tight_layout(); fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    import os
    return os.path.abspath(path)


if __name__ == "__main__":
    tab = feasibility_table([1.4, 1.0, 0.8, d_star(), 0.5, 0.35])
    for r in tab:
        print(f"  d={r['d']:.3f}: y_sing={r['y_sing']:.3f} -> "
              f"{'FEASIBLE' if r['feasible'] else 'INFEASIBLE'}")
    res = fem_verify()
    print(f"\nFEM verify (d={res['d']}, gap {res['gap_h']}): "
          f"By rel err {res['By_rel_err']*100:.3f}%, phi L2 {res['phi_L2_err']:.1e}, "
          f"pole match {res['pole_match_err']:.1e}")
    vm = vonmises_chart()
    print(f"\nvon Mises chart (d={vm['d']}, s0={vm['s0']:.3f}, {vm['nX']}x{vm['nS']}): "
          f"PDE resid {vm['pde_resid']:.2e}, Jacobian y_s>0 {vm['single_valued']} "
          f"(min {vm['jac_min']:.3f}), pole match {vm['pole_match']:.1e}")
