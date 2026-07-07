r"""2D LINEAR hodograph model + end-shaping DESIGN of a BENDING-MAGNET END (s-y plane).

Imported by hodograph_bending_sy.ipynb and by the golden
tests/feec/test_hodograph_bending_sy.py.  s = beam direction, y = gap (half-model,
mid-plane symmetry at y=0, half-gap h).  Two parts:

PART A -- the fringe model + feasibility (hodograph INVERSE, box-free).
  Forward 2D s-y FEM is an OPEN-boundary problem: the flux return is out-of-plane (the x-y
  yoke), so the field escapes and the effective-length integral is LOG-DIVERGENT (verified:
  a phi=0 air box gives a box-DEPENDENT fringe).  The inverse formulation sidesteps this:
  demand the longitudinal fringe and read off the field, box-free.
    B_y(s,0) = g(s) = (1 - tanh((s-s0)/d)) / 2         (flat B0 inside, -> 0 outside; edge d)
    g(w) = (1 - tanh((w-s0)/d)) / 2,  w=s+iy  =>  B_y=Re g,  B_s=Im g.
    phi(s,y) = -1/2 [ y - d*atan2(sinh a sin b, cosh a cos b) ],  a=(s-s0)/d, b=y/d,
  with -d phi/dy|_{y=0} = -g  =>  B_y(s,0)=g(s),  phi(s,0)=0.
  Nearest singularity y_sing = pi d/2  =>  realizable by an iron pole at gap h iff y_sing>h:
    d > d* = (2/pi) h ~ 0.64 h    -- the LONGITUDINAL fringe (Enge edge) is no sharper than
    ~0.64 x gap (caps the effective-field-boundary sharpness / edge focusing).  A
    Cauchy-Kovalevskaya analyticity bound, NOT a limit line (passive iron stays elliptic).
  NOTE the pole END is NOT the equipotential of a demanded fringe -- {phi=-h} curls up INTO
  the demand's singularity (field diverges).  The pole FACE (deep inside) is an equipotential;
  the pole END is a FREE termination -> designed in PART B.

PART B -- end-shaping DESIGN (forward FEM optimization).
  A SQUARE pole end has a reentrant 270-deg air corner -> |B| ~ r^{-1/3} (a saturation hot
  spot).  Optimize the end CHAMFER (a forward quarter-ellipse, a=chamfer length, b=rise) to
  MINIMIZE the peak pole-face |B|.  The optimum drives peak -> B0 (no enhancement), numerically
  RECOVERING the classic Rogowski electrode without conformal algebra.  The peak pole-face field
  is a LOCAL quantity -> box/mesh-CONVERGENT (unlike the log-divergent effective length).
"""
import numpy as np

S0, H, D = 0.0, 1.0, 0.9        # end location, gap, feasible fringe width
YT = 1.15                       # FEM rectangle height (< y_sing = pi D/2 = 1.414)


# ------------------------------------------------------------------ PART A: fringe model
def g(s, d=D, s0=S0):
    s = np.asarray(s, dtype=float)
    return 0.5 * (1.0 - np.tanh((s - s0) / d))


def field_complex(s, y, d=D, s0=S0):
    w = np.asarray(s, dtype=float) + 1j * np.asarray(y, dtype=float)
    gw = 0.5 * (1.0 - np.tanh((w - s0) / d))
    return gw.real, gw.imag                              # (B_y, B_s)


def phi_np(s, y, d=D, s0=S0):
    s = np.asarray(s, dtype=float); y = np.asarray(y, dtype=float)
    a, b = (s - s0) / d, y / d
    return -0.5 * (y - d * np.arctan2(np.sinh(a) * np.sin(b), np.cosh(a) * np.cos(b)))


def y_sing(d):
    return np.pi * d / 2.0


def d_star(h=H):
    return 2.0 * h / np.pi


def feasibility_table(ds, h=H):
    return [{"d": float(d), "y_sing": float(y_sing(d)), "feasible": bool(y_sing(d) > h)}
            for d in ds]


def pole_profile(sg, phi0=None, d=D, s0=S0, ytop=None):
    phi0 = -H if phi0 is None else phi0
    ytop = YT if ytop is None else ytop
    yy = np.linspace(0.0, ytop, 3000)
    yp = np.full(len(sg), np.nan)
    for j, s in enumerate(sg):
        f = phi_np(np.full_like(yy, s), yy, d, s0) - phi0
        k = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
        if len(k):
            i = k[0]; tt = -f[i] / (f[i + 1] - f[i]); yp[j] = yy[i] + tt * (yy[1] - yy[0])
    return yp


def fem_verify_sy(d=D, s0=S0, order=4, maxh=0.05, w=3.5):
    """Box-free manufactured-solution FEM check on a rectangle below the singularity."""
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry
    from ngsolve import x as Sx, y as Yc, sinh, cosh, sin, cos, atan2
    phi0 = -H
    a, b = (Sx - s0) / d, Yc / d
    phi_cf = -0.5 * (Yc - d * atan2(sinh(a) * sin(b), cosh(a) * cos(b)))
    face = WorkPlane().MoveTo(s0 - w, 0).Rectangle(2 * w, YT).Face()
    for e in face.edges:
        e.name = "bnd"
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=maxh))
        fes = ng.H1(mesh, order=order, dirichlet="bnd")
        u, v = fes.TnT()
        A = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); A.Assemble()
        gfu = ng.GridFunction(fes)
        gfu.Set(phi_cf, definedon=mesh.Boundaries("bnd"))
        r = gfu.vec.CreateVector(); r.data = -A.mat * gfu.vec
        gfu.vec.data += A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
        grad = ng.grad(gfu)
        ssv = np.linspace(s0 - 2.8, s0 + 2.8, 141)
        By = np.array([-grad(mesh(float(s), 0.02))[1] for s in ssv])
        phi_l2 = float(ng.sqrt(ng.Integrate((gfu - phi_cf) ** 2, mesh)
                               / ng.Integrate(phi_cf ** 2, mesh)))
        sp = np.linspace(s0 - 3.0, s0 - 0.6, 61)
        yy = np.linspace(0.02, YT - 0.02, 300)
        yp_fem = np.full(len(sp), np.nan)
        for j, sv in enumerate(sp):
            col = np.array([gfu(mesh(float(sv), float(tt))) for tt in yy]) - phi0
            k = np.where(np.sign(col[:-1]) != np.sign(col[1:]))[0]
            if len(k):
                i = k[0]; tt = -col[i] / (col[i + 1] - col[i]); yp_fem[j] = yy[i] + tt * (yy[1] - yy[0])
    By_t = g(ssv, d, s0)
    rel = float(np.max(np.abs(By - By_t)) / np.max(By_t))
    yp_an = pole_profile(sp, phi0, d, s0)
    m = np.isfinite(yp_fem) & np.isfinite(yp_an)
    pole_err = float(np.max(np.abs(yp_fem[m] - yp_an[m]))) if m.any() else float("nan")
    return {
        "d": float(d), "s0": float(s0), "gap_h": float(H), "phi0": phi0,
        "d_star": float(d_star(H)), "feasible": bool(y_sing(d) > H),
        "By_rel_err": rel, "phi_L2_err": phi_l2, "pole_match_err": pole_err,
        "pole_gap_deep": float(yp_an[0]),
        "_ssv": ssv, "_By_fem": By, "_By_target": By_t,
        "_sp": sp, "_yp_fem": yp_fem, "_yp_an": yp_an,
    }


def figure_fringe(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, RED, GREEN, AMBER = "#1f6feb", "#d1495b", "#2e8b57", "#e0851e"
    d, s0 = res["d"], res["s0"]
    fig, ax = plt.subplots(1, 3, figsize=(16.4, 4.7), dpi=140)
    dd = np.linspace(0.15, 1.6, 300); ds = d_star(H)
    ax[0].plot(dd, y_sing(dd), color=BLUE, lw=2.6, label=r"blow-up height $\pi d/2$")
    ax[0].axhline(H, color="#41506a", ls="--", lw=1.4, label="gap h")
    ax[0].axvspan(0.15, ds, color=RED, alpha=0.10)
    ax[0].plot([ds], [H], marker="*", ms=18, color=AMBER, mec="#1a2230", zorder=5)
    ax[0].text(0.30, 2.05, "INFEASIBLE", color=RED, fontsize=12, fontweight="bold")
    ax[0].text(1.05, 0.4, "feasible", color=GREEN, fontsize=12, fontweight="bold")
    ax[0].annotate(fr"$d^*=(2/\pi)h\approx{ds:.2f}$", (ds, H), (ds + 0.14, 1.6), fontsize=11.5,
                   arrowprops=dict(arrowstyle="->", lw=1.1))
    ax[0].set_xlabel("fringe (edge) width d"); ax[0].set_ylabel("height in the gap")
    ax[0].set_xlim(0.15, 1.6); ax[0].set_ylim(0, 2.5); ax[0].legend(fontsize=9.5, loc="lower right")
    ax[0].set_title("Feasibility: fringe no sharper than (2/pi) x gap")
    sp, yp_an, yp_fem = res["_sp"], res["_yp_an"], res["_yp_fem"]
    mm = np.isfinite(yp_an)
    ax[1].plot(sp[mm], yp_an[mm], color=BLUE, lw=3.0, label="pole FACE read off")
    ax[1].plot(sp, yp_fem, color=AMBER, lw=1.4, ls="--", label="FEM equipotential")
    ax[1].plot([s0 - 3, s0 + 2], [0, 0], color="#1a2230", lw=1.6)
    ax[1].axhline(H, color="#9aa6b6", ls=":", lw=1.0)
    ax[1].axvline(s0, color=RED, ls=":", lw=1.4, label="field edge s0")
    ax[1].text(s0 - 2.7, 0.35, "gap", color="#41506a", fontsize=10)
    ax[1].set_xlabel("s (beam direction)"); ax[1].set_ylabel("y (gap)")
    ax[1].set_xlim(s0 - 3, s0 + 1.5); ax[1].set_ylim(0, 1.3)
    ax[1].legend(fontsize=8.8, loc="upper left")
    ax[1].set_title(f"Pole FACE = equipotential; FEM match {res['pole_match_err']:.0e}")
    ax[2].plot(res["_ssv"], res["_By_target"], color=BLUE, lw=3.0, label="demanded g(s)")
    ax[2].plot(res["_ssv"], res["_By_fem"], color=RED, lw=1.3, ls="--", label="FEM B_y(s,0)")
    ax[2].axvline(s0, color="#1a2230", ls=":", lw=1.1)
    ax[2].set_xlabel("s (beam direction)"); ax[2].set_ylabel(r"$B_y(s,0)/B_0$")
    ax[2].set_xlim(s0 - 3, s0 + 3); ax[2].set_ylim(-0.05, 1.15)
    ax[2].legend(fontsize=9.5, loc="upper right")
    ax[2].set_title(f"Longitudinal fringe (FEM reproduces g to {res['By_rel_err']*100:.2f}%)")
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    import os
    return os.path.abspath(path)


# ------------------------------------------------------------------ PART B: end-shaping design
def ellipse_profile(a, b, n=60):
    """Forward quarter-ellipse chamfer: flat pole face for s<=0, then rise (0,h)->(a,h+b);
    horizontal tangent at the flat-face junction (no cusp), vertical tangent at the top."""
    t = np.linspace(0, np.pi / 2, n)
    return list(zip(a * np.sin(t), H + b * (1 - np.cos(t))))


def solve_endpeak(profile, Ls=8.0, Lr=8.0, Ly=6.0, maxh=0.05, order=3, want_field=False):
    """Forward FEM (ideal-iron pole with the given chamfer): return B0 and the peak pole-face
    |B|/B0 (a LOCAL, box/mesh-convergent quantity)."""
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry
    from ngsolve import y as Yc, grad, dx
    P = profile; s_top = P[-1][0]
    box = WorkPlane().MoveTo(-Ls, 0.0).Rectangle(Ls + Lr, Ly).Face()
    wp = WorkPlane().MoveTo(-Ls, H).LineTo(float(P[0][0]), float(P[0][1]))
    for (px, py) in P[1:]:
        wp.LineTo(float(px), float(py))
    wp.LineTo(float(s_top), Ly).LineTo(-Ls, Ly).Close()
    air = box - wp.Face()
    for e in air.edges:
        cx, cy = e.center[0], e.center[1]
        if abs(cy) < 1e-6:
            e.name = "mid"
        elif abs(cx + Ls) < 1e-6 and cy < H - 1e-9:
            e.name = "inlet"
        elif abs(cx - Lr) < 1e-6 or (abs(cy - Ly) < 1e-6 and cx > s_top + 1e-6):
            e.name = "far"
        else:
            e.name = "pole"
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(air, dim=2).GenerateMesh(maxh=maxh))
        fes = ng.H1(mesh, order=order, dirichlet="mid|inlet|pole|far")
        u, v = fes.TnT()
        A = ng.BilinearForm(grad(u) * grad(v) * dx); A.Assemble()
        gD = mesh.BoundaryCF({"mid": 0.0, "inlet": -Yc, "pole": -H, "far": 0.0}, default=0.0)
        gfu = ng.GridFunction(fes)
        gfu.Set(gD, definedon=mesh.Boundaries("mid|inlet|pole|far"))
        r = gfu.vec.CreateVector(); r.data = -A.mat * gfu.vec
        gfu.vec.data += A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
        gB = grad(gfu)
        B0 = -gB(mesh(-Ls + 1.0, 0.5 * H))[1]
        ca = (s_top + 2.0, 0.5 * H)
        peak = 0.0
        for (px, py) in ellipse_profile_pts_of(P):
            dxx, dyy = ca[0] - px, ca[1] - py
            nn = np.hypot(dxx, dyy)
            try:
                vv = gB(mesh(px + 0.012 * dxx / nn, py + 0.012 * dyy / nn))
                peak = max(peak, float(np.hypot(vv[0], vv[1])) / B0)
            except Exception:
                pass
    if want_field:
        return float(B0), float(peak), gfu, mesh, (Ls, Lr, Ly, s_top)
    return float(B0), float(peak)


def ellipse_profile_pts_of(P):
    return [(float(px), float(py)) for (px, py) in P]


DESIGN_B = 1.0                                        # optimal rise ~ gap h
DESIGN_AVALS = [0.4, 0.8, 1.4, 2.2, 3.2, 4.4]
DESIGN_OPT = (3.2, 1.0)


def design_curve(avals=None, b=DESIGN_B, maxh=0.06):
    """Peak pole-face |B|/B0 vs chamfer length a (rise b) -> monotone toward the Rogowski B0."""
    avals = DESIGN_AVALS if avals is None else avals
    peaks = []
    for a in avals:
        _, pk = solve_endpeak(ellipse_profile(a, b, 90), maxh=maxh)
        peaks.append(float(pk))
    return {"avals": list(map(float, avals)), "b": float(b), "peaks": peaks,
            "best_peak": float(min(peaks)), "best_a": float(avals[int(np.argmin(peaks))]),
            "rogowski_limit": 1.0}


def figure_design(curve, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, RED, GREEN, AMBER = "#1f6feb", "#d1495b", "#2e8b57", "#e0851e"
    shapes = [("near-square R=0.1", ellipse_profile(0.1, 0.1, 160), RED),
              ("round R=0.5", ellipse_profile(0.5, 0.5, 160), AMBER),
              (f"optimized a={DESIGN_OPT[0]} b={DESIGN_OPT[1]}", ellipse_profile(*DESIGN_OPT, 200), BLUE)]
    surf = []
    optfield = None
    for name, prof, c in shapes:
        B0, pk, gfu, mesh, geom = solve_endpeak(prof, maxh=0.045, want_field=True)
        gB = __import__("ngsolve").grad(gfu)
        ca = (prof[-1][0] + 2.0, 0.5 * H)
        vals = []
        for (px, py) in prof:
            dxx, dyy = ca[0] - px, ca[1] - py; nn = np.hypot(dxx, dyy)
            try:
                vv = gB(mesh(px + 0.012 * dxx / nn, py + 0.012 * dyy / nn))
                vals.append(float(np.hypot(vv[0], vv[1])) / B0)
            except Exception:
                vals.append(np.nan)
        surf.append((name, np.array(vals), c))
        if name.startswith("optimized"):
            optfield = (gfu, mesh, geom, B0)

    fig, ax = plt.subplots(1, 3, figsize=(16.6, 4.7), dpi=140)
    gfu, mesh, (Ls, Lr, Ly, s_top), B0 = optfield
    gB = __import__("ngsolve").grad(gfu)
    a_opt, b_opt = DESIGN_OPT
    sg = np.linspace(-1.0, 5.0, 260); yg = np.linspace(0.001, 3.2, 170)
    S, Y = np.meshgrid(sg, yg)
    BM = np.full(S.shape, np.nan)

    def ytop_of(s):
        if s < 0:
            return H
        if s <= a_opt:
            return H + b_opt * (1 - np.cos(np.arcsin(np.clip(s / a_opt, 0, 1))))
        return Ly
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            s, y = S[i, j], Y[i, j]
            if y < ytop_of(s) - 1e-3:
                try:
                    vv = gB(mesh(float(s), float(y))); BM[i, j] = float(np.hypot(vv[0], vv[1])) / B0
                except Exception:
                    pass
    pc = ax[0].pcolormesh(S, Y, np.clip(BM, 0, 1.6), cmap="magma", shading="auto")
    pp = np.array(ellipse_profile(a_opt, b_opt, 300))
    ax[0].plot(np.r_[-1, 0, pp[:, 0]], np.r_[H, H, pp[:, 1]], color="cyan", lw=2.4)
    ax[0].plot([-1, 5], [0, 0], color="w", lw=1.4)
    ax[0].text(-0.8, 0.4, "gap B0", color="w", fontsize=9)
    ax[0].text(1.4, 2.3, "iron pole", color="#cdd5df", fontsize=10)
    ax[0].set_xlabel("s"); ax[0].set_ylabel("y"); ax[0].set_xlim(-1, 5); ax[0].set_ylim(0, 3.2)
    ax[0].set_title("Optimized end: |B|/B0 (no hot spot)")
    fig.colorbar(pc, ax=ax[0], fraction=0.046)
    for name, vals, c in surf:
        xs = np.linspace(0, 1, len(vals))
        ax[1].plot(xs, vals, color=c, lw=2.2, label=f"{name} (peak {np.nanmax(vals):.2f})")
    ax[1].axhline(1.0, color="#1a2230", ls=":", lw=1.2, label="B0 (no enhancement)")
    ax[1].set_xlabel("normalized position along pole end (flat -> tip)")
    ax[1].set_ylabel(r"$|B|/B_0$ on pole surface")
    ax[1].set_ylim(0.6, 2.6); ax[1].legend(fontsize=8.4, loc="upper left")
    ax[1].set_title("Pole-surface field: square spikes, optimized ~ B0")
    ax[2].plot(curve["avals"], curve["peaks"], "o-", color=GREEN, lw=2.2, ms=6)
    ax[2].axhline(1.0, color=RED, ls="--", lw=1.4, label="Rogowski limit B0")
    ax[2].annotate("square end -> inf", (curve["avals"][0], curve["peaks"][0]),
                   (curve["avals"][0] + 0.3, 2.0), fontsize=9,
                   arrowprops=dict(arrowstyle="->", lw=1.0))
    ax[2].set_xlabel("chamfer length a  [gap h]"); ax[2].set_ylabel("peak |B|/B0")
    ax[2].set_ylim(0.95, 2.2); ax[2].legend(fontsize=9, loc="upper right")
    ax[2].set_title("Design curve: longer/gentler end -> B0")
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    import os
    return os.path.abspath(path)


if __name__ == "__main__":
    for r in feasibility_table([1.4, 1.0, 0.8, d_star(), 0.5, 0.35]):
        print(f"  d={r['d']:.3f}: y_sing={r['y_sing']:.3f} -> "
              f"{'FEASIBLE' if r['feasible'] else 'INFEASIBLE'}")
    res = fem_verify_sy()
    print(f"FEM fringe: By rel err {res['By_rel_err']*100:.3f}%, phi L2 {res['phi_L2_err']:.1e}, "
          f"pole match {res['pole_match_err']:.1e}")
    cur = design_curve()
    print(f"design curve peaks {['%.3f' % p for p in cur['peaks']]}; "
          f"best {cur['best_peak']:.3f} at a={cur['best_a']} (Rogowski 1.0)")
