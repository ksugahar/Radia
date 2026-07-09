r"""Vertical EDGE FOCUSING of a tilted dipole end -- measured by PARTICLE TRACKING.

Imported by edge_focusing_tracking.ipynb and by the golden
tests/feec/test_edge_focusing_tracking.py.  Companion of hodograph_bending_sy.py:
that unit shapes the s-y (longitudinal) end so the pole face stays at B0; THIS unit
answers the orthogonal question -- what a HORIZONTAL edge tilt (angle beta, a rotation
of the pole face about the vertical axis, in the x-s bend plane) does to the VERTICAL
optics.

Physics (hard-edge accelerator optics).  A dipole edge tilted by beta acts as a thin
vertical lens whose strength is
    | 1 / f_z | = tan(beta) / rho ,    rho = p / (q B0)   (the bend radius).
The MAGNITUDE tan(beta)/rho is the convention-independent, delicate quantity; the SIGN
is orientation-dependent (entrance vs exit edge, and which way the wedge opens).  For a
rectangular magnet BOTH edges DEFOCUS vertically; the entrance-edge orientation modeled
here (a genuinely curl-free fringe, see the sign note) FOCUSES, giving 1/f_z = +tan/rho.

Why tracking, not a field-EFB slope.  The vertical edge focusing is a SECOND-ORDER,
off-mid-plane property of the fringe; the effective-field-boundary (EFB) slope of the
mid-plane |B| CANNOT recover it (it is wrong-sign / blows up on a compact dipole; see
memory/edge_focusing_efb_slope_negative).  The correct measurement is the linearized
vertical Hill integral along the reference orbit:
    1 / f_z = (q/p) INT ( u_y dB_x/dz - u_x dB_y/dz )|_{z=0} ds ,
with u the mid-plane orbit tangent.  This module (a) builds a genuinely MAXWELLIAN
tilted hard-edge fringe (curl-free AND div-free -- see the sign note below), (b) tracks
the reference orbit + evaluates that integral, and (c) shows it reproduces tan(beta)/rho
to <1%, converging to the hard-edge law as the fringe width w -> 0, and collapsing across
rho.  That validates the METHOD on a field whose answer is known in closed form.

Sign note (why a genuinely curl-free fringe matters).  Writing the tilted mid-plane
profile B_z(s,0) = B0 g(s) with edge-normal s = (y - y_edge) cos b + x sin b, the ONLY
vacuum (curl-free + div-free) linear-in-z continuation is
    B_s = + B0 z g'(s) ,   B_z = B0 ( g(s) - 1/2 z^2 g''(s) ) ,
so dB_x/dz = +B0 g'(s) sin b, dB_y/dz = +B0 g'(s) cos b.  A div-free-BUT-NOT-curl-free
choice (B_s = -B0 z g'(s)) is a different field with a spurious current sheet at the edge
and FLIPS the focusing sign -- a real trap.  We use the curl-free field; the magnitude
tan(beta)/rho is convention-independent, which is what the tracker measures.
"""
import numpy as np

B0 = 1.0            # mid-plane flat-top field (T, normalized)
W = 0.02            # fringe (edge) width used for the reference demonstration
Y_EDGE = 0.0        # edge crossing on the reference orbit
RHO_REF = 1.0       # reference bend radius for the beta sweep
BETAS_DEG = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]


# ------------------------------------------------------------------ Maxwellian tilted edge
def edge_field(beta, b0=B0, w=W, y_edge=Y_EDGE):
    """Return a callable r=[x,y,z] -> [Bx,By,Bz]: a genuinely Maxwellian (curl-free +
    div-free) hard-edge fringe, edge tilted by beta about the vertical axis."""
    sb, cb = np.sin(beta), np.cos(beta)

    def field(r):
        x, y, z = r
        s = (y - y_edge) * cb + x * sb                 # edge-normal coordinate
        t = np.tanh(s / w)
        gp = 0.5 * (1.0 - t * t) / w                   # g'(s), g = 1/2 (1 + tanh)
        gpp = -t * (1.0 - t * t) / (w * w)             # g''(s)
        Bs = b0 * z * gp                               # curl-free continuation (sign!)
        Bz = b0 * (0.5 * (1.0 + t) - 0.5 * z * z * gpp)
        return np.array([Bs * sb, Bs * cb, Bz])

    return field


# ------------------------------------------------------------------ tracker + Hill integral
def edge_focus_integral(field, rho, b0=B0, y0=-0.35, y1=0.35, ds=2.5e-4, dz=2.0e-3):
    """Vertical edge focusing 1/f_z along the mid-plane reference orbit (RK4), via the
    linearized Hill integral 1/f_z = (q/p) INT (u_y dB_x/dz - u_x dB_y/dz)|_{z=0} ds.
    dB_x/dz, dB_y/dz by central z-difference (exact for the linear-in-z Maxwellian field).
    Also returns the orbit exit state for diagnostics."""
    qop = 1.0 / (b0 * rho)                             # q/p = 1/(B0 rho)

    def bz0(r):
        return np.array([0.0, 0.0, field([r[0], r[1], 0.0])[2]])  # mid-plane bend field

    r = np.array([0.0, y0, 0.0])
    u = np.array([0.0, 1.0, 0.0])
    K = 0.0
    n = 0
    while r[1] < y1:
        f = lambda rr, uu: qop * np.cross(uu, bz0(rr))
        k1 = f(r, u)
        k2 = f(r + 0.5 * ds * u, u + 0.5 * ds * k1)
        k3 = f(r + 0.5 * ds * (u + 0.5 * ds * k1), u + 0.5 * ds * k2)
        k4 = f(r + ds * (u + 0.5 * ds * k2), u + ds * k3)
        r = r + ds * (u + ds / 6.0 * (k1 + 2 * k2 + 2 * k3))
        u = u + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        u = u / np.linalg.norm(u)
        dBx = (field([r[0], r[1], dz])[0] - field([r[0], r[1], -dz])[0]) / (2 * dz)
        dBy = (field([r[0], r[1], dz])[1] - field([r[0], r[1], -dz])[1]) / (2 * dz)
        K += qop * (u[1] * dBx - u[0] * dBy) * ds
        n += 1
        if n > 20000:
            raise RuntimeError("edge_focus_integral: step cap")
    return {"inv_fz": float(K), "x_exit": float(r[0]), "ux_exit": float(u[0])}


def hard_edge_law(beta_deg, rho):
    """Thin-lens hard-edge vertical edge focusing for THIS (curl-free entrance-edge)
    orientation: 1/f_z = +(tan beta)/rho.  The invariant is the magnitude tan(beta)/rho;
    the sign is orientation-dependent (a rectangular-magnet edge would give -tan/rho)."""
    return np.tan(np.radians(beta_deg)) / rho


# ------------------------------------------------------------------ studies
def sweep_beta(betas_deg=None, rho=RHO_REF, w=W):
    betas_deg = BETAS_DEG if betas_deg is None else betas_deg
    rows = []
    for bd in betas_deg:
        r = edge_focus_integral(edge_field(np.radians(bd), w=w), rho)
        rows.append({"beta_deg": float(bd), "inv_fz": r["inv_fz"],
                     "hard_edge": float(hard_edge_law(bd, rho)),
                     "x_exit": r["x_exit"]})
    return rows


def w_convergence(ws=None, beta_deg=20.0, rho=RHO_REF):
    """Fit slope c(w) of 1/f_z vs -tan(beta)/rho over a small beta pencil; c -> 1 as w -> 0
    demonstrates convergence to the hard-edge law."""
    ws = [0.08, 0.04, 0.02, 0.01, 0.005] if ws is None else ws
    pencil = [5.0, 10.0, 15.0, 20.0]
    out = []
    for w in ws:
        num = np.array([edge_focus_integral(edge_field(np.radians(bd), w=w), rho)["inv_fz"]
                        for bd in pencil])
        den = np.array([hard_edge_law(bd, rho) for bd in pencil])
        c = float(np.dot(num, den) / np.dot(den, den))   # least-squares slope vs the law
        out.append({"w": float(w), "slope_vs_law": c})
    return out


def rho_collapse(rhos=None, betas_deg=None, w=W):
    """1/f_z * rho vs tan(beta) collapses onto -tan(beta) across rho (rho-independence)."""
    rhos = [0.7, 1.0, 1.6] if rhos is None else rhos
    betas_deg = [0.0, 10.0, 20.0, 30.0] if betas_deg is None else betas_deg
    out = []
    for rho in rhos:
        rows = sweep_beta(betas_deg, rho=rho, w=w)
        out.append({"rho": float(rho),
                    "tan_beta": [float(np.tan(np.radians(r["beta_deg"]))) for r in rows],
                    "inv_fz_rho": [float(r["inv_fz"] * rho) for r in rows]})
    return out


def summarize(sweep, wconv):
    """Max relative error of the tracked 1/f_z vs the hard-edge law over the sweep
    (beta=0 excluded: the law is 0 there), and the finest-w slope."""
    rel = [abs(r["inv_fz"] - r["hard_edge"]) / abs(r["hard_edge"])
           for r in sweep if abs(r["hard_edge"]) > 1e-9]
    return {"max_rel_err_vs_law": float(max(rel)),
            "beta0_baseline": float(next(r["inv_fz"] for r in sweep if r["beta_deg"] == 0.0)),
            "finest_w": float(wconv[-1]["w"]),
            "finest_w_slope": float(wconv[-1]["slope_vs_law"])}


# ------------------------------------------------------------------ figure
def figure_edge_focus(sweep, wconv, rcol, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, RED, GREEN, AMBER = "#1f6feb", "#d1495b", "#2e8b57", "#e0851e"
    fig, ax = plt.subplots(1, 3, figsize=(16.6, 4.7), dpi=140)

    tb = np.array([np.tan(np.radians(r["beta_deg"])) for r in sweep])
    meas = np.array([r["inv_fz"] for r in sweep])
    law = np.array([r["hard_edge"] for r in sweep])
    ax[0].plot(tb, law, color=BLUE, lw=3.0, label=r"hard edge $\tan\beta/\rho$")
    ax[0].plot(tb, meas, "o", color=RED, ms=7, label="tracked (Hill integral)")
    ax[0].set_xlabel(r"$\tan\beta$"); ax[0].set_ylabel(r"$1/f_z$  [1/m]")
    ax[0].legend(fontsize=9.5, loc="upper right")
    ax[0].set_title(r"Edge focusing vs tilt ($\rho=%.1f$ m)" % RHO_REF)

    ws = np.array([d["w"] for d in wconv]); cs = np.array([d["slope_vs_law"] for d in wconv])
    ax[1].semilogx(ws, cs, "o-", color=GREEN, lw=2.2, ms=6)
    ax[1].axhline(1.0, color=RED, ls="--", lw=1.4, label="hard-edge limit c=1")
    ax[1].set_xlabel("fringe width w  [m]"); ax[1].set_ylabel(r"slope $c(w)$ of $1/f_z$ vs law")
    ax[1].legend(fontsize=9.5, loc="lower left")
    ax[1].set_title(r"Converges to the law as $w\to0$")
    ax[1].invert_xaxis()

    for d, c in zip(rcol, (BLUE, AMBER, GREEN)):
        ax[2].plot(d["tan_beta"], d["inv_fz_rho"], "o-", color=c, lw=2.0, ms=5,
                   label=fr"$\rho={d['rho']:.1f}$ m")
    xx = np.linspace(0, max(rcol[0]["tan_beta"]), 50)
    ax[2].plot(xx, xx, color="#1a2230", ls=":", lw=1.6, label=r"$\tan\beta$")
    ax[2].set_xlabel(r"$\tan\beta$"); ax[2].set_ylabel(r"$\rho \cdot 1/f_z$")
    ax[2].legend(fontsize=9.0, loc="upper right")
    ax[2].set_title(r"$\rho\,/f_z$ collapses (edge is $\rho$-scaled)")

    fig.tight_layout(); fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    import os
    return os.path.abspath(path)


if __name__ == "__main__":
    sw = sweep_beta()
    print("beta[deg]   tracked 1/f_z    hard-edge -tan/rho    rel.err")
    for r in sw:
        rel = (abs(r["inv_fz"] - r["hard_edge"]) / abs(r["hard_edge"])
               if abs(r["hard_edge"]) > 1e-9 else float("nan"))
        print(f"  {r['beta_deg']:5.1f}   {r['inv_fz']:+.5f}      {r['hard_edge']:+.5f}"
              f"       {rel*100:6.2f}%")
    wc = w_convergence()
    print("\nw-convergence (slope of 1/f_z vs the hard-edge law):")
    for d in wc:
        print(f"  w={d['w']:.3f}: slope={d['slope_vs_law']:.4f}")
    s = summarize(sw, wc)
    print(f"\nmax rel err vs law {s['max_rel_err_vs_law']*100:.2f}%, "
          f"beta=0 baseline {s['beta0_baseline']:.2e}, "
          f"finest-w slope {s['finest_w_slope']:.4f}")
