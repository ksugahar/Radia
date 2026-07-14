r"""Vertical EDGE FOCUSING of a tilted dipole end -- measured by PARTICLE TRACKING.

Imported by edge_focusing_tracking.ipynb and by the golden
validation_test/feec/test_edge_focusing_tracking.py.  Companion of hodograph_bending_sy.py:
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


def psi_enge(beta_rad, K1g, rho):
    """Classical Enge fringe correction psi = (K1g/rho)(1+sin^2 beta)/cos beta (Brown,
    SLAC-75 first order), with K1g = INT g(1-g) ds along the edge NORMAL (units m; this
    is K1*gap with the gap convention divided out).  For the tanh fringe K1g = w/2."""
    return (K1g / rho) * (1.0 + np.sin(beta_rad) ** 2) / np.cos(beta_rad)


def scoff_law(beta_deg, rho, K1g):
    """SCOFF + Enge fringe-corrected vertical edge focusing for this orientation:
    1/f_z = tan(beta - psi)/rho.  At beta=0 this is -K1g/rho^2 -- exactly the tracked
    finite-fringe baseline (-w/(2 rho^2) for the tanh fringe)."""
    b = np.radians(beta_deg)
    return np.tan(b - psi_enge(b, K1g, rho)) / rho


# ------------------------------------------------------------------ studies
def sweep_beta(betas_deg=None, rho=RHO_REF, w=W):
    betas_deg = BETAS_DEG if betas_deg is None else betas_deg
    rows = []
    for bd in betas_deg:
        r = edge_focus_integral(edge_field(np.radians(bd), w=w), rho)
        rows.append({"beta_deg": float(bd), "inv_fz": r["inv_fz"],
                     "hard_edge": float(hard_edge_law(bd, rho)),
                     "enge": float(scoff_law(bd, rho, 0.5 * w)),
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
    """Max relative error of the tracked 1/f_z vs the hard-edge law AND vs the Enge
    fringe-corrected law over the sweep (beta=0 excluded for the hard-edge ratio: the
    law is 0 there), and the finest-w slope."""
    rel = [abs(r["inv_fz"] - r["hard_edge"]) / abs(r["hard_edge"])
           for r in sweep if abs(r["hard_edge"]) > 1e-9]
    rel_e = [abs(r["inv_fz"] - r["enge"]) / abs(r["enge"])
             for r in sweep if abs(r["enge"]) > 1e-9]
    return {"max_rel_err_vs_law": float(max(rel)),
            "max_rel_err_vs_enge": float(max(rel_e)),
            "beta0_baseline": float(next(r["inv_fz"] for r in sweep if r["beta_deg"] == 0.0)),
            "finest_w": float(wconv[-1]["w"]),
            "finest_w_slope": float(wconv[-1]["slope_vs_law"])}


# ================================================================== PART B: FEM (SCOFF/Enge)
# 3D FEM of a PARALLELOGRAM dipole (both edges tilted by beta), measured with the SAME
# tracker: symmetric (closed-orbit style) traversal + window decomposition + an x-uniform
# model built from the FEM's own mid-line fringe profile.  The FEM entrance-edge kick is
# compared window-by-window against the model and the SCOFF/Enge closed form
# tan(beta_eff - psi)/rho with the FEM-measured K1g.
#
# TWO independent field engines feed the same measurement chain (swap via the
# fem_scoff_study(solve_midplane=...) argument):
#   * fem_solve_midplane  -- reduced-Omega NGSolve FEM (H = Hs - grad(Omega), air box,
#     order-2 H1);
#   * hdiv_solve_midplane -- FEEC HDiv-VIM (radia.vim): iron-only tet mesh, NO air
#     discretization, exact open boundary, batch rad.Fld analytic map.
# Cross-check (2026-07-13): dK_in agrees to 0.8% at matched edge-mesh density
# (absolute-dK scatter across engines/meshes ~+-3%), and the deficit vs the
# x-uniform model persists in every configuration (dK_in/model = 0.92-0.95) --
# which REASSIGNS it as REAL 3D physics: the local iso-field tilt near x=0 is only
# ~0.95-0.96 of the geometric tan(beta), measured directly on both engines' maps.
# See edge_focusing_fem_results.json `hdiv_vim_cross_check`.
#
# The parallelogram testbed (the standard spectrometer-dipole configuration) is chosen
# for two reasons learned the hard way:
#   1. the COIL must follow the pole contour -- a straight coil front across a tilted
#      iron edge misaligns by up to x_side*tan(beta) ~ 29 mm over a ~60 mm fringe and
#      the iso-field tilt lands far below the geometric edge angle (measured dK deficit
#      0.55x); a rigid whole-coil rotation instead breaks the MMF linkage topology
#      (conductor escaping past the return leg collapsed B0 by 3x).  The sheared
#      rounded-parallelogram loop follows both edges at a fixed normal offset.
#   2. parallelogram poles + centered legs + the C2-symmetric coil pair make the WHOLE
#      magnet exactly C2-symmetric (180-deg rotation about z): Bz(x,y,0)=Bz(-x,-y,0),
#      so the C2-odd part of the sampled map is pure solver error and is removed
#      exactly by map symmetrization at BOTH edge angles.
F_GAP = 0.040            # pole gap (z)
F_POLE_W = 0.12          # pole half-width (x) -- wide vs the gap
F_T_LEG = 0.03           # return legs: flux clamp -> short fringe, iron-circuit-uniform field
F_LEG_HY = 0.05          # legs y in (-F_LEG_HY, +F_LEG_HY): centered -> C2 symmetry
F_Z_OUT = 0.100
F_L_BEAM = 0.20          # long magnet -> separated fringes, genuine flat top
F_AIR = 0.45
F_COIL_W, F_COIL_H = 0.015, 0.025
F_COIL_XSIDE = 0.08      # coil side straights at x = +/-F_COIL_XSIDE (thread the iron circuit)
F_COIL_RC = 0.02         # corner-arc radius of the rounded-parallelogram loop
F_COIL_NOFF = 0.02       # coil bar offset from the iron edge along the edge NORMAL
F_NI = 10000.0
F_MU0 = 4.0 * np.pi * 1e-7


def _fem_coil_path(z_sign, beta_rad):
    """ONE rounded-parallelogram coil loop at z = z_sign*(F_GAP/2 + F_COIL_H/2).

    Front/back bars parallel to the beta-tilted iron edges at F_COIL_NOFF normal
    offset, side straights at x = +/-F_COIL_XSIDE, F_COIL_RC corner arcs.  The
    path is C2-symmetric about the origin and closes to machine precision
    (corner tangent-trim arithmetic; verified closure gap ~3e-17 m).
    """
    from radia.coil_builder import CoilBuilder
    z_coil = z_sign * (F_GAP / 2 + F_COIL_H / 2)  # winding at the pole-face level (a
    # fully-buried winding shunts its MMF in local iron loops)
    tb, cbeta = np.tan(beta_rad), np.cos(beta_rad)
    bdeg = np.degrees(beta_rad)
    y_off = F_L_BEAM / 2 + F_COIL_NOFF / cbeta
    t_a = F_COIL_RC * np.tan(np.radians((90 - bdeg) / 2))
    t_o = F_COIL_RC * np.tan(np.radians((90 + bdeg) / 2))
    L_side = 2 * y_off - t_o - t_a
    L_bar = 2 * F_COIL_XSIDE / cbeta - t_a - t_o
    return (CoilBuilder(current=F_NI)
            .set_start([F_COIL_XSIDE, -y_off - F_COIL_XSIDE * tb + t_o, z_coil])
            .set_cross_section(width=F_COIL_W, height=F_COIL_H)
            .add_straight(L_side).add_arc(radius=F_COIL_RC, arc_angle=90 - bdeg)
            .add_straight(L_bar).add_arc(radius=F_COIL_RC, arc_angle=90 + bdeg)
            .add_straight(L_side).add_arc(radius=F_COIL_RC, arc_angle=90 - bdeg)
            .add_straight(L_bar).add_arc(radius=F_COIL_RC, arc_angle=90 + bdeg))


def fem_build_coil(beta_rad=0.0):
    """Coil pair for the parallelogram magnet.

    BOTH loops are built EXPLICITLY with the same CCW traversal and +NI current
    (a z-plane pair adds Bz on the mid-plane).  Explicit construction keeps this
    validation geometry independent of the CoilBuilder.mirror() implementation;
    mirror() itself is covered by tests/test_coil_builder_mirror.py.  Do not wrap
    the container with rad.TrfOrnt: rad.RadiaField on a
        TrfOrnt-wrapped container crashes with 0xC0000374 heap corruption
    (plain containers assemble fine).
    """
    import radia as rad
    rad.UtiDelAll()
    up = _fem_coil_path(+1, beta_rad)
    lo = _fem_coil_path(-1, beta_rad)
    return rad.ObjCnt(up.to_radia() + lo.to_radia())


def _fem_footprint(z0, z1, beta):
    from netgen.occ import WorkPlane, Axes, Pnt, Z, X
    hL = F_L_BEAM / 2
    tb = np.tan(beta)
    wp = WorkPlane(Axes(Pnt(0, 0, z0), n=Z, h=X))
    wp.MoveTo(-F_POLE_W, -hL + F_POLE_W * tb)
    wp.LineTo(F_POLE_W, -hL - F_POLE_W * tb)   # tilted ENTRANCE (-y) edge, pivot at x=0
    wp.LineTo(F_POLE_W, hL - F_POLE_W * tb)    # PARALLELOGRAM: exit edge tilted the same
    wp.LineTo(-F_POLE_W, hL + F_POLE_W * tb)   # way -> the whole magnet is C2-symmetric
    wp.Close()
    return wp.Face().Extrude(z1 - z0)


def fem_build_mesh(beta_deg, maxh_air=0.035, maxh_iron=0.014):
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    beta = np.radians(beta_deg)
    top = _fem_footprint(F_GAP / 2, F_Z_OUT, beta)
    bot = _fem_footprint(-F_Z_OUT, -F_GAP / 2, beta)
    leg_l = Box(Pnt(-F_POLE_W, -F_LEG_HY, -F_GAP / 2), Pnt(-F_POLE_W + F_T_LEG, F_LEG_HY, F_GAP / 2))
    leg_r = Box(Pnt(F_POLE_W - F_T_LEG, -F_LEG_HY, -F_GAP / 2), Pnt(F_POLE_W, F_LEG_HY, F_GAP / 2))
    iron = top + bot + leg_l + leg_r            # centered legs -> C2; clear of both edges
    iron.mat("iron"); iron.maxh = maxh_iron
    # graded refinement at the gap-facing pole-face EDGE LINES (entrance/exit, top/bot):
    # the mu->inf corner singularity there drives a smooth mesh-realization-dependent
    # odd-in-x error field on the mid-plane
    for e in iron.edges:
        c = e.center
        if abs(abs(c.z) - F_GAP / 2) < 1e-4 and abs(c.x) < 0.11 and abs(c.y) > 0.05:
            e.maxh = 0.004
    # fine air box around the beam: the tracked mid-plane gradients live on AIR elements;
    # coarse air (maxh ~ gap) leaves mesh-asymmetry dBz/dx(x=0) ~ 0.1-0.7 T/m (measured)
    fine = Box(Pnt(-0.05, -0.27, -0.03), Pnt(0.05, 0.27, 0.03)) - iron
    fine.mat("air"); fine.maxh = 0.008
    air = Box(Pnt(-F_AIR, -F_AIR, -F_AIR), Pnt(F_AIR, F_AIR, F_AIR)) - iron - fine
    air.mat("air")
    for f in air.faces:
        c = f.center
        if max(abs(c.x), abs(c.y), abs(c.z)) > 0.9 * F_AIR:
            f.name = "outer"
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Glue([air, fine, iron])).GenerateMesh(maxh=maxh_air))
    return mesh


def fem_solve_midplane(beta_deg, mu_r=1000.0, order=2, nx=81, ny=401,
                       xmax=0.05, ymax=0.26):
    """Reduced-Omega solve (H = Hs - grad(Omega), Biot-Savart coil source), then sample the
    MID-PLANE Bz(x,y) map -- by paraxial theory the vertical dynamics depends only on it."""
    import ngsolve as ng
    from ngsolve import H1, BilinearForm, LinearForm, GridFunction, grad, dx, TaskManager
    import radia as rad
    mesh = fem_build_mesh(beta_deg)
    coils = fem_build_coil(float(np.radians(beta_deg)))
    Hs = rad.RadiaField(coils, "h")
    mu = mesh.MaterialCF({"iron": mu_r * F_MU0}, default=F_MU0)
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    f = LinearForm(fes); f += mu * Hs * grad(v) * dx; f.Assemble()   # serial (RadiaField)
    with TaskManager():
        a = BilinearForm(fes); a += mu * grad(u) * grad(v) * dx; a.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    B = mu * (Hs - grad(gfu))
    xs = np.linspace(-xmax, xmax, nx)
    ys = np.linspace(-ymax, ymax, ny)
    Bz = np.zeros((nx, ny))
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            try:
                Bz[i, j] = float(B(mesh(float(x), float(y), 0.0))[2])
            except Exception:
                pass
    ne = int(mesh.ne)
    rad.UtiDelAll()
    return xs, ys, Bz, ne


def fem_profile(xs, ys, Bz, beta_deg):
    """From the x=0 mid-line: flat-top B0, flat window, K1g of each fringe (edge-normal),
    hard-edge-equivalent EFB positions."""
    i0 = int(np.argmin(np.abs(xs)))
    P = Bz[i0, :]
    B0 = float(np.mean(P[np.abs(ys) < 0.015]))
    g = P / B0
    flat = np.where(g > 0.995)[0]
    y_a, y_b = float(ys[flat[0]]), float(ys[flat[-1]])
    cb = np.cos(np.radians(beta_deg))
    dy = ys[1] - ys[0]
    gin = np.clip(np.where(ys <= 0.0, g, 1.0), 0.0, None)
    gout = np.clip(np.where(ys >= 0.0, g, 1.0), 0.0, None)
    K1g_in = float(np.sum(gin * (1 - gin)) * dy * cb)      # ds = dy cos(beta) on x=0
    K1g_out = float(np.sum(gout * (1 - gout)) * dy * cb)   # exit edge tilted too (C2 magnet)
    yE_in = y_a - float(np.sum(gin[ys <= y_a]) * dy)
    yE_out = y_b + float(np.sum(gout[ys >= y_b]) * dy)
    return {"B0": B0, "y_flat": (y_a, y_b), "K1g_in": K1g_in, "K1g_out": K1g_out,
            "yE_in": float(yE_in), "yE_out": float(yE_out), "g_line": g, "ys": ys}


def fem_model_Bz(prof, beta_deg, xs, ys):
    """X-UNIFORM model of the same magnet built from the MEASURED mid-line profile:
    Bz(x,y) = B0 Gin(s_in) Gout(s_out) with BOTH factors re-parameterized to the tilted
    edge-normal coordinates of the parallelogram (s_in = (y+hL)cb + x sb,
    s_out = (y-hL)cb + x sb).  Same tilts, same thick fringes, same two edges -- but no
    transverse structure, so FEM - model isolates the x-nonuniformity contamination."""
    from scipy.interpolate import interp1d
    g, yline = prof["g_line"], prof["ys"]
    beta = np.radians(beta_deg)
    cb, sb = np.cos(beta), np.sin(beta)
    ymid = 0.5 * (prof["y_flat"][0] + prof["y_flat"][1])
    gin_y = np.where(yline <= ymid, g, 1.0)
    gout_y = np.where(yline >= ymid, g, 1.0)
    hL = F_L_BEAM / 2
    Gin = interp1d((yline + hL) * cb, gin_y, bounds_error=False,
                   fill_value=(float(gin_y[0]), 1.0))
    Gout = interp1d((yline - hL) * cb, gout_y, bounds_error=False,
                    fill_value=(1.0, float(gout_y[-1])))
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return prof["B0"] * Gin((Y + hL) * cb + X * sb) * Gout((Y - hL) * cb + X * sb)


def track_midplane(xs, ys, Bz, B0, rho, prof, y0=-0.24, y1=0.24, ds=5e-4,
                   x0=0.0, a0=0.0):
    """Reference-orbit RK4 on the mid-plane Bz map + the linearized vertical system:
    G = kappa (u_y dBz/dx - u_x dBz/dy), 2x2 matrix M_z, and the kick integral K
    decomposed into entrance / flat / exit windows and into T1 (transverse-gradient)
    and T2 (edge-crossing) terms."""
    from scipy.interpolate import RegularGridInterpolator
    dBzdx = np.gradient(Bz, xs, axis=0)
    dBzdy = np.gradient(Bz, ys, axis=1)
    fB = RegularGridInterpolator((xs, ys), Bz, bounds_error=False, fill_value=0.0)
    fX = RegularGridInterpolator((xs, ys), dBzdx, bounds_error=False, fill_value=0.0)
    fY = RegularGridInterpolator((xs, ys), dBzdy, bounds_error=False, fill_value=0.0)
    kap = 1.0 / (B0 * rho)
    y_a, y_b = prof["y_flat"]
    r = np.array([x0, y0]); u = np.array([np.sin(a0), np.cos(a0)])
    M = np.eye(2)
    K = {"in": 0.0, "flat": 0.0, "out": 0.0}
    T1 = {"in": 0.0, "flat": 0.0, "out": 0.0}
    T2 = {"in": 0.0, "flat": 0.0, "out": 0.0}
    alpha_at_EFB = None
    x_mid = None
    x_absmax = 0.0
    n = 0
    while r[1] < y1:
        b = lambda rr: float(fB([rr])[0])
        f = lambda rr, uu: kap * np.array([uu[1] * b(rr), -uu[0] * b(rr)])
        k1 = f(r, u); k2 = f(r + 0.5 * ds * u, u + 0.5 * ds * k1)
        k3 = f(r + 0.5 * ds * (u + 0.5 * ds * k1), u + 0.5 * ds * k2)
        k4 = f(r + ds * (u + 0.5 * ds * k2), u + ds * k3)
        r = r + ds * (u + ds / 6 * (k1 + 2 * k2 + 2 * k3))
        u = u + ds / 6 * (k1 + 2 * k2 + 2 * k3 + k4); u = u / np.linalg.norm(u)
        t1 = kap * u[1] * float(fX([r])[0])
        t2 = -kap * u[0] * float(fY([r])[0])
        G = t1 + t2
        w = "in" if r[1] < y_a else ("out" if r[1] > y_b else "flat")
        K[w] += G * ds; T1[w] += t1 * ds; T2[w] += t2 * ds
        M = M + ds * np.array([[0.0, 1.0], [-G, 0.0]]) @ M
        if alpha_at_EFB is None and r[1] >= prof["yE_in"]:
            alpha_at_EFB = np.arctan2(u[0], u[1])
        if x_mid is None and r[1] >= 0.0:
            x_mid = float(r[0])
        x_absmax = max(x_absmax, abs(float(r[0])))
        n += 1
        if n > 20000:
            raise RuntimeError("track_midplane: step cap")
    if abs(r[0]) > xs[-1] * 0.98:
        raise RuntimeError(f"orbit left the sampled x-window: x={r[0]:.4f}")
    return {"K": K, "T1": T1, "T2": T2,
            "K_total": K["in"] + K["flat"] + K["out"], "M21": float(M[1, 0]),
            "alpha_EFB_deg": float(np.degrees(alpha_at_EFB or 0.0)),
            "x_mid": float(x_mid if x_mid is not None else r[0]),
            "x_absmax": float(x_absmax), "x_exit": float(r[0]),
            "theta_exit_deg": float(np.degrees(np.arctan2(u[0], u[1])))}


def track_symmetric(xs, ys, Bz, B0, rho, prof, ds=5e-4, ds_shoot=2e-3):
    """CLOSED-ORBIT style symmetric traversal: launch angle -theta/2 and offset so the
    orbit crosses the magnet centre near x=0 -> max|x| ~ sagitta, and the odd-in-x
    transverse gradient dBz/dx cancels between the two halves (suppresses T1)."""
    p1 = track_midplane(xs, ys, Bz, B0, rho, prof, ds=ds_shoot)
    a0 = -0.5 * np.radians(p1["theta_exit_deg"])
    p2 = track_midplane(xs, ys, Bz, B0, rho, prof, ds=ds_shoot, a0=a0)
    return track_midplane(xs, ys, Bz, B0, rho, prof, ds=ds, a0=a0, x0=-p2["x_mid"])


def fem_scoff_study(betas_deg=(0.0, 20.0), rho=5.0, solve_midplane=None):
    """The full SCOFF/Enge + closed-orbit FEM measurement: for each beta, solve the FEM,
    C2-symmetrize the map, measure the fringe (B0, K1g, EFB), track FEM and x-uniform
    model symmetrically, and compare the entrance-window kick against the closed form
    tan(beta_eff - psi)/rho.

    ``solve_midplane`` swaps the FIELD ENGINE only (default: the reduced-Omega
    ``fem_solve_midplane``; pass ``hdiv_solve_midplane`` for the HDiv-VIM twin).
    Everything downstream -- C2 symmetrization, profile, model, tracker, closed
    form -- is engine-independent, which is what makes the cross-check meaningful."""
    if solve_midplane is None:
        solve_midplane = fem_solve_midplane
    cases = []
    for bd in betas_deg:
        xs, ys, Bz, ne = solve_midplane(bd)
        # C2 map symmetrization: the parallelogram magnet satisfies
        # Bz(x,y,0) = Bz(-x,-y,0) EXACTLY, so the C2-odd part of the sampled map
        # is pure solver/mesh error (grids are symmetric linspaces).
        Bz = 0.5 * (Bz + Bz[::-1, ::-1])
        prof = fem_profile(xs, ys, Bz, bd)
        fem = track_symmetric(xs, ys, Bz, prof["B0"], rho, prof)
        mod = track_symmetric(xs, ys, fem_model_Bz(prof, bd, xs, ys),
                              prof["B0"], rho, prof)
        b_eff = np.radians(bd) - np.radians(fem["alpha_EFB_deg"])
        closed = float(np.tan(b_eff - psi_enge(b_eff, prof["K1g_in"], rho)) / rho)
        cases.append({"beta_deg": float(bd), "ne": ne, "rho": rho,
                      "B0": prof["B0"], "K1g_in": prof["K1g_in"],
                      "K1g_out": prof["K1g_out"], "y_flat": list(prof["y_flat"]),
                      "fem": fem, "model": mod, "closed_form_entrance": closed,
                      "beta_eff_deg": float(np.degrees(b_eff))})
    res = {"rho": rho, "cases": cases}
    if len(cases) >= 2 and cases[0]["beta_deg"] == 0.0:
        c0 = cases[0]
        for c in cases[1:]:
            c["dK_in_fem"] = c["fem"]["K"]["in"] - c0["fem"]["K"]["in"]
            c["dK_in_model"] = c["model"]["K"]["in"] - c0["model"]["K"]["in"]
            c["dK_in_closed"] = c["closed_form_entrance"] - c0["closed_form_entrance"]
    return res


# --------------------------------------------------- HDiv-VIM field engine (cross-check)
def hdiv_build_iron_mesh(beta_deg, maxh_iron=0.014, edge_maxh=0.004, face_maxh=None):
    """IRON-ONLY tet mesh for the HDiv-VIM engine: the same parallelogram poles +
    centered legs as ``fem_build_mesh``, but NO air box -- the HDiv-VIM needs no air
    discretization (the mid-plane field is evaluated by analytic surface-charge
    integrals of the solved per-element M, exact open boundary).

    ``face_maxh`` refines the GAP-FACING pole faces (|z| = F_GAP/2).  The mid-plane
    field is dominated by the surface sources on those faces at only F_GAP/2 = 20 mm
    standoff, so bulk-size (11-14 mm) elements there leave a piecewise-constant-M
    ripple on the map (flat-top-edge bumps g~1.02-1.04, K1g biased low by ~12%).
    face_maxh=0.006 (~standoff/3) removes the ripple (g_max 1.0001, K1g matches the
    reduced-Omega profile to 0.6%).  Do NOT push it much finer at fixed maxh_iron:
    a face/bulk size contrast of ~3.5x (face 4 mm vs bulk 14 mm) drives the
    mass-Riesz CG past its 4000-iteration cap (measured 2026-07-13)."""
    import ngsolve as ng
    from netgen.occ import OCCGeometry, Box, Pnt
    beta = np.radians(beta_deg)
    top = _fem_footprint(F_GAP / 2, F_Z_OUT, beta)
    bot = _fem_footprint(-F_Z_OUT, -F_GAP / 2, beta)
    leg_l = Box(Pnt(-F_POLE_W, -F_LEG_HY, -F_GAP / 2),
                Pnt(-F_POLE_W + F_T_LEG, F_LEG_HY, F_GAP / 2))
    leg_r = Box(Pnt(F_POLE_W - F_T_LEG, -F_LEG_HY, -F_GAP / 2),
                Pnt(F_POLE_W, F_LEG_HY, F_GAP / 2))
    iron = top + bot + leg_l + leg_r
    iron.mat("iron")
    iron.maxh = maxh_iron
    for e in iron.edges:      # same graded refinement at the gap-facing pole-face edges
        c = e.center
        if abs(abs(c.z) - F_GAP / 2) < 1e-4 and abs(c.x) < 0.11 and abs(c.y) > 0.05:
            e.maxh = edge_maxh
    if face_maxh is not None:
        for f in iron.faces:
            c = f.center
            if abs(abs(c.z) - F_GAP / 2) < 1e-4:
                f.maxh = face_maxh
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(iron).GenerateMesh(maxh=maxh_iron))
    return mesh


def hdiv_solve_midplane(beta_deg, mu_r=1000.0, maxh_iron=0.014, edge_maxh=0.004,
                        face_maxh=None, nx=81, ny=401, xmax=0.05, ymax=0.26):
    """FEEC HDiv-VIM twin of ``fem_solve_midplane``: ``radia.vim.MeshSoftIron`` on the
    iron-only mesh + the same explicit CoilBuilder pair, ``rad.Solve`` auto-dispatch
    (RT1), then one batch ``rad.Fld`` mid-plane map (all points in the gap/air, so the
    analytic integrals are exact for the solved piecewise-constant M).  This is fully
    engine-independent of the reduced-Omega solve, which makes it a useful cross-check
    for separating the model deficit from discretization error (2026-07-13)."""
    import radia as rad
    from radia import vim
    rad.UtiDelAll()                                   # also clears the HDiv registry
    coils = fem_build_coil(float(np.radians(beta_deg)))
    mesh = hdiv_build_iron_mesh(beta_deg, maxh_iron, edge_maxh, face_maxh)
    iron = vim.MeshSoftIron(mesh, mu_r=mu_r)
    top = rad.ObjCnt([iron, coils])
    rad.Solve(top)                                    # auto -> FEEC HDiv-VIM (RT1)
    xs = np.linspace(-xmax, xmax, nx)
    ys = np.linspace(-ymax, ymax, ny)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
    Bz = np.asarray(rad.Fld(top, "b", pts))[:, 2].reshape(X.shape)
    ne = int(mesh.ne)
    rad.UtiDelAll()
    return xs, ys, Bz, ne


def hdiv_scoff_study(betas_deg=(0.0, 20.0), rho=5.0, maxh_iron=0.014, edge_maxh=0.004,
                     face_maxh=None):
    """``fem_scoff_study`` with the field engine swapped to the FEEC HDiv-VIM.

    Cross-check result (2026-07-13, committed in edge_focusing_fem_results.json
    `hdiv_vim_cross_check`): dK_in agrees with the reduced-Omega engine to 0.8% at
    matched edge-mesh density (absolute-dK scatter across engines/meshes ~+-3%),
    and dK_in/model = 0.92-0.95 in EVERY configuration -- which REASSIGNS the
    -5..-7% model deficit from "iron-mesh discretization error" to REAL 3D physics:
    the local iso-field tilt of the entrance fringe near x=0 is only ~0.95-0.96 of
    the geometric tan(beta) (corner arcs + side bars + finite pole width), measured
    directly on both engines' maps, so the x-uniform tilted-fringe model (and
    hard-edge/SCOFF bookkeeping with the geometric beta) overpredicts the edge
    focusing by ~5% for this geometry.

    Historical mesh-dependence diagnosis (separation runs, 2026-07-13, recorded
    in the JSON's `mesh_dependence_diagnosis`): the then-used write-back field had
    a piecewise-constant-M ripple over bulk-size GAP-FACING pole-face elements at
    20 mm standoff -- NOT the edge lines.  The current production `rad.Fld`
    redirects solved HDiv objects to the full RT1 C++ evaluator, so this collapse
    is no longer its field path.  In the historical run, face_maxh=0.006
    (~standoff/3) cures it (g_max 1.0245 -> 1.0001; K1g 7.8 -> 8.80 mm, matching
    reduced-Omega 8.85 to 0.6%); face 4 mm at bulk 14 mm exceeds the mass-Riesz CG
    limit (4000-iter non-convergence).  Residual: B0 scatters +-3.5% across all
    configurations (0.425..0.456 T) while the global demag stays constant to 5
    digits -- a near-field-evaluation sensitivity that first-order-cancels in the
    B0-normalized dK measurement; do not quote the HDiv B0 absolute value."""
    return fem_scoff_study(betas_deg, rho,
                           solve_midplane=lambda bd: hdiv_solve_midplane(
                               bd, maxh_iron=maxh_iron, edge_maxh=edge_maxh,
                               face_maxh=face_maxh))


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
    enge = np.array([r["enge"] for r in sweep])
    ax[0].plot(tb, law, color=BLUE, lw=3.0, label=r"hard edge $\tan\beta/\rho$")
    ax[0].plot(tb, enge, color="#1a2230", ls="--", lw=1.6,
               label=r"Enge $\tan(\beta-\psi)/\rho$")
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
    print("beta[deg]   tracked 1/f_z    hard tan/rho    enge tan(b-psi)/rho   rel(enge)")
    for r in sw:
        rel = (abs(r["inv_fz"] - r["enge"]) / abs(r["enge"])
               if abs(r["enge"]) > 1e-9 else float("nan"))
        print(f"  {r['beta_deg']:5.1f}   {r['inv_fz']:+.5f}      {r['hard_edge']:+.5f}"
              f"        {r['enge']:+.5f}          {rel*100:6.2f}%")
    wc = w_convergence()
    print("\nw-convergence (slope of 1/f_z vs the hard-edge law):")
    for d in wc:
        print(f"  w={d['w']:.3f}: slope={d['slope_vs_law']:.4f}")
    s = summarize(sw, wc)
    print(f"\nmax rel err: vs hard edge {s['max_rel_err_vs_law']*100:.2f}%, "
          f"vs Enge {s['max_rel_err_vs_enge']*100:.2f}%; "
          f"beta=0 baseline {s['beta0_baseline']:.2e} (= -K1g/rho^2), "
          f"finest-w slope {s['finest_w_slope']:.4f}")
