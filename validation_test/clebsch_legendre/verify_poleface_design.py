"""Saturable pole-face design: free boundary with a KNOWN hodograph image.

A high-field dipole pole must deliver a uniform gap field B0 over the
good-field region while the pole iron saturates.  In the hodograph (B, theta)
plane the face is not an unknown boundary: interface continuity against a
uniform-B0 air gap pins its image to a known curve Gamma(alpha) (alpha =
local face tilt), exactly as Kirchhoff's free-streamline theory pins the
free jet boundary to q = const.  On Gamma the exact oblique condition
dA/dPsi = mu0 cot(alpha) holds and enters the weak form as a tangential-
derivative boundary term.  One LINEAR hodograph field solve followed by
linear coordinate recovery yields the face shape; two formulation identities
are proven and locked here:

    J = -(b A_theta^2 + a A_B^2)/(q B) <= 0   (folding impossible), and
    dA/ds = -B0 cos(alpha) pointwise           (uniform-B0 support exact),

plus the exact width identity (transverse extent == PHI/B0, any solve).

Stage 1 (design battery, no nonlinear shape iteration):
  - B0 = 1.5 T control: the solve COLLAPSES to the flat pole (90 % of the
    width lands below alpha ~ 1e-3) -- below the knee the fringe-free
    optimum IS the flat face and the formulation returns it.
  - B0 = 2.0 T design (face iron mu_r ~ 109): a nontrivial 0.33 mm
    shim-like bump toward the gap; identities verified to 5e-5..1e-3.

Stage 2 (independent nonlinear FEM, half H-dipole with coil, yoke, return
leg; midplane/axis symmetry; SAME termination arc + flank for both poles):
  - designed vs FLAT face, each excited to By(0) = 2.000 T:
    good-field flatness 1.1e-3 vs 4.7e-3 (x ~ 4 better).
  - low-field control (linear iron) decomposes the mechanism honestly:
    the benefit is delivered as STATIC shim-like geometry (~3.8e-3);
    the saturation-DIFFERENTIAL share is only ~2e-4.  The designed pole
    overshoots at low field (field-specific optimum, as a saturable
    design should).

Run:  python verify_poleface_design.py
Writes results_poleface_design.json + poleface_design_fem.png (committed).
Correctness validation only -- no timing claims.
"""
import datetime
import json
import math
import os
import platform
import sys
import time

import ngsolve
import numpy as np
from scipy.optimize import brentq
from ngsolve import (
    BND, BilinearForm, CoefficientFunction, GridFunction, H1, IfPos,
    InnerProduct, Integrate, LinearForm, Mesh, SetNumThreads, TaskManager,
    cos, ds, dx, grad, sin, specialcf, sqrt, x, y,
)
from netgen.geom2d import SplineGeometry

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_ipm_bridge_free_boundary import (            # noqa: E402
    MU0, recover_potential, polygon_area, _self_intersections, log,
)

# ---------------- material: hard-tail steel model ----------------
MUR0 = 7000.0
BK = 1.0
NEXP = 6


def mu_r_of(B):
    return 1.0 + (MUR0 - 1.0) / (1.0 + (B / BK) ** NEXP)


def mu_s_of(B):
    return MU0 * mu_r_of(B)


def dmu_s_dB_of(B):
    return MU0 * (-(MUR0 - 1.0) * NEXP * B ** (NEXP - 1) / BK ** NEXP
                  / (1.0 + (B / BK) ** NEXP) ** 2)


def mu_d_of(B):
    return mu_s_of(B) ** 2 / (mu_s_of(B) - B * dmu_s_dB_of(B))


def nu_iron_np(B):
    return 1.0 / mu_s_of(B)


# ---------------- design case (globals set by run_design) ----------------
B0 = 2.0                       # gap spec [T]
ALPHA_MAX = 0.05               # face tilt at the designed edge [rad]
W_HALF = 0.055                 # designed half-width [m] (overhang ~ gap)
PHI = B0 * W_HALF              # tube flux [Wb/m]; tube scales EXACTLY w/ PHI
B_E = 1.35                     # entry loading [T]
TH_E = 0.35                    # field angle at the outer entry corner [rad]
GAP = 0.025                    # half-gap for the embedding [m]
DELTA = 2.0e-3                 # probe inset into the domain [(B,theta) units]
TRIM = 4


# ---------------- Gamma(alpha): the known face image ----------------
def gamma_B_of_alpha(al):
    """Iron-side |B| at face tilt al (unique transversal root)."""
    if al == 0.0:
        return B0

    def G(B):
        m = mu_r_of(B)
        return B * B - (B0 * math.cos(al)) ** 2 - (m * B0 * math.sin(al)) ** 2

    return brentq(G, B0 * 0.999, 8.0, xtol=1e-13, rtol=8.9e-16)


def theta_of_alpha(al, Broot):
    return al + math.atan(mu_r_of(Broot) * math.tan(al))


def build_gamma_table(al):
    al = np.asarray(al, dtype=float)
    fB = np.array([gamma_B_of_alpha(a) for a in al])
    fT = np.array([theta_of_alpha(a, b) for a, b in zip(al, fB)])
    return fB, fT


# ---------------- hodograph geometry ----------------
def build_hodo_mesh(report, maxh=0.03):
    al_geo = np.concatenate([[0.0], np.geomspace(3e-5, ALPHA_MAX, 45)])
    fB, fT = build_gamma_table(al_geo)
    B_edge, th_max = float(fB[-1]), float(fT[-1])
    log(f"design B0={B0}: Gamma edge = ({B_edge:.4f} T, "
        f"{math.degrees(th_max):.2f} deg)")

    tw = np.linspace(0.0, 1.0, 41)
    wB = B_edge + (B_E - B_edge) * tw
    wT = TH_E + (th_max - TH_E) * np.cos(0.5 * math.pi * tw)

    sym_pts = [(float(b), 0.0) for b in np.linspace(B_E, B0, 13)]
    face_pts = [(float(b), float(t)) for b, t in zip(fB, fT)]
    wall_pts = [(float(b), float(t)) for b, t in zip(wB, wT)]
    entry_pts = [(B_E, float(t)) for t in np.linspace(TH_E, 0.0, 9)]

    geo = SplineGeometry()
    allp, tags = [], []
    for pts_, bc in ((sym_pts, "sym"), (face_pts, "face"),
                     (wall_pts, "wall"), (entry_pts, "entry")):
        for i in range(len(pts_) - 1):
            allp.append(pts_[i])
            tags.append(bc)
    ids = [geo.AppendPoint(*p) for p in allp]
    for i, bc in enumerate(tags):
        geo.Append(["line", ids[i], ids[(i + 1) % len(ids)]], bc=bc,
                   leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    report.setdefault("gamma", {})[f"B0_{B0:g}"] = {
        "alpha": [float(a) for a in al_geo],
        "B_T": [float(b) for b in fB],
        "theta_rad": [float(t) for t in fT],
        "B_edge_T": B_edge, "theta_max_rad": th_max,
    }
    tables = {"al_geo": al_geo, "fB": fB, "fT": fT,
              "B_edge": B_edge, "th_max": th_max,
              "wall": np.column_stack([wB, wT]),
              "entry": np.array(entry_pts), "sym": np.array(sym_pts)}
    return mesh, tables


def piecewise_linear_cf(ts, vs, var):
    """PW-linear CF of `var` through (ts, vs); linear-from-origin below
    ts[0], clamped above ts[-1].  Robust at the alpha->0 corner (a log-fit
    extrapolates dangerously there)."""
    cf = CoefficientFunction(float(vs[-1]))
    for i in range(len(ts) - 2, -1, -1):
        t0, t1 = float(ts[i]), float(ts[i + 1])
        v0, v1 = float(vs[i]), float(vs[i + 1])
        seg = v0 + (v1 - v0) * (var - t0) / (t1 - t0)
        cf = IfPos(var - t1, cf, seg)
    return IfPos(var - float(ts[0]), cf, float(vs[0]) * var / float(ts[0]))


# ---------------- oblique-BC design solve ----------------
def solve_design(mesh, tables):
    aB = x * mu_d_of(x) / mu_s_of(x) ** 2
    bB = 1.0 / (mu_s_of(x) * x)
    qB = x / mu_s_of(x)

    al_cf = np.geomspace(1e-5, ALPHA_MAX, 60)
    cB, cT = build_gamma_table(al_cf)
    tan_alpha_cf = piecewise_linear_cf(cT, np.tan(al_cf), y)

    # tangential orientation probe: integral of tau_B over "face" must equal
    # +/-(B_edge - B0); the sign fixes sgn_t so sgn_t*tau is the CCW
    # (domain-on-the-left) direction of the conormal identity dPsi/ds.
    tau = specialcf.tangential(2)
    face_reg = mesh.Boundaries("face")
    probe = Integrate(tau[0], mesh, definedon=face_reg)
    dB_ccw = tables["B_edge"] - B0
    if abs(abs(probe) - abs(dB_ccw)) > 1e-8 + 1e-6 * abs(dB_ccw):
        raise RuntimeError(f"tangential probe mismatch: {probe} vs {dB_ccw}")
    sgn_t = 1.0 if probe * dB_ccw > 0 else -1.0

    fes = H1(mesh, order=3, dirichlet="sym|wall")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (aB * grad(u)[0] * grad(v)[0] + bB * grad(u)[1] * grad(v)[1]) * dx
    # weak form: int(a A_B v_B + b A_t v_t) - int_face v dPsi/ds = 0 with
    # dPsi/ds = (tan(alpha)/mu0) dA/ds on Gamma (s = CCW arc direction)
    a += (-sgn_t / MU0) * tan_alpha_cf * (grad(u).Trace() * tau) * v \
        * ds(definedon=face_reg)
    a.Assemble()

    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"wall": PHI}, default=0.0), BND)
    res = gfA.vec.CreateVector()
    res.data = -a.mat * gfA.vec
    gfA.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * res
    return gfA, (aB, bB, qB)


# ---------------- recovery + self-check battery ----------------
def _inset(cs, delta):
    cs = np.asarray(cs, dtype=float)
    t = np.gradient(cs, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    return cs + delta * np.column_stack([-t[:, 1], t[:, 0]])  # inward (CCW)


def recover_and_check(mesh, gfA, coefs, tables, chk, battery=True):
    aB, bB, qB = coefs
    A_B, A_t = grad(gfA)[0], grad(gfA)[1]
    Psi_B, Psi_t = -bB * A_t, aB * A_B
    c, s = cos(y), sin(y)
    Fx = CoefficientFunction(((Psi_B / qB) * c - (A_B / x) * s,
                              (Psi_t / qB) * c - (A_t / x) * s))
    Fy = CoefficientFunction(((Psi_B / qB) * s + (A_B / x) * c,
                              (Psi_t / qB) * s + (A_t / x) * c))
    gfx = recover_potential(mesh, Fx)
    gfy = recover_potential(mesh, Fy)
    Jcf = grad(gfx)[0] * grad(gfy)[1] - grad(gfx)[1] * grad(gfy)[0]

    def rsamp(cs):
        return np.array([(gfx(mesh(float(b), float(t))),
                          gfy(mesh(float(b), float(t)))) for b, t in cs])

    def asamp(cs):
        return np.array([gfA(mesh(float(b), float(t))) for b, t in cs])

    # face samples (capped below alpha_max: an inset at the top corner
    # would exit through the wall boundary)
    al_chk = np.geomspace(5e-5, ALPHA_MAX * 0.985, 120)
    hB, hT = build_gamma_table(al_chk)
    face_in = _inset(np.column_stack([hB, hT]), DELTA)
    r_face = rsamp(face_in)
    A_face = asamp(face_in)

    # station distribution: where the physical width lands in alpha (the
    # free-boundary unknown).  Below the knee this collapses to the flat
    # pole (90 % of the width below alpha ~ 1e-3).
    y_tr = np.abs(r_face[:, 1] - r_face[0, 1])
    stations = {}
    for frac in (0.25, 0.5, 0.75, 0.9):
        i = min(int(np.searchsorted(y_tr, frac * y_tr[-1])), len(al_chk) - 1)
        stations[frac] = float(al_chk[i])
    print(f"  [stations B0={B0:g}] width quartiles at alpha = "
          + ", ".join(f"{f:.0%}:{a:.2e}" for f, a in stations.items()))
    chk["stations_alpha"] = {str(k): v for k, v in stations.items()}

    curves = {"face": r_face, "alpha_chk": al_chk, "A_face": A_face}
    if not battery:
        return curves, stations

    # (1) Jacobian: analytically J = -(b A_t^2 + a A_B^2)/(qB) <= 0, so the
    # dominant sign must be negative; opposite signs are admissible only at
    # degeneracies (grad A -> 0 at the A=0/A=0 face-center corner).
    from matplotlib.path import Path
    outline_h = np.vstack([tables["sym"], np.column_stack([hB, hT]),
                           tables["wall"], tables["entry"]])
    pth = Path(outline_h)
    gB, gT = np.meshgrid(np.linspace(B_E, tables["B_edge"], 70),
                         np.linspace(0.0, tables["th_max"], 70))
    pts = np.column_stack([gB.ravel(), gT.ravel()])
    keep = pth.contains_points(pts, radius=-2.5 * DELTA)
    Jv = np.array([Jcf(mesh(float(b), float(t))) for b, t in pts[keep]])
    medJ = float(np.median(np.abs(Jv)))
    wrong = Jv[np.sign(Jv) != np.sign(np.median(Jv))]
    chk["J_dominant_negative"] = bool(np.median(Jv) < 0)
    chk["J_wrong_sign_maxabs_over_med"] = float(
        np.max(np.abs(wrong)) / medJ) if len(wrong) else 0.0
    print(f"  [check] J dominant negative={chk['J_dominant_negative']} "
          f"wrong-signed max|J|/med = {chk['J_wrong_sign_maxabs_over_med']:.1e}")

    # (2)+(3) face identities on resolved segments
    al_mid = 0.5 * (al_chk[:-1] + al_chk[1:])
    seg = np.diff(r_face, axis=0)
    ds_ = np.linalg.norm(seg, axis=1)
    ok = ds_ > 2e-5
    sl = np.zeros(len(ds_), dtype=bool)
    sl[TRIM:len(ds_) - TRIM] = True
    m_all = ok & sl
    if not m_all.any():
        raise RuntimeError("no resolved face segments at all")
    dAds = np.where(ok, np.diff(A_face) / np.where(ok, ds_, 1.0), np.nan)
    bn_err = np.abs(np.abs(dAds) - B0 * np.cos(al_mid)) / B0
    chk["face_Bn_rel_err_med"] = float(np.nanmedian(bn_err[m_all]))
    chk["face_Bn_rel_err_max"] = float(np.nanmax(bn_err[m_all]))
    print(f"  [check] |dA/ds| == B0 cos(alpha): rel err med "
          f"{chk['face_Bn_rel_err_med']:.2e} max {chk['face_Bn_rel_err_max']:.2e}")

    # (4) width identity: transverse extent == (A range)/B0 exactly
    tr_meas = r_face[-1 - TRIM, 1] - r_face[TRIM, 1]
    tr_anal = (A_face[-1 - TRIM] - A_face[TRIM]) / B0
    chk["width_meas_mm"] = float(1e3 * abs(tr_meas))
    chk["width_rel_err"] = float(abs(abs(tr_meas) - abs(tr_anal))
                                 / abs(tr_anal))
    print(f"  [check] width identity: {chk['width_meas_mm']:.4f} mm "
          f"(rel {chk['width_rel_err']:.2e})")

    # (5) entry image perpendicular to the local field
    entry_in = _inset(tables["entry"], DELTA)[1:-1]
    r_entry = rsamp(entry_in)
    th_e = tables["entry"][1:-1, 1]
    th_mid = 0.5 * (th_e[:-1] + th_e[1:])
    dseg = np.diff(r_entry, axis=0)
    dev_e = (np.degrees(np.arctan2(dseg[:, 1], dseg[:, 0]))
             - (90.0 + np.degrees(th_mid)) + 90.0) % 180.0 - 90.0
    chk["entry_perp_dev_deg_max"] = float(np.max(np.abs(dev_e)))

    # (6) sag + depth (physical outputs)
    x_ax = r_face[:, 0] - r_face[0, 0]
    chk["sag_at_edge_mm"] = float(1e3 * abs(x_ax[-1 - TRIM]))
    chk["depth_entry_mm"] = float(1e3 * abs(np.mean(r_entry[:, 0])
                                            - r_face[0, 0]))
    print(f"  [check] sag at edge {chk['sag_at_edge_mm']:.3f} mm; depth to "
          f"entry {chk['depth_entry_mm']:.1f} mm; entry perp "
          f"{chk['entry_perp_dev_deg_max']:.2f} deg")
    curves["wall"] = rsamp(_inset(tables["wall"], DELTA)[2:-2])
    curves["entry"] = r_entry
    curves["sym"] = rsamp(_inset(tables["sym"], DELTA))
    return curves, stations


def run_design(report, b0, w_half, battery=True):
    global B0, W_HALF, PHI
    B0, W_HALF, PHI = b0, w_half, b0 * w_half
    chk = {}
    mesh, tables = build_hodo_mesh(report)
    gfA, coefs = solve_design(mesh, tables)
    curves, stations = recover_and_check(mesh, gfA, coefs, tables, chk,
                                         battery=battery)
    report.setdefault("design_checks", {})[f"B0_{b0:g}"] = chk
    return curves, stations, tables


def embed_face(curves):
    """Design frame -> magnet frame, face center -> (0, GAP).  The tube sits
    at x_d > fc (entry deep) and the gap at x_d < fc: pole axis x_d -> +y_m
    (iron above the face, midplane at y_m = 0), transverse y_d -> x_m."""
    fc = curves["face"][0]
    out = {}
    for k in ("face", "wall", "entry", "sym"):
        out[k] = np.column_stack([curves[k][:, 1] - fc[1],
                                  (curves[k][:, 0] - fc[0]) + GAP])
    return out


# ================= independent nonlinear FEM (half H-dipole) =============
H_POLE = 0.045
P_ROOT = GAP + H_POLE
T_YOKE = 0.077
Y_YOKE = P_ROOT + T_YOKE
X_LEG_IN, X_LEG_OUT = 0.140, 0.210
X_BOX, Y_BOX = 0.280, 0.240
R_T = 0.025
PHI_END = math.radians(70.0)
COIL = (0.098, 0.134, 0.008, 0.060)
MAXH_FACE = 1.5e-3
MAXH_GAPB = 2.5e-3
MAXH_GLOB = 0.010
ORDER = 3
PROBE_Y = 5.0e-4
X_GFR = 0.030
NI_LIN = 1.2e4
NU0 = 1.0 / MU0
XS = np.linspace(1e-4, 0.058, 117)


def face_polyline(face_emb, designed):
    xs = face_emb[:, 0].copy()
    ys = face_emb[:, 1].copy()
    if not designed:
        ys = np.full_like(xs, GAP)
    if xs[0] > 1e-6:
        xs = np.concatenate([[0.0], xs])
        ys = np.concatenate([[ys[0]], ys])
    # dedup: the hodograph corner degeneracy maps many small-alpha stations
    # to nearly the same physical point (1e-12 m segments crash netgen)
    pts = np.column_stack([xs, ys])
    keep = [0]
    for i in range(1, len(pts)):
        if np.hypot(*(pts[i] - pts[keep[-1]])) > 1e-4:
            keep.append(i)
    if keep[-1] != len(pts) - 1:
        keep.append(len(pts) - 1)
    return pts[keep]


def termination(face):
    """Tangent-matched roll-off arc + 70 deg flank up to the root plane."""
    p_e = face[-1]
    t_e = face[-1] - face[-2]
    alpha = math.atan2(-t_e[1], t_e[0])
    cx = p_e[0] + R_T * math.sin(alpha)
    cy = p_e[1] + R_T * math.cos(alpha)
    ph = np.linspace(-alpha, PHI_END, 30)[1:]
    arc = np.column_stack([cx + R_T * np.sin(ph), cy - R_T * np.cos(ph)])
    a_end = arc[-1]
    x_rt = a_end[0] + (P_ROOT - a_end[1]) / math.tan(PHI_END)
    return np.vstack([arc, [[x_rt, P_ROOT]]])


def build_magnet_mesh(face_emb, designed):
    face = face_polyline(face_emb, designed)
    term = termination(face)
    x_rt = term[-1][0]
    geo = SplineGeometry()
    P = {}

    def pt(x_, y_):
        key = (round(x_, 9), round(y_, 9))
        if key not in P:
            P[key] = geo.AppendPoint(x_, y_)
        return P[key]

    def seg(p0, p1, bc, l, r, maxh=None):
        kw = dict(bc=bc, leftdomain=l, rightdomain=r)
        if maxh is not None:
            kw["maxh"] = maxh
        geo.Append(["line", pt(*p0), pt(*p1)], **kw)

    seg((0, 0), (X_LEG_IN, 0), "mid", 1, 0, MAXH_GAPB)
    seg((X_LEG_IN, 0), (X_LEG_OUT, 0), "mid", 2, 0)
    seg((X_LEG_OUT, 0), (X_BOX, 0), "mid", 1, 0)
    seg((X_BOX, 0), (X_BOX, Y_BOX), "outer", 1, 0)
    seg((X_BOX, Y_BOX), (0, Y_BOX), "outer", 1, 0)
    seg((0, Y_BOX), (0, Y_YOKE), "outer", 1, 0)
    seg((0, Y_YOKE), (0, GAP), "outer", 2, 0)
    seg((0, GAP), (0, 0), "outer", 1, 0, MAXH_GAPB)
    fpts = [tuple(p) for p in face] + [tuple(p) for p in term]
    for i in range(len(fpts) - 1):
        bc = "face" if i < len(face) - 1 else "term"
        mh = MAXH_FACE if i < len(face) - 1 else None
        seg(fpts[i], fpts[i + 1], bc, 2, 1, mh)
    seg((x_rt, P_ROOT), (X_LEG_IN, P_ROOT), "ironbnd", 2, 1)
    seg((X_LEG_IN, P_ROOT), (X_LEG_IN, 0), "ironbnd", 2, 1)
    seg((X_LEG_OUT, 0), (X_LEG_OUT, Y_YOKE), "ironbnd", 2, 1)
    seg((X_LEG_OUT, Y_YOKE), (0, Y_YOKE), "ironbnd", 2, 1)
    x0, x1, y0, y1 = COIL
    seg((x0, y0), (x1, y0), "coilbnd", 3, 1)
    seg((x1, y0), (x1, y1), "coilbnd", 3, 1)
    seg((x1, y1), (x0, y1), "coilbnd", 3, 1)
    seg((x0, y1), (x0, y0), "coilbnd", 3, 1)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "iron")
    geo.SetMaterial(3, "coil")
    mesh = Mesh(geo.GenerateMesh(maxh=MAXH_GLOB))
    print(f"  [mesh] {'designed' if designed else 'flat'}: {mesh.ne} el, "
          f"face x_end {1e3*face[-1][0]:.2f} mm, "
          f"bump {1e3*(GAP-face[-1][1]):+.3f} mm")
    return mesh


def solve_magnet(mesh, NI, gfA=None, tol=3.0e-6, maxit=600):
    """Adaptive-damping Picard (omega halves on stall, floor 0.08)."""
    fes = H1(mesh, order=ORDER, dirichlet="outer")
    u, v = fes.TnT()
    if gfA is None or gfA.space.ndof != fes.ndof:
        gfA = GridFunction(fes)
    x0, x1, y0, y1 = COIL
    f = LinearForm(fes)
    f += (NI / ((x1 - x0) * (y1 - y0))) * v * dx("coil")
    f.Assemble()
    omega, hist = 0.35, []
    gfw = GridFunction(fes)
    for it in range(maxit):
        Bmag = sqrt(InnerProduct(grad(gfA), grad(gfA)) + 1e-12)
        nu = mesh.MaterialCF({"iron": nu_iron_np(Bmag)}, default=NU0)
        a = BilinearForm(fes)
        a += nu * InnerProduct(grad(u), grad(v)) * dx
        a.Assemble()
        gfw.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                     inverse="pardiso") * f.vec
        dn = np.linalg.norm(gfw.vec.FV().NumPy() - gfA.vec.FV().NumPy())
        rel = dn / max(np.linalg.norm(gfw.vec.FV().NumPy()), 1e-30)
        hist.append(rel)
        if rel < tol:
            gfA.vec.data = gfw.vec
            return gfA, it + 1
        if len(hist) > 8 and rel > 0.98 * hist[-8] and omega > 0.08:
            omega *= 0.5
        gfA.vec.FV().NumPy()[:] = (gfA.vec.FV().NumPy()
                                   + omega * (gfw.vec.FV().NumPy()
                                              - gfA.vec.FV().NumPy()))
    raise RuntimeError(f"Picard stalled: rel={rel:.2e} after {maxit} its")


def extract_profiles(mesh, gfA):
    g = grad(gfA)
    By0 = -g(mesh(1e-4, PROBE_Y))[0]
    sgn = -1.0 if By0 < 0 else 1.0
    By = sgn * np.array([-g(mesh(float(xx), PROBE_Y))[0] for xx in XS])
    xf = np.linspace(1e-4, 0.050, 101)
    Byf = sgn * np.array([-g(mesh(float(xx), GAP - 1.5e-3))[0] for xx in xf])
    Bi = np.array([math.hypot(*g(mesh(float(xx), GAP + 1.0e-3)))
                   for xx in xf])
    gfr = XS <= X_GFR + 1e-9
    flat = By[gfr] / By[0] - 1.0
    return {"By0_T": float(abs(By0)),
            "x_m": [float(v) for v in XS],
            "By_T": [float(v) for v in By],
            "flatness_max_gfr": float(np.max(np.abs(flat))),
            "x_face_m": [float(v) for v in xf],
            "By_faceline_T": [float(v) for v in Byf],
            "iron_B_max_near_face_T": float(np.max(Bi))}


def run_pole(face_emb, designed, report, NI_start):
    tag = "designed" if designed else "flat"
    mesh = build_magnet_mesh(face_emb, designed)
    NI, gfA = NI_start, None
    for k in range(7):
        t0 = time.time()
        with TaskManager():
            gfA, its = solve_magnet(mesh, NI, gfA)
        By0 = -grad(gfA)(mesh(1e-4, PROBE_Y))[0]
        print(f"  [{tag}] NI={NI:9.1f} A -> By(0)={By0:+.5f} T "
              f"({its} its, {time.time()-t0:.0f} s)", flush=True)
        if abs(abs(By0) - B0) < 2e-3:
            break
        NI = NI * (B0 / abs(By0)) ** 1.2
    else:
        raise RuntimeError(f"{tag}: NI search did not converge")
    res = extract_profiles(mesh, gfA)
    res.update({"NI_A": float(NI), "ne": int(mesh.ne)})
    print(f"  [{tag}] flatness max |dB/B| (GFR) {res['flatness_max_gfr']:.3e}"
          f"; iron near face {res['iron_B_max_near_face_T']:.2f} T")
    report[tag] = res
    return NI, mesh


# ---------------- golden bands (locked 2026-07-29) ----------------
GOLD = {
    "collapse_station90_max": 1.5e-3,       # measured 6.8e-4 at B0=1.5
    "design_station90": (2.0e-3, 7.0e-3),   # measured 3.63e-3 at B0=2.0
    "width_rel_err_max": 5.0e-4,            # measured 4.8e-5
    "bn_med_max": 5.0e-3,                   # measured 9.8e-4
    "j_wrong_over_med_max": 1.0e-3,         # measured 0.0
    "sag_mm": (0.28, 0.40),                 # measured 0.334
    "depth_mm": (50.0, 72.0),               # measured 60.7
    "entry_perp_deg_max": 1.0,              # measured 0.24
    "fem_flat_designed": (5.0e-4, 2.0e-3),  # measured 1.12e-3
    "fem_flat_flat": (3.0e-3, 6.5e-3),      # measured 4.73e-3
    "fem_ratio_min": 2.5,                   # measured 4.23
    "fem_lin_designed_max": 2.0e-3,         # measured 7.8e-4 (overshoot)
    "fem_sat_share_max": 1.0e-3,            # measured 2.06e-4
    "fem_geom_share": (2.5e-3, 5.5e-3),     # measured 3.82e-3
    "ni_premium": (0.002, 0.035),           # measured 0.0144
    "iron_face": (2.0, 2.6),                # measured 2.25-2.26
}


def _band(name, val, lo, hi):
    if not (lo <= val <= hi):
        raise AssertionError(f"golden {name}: {val:.4e} outside "
                             f"[{lo:.3e}, {hi:.3e}]")
    print(f"  [golden] {name} = {val:.4e} in [{lo:.3e}, {hi:.3e}]  OK")


def main():
    SetNumThreads(4)
    here = os.path.dirname(os.path.abspath(__file__))
    report = {
        "schema": "radia.validation.clebsch-poleface-design.v1",
        "case": {
            "material": f"mu_r(B)=1+({MUR0:.0f}-1)/(1+(B/{BK})^{NEXP})",
            "design": {
                "B0_T": 2.0,
                "half_width_m": 0.055,
                "alpha_max_rad": ALPHA_MAX,
                "B_entry_T": B_E,
                "theta_entry_rad": TH_E,
            },
            "collapse_control": {"B0_T": 1.5, "half_width_m": 0.040},
            "fem": {
                "half_gap_m": GAP,
                "gfr_halfwidth_m": X_GFR,
                "order": ORDER,
                "H_POLE": H_POLE,
                "T_YOKE": T_YOKE,
                "R_T": R_T,
                "PHI_END_deg": math.degrees(PHI_END),
                "COIL": COIL,
                "NI_LIN_A": NI_LIN,
            },
        },
    }
    with TaskManager():
        print("step (1) collapse control at B0 = 1.5 T (below the knee)")
        _, st15, _ = run_design(report, 1.5, 0.040, battery=False)
        _band("collapse_station90", st15[0.9], 0.0,
              GOLD["collapse_station90_max"])

        print("step (2) design at B0 = 2.0 T (oblique BC + full battery)")
        curves, st20, tables = run_design(report, 2.0, 0.055, battery=True)
        chk = report["design_checks"]["B0_2"]
        _band("design_station90", st20[0.9], *GOLD["design_station90"])
        _band("width_rel_err", chk["width_rel_err"], 0.0,
              GOLD["width_rel_err_max"])
        _band("face_Bn_med", chk["face_Bn_rel_err_med"], 0.0,
              GOLD["bn_med_max"])
        _band("J_wrong_over_med", chk["J_wrong_sign_maxabs_over_med"], 0.0,
              GOLD["j_wrong_over_med_max"])
        if not chk["J_dominant_negative"]:
            raise AssertionError("J dominant sign is not negative")
        _band("sag_mm", chk["sag_at_edge_mm"], *GOLD["sag_mm"])
        _band("depth_mm", chk["depth_entry_mm"], *GOLD["depth_mm"])
        _band("entry_perp_deg", chk["entry_perp_dev_deg_max"], 0.0,
              GOLD["entry_perp_deg_max"])
    emb = embed_face(curves)

    print("step (3) independent nonlinear FEM: designed vs flat pole")
    NI0 = B0 / MU0 * GAP * 1.05
    NI_d, mesh_d = run_pole(emb["face"], True, report, NI0)
    NI_f, mesh_f = run_pole(emb["face"], False, report, NI0)

    print("step (4) low-field control (linear iron)")
    for tag, mesh in (("designed", mesh_d), ("flat", mesh_f)):
        with TaskManager():
            gfL, itsL = solve_magnet(mesh, NI_LIN)
        res = extract_profiles(mesh, gfL)
        res["NI_A"] = NI_LIN
        report[tag + "_lin"] = res
        print(f"  [{tag}_lin] By(0)={res['By0_T']:.4f} T, flatness "
              f"{res['flatness_max_gfr']:.3e} ({itsL} its)")

    d, fl = report["designed"], report["flat"]
    gfr = XS <= X_GFR + 1e-9

    def fcurve(r):
        By = np.array(r["By_T"])
        return By[gfr] / By[0] - 1.0

    dF_tot = fcurve(d) - fcurve(fl)
    dF_geo = fcurve(report["designed_lin"]) - fcurve(report["flat_lin"])
    dF_sat = dF_tot - dF_geo
    ratio = fl["flatness_max_gfr"] / max(d["flatness_max_gfr"], 1e-30)
    report["verdict"] = {
        "flatness_designed": d["flatness_max_gfr"],
        "flatness_flat": fl["flatness_max_gfr"],
        "improvement_factor": float(ratio),
        "geometry_share_max": float(np.max(np.abs(dF_geo))),
        "saturation_share_max": float(np.max(np.abs(dF_sat))),
        "ni_premium": float(NI_d / NI_f - 1.0),
    }
    print(f"\n== flatness designed {d['flatness_max_gfr']:.3e} vs flat "
          f"{fl['flatness_max_gfr']:.3e} (x{ratio:.2f}); bump effect = "
          f"geometry {np.max(np.abs(dF_geo)):.3e} + saturation-diff "
          f"{np.max(np.abs(dF_sat)):.3e}")

    _band("fem_flat_designed", d["flatness_max_gfr"],
          *GOLD["fem_flat_designed"])
    _band("fem_flat_flat", fl["flatness_max_gfr"], *GOLD["fem_flat_flat"])
    _band("fem_ratio", ratio, GOLD["fem_ratio_min"], 1e9)
    _band("fem_lin_designed", report["designed_lin"]["flatness_max_gfr"],
          0.0, GOLD["fem_lin_designed_max"])
    _band("fem_sat_share", float(np.max(np.abs(dF_sat))), 0.0,
          GOLD["fem_sat_share_max"])
    _band("fem_geom_share", float(np.max(np.abs(dF_geo))),
          *GOLD["fem_geom_share"])
    _band("ni_premium", NI_d / NI_f - 1.0, *GOLD["ni_premium"])
    _band("iron_face_designed", d["iron_B_max_near_face_T"],
          *GOLD["iron_face"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    ax = axs[0, 0]
    g20 = report["gamma"]["B0_2"]
    ax.plot(g20["B_T"], np.degrees(g20["theta_rad"]), "r.-", ms=3,
            label="Gamma(alpha) face")
    ax.plot(tables["wall"][:, 0], np.degrees(tables["wall"][:, 1]), "b-",
            label="wall (A=Phi)")
    ax.plot(tables["entry"][:, 0], np.degrees(tables["entry"][:, 1]), "g-",
            label="entry (natural)")
    ax.plot(tables["sym"][:, 0], np.degrees(tables["sym"][:, 1]), "-",
            color="0.4", label="sym (A=0)")
    ax.set_xlabel("B [T]")
    ax.set_ylabel("theta [deg]")
    ax.legend(fontsize=7)
    ax = axs[0, 1]
    for k, col in (("face", "tab:red"), ("wall", "tab:blue"),
                   ("entry", "tab:green"), ("sym", "0.4")):
        cv = emb[k]
        ax.plot(1e3 * cv[:, 0], 1e3 * cv[:, 1], color=col, lw=1.3, label=k)
        ax.plot(-1e3 * cv[:, 0], 1e3 * cv[:, 1], color=col, lw=0.6,
                alpha=0.5)
    ax.axhline(0.0, color="k", lw=0.6, ls="--")
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.legend(fontsize=7)
    ax = axs[0, 2]
    fm = emb["face"]
    ax.plot(1e3 * fm[:, 0], 1e3 * (GAP - fm[:, 1]), "r-", lw=1.5)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("bump toward the gap [mm]")
    ax = axs[1, 0]
    for tag, col, lsty in (("designed", "tab:red", "-"),
                           ("flat", "tab:blue", "-"),
                           ("designed_lin", "tab:red", "--"),
                           ("flat_lin", "tab:blue", "--")):
        r = report[tag]
        By = np.array(r["By_T"])
        ax.plot(1e3 * np.array(r["x_m"]), 1e4 * (By / By[0] - 1.0),
                color=col, ls=lsty, lw=1.5 if lsty == "-" else 1.0,
                label=tag)
    ax.axvline(1e3 * X_GFR, color="k", lw=0.7, ls="--")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("(By/By0 - 1) x 1e4")
    ax.legend(fontsize=7)
    ax = axs[1, 1]
    for tag, col in (("designed", "tab:red"), ("flat", "tab:blue")):
        r = report[tag]
        ax.plot(1e3 * np.array(r["x_face_m"]), np.array(r["By_faceline_T"]),
                color=col, lw=1.5, label=tag)
    ax.axhline(B0, color="k", lw=0.7, ls="--")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("By below the face [T]")
    ax.legend(fontsize=7)
    ax = axs[1, 2]
    xg = 1e3 * XS[gfr]
    ax.plot(xg, 1e4 * dF_tot, "k-", lw=1.5, label="designed - flat (2 T)")
    ax.plot(xg, 1e4 * dF_geo, "g--", lw=1.2, label="geometry share (lin)")
    ax.plot(xg, 1e4 * dF_sat, "m-", lw=1.2, label="saturation share")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("dF x 1e4")
    ax.legend(fontsize=7)
    fig.tight_layout()
    png = os.path.join(here, "poleface_design_fem.png")
    fig.savefig(png, dpi=140)
    print(f"figure -> {png}")

    report["meta"] = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "ngsolve_version": ngsolve.__version__,
        "numpy_version": np.__version__,
        "purpose": "correctness validation only (no timing claims)",
    }
    out = os.path.join(here, "results_poleface_design.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"results -> {out}")
    print("ALL GOLDEN CHECKS PASSED")


if __name__ == "__main__":
    main()
