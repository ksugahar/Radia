"""Rung 1.5: the COLLECTING SynRM channel -- where the hodograph earns its keep.

Rung 1 established that the pure turn (flux enters only through the end faces)
is exactly solvable for any material law (H = C/r), so the hodograph was
validated there but not needed.  The real SynRM channel is different: flux
enters DISTRIBUTED along the gap-side face and accumulates.  That kills the
azimuthal-H structure -- no quadrature exists -- and the design problem
becomes genuinely two-dimensional in the hodograph.

New boundary-condition class exercised here: the entry face is prescribed in
the hodograph as the segment { B = B_e, theta in [0, theta_c] } carrying the
DIRICHLET RAMP A(theta) = Phi * theta / theta_c (uniform collection per unit
field angle).  Wall assignment is forced by a non-crossing argument in the
(B, theta) plane: the CAP wall must be the LONG wall attached at the theta=0
corner of the entry face (A = 0); the opposite barrier wall (A = Phi) is the
short wall attached at theta = theta_c.  The reversed assignment makes the
two wall images cross and the domain pinch.

Domain (B horizontal, theta vertical; all walls straight segments there):

    entry   : B = B_e,            theta in [0, theta_c],  A = Phi*theta/theta_c
    cap wall: A = 0,   B ramps B_e -> B_cap over [0, theta_r], flat to theta_m
    low wall: A = Phi, B ramps B_e -> B_out1 over [theta_c, theta_b],
              flat over [theta_b, theta_m]
    exit    : theta = theta_m, natural (Psi = const, the q-axis symmetry face)

Verification: mesh the recovered outline and solve the independent nonlinear
FEM with the SAME distributed entry (the design's A distribution along the
recovered face, imposed as Dirichlet data), walls at their design A values,
natural exit.  Compare wall |B| against the prescribed profiles and the MMF.

Baseline: a compass-drawn collector THROUGH THE SAME ENTRY FACE, same exit
midpoint/width/direction, walls = circular arcs -- what an engineer sketches.
Same flux, same entry data.  Its peak wall |B| vs the cap is the payoff metric.

Golden bands asserted at the end of the run (2026-07-28 baseline, LAB):
  orientation            : J single-signed on every sampled wall point
  entry Dirichlet data   : A(t) chord-polynomial fit residual < 1e-3 of Phi
  cap wall vs profile    : mean < 1.5 % (full profile incl. the ramp)
  flat-cap region        : max < 3.0 %
  MMF                    : rel < 0.5 % x (mu_s/mu_d at the low-wall level) --
                           MMF is an H-quantity, and on the saturating curve
                           dH/H = (mu_s/mu_d) dB/B ~ 4.3x amplification here

Measured at that baseline (h/16): flat-cap region mean 0.032 % / max 0.068 %
with peak 1.900 T; full cap profile mean 0.20 %; low wall mean 0.34 %; MMF
design 1.658 A (prescription-direct; psi-projection route agrees to 0.15 %)
vs FEM 1.683 A = the 0.35 % |B| agreement seen through the 4.26x saturation
slope.  Entry density dA/ds = 0.868..0.884 T (near-uniform, ~ B_e).  The
compass baseline through the SAME face and exit uses 3.29 mm^2 against the
design's 1.84 mm^2 (+78 %) and peaks at 1.371 T -- 28 % of the cap unused:
with fixed terminals, circular walls cannot follow the accumulating flux.

Design rules learned here (each cost one failed iteration; see the
hodograph-wall-cusp bug pattern):
  1. Wall assignment is forced by non-crossing in (B, theta): the cap wall is
     the LONG wall attached at the theta = 0 corner of the entry face.
  2. Keep rho_local = H(B_capwall)/H(B_lowwall) <= ~5 at every theta (the
     Rung-1 chart reused as a LOCAL rule) or the wall demands a ~50 um
     turning radius and cusps.
  3. Ramps must be C1 -- wall-advance speed is |Psi_theta + Psi_B B'|/q, so a
     B' discontinuity puts a cusp exactly at the kink.  sin(pi t/2T) ramps,
     feature angles staggered.

Run:  python verify_synrm_collector_design.py
Writes results_synrm_collector_design.json and synrm_collector_outlines.png
next to this file (committed).
"""
import datetime
import json
import math
import os
import platform
import sys

import numpy as np
from scipy.spatial import cKDTree
from ngsolve import (
    BND, BilinearForm, CoefficientFunction, GridFunction, H1, Mesh,
    SetNumThreads, TaskManager, cos, dx, grad, sin, sqrt, x, y,
)
from netgen.geom2d import SplineGeometry

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_ipm_bridge_free_boundary import (            # noqa: E402
    B_KNEE, DFLUX, log, mu_d_of, mu_s_of, nu_cf_of,
    polygon_area, recover_potential, _self_intersections,
)

# ---------------- design case ----------------
# Wall profiles are chosen so the LOCAL field contrast between the walls,
# rho_local(theta) = H(B_capwall)/H(B_lowwall), stays <= ~5 everywhere -- the
# Rung-1 chart's flat-optimum region.  The first attempt (cap ramp done by
# 30 deg while the low wall still sat at B_e) drove rho_local to ~17, which
# demands a ~50 um local turning radius: the cap wall curled into a cusp and
# the outline self-intersected.  The contrast bound is a REUSABLE design rule.
B_E = 0.90                        # iron-side |B| along the entry face
B_CAP = B_KNEE                    # 1.90 T
B_OUT1 = 1.45                     # low wall level at the exit (~rho=3.3 contrast)
TH_C = math.radians(30.0)         # collection span (field angle)
TH_R = math.radians(65.0)         # cap wall reaches the cap here
TH_LB = math.radians(55.0)        # low wall ramp ends (starts at TH_C)
TH_B = math.radians(100.0)        # end of the designed body (before margin)
TH_M = math.radians(120.0)        # exit (incl. 20 deg lead-out margin)
PHI = DFLUX                       # 1.1e-3 Wb per metre of stack
NS = 201                          # samples per boundary piece
EPS = 1e-9
EPSB = 5e-4                       # B-inset > polyline sagitta (~2e-4)
CTRIM = 0.008                     # rad; skip the zero-width corner wedge tip

# Wall-advance speed along an A = const wall is |Psi_theta + Psi_B * B'|/q, so
# a DISCONTINUOUS ramp slope B' makes the wall speed jump and the recovered
# wall curls into a cusp right at the kink (observed at the exact kink angle,
# twice).  The ramps therefore use sin(pi t / 2T): finite slope at the start
# (a clean wedge against the entry edge) and C1-flat into the plateau.  The
# two ramp ends are STAGGERED (65 vs 55 deg) so no two profile features
# coincide at one angle.


def _sramp(u):
    return np.sin(0.5 * math.pi * np.clip(u, 0.0, 1.0))


def cap_B(t):
    t = np.asarray(t, dtype=float)
    return B_E + (B_CAP - B_E) * _sramp(t / TH_R)


def low_B(t):
    t = np.asarray(t, dtype=float)
    return B_E + (B_OUT1 - B_E) * _sramp((t - TH_C) / (TH_LB - TH_C))


# ---------------- step (1): the hodograph design ----------------
def design_collector(report, maxh=0.02):
    geo = SplineGeometry()
    # boundary polylines in (B, theta); the smooth ramps are sampled densely
    tr = np.linspace(0.0, TH_R, 41)
    cap_edge = [(float(cap_B(t)), float(t)) for t in tr] + [(B_CAP, TH_M)]
    tl = np.linspace(TH_LB, TH_C, 41)
    low_edge = [(B_OUT1, TH_M)] + [(float(low_B(t)), float(t)) for t in tl]
    pieces = [(cap_edge, "inner"),
              ([(B_CAP, TH_M), (B_OUT1, TH_M)], "outlet"),
              (low_edge, "outer"),
              ([(B_E, TH_C), (B_E, 0.0)], "entry")]
    allpts, alltags = [], []
    for pts_, bc in pieces:
        for i in range(len(pts_) - 1):
            allpts.append(pts_[i])
            alltags.append(bc)
    ids = [geo.AppendPoint(*p) for p in allpts]
    for i, bc in enumerate(alltags):
        geo.Append(["line", ids[i], ids[(i + 1) % len(ids)]], bc=bc,
                   leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    log(f"design: hodograph mesh {mesh.ne} elements")

    aB = x * mu_d_of(x) / mu_s_of(x) ** 2
    bB = 1.0 / (mu_s_of(x) * x)
    qB = x / mu_s_of(x)

    fes = H1(mesh, order=3, dirichlet="inner|outer|entry")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (aB * grad(u)[0] * grad(v)[0] + bB * grad(u)[1] * grad(v)[1]) * dx
    a.Assemble()
    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"inner": 0.0, "outer": PHI,
                             "entry": PHI * y / TH_C}, default=0.0), BND)
    res = gfA.vec.CreateVector()
    res.data = -a.mat * gfA.vec
    gfA.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res

    A_B, A_t = grad(gfA)[0], grad(gfA)[1]
    Psi_B, Psi_t = -bB * A_t, aB * A_B
    gfPsi = recover_potential(mesh, CoefficientFunction((Psi_B, Psi_t)))
    c, s = cos(y), sin(y)
    Fx = CoefficientFunction(((Psi_B / qB) * c - (A_B / x) * s,
                              (Psi_t / qB) * c - (A_t / x) * s))
    Fy = CoefficientFunction(((Psi_B / qB) * s + (A_B / x) * c,
                              (Psi_t / qB) * s + (A_t / x) * c))
    gfx = recover_potential(mesh, Fx)
    gfy = recover_potential(mesh, Fy)
    Jcf = grad(gfx)[0] * grad(gfy)[1] - grad(gfx)[1] * grad(gfy)[0]
    log("design: coordinates recovered")

    def sample(curve_pts):
        out = []
        for bb, tt in curve_pts:
            p = mesh(float(bb), float(tt))
            out.append((gfx(p), gfy(p)))
        return np.array(out)

    ths_cap = np.linspace(CTRIM, TH_M - EPS, NS)
    cap_pts = sample([(cap_B(t) - EPSB, t) for t in ths_cap])
    ths_low = np.linspace(TH_C + EPS, TH_M - EPS, NS)
    low_pts = sample([(low_B(t) + EPSB, t) for t in ths_low])
    ths_ent = np.linspace(CTRIM, TH_C - EPS, NS)
    ent_pts = sample([(B_E + EPSB, t) for t in ths_ent])
    ent_A = PHI * ths_ent / TH_C
    bs_exit = np.linspace(B_OUT1 + EPS, B_CAP - EPS, 80)
    exit_pts = sample([(bb, TH_M - EPS) for bb in bs_exit])

    # orientation monitor on all walls
    Jv = []
    for bb, tt in ([(cap_B(t) - EPSB, t) for t in ths_cap]
                   + [(low_B(t) + EPSB, t) for t in ths_low]
                   + [(B_E + EPSB, t) for t in ths_ent]):
        Jv.append(Jcf(mesh(float(bb), float(tt))))
    Jv = np.array(Jv)

    # design MMF along the low wall (flux line: Delta Psi = int H.dl)
    mmf = float(gfPsi(mesh(low_B(TH_M - EPS) + EPSB, TH_M - EPS))
                - gfPsi(mesh(B_E + EPSB, TH_C + EPS)))
    psi_low = np.array([gfPsi(mesh(float(low_B(t)) + EPSB, float(t)))
                        for t in ths_low])
    # Direct design MMF: the low wall is a flux line, so the tangential |H|
    # ON it is exactly q(B_low(theta)) from the prescription -- no potential
    # recovery involved (gfPsi is a Galerkin projection of a field that is a
    # gradient only for the exact solution, so it carries O(FE) error).
    _mid_t = 0.5 * (ths_low[:-1] + ths_low[1:])
    _q_mid = np.array([b / mu_s_of(b) for b in low_B(_mid_t)])
    _seg = np.linalg.norm(np.diff(low_pts, axis=0), axis=1)
    mmf_q_cum = np.concatenate([[0.0], np.cumsum(_q_mid * _seg)])

    # entry-face physics: arc length, local entry density dA/ds, face angle
    seg = np.linalg.norm(np.diff(ent_pts, axis=0), axis=1)
    s_ent = np.concatenate([[0.0], np.cumsum(seg)])
    dAds = np.gradient(ent_A, s_ent)
    face_len = float(s_ent[-1])

    loop = np.vstack([cap_pts, exit_pts[::-1][1:], low_pts[::-1][1:],
                      ent_pts[::-1][1:-1]])
    area = abs(polygon_area(loop[:, 0], loop[:, 1]))
    d = {
        "J_single_sign": bool(np.all(Jv < 0) or np.all(Jv > 0)),
        "min_absJ": float(np.min(np.abs(Jv))),
        "median_absJ": float(np.median(np.abs(Jv))),
        "mmf_design_A": mmf,
        "iron_area_mm2": 1e6 * area,
        "entry_face_len_mm": 1e3 * face_len,
        "entry_dAds_T": [float(dAds.min()), float(dAds.mean()),
                         float(dAds.max())],
        "exit_width_mm": 1e3 * float(np.linalg.norm(cap_pts[-1] - low_pts[-1])),
    }
    print(f"  [design] J single sign={d['J_single_sign']}  "
          f"min|J|={d['min_absJ']:.3e} (median {d['median_absJ']:.3e})")
    print(f"  [design] entry face {d['entry_face_len_mm']:.4f} mm, "
          f"entry density dA/ds = {dAds.min():.3f}/{dAds.mean():.3f}/"
          f"{dAds.max():.3f} T (min/mean/max)")
    print(f"  [design] exit width {d['exit_width_mm']:.4f} mm, iron "
          f"{d['iron_area_mm2']:.4f} mm^2, MMF {mmf:.4f} A")
    report["design"] = d
    return {"cap": cap_pts, "low": low_pts, "ent": ent_pts, "exit": exit_pts,
            "ths_cap": ths_cap, "ths_low": ths_low, "ent_A": ent_A,
            "s_ent": s_ent, "psi_low": psi_low, "mmf_q_cum": mmf_q_cum}


# ---------------- entry Dirichlet data as a closed-form CF ----------------
def entry_A_cf(ent_pts, ent_A, tol=1e-3):
    """A along the face as a polynomial in the chord coordinate (fail loud)."""
    p0, p1 = ent_pts[0], ent_pts[-1]
    tvec = (p1 - p0) / np.linalg.norm(p1 - p0)
    t = (ent_pts - p0) @ tvec
    if not np.all(np.diff(t) > 0):
        raise RuntimeError("entry face folds against its chord; the chord "
                           "parameterization is invalid")
    coef = np.polyfit(t, ent_A, 5)
    fit = np.polyval(coef, t)
    resid = np.abs(fit - ent_A).max() / PHI
    if resid > tol:
        raise RuntimeError(f"entry A(t) poly fit residual {resid:.2e} > {tol}")
    tcf = (x - float(p0[0])) * float(tvec[0]) + (y - float(p0[1])) * float(tvec[1])
    acf = CoefficientFunction(0.0)
    for ck in coef:
        acf = acf * tcf + float(ck)
    return acf, resid


# ---------------- step (2): independent nonlinear FEM ----------------
def forward_collector(curves, acf, maxh, label):
    cap_pts, low_pts, ent_pts, exit_pts = (curves["cap"], curves["low"],
                                           curves["ent"], curves["exit"])
    # the cap wall ends at the CAP exit corner, so the exit face must be
    # traversed cap-corner -> low-corner, i.e. REVERSED relative to sampling.
    # Guard the circuit: a piece concatenated in the wrong direction produces
    # dozens of phantom "self-intersections" that mimic a genuine design fold
    # (that misdiagnosis cost three debugging rounds).
    seq = [cap_pts, exit_pts[::-1], low_pts[::-1], ent_pts[::-1]]
    span = np.linalg.norm(np.vstack(seq).max(axis=0)
                          - np.vstack(seq).min(axis=0))
    for k in range(len(seq)):
        gap = np.linalg.norm(seq[k][-1] - seq[(k + 1) % len(seq)][0])
        if gap > 0.1 * span:
            raise RuntimeError(
                f"{label}: outline pieces {k} and {(k + 1) % len(seq)} do not "
                f"connect (seam gap {gap:.3e} vs extent {span:.3e}) -- a "
                f"piece is traversed in the wrong direction")
    parts = [(cap_pts, "wall_cap"), (exit_pts[::-1][1:], "outlet"),
             (low_pts[::-1][1:], "wall_low"), (ent_pts[::-1][1:-1], "entry")]
    loop = np.vstack([p for p, _ in parts])
    tags = [t for p, t in parts for _ in range(len(p))]
    tol = 0.05 * maxh
    keep = [0]
    for i in range(1, len(loop)):
        if np.linalg.norm(loop[i] - loop[keep[-1]]) > tol:
            keep.append(i)
    loop = loop[keep]
    tags = [tags[i] for i in keep]
    n = len(loop)
    hits = _self_intersections(loop)
    if hits:
        raise RuntimeError(f"{label}: outline self-intersects at {hits[:5]} "
                           f"({len(hits)} total)")
    if polygon_area(loop[:, 0], loop[:, 1]) < 0:
        loop = loop[::-1]
        tags = [tags[(n - 2 - i) % n] for i in range(n)]

    geo = SplineGeometry()
    ids = [geo.AppendPoint(float(px), float(py)) for px, py in loop]
    for i in range(n):
        geo.Append(["line", ids[i], ids[(i + 1) % n]], bc=tags[i],
                   leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    log(f"{label}: mesh {mesh.ne} elements, {mesh.nv} vertices")

    fes = H1(mesh, order=3, dirichlet="wall_cap|wall_low|entry")
    u, v = fes.TnT()
    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"wall_cap": 0.0, "wall_low": PHI, "entry": acf},
                            default=0.0), BND)
    gfPrev, gfTrial = GridFunction(fes), GridFunction(fes)
    eps = 1e-6
    omega = 0.35
    hist = []
    for _ in range(600):
        gfPrev.vec.data = gfA.vec
        a = BilinearForm(fes)
        a += nu_cf_of(sqrt(grad(gfA) * grad(gfA) + eps ** 2)) \
            * grad(u) * grad(v) * dx
        a.Assemble()
        gfTrial.vec.data = gfA.vec
        r = gfA.vec.CreateVector()
        r.data = -a.mat * gfTrial.vec
        gfTrial.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                          inverse="sparsecholesky") * r
        gfA.vec.data = gfPrev.vec + omega * (gfTrial.vec - gfPrev.vec)
        dvec = gfPrev.vec.CreateVector()
        dvec.data = gfA.vec - gfPrev.vec
        hist.append(dvec.Norm() / max(gfA.vec.Norm(), 1e-30))
        if hist[-1] < 1e-9:
            break
    log(f"{label}: Picard {len(hist)} iterations, final rel step {hist[-1]:.2e}")
    if hist[-1] > 1e-7:
        raise RuntimeError(f"{label}: Picard did not converge "
                           f"(rel step {hist[-1]:.2e})")

    # wall probes: inset along the INWARD wall normal (side chosen by a
    # point-in-polygon test; a nearest-opposite-point chord can exit the
    # domain on the concave side of a strongly curved wall)
    from matplotlib.path import Path
    poly = Path(loop)
    Bcf = sqrt(grad(gfA) * grad(gfA))

    def inset_points(pts, gap_ref, frac):
        tree = cKDTree(gap_ref)
        dist, _ = tree.query(pts)
        tang = np.gradient(pts, axis=0)
        tang /= np.maximum(np.linalg.norm(tang, axis=1)[:, None], 1e-30)
        nrm = np.c_[-tang[:, 1], tang[:, 0]]
        out = []
        for k, (px, py) in enumerate(pts):
            step = frac * dist[k]
            for sgn in (1.0, -1.0):
                q = np.array([px, py]) + sgn * step * nrm[k]
                if poly.contains_point(q):
                    out.append(q)
                    break
            else:
                raise RuntimeError(f"{label}: no inward normal at "
                                   f"({px:.6e}, {py:.6e})")
        return np.array(out)

    def probe(cf, pts, what):
        vals = []
        for qx, qy in pts:
            try:
                vals.append(cf(mesh(float(qx), float(qy))))
            except Exception as exc:                     # noqa: BLE001
                raise RuntimeError(f"{label}: {what} probe ({qx:.6e}, "
                                   f"{qy:.6e}) outside the mesh ({exc})") from exc
        return np.array(vals)

    PT = 6                      # corner trim: probes stay off the chamfers
    b_cap_fem = probe(Bcf, inset_points(cap_pts[PT:-PT],
                                        np.vstack([low_pts, ent_pts]), 0.005),
                      "cap wall")
    b_low_fem = probe(Bcf, inset_points(low_pts[PT:-PT],
                                        np.vstack([cap_pts, ent_pts]), 0.005),
                      "low wall")

    # MMF along the low wall (flux line): int H.dl over the trimmed path.
    # H varies ACROSS the channel (toward the cap wall it triples), so a probe
    # a fixed fraction inside reads high; integrate at two insets and
    # extrapolate linearly to the wall.
    Bvec = CoefficientFunction((grad(gfA)[1], -grad(gfA)[0]))
    Hvec = nu_cf_of(sqrt(grad(gfA) * grad(gfA) + eps ** 2)) * Bvec
    lp = low_pts[PT:-PT]
    mids = 0.5 * (lp[:-1] + lp[1:])
    segs = lp[1:] - lp[:-1]

    def mmf_at(frac):
        hq = inset_points(mids, np.vstack([cap_pts, ent_pts]), frac)
        tot = 0.0
        for i in range(len(segs)):
            hv = Hvec(mesh(float(hq[i][0]), float(hq[i][1])))
            tot += hv[0] * segs[i][0] + hv[1] * segs[i][1]
        return abs(float(tot))

    m1, m2 = mmf_at(0.01), mmf_at(0.005)
    mmf = 2.0 * m2 - m1                     # linear extrapolation to the wall
    return {"n_elements": mesh.ne, "picard_iterations": len(hist),
            "mmf_fem_at_1pct_A": m1, "mmf_fem_at_05pct_A": m2,
            "mmf_fem_A": mmf}, b_cap_fem, b_low_fem


# ---------------- step (3): the compass-drawn baseline ----------------
def arc_through(P, Q, t_end, nsamp):
    """Circular arc from P to Q whose tangent AT Q is t_end (unit)."""
    n_end = np.array([-t_end[1], t_end[0]])       # normal at Q
    d = P - Q
    denom = 2.0 * float(d @ n_end)
    if abs(denom) < 1e-12 * np.linalg.norm(d):
        tt = np.linspace(1.0, 0.0, nsamp)[:, None]
        return Q + tt * d                          # degenerate: straight P->Q
    R = float(d @ d) / denom                       # signed: centre C = Q + R n
    C = Q + R * n_end
    rad = abs(R)
    a0 = math.atan2(P[1] - C[1], P[0] - C[0])
    a1 = math.atan2(Q[1] - C[1], Q[0] - C[0])
    da = a1 - a0
    while da > math.pi:
        da -= 2 * math.pi
    while da < -math.pi:
        da += 2 * math.pi
    angs = a0 + da * np.linspace(0.0, 1.0, nsamp)
    return C + rad * np.c_[np.cos(angs), np.sin(angs)]


def naive_collector(curves):
    """Same entry face, same exit midpoint/width/direction; arc walls."""
    cap_pts, low_pts, ent_pts, exit_pts = (curves["cap"], curves["low"],
                                           curves["ent"], curves["exit"])
    E_cap, E_low = ent_pts[0], ent_pts[-1]         # face ends (A=0 / A=Phi)
    X_cap, X_low = cap_pts[-1], low_pts[-1]        # design exit corners
    t_exit = X_cap - X_low
    t_exit = np.array([-t_exit[1], t_exit[0]])
    t_exit /= np.linalg.norm(t_exit)
    # orient the exit tangent along the channel (away from the entry)
    if float(t_exit @ (0.5 * (X_cap + X_low) - 0.5 * (E_cap + E_low))) < 0:
        t_exit = -t_exit
    ncap = arc_through(E_cap, X_cap, t_exit, len(cap_pts))
    nlow = arc_through(E_low, X_low, t_exit, len(low_pts))
    return {"cap": ncap, "low": nlow, "ent": ent_pts, "exit": exit_pts,
            "ths_cap": curves["ths_cap"], "ths_low": curves["ths_low"]}


# ---------------- main ----------------
def main():
    SetNumThreads(4)
    failures = []
    report = {"case": {
        "B_entry_T": B_E, "B_cap_T": B_CAP, "B_out1_T": B_OUT1,
        "collection_deg": math.degrees(TH_C),
        "cap_ramp_end_deg": math.degrees(TH_R),
        "low_ramp_end_deg": math.degrees(TH_LB),
        "body_end_deg": math.degrees(TH_B), "exit_deg": math.degrees(TH_M),
        "design_rule": "keep rho_local = H(B_cap_wall)/H(B_low_wall) <= ~5",
        "flux_Wb_per_m": PHI,
        "material": "same representative curve as the promoted bridge driver",
    }}
    with TaskManager():
        print("step (1) hodograph design of the collecting channel")
        curves = design_collector(report)
        if not report["design"]["J_single_sign"]:
            failures.append("design: inverse map folds (J changes sign)")

        acf, resid = entry_A_cf(curves["ent"], curves["ent_A"])
        report["entry_fit_residual"] = resid
        print(f"  [entry] A(t) poly fit residual {resid:.2e} (of Phi)")

        w_exit = float(np.linalg.norm(curves["cap"][-1] - curves["low"][-1]))
        print("\nstep (2) independent nonlinear FEM on the designed outline")
        report["verify"] = {}
        for div, lab in ((8.0, "design_h8"), (16.0, "design_h16")):
            res, bc_, bl_ = forward_collector(curves, acf, w_exit / div, lab)
            # probes are corner-trimmed inside forward_collector (PT = 6)
            thc = curves["ths_cap"][6:-6]
            thl = curves["ths_low"][6:-6]
            tgt_c = cap_B(thc)
            tgt_l = low_B(thl)
            ec = np.abs(bc_ - tgt_c) / tgt_c
            el = np.abs(bl_ - tgt_l) / tgt_l
            flat = (thc >= TH_R + math.radians(5.0)) & (thc <= TH_B)
            efc = np.abs(bc_[flat] - B_CAP) / B_CAP
            res.update({
                "cap_wall_rel_err_mean": float(ec.mean()),
                "cap_wall_rel_err_max": float(ec.max()),
                "cap_flat_rel_err_mean": float(efc.mean()),
                "cap_flat_rel_err_max": float(efc.max()),
                "cap_flat_peak_T": float(bc_[flat].max()),
                "low_wall_rel_err_mean": float(el.mean()),
                "low_wall_rel_err_max": float(el.max()),
                "cap_B_fem_T": bc_.tolist(), "low_B_fem_T": bl_.tolist(),
            })
            print(f"  [{lab}] cap wall vs prescribed profile: mean "
                  f"{ec.mean()*100:.3f}% max {ec.max()*100:.3f}%   "
                  f"flat-cap region: mean {efc.mean()*100:.3f}% "
                  f"max {efc.max()*100:.3f}% peak {bc_[flat].max():.3f} T")
            print(f"  [{lab}] low wall: mean {el.mean()*100:.3f}% "
                  f"max {el.max()*100:.3f}%   MMF {res['mmf_fem_A']:.4f} A")
            report["verify"][lab] = res

        # design MMF over the SAME trimmed low-wall path as the FEM integral,
        # from the prescription directly (int q(B_low) ds); the psi-projection
        # value is kept in the JSON for reference.  NOTE the sensitivity: MMF
        # is an H-quantity, and on the saturating curve dH/H = (mu_s/mu_d) *
        # dB/B -- at the low-wall level 1.45 T that factor is ~4.3, so the
        # observed ~0.35% wall-|B| agreement corresponds to ~1.5% in MMF.
        # The MMF band is therefore the |B| band times this amplification.
        mmf_d = abs(float(curves["mmf_q_cum"][-7] - curves["mmf_q_cum"][6]))
        report["design"]["mmf_trimmed_path_A"] = mmf_d
        report["design"]["mmf_trimmed_path_psi_A"] = abs(
            float(curves["psi_low"][-7] - curves["psi_low"][6]))
        amp = float((mu_s_of(B_OUT1) / mu_d_of(B_OUT1)))
        report["design"]["mmf_sensitivity_mus_over_mud"] = amp
        for lab, v in report["verify"].items():
            rel = abs(v["mmf_fem_A"] - mmf_d) / mmf_d
            v["mmf_rel_diff"] = rel
            if v["cap_wall_rel_err_mean"] > 0.015:
                failures.append(f"{lab}: cap wall mean err "
                                f"{v['cap_wall_rel_err_mean']:.4f} > 1.5%")
            if v["cap_flat_rel_err_max"] > 0.03:
                failures.append(f"{lab}: flat-cap max err "
                                f"{v['cap_flat_rel_err_max']:.4f} > 3.0%")
            if rel > 0.005 * amp:       # |B| band 0.5% x saturation slope
                failures.append(f"{lab}: MMF design {mmf_d:.4f} vs FEM "
                                f"{v['mmf_fem_A']:.4f} A (rel {rel:.4f} > "
                                f"{0.005 * amp:.4f})")

        print("\nstep (3) compass-drawn baseline (same face, same exit, "
              "arc walls)")
        ncurves = naive_collector(curves)
        nloop = np.vstack([ncurves["cap"], ncurves["exit"][::-1][1:],
                           ncurves["low"][::-1][1:], ncurves["ent"][::-1][1:-1]])
        n_area = abs(polygon_area(nloop[:, 0], nloop[:, 1]))
        nres, nbc, nbl = forward_collector(ncurves, acf, w_exit / 16.0,
                                           "naive_h16")
        peak_n = float(max(nbc.max(), nbl.max()))
        nres.update({
            "iron_area_mm2": 1e6 * n_area,
            "peak_wall_B_T": peak_n,
            "cap_overshoot_pct": 100.0 * (peak_n - B_CAP) / B_CAP,
            "cap_B_fem_T": nbc.tolist(), "low_B_fem_T": nbl.tolist(),
        })
        print(f"  [naive] iron {1e6*n_area:.4f} mm^2 (design "
              f"{report['design']['iron_area_mm2']:.4f}), peak wall |B| = "
              f"{peak_n:.3f} T -> cap {B_CAP:.2f} T "
              f"{'EXCEEDED' if peak_n > B_CAP else 'respected'} "
              f"({nres['cap_overshoot_pct']:+.2f}%)")
        report["naive"] = nres

    # ---------------- outline picture ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    for a_, cv, ttl in ((ax[0], curves, "hodograph design"),
                        (ax[1], ncurves, "compass-drawn baseline")):
        for k, col in (("cap", "tab:red"), ("low", "tab:blue"),
                       ("ent", "tab:green"), ("exit", "0.4")):
            a_.plot(1e3 * cv[k][:, 0], 1e3 * cv[k][:, 1], col, lw=1.4)
        a_.set_title(ttl)
        a_.set_aspect("equal")
        a_.set_xlabel("x [mm]")
    ax[0].set_ylabel("y [mm]")
    png = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "synrm_collector_outlines.png")
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    print(f"outlines -> {png}")

    report["meta"] = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "hostname": platform.node(), "python_version": platform.python_version(),
        "purpose": "correctness validation only (no timing claims)",
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results_synrm_collector_design.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"results -> {out}")
    if failures:
        for f_ in failures:
            print("CHECK FAIL:", f_)
        raise SystemExit(1)
    print("all machinery checks passed")


if __name__ == "__main__":
    main()
