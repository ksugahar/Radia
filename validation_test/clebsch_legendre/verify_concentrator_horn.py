"""Flux-concentrator horn: hodograph design -> shape -> nonlinear FEM check.

A planar sensor flux concentrator collects ambient flux over a wide face and
funnels it into a narrow tip; the sensor sits in the gap between two facing
tips.  The design constraint is a CAP: at the maximum rated ambient field,
|B| in the iron must stay below the linearity limit B_max everywhere.  The
engineering objective is a short horn delivering a given geometric gain
under that cap.

Hodograph formulation (the cleanest domain in this program): half-horn,
mirror symmetry about the axis.  In the (B, theta) plane the domain is a
curvilinear quadrilateral, ALL-Dirichlet in A:
    centerline: theta = 0,  A = 0,        B in [B_face, B_tip]
    face      : B = B_face, A = ramp 0 -> Phi/2, theta in [0, TH_C]
    wall      : A = Phi/2,  prescribed C1 profile (B_face,TH_C)->(B_tip,TH_T)
    tip       : B = B_tip,  A = ramp 0 -> Phi/2, theta in [0, TH_T]
One LINEAR hodograph field solve plus linear coordinate recovery yields the
wall curve, which IS the horn profile.  The wall tilt theta is the local taper
angle, so the prescription "theta unwinds from
TH_C to TH_T while B climbs to the cap" is literally "taper aggressively
while the iron is cheap, straighten out as the cap approaches".

Verification (locked): mirror the recovered wall, build the FACING PAIR with
a sensor gap (mid-plane homogeneous-Neumann mirror), immerse in a uniform
ambient (Dirichlet A = B0*y on the far box, ladder-checked), solve the
independent nonlinear FEM over an ambient sweep, and compare against a
STRAIGHT TAPER with the same face width, tip width, and length:

    gain(low)      horn 10.7 vs straight 10.2  (about +5 % at the SAME
                   footprint -- air-path-dominated, so this is real signal)
    iron peak at   horn ~0.96 T vs straight ~0.99 T (the designed wall
    rated ambient  spreads the load; the straight taper concentrates at
                   the tip corner)
    ladder         gap field shifts < 1 % at a 1.6x box

Honest framing: +5 % gain at identical footprint directly multiplies sensor
sensitivity and comes without nonlinear shape iteration, with the internal
|B| controlled by construction; a hand-tuned bulged spline could plausibly
match the shape -- the hodograph's edge is the constructive cap, not
exclusivity.

Run:  python verify_concentrator_horn.py
Writes results_concentrator_horn.json + concentrator_horn_verify.png.
Correctness validation only -- no timing claims.
"""
import datetime
import json
import math
import os
import platform
import sys

import ngsolve
import numpy as np
from ngsolve import (
    BND, BilinearForm, CoefficientFunction, GridFunction, H1, Mesh,
    SetNumThreads, TaskManager, cos, dx, grad, sin, sqrt, x, y,
)
from netgen.geom2d import SplineGeometry

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_ipm_bridge_free_boundary import (            # noqa: E402
    MU0, recover_potential, polygon_area, _self_intersections, log,
)

# ---------------- material: permalloy-like ----------------
MUR0 = 2.0e4
BK = 0.75
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


def nu_cf_of(Bmag):
    return 1.0 / (MU0 * (1.0 + (MUR0 - 1.0) / (1.0 + (Bmag / BK) ** NEXP)))


# ---------------- design case ----------------
B_F = 0.12                     # iron-side |B| on the collection face
B_T = 1.00                     # tip level = the linearity cap
TH_C = math.radians(30.0)      # wall tilt at the face (taper half-angle)
TH_T = math.radians(4.0)       # residual tilt at the tip
PHI2 = 6.0e-5                  # Phi/2 [Wb/m]: face half-width ~ 0.5 mm
NS = 161
EPS = 1e-9
EPSB = 4e-4

GAP_HALF = 20e-6               # sensor half-gap (mid-plane is the mirror)
BOX_X = 12e-3
BOX_Y = 10e-3
CAP = B_T


def _sramp(u):
    return np.sin(0.5 * math.pi * np.clip(u, 0.0, 1.0))


def wall_curve(t):
    """C1 wall profile in (B, theta): B climbs face->tip, tilt unwinds."""
    t = np.asarray(t, dtype=float)
    Bv = B_F + (B_T - B_F) * _sramp(t)
    th = TH_T + (TH_C - TH_T) * np.cos(0.5 * math.pi * np.clip(t, 0, 1))
    return Bv, th


# ---------------- design solve ----------------
def design_horn(report, maxh=0.02):
    geo = SplineGeometry()
    tw = np.linspace(0.0, 1.0, 61)
    wB, wT = wall_curve(tw)
    # counterclockwise circuit (a clockwise loop makes leftdomain=1 the
    # unbounded exterior and GenerateMesh grinds forever)
    pieces = [([(B_F, 0.0), (B_T, 0.0)], "center"),
              ([(B_T, 0.0), (B_T, TH_T)], "tip"),
              ([(float(b), float(t)) for b, t in
                zip(wB[::-1], wT[::-1])], "wall"),
              ([(B_F, TH_C), (B_F, 0.0)], "face")]
    allp, tags = [], []
    for pts_, bc in pieces:
        for i in range(len(pts_) - 1):
            allp.append(pts_[i])
            tags.append(bc)
    ids = [geo.AppendPoint(*p) for p in allp]
    for i, bc in enumerate(tags):
        geo.Append(["line", ids[i], ids[(i + 1) % len(ids)]], bc=bc,
                   leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    log(f"design: hodograph mesh {mesh.ne} elements")

    aB = x * mu_d_of(x) / mu_s_of(x) ** 2
    bB = 1.0 / (mu_s_of(x) * x)
    qB = x / mu_s_of(x)

    fes = H1(mesh, order=3, dirichlet="face|wall|tip|center")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (aB * grad(u)[0] * grad(v)[0] + bB * grad(u)[1] * grad(v)[1]) * dx
    a.Assemble()
    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"center": 0.0, "wall": PHI2,
                             "face": PHI2 * y / TH_C,
                             "tip": PHI2 * y / TH_T}, default=0.0), BND)
    res = gfA.vec.CreateVector()
    res.data = -a.mat * gfA.vec
    gfA.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                  inverse="sparsecholesky") * res

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
    log("design: coordinates recovered")

    def sample(cs):
        return np.array([(gfx(mesh(float(b), float(t))),
                          gfy(mesh(float(b), float(t)))) for b, t in cs])

    tw2 = np.linspace(0.0, 1.0, NS)
    wB2, wT2 = wall_curve(tw2)
    # inset the wall samples INTO the domain (smaller B and smaller theta)
    wall = sample([(min(max(b - EPSB, B_F + EPSB), B_T - EPSB),
                    max(t - 5e-3, EPS)) for b, t in zip(wB2, wT2)])
    cent = sample([(b, EPS) for b in np.linspace(B_F + EPSB, B_T - EPSB, NS)])
    # face/tip theta insets MATCH the wall's 5e-3 inset so the assembled
    # outline's corner points coincide (mismatched insets zigzag and cross)
    face = sample([(B_F + EPSB, t) for t in
                   np.linspace(EPS, TH_C - 5e-3, 81)])
    tipf = sample([(B_T - EPSB, t) for t in
                   np.linspace(EPS, TH_T - 5e-3, 41)])

    Jv = np.array([Jcf(mesh(float(b), max(float(t) - 5e-3, EPS)))
                   for b, t in zip(wB2, wT2)]
                  + [Jcf(mesh(float(b), EPS + 1e-4))
                     for b in np.linspace(B_F + EPSB, B_T - EPSB, 60)])
    single = bool(np.all(Jv > 0) or np.all(Jv < 0))

    length = float(np.linalg.norm(cent[-1] - cent[0]))
    w_face = 2.0 * float(np.linalg.norm(face[-1] - face[0]))
    w_tip = 2.0 * float(np.linalg.norm(tipf[-1] - tipf[0]))
    d = {"J_single_sign": single, "min_absJ": float(np.min(np.abs(Jv))),
         "length_mm": 1e3 * length, "face_width_mm": 1e3 * w_face,
         "tip_width_mm": 1e3 * w_tip,
         "geometric_gain": w_face / w_tip}
    print(f"  [design] J single sign={single}; horn length "
          f"{1e3*length:.3f} mm, face {1e3*w_face:.3f} -> tip "
          f"{1e3*w_tip:.3f} mm (gain {w_face/w_tip:.2f})")
    report["design"] = d
    return {"wall": wall, "cent": cent, "face": face, "tip": tipf}


# ---------------- FEM verification ----------------
def assemble_outline(curves):
    """Full horn outline from the half-design: mirror across the centerline
    (theta = 0 flux line: exactly straight)."""
    y_c = float(curves["cent"][:, 1].mean())

    def mir(c):
        out = c.copy()
        out[:, 1] = 2 * y_c - out[:, 1]
        return out

    face, wall, tip = curves["face"], curves["wall"], curves["tip"]
    outline = np.vstack([
        face,                       # centre of face -> upper corner
        wall[1:],                   # upper corner -> tip upper corner
        tip[::-1][1:],              # tip upper corner -> tip centre
        mir(tip)[1:],               # -> tip lower corner
        mir(wall[::-1])[1:],        # -> face lower corner
        mir(face)[::-1][1:-1],      # -> back toward the face centre
    ])
    keep = [0]
    for k in range(1, len(outline)):
        if np.linalg.norm(outline[k] - outline[keep[-1]]) > 2e-6:
            keep.append(k)
    outline = outline[keep]
    span = np.linalg.norm(outline.max(axis=0) - outline.min(axis=0))
    for k in range(1, len(outline)):
        if np.linalg.norm(outline[k] - outline[k - 1]) > 0.15 * span:
            raise RuntimeError(f"horn outline seam gap at {k}")
    return outline, y_c


def straight_outline(curves):
    """Straight taper: same face corners, tip corners, flat faces."""
    y_c = float(curves["cent"][:, 1].mean())
    fc = curves["face"][-1]
    tc = curves["wall"][-1]
    xf = float(curves["face"][:, 0].mean())
    xt = float(curves["tip"][:, 0].mean())
    fu = np.array([xf, fc[1]])
    tu = np.array([xt, tc[1]])
    li = np.linspace(0, 1, 60)[:, None]
    face_u = np.array([xf, y_c]) + li * (fu - np.array([xf, y_c]))
    wall_u = fu + li * (tu - fu)
    tip_u = tu + li * (np.array([xt, y_c]) - tu)

    def mir(c):
        out = c.copy()
        out[:, 1] = 2 * y_c - out[:, 1]
        return out

    outline = np.vstack([face_u, wall_u[1:], tip_u[1:],
                         mir(tip_u)[::-1][1:], mir(wall_u)[::-1][1:],
                         mir(face_u)[::-1][1:-1]])
    return outline, y_c


def build_box_with_iron(outline, y_c):
    """Place the horn with its tip at x = -GAP_HALF (mid-plane at x = 0),
    axis on y = 0; box with a graded mid-plane edge."""
    xs = outline[:, 0]
    left_w = np.ptp(outline[xs < xs.mean()][:, 1])
    right_w = np.ptp(outline[xs >= xs.mean()][:, 1])
    if left_w < right_w:            # tip on the left: flip x
        outline = outline.copy()
        outline[:, 0] = -outline[:, 0]
    shift = np.array([-GAP_HALF - outline[:, 0].max(), -y_c])
    poly = outline + shift

    geo = SplineGeometry()
    yg = np.unique(np.concatenate([
        np.linspace(-BOX_Y, -0.4e-3, 12), np.linspace(-0.4e-3, 0.4e-3, 33),
        np.linspace(0.4e-3, BOX_Y, 12)]))
    pts_r = [(0.0, float(v)) for v in yg]
    bl = geo.AppendPoint(-BOX_X, -BOX_Y)
    br = geo.AppendPoint(0.0, -BOX_Y)
    rid = [geo.AppendPoint(*p) for p in pts_r[1:-1]]
    tr = geo.AppendPoint(0.0, BOX_Y)
    tl = geo.AppendPoint(-BOX_X, BOX_Y)
    geo.Append(["line", bl, br], bc="amb", leftdomain=1, rightdomain=0)
    chain = [br] + rid + [tr]
    for i in range(len(chain) - 1):
        geo.Append(["line", chain[i], chain[i + 1]], bc="mid",
                   leftdomain=1, rightdomain=0)
    geo.Append(["line", tr, tl], bc="amb", leftdomain=1, rightdomain=0)
    geo.Append(["line", tl, bl], bc="amb", leftdomain=1, rightdomain=0)
    if polygon_area(poly[:, 0] * 1.0, poly[:, 1] * 1.0) < 0:
        poly = poly[::-1]
    hits = _self_intersections(poly)
    if hits:
        raise RuntimeError(f"iron outline self-intersects: {hits[:4]}")
    ids = [geo.AppendPoint(float(px), float(py)) for px, py in poly]
    for i in range(len(ids)):
        geo.Append(["line", ids[i], ids[(i + 1) % len(ids)]], bc="iron",
                   leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "iron")
    geo.SetDomainMaxH(1, 1.2e-3)
    geo.SetDomainMaxH(2, 0.03e-3)
    mesh = Mesh(geo.GenerateMesh(maxh=2e-3))
    return mesh, poly


def solve_ambient(mesh, B0, label, maxit=600, omega=0.4):
    fes = H1(mesh, order=2, dirichlet="amb")
    u, v = fes.TnT()
    gfA = GridFunction(fes)
    gfA.Set(B0 * y, BND, definedon=mesh.Boundaries("amb"))
    gfP, gfT = GridFunction(fes), GridFunction(fes)
    hist = []
    for it in range(maxit):
        gfP.vec.data = gfA.vec
        Bmag = sqrt(grad(gfA) * grad(gfA) + 1e-12)
        nu_cf = mesh.MaterialCF({"iron": nu_cf_of(Bmag)}, default=1.0 / MU0)
        a = BilinearForm(fes)
        a += nu_cf * grad(u) * grad(v) * dx
        a.Assemble()
        r = gfA.vec.CreateVector()
        r.data = -a.mat * gfA.vec
        gfT.vec.data = gfA.vec
        gfT.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * r
        gfA.vec.data = gfP.vec + omega * (gfT.vec - gfP.vec)
        d = gfP.vec.CreateVector()
        d.data = gfA.vec - gfP.vec
        rel = d.Norm() / max(gfA.vec.Norm(), 1e-30)
        hist.append(rel)
        if rel < 1e-9:
            break
        if len(hist) > 8 and rel > 0.98 * hist[-8] and omega > 0.08:
            omega *= 0.5
    if hist[-1] > 1e-7:
        raise RuntimeError(f"{label}: Picard stalled at {hist[-1]:.2e}")
    log(f"{label}: {len(hist)} its")
    return gfA


def measure(mesh, gfA, poly, B0):
    from matplotlib.path import Path
    Bx = grad(gfA)[1]                       # B_x = dA/dy
    Bmag = sqrt(grad(gfA) * grad(gfA))
    b_gap = float(Bx(mesh(-0.5 * GAP_HALF, 0.0)))
    xs = np.linspace(poly[:, 0].min(), poly[:, 0].max(), 160)
    ys = np.linspace(poly[:, 1].min(), poly[:, 1].max(), 90)
    pth = Path(poly)
    pk = 0.0
    sampled = 0
    failures = []
    ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
    for px in xs:
        for py in ys:
            if pth.contains_point((px, py), radius=-1e-5):
                try:
                    if ind(mesh(float(px), float(py))) > 0.5:
                        pk = max(pk, float(Bmag(mesh(float(px), float(py)))))
                        sampled += 1
                except Exception as exc:        # noqa: BLE001
                    failures.append((float(px), float(py), str(exc)))
    if failures:
        px, py, msg = failures[0]
        raise RuntimeError(
            f"iron peak sampling failed at ({px:.6e}, {py:.6e}): {msg} "
            f"({len(failures)} failed point(s))")
    if sampled == 0 or not math.isfinite(pk) or pk <= 0.0:
        raise RuntimeError(
            f"iron peak sampling produced no valid samples at B0={B0:.6e}")
    if not math.isfinite(b_gap):
        raise RuntimeError(f"gap-field sampling is not finite at B0={B0:.6e}")
    return b_gap, pk


# ---------------- golden bands (locked 2026-07-29) ----------------
GOLD = {
    "geometric_gain": (8.0, 9.6),           # measured 8.79
    "length_mm": (0.80, 1.05),              # measured 0.916
    "B0_rated_mT": (55.0, 100.0),           # measured 76.1
    "gain_low_horn": (10.0, 11.4),          # measured 10.72
    "gain_low_straight": (9.5, 10.9),       # measured 10.20
    "gain_advantage": (0.02, 0.09),         # measured +5.1 %
    "peak_rated_horn": (0.90, 1.03),        # measured 0.964
    "peak_margin_min": 0.005,               # straight-horn, measured 0.030
    "ladder_max": 0.01,                     # measured 0.0033
    "gain3x_deficit_max": 0.05,             # horn 9.50 vs straight 9.40
}


def _band(name, val, lo, hi):
    if not (lo <= val <= hi):
        raise AssertionError(f"golden {name}: {val:.4e} outside "
                             f"[{lo:.3e}, {hi:.3e}]")
    print(f"  [golden] {name} = {val:.4e} in [{lo:.3e}, {hi:.3e}]  OK")


def main():
    global BOX_X, BOX_Y
    SetNumThreads(4)
    here = os.path.dirname(os.path.abspath(__file__))
    report = {
        "schema": "radia.validation.clebsch-concentrator-horn.v1",
        "case": {
            "gap_um": 2e6 * GAP_HALF,
            "box_mm": [1e3 * BOX_X, 1e3 * BOX_Y],
            "cap_T": CAP,
            "material": f"permalloy-like MUR0={MUR0:.0f} BK={BK} NEXP={NEXP}",
            "B_face_T": B_F,
            "B_tip_cap_T": B_T,
            "taper_deg_face_to_tip": [math.degrees(TH_C), math.degrees(TH_T)],
            "phi_half_Wb_per_m": PHI2,
        },
    }
    with TaskManager():
        print("step (1) hodograph design of the horn")
        curves = design_horn(report)
        if not report["design"]["J_single_sign"]:
            raise AssertionError("design Jacobian changes sign (fold)")
        _band("geometric_gain", report["design"]["geometric_gain"],
              *GOLD["geometric_gain"])
        _band("length_mm", report["design"]["length_mm"], *GOLD["length_mm"])

        out_h, yc_h = assemble_outline(curves)
        out_s, yc_s = straight_outline(curves)
        mesh_h, poly_h = build_box_with_iron(out_h, yc_h)
        mesh_s, poly_s = build_box_with_iron(out_s, yc_s)
        area_h = abs(polygon_area(poly_h[:, 0], poly_h[:, 1]))
        area_s = abs(polygon_area(poly_s[:, 0], poly_s[:, 1]))
        print(f"  meshes: horn {mesh_h.ne} el, straight {mesh_s.ne} el; "
              f"iron areas {1e6*area_h:.4f} / {1e6*area_s:.4f} mm^2")
        report["iron_area_mm2"] = {"horn": 1e6 * area_h,
                                   "straight": 1e6 * area_s}

        print("step (2) rated ambient search (horn iron peak -> cap)")
        B_lo = B_hi = None
        B0 = 2.0e-3
        for _ in range(14):
            g = solve_ambient(mesh_h, B0, f"horn_B{B0:.2e}")
            bg, pk = measure(mesh_h, g, poly_h, B0)
            print(f"  [rated?] B0={1e3*B0:.3f} mT: iron peak {pk:.3f} T, "
                  f"gap {1e3*bg:.2f} mT", flush=True)
            if pk < 0.96 * CAP:
                B_lo = B0
                B0 = B0 * 2.0 if B_hi is None else math.sqrt(B0 * B_hi)
            elif pk > 1.02 * CAP:
                B_hi = B0
                B0 = B0 / 2.0 if B_lo is None else math.sqrt(B0 * B_lo)
            else:
                break
        else:
            raise RuntimeError("rated ambient search did not converge")
        B0_rated = B0
        report["B0_rated_mT"] = 1e3 * B0_rated
        _band("B0_rated_mT", 1e3 * B0_rated, *GOLD["B0_rated_mT"])

        print("step (3) ambient sweep: gain + linearity + peak, both shapes")
        sweep = [0.1, 0.3, 0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.6, 2.0,
                 2.5, 3.0]
        results = {}
        for name, mesh_x, poly_x in (("horn", mesh_h, poly_h),
                                     ("straight", mesh_s, poly_s)):
            rows = []
            for m in sweep:
                B0m = m * B0_rated
                g = solve_ambient(mesh_x, B0m, f"{name}_x{m:g}")
                bg, pk = measure(mesh_x, g, poly_x, B0m)
                rows.append({"mult": m, "B0_mT": 1e3 * B0m,
                             "B_gap_mT": 1e3 * bg, "gain": bg / B0m,
                             "iron_peak_T": pk})
                print(f"  [{name} x{m:g}] gain {bg/B0m:7.3f}  iron peak "
                      f"{pk:.3f} T", flush=True)
            g0 = rows[0]["gain"]
            for r in rows:
                r["linearity_pct"] = 100.0 * (r["gain"] / g0 - 1.0)
            results[name] = rows
        report["sweep"] = results

        print("step (4) truncation ladder: rated case at a 1.6x box")
        bx, by = BOX_X, BOX_Y
        BOX_X, BOX_Y = 1.6 * bx, 1.6 * by
        mesh_L, poly_L = build_box_with_iron(out_h, yc_h)
        gL = solve_ambient(mesh_L, B0_rated, "horn_bigbox")
        bgL, pkL = measure(mesh_L, gL, poly_L, B0_rated)
        BOX_X, BOX_Y = bx, by
        g_small = [r for r in results["horn"] if r["mult"] == 1.0][0]
        lad = abs(bgL * 1e3 - g_small["B_gap_mT"]) / g_small["B_gap_mT"]
        report["ladder"] = {"B_gap_mT_box1": g_small["B_gap_mT"],
                            "B_gap_mT_box1p6": 1e3 * bgL,
                            "rel_shift": lad}
        print(f"  [ladder] gap field shifts {100*lad:.2f}% at a 1.6x box")

    hr = {r["mult"]: r for r in results["horn"]}
    sr = {r["mult"]: r for r in results["straight"]}
    report["verdict"] = {
        "gain_low": {"horn": hr[0.1]["gain"], "straight": sr[0.1]["gain"]},
        "gain_advantage": hr[0.1]["gain"] / sr[0.1]["gain"] - 1.0,
        "iron_peak_at_rated_T": {"horn": hr[1.0]["iron_peak_T"],
                                 "straight": sr[1.0]["iron_peak_T"]},
        "gain_at_3x": {"horn": hr[3.0]["gain"],
                       "straight": sr[3.0]["gain"]},
    }
    print(f"\n== gain(low): horn {hr[0.1]['gain']:.3f} vs straight "
          f"{sr[0.1]['gain']:.3f} ({100*report['verdict']['gain_advantage']:+.1f}%)")
    print(f"== iron peak at rated: horn {hr[1.0]['iron_peak_T']:.3f} T vs "
          f"straight {sr[1.0]['iron_peak_T']:.3f} T (cap {CAP})")

    _band("gain_low_horn", hr[0.1]["gain"], *GOLD["gain_low_horn"])
    _band("gain_low_straight", sr[0.1]["gain"], *GOLD["gain_low_straight"])
    _band("gain_advantage", report["verdict"]["gain_advantage"],
          *GOLD["gain_advantage"])
    _band("peak_rated_horn", hr[1.0]["iron_peak_T"],
          *GOLD["peak_rated_horn"])
    _band("peak_margin", sr[1.0]["iron_peak_T"] - hr[1.0]["iron_peak_T"],
          GOLD["peak_margin_min"], 1.0)
    _band("ladder", lad, 0.0, GOLD["ladder_max"])
    _band("gain3x_horn_vs_straight",
          sr[3.0]["gain"] - hr[3.0]["gain"], -10.0,
          GOLD["gain3x_deficit_max"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    y_c = float(curves["cent"][:, 1].mean())
    for k, col, lb in (("wall", "tab:red", "wall (A=Phi/2)"),
                       ("cent", "0.4", "centerline"),
                       ("face", "tab:green", "collection face"),
                       ("tip", "tab:blue", "tip face")):
        cv = curves[k]
        ax[0].plot(1e3 * cv[:, 0], 1e3 * cv[:, 1], color=col, lw=1.4,
                   label=lb)
        ax[0].plot(1e3 * cv[:, 0], 1e3 * (2 * y_c - cv[:, 1]), color=col,
                   lw=0.7, alpha=0.5)
    ax[0].set_aspect("equal")
    ax[0].legend(fontsize=7)
    ax[0].set_xlabel("x [mm]")
    ax[0].set_ylabel("y [mm]")
    for name, col in (("horn", "tab:red"), ("straight", "tab:blue")):
        rows = results[name]
        ax[1].plot([r["B0_mT"] for r in rows], [r["gain"] for r in rows],
                   "o-", color=col, label=name, ms=3)
        ax[2].plot([r["B0_mT"] for r in rows],
                   [r["iron_peak_T"] for r in rows], "o-", color=col,
                   label=name, ms=3)
    ax[1].set_xlabel("ambient B0 [mT]")
    ax[1].set_ylabel("gain")
    ax[1].legend(fontsize=8)
    ax[2].axhline(CAP, color="0.4", ls="--", lw=0.8)
    ax[2].set_xlabel("ambient B0 [mT]")
    ax[2].set_ylabel("iron peak |B| [T]")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    png = os.path.join(here, "concentrator_horn_verify.png")
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
    out = os.path.join(here, "results_concentrator_horn.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"results -> {out}")
    print("ALL GOLDEN CHECKS PASSED")


if __name__ == "__main__":
    main()
