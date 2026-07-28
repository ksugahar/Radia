"""IPM rotor bridge: free-boundary hodograph design -> shape -> nonlinear FEM check.

Engineering question
--------------------
Flux leaks around the tip of an IPM flux barrier through a thin iron bridge.
The bridge is a FLUX FUNNEL that turns: flux is collected on the wide side and
concentrated into the narrow throat.  The hard constraint is a CAP:

    |B| must nowhere exceed the knee B_knee

because past the knee the iron stops limiting the leakage predictably and the
local loss density explodes.  The cheapest iron is the iron that sits EXACTLY
at the cap on its most-loaded surface -- no hot spot, no wasted margin.

In physical space that surface is a FREE BOUNDARY (an unknown curve on which a
field condition holds), and finding it needs a nonlinear shape iteration.  In
the field-plane (Chaplygin) hodograph the same condition is a COORDINATE LINE

    B = B_knee

so the free boundary becomes a fixed Dirichlet edge and the whole design is
ONE LINEAR SOLVE.  That is the claim under test.

Design spec (hodograph domain)
    inner (barrier-side) wall  B = B_knee              CONSTANT   <- the cap
    outer (rotor-OD-side) wall B = B_out(theta)        ramps up   <- funnel closes
    flux per unit stack length dA = Phi_leak
    lead-in / lead-out of MARGIN degrees at each end hold B_out flat, so the
    designed BODY (0 <= theta <= Theta) is not contaminated by the terminals.
    A real bridge merges into the core the same way; an abrupt terminal is a
    modelling artefact, not a design feature.

Scale-freedom (used, not assumed): the hodograph equation is LINEAR in A and
Psi is linear in A, so the recovered geometry scales EXACTLY linearly with
Phi_leak.  One solve gives the whole flux family; the run checks this.

Forward verification
    mesh the designed outline, solve  div( nu(|grad A|) grad A ) = 0  with the
    same Dirichlet flux values, and compare |B| on the walls with the spec.

Naive baseline
    what an engineer draws instead: circular centreline through the same two
    body end mid-points, width tapered linearly between the same widths, the
    same lead-in/lead-out, the same flux.  Its inner-wall peak |B| is the
    number that matters -- if it overshoots B_knee the cap is violated.

Formulation (B-radial A-form; identical to the verified 90-degree bend case in
validation_test/clebsch_legendre/verify_chaplygin_bend_design.py)
    d/dB( a A_B ) + b A_thth = 0,   a = B mu_d / mu_s^2,   b = 1/(mu_s B)
    Psi_B = -b A_th,                Psi_th = a A_B
    dr    = (dPsi/q) e_H + (dA/B) e_perp,   q = B/mu_s = |H|

Golden bands asserted at the end of the run (2026-07-28 baseline, LAB):
  constant-mu sanity    : both designed walls are an exact annulus (dev < 1e-6)
  flux scale-freedom    : halving Phi halves every wall coordinate (dev < 1e-6)
  orientation           : J keeps one sign on every design (no folding)
  inner wall vs the cap : mean < 1.0 %, max < 2.0 % over the designed body,
                          both mesh resolutions
  outer wall vs the ramp: mean < 1.5 %
  MMF                   : |design - FEM| / design < 1.0 %
  naive baseline        : must overshoot the cap by > 1 %, else the comparison
                          has no content as posed

Measured at that baseline: designed body inner wall 1.875..1.902 T against a
1.900 T cap (overshoot +0.11 %, mean err 0.43 %), outer wall 0.908..1.694 T
against the 0.900..1.700 T ramp (mean err 0.13 %), MMF 2.4447 A design vs
2.4476 A FEM (0.12 %), h/8 and h/16 agreeing to three digits.  The naive
baseline spans 1.326..2.157 T on the same wall: cap overshot by 13.5 %, spread
43.7 % against 1.4 %, and 6.0 % more iron for the same flux.

The lead-in / lead-out is NOT cosmetic.  With the terminals attached directly
to the body the inlet corner contaminates the first ~11 degrees of the inner
wall with up to 9.1 % error, and it is mesh-INDEPENDENT (9.07 % at h/8, 9.00 %
at h/16), i.e. a real terminal effect and not discretisation.  The near-
degenerate Jacobian is localised there too (min |J| 2.0e-8 against a 1.2e-6
median, at the terminal corner, outside the body).

Run:  python verify_ipm_bridge_free_boundary.py
Writes results_ipm_bridge_free_boundary.json next to this file (committed).
"""
import datetime
import json
import math
import os
import platform
import time

import numpy as np
from netgen.geom2d import SplineGeometry
from ngsolve import (
    BND, BitArray, BilinearForm, CoefficientFunction, GridFunction, H1,
    LinearForm, Mesh, TaskManager, cos, dx, grad, sin, sqrt, x, y, SetNumThreads,
)

_T0 = time.time()


def log(msg):
    print(f"[{time.time() - _T0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- material
# Representative non-oriented silicon steel.  NOT a datasheet fit: what is
# under test is that the hodograph handles a first-order saturating curve, so
# the model only has to be qualitatively right.  (B, H, mu_r, mu_d) samples go
# in the JSON so the curve actually used is auditable.
MU0 = 4.0e-7 * math.pi
MUR0 = 7000.0
BK = 1.0
NEXP = 4


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


# ---------------------------------------------------------------- design case
THETA = 0.5 * math.pi           # designed body: flux turns 90 deg
MARGIN = math.radians(20.0)     # lead-in / lead-out, flat B_out
B_KNEE = 1.90                   # the cap: barrier-side surface sits here
B_OUT0 = 0.90                   # outer wall at the wide inlet
B_OUT1 = 1.70                   # outer wall at the narrow throat
DFLUX = 1.1e-3                  # Wb per metre of stack
NSAMP = 175                     # wall samples over the FULL span
TRIM = 6                        # samples dropped at each terminal corner


def b_out_of(th, b1=B_OUT1):
    t = np.clip(np.asarray(th, dtype=float) / THETA, 0.0, 1.0)
    return B_OUT0 + (b1 - B_OUT0) * t


# ---------------------------------------------------------------- helpers
def recover_potential(mesh, F, order=3):
    """Return u with grad(u) = F.  Gauge: pin one DOF (a NumberSpace row would
    couple every DOF and densify the direct factorisation)."""
    fes = H1(mesh, order=order)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    f = LinearForm(fes)
    f += (F * grad(v)) * dx
    a.Assemble()
    f.Assemble()
    free = BitArray(fes.FreeDofs())
    free[0] = False
    gf = GridFunction(fes)
    gf.vec[:] = 0.0
    gf.vec.data += a.mat.Inverse(free, inverse="sparsecholesky") * f.vec
    return gf


def polygon_area(px, py):
    return 0.5 * float(np.sum(px * np.roll(py, -1) - np.roll(px, -1) * py))


def band_area(inner, outer):
    """Signed area of the strip bounded by the two sampled wall curves."""
    loop = np.vstack([outer, inner[::-1]])
    return abs(polygon_area(loop[:, 0], loop[:, 1]))


def _self_intersections(loop):
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    n = len(loop)
    hits = []
    for i in range(n):
        p, p2 = loop[i], loop[(i + 1) % n]
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            q, q2 = loop[j], loop[(j + 1) % n]
            d1, d2 = cross(q, q2, p), cross(q, q2, p2)
            d3, d4 = cross(p, p2, q), cross(p, p2, q2)
            if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
                hits.append((i, j))
    return hits


# ---------------------------------------------------------------- step (1)
def design(report, dflux=DFLUX, mu_r_const=None, ramp=True, tag=None,
           maxh=0.025):
    """One linear hodograph solve -> physical wall curves (full span)."""
    if tag is None:
        tag = ("linear" if mu_r_const else "saturable") + ("" if ramp else "_flat")
    b1 = B_OUT1 if ramp else B_OUT0
    t0, t1 = -MARGIN, THETA + MARGIN

    def mus(B):
        return MU0 * mu_r_const if mu_r_const else mu_s_of(B)

    def mud(B):
        return MU0 * mu_r_const if mu_r_const else mu_d_of(B)

    geo = SplineGeometry()
    corners = [(B_OUT0, t0), (B_KNEE, t0), (B_KNEE, t1), (b1, t1),
               (b1, THETA), (B_OUT0, 0.0)]
    edges = ["inlet", "inner", "outlet", "outer", "outer", "outer"]
    if not ramp:                       # flat: the two outer kinks are collinear
        corners = [(B_OUT0, t0), (B_KNEE, t0), (B_KNEE, t1), (B_OUT0, t1)]
        edges = ["inlet", "inner", "outlet", "outer"]
    p = [geo.AppendPoint(*c) for c in corners]
    for i, bc in enumerate(edges):
        geo.Append(["line", p[i], p[(i + 1) % len(p)]], bc=bc,
                   leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    log(f"design/{tag}: hodograph mesh {mesh.ne} elements")

    aB = x * mud(x) / mus(x) ** 2
    bB = 1.0 / (mus(x) * x)
    qB = x / mus(x)

    fes = H1(mesh, order=3, dirichlet="inner|outer")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (aB * grad(u)[0] * grad(v)[0] + bB * grad(u)[1] * grad(v)[1]) * dx
    a.Assemble()
    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"outer": 0.0, "inner": dflux}, default=0.0), BND)
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
    log(f"design/{tag}: coordinates recovered")

    ths = np.linspace(t0, t1, NSAMP)
    eps = 1e-9
    inner, outer, Jv = [], [], []
    for th in ths:
        tc = min(max(th, t0 + eps), t1 - eps)
        pi_ = mesh(B_KNEE - eps, tc)
        po_ = mesh(float(b_out_of(th, b1)) + eps, tc)
        inner.append((gfx(pi_), gfy(pi_)))
        outer.append((gfx(po_), gfy(po_)))
        Jv.append((Jcf(pi_), B_KNEE, th))
        Jv.append((Jcf(po_), float(b_out_of(th, b1)), th))
    inner = np.array(inner)
    outer = np.array(outer)

    bs0 = np.linspace(B_OUT0 + eps, B_KNEE - eps, 60)
    bs1 = np.linspace(b1 + eps, B_KNEE - eps, 60)
    inlet = np.array([(gfx(mesh(b, t0 + eps)), gfy(mesh(b, t0 + eps))) for b in bs0])
    outlet = np.array([(gfx(mesh(b, t1 - eps)), gfy(mesh(b, t1 - eps))) for b in bs1])

    mmf = float(gfPsi(mesh(0.5 * (b1 + B_KNEE), t1 - eps))
                - gfPsi(mesh(0.5 * (B_OUT0 + B_KNEE), t0 + eps)))

    body = (ths >= -1e-12) & (ths <= THETA + 1e-12)
    Jarr = np.array([j for j, _, _ in Jv])
    kmin = int(np.argmin(np.abs(Jarr)))
    w0 = float(np.linalg.norm(inner[body][0] - outer[body][0]))
    w1 = float(np.linalg.norm(inner[body][-1] - outer[body][-1]))
    out = {
        "tag": tag, "dflux_Wb_per_m": dflux,
        "margin_deg": math.degrees(MARGIN),
        "J_single_sign": bool(np.all(Jarr < 0) or np.all(Jarr > 0)),
        "min_absJ": float(np.min(np.abs(Jarr))),
        "min_absJ_at_B_theta_deg": [Jv[kmin][1], math.degrees(Jv[kmin][2])],
        "median_absJ": float(np.median(np.abs(Jarr))),
        "mmf_design_A": mmf,
        "body_inlet_width_mm": 1e3 * w0, "body_throat_width_mm": 1e3 * w1,
        "body_iron_area_mm2": 1e6 * band_area(inner[body], outer[body]),
        "total_length_mm": 1e3 * float(np.sum(np.linalg.norm(
            np.diff(0.5 * (inner + outer), axis=0), axis=1))),
    }
    print(f"  [design/{tag}] J single sign={out['J_single_sign']}  "
          f"min|J|={out['min_absJ']:.3e} at B={Jv[kmin][1]:.2f} T, "
          f"theta={math.degrees(Jv[kmin][2]):.1f} deg (median {out['median_absJ']:.3e})")
    print(f"  [design/{tag}] body inlet={1e3*w0:.3f} mm  throat={1e3*w1:.3f} mm"
          f"  body iron={out['body_iron_area_mm2']:.4f} mm^2  MMF={mmf:.3f} A")
    report.setdefault("design", {})[tag] = out
    return inner, outer, inlet, outlet, ths, body


def annulus_sanity(inner, outer, report):
    """Constant-mu + both walls at constant B must come out an exact annulus."""
    def fit_circle(pts):
        A = np.c_[2 * pts[:, 0], 2 * pts[:, 1], np.ones(len(pts))]
        bb = (pts ** 2).sum(axis=1)
        sol, *_ = np.linalg.lstsq(A, bb, rcond=None)
        cx, cy = sol[0], sol[1]
        r = math.sqrt(sol[2] + cx * cx + cy * cy)
        dev = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - r).max() / r
        return (cx, cy), r, dev
    ci, ri, di = fit_circle(inner)
    co, ro, do = fit_circle(outer)
    res = {
        "inner_radius_mm": 1e3 * ri, "outer_radius_mm": 1e3 * ro,
        "inner_circularity_dev": di, "outer_circularity_dev": do,
        "radius_ratio": ro / ri,
        "expected_ratio_Bknee_over_Bout": B_KNEE / B_OUT0,
        "ratio_rel_err": abs(ro / ri - B_KNEE / B_OUT0) / (B_KNEE / B_OUT0),
        "center_offset_mm": 1e3 * float(np.hypot(ci[0] - co[0], ci[1] - co[1])),
    }
    print(f"  [sanity] constant-mu -> annulus r_in={1e3*ri:.3f} r_out={1e3*ro:.3f} mm"
          f"  circ dev {di:.2e}/{do:.2e}  ratio err {res['ratio_rel_err']:.2e}")
    report["sanity_constant_mu_annulus"] = res


# ---------------------------------------------------------------- naive shape
def naive_funnel(inner, outer, ths, body):
    """What an engineer draws instead: a circular centreline through the same
    two BODY end mid-points, turning by the same angle, width tapered linearly
    between the same two widths, with the same flat lead-in / lead-out."""
    ib = np.flatnonzero(body)
    i0, i1 = ib[0], ib[-1]
    m0 = 0.5 * (inner[i0] + outer[i0])
    m1 = 0.5 * (inner[i1] + outer[i1])
    w0 = float(np.linalg.norm(inner[i0] - outer[i0]))
    w1 = float(np.linalg.norm(inner[i1] - outer[i1]))

    chord = m1 - m0
    cl = float(np.linalg.norm(chord))
    R = cl / (2.0 * math.sin(0.5 * THETA))
    hh = R * math.cos(0.5 * THETA)
    mid = 0.5 * (m0 + m1)
    nrm = np.array([-chord[1], chord[0]]) / cl
    # the arc bulges AWAY from its centre, so the designed centreline's
    # mid-point says which side of the chord the centre is NOT on
    k = (i0 + i1) // 2
    side = np.sign(np.dot(0.5 * (inner[k] + outer[k]) - mid, nrm)) or 1.0
    ctr = mid - side * hh * nrm

    a0 = math.atan2(m0[1] - ctr[1], m0[0] - ctr[0])
    a1 = math.atan2(m1[1] - ctr[1], m1[0] - ctr[0])
    while a1 - a0 > math.pi:
        a1 -= 2 * math.pi
    while a1 - a0 < -math.pi:
        a1 += 2 * math.pi
    # map the design's theta samples onto the arc, margins included
    sgn = 1.0 if a1 > a0 else -1.0
    angs = a0 + sgn * (ths - ths[np.flatnonzero(body)[0]])
    half = 0.5 * (w0 + (w1 - w0) * np.clip(np.asarray(ths) / THETA, 0.0, 1.0))

    ur = np.c_[np.cos(angs), np.sin(angs)]
    inn = ctr + (R - half)[:, None] * ur      # inside of the turn = high |B|
    out = ctr + (R + half)[:, None] * ur
    ss = np.linspace(0.0, 1.0, 60)[:, None]   # end faces: outer -> inner
    inlet = out[0] + (inn[0] - out[0]) * ss
    outlet = out[-1] + (inn[-1] - out[-1]) * ss
    return inn, out, inlet, outlet, {
        "R_mm": 1e3 * R, "centre_m": ctr.tolist(),
        "body_turn_deg": abs(math.degrees(a1 - a0)),
        "body_inlet_width_mm": 1e3 * w0, "body_throat_width_mm": 1e3 * w1,
        "body_iron_area_mm2": 1e6 * band_area(inn[body], out[body]),
    }


# ---------------------------------------------------------------- step (2)
def forward_solve(inner, outer, inlet, outlet, maxh, label, dflux=DFLUX):
    """Mesh the outline and solve the nonlinear physical problem on it."""
    parts = [(outer, "wall_out"), (outlet[1:], "outlet"),
             (inner[::-1][1:], "wall_in"), (inlet[::-1][1:-1], "inlet")]
    loop = np.vstack([p for p, _ in parts])
    tags = [t for p, t in parts for _ in range(len(p))]
    # corners are sampled twice (once by each adjacent curve); the leftover
    # sub-micron segments stall the mesher, so dedup at a mesh-relative scale
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
    pids = [geo.AppendPoint(float(px), float(py)) for px, py in loop]
    for i in range(n):
        geo.Append(["line", pids[i], pids[(i + 1) % n]], bc=tags[i],
                   leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    log(f"{label}: mesh {mesh.ne} elements, {mesh.nv} vertices")

    fes = H1(mesh, order=3, dirichlet="wall_in|wall_out")
    u, v = fes.TnT()
    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"wall_out": 0.0, "wall_in": dflux}, default=0.0), BND)
    gfPrev, gfTrial = GridFunction(fes), GridFunction(fes)
    eps = 1e-6
    omega = 0.35                 # under-relaxation; plain Picard oscillates
    hist = []
    for _ in range(600):
        gfPrev.vec.data = gfA.vec
        a = BilinearForm(fes)
        a += nu_cf_of(sqrt(grad(gfA) * grad(gfA) + eps ** 2)) * grad(u) * grad(v) * dx
        a.Assemble()
        gfTrial.vec.data = gfA.vec
        res = gfA.vec.CreateVector()
        res.data = -a.mat * gfTrial.vec
        gfTrial.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                          inverse="sparsecholesky") * res
        gfA.vec.data = gfPrev.vec + omega * (gfTrial.vec - gfPrev.vec)
        d = gfPrev.vec.CreateVector()
        d.data = gfA.vec - gfPrev.vec
        hist.append(d.Norm() / max(gfA.vec.Norm(), 1e-30))
        if hist[-1] < 1e-9:
            break
    log(f"{label}: Picard {len(hist)} iterations, final rel step {hist[-1]:.2e}")
    if hist[-1] > 1e-7:
        raise RuntimeError(f"{label}: Picard did not converge "
                           f"(rel step {hist[-1]:.2e} after {len(hist)} its)")

    # Probe the walls a hair INSIDE the domain: the sample points sit exactly on
    # the (polygonal) boundary, where NGSolve's point search is not reliable.
    # 0.5 % of the local channel width biases |B| by <0.3 % of B_knee here, far
    # below the bands being checked.  The first/last TRIM samples are dropped
    # outright: a transverse inset cannot move a corner sample off the end face.
    span = np.linalg.norm(inner - outer, axis=1)[:, None]
    dirn = (outer - inner) / np.maximum(span, 1e-30)
    pin = (inner + 0.005 * span * dirn)[TRIM:-TRIM]
    pout = (outer - 0.005 * span * dirn)[TRIM:-TRIM]

    def probe(cf, pts, what):
        vals = []
        for px, py in pts:
            try:
                vals.append(cf(mesh(float(px), float(py))))
            except Exception as exc:                      # noqa: BLE001
                raise RuntimeError(
                    f"{label}: {what} probe ({px:.6e}, {py:.6e}) is not inside "
                    f"the meshed outline ({exc})") from exc
        return np.array(vals)

    Bcf = sqrt(grad(gfA) * grad(gfA))
    b_in_fem = probe(Bcf, pin, "inner wall")
    b_out_fem = probe(Bcf, pout, "outer wall")

    Bvec = CoefficientFunction((grad(gfA)[1], -grad(gfA)[0]))
    Hvec = nu_cf_of(sqrt(grad(gfA) * grad(gfA) + eps ** 2)) * Bvec
    mid = 0.5 * (inner + outer)
    mmf = 0.0
    for i in range(len(mid) - 1):
        seg = mid[i + 1] - mid[i]
        pm = 0.5 * (mid[i] + mid[i + 1])
        hv = Hvec(mesh(pm[0], pm[1]))
        mmf += hv[0] * seg[0] + hv[1] * seg[1]
    return {
        "n_elements": mesh.ne, "picard_iterations": len(hist),
        "wall_probe_inset_frac_of_width": 0.005, "trim_samples_each_end": TRIM,
        "inner_B_fem_T": b_in_fem.tolist(), "outer_B_fem_T": b_out_fem.tolist(),
        "mmf_fem_A": abs(float(mmf)),
    }, b_in_fem, b_out_fem


def summarise(bi, bo, thc, body_c, label, tag_out):
    e_in = np.abs(bi[body_c] - B_KNEE) / B_KNEE
    tgt_o = b_out_of(thc[body_c])
    e_out = np.abs(bo[body_c] - tgt_o) / tgt_o
    d = {
        "inner_wall_rel_err_mean": float(e_in.mean()),
        "inner_wall_rel_err_max": float(e_in.max()),
        "outer_wall_rel_err_mean": float(e_out.mean()),
        "outer_wall_rel_err_max": float(e_out.max()),
        "inner_peak_B_T": float(bi[body_c].max()),
        "inner_min_B_T": float(bi[body_c].min()),
        "cap_overshoot_pct": 100.0 * (float(bi[body_c].max()) - B_KNEE) / B_KNEE,
        "inner_spread_pct": 100.0 * (float(bi[body_c].max())
                                     - float(bi[body_c].min())) / B_KNEE,
    }
    print(f"  [{label}] body inner |B| = {d['inner_min_B_T']:.3f}..{d['inner_peak_B_T']:.3f} T"
          f"  (cap {B_KNEE:.2f} T, overshoot {d['cap_overshoot_pct']:+.2f}%,"
          f" spread {d['inner_spread_pct']:.2f}%)")
    print(f"  [{label}] rel err inner mean={e_in.mean()*100:.3f}% max={e_in.max()*100:.3f}%"
          f"   outer mean={e_out.mean()*100:.3f}% max={e_out.max()*100:.3f}%")
    return d


def main():
    SetNumThreads(4)
    bh = [0.5, 1.0, 1.43, 1.7, 1.9, 2.1]
    report = {"case": {
        "material": f"mu_r(B)=1+({MUR0}-1)/(1+(B/{BK})^{NEXP}) "
                    f"[representative, not a datasheet fit]",
        "material_samples": [
            {"B_T": b, "H_A_per_m": b / mu_s_of(b), "mu_r_secant": mu_r_of(b),
             "mu_r_differential": mu_d_of(b) / MU0,
             "anisotropy_mud_over_mus": mu_d_of(b) / mu_s_of(b)} for b in bh],
        "body_turn_deg": math.degrees(THETA), "margin_deg": math.degrees(MARGIN),
        "B_knee_cap_T": B_KNEE, "B_outer_T": [B_OUT0, B_OUT1],
        "flux_Wb_per_m": DFLUX,
    }}
    for s in report["case"]["material_samples"]:
        print(f"  [material] B={s['B_T']:.2f} T  H={s['H_A_per_m']:8.1f} A/m  "
              f"mu_r,s={s['mu_r_secant']:7.1f}  mu_r,d={s['mu_r_differential']:6.1f}"
              f"  mu_d/mu_s={s['anisotropy_mud_over_mus']:.4f}")

    with TaskManager():
        print("\nstep (0) machinery sanity: constant mu, both walls at constant B")
        li, lo, _, _, _, _ = design(report, mu_r_const=200.0, ramp=False)
        annulus_sanity(li, lo, report)

        print("\nstep (1) design: saturable funnel, inner wall pinned at the knee")
        inner, outer, inlet, outlet, ths, body = design(report)

        print("\nstep (1b) flux scale-freedom: same field spec at half the flux")
        hi, ho, _, _, _, _ = design(report, dflux=0.5 * DFLUX,
                                    tag="saturable_halfflux")
        sc = np.concatenate([
            np.linalg.norm(hi, axis=1) / np.maximum(np.linalg.norm(inner, axis=1), 1e-30),
            np.linalg.norm(ho, axis=1) / np.maximum(np.linalg.norm(outer, axis=1), 1e-30)])
        report["flux_scale_freedom"] = {
            "expected_ratio": 0.5, "measured_mean": float(sc.mean()),
            "measured_max_dev": float(np.abs(sc - 0.5).max())}
        print(f"  [scale] geometry ratio at half flux: mean={sc.mean():.8f} "
              f"(expect 0.5), max dev={np.abs(sc-0.5).max():.2e}")

        thc = ths[TRIM:-TRIM]
        body_c = body[TRIM:-TRIM]
        wt = float(np.linalg.norm(inner[body][-1] - outer[body][-1]))

        print("\nstep (2) forward nonlinear FEM on the designed shape")
        report["verify"] = {}
        for div in (8.0, 16.0):
            lab = f"design_h{div:g}"
            res, bi, bo = forward_solve(inner, outer, inlet, outlet, wt / div, lab)
            res.update(summarise(bi, bo, thc, body_c, lab, res))
            report["verify"][lab] = res

        print("\nstep (3) naive baseline: circular centreline + linear width taper")
        ni, no, nin, nout, ninfo = naive_funnel(inner, outer, ths, body)
        nres, nbi, nbo = forward_solve(ni, no, nin, nout, wt / 16.0, "naive_h16")
        nres["arc"] = ninfo
        nres.update(summarise(nbi, nbo, thc, body_c, "naive_h16", nres))
        print(f"  [naive] arc R={ninfo['R_mm']:.3f} mm, body turn="
              f"{ninfo['body_turn_deg']:.1f} deg, widths "
              f"{ninfo['body_inlet_width_mm']:.3f}->"
              f"{ninfo['body_throat_width_mm']:.3f} mm, body iron="
              f"{ninfo['body_iron_area_mm2']:.4f} mm^2")
        report["naive"] = nres

    d = report["design"]["saturable"]
    n = report["naive"]["arc"]
    report["comparison"] = {
        "body_iron_area_mm2": {"hodograph": d["body_iron_area_mm2"],
                               "naive": n["body_iron_area_mm2"],
                               "naive_excess_pct": 100.0 * (n["body_iron_area_mm2"]
                                                            - d["body_iron_area_mm2"])
                               / d["body_iron_area_mm2"]},
        "cap_overshoot_pct": {"hodograph": report["verify"]["design_h16"]["cap_overshoot_pct"],
                              "naive": report["naive"]["cap_overshoot_pct"]},
        "inner_spread_pct": {"hodograph": report["verify"]["design_h16"]["inner_spread_pct"],
                             "naive": report["naive"]["inner_spread_pct"]},
    }

    # ------------------------------------------------- checks (fail loud)
    failures = []
    san = report["sanity_constant_mu_annulus"]
    if san["inner_circularity_dev"] > 1e-6 or san["outer_circularity_dev"] > 1e-6:
        failures.append(f"constant-mu walls are not circles: {san}")
    if san["ratio_rel_err"] > 1e-6:
        failures.append(f"constant-mu radius ratio off: {san['ratio_rel_err']:.3e}")
    for tag, dd in report["design"].items():
        if not dd["J_single_sign"]:
            failures.append(f"design {tag}: inverse map folds")
    fs = report["flux_scale_freedom"]
    if fs["measured_max_dev"] > 1e-6:
        failures.append(f"flux scale-freedom broken: dev {fs['measured_max_dev']:.2e}")
    mmf_d = abs(report["design"]["saturable"]["mmf_design_A"])
    for label, v in report["verify"].items():
        if v["inner_wall_rel_err_mean"] > 0.010:
            failures.append(f"{label}: inner wall mean err "
                            f"{v['inner_wall_rel_err_mean']:.4f} > 1.0%")
        if v["inner_wall_rel_err_max"] > 0.020:
            failures.append(f"{label}: inner wall max err "
                            f"{v['inner_wall_rel_err_max']:.4f} > 2.0%")
        if v["outer_wall_rel_err_mean"] > 0.015:
            failures.append(f"{label}: outer wall mean err "
                            f"{v['outer_wall_rel_err_mean']:.4f} > 1.5%")
        rel = abs(v["mmf_fem_A"] - mmf_d) / mmf_d
        v["mmf_rel_diff"] = rel
        if rel > 0.010:
            failures.append(f"{label}: MMF design {mmf_d:.4f} A vs FEM "
                            f"{v['mmf_fem_A']:.4f} A (rel {rel:.4f} > 1.0%)")
    if report["naive"]["cap_overshoot_pct"] <= 1.0:
        failures.append("naive baseline does NOT violate the cap -- the "
                        "comparison has no content as posed")
    report["golden"] = {"passed": not failures, "failures": failures}
    report["meta"] = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "hostname": platform.node(), "python_version": platform.python_version(),
        "purpose": "correctness validation only (no timing claims)",
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results_ipm_bridge_free_boundary.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"\nresults -> {out}")
    cmp_ = report["comparison"]
    print(f"  hodograph vs naive: cap overshoot "
          f"{cmp_['cap_overshoot_pct']['hodograph']:+.2f}% vs "
          f"{cmp_['cap_overshoot_pct']['naive']:+.2f}%;  inner-wall spread "
          f"{cmp_['inner_spread_pct']['hodograph']:.2f}% vs "
          f"{cmp_['inner_spread_pct']['naive']:.2f}%;  naive uses "
          f"{cmp_['body_iron_area_mm2']['naive_excess_pct']:+.2f}% more iron")
    if failures:
        for f_ in failures:
            print("CHECK FAIL:", f_)
        raise SystemExit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
