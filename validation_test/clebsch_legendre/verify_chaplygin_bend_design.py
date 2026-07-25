"""Saturable 90-degree bend: hodograph design -> physical shape -> nonlinear FEM check.

End-to-end verification that the field-plane (Chaplygin) hodograph DESIGN
direction works for a saturable material: one linear solve in hodograph
coordinates produces a physical shape, and an independent nonlinear FEM on
that shape reproduces the design spec.

Step (1) physically grounded design case
  material   mu_r(B) = 1 + (mu_r0-1)/(1+(B/Bk)^2)
             (same saturating model as the flux-guide comparison in
             docs/clebsch_hodograph/demos/chaplygin_hodograph_2d.py)
  spec       90 deg turn, walls are flux lines (A = const)
             outer wall field  B_out            = 1.00 T (constant)
             inner wall field  B_in(theta)      = 1.30 -> 1.75 T (tapered)
             flux per unit depth  dA            = 0.05 Wb/m
  design domain in hodograph space is the curvilinear trapezoid
             Omega_h = { (B, theta) : 0<=theta<=Theta,  B_out <= B <= B_in(theta) }

Step (2) forward verification
  the designed wall curves become a physical 2-D region; solve the nonlinear
  magnetostatic problem  div( nu(|grad A|) grad A ) = 0  there with the SAME
  Dirichlet flux values, and compare |B| on the walls against the design spec.

Formulation (B-radial A-form; every relation is locked symbolically by
site_builder verify_geometry05_exact_solution.wls and re-derived in the
companion verify_clebsch_legendre_transform.py family):
  d/dB( a A_B ) + b A_thth = 0,   a = B mu_d / mu_s^2,   b = 1/(mu_s B)
  Psi_B = -b A_th,                Psi_th = a A_B          (conjugate MMF potential)
  dr    = (dPsi/q) e_H + (dA/B) e_perp,   q = B/mu_s = |H|

Golden bands asserted at the end of the run (2026-07-23 baseline, LAB):
  constant-mu sanity   : designed walls are an exact annulus (dev < 1e-6)
  wall |B| (5..85 deg) : mean rel err < 1.0 % on the inner wall,
                         < 1.5 % on the outer wall, both mesh resolutions
  MMF                  : |design - FEM| / design < 0.5 %
  orientation          : J keeps one sign on every design (no folding)

Run:  python verify_chaplygin_bend_design.py
Writes results_chaplygin_bend_design.json next to this file (committed).
"""
import json
import math
import os
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

MU0 = 4.0e-7 * math.pi
MUR0 = 200.0
BK = 1.0

THETA = 0.5 * math.pi
B_OUT = 1.00
B_IN0 = 1.30
B_IN1 = 1.75
DFLUX = 0.05          # Wb/m, flux per unit depth between the two walls
NSAMP = 121           # samples along each wall


# ------------------------------------------------------------------ material
def mu_r_of(B):
    return 1.0 + (MUR0 - 1.0) / (1.0 + (B / BK) ** 2)


def mu_s_of(B):
    return MU0 * mu_r_of(B)


def dmu_s_dB_of(B):
    return MU0 * (-(MUR0 - 1.0) * 2.0 * B / BK**2 / (1.0 + (B / BK) ** 2) ** 2)


def mu_d_of(B):
    return mu_s_of(B) ** 2 / (mu_s_of(B) - B * dmu_s_dB_of(B))


def b_in_of(th):
    return B_IN0 + (B_IN1 - B_IN0) * th / THETA


# ------------------------------------------------------------------ helpers
def recover_potential(mesh, F, order=3):
    """Return u with grad(u) = F (Galerkin projection, gauge: u = 0 at one DOF).

    The gauge is fixed by pinning a single DOF rather than by a NumberSpace
    constraint: a NumberSpace row couples every DOF and destroys sparsity, which
    turns the direct factorization into a dense one.
    """
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


def _self_intersections(loop):
    """Segment pairs of the closed polyline that properly cross."""
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


# ------------------------------------------------------------------ step (1)
def design(report, mu_r_const=None, taper=True, tag=None):
    """Solve the hodograph design problem; return sampled wall curves."""
    if tag is None:
        tag = ("linear" if mu_r_const else "saturable") + ("" if taper else "_flat")
    b_in_1 = B_IN1 if taper else B_IN0

    def b_in(th):
        return B_IN0 + (b_in_1 - B_IN0) * th / THETA

    def mus(B):
        return MU0 * mu_r_const if mu_r_const else mu_s_of(B)

    def mud(B):
        return MU0 * mu_r_const if mu_r_const else mu_d_of(B)

    geo = SplineGeometry()
    p = [geo.AppendPoint(*pt) for pt in
         ((B_OUT, 0.0), (B_IN0, 0.0), (b_in_1, THETA), (B_OUT, THETA))]
    geo.Append(["line", p[0], p[1]], bc="inlet", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[1], p[2]], bc="inner", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], bc="outlet", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], bc="outer", leftdomain=1, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=0.03))
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
    gfA.Set(mesh.BoundaryCF({"outer": 0.0, "inner": DFLUX}, default=0.0), BND)
    res = gfA.vec.CreateVector()
    res.data = -a.mat * gfA.vec
    gfA.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res

    log(f"design/{tag}: A solved")
    A_B, A_t = grad(gfA)[0], grad(gfA)[1]
    Psi_B, Psi_t = -bB * A_t, aB * A_B
    gfPsi = recover_potential(mesh, CoefficientFunction((Psi_B, Psi_t)))
    log(f"design/{tag}: Psi recovered")

    c, s = cos(y), sin(y)
    Fx = CoefficientFunction(((Psi_B / qB) * c - (A_B / x) * s,
                              (Psi_t / qB) * c - (A_t / x) * s))
    Fy = CoefficientFunction(((Psi_B / qB) * s + (A_B / x) * c,
                              (Psi_t / qB) * s + (A_t / x) * c))
    gfx = recover_potential(mesh, Fx)
    gfy = recover_potential(mesh, Fy)
    log(f"design/{tag}: coordinates recovered")

    # Jacobian of the inverse map (must keep one sign: no folding)
    Jcf = (grad(gfx)[0] * grad(gfy)[1] - grad(gfx)[1] * grad(gfy)[0])

    ths = np.linspace(0.0, THETA, NSAMP)
    eps = 1e-9
    inner, outer, Jvals = [], [], []
    for th in ths:
        bi = min(b_in(th), b_in_1) - eps
        pi_ = mesh(max(B_OUT + eps, bi), min(max(th, eps), THETA - eps))
        po_ = mesh(B_OUT + eps, min(max(th, eps), THETA - eps))
        inner.append((gfx(pi_), gfy(pi_)))
        outer.append((gfx(po_), gfy(po_)))
        Jvals.append(Jcf(pi_))
        Jvals.append(Jcf(po_))
    inner = np.array(inner)
    outer = np.array(outer)

    # end faces (Psi = const terminals)
    bs0 = np.linspace(B_OUT + eps, B_IN0 - eps, 60)
    bs1 = np.linspace(B_OUT + eps, b_in_1 - eps, 60)
    inlet = np.array([(gfx(mesh(b, eps)), gfy(mesh(b, eps))) for b in bs0])
    outlet = np.array([(gfx(mesh(b, THETA - eps)), gfy(mesh(b, THETA - eps)))
                       for b in bs1])

    mmf = float(gfPsi(mesh(0.5 * (B_OUT + b_in_1), THETA - eps))
                - gfPsi(mesh(0.5 * (B_OUT + B_IN0), eps)))
    Jarr = np.array(Jvals)
    out = {
        "tag": tag,
        "J_single_sign": bool(np.all(Jarr < 0) or np.all(Jarr > 0)),
        "min_absJ": float(np.min(np.abs(Jarr))),
        "mmf_design_A": mmf,
        "inner_wall_m": inner.tolist(),
        "outer_wall_m": outer.tolist(),
    }
    width0 = float(np.linalg.norm(inner[0] - outer[0]))
    width1 = float(np.linalg.norm(inner[-1] - outer[-1]))
    out["inlet_width_mm"] = 1e3 * width0
    out["outlet_width_mm"] = 1e3 * width1
    print(f"  [design/{tag}] J single sign={out['J_single_sign']}  "
          f"min|J|={out['min_absJ']:.3e}")
    print(f"  [design/{tag}] inlet width={1e3*width0:.2f} mm  "
          f"outlet width={1e3*width1:.2f} mm  MMF={mmf:.1f} A")
    report.setdefault("design", {})[tag] = out
    return inner, outer, inlet, outlet, ths


def annulus_sanity(inner, outer, report):
    """Constant-mu design must come out as an exact circular annulus."""
    def fit_circle(pts):
        A = np.c_[2 * pts[:, 0], 2 * pts[:, 1], np.ones(len(pts))]
        bb = (pts**2).sum(axis=1)
        sol, *_ = np.linalg.lstsq(A, bb, rcond=None)
        cx, cy = sol[0], sol[1]
        r = math.sqrt(sol[2] + cx * cx + cy * cy)
        dev = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - r).max() / r
        return (cx, cy), r, dev
    ci, ri, di = fit_circle(inner)
    co, ro, do = fit_circle(outer)
    ratio = ro / ri
    expect = B_IN0 / B_OUT          # constant mu: B ~ 1/r  =>  r_out/r_in = B_in/B_out
    res = {
        "inner_radius_mm": 1e3 * ri, "outer_radius_mm": 1e3 * ro,
        "inner_circularity_dev": di, "outer_circularity_dev": do,
        "radius_ratio": ratio, "expected_ratio_Bin_over_Bout": expect,
        "ratio_rel_err": abs(ratio - expect) / expect,
        "center_offset_mm": 1e3 * float(np.hypot(ci[0] - co[0], ci[1] - co[1])),
    }
    print(f"  [sanity] constant-mu design -> annulus: r_in={1e3*ri:.2f} mm, "
          f"r_out={1e3*ro:.2f} mm")
    print(f"  [sanity] circularity dev: inner={di:.2e} outer={do:.2e}; "
          f"r_out/r_in={ratio:.6f} vs B_in/B_out={expect:.6f} "
          f"(rel err {res['ratio_rel_err']:.2e})")
    report["sanity_constant_mu_annulus"] = res


# ------------------------------------------------------------------ step (2)
def forward_verify(inner, outer, inlet, outlet, ths, report, maxh=None,
                   label=None, b_in_target=None):
    """Mesh the designed shape and solve the nonlinear physical problem."""
    # Corners are sampled twice (once by each adjacent curve); dropping the
    # shared endpoints is what keeps the closed loop free of the degenerate
    # crossing segments that stall the mesher.
    parts = [(outer, "wall_out"),
             (outlet[1:], "outlet"),
             (inner[::-1][1:], "wall_in"),
             (inlet[::-1][1:-1], "inlet")]
    loop = np.vstack([p for p, _ in parts])
    tags = [t for p, t in parts for _ in range(len(p))]
    # drop near-duplicate consecutive points (tolerance: 10 um)
    keep = [0]
    for i in range(1, len(loop)):
        if np.linalg.norm(loop[i] - loop[keep[-1]]) > 1e-5:
            keep.append(i)
    loop = loop[keep]
    tags = [tags[i] for i in keep]
    n = len(loop)
    hits = _self_intersections(loop)
    if hits:
        raise RuntimeError(
            f"designed outline self-intersects at segment pairs {hits[:5]} "
            f"({len(hits)} total); the design is not manufacturable as posed")
    if polygon_area(loop[:, 0], loop[:, 1]) < 0:
        # segment i of the reversed loop is the reverse of original segment n-2-i
        loop = loop[::-1]
        tags = [tags[(n - 2 - i) % n] for i in range(n)]
        order_flag = "reversed"
    else:
        order_flag = "forward"

    geo = SplineGeometry()
    pids = [geo.AppendPoint(float(px), float(py)) for px, py in loop]
    for i in range(n):
        geo.Append(["line", pids[i], pids[(i + 1) % n]], bc=tags[i],
                   leftdomain=1, rightdomain=0)
    width = np.linalg.norm(inner[0] - outer[0])
    mesh = Mesh(geo.GenerateMesh(maxh=maxh or width / 8.0))
    log(f"verify: mesh {mesh.ne} elements, {mesh.nv} vertices "
        f"({order_flag} loop), boundaries={set(mesh.GetBoundaries())}")

    fes = H1(mesh, order=3, dirichlet="wall_in|wall_out")
    u, v = fes.TnT()
    gfA = GridFunction(fes)
    gfA.Set(mesh.BoundaryCF({"wall_out": 0.0, "wall_in": DFLUX}, default=0.0), BND)
    gfPrev = GridFunction(fes)
    gfTrial = GridFunction(fes)
    eps = 1e-6
    omega = 0.35           # under-relaxation; plain Picard oscillates here
    hist = []
    for _ in range(400):
        gfPrev.vec.data = gfA.vec
        Bmag = sqrt(grad(gfA) * grad(gfA) + eps**2)
        nu = 1.0 / (MU0 * (1.0 + (MUR0 - 1.0) / (1.0 + (Bmag / BK) ** 2)))
        a = BilinearForm(fes)
        a += nu * grad(u) * grad(v) * dx
        a.Assemble()
        gfTrial.vec.data = gfA.vec
        res = gfA.vec.CreateVector()
        res.data = -a.mat * gfTrial.vec
        gfTrial.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                          inverse="sparsecholesky") * res
        gfA.vec.data = gfPrev.vec + omega * (gfTrial.vec - gfPrev.vec)
        d = gfPrev.vec.CreateVector()
        d.data = gfA.vec - gfPrev.vec
        rel = d.Norm() / max(gfA.vec.Norm(), 1e-30)
        hist.append(rel)
        if rel < 1e-9:
            break
    log(f"verify: Picard {len(hist)} iterations, final rel step {hist[-1]:.2e}")
    if hist[-1] > 1e-7:
        raise RuntimeError(
            f"forward Picard did not converge (rel step {hist[-1]:.2e} after "
            f"{len(hist)} iterations); the verification number would be meaningless")

    # independent global check: MMF = integral of H.dl along the guide mid-line
    Bvec = CoefficientFunction((grad(gfA)[1], -grad(gfA)[0]))
    Bmag_cf = sqrt(grad(gfA) * grad(gfA) + eps**2)
    nu_cf = 1.0 / (MU0 * (1.0 + (MUR0 - 1.0) / (1.0 + (Bmag_cf / BK) ** 2)))
    Hvec = nu_cf * Bvec
    mid = 0.5 * (inner + outer)
    mmf_fem = 0.0
    for i in range(len(mid) - 1):
        seg = mid[i + 1] - mid[i]
        pm = 0.5 * (mid[i] + mid[i + 1])
        hv = Hvec(mesh(pm[0], pm[1]))
        mmf_fem += hv[0] * seg[0] + hv[1] * seg[1]
    mmf_fem = abs(float(mmf_fem))

    Bcf = sqrt(grad(gfA) * grad(gfA))
    b_in_fem, b_out_fem = [], []
    for i, th in enumerate(ths):
        b_in_fem.append(Bcf(mesh(inner[i][0], inner[i][1])))
        b_out_fem.append(Bcf(mesh(outer[i][0], outer[i][1])))
    b_in_fem = np.array(b_in_fem)
    b_out_fem = np.array(b_out_fem)
    b_in_tgt = b_in_of(ths) if b_in_target is None else np.asarray(b_in_target)
    b_out_tgt = np.full_like(ths, B_OUT)

    # ignore the corner samples: |B| evaluated exactly at a polygon corner is
    # meaningless (the two Dirichlet walls meet the natural end face there)
    sl = slice(6, -6)
    e_in = np.abs(b_in_fem[sl] - b_in_tgt[sl]) / b_in_tgt[sl]
    e_out = np.abs(b_out_fem[sl] - b_out_tgt[sl]) / b_out_tgt[sl]
    res = {
        "n_elements": mesh.ne, "picard_iterations": len(hist),
        "inner_wall_rel_err_max": float(e_in.max()),
        "inner_wall_rel_err_mean": float(e_in.mean()),
        "outer_wall_rel_err_max": float(e_out.max()),
        "outer_wall_rel_err_mean": float(e_out.mean()),
        "mmf_fem_A": mmf_fem,
        "inner_B_fem_T": b_in_fem.tolist(),
        "outer_B_fem_T": b_out_fem.tolist(),
        "target_inner_B_T": b_in_tgt.tolist(),
    }
    log(f"verify: MMF  FEM={mmf_fem:.1f} A")
    print(f"  [verify] inner wall |B|: design {b_in_tgt[0]:.3f}->{b_in_tgt[-1]:.3f} T, "
          f"FEM {b_in_fem[3]:.3f}->{b_in_fem[-4]:.3f} T")
    print(f"  [verify] outer wall |B|: design {B_OUT:.3f} T constant, "
          f"FEM {b_out_fem[3]:.3f}..{b_out_fem[-4]:.3f} T")
    print(f"  [verify] rel err  inner mean={e_in.mean()*100:.3f}% max={e_in.max()*100:.3f}%"
          f"   outer mean={e_out.mean()*100:.3f}% max={e_out.max()*100:.3f}%")
    report.setdefault("verify", {})[label or f"maxh_{maxh}"] = res
    return res


def main():
    SetNumThreads(4)
    report = {"case": {
        "material": "mu_r(B)=1+(mu_r0-1)/(1+(B/Bk)^2)", "mu_r0": MUR0, "Bk_T": BK,
        "turn_deg": 90.0, "B_outer_T": B_OUT,
        "B_inner_T": [B_IN0, B_IN1], "flux_Wb_per_m": DFLUX,
    }}
    with TaskManager():
        print("step (1) design -- constant-mu sanity check (flat inner wall)")
        li, lo, _, _, _ = design(report, mu_r_const=100.0, taper=False)
        annulus_sanity(li, lo, report)

        # bisection: does the saturable path work at all without the taper?
        print("step (1) design -- saturable, FLAT inner wall (bisection case)")
        fi, fo, fin, fout, fths = design(report, taper=False)
        print("step (2) forward verification -- saturable flat")
        wf = float(np.linalg.norm(fi[0] - fo[0]))
        forward_verify(fi, fo, fin, fout, fths, report, maxh=wf / 8.0,
                       label="saturable_flat", b_in_target=np.full_like(fths, B_IN0))

        print("step (1) design -- saturable, tapered inner wall")
        inner, outer, inlet, outlet, ths = design(report)
        print("step (2) forward verification -- saturable tapered")
        width = float(np.linalg.norm(inner[0] - outer[0]))
        for div in (8.0, 16.0):
            forward_verify(inner, outer, inlet, outlet, ths, report,
                           maxh=width / div, label=f"saturable_taper_h{div:g}")
            log(f"verify: resolution width/{div:g} done")

    # ---------------- golden bands (fail loud) ----------------
    failures = []
    san = report["sanity_constant_mu_annulus"]
    if san["inner_circularity_dev"] > 1e-6 or san["outer_circularity_dev"] > 1e-6:
        failures.append(f"constant-mu walls are not circles: {san}")
    if san["ratio_rel_err"] > 1e-6:
        failures.append(f"constant-mu radius ratio off: {san['ratio_rel_err']:.3e}")
    for tag, d in report["design"].items():
        if not d["J_single_sign"]:
            failures.append(f"design {tag}: inverse map folds (J changes sign)")
    for label, v in report["verify"].items():
        if v["inner_wall_rel_err_mean"] > 0.010:
            failures.append(f"{label}: inner wall mean err "
                            f"{v['inner_wall_rel_err_mean']:.4f} > 1.0%")
        if v["outer_wall_rel_err_mean"] > 0.015:
            failures.append(f"{label}: outer wall mean err "
                            f"{v['outer_wall_rel_err_mean']:.4f} > 1.5%")
    for label, design_tag in (("saturable_flat", "saturable_flat"),
                              ("saturable_taper_h8", "saturable"),
                              ("saturable_taper_h16", "saturable")):
        mmf_d = report["design"][design_tag]["mmf_design_A"]
        mmf_f = report["verify"][label]["mmf_fem_A"]
        rel = abs(mmf_f - abs(mmf_d)) / abs(mmf_d)
        report["verify"][label]["mmf_rel_diff"] = rel
        if rel > 0.005:
            failures.append(f"{label}: MMF design {mmf_d:.1f} A vs FEM "
                            f"{mmf_f:.1f} A (rel {rel:.4f} > 0.5%)")
    report["golden"] = {"passed": not failures, "failures": failures}

    import datetime
    import platform
    report["meta"] = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "purpose": "correctness validation only (no timing claims)",
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results_chaplygin_bend_design.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"results -> {out}")
    if failures:
        for f_ in failures:
            print("GOLDEN FAIL:", f_)
        raise SystemExit(1)
    print("all golden bands passed")


if __name__ == "__main__":
    main()
