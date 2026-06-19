# -*- coding: utf-8 -*-
# DEMO 1 (verified, verdict=confirms): hp vs corner singularity.
# L-shape 270-deg re-entrant corner, exact u = r^(2/3) sin(2 th/3), th in [0,3pi/2].
# Measure H1-seminorm error vs DOF: h-version (p=1) ~ N^-1/3, p-version ~ N^-2/3 (both algebraic).
import numpy as np
import ngsolve as ng
from netgen.geom2d import SplineGeometry


def make_lshape(maxh):
    geo = SplineGeometry()
    pts = [(0, 0), (1, 0), (1, 1), (-1, 1), (-1, -1), (0, -1)]
    pnums = [geo.AppendPoint(*p) for p in pts]
    for i in range(len(pnums)):
        geo.Append(["line", pnums[i], pnums[(i + 1) % len(pnums)]], bc="dir")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def exact_fields():
    x, y = ng.x, ng.y
    r = ng.sqrt(x * x + y * y)
    th0 = ng.atan2(y, x)
    th = ng.IfPos(-th0, th0 + 2 * np.pi, th0)  # shift to [0, 3pi/2)
    uex = r ** (2.0 / 3.0) * ng.sin(2.0 * th / 3.0)
    u_r = (2.0 / 3.0) * r ** (-1.0 / 3.0) * ng.sin(2.0 * th / 3.0)
    u_th = (2.0 / 3.0) * r ** (2.0 / 3.0) * ng.cos(2.0 * th / 3.0)
    ux = ng.cos(th) * u_r - (ng.sin(th) / r) * u_th
    uy = ng.sin(th) * u_r + (ng.cos(th) / r) * u_th
    return uex, ng.CoefficientFunction((ux, uy))


def solve(mesh, order):
    uex, graduex = exact_fields()
    fes = ng.H1(mesh, order=order, dirichlet="dir")
    gfu = ng.GridFunction(fes)
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    f = ng.LinearForm(fes); f.Assemble()
    gfu.Set(uex, ng.BND)
    res = f.vec.CreateVector(); res.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * res
    diff = ng.grad(gfu) - graduex
    return np.sqrt(ng.Integrate(diff * diff, mesh, order=2 * order + 8)), fes.ndof


def fit_alpha(Ns, errs):
    lN = np.log(np.array(Ns, float)); le = np.log(np.array(errs, float))
    A = np.vstack([np.ones_like(lN), -lN]).T
    return np.linalg.lstsq(A, le, rcond=None)[0][1]


if __name__ == "__main__":
    print("=== (A) h-refinement, order p=1 ===")
    hN, hE = [], []
    for maxh in [0.4, 0.2, 0.1, 0.05, 0.025]:
        e, n = solve(make_lshape(maxh), 1); hN.append(n); hE.append(e)
        print(f"  maxh={maxh:7.4f}  N={n:7d}  err={e:.6e}")
    print("=== (B) p-refinement, fixed coarse mesh maxh~0.3 ===")
    mesh = make_lshape(0.3); pN, pE = [], []
    for p in [1, 2, 3, 4, 5, 6]:
        e, n = solve(mesh, p); pN.append(n); pE.append(e)
        print(f"  order={p}  N={n:7d}  err={e:.6e}")
    ah = fit_alpha(hN, hE); ap = fit_alpha(pN, pE)
    print(f"alpha_h={ah:.4f} (~1/3)  alpha_p={ap:.4f} (~2/3)  ratio alpha_p/alpha_h={ap/ah:.2f} (asymptotic 2)")

    # --- verification (backs manuscript reentrant-corner exponents) ---
    N_FAIL = 0
    def check(name, cond, detail=""):
        global N_FAIL
        if not cond: N_FAIL += 1
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")
    check("h-version rate alpha_h ~ 1/3 (270deg reentrant corner, lambda=2/3, p=1)",
          0.30 < ah < 0.40, f"alpha_h = {ah:.4f}")
    check("p-version rate alpha_p ~ 2/3", 0.60 < ap < 0.72, f"alpha_p = {ap:.4f}")
    check("p-version ~ 2x h-version (measured ~1.85, asymptotic 2; NOT exactly 2)",
          1.6 < ap / ah < 2.1, f"alpha_p/alpha_h = {ap/ah:.2f}")
    assert N_FAIL == 0, f"{N_FAIL} checks failed"
