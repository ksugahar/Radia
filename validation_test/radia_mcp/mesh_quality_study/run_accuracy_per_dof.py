"""Which mesh is more ACCURATE for the same cost (dof)? tet vs tet vs hex.

Promoted from C:/temp/mesh_quality_study (2026-08-06) with its committed
results JSON (Data Persistence Policy). Re-run with
`python run_accuracy_per_dof.py` (requires Cubit + netgen + build123d;
scratch meshes land in artifacts/, gitignored).

The equal-budget pairs in results_solver_impact.json give one point per
mesher, which cannot support an accuracy-per-dof ranking: at equal
ELEMENT budget the meshers land at different dof counts. And the hex
routes were never solved on at all -- study 1 only scored their minSICN
and element count, which says nothing about accuracy.

This script measures the actual convergence curve -- relative L2 / H1
error vs ndof -- for every route on the same STEP, over ~1 decade of
ndof, and compares the curves at a MATCHED ndof.

Manufactured solution, per geometry:
    u_ex = sin(a x) sin(b y) exp(c z)  ->  -lap u_ex = (a^2+b^2-c^2) u_ex
H1 order 1, Dirichlet u_ex on all boundaries. (a,b,c) are chosen so the
solution genuinely varies over that geometry's extent -- in particular
the thin plate gets a real through-thickness variation, which is where
hex is supposed to earn its keep.
"""
import json
import math
import os
import platform
import sys
import time
from datetime import datetime

sys.path.insert(0, r"S:\Radia\01_GitHub\packages\radia-mcp\src")

from build123d import Box, Sphere, export_step
from radia_mcp.build123d.archetypes import c_core
from radia_mcp.cubit.server import _run_batch

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "artifacts")
os.makedirs(OUT, exist_ok=True)


def solve_poisson(mesh, abc):
    from ngsolve import (BND, BilinearForm, CGSolver, CoefficientFunction,
                         GridFunction, H1, Integrate, InnerProduct, LinearForm,
                         Preconditioner, cos, dx, exp, grad, sin, sqrt, x, y, z)
    a_, b_, c_ = abc
    uex = sin(a_ * x) * sin(b_ * y) * exp(c_ * z)
    guex = CoefficientFunction((
        a_ * cos(a_ * x) * sin(b_ * y) * exp(c_ * z),
        b_ * sin(a_ * x) * cos(b_ * y) * exp(c_ * z),
        c_ * sin(a_ * x) * sin(b_ * y) * exp(c_ * z)))
    lam = a_ * a_ + b_ * b_ - c_ * c_          # -lap u = lam * u

    fes = H1(mesh, order=1, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    f = LinearForm(fes)
    f += lam * uex * v * dx
    jac = Preconditioner(a, "local")
    a.Assemble()
    f.Assemble()
    free = sum(fes.FreeDofs())
    if free == 0:
        raise RuntimeError(
            f"degenerate case: ndof={fes.ndof} but 0 free dofs -- every node "
            f"is on a Dirichlet boundary, so nothing is solved and the "
            f"'error' is pure interpolation error. Refine so the geometry "
            f"has interior nodes.")
    gfu = GridFunction(fes)
    gfu.Set(uex, BND)
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    inv = CGSolver(a.mat, jac.mat, precision=1e-12, maxsteps=200000,
                   printrates=False)
    gfu.vec.data += inv * res
    l2 = sqrt(Integrate((gfu - uex) ** 2, mesh, order=4))
    h1 = sqrt(Integrate(InnerProduct(grad(gfu) - guex, grad(gfu) - guex),
                        mesh, order=4))
    return {"ndof": fes.ndof, "free_dofs": free, "interior_frac": free / fes.ndof,
            "ne": mesh.ne, "cg_iters": inv.GetSteps(),
            "l2_rel": l2 / sqrt(Integrate(uex ** 2, mesh, order=4)),
            "h1_rel": h1 / sqrt(Integrate(InnerProduct(guex, guex),
                                          mesh, order=4))}


def mesh_cubit(step, size, tag, hex_mode):
    from ngsolve import Mesh
    vol = os.path.join(OUT, f"ac_{tag}.vol")
    cmds = [] if hex_mode else ["volume all scheme tetmesh"]
    cmds += [f"volume all size {size}", "mesh volume all",
             "block 1 add volume all", 'block 1 name "mesh"',
             f'export netgen "{vol.replace(os.sep, "/")}" order 1']
    r = _run_batch(step, cmds, timeout_s=900)
    assert r["status"] == "ok", r
    return Mesh(vol)


def mesh_netgen(step, maxh):
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    return Mesh(OCCGeometry(step).GenerateMesh(maxh=maxh))


MIN_INTERIOR_FRAC = 0.20


def fit(pts, key):
    """least-squares fit of log(err) = intercept + slope*log(ndof).

    Points whose mesh is almost all boundary are EXCLUDED: there the
    discrete solution is dominated by the Dirichlet interpolant, the
    error does not respond to refinement, and including them silently
    corrupts the fitted rate.
    """
    pts = [p for p in pts if p["interior_frac"] >= MIN_INTERIOR_FRAC]
    if len(pts) < 3:
        return {"slope": None, "intercept": None, "n_points": len(pts),
                "note": "too few non-degenerate points to fit"}
    xs = [math.log(p["ndof"]) for p in pts]
    ys = [math.log(p[key]) for p in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    den = sum((t - mx) ** 2 for t in xs)
    slope = sum((s - mx) * (t - my) for s, t in zip(xs, ys)) / den
    return {"slope": round(slope, 3), "intercept": round(my - slope * mx, 3),
            "n_points": n, "ndof_range": [min(p["ndof"] for p in pts),
                                          max(p["ndof"] for p in pts)]}


def main():
    t0 = time.time()
    cc = os.path.join(OUT, "c_core.step")
    sp = os.path.join(OUT, "sphere.step")
    tp = os.path.join(OUT, "thin_plate.step")
    export_step(c_core(width=80, height=60, depth=25, leg=15, gap=8), cc)
    export_step(Sphere(1.0), sp)
    export_step(Box(8, 8, 1), tp)

    # name, step, (a,b,c), sizes, routes
    plan = [
        ("sphere", sp, (1.0, 1.0, 1.0),
         [0.35, 0.28, 0.22, 0.18, 0.14, 0.11],
         ("cubit_tet", "netgen")),
        ("c_core", cc, (0.02, 0.02, 0.02),
         [8.0, 6.0, 5.0, 4.0, 3.2, 2.6],
         ("cubit_tet", "netgen")),
        # thin plate 8x8x1: sizes chosen so the sweep RESOLVES the
        # thickness (2..5 elements through it). A coarser sweep on a
        # 50x50x1 plate produced meshes whose every node was on a
        # Dirichlet face -- zero free dofs, no solve, a flat "error"
        # curve of pure interpolation error. solve_poisson() now raises
        # on that case instead of reporting it as a result.
        ("thin_plate", tp, (0.4, 0.4, 1.5),
         [0.5, 0.4, 0.32, 0.25, 0.2],
         ("cubit_tet", "cubit_hex", "netgen")),
    ]

    from ngsolve import TaskManager
    results = {"timestamp": datetime.now().isoformat(),
               "hostname": platform.node(),
               "problem": ("H1 order-1 Poisson, manufactured "
                           "u=sin(ax)sin(by)exp(cz); relative L2/H1 error "
                           "vs ndof for every route on the same STEP"),
               "geometries": {}}

    with TaskManager():
        for name, step, abc, sizes, routes in plan:
            g = {"abc": abc, "sizes": sizes}
            for route in routes:
                rows = []
                for i, s in enumerate(sizes):
                    if route == "netgen":
                        m = mesh_netgen(step, s)
                    else:
                        m = mesh_cubit(step, s, f"{name}_{route}_{i}",
                                       route == "cubit_hex")
                    row = solve_poisson(m, abc)
                    row["h"] = s
                    rows.append(row)
                    print(f"{name:10s} {route:9s} h={s:6.3f} "
                          f"ne={row['ne']:7d} ndof={row['ndof']:7d} "
                          f"int={row['interior_frac']*100:4.0f}% "
                          f"L2={row['l2_rel']*100:8.4f}% "
                          f"H1={row['h1_rel']*100:8.4f}% "
                          f"CG={row['cg_iters']:5d}", flush=True)
                g[route] = rows
                for key in ("l2_rel", "h1_rel"):
                    g[f"{route}_{key}_fit"] = fit(rows, key)

            # compare every route at a common ndof (geometric mid of overlap)
            lo = max(min(p["ndof"] for p in g[r]) for r in routes)
            hi = min(max(p["ndof"] for p in g[r]) for r in routes)
            mid = math.sqrt(lo * hi)
            for key in ("l2_rel", "h1_rel"):
                vals = {}
                for r in routes:
                    fr = g[f"{r}_{key}_fit"]
                    if fr.get("slope") is None:
                        continue
                    vals[r] = math.exp(fr["intercept"]
                                       + fr["slope"] * math.log(mid))
                if "cubit_tet" not in vals:
                    g[f"{key}_at_matched_ndof"] = {"note": "no valid fit"}
                    continue
                base = vals["cubit_tet"]
                g[f"{key}_at_matched_ndof"] = {
                    "ndof": round(mid),
                    "values": vals,
                    "vs_cubit_tet": {r: round(v / base, 3)
                                     for r, v in vals.items()}}
                txt = "  ".join(f"{r}={v*100:.4f}%" for r, v in vals.items())
                print(f"{name:10s} {key} @ ndof={round(mid):6d}: {txt}",
                      flush=True)
            for r in routes:
                fl, fh = g[f"{r}_l2_rel_fit"], g[f"{r}_h1_rel_fit"]
                if fl.get("slope") is None:
                    print(f"{name:10s} {r:9s} rate: {fl['note']}", flush=True)
                    continue
                print(f"{name:10s} {r:9s} rate L2 {fl['slope']:+.3f}"
                      f"  H1 {fh['slope']:+.3f} (err ~ ndof^s, "
                      f"{fl['n_points']} pts, ndof {fl['ndof_range'][0]}-"
                      f"{fl['ndof_range'][1]})", flush=True)
            results["geometries"][name] = g
            print(flush=True)

    out = os.path.join(_HERE, "results_accuracy_per_dof.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"saved {out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
