"""Does the H1 mesh verdict survive on VECTOR elements (HCurl / HDiv)?

Promoted from C:/temp/mesh_quality_study (2026-08-07) with its committed
results JSON (Data Persistence Policy). Re-run with
`python run_vector_elements.py` (requires Cubit + netgen + gmsh +
build123d).

run_accuracy_per_dof.py answered "which mesh is better" for an isotropic
scalar Poisson problem at H1 order 1, and its README explicitly refuses
to carry that verdict to HCurl/HDiv -- which is precisely what this lab
solves (HDiv-VIM, eddy currents). This closes that hole.

Two manufactured problems, both posed identically on every mesh:

  HCurl (curl-curl + mass, the eddy-current shape):
      u_ex = (sin ky, sin kz, sin kx)   ->   curl curl u_ex = k^2 u_ex
      solve  curl curl u + u = (k^2+1) u_ex,  tangential trace from u_ex

  HDiv (div-div + mass, the flux-space shape):
      u_ex = (sin kx, sin ky, sin kz)   ->   -grad div u_ex = k^2 u_ex
      solve  -grad div u + u = (k^2+1) u_ex,  normal trace from u_ex

Lowest order in each space (Whitney edge / RT face), Jacobi-CG, and the
residual is CHECKED -- a run that hits maxsteps is recorded as not
converged and excluded from fits rather than reported as a result.

Observables per route: ndof, free-dof fraction, CG iterations, relative
L2 error, and relative error in the natural (curl / div) seminorm.
"""
import json
import math
import os
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "packages" / "radia-mcp" / "src"))

from build123d import Box, Sphere, export_step
from radia_mcp.build123d.archetypes import c_core
from radia_mcp.cubit.server import _run_batch
from radia_mcp.gmsh.msh_inspect import mesh_quality

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "artifacts")
os.makedirs(OUT, exist_ok=True)
MAXSTEPS = 30000
MIN_FREE_FRAC = 0.20


def _cf(*comps):
    from ngsolve import CoefficientFunction
    return CoefficientFunction(tuple(comps))


def solve_hcurl(mesh, k):
    from ngsolve import (BND, BilinearForm, CGSolver, GridFunction, HCurl,
                         Integrate, InnerProduct, LinearForm, Preconditioner,
                         cos, curl, dx, sin, sqrt, x, y, z)
    uex = _cf(sin(k * y), sin(k * z), sin(k * x))
    cuex = _cf(-k * cos(k * z), -k * cos(k * x), -k * cos(k * y))
    fes = HCurl(mesh, order=0, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += (InnerProduct(curl(u), curl(v)) + InnerProduct(u, v)) * dx
    f = LinearForm(fes)
    f += (k * k + 1.0) * InnerProduct(uex, v) * dx
    return _solve(fes, a, f, mesh, uex, cuex, curl, "curl")


def solve_hdiv(mesh, k):
    from ngsolve import (BilinearForm, CGSolver, GridFunction, HDiv, Integrate,
                         InnerProduct, LinearForm, Preconditioner, cos, div,
                         dx, sin, sqrt, x, y, z)
    uex = _cf(sin(k * x), sin(k * y), sin(k * z))
    duex = k * (cos(k * x) + cos(k * y) + cos(k * z))
    fes = HDiv(mesh, order=0, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += (div(u) * div(v) + InnerProduct(u, v)) * dx
    f = LinearForm(fes)
    f += (k * k + 1.0) * InnerProduct(uex, v) * dx
    return _solve(fes, a, f, mesh, uex, duex, div, "div")


def _solve(fes, a, f, mesh, uex, dex, deriv, dname):
    from ngsolve import (BND, CGSolver, GridFunction, Integrate, InnerProduct,
                         Preconditioner, sqrt)
    jac = Preconditioner(a, "local")
    a.Assemble()
    f.Assemble()
    free = sum(fes.FreeDofs())
    if free == 0:
        raise RuntimeError(f"degenerate: ndof={fes.ndof} but 0 free dofs")
    gfu = GridFunction(fes)
    gfu.Set(uex, BND)
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    r0 = res.Norm()
    inv = CGSolver(a.mat, jac.mat, precision=1e-10, maxsteps=MAXSTEPS,
                   printrates=False)
    gfu.vec.data += inv * res
    # verify the solve rather than trusting the iteration count
    chk = f.vec.CreateVector()
    chk.data = f.vec - a.mat * gfu.vec
    for i in range(len(chk)):
        if not fes.FreeDofs()[i]:
            chk[i] = 0.0
    rel_res = chk.Norm() / max(r0, 1e-300)
    converged = inv.GetSteps() < MAXSTEPS and rel_res < 1e-6

    def _err(cf_diff):
        return sqrt(Integrate(InnerProduct(cf_diff, cf_diff), mesh, order=4))

    l2 = _err(gfu - uex)
    l2n = _err(uex)
    d = _err(deriv(gfu) - dex)
    dn = _err(dex)
    return {"ndof": fes.ndof, "free_dofs": free,
            "free_frac": free / fes.ndof, "ne": mesh.ne,
            "cg_iters": inv.GetSteps(), "rel_residual": rel_res,
            "converged": bool(converged),
            "l2_rel": l2 / l2n, f"{dname}_rel": d / dn}


def quality_min(msh):
    q = mesh_quality(msh)
    bts = q.get("by_type") or []
    if not bts:
        return None
    return round(min(b["min_quality"] for b in bts), 4)


def mesh_cubit(step, size, tag, hex_mode):
    from ngsolve import Mesh
    vol = os.path.join(OUT, f"ve_{tag}.vol")
    msh = os.path.join(OUT, f"ve_{tag}.msh")
    cmds = [] if hex_mode else ["volume all scheme tetmesh"]
    cmds += [f"volume all size {size}", "mesh volume all",
             "block 1 add volume all", 'block 1 name "mesh"',
             f'export netgen "{vol.replace(os.sep, "/")}" order 1',
             f'export gmsh "{msh.replace(os.sep, "/")}" overwrite']
    r = _run_batch(step, cmds, timeout_s=900)
    assert r["status"] == "ok", r
    return Mesh(vol), quality_min(msh)


def mesh_netgen(step, maxh, tag):
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    from radia_mcp.cubit.server import _netgen_mesh_to_msh
    from pathlib import Path
    msh = Path(OUT) / f"ve_{tag}_ng.msh"
    _netgen_mesh_to_msh(Path(step), maxh, msh)
    return Mesh(OCCGeometry(step).GenerateMesh(maxh=maxh)), quality_min(msh)


def fit(pts, key):
    pts = [p for p in pts
           if p.get("converged") and p["free_frac"] >= MIN_FREE_FRAC]
    if len(pts) < 3:
        return {"slope": None, "n_points": len(pts),
                "note": "too few converged non-degenerate points"}
    xs = [math.log(p["ndof"]) for p in pts]
    ys = [math.log(p[key]) for p in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    slope = (sum((s - mx) * (t - my) for s, t in zip(xs, ys))
             / sum((t - mx) ** 2 for t in xs))
    return {"slope": round(slope, 3), "intercept": round(my - slope * mx, 3),
            "n_points": n}


def main():
    t0 = time.time()
    cc = os.path.join(OUT, "c_core.step")
    sp = os.path.join(OUT, "sphere.step")
    tp = os.path.join(OUT, "thin_plate.step")
    export_step(c_core(width=80, height=60, depth=25, leg=15, gap=8), cc)
    export_step(Sphere(1.0), sp)
    export_step(Box(8, 8, 1), tp)

    plan = [
        ("sphere", sp, 1.0, [0.35, 0.28, 0.22, 0.18],
         ("cubit_tet", "netgen")),
        ("c_core", cc, 0.02, [8.0, 6.0, 5.0, 4.0],
         ("cubit_tet", "netgen")),
        ("thin_plate", tp, 0.4, [0.5, 0.4, 0.32, 0.25],
         ("cubit_tet", "cubit_hex", "netgen")),
    ]

    from ngsolve import TaskManager
    results = {"timestamp": datetime.now().isoformat(),
               "hostname": platform.node(),
               "problem": ("HCurl curl-curl+mass and HDiv div-div+mass, "
                           "lowest order, manufactured solutions; "
                           "Jacobi-CG with a verified residual"),
               "maxsteps": MAXSTEPS, "spaces": {}}

    with TaskManager():
        for space, solver, dname in (("hcurl", solve_hcurl, "curl"),
                                     ("hdiv", solve_hdiv, "div")):
            geo_out = {}
            for name, step, k, sizes, routes in plan:
                g = {"k": k}
                for route in routes:
                    rows = []
                    for i, s in enumerate(sizes):
                        tag = f"{name}_{route}_{i}"
                        if route == "netgen":
                            m, qmin = mesh_netgen(step, s, tag)
                        else:
                            m, qmin = mesh_cubit(step, s, tag,
                                                 route == "cubit_hex")
                        row = solver(m, k)
                        row["h"] = s
                        row["min_quality"] = qmin
                        rows.append(row)
                        print(f"{space:5s} {name:10s} {route:9s} h={s:6.3f} "
                              f"ne={row['ne']:7d} ndof={row['ndof']:7d} "
                              f"free={row['free_frac']*100:3.0f}% "
                              f"min={qmin} "
                              f"L2={row['l2_rel']*100:8.3f}% "
                              f"{dname}={row[dname+'_rel']*100:8.3f}% "
                              f"CG={row['cg_iters']:6d}"
                              f"{'' if row['converged'] else '  NOT-CONVERGED'}",
                              flush=True)
                    g[route] = rows
                    for key in ("l2_rel", f"{dname}_rel"):
                        g[f"{route}_{key}_fit"] = fit(rows, key)
                geo_out[name] = g
                print(flush=True)
            results["spaces"][space] = geo_out

    out = os.path.join(_HERE, "results_vector_elements.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"saved {out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
