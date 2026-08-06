"""Does mesh quality actually change SOLVER behaviour?

Promoted from C:/temp/mesh_quality_study (2026-08-06) with its committed
results JSON (Data Persistence Policy). Re-run with
`python run_solver_impact.py` (requires Cubit + netgen + gmsh +
build123d; scratch meshes land in artifacts/, gitignored).

run_study.py ranked the meshers by minSICN. This script asks whether that
ranking has any consequence, using a manufactured-solution Poisson problem
posed identically on every mesh:

    u_ex = sin(kx) sin(ky) exp(kz)   ->   -lap u_ex = k^2 u_ex

(k chosen per geometry so k*L ~ 1, keeping exp() well scaled), H1 order 1,
Dirichlet u_ex on all boundaries. Two deterministic observables, no timing:

  * CG iterations to a fixed relative precision, unpreconditioned
    (Projector on the free dofs) AND Jacobi-preconditioned -> CONVERGENCE
  * relative L2 / H1 error against u_ex                     -> ACCURACY

Three pairs:
  * sphere  equal-budget  -- cubit tetmesh 0.3 vs netgen maxh 0.2006
  * c_core  equal-budget  -- cubit tetmesh 4.0 vs netgen maxh 3.1024
  * c_core  netgen valley -- SAME mesher, refining, min quality falling
    0.51 -> 0.37 -> 0.27: does the worst element show up at all?

Cubit meshes travel the PRODUCTION route (.vol -> NGSolve).

Quality-class run (correctness, not timing) -- LAB execution allowed.
"""
import json
import os
import platform
import time
from datetime import datetime

from build123d import Sphere, export_step
from radia_mcp.build123d.archetypes import c_core
from radia_mcp.cubit.server import _run_batch
from radia_mcp.gmsh.msh_inspect import mesh_quality

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "artifacts")
os.makedirs(OUT, exist_ok=True)


def quality_of(mesh, tag):
    from radia.gmsh_post_export import GmshPostExport
    msh = os.path.join(OUT, f"si_{tag}.msh")
    GmshPostExport(mesh).write(msh)
    q = mesh_quality(msh)
    bts = q["by_type"]
    n = sum(bt["n_elements"] for bt in bts)
    return {"n_elements": n,
            "min_quality": round(min(bt["min_quality"] for bt in bts), 4),
            "mean_quality": round(sum(bt["mean_quality"] * bt["n_elements"]
                                      for bt in bts) / n, 4),
            "negative": sum(bt["negative"] for bt in bts)}


def solve_poisson(mesh, k):
    """Manufactured-solution Poisson; CG iteration counts + errors."""
    from ngsolve import (BND, BilinearForm, CGSolver, CoefficientFunction,
                         GridFunction, H1, Integrate, InnerProduct, LinearForm,
                         Preconditioner, Projector, cos, dx, exp, grad, sin,
                         sqrt, x, y, z)

    uex = sin(k * x) * sin(k * y) * exp(k * z)
    guex = CoefficientFunction((
        k * cos(k * x) * sin(k * y) * exp(k * z),
        k * sin(k * x) * cos(k * y) * exp(k * z),
        k * sin(k * x) * sin(k * y) * exp(k * z)))
    rhs = k * k * uex                       # -lap u_ex = k^2 u_ex

    fes = H1(mesh, order=1, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    f = LinearForm(fes)
    f += rhs * v * dx
    jac = Preconditioner(a, "local")        # Jacobi, honours FreeDofs
    a.Assemble()
    f.Assemble()

    out = {"ndof": fes.ndof, "free_dofs": sum(fes.FreeDofs())}
    for label, pre in (("cg_plain", Projector(fes.FreeDofs(), True)),
                       ("cg_jacobi", jac.mat)):
        gfu = GridFunction(fes)
        gfu.Set(uex, BND)                   # inhomogeneous Dirichlet lift
        res = f.vec.CreateVector()
        res.data = f.vec - a.mat * gfu.vec
        inv = CGSolver(a.mat, pre, precision=1e-10, maxsteps=100000,
                       printrates=False)
        gfu.vec.data += inv * res
        out[label + "_iters"] = inv.GetSteps()
        if label == "cg_jacobi":
            l2 = sqrt(Integrate((gfu - uex) ** 2, mesh, order=4))
            h1 = sqrt(Integrate(
                InnerProduct(grad(gfu) - guex, grad(gfu) - guex),
                mesh, order=4))
            out["l2_error"] = l2
            out["h1_error"] = h1
            out["l2_rel"] = l2 / sqrt(Integrate(uex ** 2, mesh, order=4))
            out["h1_rel"] = h1 / sqrt(Integrate(InnerProduct(guex, guex),
                                                mesh, order=4))
    return out


def case_cubit(step, size, tag, k):
    from ngsolve import Mesh
    vol = os.path.join(OUT, f"si_{tag}.vol")
    r = _run_batch(step, ["volume all scheme tetmesh",
                          f"volume all size {size}", "mesh volume all",
                          "block 1 add volume all", 'block 1 name "mesh"',
                          f'export netgen "{vol.replace(os.sep, "/")}" order 1'],
                   timeout_s=600)
    assert r["status"] == "ok", r
    mesh = Mesh(vol)
    row = {"mesher": "cubit_tet", "h": size}
    row.update(quality_of(mesh, tag))
    row.update(solve_poisson(mesh, k))
    return row


def case_netgen(step, maxh, tag, k):
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    mesh = Mesh(OCCGeometry(step).GenerateMesh(maxh=maxh))
    row = {"mesher": "netgen", "h": maxh}
    row.update(quality_of(mesh, tag))
    row.update(solve_poisson(mesh, k))
    return row


def main():
    t0 = time.time()
    cc = os.path.join(OUT, "c_core.step")
    sp = os.path.join(OUT, "sphere.step")
    export_step(c_core(width=80, height=60, depth=25, leg=15, gap=8), cc)
    export_step(Sphere(1.0), sp)

    from ngsolve import TaskManager
    results = {"timestamp": datetime.now().isoformat(),
               "hostname": platform.node(),
               "problem": ("H1 order-1 Poisson, manufactured "
                           "u=sin(kx)sin(ky)exp(kz), Dirichlet on all "
                           "boundaries; CG precision 1e-10, unpreconditioned "
                           "(Projector) and Jacobi"),
               "referee": "gmsh minSICN",
               "cases": []}

    # k chosen so k*L ~ 1 (sphere L=2, c_core L=80)
    pairs = [
        ("sphere_equal_budget", 1.0,
         [("cubit", sp, 0.3), ("netgen", sp, 0.2006)]),
        ("c_core_equal_budget", 0.02,
         [("cubit", cc, 4.0), ("netgen", cc, 3.1024)]),
        ("c_core_netgen_valley", 0.02,
         [("netgen", cc, 8.0435), ("netgen", cc, 6.4130),
          ("netgen", cc, 6.0870)]),
    ]

    with TaskManager():
        for pair_name, k, legs in pairs:
            entry = {"pair": pair_name, "k": k, "legs": []}
            for i, (kind, step, size) in enumerate(legs):
                tag = f"{pair_name}_{i}_{kind}"
                row = (case_cubit(step, size, tag, k) if kind == "cubit"
                       else case_netgen(step, size, tag, k))
                entry["legs"].append(row)
                print(f"[{pair_name}] {row['mesher']:9s} "
                      f"h={row['h']:7.4f} n={row['n_elements']:6d} "
                      f"min={row['min_quality']:.3f} "
                      f"ndof={row['ndof']:6d} "
                      f"CG={row['cg_plain_iters']:6d} "
                      f"CG+Jac={row['cg_jacobi_iters']:6d} "
                      f"L2rel={row['l2_rel']*100:6.2f}% "
                      f"H1rel={row['h1_rel']*100:6.2f}%", flush=True)
            results["cases"].append(entry)

    out = os.path.join(_HERE, "results_solver_impact.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"saved {out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
