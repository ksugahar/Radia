"""Control: is hex's huge L2 win on the thin plate just AXIS ALIGNMENT?

Promoted from C:/temp/mesh_quality_study (2026-08-07) with its committed
results JSON (Data Persistence Policy). Re-run with
`python run_rotation_control.py` (requires Cubit + netgen + build123d).

vector_elements.py found cubit_hex beating cubit_tet by ~30x in the HCurl
L2 error on an axis-aligned plate -- while the curl seminorm showed only
~1.4x and the HDiv div seminorm showed hex *behind*. That asymmetry is
the signature of a superconvergence artifact: the manufactured field is
separable and axis-aligned, and so is the hex lattice.

Rigid rotation commutes with curl and div, so for a proper rotation R

    u'(x) = R u(R^T x)   =>   curl curl u' = k^2 u'   (HCurl case)
                              -grad div u' = k^2 u'   (HDiv case)

still holds exactly. Rotating the FIELD (not the mesh) therefore breaks
the alignment while keeping the same geometry, the same meshes, and the
same exact-solution machinery -- so any change in the hex/tet ratio is
attributable to alignment alone.
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

from build123d import Box, export_step
from radia_mcp.cubit.server import _run_batch

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "artifacts")
os.makedirs(OUT, exist_ok=True)
MAXSTEPS = 30000

# proper rotation: 30 deg about z, then 20 deg about y (no reflection)
_cz, _sz = math.cos(math.radians(30)), math.sin(math.radians(30))
_cy, _sy = math.cos(math.radians(20)), math.sin(math.radians(20))
RZ = [[_cz, -_sz, 0.0], [_sz, _cz, 0.0], [0.0, 0.0, 1.0]]
RY = [[_cy, 0.0, _sy], [0.0, 1.0, 0.0], [-_sy, 0.0, _cy]]
R = [[sum(RY[i][k] * RZ[k][j] for k in range(3)) for j in range(3)]
     for i in range(3)]


def _rot_apply(R, vec):
    from ngsolve import CoefficientFunction
    return CoefficientFunction(tuple(
        sum(R[i][j] * vec[j] for j in range(3)) for i in range(3)))


def _rt_coords(R):
    """(R^T x) as three scalar CFs."""
    from ngsolve import x, y, z
    X = (x, y, z)
    return [sum(R[j][i] * X[j] for j in range(3)) for i in range(3)]


def fields(space, k, rotate):
    """Return (u_ex, deriv_ex) for the given space, optionally rotated."""
    from ngsolve import CoefficientFunction, cos, sin, x, y, z
    Rm = R if rotate else [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    X, Y, Z = _rt_coords(Rm)
    if space == "hcurl":
        u = [sin(k * Y), sin(k * Z), sin(k * X)]
        d = [-k * cos(k * Z), -k * cos(k * X), -k * cos(k * Y)]
        return _rot_apply(Rm, u), _rot_apply(Rm, d)     # curl is a vector
    u = [sin(k * X), sin(k * Y), sin(k * Z)]
    div = k * (cos(k * X) + cos(k * Y) + cos(k * Z))    # div is a scalar
    return _rot_apply(Rm, u), CoefficientFunction(div)


def solve(mesh, space, k, rotate):
    from ngsolve import (BND, BilinearForm, CGSolver, GridFunction, HCurl,
                         HDiv, Integrate, InnerProduct, LinearForm,
                         Preconditioner, curl, div, dx, sqrt)
    uex, dex = fields(space, k, rotate)
    if space == "hcurl":
        fes = HCurl(mesh, order=0, dirichlet=".*")
        deriv = curl
    else:
        fes = HDiv(mesh, order=0, dirichlet=".*")
        deriv = div
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += (InnerProduct(deriv(u), deriv(v)) + InnerProduct(u, v)) * dx
    f = LinearForm(fes)
    f += (k * k + 1.0) * InnerProduct(uex, v) * dx
    jac = Preconditioner(a, "local")
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(uex, BND)
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    r0 = res.Norm()
    inv = CGSolver(a.mat, jac.mat, precision=1e-10, maxsteps=MAXSTEPS,
                   printrates=False)
    gfu.vec.data += inv * res
    chk = f.vec.CreateVector()
    chk.data = f.vec - a.mat * gfu.vec
    for i in range(len(chk)):
        if not fes.FreeDofs()[i]:
            chk[i] = 0.0
    rel_res = chk.Norm() / max(r0, 1e-300)

    def _n(cf):
        return sqrt(Integrate(InnerProduct(cf, cf), mesh, order=4))

    return {"ndof": fes.ndof, "ne": mesh.ne, "cg_iters": inv.GetSteps(),
            "rel_residual": rel_res,
            "converged": inv.GetSteps() < MAXSTEPS and rel_res < 1e-6,
            "l2_rel": _n(gfu - uex) / _n(uex),
            "d_rel": _n(deriv(gfu) - dex) / _n(dex)}


def mesh_cubit(step, size, tag, hex_mode):
    from ngsolve import Mesh
    vol = os.path.join(OUT, f"rc_{tag}.vol")
    cmds = [] if hex_mode else ["volume all scheme tetmesh"]
    cmds += [f"volume all size {size}", "mesh volume all",
             "block 1 add volume all", 'block 1 name "mesh"',
             f'export netgen "{vol.replace(os.sep, "/")}" order 1']
    r = _run_batch(step, cmds, timeout_s=900)
    assert r["status"] == "ok", r
    return Mesh(vol)


def main():
    t0 = time.time()
    tp = os.path.join(OUT, "thin_plate.step")
    export_step(Box(8, 8, 1), tp)
    sizes = [0.5, 0.4, 0.32, 0.25]
    k = 0.4

    from ngsolve import TaskManager
    out = {"timestamp": datetime.now().isoformat(),
           "hostname": platform.node(),
           "rotation_deg": {"z": 30, "y": 20},
           "note": ("same geometry, same meshes, same exact-solution "
                    "identity; only the field's orientation changes"),
           "cases": []}

    with TaskManager():
        meshes = {r: [mesh_cubit(tp, s, f"{r}_{i}", r == "cubit_hex")
                      for i, s in enumerate(sizes)]
                  for r in ("cubit_tet", "cubit_hex")}
        for space in ("hcurl", "hdiv"):
            for rotate in (False, True):
                for route in ("cubit_tet", "cubit_hex"):
                    for m, s in zip(meshes[route], sizes):
                        row = solve(m, space, k, rotate)
                        row.update({"space": space, "rotated": rotate,
                                    "route": route, "h": s})
                        out["cases"].append(row)
                        print(f"{space:5s} rot={str(rotate):5s} {route:9s} "
                              f"h={s:5.2f} ndof={row['ndof']:6d} "
                              f"L2={row['l2_rel']*100:8.4f}% "
                              f"d={row['d_rel']*100:8.4f}% "
                              f"CG={row['cg_iters']:5d}"
                              f"{'' if row['converged'] else ' NOT-CONV'}",
                              flush=True)
                print(flush=True)

    # hex/tet ratio at each size, rotated vs not
    print("=== hex/tet ratio (same size, so nearly matched geometry) ===")
    summ = []
    for space in ("hcurl", "hdiv"):
        for rotate in (False, True):
            for s in sizes:
                def g(route):
                    return next(c for c in out["cases"]
                                if c["space"] == space and c["rotated"] == rotate
                                and c["route"] == route and c["h"] == s)
                t, h = g("cubit_tet"), g("cubit_hex")
                row = {"space": space, "rotated": rotate, "h": s,
                       "l2_hex_over_tet": h["l2_rel"] / t["l2_rel"],
                       "d_hex_over_tet": h["d_rel"] / t["d_rel"],
                       "cg_hex_over_tet": h["cg_iters"] / t["cg_iters"],
                       "ndof_hex_over_tet": h["ndof"] / t["ndof"]}
                summ.append(row)
                print(f"{space:5s} rot={str(rotate):5s} h={s:5.2f}  "
                      f"L2 {row['l2_hex_over_tet']:8.4f}x  "
                      f"deriv {row['d_hex_over_tet']:7.3f}x  "
                      f"CG {row['cg_hex_over_tet']:6.3f}x  "
                      f"(ndof {row['ndof_hex_over_tet']:.2f}x)")
    out["hex_over_tet"] = summ
    p = os.path.join(_HERE, "results_rotation_control.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"saved {p} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
