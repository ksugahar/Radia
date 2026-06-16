"""HDiv-VIM demag on HEX and WEDGE meshes via the dense analytic POLYTOPE charge Gram.

The tet/triangle analytic charge Gram (phi_tet / tri_potential, the exact Wilton potentials) generalizes
to ANY flat-faced convex cell + quad face with no new singular quadrature: a cell's Newtonian potential
is the divergence-theorem sum over its convex-hull triangular faces of the SAME Wilton triangle
potential, and a quad face is two flat triangles (radia.vim._core: _polytope_potential / _cell_hull_tris
/ _face_subtris).  So `radia.vim.hdiv_demag_solve` solves hex AND wedge soft iron, not just tet.

This demonstrates the physics gate: a CUBE meshed three ways -- hexes, wedges (prisms), tets -- has the
same demag factor (~1/3, isotropic uniform-M body) and the same solved volume-averaged magnetization in
a uniform applied field, to <1% across element types.  The scalable C++ charge-Gram H-matrix is
tetrahedron-only; hex/wedge take the dense analytic path here (correct, O(N^2)).

Run:  python hdiv_demag_hex_wedge.py   (writes hdiv_demag_hex_wedge.json next to this script)
"""
import json
import math
import os

import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from netgen.occ import Box, OCCGeometry

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))
import radia.vim as vim  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
H0 = 1000.0          # applied H (A/m)
L = 0.02             # cube edge (m)
MU_R = 100.0


def _solve(mesh, tag):
    with ng.TaskManager():
        res = vim.hdiv_demag_solve(mesh, mu_r=MU_R, H_ext=ng.CoefficientFunction((0, 0, H0)))
    out = dict(tag=tag, n_el=int(res["n_el"]), ndof=int(res["ndof"]), iters=int(res["iters"]),
               demag=float(res["demag"]), M_avg_z=float(res["M_avg"][2]))
    print(f"[{tag}] n_el={out['n_el']:5d} ndof={out['ndof']:5d} iters={out['iters']:3d} "
          f"demag={out['demag']:.4f} M_avg_z={out['M_avg_z']:.1f}")
    return out


def main():
    mp = lambda x, y, z: (L * x, L * y, L * z)  # noqa: E731
    results = []
    results.append(_solve(MakeStructured3DMesh(hexes=True, nx=5, ny=5, nz=5, mapping=mp), "hex 5^3"))
    results.append(_solve(MakeStructured3DMesh(hexes=False, prism=True, nx=5, ny=5, nz=5, mapping=mp),
                          "wedge 5^3"))
    with ng.TaskManager():
        tetm = ng.Mesh(OCCGeometry(Box((0, 0, 0), (L, L, L))).GenerateMesh(maxh=L / 5))
    results.append(_solve(tetm, "tet h=L/5"))

    mz = {r["tag"]: r["M_avg_z"] for r in results}
    ref = mz["tet h=L/5"]
    rel = {t: abs(mz[t] - ref) / abs(ref) for t in mz}
    print("element-type agreement of M_avg_z vs tet:")
    for t in ("hex 5^3", "wedge 5^3"):
        print(f"  {t}: rel = {rel[t]:.2e}")

    data = {
        "description": "HDiv-VIM demag on a cube meshed as hex / wedge / tet (polytope analytic Gram); "
                       "demag factor ~1/3 and M_avg element-type independent.",
        "applied_H_Am": H0, "cube_edge_m": L, "mu_r": MU_R,
        "results": results,
        "M_avg_z_rel_vs_tet": rel,
    }
    path = os.path.join(HERE, "hdiv_demag_hex_wedge.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print("saved", path)


if __name__ == "__main__":
    main()
