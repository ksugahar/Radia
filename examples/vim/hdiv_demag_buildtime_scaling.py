"""hdiv_demag_buildtime_scaling.py -- HDiv-VIM analytic charge-Gram H-matrix BUILD time vs N.

The wall-clock head-to-head vs yano-type (Radia MMM/MSC) on the C-type electromagnet has two phases:
the one-time H-matrix BUILD and the iterative SOLVE.  The SOLVE is the clear HDiv-VIM win (5-6 Newton
iters vs yano's 214 -- loops field-null by construction).  This benchmark measures the BUILD, the
dominant remaining unknown ("is the H-matrix build slow?"; the saved yano-type C-type build = 582 s @
165600 DOF, examples/c_type_electromagnet/.../hacapk/165600DOF.json).

KELVIN-LESS by construction: the HDiv-VIM is a VOLUME INTEGRAL method (like MMM/MSC) -- the 1/r charge
Gram handles the OPEN boundary analytically, so only the IRON is meshed (here a cube as the iron proxy):
NO air domain, NO Kelvin transformation (that is a FEM-only need).

RESULT (cube, eps=1e-4, leaf=32; charges extracted straight from the tet mesh, no dense build_demag):
  n_charge   t_Hbuild   compr
       281       0.19s    1.00
       901       1.18s    0.88
      1541       3.80s    0.70
      4225      13.2s     0.28
      7278      23.9s     0.21

The build scales ~N^1.1-1.3 at large N (ACA compression kicking in: compr 1.0 -> 0.21).  Extrapolated
to C-type scale (n_charge ~30k-150k) -> ~150-1200 s, the SAME ORDER as yano's 582 s build (comparable,
not a clear win).  CAUSE: the entry is ALWAYS-analytic (every pair uses the expensive PhiTet/TriPotential
~40-80 calls).  LEVER: a near/far split (cheap centroid-monopole for ACA-far pairs, analytic only for
near) would make the build a clear win.  Matvec is fast (8.5 ms @ 7278, O(N log N)) -> the solve phase
(~5-6 Newton iters) is ~tens of seconds vs yano's 1953 s.  Total estimate: HDiv ~2-6x faster, BUILD-LIMITED.
"""
import json
import os
import time

import numpy as np
import ngsolve as ng
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import radia._radia_pybind as _rp

HERE = os.path.dirname(os.path.abspath(__file__))


def _cube(h):
    g = CSGeometry(); g.Add(OrthoBrick(Pnt(-.5, -.5, -.5), Pnt(.5, .5, .5)))
    return ng.Mesh(g.GenerateMesh(maxh=h))


def _charges(mesh):
    el_V = np.array([[list(mesh[v].point) for v in el.vertices] for el in mesh.Elements(ng.VOL)])
    bf_V = np.array([[list(mesh[v].point) for v in el.vertices] for el in mesh.Elements(ng.BND)])
    return el_V.ravel(), bf_V.ravel(), len(el_V), len(bf_V)


def run(hs=(0.30, 0.22, 0.16, 0.12, 0.09)):
    rows = []
    for h in hs:
        mesh = _cube(h)
        cv, fv, n_el, n_bf = _charges(mesh)
        n = n_el + n_bf
        with ng.TaskManager():
            t0 = time.perf_counter()
            H = _rp._ChargeGramHMatrix(cell_verts=list(cv), face_verts=list(fv), n_el=n_el,
                                       eps=1e-4, leaf=32, eta=2.0)
            t_build = time.perf_counter() - t0
        st = H.stats()
        rng = np.random.default_rng(0); q = rng.standard_normal(n)
        t0 = time.perf_counter(); H.matvec(q.tolist()); t_mv = time.perf_counter() - t0
        rows.append(dict(maxh=h, n_el=int(n_el), n_charge=int(n), t_Hbuild_s=float(t_build),
                         compression=float(st["compression"]), memory_mb=float(st["memory_mb"]),
                         matvec_ms=float(t_mv * 1e3)))
        print(f"  n_charge={n:>6}  t_Hbuild={t_build:>7.2f}s  compr={st['compression']:.3f}  mv={t_mv*1e3:.1f}ms")
    return rows


if __name__ == "__main__":
    rows = run()
    with open(os.path.join(HERE, "hdiv_demag_buildtime_scaling.json"), "w") as f:
        json.dump({"kernel": "analytic charge-Gram (HDiv-VIM, Kelvin-less)", "body": "cube (iron proxy)",
                   "eps": 1e-4, "yano_ctype_build_s": 582.0, "yano_ctype_dof": 165600, "rows": rows}, f, indent=2)
    print("saved hdiv_demag_buildtime_scaling.json")
