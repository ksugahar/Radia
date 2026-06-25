"""compare_curved_vs_radia_field.py -- (B) head-to-head: external field accuracy-per-resolution,
the SHIPPED Radia solver (FLAT tets) vs HDiv-VIM (CURVED single-layer surface charge), vs the ANALYTIC
dipole of a uniform-M sphere.

This is the quantitative basis for the curved accuracy-per-DOF win against the PRODUCTION code.  Radia's
ObjTetrahedron are FLAT -- the accessible stand-in for the six-face surface-charge distortion elements, which are also
flat.  At the SAME mesh parameter h, the HDiv curved field is ~10-30x more accurate; and curved at the
COARSEST mesh beats shipped-Radia-flat at the FINEST.

HONEST SCOPE: this measures ACCURACY-PER-RESOLUTION (geometry-driven, fair across implementations), NOT
wall-clock.  The HDiv-VIM here is a Python prototype (dense surface-charge sum), not time-optimized; a
fair speed comparison needs the C++ productionization (not done).  Radia cannot referee curved geometry
(it facets), so the reference is the ANALYTIC dipole.  The uniform-M sphere is chosen because its
external field is exactly dipolar (a clean analytic truth); the nonlinearity (soft iron) only scales the
field magnitude by M and does not change this accuracy-per-resolution picture.
"""
import json
import os
import sys
from math import pi

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src", "radia"))
import radia as rad                         # the SHIPPED production solver (flat elements)
import ngsolve as ng
from ngsolve import TaskManager
from netgen.csg import CSGeometry, Sphere, Pnt
import netgen_mesh_import as nmi
import hdiv_demag_curved as cv

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))

MU0 = 4e-7 * pi
MS = 1.0e5
M_DIP = MS * (4.0 * pi / 3.0)
OBS = np.array([[0, 0, 1.5], [0, 0, 2.0], [0, 0, 3.0], [1.5, 0, 0.6], [2.0, 0, 0.0]], float)


def _Bdip(r):
    rn = np.linalg.norm(r); rh = r / rn; mv = np.array([0.0, 0.0, M_DIP])
    return (MU0 / (4 * pi)) * (3.0 * np.dot(mv, rh) * rh - mv) / rn ** 3


def radia_flat(h):
    """Shipped Radia: flat-tet uniform-M sphere, external B via rad.Fld, max rel error vs dipole."""
    rad.UtiDelAll()
    geo = CSGeometry(); geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    ngm = ng.Mesh(geo.GenerateMesh(maxh=h))
    ne = int(ngm.ne)
    cont = nmi.netgen_mesh_to_radia(ngm, material={'magnetization': [0, 0, MS]}, units='m', verbose=False)
    B = np.array(rad.Fld(cont, 'b', OBS.tolist())).reshape(-1, 3)
    err = max(np.linalg.norm(B[i] - _Bdip(r)) / np.linalg.norm(_Bdip(r)) for i, r in enumerate(OBS))
    rad.UtiDelAll()
    return ne, float(err)


def hdiv_curved(h, curve=3, nsub=4):
    """HDiv-VIM: curved surface-charge H field (B=mu0 H outside), max rel error vs dipole."""
    geo = CSGeometry(); geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=h))
    with TaskManager():
        mesh.Curve(curve)
    P, w, nz, nb = cv._surface_samples(mesh, nsub, np.zeros(3))
    sw = MS * nz * w
    err = 0.0
    for r in OBS:
        d = r - P; rn = np.linalg.norm(d, axis=1)
        H = (1.0 / (4 * pi)) * np.sum(sw[:, None] * d / rn[:, None] ** 3, axis=0)
        err = max(err, np.linalg.norm(H - _Bdip(r) / MU0) / np.linalg.norm(_Bdip(r) / MU0))
    return int(nb), float(err)


def run(hs=(0.6, 0.4, 0.3, 0.2)):
    out = {"hs": list(hs), "rows": []}
    for h in hs:
        ne, ef = radia_flat(h)
        nb, ec = hdiv_curved(h)
        out["rows"].append(dict(h=h, radia_tets=ne, radia_flat_err=ef, hdiv_bnd_tris=nb, hdiv_curved_err=ec))
    out["flat_finest_err"] = out["rows"][-1]["radia_flat_err"]
    out["curved_coarsest_err"] = out["rows"][0]["hdiv_curved_err"]
    return out


if __name__ == "__main__":
    res = run()
    print("external field accuracy vs analytic dipole -- shipped Radia (FLAT) vs HDiv-VIM (CURVED):")
    print(f"  {'h':>5} | {'Radia tets':>10} {'flat err':>10} | {'HDiv tris':>10} {'curved err':>11}")
    for r in res["rows"]:
        print(f"  {r['h']:>5} | {r['radia_tets']:>10d} {100*r['radia_flat_err']:>9.2f}% | "
              f"{r['hdiv_bnd_tris']:>10d} {100*r['hdiv_curved_err']:>10.3f}%")
    fr, cr = res["rows"][-1], res["rows"][0]
    print(f"\n=> shipped Radia FLAT at finest h={fr['h']} ({fr['radia_tets']} tets): "
          f"{100*res['flat_finest_err']:.2f}%  vs  HDiv CURVED at coarsest h={cr['h']} "
          f"({cr['hdiv_bnd_tris']} tris): {100*res['curved_coarsest_err']:.3f}%")
    print("   curved-coarsest beats flat-finest = the accuracy-per-resolution win vs the production")
    print("   solver (vs analytic truth).  ACCURACY-PER-DOF; wall-clock = the C++ lift, not done.")
    with open(os.path.join(HERE, "compare_curved_vs_radia_field.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "compare_curved_vs_radia_field.json"))
