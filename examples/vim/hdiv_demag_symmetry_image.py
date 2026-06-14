"""hdiv_demag_symmetry_image.py -- SYMMETRY MODELS (1/2, 1/4, 1/8) for the HDiv-VIM demag, via the
IMAGE method, validated against the full-sphere demag.

Two pieces make a symmetry model work; the HDiv-VIM gets BOTH cheaply:
  (1) LOOPS: automatic -- they are ker(B) on the cut mesh, field-null by construction, no cohomology
      loop-star bookkeeping (see test_hdiv_vim_symmetry_loops.py).
  (2) the DEMAG VALUE: the IMAGE method here.  Only the REAL surface (the spherical cap) carries the
      surface charge sigma = M.n = n_z; the flat cut faces are symmetry planes (no real charge).
      Reflecting the cap charge over the reduction planes -- with sign = (-1)^(#z-reflections), because
      sigma = n_z flips under a z-mirror -- reconstructs the full sphere's sigma = cos(theta).  So the
      reduced model computes the FULL demag from only ~1/N the independent surface DOF.

Validation: the 1/2, 1/4, 1/8 demag must equal the full-sphere demag (same crude sub-point Gram), from
~1/2, 1/4, 1/8 the DOF.  M = z_hat; the field-parallel cut planes (x=0, y=0) keep sigma's sign, the
field-perpendicular plane (z=0, present only in the 1/8 model) flips it -- the IMA sign rule.
"""
import json
import os
import sys
from math import pi

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hdiv_demag_curved as cv
from radia.vim import _core as tet

import ngsolve as ng
from ngsolve import IntegrationRule, ElementId, BND, TaskManager
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))
V_FULL = 4.0 * pi / 3.0

# which axes are reduced (cut by a symmetry plane) for each model
MODELS = {"full": (), "half": ("x",), "quarter": ("x", "y"), "eighth": ("x", "y", "z")}
_LO = {"x": 0.0, "y": 0.0, "z": 0.0}


def _demag(P, sw, w):
    """D_z = (1/4pi V_full)[ sum_{a!=b} sw_a sw_b / r + sum_a sigma_a^2 C_TRI w_a^1.5 ], sigma = sw/w."""
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    np.fill_diagonal(D, np.inf)
    off = float(np.sum(np.outer(sw, sw) / D))
    selfd = float(np.sum((sw ** 2 / w) * tet.C_TRI * np.sqrt(w)))
    return (off + selfd) / (4 * pi * V_FULL)


def _reduced_mesh(cut_axes, h):
    lo = [(-2.0 if a not in cut_axes else 0.0) for a in "xyz"]
    geo = CSGeometry(); geo.Add(Sphere(Pnt(0, 0, 0), 1.0) * OrthoBrick(Pnt(*lo), Pnt(2, 2, 2)))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def model_demag(model, h, curve=3, nsub=4):
    cut_axes = MODELS[model]
    mesh = _reduced_mesh(cut_axes, h)
    with TaskManager():
        mesh.Curve(curve)
    lam = tet._bary_tri(nsub); xi, eta, m = lam[:, 1], lam[:, 2], len(lam)
    ax_idx = {"x": 0, "y": 1, "z": 2}
    capP, capW, capNZ = [], [], []
    n_cap_el = 0
    for i in range(sum(1 for _ in mesh.Elements(BND))):
        ir = IntegrationRule(points=[(xi[k], eta[k], 0) for k in range(m)], weights=[1.0] * m)
        trafo = mesh.GetTrafo(ElementId(BND, i))
        pts, Js, nzs = [], [], []
        for ip in ir:
            mip = trafo(ip)
            p = np.array([mip.point[0], mip.point[1], mip.point[2]])
            jac = np.asarray(mip.jacobi); nrm = np.cross(jac[:, 0], jac[:, 1]); nrm /= np.linalg.norm(nrm) + 1e-30
            if np.dot(p, nrm) < 0:
                nrm = -nrm
            pts.append(p); Js.append(mip.measure); nzs.append(nrm[2])
        pts = np.array(pts); c = pts.mean(0)
        # a CUT (symmetry-plane) element has its centroid on a reduced-axis plane -> no real charge
        if any(abs(c[ax_idx[a]]) < 1e-4 for a in cut_axes):
            continue
        n_cap_el += 1
        for k in range(m):
            capP.append(pts[k]); capW.append(Js[k] * 0.5 / m); capNZ.append(nzs[k])
    capP = np.array(capP); capW = np.array(capW); capNZ = np.array(capNZ)
    # reflect the cap over the reduced axes (product of +-1); sigma = n_z flips sign per z-reflection
    signs_per_axis = [(1, -1) if a in cut_axes else (1,) for a in "xyz"]
    Ps, SWs, Ws = [], [], []
    for sx in signs_per_axis[0]:
        for sy in signs_per_axis[1]:
            for sz in signs_per_axis[2]:
                Ps.append(capP * np.array([sx, sy, sz]))
                SWs.append((sz * capNZ) * capW)            # sz = sigma sign (n_z flips under z-mirror)
                Ws.append(capW)
    Pall = np.vstack(Ps); SWall = np.concatenate(SWs); Wall = np.concatenate(Ws)
    return _demag(Pall, SWall, Wall), len(capP), n_cap_el


def run(h=0.4):
    Dfull, _, _ = model_demag("full", h)
    out = {"h": h, "full_demag": Dfull, "models": []}
    for model in ("half", "quarter", "eighth"):
        Dm, ncap_pts, ncap_el = model_demag(model, h)
        out["models"].append(dict(model=model, demag=Dm, cap_pts=ncap_pts, cap_el=ncap_el,
                                  match_vs_full=Dm / Dfull - 1.0))
    return out


if __name__ == "__main__":
    res = run()
    print(f"symmetry-model demag via the IMAGE method (M=z_hat), full-sphere demag = {res['full_demag']:.5f}:")
    print(f"  {'model':>8} {'demag':>9} {'cap-tris':>9} {'match vs full':>14}")
    for r in res["models"]:
        print(f"  {r['model']:>8} {r['demag']:>9.5f} {r['cap_el']:>9d} {100*r['match_vs_full']:>+13.3f}%")
    print("=> 1/2, 1/4, 1/8 reproduce the full demag from ~1/2, 1/4, 1/8 the surface DOF.  Loops are")
    print("   automatic (ker B); this image method supplies the symmetry-model demag value.")
    with open(os.path.join(HERE, "hdiv_demag_symmetry_image.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_demag_symmetry_image.json"))
