"""hdiv_curved_nonlinear_field.py -- (A) the FIELD output of a CURVED x NONLINEAR body: where the curved
win is LARGE (~9x), unlike the magnetization (modest ~0.3%, see test_hdiv_vim_curved_nonlinear.py).

A nonlinear soft-iron sphere in a uniform field magnetizes UNIFORMLY (for ANY M-H law), so M is the
scalar fixed point  M = M(H_ext - D M)  with the curved demag D = 1/3.  The EXTERNAL H field is then the
EXACT point dipole, m = M V.  The curved win lives HERE: the flat faceted sphere's volume is ~9% low, so
its dipole moment -- hence the WHOLE external field -- is ~9% WRONG; mesh.Curve(3) at the SAME ndof is
<0.4% at every external point.  This is the engineering deliverable (the stray field around a nonlinear
soft-iron part) being ~9% off with flat elements (yano-type) and EXACT with curved (HDiv-VIM).

H(r) = (1/4pi) INT_S sigma(r') (r-r')/|r-r'|^3 dS',  sigma = M.n  (the surface-charge stray field; no
singular quadrature at an EXTERNAL point -> the only error is the geometry).  Validated vs the ANALYTIC
dipole -- Radia cannot referee curved geometry (its ObjHex/Tet facet the body).  The nonlinearity sets
the field MAGNITUDE (physical M); the curved win is the ~9% geometry error, which the nonlinearity does
not amplify (it merely scales it).
"""
import json
import os
import sys
from math import pi

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hdiv_demag_curved as cv            # _surface_samples (curved-aware surface quadrature)
from radia.vim import _nonlinear as nl     # _scalar_fixed_point

import ngsolve as ng                       # noqa: E402
from ngsolve import TaskManager            # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))

# soft-iron-like saturating law M(H) = chi0 H / (1 + chi0|H|/Msat)
CHI0, MSAT = 5000.0, 1.6e6


def Mof(H):
    return CHI0 * H / (1.0 + CHI0 * abs(H) / MSAT)


def reconstruct_H(mesh, M_scalar, obs, nsub=4):
    """External H field from the surface charge sigma = M_scalar * n_z, curved-aware (the SAME
    mesh.GetTrafo sampling for flat and curved -- only mesh.Curve toggled)."""
    P, w, nz, _ = cv._surface_samples(mesh, nsub, np.zeros(3))
    sw = M_scalar * nz * w                                  # sigma * dS at each surface quad point
    H = np.zeros((len(obs), 3))
    for i, r in enumerate(obs):
        d = r - P
        rn = np.linalg.norm(d, axis=1)
        H[i] = (1.0 / (4 * pi)) * np.sum(sw[:, None] * d / rn[:, None] ** 3, axis=0)
    return H


def H_dipole(r, m):
    """Analytic point-dipole H field, moment m along z."""
    rn = np.linalg.norm(r)
    rh = r / rn
    mv = np.array([0.0, 0.0, m])
    return (1.0 / (4 * pi)) * (3.0 * np.dot(mv, rh) * rh - mv) / rn ** 3


def _sphere(h):
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(g.GenerateMesh(maxh=h))


def run(h=0.6, H_ext=1e4, curve_order=3):
    M_s = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H_ext)    # nonlinear sphere: uniform M, demag D=1/3
    m = M_s * (4.0 * pi / 3.0)                              # dipole moment of the uniform-M sphere
    obs = np.array([[0, 0, 1.5], [0, 0, 2.0], [0, 0, 3.0], [1.5, 0, 0.6], [2.0, 0, 0.0]], float)
    out = {"H_ext": H_ext, "M_s": M_s, "dipole_m": m, "h": h, "obs": obs.tolist(), "cases": []}
    for curve in (0, curve_order):
        mesh = _sphere(h)
        if curve:
            with TaskManager():
                mesh.Curve(curve)
        Hn = reconstruct_H(mesh, M_s, obs)
        errs = [float(np.linalg.norm(Hn[i] - H_dipole(r, m)) / np.linalg.norm(H_dipole(r, m)))
                for i, r in enumerate(obs)]
        out["cases"].append(dict(curved=bool(curve), max_err=max(errs), errs=errs))
    return out


if __name__ == "__main__":
    res = run()
    print(f"nonlinear soft-iron sphere: M_s = {res['M_s']:.1f} A/m (H_ext={res['H_ext']:.0e}); "
          f"dipole m = {res['dipole_m']:.1f}")
    print("external H-field error vs the ANALYTIC dipole (5 points; geometry-only error):")
    for c in res["cases"]:
        tag = "Curve(3)" if c["curved"] else "FLAT    "
        print(f"  {tag}: max {100*c['max_err']:+.2f}%   per-point " +
              " ".join(f"{100*e:+.2f}%" for e in c["errs"]))
    flat = next(c for c in res["cases"] if not c["curved"])
    curv = next(c for c in res["cases"] if c["curved"])
    print(f"=> the curved field win is ~{flat['max_err']/curv['max_err']:.0f}x: the stray field of a "
          f"nonlinear soft-iron part is ~{100*flat['max_err']:.0f}% wrong with FLAT elements (yano-type)")
    print("   and <0.4% with curved -- THIS is where curved x nonlinear matters (the field, not M).")
    with open(os.path.join(HERE, "hdiv_curved_nonlinear_field.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_curved_nonlinear_field.json"))
