"""hdiv_demag_tet_gram_study.py -- resolve the TRUE sphere demag factor of the HDiv-type VIM (the
verify-first foundation for production step #3, the Wilton/exact Gram).

FINDING that motivates this study: for a UNIFORM magnetization M_z, div M = 0, so the RT0 volume
(cell) charges are ~0 and the demag is ENTIRELY a BOUNDARY FACE-FACE (surface charge sigma = M.n)
Gram phenomenon.  The shipped Gram uses CENTROID-MONOPOLE for off-diagonal entries, which is
INACCURATE for adjacent (near) boundary faces -- and the demag VALUE depends on it.

This study replaces the boundary face-face off-diagonal block with SUB-POINT quadrature at increasing
nsub (the diagonal self-energy already uses sub-points).  nsub=1 (one centroid sub-point per face)
reduces EXACTLY to the centroid-monopole -> a built-in consistency check; nsub->inf converges to the
exact Coulomb Gram.  So demag(nsub) shows the Gram-convergence at fixed mesh, and demag(h) the mesh
(discretization) convergence.  Together they answer: what is the TRUE (Gram-exact) sphere demag, and
does it converge to the analytic 1/3?

Run: python hdiv_demag_tet_gram_study.py  (writes hdiv_demag_tet_gram_study.json).
"""
import json
import os
from math import pi

import numpy as np

import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt

from radia.hdiv_vim import _core as tet   # build_demag, demag_factor, _bary_tri

HERE = os.path.dirname(os.path.abspath(__file__))
INV4PI = 1.0 / (4.0 * pi)


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():                               # caller wraps (build_demag is a helper, no internal TM)
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def _tri_area(V):
    return 0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0]))


def demag_facegram_nsub(mesh, nsub):
    """Demag factor with the boundary face-face off-diagonal Gram block computed by nsub-sub-point
    quadrature (diagonal self-energy unchanged).  nsub=1 == centroid-monopole."""
    d = tet.build_demag(mesh, nsub=4)            # baseline G (monopole off-diag + sub-point self)
    bf_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.BND)]
    n_bf = len(bf_V)
    n_cell = d["n_charge"] - n_bf
    G = d["G"].copy()
    # boundary-face sub-points + weights (barycentric); nsub=1 -> single centroid (== monopole)
    lam = tet._bary_tri(nsub)
    SP = [lam @ V for V in bf_V]
    SW = [np.full(len(lam), _tri_area(V) / len(lam)) for V in bf_V]
    # overwrite the boundary face-face OFF-diagonal block with the sub-point Coulomb sum
    for a in range(n_bf):
        ca = n_cell + a
        for b in range(a + 1, n_bf):
            cb = n_cell + b
            D = np.linalg.norm(SP[a][:, None, :] - SP[b][None, :, :], axis=2)
            g = float(np.sum(np.outer(SW[a], SW[b]) * INV4PI / D))
            G[ca, cb] = g
            G[cb, ca] = g
    N = d["B"].T @ G @ d["B"]
    m = d["m_unit"]
    return float((m @ N @ m) / (m @ d["M_mass"] @ m)), n_bf


def run():
    out = {"note": "sphere demag (analytic 1/3=0.3333); nsub=1 == centroid-monopole baseline",
           "rows": []}
    print("sphere demag factor (analytic = 1/3 = 0.33333)")
    print(f"{'h':>6} {'n_bf':>5} | " + " ".join(f"nsub={n:<2d}" for n in (1, 2, 4, 8)))
    for h in (0.5, 0.35, 0.25):
        mesh = _sphere(h)
        vals = {}
        n_bf = 0
        for nsub in (1, 2, 4, 8):
            Dz, n_bf = demag_facegram_nsub(mesh, nsub)
            vals[nsub] = Dz
        print(f"{h:6.2f} {n_bf:5d} | " + " ".join(f"{vals[n]:.4f} " for n in (1, 2, 4, 8)))
        out["rows"].append({"h": h, "n_bf": n_bf, "demag": {str(n): vals[n] for n in (1, 2, 4, 8)}})
    with open(os.path.join(HERE, "hdiv_demag_tet_gram_study.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_demag_tet_gram_study.json"))
    return out


if __name__ == "__main__":
    run()
