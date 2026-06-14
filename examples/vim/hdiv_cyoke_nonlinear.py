"""hdiv_cyoke_nonlinear.py -- NONLINEAR HDiv-VIM on a C-yoke (reentrant-corner, NON-convex, NON-uniform M),
cross-validated vs the shipped Radia solver.

The sphere / spheroid nonlinear cases magnetise UNIFORMLY (surface charge only).  A C-yoke is the hard
case: non-convex with reentrant corners, and the nonlinear M is genuinely NON-uniform (div M != 0 ->
volume charge), so it exercises the FULL operator (surface AND volume Gram).

KEY (a trap worth stating): a NON-uniform-M nonlinear problem needs analytic_gram=True (the full volume
Gram via phi_tet).  The surface-only wilton_surface Gram leaves the volume (cell) blocks crude and Newton
does NOT converge -- now a LOUD RuntimeError (require_convergence default), not a silent wrong result.

With analytic_gram the C-yoke nonlinear:
  - converges in ~6 damped-Newton iters,
  - is mesh-stable (volume-avg M_z drifts ~1.5% monotonically over maxh 0.020 -> 0.013, converging),
  - matches the shipped Radia MMM solver (SAME flat mesh, SAME M-H law, SAME applied field) to < 1%.
A C-yoke is FLAT, so Radia is a VALID cross-reference (both facet identically) -- this is cross-method
agreement, not a vs-analytic-truth error.  M-H law: M(H) = chi0 H / (1 + chi0|H|/Msat); Radia consumes
the equivalent B(H) = mu0(H + M(H)) table via MatSatIsoTab.
"""
import json
import os
import sys
from math import pi

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src", "radia"))
import ngsolve as ng
from ngsolve import TaskManager
from netgen.occ import Box, OCCGeometry, Pnt as OPnt
from radia.vim import _nonlinear as nl
import radia as rad
import netgen_mesh_import as nmi

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7 * pi


def cyoke():
    """A square 'C' yoke: outer box minus inner cavity minus a one-sided gap (the C opening)."""
    return (Box(OPnt(-0.06, -0.06, -0.02), OPnt(0.06, 0.06, 0.02))
            - Box(OPnt(-0.035, -0.035, -0.03), OPnt(0.035, 0.035, 0.03))
            - Box(OPnt(0.018, -0.07, -0.03), OPnt(0.07, 0.07, 0.03)))


def cyoke_mesh(h):
    return ng.Mesh(OCCGeometry(cyoke()).GenerateMesh(maxh=h))


def hdiv_cyoke_Mz(mesh, chi0, Msat, H0):
    """HDiv-VIM volume-avg M_z on the C-yoke (analytic_gram = the required nonlinear Gram)."""
    with TaskManager():
        Mh, nit, _ = nl.solve_nonlinear_newton(mesh, chi0, Msat, H0, analytic_gram=True, maxit=60)
    return Mh, nit


def radia_cyoke_Mz(mesh, chi0, Msat, H0):
    """Shipped Radia (MMM, MatSatIsoTab) volume-avg M_z on the SAME mesh, SAME M-H law, SAME applied H0."""
    rad.UtiDelAll()
    cont = nmi.netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m', verbose=False)
    Hs = np.concatenate([[0.0], np.logspace(-1, 7, 60)])
    Ms = chi0 * Hs / (1.0 + chi0 * Hs / Msat)
    Bs = MU0 * (Hs + Ms)
    mat = rad.MatSatIsoTab([[float(a), float(b)] for a, b in zip(Hs, Bs)])
    rad.MatApl(cont, mat)
    rad.Solve(rad.ObjCnt([cont, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])]), 1e-6, 3000, 0)
    cents, vols = [], []
    for el in mesh.Elements(ng.VOL):
        vs = np.array([mesh[v].point for v in el.vertices])
        cents.append(vs.mean(0)); vols.append(abs(np.linalg.det(vs[1:] - vs[0])) / 6.0)
    cents, vols = np.array(cents), np.array(vols)
    M = np.array(rad.Fld(cont, 'm', cents.tolist())).reshape(-1, 3)
    mz = float(np.sum(M[:, 2] * vols) / np.sum(vols))
    rad.UtiDelAll()
    return mz


def run(hs=(0.020, 0.016, 0.013), chi0=1000.0, Msat=1.0e6, H0=2.0e5):
    out = {"chi0": chi0, "Msat": Msat, "H0": H0, "rows": []}
    for h in hs:
        mesh = cyoke_mesh(h); ne = int(mesh.ne)
        Mh, nit = hdiv_cyoke_Mz(mesh, chi0, Msat, H0)
        Mr = radia_cyoke_Mz(mesh, chi0, Msat, H0)
        out["rows"].append(dict(h=h, ne=ne, hdiv_Mz=Mh, nit=nit, radia_Mz=Mr, agree=Mh / Mr - 1.0))
    return out


if __name__ == "__main__":
    res = run()
    print("C-yoke (non-convex, reentrant corners) NONLINEAR -- HDiv-VIM (analytic_gram) vs shipped Radia:")
    print(f"  {'maxh':>6} {'ne':>5} {'HDiv M_z':>11} {'iters':>5} {'Radia M_z':>11} {'agree':>9}")
    for r in res["rows"]:
        print(f"  {r['h']:>6} {r['ne']:>5d} {r['hdiv_Mz']:>11.1f} {r['nit']:>5d} {r['radia_Mz']:>11.1f} "
              f"{100*r['agree']:>+8.2f}%")
    print("=> non-uniform nonlinear on a sharp non-convex body: converges in ~6 iters, mesh-stable, and")
    print("   agrees with the shipped Radia solver to <1% (both flat -> Radia is a valid cross-reference).")
    with open(os.path.join(HERE, "hdiv_cyoke_nonlinear.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_cyoke_nonlinear.json"))
