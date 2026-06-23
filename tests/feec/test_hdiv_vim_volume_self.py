"""Golden: the order-0 (RT0) ANALYTIC charge-Gram VOLUME self-energy

    G[a][a] = (1/4pi) INT_T INT_T 1/|x-y| dy dx

matches the converged value of that integral.

REGRESSION GUARD for the 2026-06-12 fix.  The order-0 analytic ctor's CELL outer-quadrature was a crude
10-point equal-weight _bary_tet(3) rule that under-integrated the volume self-energy by ~6.5%.  It was
INVISIBLE to every other HDiv-VIM golden because they all use UNIFORM M -> div M = 0 -> the volume charge
(and hence the volume self) never participates; the demag is surface-charge only.  This test reads the RAW
volume-self Gram diagonal of the C++ charge-Gram H-matrix directly, so div M != 0 is exercised and the
defect cannot hide.

The reference is the SAME integral with an exact analytic inner (the Wilton uniform-tet potential, built
here from the kept tri_potential via the divergence theorem) and a FINE Gauss-Duffy outer rule -- it
self-converges, so it is an independent gold standard, not a re-statement of the production rule.
"""
import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Box, OCCGeometry, Pnt

import radia._radia_pybind as rp
from radia.vim._core import tri_potential, _gauss_duffy_tet


def _tet_vol(V):
    return abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]))) / 6.0


def _phi_tet(V, P):
    """Exact Newtonian potential INT_tet 1/|P-r'| dV' of a uniform tet via the divergence theorem
    (nabla'^2 R = 2/R -> (1/2) sum_{4 faces} d_face INT_face 1/R dA'), the inner face integral = the
    kept exact Wilton triangle potential tri_potential.  P: (M,3) -> (M,)."""
    V = np.asarray(V, float)
    P = np.atleast_2d(np.asarray(P, float))
    cen = V.mean(0)
    tot = np.zeros(len(P))
    for f in ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)):
        Fv = V[list(f)]
        n = np.cross(Fv[1] - Fv[0], Fv[2] - Fv[0])
        n = n / np.linalg.norm(n)
        if np.dot(Fv.mean(0) - cen, n) < 0:
            n = -n                                          # outward face normal
        d = (Fv[0] - P) @ n                                 # (M,) signed distance P -> face plane
        tot += d * tri_potential(Fv, P)
    return 0.5 * tot


def _gold_self(V):
    """(1/4pi) INT_T phi(V, x) dx via a FINE (o=12) Gauss-Duffy outer rule; phi is the EXACT Wilton
    uniform-tet potential (_phi_tet), so this self-converges to the true volume self-energy."""
    V = np.asarray(V, float)
    refP, refW = _gauss_duffy_tet(12)
    P = V[0] + refP @ (V[1:] - V[0])
    return float(np.sum(refW * 6.0 * _tet_vol(V) * _phi_tet(V, P)) / (4.0 * np.pi))


def test_order0_volume_self_energy_matches_gold():
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.6))
        vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
        bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
        el_V = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in vels]
        bf_V = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in bels]

        G = rp._ChargeGramHMatrix(
            cell_verts=np.concatenate([V.ravel() for V in el_V]).tolist(),
            face_verts=np.concatenate([V.ravel() for V in bf_V]).tolist(),
            n_el=len(el_V), eps=1e-9, leaf=32, eta=2.0, near_factor=1e30,
        )
        n = G.ndof()
        worst = 0.0
        for a in range(len(el_V)):
            e = [0.0] * n
            e[a] = 1.0
            Gaa = G.matvec(e)[a]                       # the raw volume-self diagonal G[a][a]
            gold = _gold_self(el_V[a])
            rel = Gaa / gold - 1.0
            worst = max(worst, abs(rel))
            # band 0.5%: the fixed 4-pt Gauss-Duffy rule lands ~0.01-0.2%; the old buggy rule was ~-6.5%
            assert abs(rel) < 5e-3, (
                f"cell {a}: order-0 volume self {Gaa:.6e} vs gold {gold:.6e} = {100*rel:+.3f}% "
                f"(was ~-6.5% before the 2026-06-12 _bary_tet(3) -> Gauss-Duffy fix)"
            )
        assert worst < 5e-3
