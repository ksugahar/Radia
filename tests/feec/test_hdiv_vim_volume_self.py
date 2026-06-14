"""Golden: the order-0 (RT0) ANALYTIC charge-Gram VOLUME self-energy

    G[a][a] = (1/4pi) INT_T INT_T 1/|x-y| dy dx

matches the converged value of that integral.

REGRESSION GUARD for the 2026-06-12 fix.  The order-0 analytic ctor's CELL outer-quadrature was a crude
10-point equal-weight _bary_tet(3) rule that under-integrated the volume self-energy by ~6.5%.  It was
INVISIBLE to every other HDiv-VIM golden because they all use UNIFORM M -> div M = 0 -> the volume charge
(and hence the volume self) never participates; the demag is surface-charge only.  AND the Python dense
reference (analytic_charge_gram) shared the same crude rule, so the "validated against dense" cross-check was
blind in the same direction.  This test reads the RAW volume-self Gram diagonal directly, so div M != 0 is
exercised and the defect cannot hide.

The reference is the SAME integral with an exact analytic inner (phi_tet, the Wilton uniform-tet potential)
and a FINE Gauss-Duffy outer rule -- it self-converges, so it is an independent gold standard, not a
re-statement of the production rule.
"""
import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Box, OCCGeometry, Pnt

import radia._radia_pybind as rp
from radia.vim._core import phi_tet, _gauss_duffy_tet


def _tet_vol(V):
    return abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]))) / 6.0


def _gold_self(V):
    """(1/4pi) INT_T phi_tet(V, x) dx via a FINE (o=12) Gauss-Duffy outer rule; phi_tet is the EXACT
    Wilton uniform-tet potential, so this self-converges to the true volume self-energy."""
    V = np.asarray(V, float)
    refP, refW = _gauss_duffy_tet(12)
    P = V[0] + refP @ (V[1:] - V[0])
    return float(np.sum(refW * 6.0 * _tet_vol(V) * phi_tet(V, P)) / (4.0 * np.pi))


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
