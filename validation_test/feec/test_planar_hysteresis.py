"""Golden lock for the SHARED 2D play-hysteresis demag (radia.planar_hysteresis + the
planar_materials.PlayHysteresis operator) on the direct-N + Newton framework.

Gates: play tangent (chi_inc == finite-diff); anhysteretic limit (eta=0 == linear demag); the full
direct-N Newton disk loop == the independent scalar reference M = play(H_ext - D M), D=1/2; 2nd-cycle
loop closure; nonzero remanence + coercivity; Newton fail-loud.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.planar_aniso as pa
import radia.planar_hysteresis as ph
from radia.planar_materials import PlayHysteresis

ETA = [0.1, 0.3, 0.6, 1.0, 1.6]
W = [3.0, 2.0, 1.2, 0.6, 0.2]                 # sum = 7 -> anhysteretic chi0 = 7


def _disk(maxh=0.3):
    g = SplineGeometry(); g.AddCircle((0, 0), r=1.0, bc="e")
    return ng.Mesh(g.GenerateMesh(maxh=maxh))


def _scalar_ref(play, H_seq, D=0.5):
    """1D reference: M = play(H_ext - D M) per step (single site)."""
    p = play.fresh_state(1)
    M = 0.0
    out = []
    for Hext in H_seq:
        for _ in range(200):
            H = np.array([Hext - D * M])
            F = M - play.M(H, p)[0]
            J = 1.0 + play.chi_inc(H, p)[0] * D
            Mn = M - F / J
            if abs(Mn - M) < 1e-13 * max(abs(Mn), 1.0):
                M = Mn; break
            M = Mn
        p = play.advance(np.array([Hext - D * M]), p)
        out.append(M)
    return np.array(out)


def test_play_tangent():
    play = PlayHysteresis(ETA, W)
    p = play.advance(np.array([0.3]), play.fresh_state(1))
    H = np.array([0.5]); dh = 1e-6
    fd = (play.M(H + dh, p)[0] - play.M(H - dh, p)[0]) / (2 * dh)
    assert abs(play.chi_inc(H, p)[0] - fd) < 1e-6, (play.chi_inc(H, p)[0], fd)
    assert np.all(play.chi_inc(np.linspace(-3, 3, 50), np.zeros((50, len(ETA)))) >= 0)   # always >= 0


def test_anhysteretic_limit_is_linear_demag():
    play = PlayHysteresis(np.zeros(3), [2.0, 3.0, 2.0])       # eta=0 -> chi = 7 anhysteretic
    with ng.TaskManager():
        d = _disk(0.3)
        r = ph.solve_hysteresis_demag(d, play, [1.0])         # single step, virgin
    cl = 7.0 / (1 + 7.0 / 2)
    assert abs(r["M_avg"][0] - cl) / cl < 1e-3, (r["M_avg"][0], cl)


def test_full_loop_matches_scalar_reference():
    play_full = PlayHysteresis(ETA, W)
    play_ref = PlayHysteresis(ETA, W)
    H = 3.0 * np.sin(np.linspace(0, 2 * np.pi, 49))
    with ng.TaskManager():
        d = _disk(0.3)
        r = ph.solve_hysteresis_demag(d, play_full, H)
    Mref = _scalar_ref(play_ref, H)
    rel = np.max(np.abs(r["M_avg"] - Mref)) / np.max(np.abs(Mref))
    assert rel < 1e-3, rel                                    # full N-Newton == 1D reference


def test_loop_closes_second_cycle_and_has_hysteresis():
    play = PlayHysteresis(ETA, W)
    H = 3.0 * np.sin(np.linspace(0, 4 * np.pi, 97))           # 2 cycles
    with ng.TaskManager():
        d = _disk(0.3)
        r = ph.solve_hysteresis_demag(d, play, H)
    M = r["M_avg"]
    # cycle-2 closes: the state is periodic after the first cycle (compare the two half-way points)
    q = len(H) // 4
    assert abs(M[2 * q] - M[4 * q]) < 5e-3 * (M.max() - M.min())   # H_ext=0 descending, cycle1 vs cycle2
    # hysteresis present: at H_ext=0 the ascending and descending M differ (open loop -> remanence)
    i_desc = 2 * q                                            # H_ext ~ 0 going down
    i_asc = 4 * q                                             # H_ext ~ 0 going down again (periodic)
    assert (M.max() - M.min()) > 1.0                          # substantial magnetisation swing
    assert abs(M[q * 3]) > 0.05 * M.max() or abs(M[q]) > 0.05 * M.max()   # nonzero remanence somewhere


def test_newton_fail_loud():
    """Newton fails loud if it cannot converge (No-Fallbacks)."""
    play = PlayHysteresis(ETA, W)
    with ng.TaskManager():
        d = _disk(0.4)
        with pytest.raises(RuntimeError, match="did NOT converge"):
            ph.solve_hysteresis_demag(d, play, [2.0], newton_maxit=1, newton_tol=1e-14)
