# -*- coding: utf-8 -*-
"""Golden tests for radia.open_boundary (the exact Zs-DtN-CLN open boundary).

Locks the VERIFIED properties ported from the research demos
for the Zs-DtN-CLN open boundary:
  - the eddy DtN is exactly rational in q=sqrt(s); the Cauer-in-q ladder is exact
    at n+1 stages (NRMSE ~1e-15) for n=1..6 and well-conditioned (spread <1e3);
  - the ladder reproduces the exact symbol in general (R0, mu_sigma) units;
  - the wave (in s) and diffusion (in sqrt(s)) DtN share the SAME poles roots(theta_n);
  - the companion auxiliary-ODE rates are roots(theta_n), all Re<0 => passive/stable;
  - a Foster fit in s FLOORS + ILL-CONDITIONS (the structural contrast);
  - the sqrt(s) passive ladder has real negative poles (stable);
  - the eddy DtN is analytic/bounded in Re(s)>0 (no RHP pole => passive).
"""
import numpy as np
import pytest

import radia.open_boundary as ob

OMEGA = np.logspace(-1, 2, 60)
MODES = (1, 2, 3, 4, 5, 6)


def _nrmse(a, b):
    return float(np.sqrt(np.mean(np.abs(a - b) ** 2)) / np.sqrt(np.mean(np.abs(b) ** 2)))


@pytest.mark.parametrize("n", MODES)
def test_cauer_exact_at_n_plus_1_stages(n):
    """The Cauer-in-q ladder is EXACT at exactly n+1 stages, well-conditioned."""
    stages = ob.cauer_ladder(n)
    assert len(stages) == n + 1, f"n={n}: expected n+1 stages, got {len(stages)}"
    Zc = np.array([ob.eval_ladder(stages, 1j * w) for w in OMEGA])
    Zr = np.array([ob.eddy_dtn(n, 1j * w) for w in OMEGA])
    assert _nrmse(Zc, Zr) < 1e-10, f"n={n}: Cauer ladder not exact vs symbol"
    allc = np.abs(np.concatenate([np.asarray(s, float) for s in stages]))
    spread = float(np.max(allc) / np.min(allc[allc > 0]))
    assert spread < 1e3, f"n={n}: ladder ill-conditioned (spread {spread:.1e})"


@pytest.mark.parametrize("R0,mu_sigma", [(1.0, 1.0), (0.1, 0.5), (0.03, 4.0e-7 * np.pi * 5.8e7)])
def test_ladder_matches_symbol_general_units(R0, mu_sigma):
    """eval_ladder reproduces eddy_dtn for non-unit (R0, mu_sigma)."""
    for n in (1, 2, 3):
        stages = ob.cauer_ladder(n)
        Zc = np.array([ob.eval_ladder(stages, 1j * w, R0, mu_sigma) for w in OMEGA])
        Zr = np.array([ob.eddy_dtn(n, 1j * w, R0, mu_sigma) for w in OMEGA])
        assert _nrmse(Zc, Zr) < 1e-9, f"n={n} R0={R0} mu_sigma={mu_sigma}: ladder != symbol"


@pytest.mark.parametrize("n", MODES)
def test_wave_diffusion_share_reverse_bessel_poles(n):
    """Wave (in s) and diffusion (in q=sqrt(s)) DtN have the SAME poles roots(theta_n)."""
    A, den = ob.eddy_dtn_rational_q(n)
    qpoles = np.sort_complex(np.roots(den[::-1].copy()).astype(complex))
    rbr = np.sort_complex(ob.reverse_bessel_roots(n))
    assert np.max(np.abs(qpoles - rbr)) < 1e-9
    # diffusion DtN rational form reproduces the scipy K-Bessel reference
    for s in (0.3, 1.0, 5.0, 2.0 + 1.0j):
        q = np.sqrt(complex(s))
        val = (sum(A[k] * q ** k for k in range(len(A)))
               / sum(den[k] * q ** k for k in range(len(den))))
        assert abs(val - complex(ob.eddy_dtn(n, s))) / abs(complex(ob.eddy_dtn(n, s))) < 1e-9


@pytest.mark.parametrize("n", MODES)
def test_companion_poles_passive_stable(n):
    """Companion auxiliary-ODE rates = roots(theta_n), all Re<0 (passive/stable)."""
    poles = ob.companion_poles(n)
    assert len(poles) == n
    assert np.all(poles.real < 0.0), f"n={n}: a companion pole has Re>=0 (unstable)"


@pytest.mark.parametrize("n", (1, 2, 3))
def test_foster_in_s_floors_and_illconditions(n):
    """The structural contrast: a Foster fit in s cannot match the sqrt(s) DtN --
    it FLOORS (never exact) and ILL-CONDITIONS, unlike the exact Cauer-in-q."""
    def foster(M):
        s = 1j * OMEGA
        cols = [s / (s + pj) for pj in np.logspace(np.log10(OMEGA[0]), np.log10(OMEGA[-1]), M)]
        Amat = np.column_stack(cols + [np.ones_like(s)])
        rhs = np.array([ob.eddy_dtn(n, sv) for sv in s], complex)
        coef, *_ = np.linalg.lstsq(np.vstack([Amat.real, Amat.imag]),
                                   np.concatenate([rhs.real, rhs.imag]), rcond=None)
        nz = np.abs(coef[np.abs(coef) > 0])
        spread = float(np.max(np.abs(coef)) / (np.min(nz) if nz.size else 1.0))
        return _nrmse(Amat @ coef, rhs), spread

    e32, sp32 = foster(32)
    assert e32 > 1e-4, f"n={n}: Foster should floor (got {e32:.1e})"
    assert sp32 > 1e4, f"n={n}: Foster should ill-condition (spread {sp32:.1e})"


def test_sqrt_s_passive_ladder_stable():
    """The sqrt(s) memory ladder has real positive p (poles -p<0 => stable) + fits."""
    g, p, nrmse = ob.sqrt_s_passive_ladder(OMEGA, 12)
    assert np.all(p > 0.0), "sqrt(s) ladder pole not in LHP"
    assert np.all(g >= 0.0), "sqrt(s) ladder not passive (g_m<0)"
    assert nrmse < 5e-2


@pytest.mark.parametrize("n", MODES)
def test_eddy_dtn_analytic_in_right_half_plane(n):
    """G_n is finite over a Re(s)>0 grid (poles only on the non-physical sqrt(s)
    sheet) => passive, stable open boundary."""
    grid = [(sr + 1j * si) for sr in (0.05, 0.5, 2.0, 10.0) for si in (-8, -1, 0, 1, 8)]
    worst = max(abs(complex(ob.eddy_dtn(n, s))) for s in grid)
    assert worst < 1e3, f"n={n}: |G_n| blows up in Re(s)>0 (RHP pole)"
