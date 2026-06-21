# -*- coding: utf-8 -*-
r"""Energy-based (vector Play) hysteresis in the radia-ngsolve core -- physics-gated.

Locks the lab-canonical energy-based hysteresis (Henrotte 2006 / Francois-Lavet 2013 / Jacques 2018;
the same B-input vector Play model as Radia's MatPlayHysteresis / MatEnergyHysteresis, see
radia.hysteresis_io) now available in the NGSolve core (radia_ngsolve.hysteresis).  Every property is
gated against a closed form or an exact physical law -- no fitted coefficients:

 (1) loop area = the analytic single-cell form  4 a eta (Bm - eta)  (and the multi-cell sum) to
     machine precision -- the model is exact;
 (2) the B-H loop closes (return-point memory);
 (3) STEINMETZ: geometric (~1/eta) thresholds give a loop-area loss ~ Bm^2 -- the empirical
     P_hyst = k_h f Bm^2 with k_h DERIVED from the cell parameters, not fitted (the principled
     replacement for coreloss.steinmetz_loss_density's fitted k_h);
 (4) RAYLEIGH: uniform thresholds give the low-field cubic loss ~ Bm^3;
 (5) CONGRUENCY: equal-amplitude minor loops are congruent regardless of bias (the Play hallmark);
 (6) 2nd LAW: the friction dissipation dD = sum_k a_k eta_k |dp_k| is >= 0 at every step and sums
     over a cycle to exactly the loop area (energy-consistent, the whole point of "energy-based").
"""
import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.hysteresis import (PlayHysteresis, rayleigh_cells, steinmetz_cells,
                                                play_cells, dissipation_increments)

_trapz = getattr(np, "trapezoid", None) or np.trapz


def test_single_cell_loop_area_is_analytic():
    m = PlayHysteresis(eta=[0.0, 0.5], a=[100.0, 200.0])
    for Bm in (0.8, 1.2, 2.0):
        num = m.loss_per_cycle(Bm, n=2001)
        ana = m.analytic_loss_per_cycle(Bm)              # 4*200*0.5*(Bm-0.5)
        assert abs(num - ana) / ana < 1e-4, f"Bm={Bm}: {num} vs {ana}"
    # below threshold -> no loss (purely reversible)
    assert m.loss_per_cycle(0.4, n=2001) < 1e-6


def test_multicell_loop_area_is_analytic_sum():
    m = PlayHysteresis(eta=[0.0, 0.3, 0.6, 0.9, 1.2], a=[80.0, 120.0, 100.0, 80.0, 60.0])
    for Bm in (1.0, 1.5):
        num = m.loss_per_cycle(Bm, n=2001)
        ana = m.analytic_loss_per_cycle(Bm)
        assert abs(num - ana) / ana < 1e-4, f"Bm={Bm}: {num} vs {ana}"


def test_loop_closes():
    m = PlayHysteresis(eta=[0.0, 0.3, 0.6, 0.9, 1.2], a=[80.0, 120.0, 100.0, 80.0, 60.0])
    B, H = m.steady_loop(1.5, n=2001)
    assert abs(H[0] - H[-1]) < 1e-9 * (abs(H).max() + 1) and abs(B[0] - B[-1]) < 1e-12


def test_steinmetz_quadratic_loss():
    """Geometric thresholds (~1/eta density) -> hysteresis loss ~ Bm^2 (Steinmetz), k_h derived."""
    s = steinmetz_cells(eta_min=0.05, eta_max=1.4, K=40, a_each=120.0)
    Bms = np.array([0.4, 0.6, 0.8, 1.0, 1.2])
    loss = np.array([s.loss_per_cycle(b, n=4001) for b in Bms])
    beta = np.polyfit(np.log(Bms), np.log(loss), 1)[0]
    print(f"Steinmetz beta={beta:.3f} (loss ~ Bm^2)")
    assert 1.9 < beta < 2.3, f"geometric-threshold loss should be ~ Bm^2, got beta={beta:.3f}"


def test_rayleigh_cubic_loss():
    """Uniform thresholds -> low-field cubic loss ~ Bm^3 (Rayleigh)."""
    r = rayleigh_cells(eta_max=1.5, K=40, a_each=120.0)
    Bms = np.array([0.4, 0.6, 0.8, 1.0, 1.2])
    loss = np.array([r.loss_per_cycle(b, n=4001) for b in Bms])
    beta = np.polyfit(np.log(Bms), np.log(loss), 1)[0]
    print(f"Rayleigh beta={beta:.3f} (loss ~ Bm^3)")
    assert 2.85 < beta < 3.15, f"uniform-threshold loss should be ~ Bm^3, got beta={beta:.3f}"


def test_minor_loop_congruency():
    """Equal-amplitude minor loops are congruent regardless of DC bias (the Play-model property)."""
    s = steinmetz_cells(eta_min=0.05, eta_max=1.4, K=40, a_each=120.0)

    def minor_area(bias, dB, n=1201, cycles=5):
        th = np.linspace(0, 2 * math.pi, n); osc = bias + dB * np.cos(th)
        pre = np.linspace(0, bias + dB, 600)
        wave = np.concatenate([pre, np.tile(osc[:-1], cycles), osc])
        H, _ = s.bh_history(wave)
        return abs(_trapz(H[-n:], osc))

    a1, a2 = minor_area(0.3, 0.4), minor_area(0.6, 0.4)
    assert abs(a1 - a2) / max(a1, 1e-12) < 1e-6, f"minor loops not congruent: {a1} vs {a2}"


def test_second_law_dissipation_nonnegative_and_closes():
    """Friction dissipation dD = sum_k a_k eta_k |dp_k| >= 0 always, and sums to the loop area."""
    m = PlayHysteresis(eta=[0.0, 0.3, 0.6, 0.9, 1.2], a=[80.0, 120.0, 100.0, 80.0, 60.0])
    n = 2001
    th = np.linspace(0, 2 * math.pi, n)
    one = 1.5 * np.cos(th)
    wave = np.concatenate([np.tile(one[:-1], 3), one])
    dD = dissipation_increments(m, wave)
    cyc = dD[-(n - 1):]
    area = m.analytic_loss_per_cycle(1.5)
    print(f"min dD={dD.min():.2e}, sum dD/cycle={cyc.sum():.4f}, loop area={area:.4f}")
    assert dD.min() >= -1e-12, "dissipation must be non-negative (2nd law)"
    assert abs(cyc.sum() - area) / area < 1e-4, "dissipation must sum to the loop area"


def test_play_cells_standard_discretization():
    """play_cells: equal-interval thresholds eta_k=(k-1/2)Bmax/N (Matsuo/Hane Eq.3) + analytic loss."""
    Bmax, N, nu = 1.5, 10, 100.0
    m = play_cells(Bmax, N, nu)
    expect = (np.arange(1, N + 1) - 0.5) * Bmax / N
    assert np.allclose(m.eta[1:], expect) and m.eta[0] == 0.0
    num = m.loss_per_cycle(1.2, n=3001)
    assert abs(num - m.analytic_loss_per_cycle(1.2)) / m.analytic_loss_per_cycle(1.2) < 1e-3


def test_h_axis_congruency_of_shape():
    """B-input play model: minor loops of the same B-range have IDENTICAL H-shape regardless of DC bias
    (H-axis congruency -- the property the lab measured in NOES and the B-input model reproduces)."""
    s = steinmetz_cells(eta_min=0.05, eta_max=1.4, K=40, a_each=120.0)
    n, dB = 1201, 0.4

    def centered_minor(bias, cycles=6):
        th = np.linspace(0, 2 * math.pi, n); osc = bias + dB * np.cos(th)
        wave = np.concatenate([np.linspace(0, bias + dB, 600), np.tile(osc[:-1], cycles), osc])
        H, _ = s.bh_history(wave)
        Hc = H[-n:]
        return dB * np.cos(th), Hc - Hc.mean()       # B relative to bias, H relative to its mean

    B1, H1 = centered_minor(0.3)
    B2, H2 = centered_minor(0.6)                      # same B-grid (relative), different bias
    excursion = H1.max() - H1.min()
    rel = np.max(np.abs(H1 - H2)) / excursion
    print(f"H-axis congruency: max|dH| / H-excursion = {rel:.2e}")
    assert rel < 1e-6, f"B-input minor loops must be H-axis congruent across bias, got {rel:.2e}"


def main():
    test_play_cells_standard_discretization()
    test_h_axis_congruency_of_shape()
    test_single_cell_loop_area_is_analytic()
    test_multicell_loop_area_is_analytic_sum()
    test_loop_closes()
    test_steinmetz_quadratic_loss()
    test_rayleigh_cubic_loss()
    test_minor_loop_congruency()
    test_second_law_dissipation_nonnegative_and_closes()
    print("[OK] energy-based Play hysteresis in the radia-ngsolve core: analytic loop area, loop "
          "closure, Steinmetz Bm^2 / Rayleigh Bm^3 loss (no fitted k_h), congruency, 2nd-law dissipation.")


if __name__ == "__main__":
    main()
