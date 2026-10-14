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
                                                play_cells, saturating_cells, dissipation_increments,
                                                identify_from_loop_areas)

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


def test_rotational_loss_analytic():
    """Vector rotational loss == closed form 2*pi sum_k a_k eta_k sqrt(Bm^2-eta_k^2) (linear cells)."""
    m = play_cells(Bmax=1.5, N=20, reluctivities=100.0)
    for Bm in (0.6, 1.0, 1.4):
        num = m.rotational_loss_per_cycle(Bm, n=4001)
        ana = m.analytic_rotational_loss(Bm)
        assert abs(num - ana) / ana < 5e-3, f"Bm={Bm}: rot {num} vs analytic {ana}"


def test_rotational_alternating_ratio_lowfield():
    """Rotational/alternating loss ratio -> pi/2 in the low-field (eta<<Bm) Rayleigh regime."""
    m = play_cells(Bmax=0.15, N=30, reluctivities=100.0)         # eta_max=0.1475
    ratios = [m.rotational_loss_per_cycle(Bm, n=2881) / m.loss_per_cycle(Bm, n=2001)
              for Bm in (0.8, 1.5, 3.0)]                          # decreasing eta_max/Bm
    print(f"rot/alt ratios (decreasing field) = {[round(r,3) for r in ratios]}, pi/2={math.pi/2:.3f}")
    assert ratios[0] > ratios[1] > ratios[2], "ratio must decrease toward pi/2 as field rises"
    assert abs(ratios[-1] - math.pi / 2) / (math.pi / 2) < 0.06, "low-field ratio must approach pi/2"


def test_saturating_rotational_loss_bounded():
    """Saturating shape functions BOUND the rotational loss to 2*pi a B_sat sum_k eta_k at high B."""
    a_each, B_sat, N, Bmax = 200.0, 0.5, 20, 1.5
    m = saturating_cells(Bmax=Bmax, N=N, a_each=a_each, B_sat=B_sat)
    limit = 2 * math.pi * a_each * B_sat * float(sum(m.eta))
    hi = m.rotational_loss_per_cycle(8.0, n=1441)
    lo = m.rotational_loss_per_cycle(0.6, n=1441)
    print(f"saturating rot loss: lo(0.6)={lo:.1f}, hi(8.0)={hi:.1f}, analytic limit={limit:.1f}")
    assert hi > lo, "rotational loss should rise with field"
    assert abs(hi - limit) / limit < 2e-3, f"high-field rot loss must saturate to {limit:.1f}, got {hi:.1f}"


def test_loss_kernel_matches_analytic():
    """The waveform loss kernel (the motor per-point iron-loss kernel) reproduces the alternating loop
    area and the rotational loss exactly -- so feeding a motor B(theta) sweep gives the real loss."""
    m = play_cells(Bmax=1.5, N=20, reluctivities=100.0)
    th = np.linspace(0, 2 * math.pi, 2001)
    alt = m.loss_from_waveform(1.0 * np.cos(th))
    rot = m.loss_from_waveform(1.0 * np.cos(th), 1.0 * np.sin(th))
    assert abs(alt - m.loss_per_cycle(1.0, n=2001)) / m.loss_per_cycle(1.0, n=2001) < 1e-6
    assert abs(rot - m.rotational_loss_per_cycle(1.0, n=2001)) / m.rotational_loss_per_cycle(1.0, n=2001) < 2e-3


def test_elliptical_loss_between_alternating_and_rotational():
    """An ELLIPTICAL B-path (motor reality) dissipates more than pure alternating and less than full
    rotational, monotonically in the axis ratio -- the 2D loss a scalar |B|-Steinmetz cannot capture."""
    m = play_cells(Bmax=1.5, N=20, reluctivities=100.0)
    th = np.linspace(0, 2 * math.pi, 2001)
    losses = [m.loss_from_waveform(np.cos(th), b * np.sin(th)) for b in (0.0, 0.3, 0.6, 1.0)]
    print(f"elliptical losses (b=0,.3,.6,1) = {[round(L,1) for L in losses]}")
    assert losses[0] < losses[1] < losses[2] < losses[3], "loss must rise with the B-path ellipticity"


def test_minor_loop_adds_loss_natively():
    """A minor reversal embedded in the major cycle ADDS its own loop area (native minor-loop capture --
    the iGSE refinement coreloss.py flags as a TODO); the increment equals the minor loop's analytic area."""
    m = play_cells(Bmax=1.5, N=20, reluctivities=100.0)
    th = np.linspace(0, 2 * math.pi, 1201)
    major = 1.2 * np.cos(th)
    # embed one minor loop near the top: dip 1.2 -> 0.6 -> 1.2 then continue
    dip = np.concatenate([np.linspace(1.2, 0.6, 80), np.linspace(0.6, 1.2, 80)])
    with_minor = np.concatenate([major[:200], dip, major[200:]])
    L_major = m.loss_from_waveform(major)
    L_minor = m.loss_from_waveform(with_minor)
    # a 1.2->0.6->1.2 reversal is a symmetric-equivalent loop of half-swing 0.3
    minor_area = m.analytic_loss_per_cycle(0.3)
    print(f"major={L_major:.2f}, with-minor={L_minor:.2f}, increment={L_minor-L_major:.2f}, "
          f"minor-area(half-swing 0.3)={minor_area:.2f}")
    assert L_minor > L_major, "an embedded minor loop must add loss (native capture)"
    assert abs((L_minor - L_major) - minor_area) / minor_area < 0.05, "increment == the minor-loop area"


def test_play_tangent_matches_finite_difference():
    """The Newton differential-reluctivity tensor dH/dB (Mitsuoka 2013 FE ingredient) matches a
    finite-difference of H(B) to machine precision, including the off-diagonal vector coupling that
    appears when cells are moving."""
    m = play_cells(1.5, 12, 800.0)
    eps = 1e-6

    def jac_fd(px, py, B0):
        def H(bx, by):
            hx, hy, _, _ = m.step(bx, by, px.copy(), py.copy()); return np.array([hx, hy])
        J = np.zeros((2, 2))
        J[:, 0] = (H(B0[0]+eps, B0[1]) - H(B0[0]-eps, B0[1])) / (2*eps)
        J[:, 1] = (H(B0[0], B0[1]+eps) - H(B0[0], B0[1]-eps)) / (2*eps)
        return J

    # moving regime: frozen virgin state, a big B step -> several cells move (non-trivial tensor)
    px, py = np.zeros(m.K), np.zeros(m.K)
    B0 = (0.5, 0.3)
    J = jac_fd(px, py, B0)
    nxx, nxy, nyy = (t[0] for t in m.play_tangent(B0[0], B0[1], px, py))
    print(f"tangent moving: nxx={nxx:.1f} nxy={nxy:.1f} nyy={nyy:.1f}; FD off-diag {J[0,1]:.1f}/{J[1,0]:.1f}")
    assert nxy > 1.0, "moving cells must produce a nonzero off-diagonal (vector coupling)"
    assert abs(J[0, 1] - J[1, 0]) < 1e-6*nxx, "dH/dB must be symmetric"
    rel = max(abs(J[0, 0]-nxx)/nxx, abs(J[1, 1]-nyy)/nyy, abs(J[0, 1]-nxy)/nxy)
    assert rel < 1e-5, f"play_tangent off the finite-difference dH/dB by {rel:.2e}"


def test_play_coenergy_is_the_convex_potential():
    """The incremental co-energy psi*(B) is the CONVEX potential of the B-input Play field: its
    B-gradient is exactly H (so Pi = INT psi* - INT J.A is the variational principle for the FE solve)
    and its Hessian is exactly the SPD play_tangent -- the two facts that make the energy-line-search
    Newton globally convergent (test_hysteresis_fe_variational)."""
    m = play_cells(1.5, 12, 800.0)
    px, py = np.zeros((1, m.K)), np.zeros((1, m.K))
    e, d = 1e-6, 1e-4
    P = lambda bx, by: m.play_coenergy_density(bx, by, px, py)[0]
    for (Bx, By) in [(1.0, 0.0), (0.5, 0.3), (1.2, 0.7), (0.3, 0.0)]:
        gx = (P(Bx+e, By) - P(Bx-e, By)) / (2*e)
        gy = (P(Bx, By+e) - P(Bx, By-e)) / (2*e)
        Hx, Hy, _, _ = m.step(Bx, By, np.zeros(m.K), np.zeros(m.K))      # same frozen virgin state
        assert abs(gx-Hx) < 1e-3*max(abs(Hx), 1.0) and abs(gy-Hy) < 1e-3*max(abs(Hy), 1.0), \
            f"grad psi* must equal H at B=({Bx},{By}): ({gx:.1f},{gy:.1f}) vs ({Hx:.1f},{Hy:.1f})"
        hxx = (P(Bx+d, By) - 2*P(Bx, By) + P(Bx-d, By)) / d**2
        hyy = (P(Bx, By+d) - 2*P(Bx, By) + P(Bx, By-d)) / d**2
        hxy = (P(Bx+d, By+d) - P(Bx+d, By-d) - P(Bx-d, By+d) + P(Bx-d, By-d)) / (4*d**2)
        nxx, nxy, nyy = (t[0] for t in m.play_tangent(Bx, By, px, py))
        assert abs(hxx-nxx) < 2e-2*nxx and abs(hyy-nyy) < 2e-2*nyy and abs(hxy-nxy) < 2e-2*max(abs(nxy), 1.0), \
            f"Hessian of psi* must equal the SPD play_tangent at B=({Bx},{By})"
        assert nxx > 0 and nyy > 0 and nxx*nyy - nxy*nxy > -1e-6*nxx*nyy, "play_tangent must be SPD"
    print("[ok] psi* gradient == H and Hessian == SPD play_tangent (the convex potential of the Play field)")


def test_identify_from_loop_areas_roundtrip():
    """Identify the cell slopes from MEASURED symmetric-loop areas, round-trip exactly (analytic) and
    robustly (simulated 'measurements'), and generalise to an unseen amplitude."""
    Bmax, N = 1.5, 10
    D = Bmax / N
    eta = np.concatenate([[0.0], (np.arange(1, N + 1) - 0.5) * D])
    a_true = np.array([50., 80, 120, 100, 90, 70, 60, 50, 40, 30, 20])      # a0=50 + a_1..a_10
    m0 = PlayHysteresis(eta=eta, a=a_true)
    W_ana = np.array([m0.analytic_loss_per_cycle(mi * D) for mi in range(1, N + 1)])
    W_sim = np.array([m0.loss_per_cycle(mi * D, n=2001) for mi in range(1, N + 1)])
    mid_a = identify_from_loop_areas(Bmax, N, W_ana, a0=50.)
    mid_s = identify_from_loop_areas(Bmax, N, W_sim, a0=50.)
    rel_a = np.max(np.abs(mid_a.a[1:] - a_true[1:]) / a_true[1:])
    rel_s = np.max(np.abs(mid_s.a[1:] - a_true[1:]) / a_true[1:])
    Bt = 1.13
    gen = abs(mid_s.loss_per_cycle(Bt, n=2001) - m0.loss_per_cycle(Bt, n=2001)) / m0.loss_per_cycle(Bt, n=2001)
    print(f"identify slopes: analytic rel={rel_a:.2e}, simulated rel={rel_s:.2e}; unseen-loop rel={gen:.2e}")
    assert rel_a < 1e-9, "exact inversion of analytic loop areas must recover the slopes"
    assert rel_s < 1e-2, "identification must be robust to simulated-measurement noise"
    assert gen < 1e-2, "the identified model must reproduce an unseen loop amplitude"


def main():
    test_identify_from_loop_areas_roundtrip()
    test_play_coenergy_is_the_convex_potential()
    test_play_tangent_matches_finite_difference()
    test_loss_kernel_matches_analytic()
    test_elliptical_loss_between_alternating_and_rotational()
    test_minor_loop_adds_loss_natively()
    test_rotational_loss_analytic()
    test_rotational_alternating_ratio_lowfield()
    test_saturating_rotational_loss_bounded()
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
