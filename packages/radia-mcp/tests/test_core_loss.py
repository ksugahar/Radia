# -*- coding: utf-8 -*-
"""Core-loss separation post-processor (radia_ngsolve.coreloss) validated
against closed forms and the FFT harmonic path.

1. The classical-eddy term with k_c = pi^2 sigma d^2 / 6 reproduces the textbook
   lamination eddy loss P = pi^2 sigma d^2 f^2 Bpk^2 / 6 (exact).
2. Classical loss scales as f^2 and Bpk^2.
3. FFT of a pure sinusoid over one period recovers the fundamental amplitude;
   harmonic_loss_density then equals bertotti_loss_density (sinusoidal limit).
4. FFT of a 2-harmonic waveform recovers both amplitudes; classical loss sums
   per harmonic, P_class = k_c [ (f1 B1)^2 + (2 f1 B2)^2 ].
5. core_loss_in_region integrates density * area over a unit-square region.

This is the radia-ngsolve native rotor-sweep core loss -- FEMM cannot compute
DC-magnet iron loss without the .ans frequency hand-edit hack.
"""
import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.coreloss import (
    lamination_kc, classical_eddy_loss_density, bertotti_loss_density,
    steinmetz_loss_density, fft_amplitudes, harmonic_loss_density,
    core_loss_in_region)

SIGMA, D = 2.0e6, 0.5e-3          # electrical steel ~2 MS/m, 0.5 mm lamination
K_H, K_E = 150.0, 1.0             # representative hysteresis / excess coeffs


def test_classical_term_matches_closed_form():
    kc = lamination_kc(SIGMA, D)
    for Bpk, f in ((1.5, 50.0), (1.0, 400.0), (0.8, 1000.0)):
        closed = classical_eddy_loss_density(Bpk, f, SIGMA, D)
        # bertotti classical-only must equal the closed form
        berto = bertotti_loss_density(Bpk, f, k_h=0.0, k_c=kc, k_e=0.0)
        assert math.isclose(berto, closed, rel_tol=1e-12), (berto, closed)


def test_classical_scaling():
    kc = lamination_kc(SIGMA, D)
    p1 = bertotti_loss_density(1.0, 50.0, 0.0, kc)
    assert math.isclose(bertotti_loss_density(1.0, 100.0, 0.0, kc), 4 * p1, rel_tol=1e-12)  # f^2
    assert math.isclose(bertotti_loss_density(2.0, 50.0, 0.0, kc), 4 * p1, rel_tol=1e-12)  # Bpk^2


def test_fft_sinusoid_recovers_fundamental():
    N, Bm = 256, 1.3
    t = np.arange(N) / N
    wave = Bm * np.sin(2 * np.pi * t)
    A = fft_amplitudes(wave)
    assert abs(A[1] - Bm) < 1e-9
    assert abs(A[0]) < 1e-9                      # zero mean
    assert np.max(A[2:]) < 1e-9                  # no higher harmonics


def test_harmonic_reduces_to_bertotti_for_sinusoid():
    N, Bm, f1 = 256, 1.3, 50.0
    kc = lamination_kc(SIGMA, D)
    A = fft_amplitudes(Bm * np.sin(2 * np.pi * np.arange(N) / N))
    h = harmonic_loss_density(A, f1, K_H, kc, K_E)
    b = bertotti_loss_density(Bm, f1, K_H, kc, K_E)
    assert math.isclose(h["P_total"], b, rel_tol=1e-6), (h["P_total"], b)


def test_two_harmonic_classical_sum():
    N, B1, B2, f1 = 512, 1.2, 0.4, 50.0
    kc = lamination_kc(SIGMA, D)
    t = np.arange(N) / N
    wave = B1 * np.sin(2 * np.pi * t) + B2 * np.sin(2 * np.pi * 2 * t)
    A = fft_amplitudes(wave)
    assert abs(A[1] - B1) < 1e-9 and abs(A[2] - B2) < 1e-9
    h = harmonic_loss_density(A, f1, k_h=0.0, k_c=kc, k_e=0.0)
    expected = kc * ((f1 * B1) ** 2 + (2 * f1 * B2) ** 2)
    assert math.isclose(h["P_class"], expected, rel_tol=1e-6)


def test_core_loss_in_region_integrates_density():
    from ngsolve import Mesh, CoefficientFunction
    from netgen.occ import OCCGeometry, MoveTo
    sq = MoveTo(0, 0).Rectangle(1.0, 1.0).Face()
    sq.faces.name = "steel"
    mesh = Mesh(OCCGeometry(sq, dim=2).GenerateMesh(maxh=0.3))
    kc = lamination_kc(SIGMA, D)
    Bpk, f = 1.5, 50.0
    total = core_loss_in_region(CoefficientFunction(Bpk), mesh, "steel", f, K_H, kc, k_e=0.0)
    expected = bertotti_loss_density(Bpk, f, K_H, kc, 0.0) * 1.0   # area = 1 m^2
    assert math.isclose(total, expected, rel_tol=1e-9), (total, expected)


if __name__ == "__main__":
    test_classical_term_matches_closed_form()
    test_classical_scaling()
    test_fft_sinusoid_recovers_fundamental()
    test_harmonic_reduces_to_bertotti_for_sinusoid()
    test_two_harmonic_classical_sum()
    test_core_loss_in_region_integrates_density()
    print("[OK] core-loss separation validated: classical closed form, f^2/B^2 "
          "scaling, FFT harmonic path, region integration.")
