"""Validate the canonical induction-heating cylinder Joule-loss formula
against its analytic asymptotic limits."""

from __future__ import annotations

import math

import pytest

from radia.analytical_formulas.induction_heating import (
    cylinder_axial_eddy_loss,
    cylinder_axial_eddy_loss_small_ka,
    cylinder_axial_eddy_loss_thin_skin,
)


MU_0 = 4.0e-7 * math.pi


def _ka(a, omega, sigma, mu_r=1.0):
    return a * math.sqrt(omega * sigma * mu_r * MU_0)


# Cu workpiece geometry from the IH panel sample
SIGMA_CU = 5.8e7
A_CU = 0.02  # 20 mm radius
H0 = 1000.0  # A/m


# ---------------------------------------------------------------------------
# Asymptotic limits
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "ka, rtol",
    [
        (0.001, 1e-4),
        (0.01, 1e-4),
        (0.1, 1e-3),
    ],
)
def test_small_ka_asymptote(ka, rtol):
    """At small ka the full Bessel form must match the elementary
    Faraday + Joule eddy-current formula ``π σ ω² μ² H_0² a^4 / 16``.
    """
    omega = (ka / A_CU) ** 2 / (SIGMA_CU * MU_0)
    p_full = cylinder_axial_eddy_loss(H0, A_CU, omega, SIGMA_CU)
    p_small = cylinder_axial_eddy_loss_small_ka(H0, A_CU, omega, SIGMA_CU)
    assert p_full > 0, f"Joule loss must be positive, got {p_full}"
    assert p_small == pytest.approx(p_full, rel=rtol), (
        f"At ka={ka}, small-ka asymptote {p_small:.4e} disagrees with full "
        f"Bessel {p_full:.4e} by more than {rtol}"
    )


@pytest.mark.parametrize(
    "ka, rtol",
    [
        (10.0, 0.10),
        (30.0, 0.05),
        (100.0, 0.02),
    ],
)
def test_thin_skin_asymptote(ka, rtol):
    """At large ka the full Bessel form must match the planar surface-
    impedance form ``π a H_0² Re(Z_s)``.
    """
    omega = (ka / A_CU) ** 2 / (SIGMA_CU * MU_0)
    p_full = cylinder_axial_eddy_loss(H0, A_CU, omega, SIGMA_CU)
    p_thin = cylinder_axial_eddy_loss_thin_skin(H0, A_CU, omega, SIGMA_CU)
    assert p_thin == pytest.approx(p_full, rel=rtol), (
        f"At ka={ka}, thin-skin asymptote {p_thin:.4e} disagrees with full "
        f"Bessel {p_full:.4e} by more than {rtol}"
    )


# ---------------------------------------------------------------------------
# Physical sanity checks
# ---------------------------------------------------------------------------
def test_positive_joule_loss_at_all_ka():
    """Joule loss must be strictly positive for any ka > 0 (regression
    against the OCR'd Stafl §4.3 form which gave NEGATIVE losses)."""
    for ka in (0.01, 0.1, 1.0, 10.0, 100.0):
        omega = (ka / A_CU) ** 2 / (SIGMA_CU * MU_0)
        p = cylinder_axial_eddy_loss(H0, A_CU, omega, SIGMA_CU)
        assert p > 0, f"P({ka=}) = {p}; must be positive"


def test_quadratic_in_H_0():
    """Joule loss is quadratic in H_0 (single linearity check)."""
    omega = 6.28e4
    p1 = cylinder_axial_eddy_loss(1.0, A_CU, omega, SIGMA_CU)
    p3 = cylinder_axial_eddy_loss(3.0, A_CU, omega, SIGMA_CU)
    assert p3 == pytest.approx(9.0 * p1, rel=1e-13)


def test_zero_H_0_gives_zero_loss():
    omega = 6.28e4
    p = cylinder_axial_eddy_loss(0.0, A_CU, omega, SIGMA_CU)
    assert p == 0.0


def test_invalid_inputs_raise():
    omega = 6.28e4
    with pytest.raises(ValueError):
        cylinder_axial_eddy_loss(-1.0, A_CU, omega, SIGMA_CU)
    with pytest.raises(ValueError):
        cylinder_axial_eddy_loss(H0, -A_CU, omega, SIGMA_CU)
    with pytest.raises(ValueError):
        cylinder_axial_eddy_loss(H0, A_CU, -1.0, SIGMA_CU)
    with pytest.raises(ValueError):
        cylinder_axial_eddy_loss(H0, A_CU, omega, -1.0)
    with pytest.raises(ValueError):
        cylinder_axial_eddy_loss(H0, A_CU, omega, SIGMA_CU, mu_r=0.5)


# ---------------------------------------------------------------------------
# IH-panel-relevant absolute number cross-check
# ---------------------------------------------------------------------------
def test_canonical_ih_cu_workpiece_10kHz():
    """Cu workpiece, 20 mm radius, 10 kHz, NI=1000 A/m: thin-skin
    formula gives ~1.6 W/m. Used as the IH workpiece sanity benchmark.
    """
    omega = 2 * math.pi * 10_000  # 10 kHz
    p_full = cylinder_axial_eddy_loss(H0, A_CU, omega, SIGMA_CU)
    # Thin-skin estimate ~ pi a H^2 sqrt(omega mu / (2 sigma))
    delta = math.sqrt(2 / (omega * MU_0 * SIGMA_CU))
    p_thin = math.pi * A_CU * H0 ** 2 / (SIGMA_CU * delta)
    assert p_full == pytest.approx(p_thin, rel=0.05)
    # Magnitude order: 1 to 2 W/m
    assert 1.0 < p_full < 5.0, f"Cu 20mm 10kHz IH cylinder loss = {p_full} W/m, expected ~ 1-2 W/m"
