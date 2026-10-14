import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (C0, MU0,
                                               rectangular_waveguide_cutoff,
                                               rectangular_waveguide_te10_port_normalization,
                                               waveguide_wave_impedance)


A, B = 0.02286, 0.01016


def test_te10_port_normalization_power_integral():
    row = rectangular_waveguide_te10_port_normalization(10.0e9, A, B, power_w=1.0)
    fc = rectangular_waveguide_cutoff(A, B, 1, 0)
    z_te = waveguide_wave_impedance(10.0e9, fc, "TE")["Z"]

    assert row["fc"] == pytest.approx(fc)
    assert row["Z_TE_ohm"] == pytest.approx(z_te)
    assert row["sin2_area_integral_m2"] == pytest.approx(0.5 * A * B)
    assert row["H_x_peak_A_per_m"] == pytest.approx(row["E_y_peak_V_per_m"] / z_te)
    assert row["poynting_power_W"] == pytest.approx(1.0, rel=1e-14)
    assert row["poynting_abs_error_W"] < 1e-14


def test_te10_longitudinal_h_ratio_and_power_scaling():
    one = rectangular_waveguide_te10_port_normalization(10.0e9, A, B, power_w=1.0)
    four = rectangular_waveguide_te10_port_normalization(10.0e9, A, B, power_w=4.0)
    zero = rectangular_waveguide_te10_port_normalization(10.0e9, A, B, power_w=0.0)

    kc = math.pi / A
    assert one["H_z_over_H_x_peak"] == pytest.approx(kc / one["beta"])
    assert four["E_y_peak_V_per_m"] == pytest.approx(2.0 * one["E_y_peak_V_per_m"])
    assert four["H_x_peak_A_per_m"] == pytest.approx(2.0 * one["H_x_peak_A_per_m"])
    assert four["poynting_power_W"] == pytest.approx(4.0, rel=1e-14)
    assert zero["E_y_peak_V_per_m"] == 0.0
    assert zero["H_x_peak_A_per_m"] == 0.0
    assert math.isinf(zero["H_z_over_H_x_peak"])


def test_te10_port_normalization_frequency_trends():
    low = rectangular_waveguide_te10_port_normalization(8.2e9, A, B)
    mid = rectangular_waveguide_te10_port_normalization(10.0e9, A, B)
    high = rectangular_waveguide_te10_port_normalization(12.4e9, A, B)

    assert low["fc"] == pytest.approx(0.5 * C0 / A)
    assert low["v_group"] < mid["v_group"] < high["v_group"]
    assert low["Z_TE_ohm"] > mid["Z_TE_ohm"] > high["Z_TE_ohm"] > MU0 * C0
    assert low["E_y_peak_V_per_m"] > mid["E_y_peak_V_per_m"] > high["E_y_peak_V_per_m"]


def test_te10_port_normalization_validation():
    with pytest.raises(ValueError):
        rectangular_waveguide_te10_port_normalization(5.0e9, A, B)
    with pytest.raises(ValueError):
        rectangular_waveguide_te10_port_normalization(10.0e9, -A, B)
    with pytest.raises(ValueError):
        rectangular_waveguide_te10_port_normalization(10.0e9, A, B, power_w=-1.0)
