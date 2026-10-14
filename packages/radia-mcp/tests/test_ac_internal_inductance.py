r"""AC internal-inductance roll-off of a round wire (skin effect) -- regression test (#45).

L_int(omega)/L_int_dc = (4/q)[ber ber'+bei bei']/[ber'^2+bei'^2] (Kelvin, q=sqrt(2)a/delta),
the inductive twin of the skin-effect Rac/Rdc. radia's `solve_planar_eddy` on the WIRE ALONE
(A_z=0 on the surface -> Z=Vc/I is the pure internal impedance) reproduces both Rac=Re(Z) and
L_int=Im(Z)/omega. L_int -> mu0/8pi (#36 DC limit) at low freq, rolls off ~4/q high. Pure
Kelvin closed form -> tool-independent."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (internal_inductance_round_wire,
                                           round_wire_internal_impedance,
                                           skin_effect_resistance_ratio,
                                           skin_effect_internal_inductance_ratio)

MU0 = 4e-7 * math.pi
A, SIGMA, I = 1e-3, 5.8e7, 1.0


def test_kelvin_ratio_limits():
    pytest.importorskip("scipy")
    # low frequency: both ratios -> 1 (DC: uniform current, L_int = mu0/8pi)
    assert abs(skin_effect_resistance_ratio(0.01) - 1.0) < 1e-3
    assert abs(skin_effect_internal_inductance_ratio(0.01) - 1.0) < 1e-3
    # high frequency asymptotes (q = sqrt(2) a/delta): Rac/Rdc -> q/(2 sqrt2) + 1/4,
    # L_int/L_dc -> 2 sqrt2 / q  (current confined to the ~delta skin)
    q = 30.0
    assert abs(skin_effect_resistance_ratio(q) - (q / (2 * math.sqrt(2)) + 0.25)) / (q / (2 * math.sqrt(2))) < 0.02
    assert abs(skin_effect_internal_inductance_ratio(q) - 2 * math.sqrt(2) / q) / (2 * math.sqrt(2) / q) < 0.05
    # monotone: Rac rises, L_int falls with q
    assert skin_effect_resistance_ratio(1) < skin_effect_resistance_ratio(2) < skin_effect_resistance_ratio(4)
    assert skin_effect_internal_inductance_ratio(1) > skin_effect_internal_inductance_ratio(2) > skin_effect_internal_inductance_ratio(4)


def test_round_wire_internal_impedance_decomposes_into_r_and_l():
    pytest.importorskip("scipy")
    freq = 20_000.0
    res = round_wire_internal_impedance(A, SIGMA, freq)
    omega = 2.0 * math.pi * freq
    q = A * math.sqrt(omega * MU0 * SIGMA)
    rdc = 1.0 / (SIGMA * math.pi * A * A)
    lint_dc = internal_inductance_round_wire()

    assert math.isclose(res["q"], q, rel_tol=1e-12)
    assert math.isclose(res["skin_depth"], math.sqrt(2.0) * A / q, rel_tol=1e-12)
    assert math.isclose(res["Rdc_per_m"], rdc, rel_tol=1e-12)
    assert math.isclose(res["Rac_per_m"], rdc * skin_effect_resistance_ratio(q), rel_tol=1e-12)
    assert math.isclose(res["Lint_per_m"], lint_dc * skin_effect_internal_inductance_ratio(q), rel_tol=1e-12)
    assert math.isclose(res["Z_per_m"].real, res["Rac_per_m"], rel_tol=1e-12)
    assert math.isclose(res["Z_per_m"].imag, omega * res["Lint_per_m"], rel_tol=1e-12)


def test_round_wire_internal_impedance_limits_and_validation():
    pytest.importorskip("scipy")
    dc = round_wire_internal_impedance(A, SIGMA, 0.0)
    assert dc["Z_per_m"] == complex(dc["Rdc_per_m"], 0.0)
    assert math.isinf(dc["skin_depth"])
    assert math.isclose(dc["Lint_per_m"], internal_inductance_round_wire(), rel_tol=1e-12)

    low = round_wire_internal_impedance(A, SIGMA, 1e-6)
    assert math.isclose(low["Rac_per_m"], low["Rdc_per_m"], rel_tol=1e-7)
    assert math.isclose(low["Lint_per_m"], internal_inductance_round_wire(), rel_tol=1e-7)

    high = round_wire_internal_impedance(A, SIGMA, 2.0e7)
    assert high["Rac_per_m"] > 10.0 * high["Rdc_per_m"]
    assert high["Lint_per_m"] < 0.2 * internal_inductance_round_wire()

    for bad in (
        lambda: round_wire_internal_impedance(0.0, SIGMA, 1.0),
        lambda: round_wire_internal_impedance(A, 0.0, 1.0),
        lambda: round_wire_internal_impedance(A, SIGMA, -1.0),
        lambda: round_wire_internal_impedance(A, SIGMA, 1.0, mu_r=0.0),
    ):
        with pytest.raises(ValueError):
            bad()
