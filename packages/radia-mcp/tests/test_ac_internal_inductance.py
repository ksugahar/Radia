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
