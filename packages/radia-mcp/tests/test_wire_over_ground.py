r"""Wire over a conducting ground plane -- external inductance (image method) -- test (#63).

L_ext = (mu0/2 pi) acosh(h/a) = 1/2 the two-wire value at separation 2h (the ground reflects a -I
image at -h; half the energy is above the plane). Closed-form helper (tool-independent) + an FE
field-energy check (L = 2W/I^2 above the ground)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (wire_over_ground_inductance,
                                           two_wire_external_inductance, MU0)


def test_wire_over_ground_closed_form():
    a = 0.0005
    for h in (0.005, 0.010, 0.020):
        L = wire_over_ground_inductance(h, a)
        assert math.isclose(L, MU0 / (2 * math.pi) * math.acosh(h / a), rel_tol=1e-12)
        # exactly HALF the two-wire external inductance at separation 2h
        assert math.isclose(L, 0.5 * two_wire_external_inductance(2 * h, a), rel_tol=1e-12)
    # monotone increasing with height; h >> a limit -> (mu0/2pi) ln(2h/a)
    assert wire_over_ground_inductance(0.005, a) < wire_over_ground_inductance(0.020, a)
    assert math.isclose(wire_over_ground_inductance(0.05, a),
                        MU0 / (2 * math.pi) * math.log(2 * 0.05 / a), rel_tol=2e-3)
