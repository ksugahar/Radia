r"""Two-wire line TOTAL loop inductance (external + internal) -- regression test (#48).

L_loop = (mu0/pi) acosh(D/2a) + 2*(mu0/8pi) = external (#30) + both wires' internal (#36).
radia partitions the FE field energy (air -> external, wire interiors -> internal); the closed
form is the tool-independent gate. The +-I pair's energy converges, so the internal part is
exact (mu0/8pi each, radius-independent) and the external is the acosh term to ~1.8 %."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (two_wire_external_inductance,
                                           two_wire_loop_inductance,
                                           internal_inductance_round_wire)

MU0 = 4e-7 * math.pi
A, D, ROUT, I = 0.001, 0.010, 0.12, 1.0


def test_loop_inductance_formula():
    Le = two_wire_external_inductance(D, A)
    Ll = two_wire_loop_inductance(D, A)
    # loop = external + 2 * internal (mu0/8pi each)
    assert math.isclose(Ll, Le + MU0 / (4 * math.pi), rel_tol=1e-12)
    assert math.isclose(Ll, Le + 2 * internal_inductance_round_wire(), rel_tol=1e-12)
    assert Ll > Le                                   # loop exceeds external-only
    # the internal term is radius-independent (mu0/4pi)
    assert math.isclose(two_wire_loop_inductance(D, A) - two_wire_external_inductance(D, A),
                        two_wire_loop_inductance(2 * D, 2 * A) - two_wire_external_inductance(2 * D, 2 * A),
                        rel_tol=1e-12)
