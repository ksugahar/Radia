r"""Multilayer thermal resistance / die-stack junction temperature -- regression test (#50).

Heat-generating top layer -> conduction down through passive layers (series thermal R=L/k) ->
cooled base. Drop across passive layer i = flux*L_i/k_i; junction rise = q L1^2/(2k1) +
Q*(L2/k2+L3/k3), Q=q L1. radia's solve_heat_steady (region-wise k) reproduces the interface
temperatures; the closed-form thermal network is the tool-independent gate. The
electronics-cooling junction-to-ambient block."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import thermal_resistance_series

WD = 0.02
L1, L2, L3 = 0.002, 0.002, 0.002
K1, K2, K3 = 200.0, 10.0, 400.0
Q_VOL = 1.0e7
QF = Q_VOL * L1


def test_thermal_resistance_series():
    R = thermal_resistance_series([L2, L3], [K2, K3])
    assert math.isclose(R, L2 / K2 + L3 / K3, rel_tol=1e-12)
    # the low-k layer (solder, k2=10) dominates over the high-k spreader (k3=400)
    assert (L2 / K2) > 0.9 * R
    # series adds; a thicker low-k layer raises R
    assert thermal_resistance_series([2 * L2, L3], [K2, K3]) > R
