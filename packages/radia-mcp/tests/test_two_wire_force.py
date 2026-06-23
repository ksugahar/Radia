r"""Ampere's force law: force per unit length between two parallel wires -- regression test.

Two parallel wires a distance d apart carrying currents I1, I2 attract/repel with
F = mu0 I1 I2/(2 pi d) (two_wire_force_per_length); like currents attract. The analytic primitive of
the image-wire force (image_force_wire_iron = this with the image at 2h). The FE leg gets the force
as the Lorentz body force int J x B over one wire (only the OTHER wire's field nets a force; the
symmetric self-field integrates to zero). Tool-independent (Ampere's force law).
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (two_wire_force_per_length, image_force_wire_iron,
                                           MU0)

a, d, I = 0.002, 0.02, 1.0


def test_two_wire_force_closed_form():
    # Ampere's force law + the sign/relationship conventions
    F = two_wire_force_per_length(I, I, d)
    assert F == pytest.approx(MU0 * I * I / (2 * math.pi * d))      # mu0 I1 I2/(2 pi d)
    assert F == pytest.approx(1.0e-5)                               # I=1, d=20mm -> 1e-5 N/m
    assert F > 0                                                    # like currents (I1 I2>0) attract
    assert two_wire_force_per_length(I, -I, d) < 0                  # opposite currents repel (sign flip)
    assert two_wire_force_per_length(2 * I, I, d) == pytest.approx(2 * F)   # linear in each current
    assert two_wire_force_per_length(I, I, 2 * d) == pytest.approx(F / 2)   # 1/d falloff
    # image-wire force is the special case I1=I2=I at separation 2h
    h = d / 2
    assert image_force_wire_iron(I, h) == pytest.approx(two_wire_force_per_length(I, I, 2 * h))
