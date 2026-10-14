"""Electrostatic Maxwell stress / traction helpers."""

import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    EPS0,
    electrostatic_stress_tensor,
    electrostatic_traction,
    electrostatic_traction_summary,
)


def test_electrostatic_normal_field_pressure():
    field = 1.0e6
    pressure = 0.5 * EPS0 * field * field
    tensor = electrostatic_stress_tensor((0.0, 0.0, field))
    traction = electrostatic_traction((0.0, 0.0, field), (0.0, 0.0, 1.0))
    summary = electrostatic_traction_summary(
        (0.0, 0.0, field),
        (0.0, 0.0, 2.0),
        area_m2=0.25,
    )

    assert tensor[0][0] == pytest.approx(-pressure)
    assert tensor[1][1] == pytest.approx(-pressure)
    assert tensor[2][2] == pytest.approx(pressure)
    assert traction == pytest.approx([0.0, 0.0, pressure])
    assert summary["normal_traction_Pa"] == pytest.approx(pressure)
    assert summary["normal_traction_identity_Pa"] == pytest.approx(pressure)
    assert summary["force_N"] == pytest.approx([0.0, 0.0, 0.25 * pressure])


def test_electrostatic_tangential_field_is_tension():
    field = 2.0e6
    pressure = 0.5 * EPS0 * field * field
    summary = electrostatic_traction_summary((field, 0.0, 0.0), (0.0, 0.0, 1.0))

    assert summary["E_normal_V_per_m"] == pytest.approx(0.0)
    assert summary["E_tangent_V_per_m"] == pytest.approx(field)
    assert summary["normal_traction_Pa"] == pytest.approx(-pressure)
    assert summary["tangential_traction_magnitude_Pa"] == pytest.approx(0.0)


def test_electrostatic_oblique_field_decomposes():
    summary = electrostatic_traction_summary((3.0e6, 4.0e6, 0.0), (1.0, 0.0, 0.0))

    assert summary["E_normal_V_per_m"] == pytest.approx(3.0e6)
    assert summary["E_tangent_V_per_m"] == pytest.approx(4.0e6)
    assert summary["normal_traction_Pa"] == pytest.approx(0.5 * EPS0 * (9.0e12 - 16.0e12))
    assert summary["tangential_traction_Pa"] == pytest.approx([0.0, EPS0 * 12.0e12, 0.0])


def test_electrostatic_traction_rejects_bad_inputs():
    with pytest.raises(ValueError):
        electrostatic_stress_tensor((1.0,), eps=EPS0)
    with pytest.raises(ValueError):
        electrostatic_stress_tensor((1.0, 0.0), eps=0.0)
    with pytest.raises(ValueError):
        electrostatic_traction((1.0, 0.0), (0.0, 0.0))
    with pytest.raises(ValueError):
        electrostatic_traction_summary((1.0, 0.0), (1.0, 0.0), area_m2=-1.0)
