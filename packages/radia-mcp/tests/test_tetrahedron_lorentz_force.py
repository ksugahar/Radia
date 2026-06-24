"""P1 tetrahedron Lorentz body-force loads."""

import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.force import tetrahedron_lorentz_force_summary  # noqa: E402


UNIT_TET = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
]


def test_tetrahedron_lorentz_force_integrates_constant_j_cross_b():
    row = tetrahedron_lorentz_force_summary(UNIT_TET, (2.0, 0.0, 0.0), (0.0, 3.0, 0.0))

    assert row["volume_m3"] == pytest.approx(1.0 / 6.0)
    assert row["force_density_N_per_m3"] == pytest.approx([0.0, 0.0, 6.0])
    assert row["integrated_force_N"] == pytest.approx([0.0, 0.0, 1.0])
    assert row["p1_shape_function_integral_m3"] == pytest.approx(1.0 / 24.0)
    for node_force in row["nodal_force_loads_N"]:
        assert node_force == pytest.approx([0.0, 0.0, 0.25])


def test_tetrahedron_lorentz_force_orientation_independent_volume():
    row = tetrahedron_lorentz_force_summary(list(reversed(UNIT_TET)), (0.0, 4.0, 0.0), (0.0, 0.0, 5.0))

    assert row["volume_m3"] == pytest.approx(1.0 / 6.0)
    assert row["force_density_N_per_m3"] == pytest.approx([20.0, 0.0, 0.0])
    assert row["integrated_force_N"] == pytest.approx([20.0 / 6.0, 0.0, 0.0])


def test_tetrahedron_lorentz_force_rejects_bad_inputs():
    with pytest.raises(ValueError):
        tetrahedron_lorentz_force_summary(UNIT_TET[:3], (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    with pytest.raises(ValueError):
        tetrahedron_lorentz_force_summary(UNIT_TET, (1.0, 0.0), (0.0, 1.0, 0.0))
