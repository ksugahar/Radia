"""Independent enthalpy and axisymmetric-radius correctness checks."""

from pathlib import Path
import sys

import numpy as np
import pytest

from radia.ih_thermal_material import ThermalMaterial


def test_default_table_enthalpy_integrates_internal_cp_knot():
    material = ThermalMaterial(
        rho=1.0,
        T=np.array([0.0, 100.0, 200.0]),
        k=np.array([1.0, 1.0, 1.0]),
        cp=np.array([1.0, 10.0, 1.0]),
    )

    # Independent piecewise-linear integral:
    # 100*(1+10)/2 + 100*(10+1)/2 = 1100 J/kg.
    assert float(material.enthalpy(200.0)) == pytest.approx(1100.0)


def test_constant_source_shortcut_requires_a_truly_constant_table():
    material = ThermalMaterial(
        rho=1.0,
        T=np.array([0.0, 100.0, 200.0]),
        k=np.array([1.0, 1.0, 1.0]),
        cp=np.array([1.0, 10.0, 1.0]),
        source="constant",
    )
    assert float(material.enthalpy(200.0)) == pytest.approx(1100.0)


PANELS = Path(__file__).resolve().parents[1] / "src" / "radia" / "panels"
if str(PANELS) not in sys.path:
    sys.path.insert(0, str(PANELS))
import calc_heat_axisym  # noqa: E402


class _Vertex:
    def __init__(self, radius):
        self.point = (radius, 0.0, 0.0)


class _Mesh:
    def __init__(self, radii):
        self.dim = 2
        self.vertices = [_Vertex(radius) for radius in radii]


def test_axisymmetric_radius_validation_accepts_axis_and_positive_radius():
    calc_heat_axisym._validate_axisym_radii(_Mesh([0.0, 0.01, 0.02]))


@pytest.mark.parametrize("radius", [-1.0e-12, np.nan, np.inf])
def test_axisymmetric_radius_validation_rejects_invalid_weight(radius):
    with pytest.raises(ValueError, match=r"r >= 0|finite"):
        calc_heat_axisym._validate_axisym_radii(_Mesh([0.0, radius]))


def test_axisymmetric_solve_rejects_negative_radius_before_assembly():
    result = calc_heat_axisym.solve_heat_axisym(
        "unused.vol", _wp_mesh=_Mesh([0.0, -1.0e-6]),
        _write_solution=False,
    )
    assert "negative radial coordinate" in result["error"]
