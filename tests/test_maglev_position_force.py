import numpy as np
import pytest

from radia.maglev import PositionForceCurve


def test_position_force_curve_interpolates_equilibrium_and_compares():
    candidate = PositionForceCurve(
        positions_m=np.array([0.0, 1.0e-3, 2.0e-3]),
        force_N=np.array([-3.0, -2.0, -1.0]),
        name="candidate",
    )
    reference = PositionForceCurve(
        positions_m=np.array([0.0, 1.0e-3, 2.0e-3]),
        force_N=np.array([-3.0, -2.1, -1.0]),
        name="reference",
    )

    assert candidate.at(0.5e-3) == pytest.approx(-2.5)
    np.testing.assert_allclose(candidate.crossings(-1.5), [1.5e-3])
    comparison = candidate.compare(reference)
    assert comparison["sample_count"] == 3
    assert comparison["max_abs_error_N"] == pytest.approx(0.1)
    assert comparison["max_abs_error_normalized"] == pytest.approx(0.1 / 3.0)

    result = candidate.force_result_at(0.5e-3, lift_axis=2)
    assert result["schema"] == "radia.force-result/v1"
    assert result["force_N"] == [0.0, 0.0, -2.5]
    assert result["method"] == "interpolated_lorentz_force"
    assert result["position_m"] == pytest.approx(0.5e-3)
    assert result["force_curve"] == "candidate"


def test_position_force_curve_forbids_extrapolation():
    curve = PositionForceCurve(
        positions_m=np.array([0.0, 1.0]),
        force_N=np.array([0.0, 1.0]),
    )
    with pytest.raises(ValueError, match="outside"):
        curve.at(1.1)
