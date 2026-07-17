import numpy as np
import pytest

from radia.maglev import MovingHCurlCLNFamily, PositionForceCurve
from radia.vim import HCurlEddyCLNModel


def _model(resistance, inductance, port):
    return HCurlEddyCLNModel(
        resistance=np.diag(resistance),
        inductance=np.diag(inductance),
        surface_mass=np.zeros((2, 2)),
        port_rhs=np.asarray(port, dtype=float).reshape(2, 1),
        basis_names=("cycle", "bulk"),
        blocks={"bridge": (0, 1), "volume": (1, 2)},
    )


def test_moving_hcurl_cln_family_interpolates_constant_basis_passively():
    left = _model([2.0, 4.0], [3.0, 5.0], [1.0, 2.0])
    right = _model([4.0, 8.0], [5.0, 9.0], [3.0, 6.0])
    family = MovingHCurlCLNFamily(
        positions_m=np.array([0.0, 0.01]),
        models=(left, right),
    )

    middle = family.at(0.005)
    np.testing.assert_allclose(middle.resistance, np.diag([3.0, 6.0]))
    np.testing.assert_allclose(middle.inductance, np.diag([4.0, 7.0]))
    np.testing.assert_allclose(middle.port_rhs[:, 0], [2.0, 4.0])
    assert middle.diagnostics()["passive"] is True
    assert family.diagnostics() == {
        "position_samples": 2,
        "position_min_m": 0.0,
        "position_max_m": 0.01,
        "state_order": 2,
        "port_count": 1,
        "constant_basis": True,
        "all_samples_passive": True,
        "all_samples_finite_rl": True,
    }


def test_moving_hcurl_cln_family_rejects_nonconstant_coordinates():
    left = _model([2.0, 4.0], [3.0, 5.0], [1.0, 2.0])
    right = HCurlEddyCLNModel(
        resistance=np.eye(2),
        inductance=np.eye(2),
        surface_mass=np.zeros((2, 2)),
        port_rhs=np.ones((2, 1)),
        basis_names=("different", "basis"),
    )

    with pytest.raises(ValueError, match="basis names"):
        MovingHCurlCLNFamily(np.array([0.0, 1.0]), (left, right))


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


def test_position_force_curve_forbids_extrapolation():
    curve = PositionForceCurve(
        positions_m=np.array([0.0, 1.0]),
        force_N=np.array([0.0, 1.0]),
    )
    with pytest.raises(ValueError, match="outside"):
        curve.at(1.1)
