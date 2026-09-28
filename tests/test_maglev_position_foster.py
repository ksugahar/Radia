import numpy as np
import pytest

from radia.maglev import MovingHCurlFosterFamily
from radia.vim import HCurlEddyFosterModel


def _model(port, resistance=(2.0, 4.0), inductance=(3.0, 5.0)):
    return HCurlEddyFosterModel(
        resistance=np.diag(resistance),
        inductance=np.diag(inductance),
        surface_mass=np.zeros((2, 2)),
        port_rhs=np.asarray(port, dtype=float).reshape(2, 1),
        basis_names=("cycle", "bulk"),
        blocks={"bridge": (0, 1), "volume": (1, 2)},
    )


def test_moving_foster_family_interpolates_ports_with_shared_modes():
    left = _model([1.0, 2.0])
    right = _model([3.0, 6.0])
    family = MovingHCurlFosterFamily(positions_m=np.array([0.0, 0.01]), models=(left, right))

    middle = family.at(0.005)
    np.testing.assert_allclose(middle.port_rhs[:, 0], [2.0, 4.0])
    np.testing.assert_allclose(middle.decay_rates, left.decay_rates)
    s = 2j * np.pi * 50.0
    direct = np.linalg.solve(np.diag([2.0, 4.0]) + s * np.diag([3.0, 5.0]),
                             np.array([2.0, 4.0]))
    np.testing.assert_allclose(middle.solve(s, 1.0), direct, rtol=1e-12)
    assert family.diagnostics() == {
        "position_samples": 2,
        "position_min_m": 0.0,
        "position_max_m": 0.01,
        "state_order": 2,
        "port_count": 1,
        "shared_modes": True,
        "all_samples_passive": True,
        "all_samples_finite_rl": True,
    }
    with pytest.raises(ValueError, match="outside"):
        family.at(0.02)


def test_moving_foster_family_rejects_position_dependent_operators():
    left = _model([1.0, 2.0])
    right = _model([1.0, 2.0], resistance=(2.5, 4.0))
    with pytest.raises(ValueError, match="shared R/L"):
        MovingHCurlFosterFamily(np.array([0.0, 1.0]), (left, right))


def test_moving_foster_family_rejects_different_coordinates():
    left = _model([1.0, 2.0])
    right = HCurlEddyFosterModel(
        resistance=np.diag([2.0, 4.0]),
        inductance=np.diag([3.0, 5.0]),
        surface_mass=np.zeros((2, 2)),
        port_rhs=np.ones((2, 1)),
        basis_names=("different", "basis"),
    )
    with pytest.raises(ValueError, match="basis names"):
        MovingHCurlFosterFamily(np.array([0.0, 1.0]), (left, right))
