import numpy as np
import pytest


native = pytest.importorskip("radia._radia_pybind")
from radia.beam import propagate_variational_map


def _zero_jets(count):
    return (
        np.zeros((count, 6, 6), dtype=float),
        np.zeros((count, 6, 6, 6), dtype=float),
        np.zeros((count, 6, 6, 6, 6), dtype=float),
    )


def test_native_variational_map_attributes_ordered_nonlinear_regions():
    lengths = np.array([0.3, 0.4, 0.2])
    a, f2, f3 = _zero_jets(3)
    f2[0, 1, 0, 0] = 2.0
    f2[1, 2, 0, 1] = 1.5
    f2[1, 2, 1, 0] = 1.5
    f3[2, 3, 0, 0, 0] = -0.7

    report = propagate_variational_map(
        lengths,
        a,
        f2,
        f3,
        ["upstream_sextupole", "downstream_sextupole", "direct_octupole"],
        maximum_order=3,
        maximum_step_m=1.0,
    )

    assert report["backend"] == "native-cpp"
    assert report["R"].shape == (6, 6)
    assert report["T"].shape == (6, 6, 6)
    assert report["U"].shape == (6, 6, 6, 6)
    assert report["station_R"].shape == (4, 6, 6)
    assert report["region_T"].shape == (3, 6, 6, 6)
    np.testing.assert_allclose(report["T"][1, 0, 0], 0.6, atol=2e-14)
    np.testing.assert_allclose(report["U"][2, 0, 0, 0], 1.08, atol=3e-14)
    np.testing.assert_allclose(report["U"][3, 0, 0, 0], -0.14, atol=2e-14)
    np.testing.assert_array_equal(report["pair_regions"], [[0, 1]])
    np.testing.assert_allclose(
        report["pair_U_cascade"][0, 2, 0, 0, 0], 1.08, atol=3e-14
    )
    assert report["diagnostics"]["T_reconstruction_error"] < 3e-14
    assert report["diagnostics"]["U_reconstruction_error"] < 4e-14
    assert report["diagnostics"]["T_symmetry_defect"] < 1e-14
    assert report["diagnostics"]["U_symmetry_defect"] < 1e-14


def test_native_variational_map_matches_normal_quadrupole_matrix():
    strength = 1.7
    length = 0.8
    root = np.sqrt(strength)
    a, _, _ = _zero_jets(1)
    a[0, 0, 1] = 1.0
    a[0, 1, 0] = -strength
    a[0, 2, 3] = 1.0
    a[0, 3, 2] = strength

    report = propagate_variational_map(
        [length], a, maximum_order=1, maximum_step_m=0.001
    )
    phase = root * length
    expected_x = np.array(
        [
            [np.cos(phase), np.sin(phase) / root],
            [-root * np.sin(phase), np.cos(phase)],
        ]
    )
    expected_y = np.array(
        [
            [np.cosh(phase), np.sinh(phase) / root],
            [root * np.sinh(phase), np.cosh(phase)],
        ]
    )
    np.testing.assert_allclose(report["R"][:2, :2], expected_x, atol=8e-13)
    np.testing.assert_allclose(report["R"][2:4, 2:4], expected_y, atol=8e-13)
    assert report["diagnostics"]["R_composition_error"] < 1e-12


def test_native_variational_map_rejects_nonsymmetric_f2():
    a, f2, _ = _zero_jets(1)
    f2[0, 0, 0, 1] = 1.0
    with pytest.raises(ValueError, match="F2 input indices must be symmetric"):
        propagate_variational_map([1.0], a, f2, maximum_order=2)


def test_native_variational_map_keeps_substep_cascade_in_one_region():
    a, f2, _ = _zero_jets(1)
    f2[0, 1, 0, 0] = 2.0
    f2[0, 2, 0, 1] = 1.5
    f2[0, 2, 1, 0] = 1.5

    report = propagate_variational_map(
        [1.0], a, f2, maximum_order=3, maximum_step_m=0.1
    )

    np.testing.assert_allclose(report["U"][2, 0, 0, 0], 4.5, atol=2e-13)
    np.testing.assert_allclose(
        report["region_U_local_cascade"][0, 2, 0, 0, 0], 4.5, atol=2e-13
    )
    assert report["pair_regions"].shape == (0, 2)
    assert report["pair_U_cascade"].shape == (0, 6, 6, 6, 6)
    assert report["diagnostics"]["U_reconstruction_error"] < 3e-13
