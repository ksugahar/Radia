import numpy as np
import pytest

from radia import vim
import radia._radia_pybind as _rp


def _toy_material(*, gamma=0.0):
    eta = np.array([0.2, 0.6])
    peaks = [2.0e4, 3.0e4]
    tables = [
        (np.linspace(0.0, radius, 9), np.linspace(0.0, peak, 9))
        for radius, peak in zip(eta, peaks)
    ]
    return vim.EnergyStopMaterial(
        eta, tables, alpha=5.0, gamma=gamma, b_max=1.2
    )


def test_energy_stop_hard_projection_is_batched_pure_and_explicit():
    material = _toy_material()
    states = np.tile(material.state0(), (2, 1))
    B = np.array([[0.8, 0.0, 0.0], [0.0, 0.3, 0.4]])

    H_first = material.forward(B, states)
    H_second = material.forward(B, states)
    committed = material.commit(B, states)

    np.testing.assert_array_equal(H_first, H_second)
    np.testing.assert_allclose(committed[0, :3], [0.2, 0.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(committed[0, 3:6], [0.6, 0.0, 0.0], atol=1e-15)
    np.testing.assert_array_equal(committed[:, -3:], B)
    np.testing.assert_array_equal(states, np.zeros_like(states))
    assert material.state_size == 9
    assert material.nu_bound() == pytest.approx(150005.0)
    assert np.all(material.stored_energy(B, states) >= 0.0)


def test_energy_stop_variational_update_satisfies_stationarity():
    slope = 1.0e5
    gamma = 2.0e-6
    radius = 0.5
    r = np.linspace(0.0, radius, 11)
    material = vim.EnergyStopMaterial(
        [radius], [(r, slope * r)], alpha=5.0, gamma=gamma
    )
    B = np.array([[0.2, 0.0, 0.0]])
    committed = material.commit(B, material.state0()[None, :])
    s_radius = np.linalg.norm(committed[0, :3])

    expected = B[0, 0] / (1.0 + slope * gamma)
    residual = slope * s_radius + (s_radius - B[0, 0]) / gamma
    assert s_radius == pytest.approx(expected, abs=2e-15)
    assert abs(residual) < 1e-8


def test_energy_stop_closed_vector_loop_has_nonnegative_dissipation():
    material = _toy_material()
    theta = np.linspace(0.0, 2.0 * np.pi, 721)
    circle = np.column_stack([
        0.8 * np.cos(theta),
        0.8 * np.sin(theta),
        0.15 * np.sin(2.0 * theta),
    ])
    path = np.vstack([np.linspace(0.0, 1.0, 80)[:, None] * circle[0], circle])
    states = material.state0()[None, :]
    H = []
    for B in path:
        trial = B[None, :]
        H.append(material.forward(trial, states)[0])
        states = material.commit(trial, states)
    H = np.asarray(H)
    dB = np.diff(path, axis=0)
    dissipation = float(np.sum(0.5 * (H[1:] + H[:-1]) * dB))

    assert dissipation > 0.0


@pytest.mark.parametrize(
    "table",
    [
        (np.array([0.0, 0.1, 0.2]), np.array([0.0, 2.0, 1.0])),
        (np.array([0.0, 0.2, 0.1]), np.array([0.0, 1.0, 2.0])),
        (np.array([0.0, 0.1]), np.array([1.0, 2.0])),
    ],
)
def test_energy_stop_rejects_nonconvex_or_malformed_tables(table):
    with pytest.raises((ValueError, RuntimeError), match="EnergyStopMaterial"):
        vim.EnergyStopMaterial([0.2], [table])


def test_energy_stop_cpp_boundary_rejects_out_of_range_offsets():
    with pytest.raises((ValueError, RuntimeError), match="table offsets"):
        _rp._EnergyStopMaterial(
            np.array([0.2, 0.3]),
            np.array([0.0, 0.2, 0.0, 0.3]),
            np.array([0.0, 1.0, 0.0, 1.0]),
            np.array([0, -1, 4], dtype=np.int32),
            np.zeros(2), 5.0, 0.0,
        )
