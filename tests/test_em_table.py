import numpy as np
import pytest

from radia.panels.em_table import EMTable, interp_Zs, interp_qsurf


@pytest.fixture
def em_table():
    return EMTable(
        path="memory",
        H_grid=np.array([10.0, 100.0]),
        T_grid=np.array([20.0, 120.0]),
        Zs_re=np.array([[1.0, 2.0], [3.0, 4.0]]),
        Zs_im=np.array([[5.0, 6.0], [7.0, 8.0]]),
        q_surf=np.array([[11.0, 12.0], [31.0, 32.0]]),
        meta={},
    )


def test_qsurf_zero_and_signed_field_amplitudes(em_table):
    query = np.array([-100.0, -10.0, 0.0, 10.0, 100.0])
    q = interp_qsurf(em_table, query, 20.0)

    np.testing.assert_allclose(q, [31.0, 11.0, 0.0, 11.0, 31.0])
    np.testing.assert_allclose(
        interp_Zs(em_table, -query, 20.0), interp_Zs(em_table, query, 20.0)
    )


def test_qsurf_broadcasts_field_and_temperature(em_table):
    q = interp_qsurf(
        em_table,
        np.array([[-100.0], [0.0], [100.0]]),
        np.array([[20.0, 120.0]]),
    )

    assert q.shape == (3, 2)
    np.testing.assert_allclose(q, [[31.0, 32.0], [0.0, 0.0], [31.0, 32.0]])


@pytest.mark.parametrize(
    "H_t,T_celsius",
    [(np.nan, 20.0), (np.inf, 20.0), (10.0, np.nan), (10.0, np.inf)],
)
def test_em_table_interpolation_rejects_non_finite_queries(
    em_table, H_t, T_celsius
):
    with pytest.raises(ValueError, match="must be finite"):
        interp_qsurf(em_table, H_t, T_celsius)
    with pytest.raises(ValueError, match="must be finite"):
        interp_Zs(em_table, H_t, T_celsius)


def test_user_impedance_model_broadcast_and_peak_power():
    from radia.panels.em_table import SurfaceImpedanceModel
    model = SurfaceImpedanceModel(lambda H, T: (1 + 0.01 * T + 0.001 * H) * (1+2j),
                                  150000, (0, 1000), (0, 800), "synthetic-v1")
    H = np.array([[0.0], [100.0]])
    T = np.array([20.0, 100.0])
    z = model.impedance(H, T, frequency_hz=150000)
    assert z.shape == (2, 2)
    np.testing.assert_allclose(model.heat_flux(H, T, frequency_hz=150000),
                               0.5 * z.real * H**2)
    with pytest.raises(ValueError, match="frequency"):
        model.impedance(100, 20, frequency_hz=50000)
    with pytest.raises(ValueError, match="domain"):
        model.impedance(-1, 20, frequency_hz=150000)
    with pytest.raises(ValueError, match="domain"):
        model.impedance(100, 801, frequency_hz=150000)


@pytest.mark.parametrize("output", [-1+1j, complex(float('nan'), 0), [1, 2, 3]])
def test_user_impedance_rejects_invalid_output(output):
    from radia.panels.em_table import SurfaceImpedanceModel
    model = SurfaceImpedanceModel(lambda H, T: output, 1000, (0, 100), (0, 100), "bad")
    with pytest.raises(ValueError):
        model.impedance([1, 2], 20, frequency_hz=1000)


def test_user_table_impedance_is_strict_and_energy_consistent(em_table):
    from radia.panels.em_table import SurfaceImpedanceModel
    em_table.meta['frequency'] = 1000
    model = SurfaceImpedanceModel.from_table(em_table, name="synthetic-table")
    H = np.sqrt(1000.0)
    z = model.impedance(H, 70, frequency_hz=1000)
    assert z == pytest.approx(2.5+6.5j)
    assert model.heat_flux(H, 70, frequency_hz=1000) == pytest.approx(1250)
    em_table.Zs_re[:] = -1
    assert model.impedance(H, 70, frequency_hz=1000) == pytest.approx(z)
    with pytest.raises(ValueError, match="domain"):
        model.impedance(1, 70, frequency_hz=1000)
    with pytest.raises(ValueError, match="passivity"):
        SurfaceImpedanceModel.from_table(em_table, name="invalid")
