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
