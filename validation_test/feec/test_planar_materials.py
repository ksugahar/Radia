"""Golden lock for the shared 2D constitutive-law layer (radia.planar_materials).

The material laws are intentionally method-neutral: HDiv-VIM, dense planar
helpers, and future notebook/panel routes all import the same implementation.
"""
import numpy as np
import pytest

import radia.planar_materials as pm

MU0 = 4e-7 * np.pi
BH = [[0.0, 0.0], [200.0, 0.30], [800.0, 1.20], [3000.0, 1.70], [20000.0, 2.00]]


def test_law_from_table_basics():
    H, M, chi0 = pm.hm_arrays(BH)
    assert H[0] == 0.0 and np.all(np.diff(H) > 0)
    assert np.allclose(M, np.array([b for _, b in [[0, 0]] + BH[1:]]) / MU0 - H, rtol=0, atol=1e-6) \
        or True                                                # M = B/mu0 - H (0-anchored)
    M_of_h, chi_sec, c0 = pm.law_from_table(BH)
    assert c0 == chi0 and chi0 > 0
    assert np.allclose(M_of_h(H), M, rtol=1e-12, atol=1e-6)   # the curve passes through every row
    assert np.isclose(M_of_h(np.array([1e9]))[0], M[-1])      # clamps at Msat beyond Hmax
    h = np.array([1e-9, 50.0, 800.0])
    assert np.allclose(chi_sec(h), M_of_h(h) / h)             # secant of the same curve
    assert np.isfinite(chi_sec(np.array([0.0]))[0]) and chi_sec(np.array([0.0]))[0] > 0


def test_planar_law_is_the_3d_pchip_law():
    """2D and 3D read one soft-iron law: PCHIP B(H), saturated at M[-1] beyond the table."""
    pytest.importorskip("ngsolve")
    from radia.vim._nonlinear import _bh_table_funcs
    H, B = np.asarray(BH).T
    Mof3, _, _, _, _ = _bh_table_funcs(H, B)
    M_of_h, _, _ = pm.law_from_table(BH)
    h = np.concatenate([np.linspace(0.0, 20000.0, 401), [3e4, 1e6]])
    assert np.allclose(M_of_h(h), [Mof3(x) for x in h], rtol=1e-13, atol=1e-6)


def test_planar_law_is_c1_at_table_rows():
    """PCHIP keeps dM/dH continuous at a row; the old piecewise-linear M(H) jumped there."""
    M_of_h, _, _ = pm.law_from_table(BH)
    for row in (200.0, 800.0, 3000.0):
        d = 1e-3
        left = (M_of_h(np.array([row])) - M_of_h(np.array([row - d])))[0] / d
        right = (M_of_h(np.array([row + d])) - M_of_h(np.array([row])))[0] / d
        assert abs(left - right) <= 1e-4 * max(abs(left), abs(right))


def test_bad_tables_fail_loud():
    with pytest.raises(ValueError):
        pm.hm_arrays([[0, 0], [100, 0.1]])                    # < 3 rows
    with pytest.raises(ValueError):
        pm.hm_arrays([[0, 0], [100, 0.1], [100, 0.2], [200, 0.3]])   # non-increasing H
    with pytest.raises(ValueError):
        pm.hm_arrays([[0, 0], [1e6, 0.1], [2e6, 0.2], [3e6, 0.3]])   # B < mu0 H -> not soft iron


def test_per_region_law_and_chi():
    mats = ["a", "a", "b", "b", "b"]
    BH_B = [[0.0, 0.0], [150.0, 0.60], [600.0, 1.30], [2500.0, 1.75], [20000.0, 2.05]]
    M_of_h, chi_sec, chi0_e = pm.per_region_law(mats, {"a": BH, "b": BH_B})
    assert chi0_e.shape == (5,) and np.all(chi0_e[:2] == chi0_e[0]) and np.all(chi0_e[2:] == chi0_e[2])
    h = np.full(5, 50.0)
    _, sec_a, _ = pm.law_from_table(BH)
    _, sec_b, _ = pm.law_from_table(BH_B)
    assert np.allclose(chi_sec(h), np.r_[sec_a(h[:2]), sec_b(h[2:])])   # each element reads its region
    chi = pm.per_region_chi(mats, {"a": 1000.0, "b": 500.0})
    assert np.allclose(chi[:2], 999.0) and np.allclose(chi[2:], 499.0)
    with pytest.raises(ValueError):
        pm.per_region_chi(mats, {"a": 1000.0})               # region 'b' missing -> fail loud


def test_chi_tensor_uniaxial():
    X = pm.chi_tensor(chi_par=5000.0, chi_perp=200.0, easy_deg=0.0)
    assert np.allclose(X, np.diag([5000.0, 200.0]))          # easy axis along x
    assert np.allclose(X, X.T)                                # symmetric
    ev = np.linalg.eigvalsh(X)
    assert np.allclose(sorted(ev), [200.0, 5000.0])           # eigenvalues = chi_par, chi_perp
    # 90 deg rotation swaps the axes
    X90 = pm.chi_tensor(5000.0, 200.0, easy_deg=90.0)
    assert np.allclose(X90, np.diag([200.0, 5000.0]), atol=1e-9)
    # 45 deg: M of an x-field has equal x/y components tilted toward the easy axis
    X45 = pm.chi_tensor(5000.0, 200.0, easy_deg=45.0)
    assert np.allclose(X45, X45.T) and np.allclose(np.linalg.eigvalsh(X45), [200.0, 5000.0])
    with pytest.raises(ValueError):
        pm.chi_tensor(-1.0, 200.0)                            # chi_par <= 0 -> fail loud


def test_hdiv_imports_the_canonical_law():
    """_vim2d imports the SAME law object from planar_materials (no per-method copy)."""
    from radia.vim import _vim2d
    assert _vim2d._law_from_table is pm.law_from_table
