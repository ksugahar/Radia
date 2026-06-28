"""
Test material API:
- MatLin(mu_r) - isotropic linear material
- MatLin([mu_r_par, mu_r_perp], [ex, ey, ez]) - anisotropic linear material
- MatPM(Br, Hc, [mx, my, mz]) - permanent magnet with demagnetization
"""

import pytest


@pytest.mark.basic
def test_isotropic_linear_material(radia_clean):
    """Test MatLin(mu_r) for isotropic soft iron."""
    rad = radia_clean

    mu_r = 1000
    mat_iso = rad.MatLin(mu_r)
    assert isinstance(mat_iso, int)
    assert mat_iso > 0

    # Apply to a block
    block = rad.magnet_box([0, 0, 0], [10, 10, 10], [0, 0, 0])
    rad.MatApl(block, mat_iso)


@pytest.mark.basic
def test_anisotropic_linear_material(radia_clean):
    """Test MatLin([mu_r_par, mu_r_perp], [ex, ey, ez]) for anisotropic material."""
    rad = radia_clean

    mu_r_par = 1.06
    mu_r_perp = 1.17
    easy_axis = [0, 0, 1]
    mat_aniso = rad.MatLin([mu_r_par, mu_r_perp], easy_axis)
    assert isinstance(mat_aniso, int)
    assert mat_aniso > 0

    # Apply to a block
    block = rad.magnet_box([20, 0, 0], [10, 10, 10], [0, 0, 0])
    rad.MatApl(block, mat_aniso)


@pytest.mark.basic
def test_permanent_magnet_material(radia_clean):
    """Test MatPM(Br, Hc, [mx, my, mz]) for NdFeB N52."""
    rad = radia_clean

    Br = 1.43  # Tesla
    Hc = 876000  # A/m
    mag_axis = [0, 0, 1]
    mat_pm = rad.MatPM(Br, Hc, mag_axis)
    assert isinstance(mat_pm, int)
    assert mat_pm > 0

    # Apply to a block and verify material assignment succeeds
    block = rad.magnet_box([0, 0, 0], [20, 20, 10], [0, 0, 0])
    rad.MatApl(block, mat_pm)
