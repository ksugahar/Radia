"""Fast regression for radia.acoustics analytic sphere scattering.

Physical-limit goldens (self-validating, no MATLAB / ngsolve needed) plus a small
stored golden vector.  The heavy 3-way cross-check (MATLAB golden + ngsolve.bem
numerical) lives in validation_test/acoustics/test_scattering_3way.py.
"""
import numpy as np
import pytest

import radia.acoustics as ac

K = 2.0
R = 1.0


def _shell(radius, n=7):
    th = np.linspace(0.15, np.pi - 0.15, n)
    return np.c_[radius * np.sin(th), np.zeros(n), radius * np.cos(th)]


def _axis(zs):
    zs = np.asarray(zs, float)
    return np.c_[np.zeros_like(zs), np.zeros_like(zs), zs]


def test_soft_surface_pressure_vanishes():
    # sound-soft boundary condition: total pressure = 0 on r = R
    total = ac.soft_sphere_scattering(K, R, _shell(R))["total"]
    assert np.max(np.abs(total)) < 1e-9


def test_fluid_invisible_when_matched():
    # k1 = k0 and density_ratio = 1 -> acoustically invisible: total == incident
    f = ac.fluid_sphere_scattering(K, R, _axis([1.2, 1.5, 2.0, 3.0, 5.0]),
                                   interior_wavenumber=K, density_ratio=1.0)
    assert np.max(np.abs(f["total"] - f["incident"])) < 1e-5


def test_elastic_shear_zero_recovers_fluid():
    # Faran with zero shear = Anderson fluid sphere (interior wavenumber k0/cL)
    cL, rho = 2.0, 1.5
    pts = _axis([1.2, 1.5, 2.0, 3.0, 5.0])
    el = ac.elastic_sphere_scattering(K, R, pts, longitudinal_speed=cL,
                                      shear_speed=0.0, density_ratio=rho)["total"]
    fl = ac.fluid_sphere_scattering(K, R, pts, interior_wavenumber=K / cL,
                                    density_ratio=rho)["total"]
    assert np.max(np.abs(el - fl)) / np.max(np.abs(fl)) < 1e-4


def test_elastic_stiff_recovers_rigid():
    # very stiff / heavy elastic sphere -> sound-hard (rigid) limit (~1% at cL=50)
    pts = _axis([1.2, 1.5, 2.0, 3.0, 5.0])
    el = ac.elastic_sphere_scattering(K, R, pts, longitudinal_speed=50.0,
                                      shear_speed=25.0, density_ratio=50.0)["total"]
    rg = ac.rigid_sphere_scattering(K, R, pts)["total"]
    assert np.max(np.abs(el - rg)) / np.max(np.abs(rg)) < 0.03


def test_stored_golden_scattered_at_far_point():
    # locked against the machine-precision-verified MATLAB reference (~1e-14)
    p = np.array([[0.0, 0.0, 3.0]])
    np.testing.assert_allclose(
        ac.soft_sphere_scattering(K, R, p)["scattered"][0],
        -0.378514629643 + 0.421040785473j, atol=1e-9)
    np.testing.assert_allclose(
        ac.rigid_sphere_scattering(K, R, p)["scattered"][0],
        0.144643955105 + 0.214559921779j, atol=1e-9)
    np.testing.assert_allclose(
        ac.elastic_sphere_scattering(K, R, p, longitudinal_speed=2.0,
                                     shear_speed=1.0, density_ratio=1.5)["scattered"][0],
        -0.0867709948317 + 1.02182418537j, atol=1e-9)


@pytest.mark.parametrize("fn", [ac.soft_sphere_scattering,
                                ac.rigid_sphere_scattering,
                                ac.elastic_sphere_scattering])
def test_exterior_only_raises_on_interior(fn):
    # exterior scattered-field references: interior points are rejected
    with pytest.raises(ValueError):
        fn(K, R, np.array([[0.0, 0.0, 0.3]]))
