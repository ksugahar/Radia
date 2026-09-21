"""Independent volume-current checks for the explicit loft filament model."""
import numpy as np
import pytest
import radia as rad
from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile


def reference(point, order=24):
    # Independent volume Biot-Savart integration of the prescribed stream map.
    x, w = np.polynomial.legendre.leggauss(order)
    a, b, s = np.meshgrid(x / 2, x / 2, (x + 1) / 2, indexing='ij')
    weights = np.einsum('i,j,k->ijk', w / 2, w / 2, w / 2)
    source = np.stack((a * (4 + 4*s), 20*s, b * (6 + 4*s)), axis=-1)
    tangent = np.stack((4*a, np.full_like(a, 20), 4*b), axis=-1)
    delta = np.asarray(point) - source
    integrand = np.cross(tangent, delta) / np.linalg.norm(delta, axis=-1)[..., None]**3
    return 1e-7 * 100 * np.sum(integrand * weights[..., None], axis=(0, 1, 2))


def coil(current=100):
    return CoilBuilder(current).set_cross_section(4, 6).add_loft_straight(RectProfile(8, 10), 20)


@pytest.mark.parametrize('point', [[30, 10, 20], [-12, 4, 7], [2, 30, 15]])
def test_tapered_volume_reference_and_section_convergence(point):
    exact = reference(point)
    np.testing.assert_allclose(exact, reference(point, 32), rtol=1e-10, atol=1e-16)
    errors = []
    for n in (4, 8, 16, 32):
        value = rad.Fld(rad.ObjCnt(coil().to_radia_loft_filaments(n, n)), 'b', point)
        errors.append(np.linalg.norm(value - exact) / np.linalg.norm(exact))
    assert all(b < a / 3 for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 1e-4


def test_zero_taper_matches_solid_limit():
    c = CoilBuilder(100).set_cross_section(4, 6).add_straight(20)
    point = [20, 10, 15]
    exact = rad.Fld(rad.ObjCnt(c.to_radia()), 'b', point)
    value = rad.Fld(rad.ObjCnt(c.to_radia_loft_filaments(32, 32)), 'b', point)
    np.testing.assert_allclose(value, exact, rtol=1e-4, atol=1e-16)


def test_current_is_conserved_and_junction_is_shared(monkeypatch):
    c = coil().add_straight(10)
    calls = []
    monkeypatch.setattr(rad, 'ObjFlmCur', lambda p, i: calls.append((p, i)) or len(calls))
    assert len(c.to_radia_loft_filaments(3, 4)) == 12
    assert sum(i for p, i in calls) == pytest.approx(100)
    assert all(len(p) == 3 for p, i in calls)


def test_reversed_current():
    p = [30, 10, 20]
    plus = rad.Fld(rad.ObjCnt(coil().to_radia_loft_filaments(8, 8)), 'b', p)
    minus = rad.Fld(rad.ObjCnt(coil(-100).to_radia_loft_filaments(8, 8)), 'b', p)
    np.testing.assert_allclose(minus, -np.asarray(plus), rtol=1e-13, atol=1e-16)


@pytest.mark.parametrize('invalid', ['gap', 'arc', 'zero_width', 'current'])
@pytest.mark.parametrize('n', [1, 3])
def test_invalid_geometry_does_not_allocate(monkeypatch, invalid, n):
    c = coil()
    if invalid == 'gap':
        c.set_cross_section(2, 2).add_straight(10)
    elif invalid == 'arc':
        c.add_arc(20, 90)
    elif invalid == 'zero_width':
        c.segments[0].profile_end.w = 0
    else:
        c.segments[0].current = 99
    def fail(*args):
        pytest.fail('allocation before validation')
    monkeypatch.setattr(rad, 'ObjFlmCur', fail)
    with pytest.raises((ValueError, NotImplementedError)):
        c.to_radia_loft_filaments(n, n)
