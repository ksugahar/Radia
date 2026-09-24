"""Relative tolerance of the adaptive arc-current B quadrature ("PrcArc")."""
import numpy as np
import pytest
import radia as rad

RI, RO, HEIGHT, DENSITY = .035, .070, .105, 1e6
ANGLES = (0.2, 1.8)
# Exterior points only: the tensor Gauss reference is not accurate inside the
# current region, where the integrand is singular.
POINTS = [(.012, -.009, .065), (.09, .01, .08), (.05, .03, .07), (.12, .08, .02)]


@pytest.fixture(autouse=True)
def default_arc_tolerance():
    rad.FldLenRndSw('off')
    rad.FldCmpPrc('PrcArc->1e-9')
    yield
    rad.FldCmpPrc('PrcArc->1e-9')
    rad.FldLenRndSw('on')


def tensor_gauss_reference(point, order=64):
    """Independent volume integral of J e_phi x (x - x') / |x - x'|^3."""
    x, w = np.polynomial.legendre.leggauss(order)
    rr = ((RI+RO)/2+(RO-RI)/2*x)[:, None, None]
    zz = (HEIGHT/2*x)[None, :, None]
    phi = ((ANGLES[0]+ANGLES[1])/2+(ANGLES[1]-ANGLES[0])/2*x)[None, None, :]
    dx, dy, dz = point[0]-rr*np.cos(phi), point[1]-rr*np.sin(phi), point[2]-zz
    inv = (dx*dx+dy*dy+dz*dz)**-1.5
    weights = w[:, None, None]*w[None, :, None]*w[None, None, :]
    factor = 1e-7*DENSITY*(RO-RI)*HEIGHT*(ANGLES[1]-ANGLES[0])/8
    return factor*np.array([np.sum(weights*rr*np.cos(phi)*dz*inv),
                            np.sum(weights*rr*np.sin(phi)*dz*inv),
                            np.sum(weights*rr*(rr-point[0]*np.cos(phi)-point[1]*np.sin(phi))*inv)])


def arc_field(point):
    obj = rad.ObjArcCur([0, 0, 0], [RI, RO], list(ANGLES), HEIGHT, 4, 'man', 'z', DENSITY)
    return np.asarray(rad.Fld(obj, 'b', list(point)))


def test_explicit_default_is_bit_identical_to_the_default():
    rad.UtiDelAll()
    implicit = [arc_field(p) for p in POINTS]
    rad.FldCmpPrc('PrcArc->1e-9')
    explicit = [arc_field(p) for p in POINTS]
    for a, b in zip(implicit, explicit):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize('rtol', [1e-7, 1e-5, 1e-3])
def test_relaxed_tolerance_bounds_the_error_against_an_independent_reference(rtol):
    rad.FldCmpPrc(f'PrcArc->{rtol:g}')
    for point in POINTS:
        reference = tensor_gauss_reference(point)
        error = np.linalg.norm(arc_field(point)-reference)/np.linalg.norm(reference)
        # The reference itself agrees with order 40 to ~1e-8; allow the
        # requested tolerance with a modest factor for the vector norm.
        assert error <= 10*rtol + 2e-8, (point, error)


def test_setting_persists_across_model_deletion_until_reset():
    rad.FldCmpPrc('PrcArc->1e-3')
    loose = arc_field(POINTS[0])
    rad.UtiDelAll()
    assert np.array_equal(arc_field(POINTS[0]), loose)
    rad.FldCmpPrc('PrcArc->1e-9')
    assert not np.array_equal(arc_field(POINTS[0]), loose)


@pytest.mark.parametrize('value', ['1e-2', '1e-13', '0', 'nan', 'inf', '-1e-6'])
def test_out_of_range_tolerance_fails_loudly(value):
    before = arc_field(POINTS[0])
    with pytest.raises(RuntimeError):
        rad.FldCmpPrc(f'PrcArc->{value}')
    np.testing.assert_array_equal(arc_field(POINTS[0]), before)
