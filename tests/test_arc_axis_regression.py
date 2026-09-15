"""Native circular-current axis and azimuth regression, SI units."""
import numpy as np
import pytest
import radia as rad


@pytest.fixture(autouse=True)
def exact_coordinates():
    # Geometry perturbation is a separate policy, not analytic kernel error.
    rad.FldLenRndSw('off')
    yield
    rad.FldLenRndSw('on')


def reference_axis(rin, rout, height, current_density, z):
    x, w = np.polynomial.legendre.leggauss(32)
    radii = (rin+rout)/2 + (rout-rin)/2*x
    heights = height/2*x
    rr = radii[:, None]
    dz = z-heights[None, :]
    kernel = rr**2/(rr**2+dz**2)**1.5
    return 2*np.pi*1e-7*current_density*(rout-rin)*height/4*np.sum(kernel*w[:, None]*w[None, :])


@pytest.mark.parametrize('scale', [0.01, 1., 100.])
def test_full_circle_exact_axis(scale):
    ri, ro, h = np.array([.010,.013,.003])*scale
    obj = rad.ObjArcCur([0,0,0], [ri,ro], [0,2*np.pi], h, 40, 'man','z',1e7)
    for z in np.array([0.,.003,.009])*scale:
        b = np.asarray(rad.Fld(obj,'b',[0.,0.,z]))
        expected = reference_axis(ri,ro,h,1e7,z)
        assert np.linalg.norm(b[:2]) < abs(expected)*1e-12
        assert b[2] == pytest.approx(expected,rel=2e-7)


def test_full_circle_azimuth_covariance():
    obj = rad.ObjArcCur([0,0,0],[.010,.013],[0,2*np.pi],.003,40,'man','z',1e7)
    radius, z = .002,.006
    base = np.asarray(rad.Fld(obj,'b',[radius,0,z]))
    for theta in [np.pi/2,np.pi,3*np.pi/2]:
        c,s = np.cos(theta),np.sin(theta)
        b = rad.Fld(obj,'b',[radius*c,radius*s,z])
        expected = [base[0]*c-base[1]*s,base[0]*s+base[1]*c,base[2]]
        np.testing.assert_allclose(b,expected,rtol=2e-11,atol=1e-13)


@pytest.mark.parametrize('radius', [1e-12, 1e-10, 1e-8, 1e-7])
def test_near_axis_radial_field(radius):
    obj = rad.ObjArcCur([0,0,0],[.010,.013],[0,2*np.pi],.003,40,'man','z',1e7)
    z, dz = .006, 1e-6
    derivative = (reference_axis(.010,.013,.003,1e7,z+dz)
                  - reference_axis(.010,.013,.003,1e7,z-dz))/(2*dz)
    b = np.asarray(rad.Fld(obj,'b',[radius,0,z]))
    assert b[0] == pytest.approx(-radius*derivative/2,rel=2e-6)
    assert abs(b[1]) < 1e-15
