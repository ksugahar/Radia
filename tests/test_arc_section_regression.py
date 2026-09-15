"""Arc volume-current tests against independent tensor Gauss integration."""
import numpy as np
import pytest
import radia as rad


@pytest.fixture(autouse=True)
def exact_coordinates():
    rad.FldLenRndSw('off')
    yield
    rad.FldLenRndSw('on')


def volume_reference(point, ri, ro, height, angles, density, order=40):
    x, w = np.polynomial.legendre.leggauss(order)
    rr = ((ri+ro)/2+(ro-ri)/2*x)[:, None, None]
    zz = (height/2*x)[None, :, None]
    phi = ((angles[0]+angles[1])/2+(angles[1]-angles[0])/2*x)[None, None, :]
    dx, dy, dz = point[0]-rr*np.cos(phi), point[1]-rr*np.sin(phi), point[2]-zz
    inv = (dx*dx+dy*dy+dz*dz)**-1.5
    weights = w[:, None, None]*w[None, :, None]*w[None, None, :]
    factor = 1e-7*density*(ro-ri)*height*(angles[1]-angles[0])/8
    return factor*np.array([np.sum(weights*rr*np.cos(phi)*dz*inv),
                            np.sum(weights*rr*np.sin(phi)*dz*inv),
                            np.sum(weights*rr*(rr-point[0]*np.cos(phi)-point[1]*np.sin(phi))*inv)])


@pytest.mark.parametrize('angles', [(0., 1.7), (0.2, 4.8), (0., 1e-5)])
@pytest.mark.parametrize('point', [(0., 0., .08), (.012, -.009, .065), (.09, .01, .08)])
def test_thick_arc_volume_reference(angles, point):
    ri, ro, height, density = .035, .070, .105, 1e6
    expected = volume_reference(point, ri, ro, height, angles, density)
    refined = volume_reference(point, ri, ro, height, angles, density, 64)
    np.testing.assert_allclose(expected, refined, rtol=2e-8, atol=1e-13)
    obj = rad.ObjArcCur([0,0,0], [ri,ro], angles, height, 40, 'man', 'z', density)
    actual = np.array(rad.Fld(obj, 'b', point))
    assert np.linalg.norm(actual-refined) < 2e-7*np.linalg.norm(refined)


@pytest.mark.parametrize('scale', [.001, 1., 1000.])
def test_arc_scale_rotation_and_partition(scale):
    ri, ro, height = np.array([.035, .070, .105])*scale
    point = np.array([.012, -.009, .065])*scale
    def field(angles, p, density=1e6):
        obj=rad.ObjArcCur([0,0,0], [ri,ro], angles, height, 4, 'man','z',density)
        return np.asarray(rad.Fld(obj,'b',p))
    b=field([.2,4.8],point)
    parts=field([.2,1.2],point)+field([1.2,4.8],point)
    np.testing.assert_allclose(b,parts,rtol=1e-8,atol=1e-13*scale)
    theta=.6
    rotation=np.array([[np.cos(theta),-np.sin(theta),0],
                       [np.sin(theta),np.cos(theta),0],[0,0,1]])
    np.testing.assert_allclose(field([.8,5.4],rotation@point),rotation@b,
                               rtol=1e-8,atol=1e-13*scale)
    expected=volume_reference(point/scale,.035,.070,.105,[.2,4.8],1e6)*scale
    np.testing.assert_allclose(b,expected,rtol=2e-7,atol=1e-13*scale)
    np.testing.assert_allclose(field([.2,4.8],point,-1e6),-b,rtol=1e-12)


@pytest.mark.parametrize('point', [(0.,0.,0.), (.02,0.,0.), (.02,0.,.0525)])
def test_arc_bore_and_end_plane(point):
    angles=[0.,2.8]
    expected=volume_reference(point,.035,.070,.105,angles,1e6,80)
    refined=volume_reference(point,.035,.070,.105,angles,1e6,112)
    np.testing.assert_allclose(expected,refined,rtol=1e-8,atol=1e-12)
    obj=rad.ObjArcCur([0,0,0],[.035,.070],angles,.105,4,'man','z',1e6)
    np.testing.assert_allclose(rad.Fld(obj,'b',point),refined,rtol=2e-7,atol=1e-12)


@pytest.mark.parametrize('point', [(0.,0.,0.),(0.,0.,.08),(.02,0.,.08)])
def test_thick_full_circle_partition(point):
    def field(lo,hi):
        obj=rad.ObjArcCur([0,0,0],[.035,.070],[lo,hi],.105,4,'man','z',1e6)
        return np.asarray(rad.Fld(obj,'b',point))
    whole=field(0,2*np.pi)
    split=field(0,np.pi)+field(np.pi,2*np.pi)
    expected=volume_reference(point,.035,.070,.105,[0,2*np.pi],1e6,80)
    np.testing.assert_allclose(whole,expected,rtol=2e-7,atol=1e-12)
    np.testing.assert_allclose(whole,split,rtol=2e-7,atol=1e-12)


@pytest.mark.parametrize('radius,z', [(.05,0.),(.045,.02),(.02,.01)])
def test_full_circle_local_ampere_and_divergence(radius,z):
    obj=rad.ObjArcCur([0,0,0],[.035,.070],[0,2*np.pi],.105,4,'man','z',1e6)
    def field(r,z):
        return np.asarray(rad.Fld(obj,'b',[r,0,z]))
    center=field(radius,z)
    step=2e-6
    dr=(field(radius+step,z)-field(radius-step,z))/(2*step)
    dz=(field(radius,z+step)-field(radius,z-step))/(2*step)
    expected=4*np.pi*1e-7*1e6 if radius>.035 else 0.
    assert dz[0]-dr[2] == pytest.approx(expected,abs=2e-5)
    assert abs(dr[0]+center[0]/radius+dz[2])<2e-5
