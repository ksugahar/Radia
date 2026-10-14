"""Exercise the public builder through the native finite-section kernel."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
import radia as rad
from radia.coil_builder import CoilBuilder


@pytest.mark.parametrize('angle', [360., 120., -120.])
@pytest.mark.parametrize('angles', [[0., 0., 0.], [31., 47., 19.]])
def test_builder_arc_matches_rigidly_transformed_kernel(angle, angles):
    rad.FldLenRndSw('off')
    try:
        frame = Rotation.from_euler('xyz', angles, degrees=True).as_matrix()
        center = np.array([.004, -.007, .009])
        radius, width, height, current = .02, .008, .012, 120.
        builder = CoilBuilder(current).set_start(center+radius*frame[:,0], frame.T)
        builder.set_cross_section(width, height).add_arc(radius, angle)
        objects = builder.to_radia()
        source = rad.ObjCnt(objects)
        phi = np.deg2rad(angle)
        limits = [0., phi] if phi > 0 else [2*np.pi+phi, 2*np.pi]
        reference = rad.ObjArcCur([0,0,0], [radius-width/2, radius+width/2],
            limits, height, 4, 'man', 'z', np.sign(phi)*current/(width*height))
        local = np.array([[0,0,0], [.002,.003,.005], [.03,.004,.01]])
        expected = np.array([rad.Fld(reference,'b',p) for p in local])@frame.T
        actual = np.array([rad.Fld(source,'b',p) for p in local@frame.T+center])
        assert actual.shape == expected.shape
        np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=1e-13)
    finally:
        rad.FldLenRndSw('on')


@pytest.mark.parametrize('direction', [-1., 1.])
@pytest.mark.parametrize('angles', [[0.,0.,0.], [31.,47.,19.]])
def test_builder_nearly_closed_arc_against_axis_primitive(direction, angles):
    rad.FldLenRndSw('off')
    try:
        radius, width, height, current, z = .02, .008, .012, 120., .03
        frame = Rotation.from_euler('xyz', angles, degrees=True).as_matrix()
        center = np.array([.004,-.007,.009])
        phi = direction*(2*np.pi-1e-7)
        builder = CoilBuilder(current).set_start(center+radius*frame[:,0], frame.T)
        builder.set_cross_section(width,height).add_arc(radius,np.rad2deg(phi))
        source = rad.ObjCnt(builder.to_radia())
        ri, ro = radius-width/2, radius+width/2
        lo, hi = (0.,phi) if phi>0 else (2*np.pi+phi,2*np.pi)
        def primitive(v):
            return v*np.log((ro+np.hypot(ro,v))/(ri+np.hypot(ri,v)))
        def radial(v):
            return np.hypot(ro,v)-np.hypot(ri,v)
        transverse = radial(z-height/2)-radial(z+height/2)
        axial = primitive(z+height/2)-primitive(z-height/2)
        expected = 1e-7*direction*current/(width*height)*np.array([
            transverse*(np.sin(hi)-np.sin(lo)),
            transverse*(np.cos(lo)-np.cos(hi)), axial*(hi-lo)])
        actual = np.asarray(rad.Fld(source,'b',center+z*frame[:,2]))@frame
        np.testing.assert_allclose(actual,expected,rtol=2e-7,atol=1e-15)
    finally:
        rad.FldLenRndSw('on')
