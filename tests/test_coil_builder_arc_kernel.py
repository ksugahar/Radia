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
