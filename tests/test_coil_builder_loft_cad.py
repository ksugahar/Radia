import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile

pytest.importorskip('netgen.occ')


@pytest.mark.parametrize('angles', [[0,0,0],[23,41,-17]])
def test_straight_loft_cad_matches_current_tube(angles, tmp_path):
    rotation=Rotation.from_euler('xyz',angles,degrees=True).as_matrix()
    shift=np.array([3,-7,9])
    coil=CoilBuilder(100).set_start(shift,rotation.T).set_cross_section(4,6)
    coil.add_loft_straight(RectProfile(8,10),20)
    shape=coil.segments[0].to_occ_shape()
    corners=np.array([[x,y,z] for y,w,h in [(0,4,6),(20,8,10)]
                      for x in [-w/2,w/2] for z in [-h/2,h/2]])
    corners=corners@rotation.T+shift
    expected=np.array([corners.min(axis=0),corners.max(axis=0)])
    box=shape.bounding_box
    actual=np.array([[box[i][j] for j in range(3)] for i in range(2)])
    np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-5)
    assert shape.mass == pytest.approx(20*(24+20+16/3),rel=1e-10)
    from netgen.occ import OCCGeometry
    path=tmp_path/'straight_loft.step'
    shape.WriteStep(str(path))
    recovered=OCCGeometry(str(path)).shape
    assert len(recovered.solids)==1
    assert recovered.mass == pytest.approx(shape.mass,rel=1e-9)
    box=recovered.bounding_box
    np.testing.assert_allclose([[box[i][j] for j in range(3)] for i in range(2)],
                               expected,rtol=0,atol=1e-5)
