import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile, CircleProfile

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


@pytest.mark.parametrize('angles', [[0,0,0],[23,41,-17]])
def test_arc_loft_cad_volume_convergence_and_step(angles, tmp_path):
    from netgen.occ import OCCGeometry
    rotation=Rotation.from_euler('xyz',angles,degrees=True).as_matrix()
    coil=CoilBuilder(100).set_start([3,-7,9],rotation.T).set_cross_section(4,6)
    coil.add_loft_arc(RectProfile(8,10),20,90)
    segment=coil.segments[0]
    # Integral of R * angle * w(s) * h(s), centered radial section.
    exact=20*np.pi/2*(24+20+16/3)
    errors=[]
    for count in [8,16,32]:
        segment.n_sub=count
        shape=segment.to_occ_shape()
        assert len(shape.solids)==1
        errors.append(abs(shape.mass/exact-1))
    assert errors[-1]<1e-5
    assert errors[-1]<errors[0]
    path=tmp_path/'arc_loft.step'
    shape.WriteStep(str(path))
    recovered=OCCGeometry(str(path)).shape
    assert len(recovered.solids)==1
    assert recovered.mass==pytest.approx(shape.mass,rel=1e-8)


@pytest.mark.parametrize('radius,angle,count', [(4,90,20),(20,-90,20),
                                               (20,360,20),(20,90,3)])
def test_arc_loft_cad_rejects_unsupported_geometry(radius,angle,count):
    coil=CoilBuilder(100).set_cross_section(4,6)
    coil.add_loft_arc(RectProfile(8,10),radius,angle,n_sub=count)
    with pytest.raises(ValueError):
        coil.segments[0].to_occ_shape()


@pytest.mark.parametrize('curved', [False, True])
@pytest.mark.parametrize('angles', [[0,0,0],[23,41,-17]])
def test_circular_loft_volume_and_step(curved, angles, tmp_path):
    from netgen.occ import OCCGeometry
    from radia.coil_builder import LoftStraightSegment, LoftArcSegment
    rotation=Rotation.from_euler('xyz',angles,degrees=True).as_matrix()
    args=(100,np.array([3,-7,9]),rotation.T,CircleProfile(2),CircleProfile(4))
    if curved:
        segment=LoftArcSegment(*args,20,90,n_sub=32)
        length=20*np.pi/2
    else:
        segment=LoftStraightSegment(*args,20)
        length=20
    exact=length*np.pi*(4+8+16)/3
    shape=segment.to_occ_shape()
    assert len(shape.solids)==1
    assert shape.mass==pytest.approx(exact,rel=1e-5)
    if curved:
        segment.n_sub=4
        coarse=segment.to_occ_shape()
        assert abs(shape.mass-exact)<abs(coarse.mass-exact)
        segment.n_sub=64
        finer=segment.to_occ_shape()
        assert finer.mass==pytest.approx(shape.mass,rel=1e-8)
    path=tmp_path/'circle_loft.step'
    shape.WriteStep(str(path))
    recovered=OCCGeometry(str(path)).shape
    assert len(recovered.solids)==1
    assert recovered.mass==pytest.approx(shape.mass,rel=1e-8)


def test_cross_type_cad_is_not_silently_replaced():
    from radia.coil_builder import LoftStraightSegment, LoftArcSegment
    args=(100,np.zeros(3),np.eye(3),RectProfile(4,6),CircleProfile(4))
    for segment in (LoftStraightSegment(*args,20),LoftArcSegment(*args,20,90)):
        with pytest.raises(NotImplementedError):
            segment.to_occ_shape()
