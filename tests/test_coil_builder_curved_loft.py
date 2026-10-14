"""Curved stream-tube volume oracle, independent of native line segments."""
import numpy as np
import pytest
import radia as rad
from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile


def reference(order):
    x, w = np.polynomial.legendre.leggauss(order)
    a, b, s = np.meshgrid(x/2, x/2, (x+1)/2, indexing='ij')
    theta = np.pi/2*s
    radius = 20 + a*(4+4*s)
    pos = np.stack((-20+radius*np.cos(theta), radius*np.sin(theta), b*(6+4*s)), axis=-1)
    tangent = np.stack((4*a*np.cos(theta)-radius*np.pi/2*np.sin(theta),
                        4*a*np.sin(theta)+radius*np.pi/2*np.cos(theta), 4*b), axis=-1)
    delta = np.array([30, 10, 20]) - pos
    weights = np.einsum('i,j,k->ijk', w/2, w/2, w/2)
    return 1e-5*np.sum(np.cross(tangent, delta)/np.linalg.norm(delta, axis=-1)[..., None]**3
                       * weights[..., None], axis=(0,1,2))


def coil(angle=90, radius=20):
    return CoilBuilder(100).set_cross_section(4,6).add_loft_arc(RectProfile(8,10),radius,angle)


def test_joint_path_and_section_convergence():
    exact = reference(24)
    np.testing.assert_allclose(exact, reference(32), rtol=1e-11, atol=1e-16)
    errors=[]
    for n in [8,16,32,64]:
        obj=rad.ObjCnt(coil().to_radia_loft_filaments(n,n,n_arc=n))
        value=np.array(rad.Fld(obj,'b',[30,10,20]))
        errors.append(np.linalg.norm(value-exact)/np.linalg.norm(exact))
    assert all(y<x/3 for x,y in zip(errors,errors[1:]))
    assert errors[-1]<1e-4


def test_bend_exit_and_next_straight_are_continuous(monkeypatch):
    c=coil().add_straight(10)
    calls=[]
    monkeypatch.setattr(rad,'ObjFlmCur',lambda p,i: calls.append((p,i)) or len(calls))
    c.to_radia_loft_filaments(2,3,n_arc=8)
    assert len(calls)==6
    assert all(len(p)==10 for p,i in calls)
    assert sum(i for p,i in calls)==pytest.approx(100)


def test_path_refinement_with_fixed_section():
    values=[]
    for n in [8,16,32,64]:
        obj=rad.ObjCnt(coil().to_radia_loft_filaments(4,4,n_arc=n))
        values.append(np.array(rad.Fld(obj,'b',[30,10,20])))
    increments=[np.linalg.norm(b-a) for a,b in zip(values,values[1:])]
    assert all(b<a/3 for a,b in zip(increments,increments[1:]))


def test_rigid_transform_covariance():
    from scipy.spatial.transform import Rotation
    rotation=Rotation.from_euler('xyz',[23,41,-17],degrees=True).as_matrix()
    shift=np.array([3,-7,9])
    point=np.array([30,10,20])
    moved=CoilBuilder(100).set_start(shift,rotation.T).set_cross_section(4,6)
    moved.add_loft_arc(RectProfile(8,10),20,90)
    a=rad.ObjCnt(coil().to_radia_loft_filaments(8,8,16))
    b=rad.ObjCnt(moved.to_radia_loft_filaments(8,8,16))
    expected=rotation@np.array(rad.Fld(a,'b',point))
    actual=rad.Fld(b,'b',rotation@point+shift)
    np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-16)


@pytest.mark.parametrize('angle,radius',[(-90,20),(0,20),(361,20),(90,4),(90,float('nan'))])
def test_unsupported_or_degenerate_arc_rejected_before_allocation(monkeypatch,angle,radius):
    def fail(*args):
        pytest.fail('allocated invalid conductor')
    monkeypatch.setattr(rad,'ObjFlmCur',fail)
    with pytest.raises(ValueError):
        coil(angle,radius).to_radia_loft_filaments(2,2)
