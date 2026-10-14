"""Circular prescribed-current lofts checked against independent volume integrals."""
import numpy as np
import pytest
import radia as rad
from scipy.spatial.transform import Rotation
from radia.coil_builder import CoilBuilder, LoftStraightSegment, LoftArcSegment
from radia.coil_profile import CircleProfile


def make_coil(curved, rotation=None, shift=None):
    rotation=np.eye(3) if rotation is None else rotation
    shift=np.zeros(3) if shift is None else shift
    args=(100,shift,rotation.T,CircleProfile(2),CircleProfile(4))
    segment=LoftArcSegment(*args,20,90) if curved else LoftStraightSegment(*args,20)
    coil=CoilBuilder(100)
    coil.segments.append(segment)
    return coil


def reference(point, curved, order=24):
    x,w=np.polynomial.legendre.leggauss(order)
    # Integrate in physical normalized radius q, not the production sqrt(alpha) grid.
    q,b,s=np.meshgrid((x+1)/2,np.pi*(x+1),(x+1)/2,indexing='ij')
    weights=np.einsum('i,j,k->ijk',w/2,w/2,w/2)*2*q
    u=(2+2*s)*q*np.cos(b)
    v=(2+2*s)*q*np.sin(b)
    du=2*q*np.cos(b)
    dv=2*q*np.sin(b)
    if curved:
        theta=np.pi*s/2
        c,t=np.cos(theta),np.sin(theta)
        source=np.stack((-20+(20+u)*c,(20+u)*t,v),axis=-1)
        tangent=np.stack((du*c-(20+u)*t*np.pi/2,
                          du*t+(20+u)*c*np.pi/2,dv),axis=-1)
    else:
        source=np.stack((u,20*s,v),axis=-1)
        tangent=np.stack((du,np.full_like(u,20),dv),axis=-1)
    delta=np.asarray(point)-source
    integrand=np.cross(tangent,delta)/np.linalg.norm(delta,axis=-1)[...,None]**3
    return 1e-7*100*np.sum(weights[...,None]*integrand,axis=(0,1,2))


@pytest.mark.parametrize('curved',[False,True])
def test_circle_loft_independent_field_and_pose(curved):
    point=np.array([30,10,20])
    exact=reference(point,curved)
    np.testing.assert_allclose(exact,reference(point,curved,32),rtol=1e-10,atol=1e-16)
    errors=[]
    for n in (8,16,32,64):
        value=np.array(rad.Fld(rad.ObjCnt(make_coil(curved).to_radia_loft_filaments(n,n,n_arc=n)), 'b',point.tolist()))
        errors.append(np.linalg.norm(value-exact)/np.linalg.norm(exact))
    assert errors[-1]<1e-4
    assert all(b<a/3 for a,b in zip(errors,errors[1:]))
    rotation=Rotation.from_euler('xyz',[23,41,-17],degrees=True).as_matrix()
    shift=np.array([3,-7,9])
    moved=make_coil(curved,rotation,shift)
    result=rad.Fld(rad.ObjCnt(moved.to_radia_loft_filaments(64,64,n_arc=64)), 'b',(rotation@point+shift).tolist())
    np.testing.assert_allclose(result,rotation@value,rtol=1e-10,atol=1e-16)


def test_circle_join_checks_boundary_before_allocation(monkeypatch):
    coil=make_coil(False)
    coil.segments.append(LoftStraightSegment(100,np.array([0,20,0]),np.eye(3),CircleProfile(3),CircleProfile(5),10))
    monkeypatch.setattr(rad,'ObjFlmCur',lambda *args: pytest.fail('allocated before checking join'))
    with pytest.raises(ValueError,match='disconnected cross-section'):
        coil.to_radia_loft_filaments(1,1)
