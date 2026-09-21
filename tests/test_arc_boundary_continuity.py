import numpy as np
import pytest
import radia as rad
from scipy.special import ellipkm1, ellipe


def loop_reference(radius, z, order):
    nodes, weights=np.polynomial.legendre.leggauss(order)
    t,w=(nodes+1)/2,weights/2
    def intervals(lo, hi):
        edges=[lo]+([0.] if lo<0<hi else [])+[hi]
        for a,b in zip(edges,edges[1:]):
            start,end=(a,b) if abs(a)<abs(b) else (b,a)
            yield start+(end-start)*t**4,w*abs(end-start)*4*t**3
    result=np.zeros(3)
    for dr_nodes,wr in intervals(.035-radius,.070-radius):
        for dz_nodes,wz in intervals(z-.0525,z+.0525):
            dr,dz=dr_nodes[:,None],dz_nodes[None,:]
            beta=dr*dr+dz*dz
            alpha=(2*radius+dr)**2+dz*dz
            complement=beta/alpha
            K,E=ellipkm1(complement),ellipe(1-complement)
            pref=2e-7*1e6/np.sqrt(alpha)
            br=pref*dz/radius*(-K+(2*radius**2+2*radius*dr+beta)/beta*E)
            bz=pref*(K+(dr*(2*radius+dr)-dz*dz)/beta*E)
            quadrature=wr[:,None]*wz[None,:]
            result += [np.sum(br*quadrature),0,np.sum(bz*quadrature)]
    return result


@pytest.mark.parametrize('point', [[.035,0,0], [.07,0,0], [.035,0,.0525],
                                 [.07,0,.0525], [.05,0,.0525]])
def test_boundary_against_independent_elliptic_loop_formula(point):
    coarse=loop_reference(point[0],point[2],128)
    fine=loop_reference(point[0],point[2],256)
    np.testing.assert_allclose(coarse,fine,rtol=1e-8,atol=1e-12)
    rad.FldLenRndSw('off')
    try:
        obj=rad.ObjArcCur([0,0,0],[.035,.07],[0,2*np.pi],.105,40,'man','z',1e6)
        np.testing.assert_allclose(rad.Fld(obj,'b',point),fine,rtol=1e-8,atol=1e-12)
    finally:
        rad.FldLenRndSw('on')


@pytest.mark.parametrize('angle', [1.7, 2*np.pi])
@pytest.mark.parametrize('point', [[.035,0,0], [.07,0,0], [.035,0,.0525],
                                 [.07,0,.0525], [.05,0,.0525]])
def test_finite_boundary_field_matches_two_sided_limit(angle, point):
    rad.FldLenRndSw('off')
    try:
        obj=rad.ObjArcCur([0,0,0],[.035,.07],[0,angle],.105,40,'man','z',1e6)
        point=np.asarray(point,dtype=float)
        actual=np.asarray(rad.Fld(obj,'b',point))
        assert actual.shape==(3,) and np.isfinite(actual).all()
        delta=np.array([1.,1.,1.])*1e-8
        near=np.array([rad.Fld(obj,'b',point+s*delta) for s in [-1,1]])
        assert np.isfinite(near).all()
        np.testing.assert_allclose(actual,near.mean(axis=0),rtol=3e-5,atol=1e-9)
    finally:
        rad.FldLenRndSw('on')
