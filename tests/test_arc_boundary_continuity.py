import numpy as np
import pytest
import radia as rad


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
