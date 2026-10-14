import numpy as np
import pytest
import radia as rad


@pytest.mark.parametrize('coordinate', [np.nan,np.inf,-np.inf])
def test_invalid_arc_query_reports_error_and_preserves_model(coordinate):
    rad.FldLenRndSw('off')
    try:
        obj=rad.ObjArcCur([0,0,0],[.01,.013],[0,2*np.pi],.003,4,'man','z',1e7)
        expected=np.asarray(rad.Fld(obj,'b',[0,0,0]))
        with pytest.raises(RuntimeError,match='Field evaluation failed'):
            rad.Fld(obj,'b',[coordinate,0,0])
        np.testing.assert_allclose(rad.Fld(obj,'b',[0,0,0]),expected,rtol=0,atol=0)
    finally:
        rad.FldLenRndSw('on')
