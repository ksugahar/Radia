import numpy as np
import pytest
import radia as rad
from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile


def ring():
    return CoilBuilder(100).set_cross_section(4,6).add_loft_arc(RectProfile(4,6),20,360)


def test_closed_ring_matches_finite_section_analytic_field():
    # Uniform azimuthal current density integrated over r=[18,22], z=[-3,3].
    expected=2*np.pi*1e-7*(100/24)*6*(np.arcsinh(22/3)-np.arcsinh(18/3))
    errors=[]
    for n in (32,64,128):
        value=np.array(rad.Fld(rad.ObjCnt(ring().to_radia_loft_filaments(16,16,n_arc=n,require_closed=True)),
                               'b',[-20,0,0]))
        assert np.linalg.norm(value[:2])<1e-15
        errors.append(abs(value[2]/expected-1))
    assert errors[-1]<3e-4
    assert errors[-1]<errors[0]/8


def test_closed_paths_snap_without_adding_return_wire(monkeypatch):
    calls=[]
    monkeypatch.setattr(rad,'ObjFlmCur',lambda p,i: calls.append((p,i)) or len(calls))
    ring().to_radia_loft_filaments(2,3,n_arc=16,require_closed=True)
    assert len(calls)==6
    assert sum(i for p,i in calls)==pytest.approx(100)
    assert all(len(p)==17 and p[0]==p[-1] for p,i in calls)


@pytest.mark.parametrize('kind',['open','mismatched_full_turn'])
def test_invalid_closure_rejected_before_allocation(monkeypatch,kind):
    coil=CoilBuilder(100).set_cross_section(4,6)
    if kind=='open':
        coil.add_loft_straight(RectProfile(4,6),20)
    else:
        coil.add_loft_arc(RectProfile(8,10),20,360)
    monkeypatch.setattr(rad,'ObjFlmCur',lambda *a: pytest.fail('premature allocation'))
    with pytest.raises(ValueError):
        coil.to_radia_loft_filaments(1,1,require_closed=True)
