"""Native conversion must not omit or replace unsupported conductor geometry."""
import pytest
import numpy as np
import radia as rad
from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile, CircleProfile


@pytest.mark.parametrize('kind', ['loft_straight', 'loft_arc', 'circle_straight', 'circle_arc'])
def test_native_export_rejects_unsupported_geometry_before_allocation(monkeypatch, kind):
    coil=CoilBuilder(100).set_cross_section(.004,.006).add_straight(.02)
    if kind=='loft_straight':
        coil.add_loft_straight(RectProfile(.006,.008), .02)
    elif kind=='loft_arc':
        coil.add_loft_arc(RectProfile(.006,.008), .02, 90)
    elif kind=='circle_straight':
        coil.add_straight(.02, profile=CircleProfile(.002))
    else:
        coil.add_arc(.02, 90, profile=CircleProfile(.002))
    created=[]
    original = rad.ObjRecCur
    def record_allocation(*args, **kwargs):
        created.append(args)
        return original(*args, **kwargs)
    monkeypatch.setattr(rad, 'ObjRecCur', record_allocation)
    with pytest.raises(NotImplementedError, match='segment 1'):
        coil.to_radia()
    assert created==[]


@pytest.mark.parametrize('arc', [False, True])
def test_rectangular_profile_dimensions_must_match_native_geometry(monkeypatch, arc):
    coil = CoilBuilder(100).set_cross_section(.004, .006)
    profile = RectProfile(.008, .010)
    if arc:
        coil.add_arc(.02, 90, profile=profile)
    else:
        coil.add_straight(.02, profile=profile)
    def unexpected_allocation(*args, **kwargs):
        pytest.fail('native allocation before geometry validation')
    monkeypatch.setattr(rad, 'ObjRecCur', unexpected_allocation)
    monkeypatch.setattr(rad, 'ObjArcCur', unexpected_allocation)
    with pytest.raises(ValueError, match='segment 0: profile and native dimensions'):
        coil.to_radia()


def test_matching_explicit_rectangle_preserves_native_field():
    plain = CoilBuilder(100).set_cross_section(.004, .006)
    explicit = CoilBuilder(100).set_cross_section(.004, .006)
    plain.add_straight(.02).add_arc(.02, 90)
    explicit.add_straight(.02, profile=RectProfile(.004, .006)).add_arc(.02, 90)
    a, b = plain.to_radia(), explicit.to_radia()
    assert len(a) == len(b) == 2
    for obj_a, obj_b in zip(a, b):
        assert rad.Fld(obj_a, 'b', [.1, .1, .1]) == pytest.approx(
            rad.Fld(obj_b, 'b', [.1, .1, .1]), rel=1e-13, abs=1e-15)


@pytest.mark.parametrize('fault', ['length', 'radius', 'angle', 'current',
                                   'position', 'frame', 'subdivision'])
def test_invalid_tail_does_not_allocate_valid_prefix(monkeypatch, fault):
    coil=CoilBuilder(100).set_cross_section(.004,.006).add_straight(.02)
    coil.add_arc(.02,90)
    tail=coil.segments[-1]
    kwargs={}
    if fault=='length':
        coil.add_straight(.02)
        coil.segments[-1].length=-1
    elif fault=='radius':
        tail.radius=.001
    elif fault=='angle':
        tail.arc_angle=361
    elif fault=='current':
        tail.current=float('nan')
    elif fault=='position':
        tail.start_pos=np.array([float('inf'),0,0])
    elif fault=='frame':
        tail.orientation=np.diag([2.,1.,1.])
    else:
        kwargs['arc_max_segment_length']=0
    def unexpected(*args, **kwargs):
        pytest.fail('allocated native object before complete preflight')
    monkeypatch.setattr(rad,'ObjRecCur',unexpected)
    monkeypatch.setattr(rad,'ObjArcCur',unexpected)
    with pytest.raises(ValueError):
        coil.to_radia(**kwargs)
