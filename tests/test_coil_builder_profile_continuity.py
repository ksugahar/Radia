"""Cross-section state must follow each segment's exit, not its old entry."""
import numpy as np
import pytest

from radia.coil_builder import CoilBuilder
from radia.coil_profile import CircleProfile, RectProfile


@pytest.mark.parametrize('start', [RectProfile(.004, .006), CircleProfile(.003)])
@pytest.mark.parametrize('next_arc', [False, True])
def test_straight_loft_inherits_entry_and_propagates_exit(start, next_arc):
    end = RectProfile(.008, .010)
    coil = CoilBuilder(100).set_profile(start)
    coil.add_loft_straight(end, .02)
    if next_arc:
        coil.add_arc(.03, 90)
    else:
        coil.add_straight(.02)
    loft, following = coil.segments
    assert loft.profile_start is start
    assert following.profile is end
    np.testing.assert_allclose(loft.end_pos, following.start_pos, atol=1e-15)
    alpha = np.array([.1, .5, .9])
    beta = np.array([.2, .6, .8])
    np.testing.assert_allclose(loft.profile_at(1).sample_at(alpha, beta),
                               following.profile.sample_at(alpha, beta), atol=1e-15)


def test_explicit_loft_entry_does_not_override_exit_inheritance():
    entry = RectProfile(.002, .003)
    end = RectProfile(.008, .010)
    coil = CoilBuilder(100).set_profile(CircleProfile(.003))
    coil.add_loft_straight(end, .02, profile_start=entry).add_straight(.02)
    assert coil.segments[0].profile_start is entry
    assert coil.segments[1].profile is end
