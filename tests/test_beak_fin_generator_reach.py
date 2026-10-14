"""The beak's reach is a parameter; its tip radius is not.

Separating the two is what lets a sweep say whether an error belongs to the
fin's aspect or to the curvature at its tip.  If ``tip_x_mm`` moved the tip
circle as well, the geometry sweep would be varying both at once and could not
distinguish them, so that independence is checked rather than assumed.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation_test" / "induction_heating"))


def _analyse(tip_x_mm, tmp_path):
    from build123d import export_step

    from make_beak_fin_step import make_beak_fin
    from radia.peec_fin_topology import (
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    step = tmp_path / f"beak_{tip_x_mm:.1f}.step"
    export_step(make_beak_fin(12.0, tip_x_mm), str(step))
    _, cad = build_hybrid_surface_topology_from_straight_prism_step(
        step, n_peri=256, n_stations=3, lane_grading="uniform")
    return cad["section_analysis"].primary.as_dict()


@pytest.mark.parametrize("tip_x_mm", [4.0, 3.2, 2.4])
def test_the_arc_is_the_same_circle_at_every_reach(tip_x_mm):
    """The generator's contract, checked on the construction itself.

    The tip is a three-point arc through ``(centre + 0.05, +-h)`` and
    ``(tip_x, 0)``; those three points are equidistant from ``(centre, 0)``
    exactly, for any reach, so the circle is ``TIP_RADIUS`` by construction and
    not by fit.  That is what makes the reach and the curvature independent.
    """
    from make_beak_fin_step import TIP_HALF_HEIGHT, TIP_RADIUS

    centre = tip_x_mm - TIP_RADIUS
    for point in ((centre + 0.05, TIP_HALF_HEIGHT),
                  (centre + 0.05, -TIP_HALF_HEIGHT), (tip_x_mm, 0.0)):
        radius = math.hypot(point[0] - centre, point[1])
        assert radius == pytest.approx(TIP_RADIUS, rel=1e-12)


@pytest.mark.parametrize("tip_x_mm", [4.0, 3.2, 2.4])
def test_reach_moves_the_tip_and_leaves_the_root(tip_x_mm, tmp_path):
    """A shorter beak ends sooner, from the same root, rounded the same way.

    The fitted radius is checked loosely on purpose.  A shorter beak has a
    steeper flank running into the same arc, so fewer outline samples are
    unambiguously on it and the least-squares circle drifts -- 1.5% at a 2.4 mm
    reach.  That is the detector's resolution, not the geometry moving, which
    the construction test above pins exactly.
    """
    pytest.importorskip("build123d")
    from make_beak_fin_step import TIP_RADIUS

    fin = _analyse(tip_x_mm, tmp_path)
    assert fin["extremity_m"][0] == pytest.approx(tip_x_mm * 1e-3, abs=5e-6)
    assert fin["tip_radius_m"] == pytest.approx(TIP_RADIUS * 1e-3, rel=2e-2)
    # The root stays put, so the reach is the only thing that moved.
    assert fin["root_mid_m"][0] == pytest.approx(1.6e-3, abs=2e-5)


def test_the_tip_cut_follows_the_fin_rather_than_a_constant(tmp_path):
    """A fixed 3.5 mm tip cut falls past the end of a short beak.

    The sweep derives its cuts from the fin the section analysis found, which
    is what keeps the tip region non-empty when the reach changes; a constant
    would silently divide by zero, or worse, not.
    """
    pytest.importorskip("build123d")

    short = _analyse(2.4, tmp_path)
    tip_cut = short["extremity_m"][0] - short["tip_span_m"]
    assert tip_cut < 3.5e-3, "the constant cut would have been outside the fin"
    assert short["root_mid_m"][0] < tip_cut < short["extremity_m"][0]


def test_reach_outside_the_root_is_refused():
    """A beak that does not clear its own root is not a beak."""
    pytest.importorskip("build123d")
    from make_beak_fin_step import make_beak_fin

    for bad in (1.5, 1.9, 4.5):
        with pytest.raises(ValueError, match="1.6 mm root"):
            make_beak_fin(12.0, bad)


def test_the_tip_is_a_half_circle_of_the_declared_radius():
    """The generator's own arithmetic, checked where it is easy to get wrong."""
    from make_beak_fin_step import TIP_HALF_HEIGHT, TIP_RADIUS

    assert TIP_HALF_HEIGHT == pytest.approx(
        math.sqrt(TIP_RADIUS ** 2 - 0.05 ** 2), rel=1e-12)
    assert 0.0 < TIP_HALF_HEIGHT < TIP_RADIUS
