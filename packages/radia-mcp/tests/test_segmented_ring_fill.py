r"""Segmented PM ring -> continuous smeared ring geometry correction.

A 2-D axisymmetric or mean-field surrogate revolves a segment into a full ring.  For an actual ring of
discrete blocks, the first-order correction is the azimuthal fill factor: scale remanence/coercivity by
occupied angle / pitch before using the continuous-ring surrogate.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (
    segmented_ring_fill_factor,
    smeared_ring_equivalent_remanence,
)


def test_segmented_ring_fill_factor_limits():
    assert segmented_ring_fill_factor(12, 30.0) == pytest.approx(1.0)
    assert segmented_ring_fill_factor(12, 15.0) == pytest.approx(0.5)
    assert segmented_ring_fill_factor(8, 22.5) == pytest.approx(0.5)

    for n in (6, 10, 18):
        span = 0.37 * 360.0 / n
        assert segmented_ring_fill_factor(n, span) == pytest.approx(0.37)


def test_smeared_ring_equivalent_remanence_scales_source():
    res = smeared_ring_equivalent_remanence(1.2, 12, 15.0)
    assert res["fill_factor"] == pytest.approx(0.5)
    assert res["effective_remanence"] == pytest.approx(0.6)
    assert res["solid_to_segmented_volume_ratio"] == pytest.approx(2.0)

    vec = smeared_ring_equivalent_remanence((1.0, -2.0, 0.5), 10, 18.0)
    assert vec["fill_factor"] == pytest.approx(0.5)
    assert vec["effective_remanence"] == pytest.approx((0.5, -1.0, 0.25))


def test_segmented_ring_fill_preserves_thin_ring_magnet_volume():
    n, span_deg, r_mid, radial_thickness, length = 16, 12.0, 25e-3, 3e-3, 10e-3
    fill = segmented_ring_fill_factor(n, span_deg)
    solid_ring_volume = 2.0 * math.pi * r_mid * radial_thickness * length
    segmented_volume = n * (math.radians(span_deg) * r_mid) * radial_thickness * length
    assert segmented_volume / solid_ring_volume == pytest.approx(fill)


def test_segmented_ring_fill_validation():
    for bad in (
        lambda: segmented_ring_fill_factor(0, 10.0),
        lambda: segmented_ring_fill_factor(6, 0.0),
        lambda: segmented_ring_fill_factor(6, 61.0),
    ):
        with pytest.raises(ValueError):
            bad()
