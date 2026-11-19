"""Opt-in near-pair section-average mutual inductance of the PEEC kernel."""
from __future__ import annotations

import math

import numpy as np
import pytest

MU0 = 4e-7 * math.pi


def _parallel_filaments(length, d):
    """Grover: two equal, aligned, parallel filaments a distance d apart."""
    r = length / d
    return MU0 * length / (2 * math.pi) * (
        math.log(r + math.sqrt(1 + r * r)) - math.sqrt(1 + 1 / (r * r)) + 1 / r)


def _section_average_reference(length, w, h, gap, n=8):
    """Independent volume average of the mutual of two side-by-side bars
    (offset along the width axis by w + gap), n x n Gauss per section."""
    x, wt = np.polynomial.legendre.leggauss(n)
    total = 0.0
    for xi, wi in zip(0.5 * w * x, wt):
        for yi, vi in zip(0.5 * h * x, wt):
            for xj, wj in zip(w + gap + 0.5 * w * x, wt):
                for yj, vj in zip(0.5 * h * x, wt):
                    total += wi * vi * wj * vj * _parallel_filaments(
                        length, math.hypot(xj - xi, yj - yi))
    return total / 16.0


def _two_bars(builder, length, w, h, gap, *, width_dir=None):
    a0 = builder.add_node_at(0.0, 0.0, 0.0)
    a1 = builder.add_node_at(0.0, 0.0, length)
    builder.add_connected_segment(a0, a1, w, h)
    if width_dir is not None:
        builder.set_last_width_dir(*width_dir)
    b0 = builder.add_node_at(w + gap, 0.0, 0.0)
    b1 = builder.add_node_at(w + gap, 0.0, length)
    builder.add_connected_segment(b0, b1, w, h)
    if width_dir is not None:
        builder.set_last_width_dir(*width_dir)
    return np.asarray(builder.build_topology()["L"])


def test_default_and_zero_factor_keep_the_centreline_mutual():
    """Guards: the near-pair rule is opt-in; existing builders (default and
    factor 0) give bit-identical L, equal to the centreline Grover value.
    Assumes: two parallel 10 mm bars, 1 mm wide, 0.2 mm gap."""
    from radia.peec_matrices import PEECBuilder

    length, w, h, gap = 10e-3, 1e-3, 0.2e-3, 0.2e-3
    ref = _two_bars(PEECBuilder(), length, w, h, gap)
    zero = _two_bars(PEECBuilder(near_section_average=0.0), length, w, h, gap)
    assert np.array_equal(ref, zero)
    centre = _parallel_filaments(length, w + gap)
    assert abs(ref[0, 1] - centre) <= 1e-6 * centre, (ref[0, 1], centre)


def test_near_pairs_get_the_section_average_and_far_pairs_do_not():
    """Guards: with the factor on, a near pair's mutual equals an
    independent 8x8-Gauss volume average within the 3x3 rule's accuracy,
    and a pair beyond the factor keeps the centreline value.
    Assumes: width axis along x (the offset direction)."""
    from radia.peec_matrices import PEECBuilder

    length, w, h, gap = 10e-3, 1e-3, 0.2e-3, 0.2e-3
    near = _two_bars(PEECBuilder(near_section_average=2.0), length, w, h, gap,
                     width_dir=(1.0, 0.0, 0.0))
    ref = _section_average_reference(length, w, h, gap)
    centre = _parallel_filaments(length, w + gap)
    assert abs(near[0, 1] - ref) <= 2e-3 * ref, (near[0, 1], ref, centre)
    assert abs(centre - ref) > 5 * abs(near[0, 1] - ref), (near[0, 1], ref, centre)

    far_gap = 0.5  # 0.5 m apart: beyond 2 x (sum of bounding radii)
    far = _two_bars(PEECBuilder(near_section_average=2.0), length, w, h, far_gap,
                    width_dir=(1.0, 0.0, 0.0))
    far_centre = _two_bars(PEECBuilder(), length, w, h, far_gap)
    assert far[0, 1] == far_centre[0, 1]


def test_near_rule_rejects_bad_input():
    """Guards: fail-loud contract of the opt-in rule.
    Assumes: none."""
    from radia.peec_matrices import PEECBuilder

    with pytest.raises(Exception, match="factor"):
        PEECBuilder(near_section_average=float("nan"))
    with pytest.raises(Exception, match="exclusive"):
        PEECBuilder(volume_order=8, near_section_average=2.0)
    b = PEECBuilder(near_section_average=2.0)
    n0, n1 = b.add_node_at(0.0, 0.0, 0.0), b.add_node_at(0.0, 0.0, 1e-3)
    b.add_connected_segment(n0, n1, 1e-4, 1e-4)
    with pytest.raises(Exception, match="parallel"):
        b.set_last_width_dir(0.0, 0.0, 1.0)
