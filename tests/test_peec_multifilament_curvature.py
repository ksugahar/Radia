"""Multi-filament subdivision on a CURVED conductor: measured limits.

`nwinc`/`nhinc` subdivide one straight prism at a time.  The builder
sees a single segment and has no curvature information, so a sub-
filament is a parallel-offset copy that KEEPS THE PARENT LENGTH.  On a
polygonal / curved coil the outward-offset strand therefore traces a
path with the parent's perimeter instead of its own (longer) one, and
its loop inductance comes out too small -- which INVERTS the radial
current ordering.

These tests lock the two halves of that statement so the caveat stays
measurable: the partial-inductance kernel itself is correct (standalone
rings match the analytic loop inductance and grow with radius), while
the offset sub-filaments do not re-length.  Do not "fix" the second by
weakening the first.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.peec_matrices import PEECBuilder  # noqa: E402
from radia.peec_topology import PEECCircuitSolver  # noqa: E402

MU0 = 4.0e-7 * math.pi
N_SEG = 36
GAP = math.radians(8.0)
SIGMA = 5.8e7


def _ring(radius, w, h, *, nwinc=1, nhinc=1):
    b = PEECBuilder()
    ang = np.linspace(GAP / 2, 2 * math.pi - GAP / 2, N_SEG + 1)
    nd = [b.add_node_at(radius * math.cos(t), radius * math.sin(t), 0.0)
          for t in ang]
    for i in range(N_SEG):
        b.add_connected_segment(nd[i], nd[i + 1], w, h, sigma=SIGMA,
                                nwinc=nwinc, nhinc=nhinc)
    b.add_port(nd[0], nd[-1])
    return b.build_topology()


def _analytic_loop_L(radius, w, h):
    # Circular loop of rectangular cross-section, GMD radius approximation
    a = 0.2235 * (w + h)
    return MU0 * radius * (math.log(8.0 * radius / a) - 2.0)


def test_standalone_ring_inductance_matches_analytic_and_grows():
    """The partial-inductance kernel is correct: a ring built at its true
    radius reproduces mu0 R (ln(8R/a) - 2) and grows with R."""
    w, h = 4e-3, 0.8e-3
    radii = (0.02331, 0.025, 0.02651)
    Ls = []
    for r in radii:
        L = float(np.asarray(_ring(r, w, h)["L"]).sum())
        Ls.append(L)
        assert L == pytest.approx(_analytic_loop_L(r, w, h), rel=0.03)
    assert Ls[0] < Ls[1] < Ls[2]


def test_offset_subfilaments_keep_parent_length():
    """Documented limitation: nhinc offsets do NOT re-length the strand,
    so every sub-filament of a curved segment has the parent's length."""
    topo = _ring(0.025, 4e-3, 4e-3, nhinc=5)
    lengths = np.asarray(topo["segment_lengths"])
    assert lengths.shape == (5 * N_SEG,)
    # all strands identical in length despite spanning +-1.6 mm in radius
    assert lengths.std() == pytest.approx(0.0, abs=1e-15)
    centers = np.asarray(topo["segment_centers"])
    dr = np.hypot(centers[:, 0], centers[:, 1]) - 0.025
    assert dr.max() - dr.min() > 3.0e-3


def test_curved_radial_subdivision_inverts_the_ring_inductance_order():
    """Consequence: the offset rings' loop inductance DECREASES outward,
    the opposite of the standalone rings above.  This is why nwinc/nhinc
    must not be used to resolve current across the curvature direction."""
    topo = _ring(0.025, 4e-3, 4e-3, nhinc=5)
    L = np.asarray(topo["L"])
    centers = np.asarray(topo["segment_centers"])
    dr = np.round((np.hypot(centers[:, 0], centers[:, 1]) - 0.025) * 1e3, 2)
    classes = sorted(set(dr.tolist()))
    ring_L = [float(L[np.ix_(dr == r, dr == r)].sum()) for r in classes]
    assert len(classes) == 5
    # monotonically DECREASING outward -- the inverted order
    assert all(a > b for a, b in zip(ring_L, ring_L[1:]))


def test_cross_section_axis_is_not_user_selectable():
    """Probe, don't guess: which physical direction nwinc splits is chosen
    by the builder (perpendicular to the segment), not by the caller."""
    b = PEECBuilder()
    n0 = b.add_node_at(0.0, 0.0, 0.0)
    n1 = b.add_node_at(0.1, 0.0, 0.0)
    b.add_connected_segment(n0, n1, 4e-3, 2e-3, sigma=SIGMA,
                            nwinc=5, nhinc=1)
    b.add_port(n0, n1)
    c = np.asarray(b.build_topology()["segment_centers"])
    # x-directed segment: the 5 width strands are stacked along z, not y
    assert c[:, 1].std() == pytest.approx(0.0, abs=1e-15)
    assert c[:, 2].std() > 0.0

    b = PEECBuilder()
    n0 = b.add_node_at(0.0, 0.0, 0.0)
    n1 = b.add_node_at(0.0, 0.0, 0.1)
    b.add_connected_segment(n0, n1, 4e-3, 2e-3, sigma=SIGMA,
                            nwinc=5, nhinc=1)
    b.add_port(n0, n1)
    c = np.asarray(b.build_topology()["segment_centers"])
    # z-directed segment: z is unavailable, so the strands stack along y
    assert c[:, 2].std() == pytest.approx(0.0, abs=1e-15)
    assert c[:, 1].std() > 0.0


def test_straight_bundle_skin_ordering_is_physical():
    """Where the model IS valid (a straight prism), the high-frequency
    distribution is the inductance-limited one and crowds to the corners
    with a reverse interior strand."""
    b = PEECBuilder()
    n0 = b.add_node_at(0.0, 0.0, 0.0)
    n1 = b.add_node_at(0.1, 0.0, 0.0)
    b.add_connected_segment(n0, n1, 4e-3, 4e-3, sigma=SIGMA,
                            nwinc=3, nhinc=3)
    b.add_port(n0, n1)
    topo = b.build_topology()
    solver = PEECCircuitSolver(topo)
    I = np.real(solver.compute_branch_currents(1e8, [1.0]))

    L = np.asarray(topo["L"])
    limit = np.linalg.solve(L, np.ones(9))
    limit /= limit.sum()
    assert np.allclose(I, limit, atol=1e-4)   # inductance-limited
    assert I.sum() == pytest.approx(1.0)      # Kirchhoff

    c = np.asarray(topo["segment_centers"])
    off = np.abs(c[:, 1]) + np.abs(c[:, 2])
    corners = I[off > 0.9 * off.max()]
    center = I[off < 1e-12]
    assert corners.min() > 0.2                # corners crowd
    assert center.max() < 0.0                 # interior runs backwards
