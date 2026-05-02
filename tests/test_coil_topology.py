"""Unit tests for ``radia.coil_topology`` -- the unified OPEN/CLOSED
classifier for coil solids.

Covers:
  - CLOSED full torus (``ih_closed_torus_coil.step``): no caps,
    sweep_deg = 360.
  - OPEN circular cross-section (``ih_peec_inductance_coil.step``):
    2 PLANE caps + 4 TORUS lateral, sweep ~ 355 deg.
  - OPEN rect cross-section, multi-loft united (``rect_torus_lofted_united.step``):
    28 PLANE faces, 2 smallest are caps, sweep ~ 355 deg.
  - OPEN multi-turn (``3turnCoil_work_coil_m.step``): 2 PLANE + 617
    BSPLINE.

Generate_spine round-trip:
  - For OPEN, first station ~ cap A, last station ~ cap B.
  - For CLOSED, first station at theta=0, last at theta = 2pi - dtheta.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SAMPLES = os.path.join(ROOT, "src", "radia", "panels", "samples")
GOLDEN = os.path.join(ROOT, "tests", "panels", "golden")


def _have_build123d():
    try:
        from build123d import import_step  # noqa: F401
        return True
    except Exception:
        return False


def _import(step_name: str, root: str):
    from build123d import import_step
    p = os.path.join(root, step_name)
    if not os.path.isfile(p):
        pytest.skip(f"fixture missing: {p}")
    return import_step(p)


# ---------------------------------------------------------------------------
# Cap detection
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_detect_cap_closed_full_torus_returns_none():
    from radia.coil_topology import detect_cap_faces
    sld = _import("ih_closed_torus_coil.step", SAMPLES)
    caps = detect_cap_faces(sld)
    assert caps is None, (
        f"closed torus should have no caps, got {caps}")


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_detect_cap_open_circular_cross_section():
    from radia.coil_topology import detect_cap_faces
    sld = _import("ih_peec_inductance_coil.step", SAMPLES)
    caps = detect_cap_faces(sld)
    assert caps is not None
    cap_a, cap_b = caps
    A_expected = math.pi * (0.003) ** 2   # = 28.27 mm^2
    assert abs(cap_a.area - A_expected) / A_expected < 0.05
    assert abs(cap_b.area - A_expected) / A_expected < 0.05


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_detect_cap_open_rect_lofted_united():
    from radia.coil_topology import detect_cap_faces
    sld = _import("rect_torus_lofted_united.step", GOLDEN)
    caps = detect_cap_faces(sld)
    assert caps is not None
    cap_a, cap_b = caps
    A_expected = 0.008 * 0.006   # 8 x 6 mm = 4.8e-5 m^2
    assert abs(cap_a.area - A_expected) / A_expected < 0.05
    assert abs(cap_b.area - A_expected) / A_expected < 0.05


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_detect_cap_open_multi_turn_3turncoil():
    from radia.coil_topology import detect_cap_faces
    sld = _import("3turnCoil_work_coil_m.step", SAMPLES)
    caps = detect_cap_faces(sld)
    assert caps is not None


# ---------------------------------------------------------------------------
# extract_coil_topology
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_topology_closed_full_torus():
    from radia.coil_topology import extract_coil_topology
    sld = _import("ih_closed_torus_coil.step", SAMPLES)
    topo = extract_coil_topology(sld)
    assert topo.is_open is False
    assert topo.sweep_deg == 360.0
    assert topo.cap_a is None and topo.cap_b is None


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_topology_open_circular_355_sweep():
    from radia.coil_topology import extract_coil_topology
    sld = _import("ih_peec_inductance_coil.step", SAMPLES)
    topo = extract_coil_topology(sld)
    assert topo.is_open is True
    # gap = 5 deg, sweep = 360 - 5 = 355
    assert abs(topo.sweep_deg - 355.0) < 1.0, (
        f"sweep_deg = {topo.sweep_deg:.2f}, expected ~355")


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_topology_open_rect_355_sweep():
    from radia.coil_topology import extract_coil_topology
    sld = _import("rect_torus_lofted_united.step", GOLDEN)
    topo = extract_coil_topology(sld)
    assert topo.is_open is True
    assert abs(topo.sweep_deg - 355.0) < 1.0


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_topology_open_multi_turn():
    from radia.coil_topology import extract_coil_topology
    sld = _import("3turnCoil_work_coil_m.step", SAMPLES)
    topo = extract_coil_topology(sld)
    assert topo.is_open is True
    # gap was reported as ~11 deg in the survey, so sweep ~ 349.
    # Multi-turn coils have a more complex axis; the *cap* sweep
    # angle is just the in-plane gap.
    assert 340.0 < topo.sweep_deg < 360.0


# ---------------------------------------------------------------------------
# generate_spine
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_spine_closed_full_360():
    from radia.coil_topology import extract_coil_topology, generate_spine
    sld = _import("ih_closed_torus_coil.step", SAMPLES)
    topo = extract_coil_topology(sld)
    spine = generate_spine(topo, 20)
    assert spine.shape == (20, 3)
    # First station angle ~ 0
    th0 = math.atan2(spine[0, 1], spine[0, 0])
    assert abs(th0) < 1e-9
    # endpoint=False -> last station at theta = 2pi*(19/20) = 342 deg
    th_last = math.atan2(spine[-1, 1], spine[-1, 0]) % (2 * math.pi)
    assert abs(math.degrees(th_last) - 342.0) < 0.1


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_spine_open_endpoints_match_caps():
    """For OPEN, the first and last spine stations should match the
    cap-face centroid angles (within numerical tolerance)."""
    from radia.coil_topology import extract_coil_topology, generate_spine
    sld = _import("rect_torus_lofted_united.step", GOLDEN)
    topo = extract_coil_topology(sld)
    spine = generate_spine(topo, 20)
    th0 = math.atan2(spine[0, 1], spine[0, 0])
    th_last = math.atan2(spine[-1, 1], spine[-1, 0])
    # First station = cap A angle (within 1 mdeg)
    assert abs(math.degrees(th0 - topo.theta_a)) < 0.001
    # Last station = cap B angle (mod 2pi -- in case of CW long arc)
    delta = math.degrees(th_last - topo.theta_b) % 360.0
    if delta > 180.0:
        delta = 360.0 - delta
    assert delta < 0.001, (
        f"last spine station {math.degrees(th_last):.4f} deg vs "
        f"cap_b {math.degrees(topo.theta_b):.4f} deg")


@pytest.mark.skipif(not _have_build123d(), reason="build123d unavailable")
def test_spine_open_sweeps_long_arc_not_short():
    """The total angular span of the spine must be the LONG arc
    (~355 deg for a 5 deg gap), NOT the short gap arc (~5 deg)."""
    from radia.coil_topology import extract_coil_topology, generate_spine
    sld = _import("rect_torus_lofted_united.step", GOLDEN)
    topo = extract_coil_topology(sld)
    spine = generate_spine(topo, 20)
    # Sum of unsigned angular increments along spine
    total_deg = 0.0
    for i in range(len(spine) - 1):
        th_i = math.atan2(spine[i, 1], spine[i, 0])
        th_n = math.atan2(spine[i + 1, 1], spine[i + 1, 0])
        d = (th_n - th_i + math.pi) % (2 * math.pi) - math.pi
        total_deg += abs(math.degrees(d))
    assert 354.0 < total_deg < 356.0, (
        f"spine spans {total_deg:.2f} deg, expected ~355 (long arc)")
