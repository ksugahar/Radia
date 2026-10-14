"""Structural regression tests for lead-aware PEEC chain reconstruction.

Verifies that ``_filaments_from_circle_edges_per_station`` (Path 2b in
``filaments_from_step``) produces a chain that:

1. Starts at one lead's OUTER tip (= port +) and ends at the other
   lead's OUTER tip (= port -).
2. Contains no segment longer than 5x the median segment length
   (no fake U-bends through air).
3. Walks both lead bodies (full axial extent on each lead's axis).

These checks are independent of the numerical L_coil value (which is
locked by ``test_inductance_peec_vacuum_3turn_loft`` in
test_inductance_golden.py) -- they pin the TOPOLOGY of the chain.

Pre-2026-05-13 the chain on 3turnCoil_work_coil.step had:
  - last seg = 62.71mm (Lead 1 body as a single straight bar)
  - and a 25mm fake U-bend in air to a dangling lead 2 outer tip
  - Lead 2 body entirely missing from the chain
Fixed by ``_augment_chain_with_lead_bars`` in coil_from_cad.py.
"""
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_REPO, "src", "radia"))
sys.path.insert(0, os.path.join(_REPO, "src", "radia", "panels"))


SAMPLE_STEP = os.path.join(_REPO, "src", "radia", "panels", "samples",
                           "3turnCoil_work_coil.step")


@pytest.fixture(scope="module")
def topology():
    if not os.path.isfile(SAMPLE_STEP):
        pytest.skip(f"Sample STEP not available: {SAMPLE_STEP}")
    import radia  # noqa: F401  -- ensure DLL setup runs first
    from radia.coil_from_cad import filaments_from_step
    return filaments_from_step(SAMPLE_STEP, n_peri=8, sigma=5.8e7,
                                cad_units_per_meter=1.0)


def _waypoints_of_filament(filament_polyline):
    pts = [filament_polyline[0][0]]
    for (_p0, p1) in filament_polyline:
        pts.append(p1)
    return np.array(pts)


def test_chain_has_no_fake_long_segments(topology):
    """Every filament's max segment must be <= 5x its median segment.

    A fake U-bend through air or a single-segment lead body would
    show up as a max-over-median ratio of 30-50+ on the kubota
    3turnCoil geometry.  The augmented chain subdivides lead bodies
    to ~lead_length / median_helix_seg stations and drops fake jumps.
    """
    fp = topology["filament_paths"]
    for k, polyline in enumerate(fp):
        lens = np.array([np.linalg.norm(np.array(p1) - np.array(p0))
                          for (p0, p1) in polyline])
        med = float(np.median(lens))
        mx = float(lens.max())
        assert mx <= 5.0 * med, (
            f"Filament {k}: max seg {mx*1000:.2f}mm > 5x median "
            f"{med*1000:.2f}mm -- lead-aware chain regression "
            f"(fake U-bend or single-seg lead body present).")


def test_chain_starts_and_ends_at_lead_outer_tips(topology):
    """First filament's first / last waypoints must be at the
    expected lead OUTER tip X positions (around X = 129.43 mm for
    3turnCoil_work_coil.step).  The Y differs across filaments
    because of the n_peri perimeter offset (r*cos(theta), r*sin(theta))
    around the centerline, but X is on the lead axis -- so the X
    coordinate of the first/last filament waypoint equals the lead
    outer tip X (modulo perimeter offset which has zero X-component
    since the lead axis is along +/- X).
    """
    fp = topology["filament_paths"]
    waypoints_0 = _waypoints_of_filament(fp[0])
    first_x_mm = waypoints_0[0, 0] * 1000.0
    last_x_mm = waypoints_0[-1, 0] * 1000.0
    assert abs(first_x_mm - 129.43) < 0.5, (
        f"Filament 0 first waypoint X = {first_x_mm:.2f}mm, "
        f"expected ~129.43mm (lead outer tip).")
    assert abs(last_x_mm - 129.43) < 0.5, (
        f"Filament 0 last waypoint X = {last_x_mm:.2f}mm, "
        f"expected ~129.43mm (other lead outer tip).")


def test_chain_has_both_lead_bodies(topology):
    """At least N=5 waypoints must lie within 1.5*lead_radius of
    each lead's axis line at distinct along-axis parameters, for
    BOTH leads (no lead body should be missing).
    """
    fp = topology["filament_paths"]
    waypoints = _waypoints_of_filament(fp[0])

    # Two lead axes in 3turnCoil geometry: y=+12.5 (lead 2) and y=-12.5 (lead 1)
    # extending in +/- X from helix (x~68) to lead tip (x~129).
    lead_radius_with_offset = 3.15 + 3.15 + 1.0  # axis radius + perimeter offset + slack (mm)
    leads = [
        {"y_axis": 12.5, "x_min": 68.0, "x_max": 129.43, "name": "+y lead"},
        {"y_axis": -12.5, "x_min": 68.0, "x_max": 129.43, "name": "-y lead"},
    ]
    for L in leads:
        # Distance from each waypoint to this lead's axis line in mm.
        # Axis is along X at y=L["y_axis"], z=10.  Perpendicular
        # distance to axis = sqrt((y - y_axis)**2 + (z - 10)**2).
        wp_mm = waypoints * 1000.0
        perp = np.sqrt((wp_mm[:, 1] - L["y_axis"]) ** 2
                        + (wp_mm[:, 2] - 10.0) ** 2)
        on_axis = perp < lead_radius_with_offset
        in_x_range = (wp_mm[:, 0] >= L["x_min"] - 1.0) & \
                     (wp_mm[:, 0] <= L["x_max"] + 1.0)
        # waypoints on the lead body
        on_lead = on_axis & in_x_range
        assert int(on_lead.sum()) >= 5, (
            f"Lead '{L['name']}' has only {int(on_lead.sum())} waypoints "
            f"on-axis in X range [{L['x_min']:.1f}..{L['x_max']:.1f}]mm "
            f"-- lead body missing or under-sampled.")


def test_chain_endpoint_separation_matches_lead_offset(topology):
    """The chain's first and last waypoints (= port + and port -) are
    at the OUTER tips of the two leads.  They should be ~25mm apart
    (the Y offset between the two lead axes: |+12.5 - (-12.5)| = 25mm).
    """
    fp = topology["filament_paths"]
    waypoints = _waypoints_of_filament(fp[0])
    dy_mm = abs(waypoints[0, 1] - waypoints[-1, 1]) * 1000.0
    # Filament 0 is offset by perimeter placement so dy is NOT exactly
    # 25.0 -- but it should be in the 20-30 mm range.
    assert 20.0 < dy_mm < 30.0, (
        f"Filament 0 port-to-port Y separation {dy_mm:.2f}mm -- "
        f"expected ~25mm (lead axes at y=+/-12.5).")
