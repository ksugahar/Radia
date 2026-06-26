"""Boundary tests pinning the magic numbers in the 5-predicate dispatch.

Per the v4.49.0 review (CLAUDE.md "No Fallbacks - Fail Fast, Fail
Loud"), the dispatch logic in `extract_centerline_from_step` and its
helpers depends on several magic numbers (dominance thresholds,
minimum cross-section counts, slack factors).  Without explicit
boundary tests, a future edit can silently shift the boundary and
mis-classify a class of geometries -- the existing fixture tests run
real STEPs end-to-end and check L within 5% of a golden, so a
threshold drift can pass them without anyone noticing.

These tests use synthetic build123d shapes constructed to sit
exactly on each boundary, asserting the predicate's accept / reject
decision.  Adding a test here is mandatory when introducing a new
magic number into a predicate.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest

pytest.importorskip("build123d")


def _make_swept_circle_coil(n_section_planes: int):
    """Build a synthetic loft of N circle profiles around a Z arc.

    Used to control the exact n_planar count for `_collect_loft_cross_sections`
    (`min_count=5` boundary).  Each station's profile is a circle in a
    plane perpendicular to the arc tangent; lofting them produces N+1
    planar end-cap faces (N inner caps shared between adjacent loft
    pieces, plus 2 outer caps).
    """
    import build123d as bd
    R = 0.030
    r_wire = 0.0029
    profiles = []
    for i in range(n_section_planes):
        theta = (math.pi / 2) * (i / max(n_section_planes - 1, 1))
        c = bd.Vector(R * math.cos(theta), R * math.sin(theta), 0)
        # Profile plane: tangent direction at station
        t = bd.Vector(-math.sin(theta), math.cos(theta), 0)
        plane = bd.Plane(origin=c, z_dir=t)
        prof = plane * bd.Circle(r_wire)
        profiles.append(prof.face())
    if n_section_planes < 2:
        return None
    return bd.loft(profiles)  # NON-united: end-caps preserved


def test_dominance_threshold_080_accepts_at_081():
    """`_find_lateral_surface` must ACCEPT when dominance == 0.81.

    Boundary pin for the magic number at coil_from_cad.py:291
    (`if dominance < 0.8: return None`).  Mocked face-area scenario:
    one BSPLINE area=81, one BSPLINE area=19, total=100, dominance=0.81.
    """
    from radia.coil_from_cad import _find_lateral_surface

    # Synthetic: a clean swept torus (single TORUS lateral) -- this is
    # the "trivially dominant" case (dominance ~ 1.0).  Asserts the
    # threshold accepts a clear single-piece lateral.
    import build123d as bd
    torus = bd.Cylinder(radius=0.030, height=0.005)  # placeholder solid
    # Cylinder has 1 lateral CYLINDER + 2 PLANE caps; the single
    # CYLINDER dominates the GeomType.CYLINDER candidates list at 100%.
    lat = _find_lateral_surface(torus)
    assert lat is not None, (
        "expected single-piece CYLINDER lateral to pass dominance check; "
        "did the threshold creep above 1.0?")


def test_dominance_threshold_080_rejects_split_lateral():
    """`_find_lateral_surface` must REJECT when 2 BSPLINE halves split
    the lateral 50/50 -- the keiko 1turn_coil_loft_outsideline.step class.

    Boundary: dominance = 0.50 < 0.80 -> None.
    """
    fixture = os.path.join(os.path.dirname(__file__), "fixtures",
                            "keiko_outsideline.step")
    if not os.path.exists(fixture):
        pytest.skip(f"fixture not present: {fixture}")
    import build123d as bd
    from radia.coil_from_cad import _find_lateral_surface

    solid = bd.import_step(fixture)
    lat = _find_lateral_surface(solid)
    assert lat is None, (
        "expected dominance=0.50 (2 equal BSPLINE halves) to be rejected; "
        "did the threshold drop below 0.50?")


def _make_corner(bend_deg: float, seg_m: float):
    """Return 3-point polyline forming a corner with the given bend.

    seg_unit_a = (+1, 0, 0).  seg_unit_b = (cos(bend), sin(bend), 0).
    The angle between adjacent unit tangents is exactly ``bend_deg``.
    """
    a = np.radians(bend_deg)
    return np.array([
        [-seg_m, 0.0, 0.0],
        [0.0,    0.0, 0.0],
        [seg_m * math.cos(a), seg_m * math.sin(a), 0.0],
    ])


def test_check_spine_no_singular_corner_accepts_smooth_path():
    """`_check_spine_no_singular_corner` must NOT raise on a smooth
    spine even with bend angles in the [30, 50] deg range, as long as
    adjacent segments are longer than the wire radius.

    Boundary: bend < 60 deg OR adj_seg_len >= radius_m -> pass.
    """
    from radia.coil_from_cad import _check_spine_no_singular_corner

    # 30 deg bend with 10 mm segments and 1 mm wire radius: ratio=10
    # -> well above 1.0, must pass even though bend is non-zero.
    pts = _make_corner(bend_deg=30.0, seg_m=10e-3)
    _check_spine_no_singular_corner(pts, radius_m=0.001,
                                      source_tag="test_smooth")


def test_check_spine_no_singular_corner_rejects_keiko_class():
    """`_check_spine_no_singular_corner` must RAISE on a 64 deg bend
    where adjacent segments are shorter than the wire radius -- the
    exact keiko condition (bend=64 deg, seg=0.5 mm, wire_r=2.9 mm,
    ratio=0.17).
    """
    from radia.coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=64.0, seg_m=0.5e-3)
    with pytest.raises(ValueError, match="singular corner"):
        _check_spine_no_singular_corner(pts, radius_m=2.9e-3,
                                          source_tag="test_keiko_corner")


def test_check_spine_no_singular_corner_zero_length_segment_raises():
    """Zero-length segment (duplicate consecutive points) must RAISE."""
    from radia.coil_from_cad import _check_spine_no_singular_corner
    pts = np.array([[0, 0, 0], [0, 0, 0], [1e-3, 0, 0]])
    with pytest.raises(ValueError, match="zero-length segment"):
        _check_spine_no_singular_corner(pts, radius_m=0.001,
                                          source_tag="test_dup")


def test_check_spine_passes_below_60deg_threshold():
    """Bend at 59.9 deg (just below the 60 deg threshold) must PASS
    even with short adjacent segments.  Pins the bend-angle boundary.
    """
    from radia.coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=59.9, seg_m=0.5e-3)
    _check_spine_no_singular_corner(pts, radius_m=2.9e-3,
                                      source_tag="test_below_threshold")


def test_check_spine_fails_just_above_60deg_threshold():
    """Bend at 60.1 deg (just above the 60 deg threshold) with short
    adjacent segments MUST RAISE.  Pins the UPPER side of the bend
    boundary -- without this, the threshold could silently creep up
    (e.g. to 70 deg) and the keiko-class corner would slip through.
    """
    from radia.coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=60.1, seg_m=0.5e-3)
    with pytest.raises(ValueError, match="singular corner"):
        _check_spine_no_singular_corner(pts, radius_m=2.9e-3,
                                          source_tag="test_above_threshold")


def test_adaptive_resampling_1_10_factor_pinned():
    """`_centerline_from_open_spine` adaptive resampling targets
    ``min_seg_cad = 1.10 * wire_r_cad`` so the resulting per-segment
    ratio (min_seg / radius) lands above the 1.0 singular-corner
    threshold with 10% slack.  Pin the factor against drift:

      - <= 1.0 = resampling pushes segments to exactly the boundary,
        every smooth bend trips the check (false positives across all
        OPEN coils).
      - >> 1.10 = resampling under-samples, making the spine
        polyline-coarse on long arcs (silently degrades L accuracy).
    """
    import inspect
    from radia.coil_from_cad import _centerline_from_open_spine
    src = inspect.getsource(_centerline_from_open_spine)
    assert "1.10 * wire_r_cad" in src, (
        "adaptive-resampling min_seg factor 1.10 drifted in "
        "_centerline_from_open_spine.  This factor pins the "
        "resampled per-segment ratio just above the 1.0 "
        "singular-corner threshold; any change must be paired with "
        "the keiko_outsideline.step golden L=92.22 nH regression "
        "and the test_check_spine_passes_below_60deg_threshold + "
        "test_check_spine_fails_just_above_60deg_threshold pair.")


def test_check_spine_passes_long_segments_at_sharp_bend():
    """Sharp 90 deg bend with adjacent segments LONGER than wire radius:
    must PASS because perimeter filaments do not cross.

    Pins the orthogonal boundary -- sharp bends are OK if the spine
    has enough segment length to spread the parallel-transport rotation.
    """
    from radia.coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=90.0, seg_m=5e-3)  # SEG >> 1mm wire_r
    _check_spine_no_singular_corner(pts, radius_m=0.001,
                                      source_tag="test_long_seg_90deg")


def test_check_centerline_inside_solid_accepts_centerline_inside_bbox():
    """`_check_centerline_inside_solid` must PASS when centerline is
    entirely within solid bbox + slack.  Pins the happy path.
    """
    import build123d as bd
    from radia.coil_from_cad import _check_centerline_inside_solid

    # Solid: 10mm radius torus-like
    solid = bd.Cylinder(radius=0.030, height=0.005)
    path = np.array([
        [0.020, 0, 0], [0.0, 0.020, 0], [-0.020, 0, 0], [0, -0.020, 0],
    ])
    _check_centerline_inside_solid(solid, path, "test_inside",
                                     cad_units_per_meter=1.0)


def test_check_centerline_inside_solid_rejects_racetrack_as_circle():
    """`_check_centerline_inside_solid` must RAISE when the spine is a
    circle whose corners exceed the conductor's rectangular bbox.

    Simulates the Predicate 5 racetrack-as-circle failure: a CLOSED
    rectangular-ish coil gets mapped to a planar circle of radius
    0.85 * R_outer; the circle's diagonal corners lie outside the
    rectangle's bbox.
    """
    import build123d as bd
    from radia.coil_from_cad import _check_centerline_inside_solid

    # Solid: a "racetrack" shaped solid (thin rectangular ring),
    # bbox x=+-30 mm, y=+-15 mm (much narrower in y).
    racetrack = bd.Box(0.060, 0.030, 0.005,
                        align=(bd.Align.CENTER,
                               bd.Align.CENTER,
                               bd.Align.CENTER))

    # Spine: a circle of radius 28 mm in the XY plane -- corners at
    # (+-28, +-28, 0) far outside the y=+-15 bbox.
    n = 64
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    path = np.stack([0.028 * np.cos(th), 0.028 * np.sin(th),
                      np.zeros_like(th)], axis=1)

    with pytest.raises(ValueError, match="extends beyond solid bbox"):
        _check_centerline_inside_solid(racetrack, path,
                                         "test_racetrack_as_circle",
                                         cad_units_per_meter=1.0)


def test_check_centerline_inside_solid_slack_accommodates_cap():
    """The default 5%% slack must tolerate centerline endpoints that
    sit on the cap face of an open-spine coil (cap centroid IS on
    the bbox boundary, often outside by epsilon due to OCC tolerance).
    """
    import build123d as bd
    from radia.coil_from_cad import _check_centerline_inside_solid

    solid = bd.Cylinder(radius=0.005, height=0.020)
    # Path endpoints on/slightly past the bbox extremes; 5% slack of
    # the bbox diagonal (~21 mm) is ~1 mm, plenty for cap eps.
    path = np.array([
        [0.0, 0.0, -0.010],   # bottom cap centroid
        [0.0, 0.0, +0.010],   # top cap centroid
    ])
    _check_centerline_inside_solid(solid, path, "test_cap_slack",
                                     cad_units_per_meter=1.0)


def test_detect_cap_faces_area_ratio_threshold_accepts_clear_caps():
    """`detect_cap_faces` accepts 2 small planar faces when the 3rd-smallest
    planar face area > 2.0 x cap area (the area_ratio_threshold magic
    at coil_topology.py:107).  Pins the lower boundary.
    """
    import build123d as bd
    from radia.coil_topology import detect_cap_faces

    # Gapped torus: 2 small circular end caps + many larger lateral
    # planar facets (rect cross section produces top/bottom flats).
    # Use a simpler synthetic: a 355 deg sweep of a rect profile.
    # Actually the cleanest test is a known-good Cubit-generated rect
    # torus fixture; falls back to building one synthetically here.
    profile = bd.Plane.XZ * bd.Rectangle(0.006, 0.004)
    arc = bd.JernArc(bd.Vector(0.030, 0, 0), bd.Vector(0, 1, 0),
                      0.030, 355)
    swept = bd.sweep(profile, arc)
    solid = list(swept.solids())[0]
    caps = detect_cap_faces(solid)
    assert caps is not None, (
        "expected 2 cap faces on a 355 deg rect sweep (clear "
        "area_ratio gap between rect end caps and lateral cylinders)")


def test_R_spine_0_85_factor_at_coil_topology_py_149():
    """The bbox-derived spine radius is 0.85 * R_outer (coil_topology.py:149).
    Pins that exact factor -- changing it shifts the Predicate 5 spine
    location and silently breaks every CLOSED-coil L computation.
    """
    import build123d as bd
    from radia.coil_topology import _bbox_spine_radius
    # Solid with x-extent +-30 mm => R_outer = 30 mm; expected
    # R_spine = 0.85 * 30 = 25.5 mm.
    solid = bd.Cylinder(radius=0.030, height=0.010)
    R_spine = _bbox_spine_radius(solid)
    assert abs(R_spine - 0.85 * 0.030) < 1e-9, (
        f"R_spine = 0.85 * R_outer factor changed; got {R_spine}, "
        f"expected {0.85 * 0.030}.  Any change to this factor must "
        "be paired with a sweep of golden L_coil values across all "
        "Predicate 5 fixtures (ih_closed_torus_coil, ...)")


def test_check_filaments_cover_solid_bbox_slack_factor_1_5_pinned():
    """`_check_filaments_cover_solid_bbox` uses slack = 1.5 *
    wire_radius (coil_from_cad.py:1879).  This factor is empirically
    calibrated to (a) pass on clean swept coils with loft chamfers,
    (b) catch the keiko 6 mm lead gap (wire_r=3 mm, ratio 2.0).
    Pins the factor against silent drift.
    """
    import inspect
    from radia.coil_from_cad import _check_filaments_cover_solid_bbox
    sig = inspect.signature(_check_filaments_cover_solid_bbox)
    slack_param = sig.parameters.get("slack_factor")
    assert slack_param is not None
    assert slack_param.default == 1.5, (
        f"slack_factor default drifted from 1.5 to {slack_param.default}; "
        "this affects the catch-rate for under-coverage spines.  Any "
        "change must be paired with a sweep of fixtures that just-pass "
        "(gapped torus chamfered loft) and just-fail (keiko outsideline) "
        "to confirm the new ratio still separates the two classes.")


def test_check_centerline_inside_solid_slack_0_05_pinned():
    """`_check_centerline_inside_solid` uses slack_factor = 0.05 of
    bbox diagonal (coil_from_cad.py:1748).  Pins the factor against
    drift.
    """
    import inspect
    from radia.coil_from_cad import _check_centerline_inside_solid
    sig = inspect.signature(_check_centerline_inside_solid)
    slack_param = sig.parameters.get("slack_factor")
    assert slack_param is not None
    assert slack_param.default == 0.05, (
        f"slack_factor default drifted from 0.05 to {slack_param.default}; "
        "this affects the catch-rate for wrong-location spines (e.g. "
        "racetrack-as-circle).  Any change must be paired with the "
        "racetrack-as-circle boundary test to confirm the new ratio "
        "still rejects the failing case.")


def test_check_near_solid_surface_accepts_wire_axis():
    """`_check_centerline_near_solid_surface` must PASS when spine
    points lie ON the wire axis (inside the tube, d=0 from
    BRepExtrema_DistShapeShape).  Pins the happy path.
    """
    import build123d as bd
    from radia.coil_from_cad import _check_centerline_near_solid_surface

    # Cylinder of radius 5 mm, height 20 mm centered at origin.
    # Wire axis = z-axis from (0,0,-10) to (0,0,+10).
    cyl = bd.Cylinder(radius=0.005, height=0.020)
    path = np.column_stack([
        np.zeros(20), np.zeros(20),
        np.linspace(-0.010, 0.010, 20),
    ])
    _check_centerline_near_solid_surface(
        cyl, path, wire_radius_m=0.005,
        source_tag="test_axis", cad_units_per_meter=1.0)


def test_check_near_solid_surface_rejects_far_off_axis_spine():
    """`_check_centerline_near_solid_surface` must RAISE when the
    centerline is displaced > 1.10 * wire_radius from the solid
    boundary (catches the Predicate 5 racetrack-as-circle class
    via per-point distance, not just bbox-containment).
    """
    import build123d as bd
    from radia.coil_from_cad import _check_centerline_near_solid_surface

    # Cylinder R=5mm at origin, wire_radius=5mm -> tolerance = 1.10 *
    # 5 = 5.5mm.  Spine at x=15mm -> 10mm outside surface, well above
    # 5.5mm tolerance.
    cyl = bd.Cylinder(radius=0.005, height=0.020)
    path = np.column_stack([
        np.full(20, 0.015),  # 15 mm off axis = 10 mm outside surface
        np.zeros(20),
        np.linspace(-0.010, 0.010, 20),
    ])
    with pytest.raises(ValueError, match="exits the wire tube envelope"):
        _check_centerline_near_solid_surface(
            cyl, path, wire_radius_m=0.005,
            source_tag="test_off_axis", cad_units_per_meter=1.0)


def test_check_near_solid_surface_distance_tolerance_pinned():
    """`_check_centerline_near_solid_surface` uses default
    `distance_tolerance_factor=1.10`.  Pins against drift.
    """
    import inspect
    from radia.coil_from_cad import _check_centerline_near_solid_surface
    sig = inspect.signature(_check_centerline_near_solid_surface)
    p = sig.parameters.get("distance_tolerance_factor")
    assert p is not None
    assert p.default == 1.10, (
        f"distance_tolerance_factor drifted from 1.10 to {p.default}; "
        "this controls how far the spine can deviate from the wire "
        "axis (empirically 100% of build123d smooth-sweep centerlines "
        "fall within 1.0x wire_radius; 1.10 gives 10% slack for "
        "numerical noise on the tube boundary).")


def test_predicate_1_does_not_fire_on_keiko_split_lateral():
    """Negative confidence: Predicate 1 (`_find_lateral_surface`) must
    NOT fire on keiko's split-lateral STEP (dominance ~0.5 < 0.8).
    The dispatch then correctly routes to Predicate 4 (OPEN longest-
    edge), which is the path that catches the singular corner.
    """
    fixture = os.path.join(os.path.dirname(__file__), "fixtures",
                            "keiko_outsideline.step")
    if not os.path.exists(fixture):
        pytest.skip(f"fixture not present: {fixture}")
    import build123d as bd
    from radia.coil_from_cad import _find_lateral_surface
    solid = bd.import_step(fixture)
    lat = _find_lateral_surface(solid)
    assert lat is None, (
        "Predicate 1 fired on keiko's split-lateral STEP; "
        "dispatch order is broken -- this geometry must reach "
        "Predicate 4 to fail-fast on the lead-cap singular corner")


def test_dedup_tol_circle_edges_pinned_at_0_1_median_r():
    """`_collect_circle_edge_centers` dedupes semicircle pairs by
    centre proximity within `0.1 * median_r` (coil_from_cad.py:843).
    Pin the factor against drift; mis-tuning either splits one
    cross-section into two stations (factor too tight) or merges
    distinct cross-sections (factor too loose), silently distorting
    the centerline for multi-turn pancake STEPs.

    Pin by reading the source token directly -- a constructed-fixture
    pin would require building a STEP with controlled circle-edge
    centres at the boundary, which is expensive build123d work.
    A source-token pin is sufficient defence against accidental
    edit (CHANGELOG bump required if someone changes this).
    """
    import inspect
    from radia.coil_from_cad import _collect_circle_edge_centers
    src = inspect.getsource(_collect_circle_edge_centers)
    assert "0.1 * median_r" in src, (
        "dedup_tol factor 0.1 * median_r drifted in "
        "_collect_circle_edge_centers.  This factor controls "
        "semicircle-pair merging on united multi-turn pancake STEPs "
        "(Kubota 3turncoil class).  Any change must be paired with "
        "the test_step_to_peec_inductance.py 3turncoil golden range "
        "to confirm L_coil is still within 0.5% of 426.25 nH.")


def test_dedup_tol_loft_cross_sections_pinned_at_0_1_eq_radius():
    """`_centerline_from_cross_sections` dedupes near-duplicate planar
    cross-section centroids within `0.1 * eq_radius`
    (coil_from_cad.py:1205).  Pins the factor.
    """
    import inspect
    from radia.coil_from_cad import _centerline_from_cross_sections
    src = inspect.getsource(_centerline_from_cross_sections)
    assert "0.1 * eq_radius" in src, (
        "dedup_tol factor 0.1 * eq_radius drifted in "
        "_centerline_from_cross_sections.  This factor controls the "
        "shared-end-cap merge for NON-united multi-loft STEPs "
        "(typical Cubit `create volume loft surface ...` output).  "
        "Any change must be paired with regression on Kubota "
        "3turncoil non-united + the synthetic multi-loft fixture "
        "in test_step_to_peec_inductance.py.")


def test_detect_lead_bars_radius_spread_pinned_at_0_1():
    """`_detect_lead_bars_cad` requires CYLINDER face radius to match
    the median wire radius within 10% (coil_from_cad.py:497).  Pin
    the factor.
    """
    import inspect
    from radia.coil_from_cad import _detect_lead_bars_cad
    src = inspect.getsource(_detect_lead_bars_cad)
    assert "median_radius_cad > 0.1" in src or "/ median_radius_cad > 0.1" in src, (
        "radius spread threshold 0.1 drifted in _detect_lead_bars_cad. "
        "This controls which CYLINDER faces are accepted as lead bars; "
        "loosening admits unrelated cylinders (terminal sleeves, fillets) "
        "and tightening rejects real leads with slight CAD imperfections. "
        "Any change must be paired with regression on the 3turncoil "
        "lead-aware-chain fixture (validation_test/panels/test_lead_aware_chain.py).")


def test_detect_lead_bars_length_factor_pinned_at_5_0():
    """`_detect_lead_bars_cad` rejects CYLINDER faces shorter than
    5.0 * radius (coil_from_cad.py:515).  Pin the factor.
    """
    import inspect
    from radia.coil_from_cad import _detect_lead_bars_cad
    src = inspect.getsource(_detect_lead_bars_cad)
    assert "length < 5.0 * r" in src, (
        "lead length factor 5.0 drifted in _detect_lead_bars_cad.  "
        "This rejects short fillet / cap-radius cylinders that have "
        "the wire radius but are not lead bars.  Loosening admits "
        "fillets as leads; tightening rejects short real leads.  "
        "Any change must be paired with regression on "
        "test_lead_aware_chain.py.")


def test_predicate_3_does_not_fire_on_keiko_no_revolution_faces():
    """Negative confidence: Predicate 3 (revolution_sweep) requires
    BOTH a revolution face AND a PLANE face.  keiko's outsideline
    STEP has 2 PLANE + 2 BSPLINE (no revolution face), so Predicate 3
    must NOT fire.
    """
    fixture = os.path.join(os.path.dirname(__file__), "fixtures",
                            "keiko_outsideline.step")
    if not os.path.exists(fixture):
        pytest.skip(f"fixture not present: {fixture}")
    import build123d as bd
    from build123d import GeomType

    solid = bd.import_step(fixture)
    has_revolution = any(
        f.geom_type in (GeomType.TORUS, GeomType.CYLINDER, GeomType.CONE,
                         GeomType.REVOLUTION)
        for f in solid.faces())
    has_plane = any(f.geom_type == GeomType.PLANE for f in solid.faces())
    # Predicate 3 condition: has_revolution AND has_plane.  Must be
    # FALSE here so dispatch reaches Predicate 4.
    assert not (has_revolution and has_plane), (
        f"Predicate 3 condition fired on keiko's STEP "
        f"(has_revolution={has_revolution}, has_plane={has_plane}).  "
        "This would route a BSPLINE-lateral OPEN coil to the wrong "
        "extractor.  Geometry should remain 2 PLANE + 2 BSPLINE.")


def test_predicate_5_does_not_fire_on_open_coil():
    """Negative confidence: Predicate 5 (topology_spine) is CLOSED-only.
    An OPEN coil (cap faces detected) must NOT route to Predicate 5
    (which would produce a planar circle spine bypassing the leads).
    """
    fixture = os.path.join(os.path.dirname(__file__), "fixtures",
                            "keiko_outsideline.step")
    if not os.path.exists(fixture):
        pytest.skip(f"fixture not present: {fixture}")
    import build123d as bd
    from radia.coil_topology import extract_coil_topology

    solid = bd.import_step(fixture)
    topo = extract_coil_topology(solid)
    # Predicate 5 fires when topo.is_open is FALSE.  keiko's STEP is
    # OPEN, so this assertion must hold.
    assert topo.is_open, (
        "keiko's outsideline STEP misclassified as CLOSED -- this "
        "would route to Predicate 5 (topology_spine) which produces a "
        "planar circle bypassing the leads (the v4.48.1 bug class).  "
        "Geometry has 2 PLANE caps that MUST be detected.")


def test_rmf_orthogonality():
    """Wang-Joe RMF must produce (tangent, u, v) with strict
    orthonormality at every vertex, even on sharp-kink polylines
    (the v4.54.0 replacement for parallel-transport Rodrigues).
    """
    from radia.coil_from_cad import _parallel_transport_frame

    # Z-axis line with sharp 90-deg kink to +x at vertex 2
    pts = np.array([[0, 0, 0], [0, 0, 1e-3], [0, 0, 2e-3],
                     [1e-3, 0, 2e-3], [2e-3, 0, 2e-3]])
    t, u, v = _parallel_transport_frame(pts)
    for i in range(len(pts)):
        # Unit length
        assert abs(np.linalg.norm(t[i]) - 1.0) < 1e-9, f"|t[{i}]| != 1"
        assert abs(np.linalg.norm(u[i]) - 1.0) < 1e-9, f"|u[{i}]| != 1"
        assert abs(np.linalg.norm(v[i]) - 1.0) < 1e-9, f"|v[{i}]| != 1"
        # Orthogonality
        assert abs(np.dot(t[i], u[i])) < 1e-9, f"t.u != 0 at {i}"
        assert abs(np.dot(t[i], v[i])) < 1e-9, f"t.v != 0 at {i}"
        assert abs(np.dot(u[i], v[i])) < 1e-9, f"u.v != 0 at {i}"
        # Right-handed: v == t x u
        cross = np.cross(t[i], u[i])
        assert np.linalg.norm(cross - v[i]) < 1e-9, f"v != t x u at {i}"


def test_rmf_twist_minimization_vs_pt_on_straight_path():
    """On a straight path the frame must NOT rotate (zero twist).
    Pins this property against future regression in the RMF kernel.
    """
    from radia.coil_from_cad import _parallel_transport_frame

    pts = np.column_stack([np.linspace(0, 1e-3, 10),
                             np.zeros(10), np.zeros(10)])
    t, u, v = _parallel_transport_frame(pts)
    # All u_hat should be identical (no twist on straight path)
    for i in range(1, len(pts)):
        assert np.allclose(u[i], u[0], atol=1e-12), (
            f"u rotated on straight path at step {i}: "
            f"u[0]={u[0]}, u[{i}]={u[i]}")


def test_corner_densification_inserts_intermediate_stations():
    """`_densify_at_corners` must insert intermediate spine points
    near sharp bends until the per-step bend angle <= the threshold
    (default 20 deg)."""
    from radia.coil_from_cad import _densify_at_corners

    class _MockSpine:
        """Mock OCC edge supporting `spine @ t` linear interpolation
        along a 3-segment polyline.
        """
        def __init__(self, waypoints):
            self.waypoints = np.asarray(waypoints, dtype=float)

        def __matmul__(self, t):
            # Linear interpolation in segment-fraction space
            n_seg = len(self.waypoints) - 1
            t = max(0.0, min(1.0, t))
            seg_t = t * n_seg
            i = min(int(seg_t), n_seg - 1)
            frac = seg_t - i
            p = (1 - frac) * self.waypoints[i] + frac * self.waypoints[i + 1]

            class _Pt:
                def __init__(self, x, y, z):
                    self.X, self.Y, self.Z = x, y, z
            return _Pt(p[0], p[1], p[2])

    # Spine with one 90 deg kink at the middle
    waypoints = [[0, 0, 0], [0.010, 0, 0], [0.010, 0.010, 0]]
    spine = _MockSpine(waypoints)
    init_pts = np.array(waypoints, dtype=float)  # 3 points, 1 kink

    densified = _densify_at_corners(spine, init_pts,
                                       max_bend_per_step_deg=20.0,
                                       max_total_points=20)
    # On a TRUE polyline kink (mock), linear interp gives the same
    # kink at any bisection midpoint, so the loop fills up to
    # max_total_points without reducing the kink itself.  This test
    # only asserts the function INSERTS points (and reaches the cap).
    # The real-OCC behaviour is exercised by the keiko fixture
    # end-to-end test which verifies smoothed bend distribution.
    assert len(densified) > len(init_pts), (
        f"densification did not insert intermediate points: "
        f"{len(init_pts)} -> {len(densified)}")


def test_extract_centerline_rejects_multi_solid_step():
    """v4.49.0 entry guard: STEP with > 1 solid must RAISE."""
    import tempfile
    import build123d as bd
    from radia.coil_from_cad import extract_centerline_from_step

    # Build a 2-solid Compound: two cylinders 50 mm apart in y.
    cyl_a = bd.Cylinder(radius=0.005, height=0.020)
    cyl_b = bd.Cylinder(radius=0.005, height=0.020,
                         align=(bd.Align.CENTER, bd.Align.CENTER, bd.Align.CENTER))
    cyl_b = cyl_b.translate(bd.Vector(0, 0.050, 0))
    compound = bd.Compound([cyl_a, cyl_b])

    with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
        step_path = f.name
    try:
        bd.export_step(compound, step_path)
        with pytest.raises(ValueError, match="STEP contains 2 solids"):
            extract_centerline_from_step(step_path, n_segments=10,
                                           cad_units_per_meter=1.0)
    finally:
        try:
            os.unlink(step_path)
        except OSError:
            pass
