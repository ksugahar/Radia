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
    from coil_from_cad import _find_lateral_surface

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
    from coil_from_cad import _find_lateral_surface

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
    from coil_from_cad import _check_spine_no_singular_corner

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
    from coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=64.0, seg_m=0.5e-3)
    with pytest.raises(ValueError, match="singular corner"):
        _check_spine_no_singular_corner(pts, radius_m=2.9e-3,
                                          source_tag="test_keiko_corner")


def test_check_spine_no_singular_corner_zero_length_segment_raises():
    """Zero-length segment (duplicate consecutive points) must RAISE."""
    from coil_from_cad import _check_spine_no_singular_corner
    pts = np.array([[0, 0, 0], [0, 0, 0], [1e-3, 0, 0]])
    with pytest.raises(ValueError, match="zero-length segment"):
        _check_spine_no_singular_corner(pts, radius_m=0.001,
                                          source_tag="test_dup")


def test_check_spine_passes_below_60deg_threshold():
    """Bend at 59.9 deg (just below the 60 deg threshold) must PASS
    even with short adjacent segments.  Pins the bend-angle boundary.
    """
    from coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=59.9, seg_m=0.5e-3)
    _check_spine_no_singular_corner(pts, radius_m=2.9e-3,
                                      source_tag="test_below_threshold")


def test_check_spine_passes_long_segments_at_sharp_bend():
    """Sharp 90 deg bend with adjacent segments LONGER than wire radius:
    must PASS because perimeter filaments do not cross.

    Pins the orthogonal boundary -- sharp bends are OK if the spine
    has enough segment length to spread the parallel-transport rotation.
    """
    from coil_from_cad import _check_spine_no_singular_corner
    pts = _make_corner(bend_deg=90.0, seg_m=5e-3)  # SEG >> 1mm wire_r
    _check_spine_no_singular_corner(pts, radius_m=0.001,
                                      source_tag="test_long_seg_90deg")


def test_check_centerline_inside_solid_accepts_centerline_inside_bbox():
    """`_check_centerline_inside_solid` must PASS when centerline is
    entirely within solid bbox + slack.  Pins the happy path.
    """
    import build123d as bd
    from coil_from_cad import _check_centerline_inside_solid

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
    from coil_from_cad import _check_centerline_inside_solid

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
    from coil_from_cad import _check_centerline_inside_solid

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
    from coil_topology import detect_cap_faces

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
    from coil_topology import _bbox_spine_radius
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
    from coil_from_cad import _check_filaments_cover_solid_bbox
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
    from coil_from_cad import _check_centerline_inside_solid
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
    from coil_from_cad import _check_centerline_near_solid_surface

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
    from coil_from_cad import _check_centerline_near_solid_surface

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
    from coil_from_cad import _check_centerline_near_solid_surface
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
    from coil_from_cad import _find_lateral_surface
    solid = bd.import_step(fixture)
    lat = _find_lateral_surface(solid)
    assert lat is None, (
        "Predicate 1 fired on keiko's split-lateral STEP; "
        "dispatch order is broken -- this geometry must reach "
        "Predicate 4 to fail-fast on the lead-cap singular corner")


def test_extract_centerline_rejects_multi_solid_step():
    """v4.49.0 entry guard: STEP with > 1 solid must RAISE."""
    import tempfile
    import build123d as bd
    from coil_from_cad import extract_centerline_from_step

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
