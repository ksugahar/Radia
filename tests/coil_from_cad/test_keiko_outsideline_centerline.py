"""Regression test for keiko's `1turn_coil_loft_outsideline.step` (mdx).

This STEP is an "arc + 2 lead-bars" loft with a 2-BSPLINE-symmetric
lateral surface.  Geometry summary (verified 2026-05-15):

  - 1 solid, 4 faces (2 BSPLINE lateral halves + 2 PLANE caps)
  - cap_a center (-12, +50, 0) mm, cap_b center (+13, +50, 0) mm
  - bbox x=+-33 mm, y=[-33, +50] mm, z=+-3 mm
  - cross-section: circular wire radius ~ 3 mm
  - "round" portion at R ~ 30 mm, two leads extending to y = +50 mm

Historical bug (pre-v4.48.1): ``extract_centerline_from_step`` had
a try/except cascade `revolution_sweep -> topology_spine -> open_spine`.
The middle path ``_centerline_from_topology_spine`` returned a planar
arc at R = 0.85 * bbox_max = 42.5 mm.  That arc:

  1. Sat OUTSIDE the actual conductor in x (overshoot 12 mm)
  2. MISSED the lead extension at y = +50 mm (gap 6 mm)

The downstream sanity check ``_check_filaments_cover_solid_bbox``
raised "Filament path does not match the conductor solid bbox" and
told the user to either regenerate the STEP or pass an explicit
centerline JSON.  But ``_centerline_from_open_spine`` (the next
cascade step) was already present in the codebase and produced the
correct centerline -- it was simply being shadowed by the eager
topology_spine path.

Fix (v4.48.1): per the "No Fallbacks -- Fail Fast, Fail Loud" policy
(CLAUDE.md), the try/except cascade was replaced with a
classification-based single dispatch in
``extract_centerline_from_step``.  OPEN coils (cap faces detected by
``coil_topology.extract_coil_topology``) route directly to
``_centerline_from_open_spine``; CLOSED coils route directly to
``_centerline_from_topology_spine`` (which now raises if called on
an OPEN coil -- programming-error indicator, NOT a soft-fallback
signal).

In the same release, the user-facing ``path_points_m`` parameter was
removed: STEP files are now the single source of truth for the
centerline.  If STEP auto-detection cannot recover a covering
centerline, the user fixes the CAD rather than papering over the
breakage with a hand-crafted JSON.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

pytest.importorskip("build123d")


FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
KEIKO_STEP = os.path.join(FIXTURE_DIR, "keiko_outsideline.step")


def test_keiko_outsideline_centerline_covers_lead():
    """STEP auto-detect must trace the lead bar (y -> +50 mm), not bypass it."""
    from radia.coil_from_cad import extract_centerline_from_step

    assert os.path.exists(KEIKO_STEP), KEIKO_STEP
    path_m, widths_m, heights_m = extract_centerline_from_step(
        KEIKO_STEP, n_segments=100, cad_units_per_meter=1.0)

    # 1. Lead bar at y=+50 mm must be traced (gap < 1.5 * wire_radius = 4.5 mm).
    y_max = float(path_m[:, 1].max())
    assert y_max >= 0.048, (
        f"centerline missed the y=+50 mm lead bar (y_max={y_max*1e3:.1f} mm); "
        "topology_spine bbox-fallback regression")

    # 2. Centerline must NOT overshoot solid bbox in x by more than slack.
    #    Solid bbox x=+-33 mm; tolerate +-2 mm (within wire radius slack).
    x_min, x_max = float(path_m[:, 0].min()), float(path_m[:, 0].max())
    assert x_min >= -0.035 and x_max <= 0.035, (
        f"centerline overshot solid bbox in x: [{x_min*1e3:.1f},{x_max*1e3:.1f}] "
        "mm; expected within +-35 mm (bbox-radius fallback regression)")

    # 3. Equiv-circle wire radius must recover the 3 mm wire.
    mean_area = float(np.mean(widths_m * heights_m))
    r_equiv = float(np.sqrt(mean_area / np.pi))
    assert 0.0025 < r_equiv < 0.0035, (
        f"recovered wire radius {r_equiv*1e3:.2f} mm out of band [2.5, 3.5] mm")


def test_keiko_outsideline_succeeds_with_adaptive_resampling(monkeypatch):
    """v4.53.0: adaptive resampling in `_centerline_from_open_spine`
    chooses n_segments such that segment_length >= 1.10 * wire_radius,
    which prevents the singular-corner check from firing on the
    lead-arc junction.  PEEC then produces a finite L_coil for
    keiko's 1turn_coil_loft_outsideline.step.

    Previously (v4.49.0 to v4.52.0) this STEP raised
    "singular corner" because hardcoded n_segments=100 on a 208 mm
    spine gave 2.08 mm segments (< 2.9 mm wire radius) and the lead-
    arc 64 deg bend tripped the singular-corner check.  Adaptive
    resampling caps n_segments at floor(spine_length / (1.10 *
    wire_radius)) = floor(208 / 3.19) = 65, giving segments >= 3.2 mm.
    """
    from radia.coil_from_cad import filaments_from_step
    import numpy as np

    assert os.path.exists(KEIKO_STEP), KEIKO_STEP
    # Exercise the current extractor.  Historical format-1 caches predate
    # the default-walker unit boundary and are intentionally rejected.
    monkeypatch.setenv("RADIA_PEEC_CACHE_DISABLE", "1")
    topo = filaments_from_step(KEIKO_STEP, sigma=5.8e7, n_peri=16)
    L = np.asarray(topo["solver"].L)
    assert np.isfinite(L).all(), (
        "L matrix has non-finite entries even after adaptive resampling")

    Z = topo["solver"].compute_port_impedance(150_000.0)
    L_nH = float(Z.imag) / (2 * np.pi * 150_000.0) * 1e9
    # keiko's coil: 1-turn 30 mm arc + 25 mm leads.  Physical band
    # [50, 300] nH (same band as the vertex-aligned replica test).
    assert 50.0 < L_nH < 300.0, (
        f"L_coil_nH = {L_nH:.2f} out of physical band [50, 300] nH "
        "for keiko's 1-turn 30 mm-radius coil + leads (v4.53.0 should "
        "have produced a finite physical inductance)")


def test_filaments_from_step_no_path_points_kwarg():
    """`path_points_m` was removed in v4.48.0 (STEP-only centerline policy)."""
    import inspect
    from radia.coil_from_cad import filaments_from_step

    sig = inspect.signature(filaments_from_step)
    assert "path_points_m" not in sig.parameters, (
        "filaments_from_step still accepts path_points_m -- the v4.48.0 "
        "policy is STEP auto-detect as the single source of truth")


def test_vertex_aligned_replica_routes_predicate_1_uv():
    """v4.48.2: a vertex-aligned single-piece BSPLINE coil routes through
    Predicate 1 (UV-map sampling) and returns a finite L_coil.

    This is the documented FIX path for keiko's outsideline.step:
    regenerate the STEP with a clean centerline so the lateral
    surface is a single dominant BSPLINE (dominance ~ 1.00), at
    which point `_find_lateral_surface` returns the single face,
    UV sampling MUST succeed (UV-closure check is part of the
    predicate), and the filament cross-sections come directly from
    the lateral surface without parallel transport -- no degenerate
    near-coincident segment pairs.

    Geometry: a 1-turn coil with 2 leads, matching keiko's bbox
    (x=+-30 mm, y=[-30, +50] mm, wire radius 2.9 mm).  Built via
    build123d sweep(circle, smooth_spline_path).
    """
    import math
    import build123d as bd
    import tempfile

    from radia.coil_from_cad import _find_lateral_surface, filaments_from_step

    WIRE_R = 0.0029
    ARC_R = 0.030
    LEAD_X = 0.012
    LEAD_Y_TIP = 0.050
    JUNCTION_Y = ARC_R * math.sqrt(1 - (LEAD_X / ARC_R) ** 2)

    # Spine waypoints: lead_a (top -> junction), arc CCW through bottom,
    # lead_b (junction -> top).  Use Spline so the resulting sweep has
    # ONE BSPLINE lateral instead of N CYLINDER segments.
    ang_a = math.degrees(math.atan2(JUNCTION_Y, -LEAD_X))
    ang_b = math.degrees(math.atan2(JUNCTION_Y,  LEAD_X))
    sweep_deg = 360.0 - (ang_a - ang_b)

    pts = [
        (-LEAD_X, LEAD_Y_TIP, 0.0),
        (-LEAD_X, (LEAD_Y_TIP + JUNCTION_Y) / 2, 0.0),
        (-LEAD_X, JUNCTION_Y, 0.0),
    ]
    n_arc = 60
    for i in range(1, n_arc):
        t = i / n_arc
        ang = math.radians(ang_a - sweep_deg * t)
        pts.append((ARC_R * math.cos(ang), ARC_R * math.sin(ang), 0.0))
    pts.extend([
        (LEAD_X, JUNCTION_Y, 0.0),
        (LEAD_X, (LEAD_Y_TIP + JUNCTION_Y) / 2, 0.0),
        (LEAD_X, LEAD_Y_TIP, 0.0),
    ])

    spine = bd.Spline(*[bd.Vector(*p) for p in pts])
    profile_plane = bd.Plane(origin=bd.Vector(*pts[0]),
                              z_dir=bd.Vector(0, -1, 0))
    profile = profile_plane * bd.Circle(WIRE_R)
    coil = bd.sweep(profile, spine)

    with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
        step_path = f.name
    try:
        bd.export_step(coil, step_path)

        # Sanity: lateral surface is single-piece BSPLINE
        solid = bd.import_step(step_path)
        lat = _find_lateral_surface(solid)
        assert lat is not None, (
            "expected single-piece BSPLINE lateral; "
            "smooth Spline sweep produced fragmented surface")

        # PEEC build must succeed and produce finite L
        topo = filaments_from_step(step_path, sigma=5.8e7, n_peri=16)
        assert topo["source"] == "step_uv", (
            f"expected source='step_uv' (Predicate 1 UV), "
            f"got {topo['source']!r}")
        L = np.asarray(topo["solver"].L)
        assert np.isfinite(L).all(), (
            "L matrix has non-finite entries even on the smooth "
            "vertex-aligned reference geometry")

        Z = topo["solver"].compute_port_impedance(150_000.0)
        L_nH = float(Z.imag) / (2 * 3.141592653589793 * 150_000.0) * 1e9
        # 1-turn coil at R~30 mm, leads ~25 mm: physically L ~ 100-200 nH
        assert 50.0 < L_nH < 300.0, (
            f"L_coil_nH={L_nH:.2f} out of physical band [50, 300] nH "
            "for 1-turn 30 mm-radius coil + leads")
    finally:
        try:
            os.unlink(step_path)
        except OSError:
            pass
