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
    from coil_from_cad import extract_centerline_from_step

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


def test_keiko_outsideline_filaments_pass_bbox_check():
    """End-to-end: `filaments_from_step` must NOT raise the bbox-coverage check."""
    from coil_from_cad import filaments_from_step

    assert os.path.exists(KEIKO_STEP), KEIKO_STEP
    # Must not raise -- the v4.47.x fail-fast check was the symptom of
    # the dispatch bug.  After v4.48.0 the topology_spine path falls
    # through and open_spine produces a covering centerline.
    topo = filaments_from_step(KEIKO_STEP, sigma=5.8e7, n_peri=16)
    assert "filament_paths" in topo

    paths = np.asarray(topo["filament_paths"])
    pts = paths.reshape(-1, 3)
    # Filament bbox must include the lead.
    assert pts[:, 1].max() >= 0.048, (
        f"filament y_max={pts[:, 1].max()*1e3:.1f} mm < 48 mm "
        "(lead not traced)")
    # Filament bbox must NOT overshoot the arc radius.  Allow one
    # wire-radius (3 mm) slack on each side.
    assert pts[:, 0].max() < 0.037 and pts[:, 0].min() > -0.037, (
        f"filament x bbox [{pts[:, 0].min()*1e3:.1f},{pts[:, 0].max()*1e3:.1f}] "
        "mm overshoots conductor (bbox-fallback regression)")


def test_filaments_from_step_no_path_points_kwarg():
    """`path_points_m` was removed in v4.48.0 (STEP-only centerline policy)."""
    import inspect
    from coil_from_cad import filaments_from_step

    sig = inspect.signature(filaments_from_step)
    assert "path_points_m" not in sig.parameters, (
        "filaments_from_step still accepts path_points_m -- the v4.48.0 "
        "policy is STEP auto-detect as the single source of truth")
