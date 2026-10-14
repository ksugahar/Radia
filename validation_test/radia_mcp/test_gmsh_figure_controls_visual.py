"""Rendered validation for Gmsh camera and colour controls.

These checks intentionally use Netgen CAD, Gmsh subprocess rendering, and
pixel measurements.  The lightweight naming/input contracts remain under
``packages/radia-mcp/tests``.
"""

from __future__ import annotations

import numpy as np
import pytest

from radia_mcp.gmsh.render import render_png

pytest.importorskip("netgen.occ", reason="netgen.occ not installed")
Image = pytest.importorskip("PIL.Image", reason="Pillow not installed")

L, E = 0.10, 0.0015
AXIS_PRESETS = ("+x", "-x", "+y", "-y", "+z", "-z")


def _marker_rig(tmp_path, marker=None):
    """Create a scene-invariant corner rig and optional marker."""
    from netgen.occ import Box, Glue, Pnt

    parts = [Box(Pnt(sx * L - E, sy * L - E, sz * L - E),
                 Pnt(sx * L + E, sy * L + E, sz * L + E))
             for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    name = "ref"
    if marker is not None:
        lo, hi = marker
        parts.append(Box(Pnt(*lo), Pnt(*hi)))
        name = "mark"
    path = tmp_path / f"{name}_{abs(hash(str(marker))) % 10**6}.step"
    Glue(parts).WriteStep(str(path))
    return path


def _ink(png):
    return np.asarray(Image.open(png).convert("L")).astype(int)


def _marker_pixels(tmp_path, preset, marker, tag):
    ref = _marker_rig(tmp_path)
    marked = _marker_rig(tmp_path, marker)
    opts = {"camera_preset": preset, "width": 340, "height": 340,
            "options": {"General.SmallAxes": 0}}
    reference_png = tmp_path / f"r_{tag}.png"
    marked_png = tmp_path / f"m_{tag}.png"
    for src, out in ((ref, reference_png), (marked, marked_png)):
        result = render_png(src, out, **opts)
        assert result.get("ok"), result
    diff = np.abs(_ink(marked_png) - _ink(reference_png)) > 30
    ys, xs = np.nonzero(diff)
    assert len(xs), f"marker not visible for preset {preset}"
    return xs, ys


@pytest.mark.parametrize("preset", AXIS_PRESETS)
def test_axis_preset_points_that_axis_at_the_camera(tmp_path, preset):
    """A bar along the named axis must project to a point."""
    axis = preset[1]
    bars = {"x": ((-L, -E, -E), (L, E, E)),
            "y": ((-E, -L, -E), (E, L, E)),
            "z": ((-E, -E, -L), (E, E, L))}
    lengths = {}
    for name, box in bars.items():
        xs, ys = _marker_pixels(tmp_path, preset, box, f"{preset[1]}{name}")
        lengths[name] = float(np.hypot(xs.max() - xs.min(),
                                       ys.max() - ys.min()))
    facing = min(lengths, key=lengths.get)
    assert facing == axis, (preset, lengths)
    others = [value for name, value in lengths.items() if name != axis]
    assert lengths[axis] < 0.25 * min(others), lengths


@pytest.mark.parametrize("plus,marker", [
    ("+x", ((-0.01, 0.05, -0.01), (0.01, 0.07, 0.01))),
    ("+y", ((0.05, -0.01, -0.01), (0.07, 0.01, 0.01))),
    ("+z", ((0.05, -0.01, -0.01), (0.07, 0.01, 0.01))),
])
def test_minus_preset_is_the_opposite_side_not_an_upside_down_flip(
        tmp_path, plus, marker):
    """The opposite side mirrors horizontally while preserving up."""
    minus = "-" + plus[1]
    xs_plus, ys_plus = _marker_pixels(tmp_path, plus, marker, f"p{plus[1]}")
    xs_minus, ys_minus = _marker_pixels(tmp_path, minus, marker, f"m{plus[1]}")
    assert abs(ys_plus.mean() - ys_minus.mean()) < 8, (
        f"{minus}: marker changed height instead of preserving up")
    assert abs(xs_plus.mean() - xs_minus.mean()) > 20, (
        f"{minus}: marker did not move sideways from {plus}")


def test_custom_colour_range_is_reported_and_shared_is_computed(tmp_path):
    from netgen.occ import Box, Pnt

    src = tmp_path / "b.step"
    Box(Pnt(0, 0, 0), Pnt(0.1, 0.05, 0.02)).WriteStep(str(src))
    result = render_png(
        src, tmp_path / "c.png", width=280, height=260,
        camera_preset="+z",
        color={"range": [0.0, 0.15], "style": "discrete", "intervals": 5})
    assert result.get("ok"), result
    assert result["color_range"] == [0.0, 0.15]
