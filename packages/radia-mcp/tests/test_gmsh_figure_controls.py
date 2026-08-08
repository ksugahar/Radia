"""Figure-control layer: camera presets, colour scale, glyphs, clip, axes.

The gmsh OPTIONS behind these have always existed -- what was missing was
a way to reach them without knowing their names.  These tests lock both
the naming layer (unknown names raise, with the valid list) and the two
behaviours that are easy to get silently wrong:

  * the axis camera presets.  MEASURED, not taken from another package's
    convention: a scene-invariant marker rig (identical bounding markers
    fix the zoom) shows which axis projects to a point, and an off-axis
    marker distinguishes "the opposite side" (horizontal mirror, up
    preserved) from "the same side upside down".  (180, 0, 0) is the
    latter, which is why -z is (0, 180, 0).

  * a custom colour range.  gmsh autoscales EVERY view to its own
    extrema, so two panels of the same quantity are not comparable until
    RangeType=2 + CustomMin/Max pins them to one scale.
"""

from __future__ import annotations

import numpy as np
import pytest

from radia_mcp.gmsh.post_display import CAMERA_PRESETS
from radia_mcp.gmsh.render import (_build_annotations, _build_axes,
                                   _build_clip, _build_color, _build_glyphs,
                                   render_png)

pytest.importorskip("netgen.occ", reason="netgen.occ not installed")
Image = pytest.importorskip("PIL.Image", reason="Pillow not installed")

L, E = 0.10, 0.0015
AXIS_PRESETS = ("+x", "-x", "+y", "-y", "+z", "-z")


# --------------------------------------------------------------------
# naming layer: unknown names raise with the valid list
# --------------------------------------------------------------------

def test_color_names_and_range_validated():
    cfg = _build_color({"range": [0.0, 1.5], "style": "discrete",
                        "intervals": 8, "log": True, "format": "%.2f"})
    assert cfg["range"] == [0.0, 1.5]
    assert cfg["intervals_type"] == 3.0      # discrete
    assert cfg["scale_type"] == 2.0          # log
    assert _build_color({"range": "shared"})["range"] == "shared"
    with pytest.raises(ValueError, match="continuous, discrete, iso"):
        _build_color({"style": "rainbow"})
    with pytest.raises(ValueError, match="increasing"):
        _build_color({"range": [1.0, 0.0]})


def test_glyph_and_clip_and_axes_validated():
    g = _build_glyphs({"type": "arrow3d", "sampling": 12, "size_max": 60})
    assert g["vector_type"] == 4.0 and g["sampling"] == 12
    with pytest.raises(ValueError, match="glyph type"):
        _build_glyphs({"type": "spikes"})
    with pytest.raises(ValueError, match="sampling must be >= 1"):
        _build_glyphs({"sampling": 0})

    c = _build_clip([{"normal": [1, 0, 0], "offset": -0.01,
                      "apply_to": ["views"]}])
    assert c[0]["plane"] == [1.0, 0.0, 0.0, -0.01]
    assert c[0]["targets"] == ["View[0].Clip"]
    with pytest.raises(ValueError, match="clip target"):
        _build_clip([{"normal": [1, 0, 0], "apply_to": ["walls"]}])
    with pytest.raises(ValueError, match="at most 6"):
        _build_clip([{"normal": [1, 0, 0]}] * 7)

    a = _build_axes({"mode": "box", "labels": ["x [m]", "y [m]", "z [m]"],
                     "tics": [3, 3, 3]})
    assert a["mode"] == 1.0 and a["labels"][0] == "x [m]"
    assert _build_axes(True)["mode"] == 1.0
    assert _build_axes(None) is None
    with pytest.raises(ValueError, match="axes mode"):
        _build_axes({"mode": "cage"})
    with pytest.raises(ValueError, match="three entries"):
        _build_axes({"labels": ["x", "y"]})

    ann = _build_annotations(["hello", {"text": "there", "x": 5, "y": -20}])
    assert [a["text"] for a in ann] == ["hello", "there"]
    with pytest.raises(ValueError, match="needs a 'text'"):
        _build_annotations([{"x": 1}])


# --------------------------------------------------------------------
# camera presets: measured with a scene-invariant marker rig
# --------------------------------------------------------------------

def _marker_rig(tmp_path, marker=None):
    """8 tiny corner cubes (identical in every file -> identical zoom),
    optionally plus one marker box."""
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
    mk = _marker_rig(tmp_path, marker)
    opts = dict(camera_preset=preset, width=340, height=340,
                options={"General.SmallAxes": 0})
    rp = tmp_path / f"r_{tag}.png"
    mp = tmp_path / f"m_{tag}.png"
    for src, out in ((ref, rp), (mk, mp)):
        res = render_png(src, out, **opts)
        assert res.get("ok"), res
    diff = np.abs(_ink(mp) - _ink(rp)) > 30
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
    # and it is genuinely foreshortened, not merely the smallest
    others = [v for k, v in lengths.items() if k != axis]
    assert lengths[axis] < 0.25 * min(others), lengths


@pytest.mark.parametrize("plus,marker", [
    ("+x", ((-0.01, 0.05, -0.01), (0.01, 0.07, 0.01))),
    ("+y", ((0.05, -0.01, -0.01), (0.07, 0.01, 0.01))),
    ("+z", ((0.05, -0.01, -0.01), (0.07, 0.01, 0.01))),
])
def test_minus_preset_is_the_opposite_side_not_an_upside_down_flip(
        tmp_path, plus, marker):
    """Seeing the scene from the other side mirrors HORIZONTALLY and
    keeps up pointing up.  (180, 0, 0) flips vertically instead -- the
    trap this locks out."""
    minus = "-" + plus[1]
    xs_p, ys_p = _marker_pixels(tmp_path, plus, marker, f"p{plus[1]}")
    xs_m, ys_m = _marker_pixels(tmp_path, minus, marker, f"m{plus[1]}")
    assert abs(ys_p.mean() - ys_m.mean()) < 8, (
        f"{minus}: marker changed height -> vertical flip, not the "
        f"opposite side")
    assert abs(xs_p.mean() - xs_m.mean()) > 20, (
        f"{minus}: marker did not move sideways -> same view as {plus}")


def test_camera_preset_table_keeps_the_legacy_names():
    for legacy in ("z_up_xz_from_positive_y", "positive_y_oblique",
                   "front_xz", "custom"):
        assert legacy in CAMERA_PRESETS
    for preset in AXIS_PRESETS + ("iso",):
        assert preset in CAMERA_PRESETS


# --------------------------------------------------------------------
# colour range: the comparability lever
# --------------------------------------------------------------------

def test_custom_colour_range_is_reported_and_shared_is_computed(tmp_path):
    from netgen.occ import Box, Pnt

    src = tmp_path / "b.step"
    Box(Pnt(0, 0, 0), Pnt(0.1, 0.05, 0.02)).WriteStep(str(src))
    res = render_png(src, tmp_path / "c.png", width=280, height=260,
                     camera_preset="+z",
                     color={"range": [0.0, 0.15], "style": "discrete",
                            "intervals": 5})
    assert res.get("ok"), res
    assert res["color_range"] == [0.0, 0.15]
