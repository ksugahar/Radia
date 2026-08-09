"""Fast contracts for named Gmsh figure controls.

The solver-backed camera and rendered-colour checks live in
``validation_test/radia_mcp/test_gmsh_figure_controls_visual.py``.
"""

from __future__ import annotations

import pytest

from radia_mcp.gmsh.post_display import CAMERA_PRESETS
from radia_mcp.gmsh.render import (_build_annotations, _build_axes,
                                   _build_clip, _build_color, _build_glyphs)

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
    with pytest.raises(ValueError, match="increasing"):
        _build_color({"range": [0.0, float("inf")]})
    with pytest.raises(ValueError, match="intervals"):
        _build_color({"intervals": 0})
    with pytest.raises(ValueError, match="alpha"):
        _build_color({"alpha": 1.1})


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
    assert c[0]["targets"] == ["__all_views__"]
    assert _build_clip([{"normal": [1, 0, 0],
                         "apply_to": "views"}])[0]["targets"] == [
                             "__all_views__"]
    with pytest.raises(ValueError, match="clip target"):
        _build_clip([{"normal": [1, 0, 0], "apply_to": ["walls"]}])
    with pytest.raises(ValueError, match="at most 6"):
        _build_clip([{"normal": [1, 0, 0]}] * 7)
    with pytest.raises(ValueError, match="nonzero"):
        _build_clip([{"normal": [0, 0, 0]}])

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


def test_camera_preset_table_keeps_the_legacy_names():
    for legacy in ("z_up_xz_from_positive_y", "positive_y_oblique",
                   "front_xz", "custom"):
        assert legacy in CAMERA_PRESETS
    for preset in AXIS_PRESETS + ("iso",):
        assert preset in CAMERA_PRESETS
