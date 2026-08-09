"""Fast configuration/API contracts for the Gmsh display layer.

Image generation, STEP meshing, ray casting, and LIC numerical evidence live
under ``validation_test/radia_mcp``.  This package test intentionally imports
no solver.
"""
from __future__ import annotations

import inspect

import pytest

from radia_mcp.gmsh.post_display import CAMERA_PRESETS
from radia_mcp.gmsh.raster import lic, volume_raycast
from radia_mcp.gmsh.render import (
    _build_annotations,
    _build_axes,
    _build_clip,
    _build_color,
    _build_glyphs,
)


def test_color_names_and_range_validated():
    cfg = _build_color({
        "range": [0.0, 1.5], "style": "discrete",
        "intervals": 8, "log": True, "format": "%.2f",
    })
    assert cfg["range"] == [0.0, 1.5]
    assert cfg["intervals_type"] == 3.0
    assert cfg["scale_type"] == 2.0
    assert _build_color({"range": "shared"})["range"] == "shared"
    with pytest.raises(ValueError, match="continuous, discrete, iso"):
        _build_color({"style": "rainbow"})
    with pytest.raises(ValueError, match="increasing"):
        _build_color({"range": [1.0, 0.0]})


def test_glyph_clip_axes_and_annotation_options_are_validated():
    glyph = _build_glyphs({"type": "arrow3d", "sampling": 12,
                           "size_max": 60})
    assert glyph["vector_type"] == 4.0 and glyph["sampling"] == 12
    with pytest.raises(ValueError, match="glyph type"):
        _build_glyphs({"type": "spikes"})
    with pytest.raises(ValueError, match="sampling must be >= 1"):
        _build_glyphs({"sampling": 0})

    clip = _build_clip([{
        "normal": [1, 0, 0], "offset": -0.01, "apply_to": ["views"],
    }])
    assert clip[0]["plane"] == [1.0, 0.0, 0.0, -0.01]
    assert clip[0]["targets"] == ["View[0].Clip"]
    with pytest.raises(ValueError, match="clip target"):
        _build_clip([{"normal": [1, 0, 0], "apply_to": ["walls"]}])

    axes = _build_axes({
        "mode": "box", "labels": ["x [m]", "y [m]", "z [m]"],
        "tics": [3, 3, 3],
    })
    assert axes["mode"] == 1.0 and axes["labels"][0] == "x [m]"
    assert _build_axes(True)["mode"] == 1.0

    annotations = _build_annotations([
        "hello", {"text": "there", "x": 5, "y": -20},
    ])
    assert [row["text"] for row in annotations] == ["hello", "there"]


def test_camera_preset_names_remain_stable():
    for legacy in ("z_up_xz_from_positive_y", "positive_y_oblique",
                   "front_xz", "custom"):
        assert legacy in CAMERA_PRESETS
    for preset in ("+x", "-x", "+y", "-y", "+z", "-z", "iso"):
        assert preset in CAMERA_PRESETS


def test_raster_public_signatures_keep_core_controls():
    raycast = inspect.signature(volume_raycast).parameters
    streamline = inspect.signature(lic).parameters
    assert {"view", "view_dir", "image_size", "n_steps", "step_files"} \
        <= set(raycast)
    assert {"view", "plane", "resolution", "kernel", "step_files"} \
        <= set(streamline)


def test_raster_missing_input_fails_before_solver_io(tmp_path):
    missing = tmp_path / "missing.msh"
    raycast = volume_raycast(missing, grid=4)
    streamline = lic(missing, resolution=32)
    assert not raycast["ok"] and "file not found" in raycast["error"]
    assert not streamline["ok"] and "file not found" in streamline["error"]
