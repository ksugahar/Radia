"""Fast configuration/API contracts for the Gmsh display layer.

Image generation, STEP meshing, ray casting, and LIC numerical evidence live
under ``validation_test/radia_mcp``.  This package test intentionally imports
no solver.

CI VISIBILITY (2026-08-11): ``radia_mcp.gmsh.raster`` pulls in numpy, and the
minimal-dep matrix has no numpy, so conftest's ``collect_ignore`` drops any
test MODULE that imports it -- which silently hid every pure-python contract
below from CI.  It hid a real red test for two commits: the committed
assertion still expected the pre-8d64df1f0 clip target.  The raster import is
therefore local to the one test that needs it, so the numpy-free contracts
keep running everywhere.
"""
from __future__ import annotations

import inspect

import pytest

from radia_mcp.gmsh.post_display import CAMERA_PRESETS
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
    # "views" is a SENTINEL the worker expands to every loaded view
    # (render.py: `... if target == "__all_views__" else [target]`).
    # It replaced the literal "View[0].Clip" in 8d64df1f0 so that a clip
    # applies to all views, not only the first; this assertion was left
    # stale by that commit and CI could not see it (see the module
    # docstring).
    assert clip[0]["targets"] == ["__all_views__"]
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
    pytest.importorskip("numpy", reason="raster is the numpy lane")
    from radia_mcp.gmsh.raster import lic, volume_raycast

    raycast = inspect.signature(volume_raycast).parameters
    streamline = inspect.signature(lic).parameters
    assert {"view", "view_dir", "image_size", "n_steps", "step_files"} \
        <= set(raycast)
    assert {"view", "plane", "resolution", "kernel", "step_files"} \
        <= set(streamline)


def test_raster_missing_input_fails_before_solver_io(tmp_path):
    pytest.importorskip("numpy", reason="raster is the numpy lane")
    from radia_mcp.gmsh.raster import lic, volume_raycast

    missing = tmp_path / "missing.msh"
    raycast = volume_raycast(missing, grid=4)
    streamline = lic(missing, resolution=32)
    assert not raycast["ok"] and "file not found" in raycast["error"]
    assert not streamline["ok"] and "file not found" in streamline["error"]
