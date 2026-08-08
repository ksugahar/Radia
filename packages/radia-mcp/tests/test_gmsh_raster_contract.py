"""Fast argument contracts for the optional Gmsh rasterizers.

Pixel-level algorithm checks live in ``validation_test/radia_mcp`` because
they execute Gmsh, SciPy, Netgen/OCC, and image rendering.
"""

from __future__ import annotations

import pytest
from radia_mcp.gmsh.raster import lic, volume_raycast


@pytest.fixture
def placeholder(tmp_path):
    path = tmp_path / "field.msh"
    path.write_text("placeholder", encoding="ascii")
    return path


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"grid": 7}, "grid"),
        ({"grid": 129}, "grid"),
        ({"grid": 8.0}, "integer"),
        ({"image_size": 4096}, "image_size"),
        ({"n_steps": 4096}, "n_steps"),
        ({"alpha": float("nan")}, "alpha"),
        ({"alpha_power": -1.0}, "alpha_power"),
        ({"timeout_s": float("inf")}, "timeout_s"),
        ({"step_rel_size": 1.1}, "step_rel_size"),
        ({"step_rel_size": 0.001}, "step_rel_size"),
        ({"step_color": (0.0, 0.5, 2.0)}, "step_color"),
        ({"view_dir": [1.0, 0.0, float("nan")]}, "view_dir"),
        ({"value_range": [1.0, 1.0]}, "value_range"),
        ({"value_range": [0.0, float("inf")]}, "value_range"),
    ],
)
def test_volume_raycast_rejects_unsafe_arguments_before_parsing(
        placeholder, kwargs, message):
    with pytest.raises(ValueError, match=message):
        volume_raycast(placeholder, **kwargs)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"plane": "ab"}, "unknown plane"),
        ({"resolution": 63}, "resolution"),
        ({"resolution": 4096}, "resolution"),
        ({"kernel": 1}, "kernel"),
        ({"kernel": 512}, "kernel"),
        ({"seed": -1}, "seed"),
        ({"offset": float("nan")}, "offset"),
        ({"timeout_s": 0.0}, "timeout_s"),
        ({"step_rel_size": float("inf")}, "step_rel_size"),
    ],
)
def test_lic_rejects_unsafe_arguments_before_parsing(
        placeholder, kwargs, message):
    with pytest.raises(ValueError, match=message):
        lic(placeholder, **kwargs)
