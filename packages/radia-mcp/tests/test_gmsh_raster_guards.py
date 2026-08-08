"""Fast fail-loud guards found by the 2026-08-07 debug battery.

Two silent-wrongness bugs, both of the worst kind -- a plausible
picture over garbage settings:

* ``volume_raycast(value_range=[hi, lo])`` silently became
  ``[hi, hi + 1]``: the normalization inverted and every value clipped
  to nonsense while the figure still rendered.  The check now runs
  BEFORE the grid probe, so the 262k-probe resample is never paid for
  an input that was going to be refused.
* a clip plane with a ZERO normal was forwarded to gmsh, which accepts
  it and clips nothing/everything without a word.

These tests are deliberately light (no gmsh subprocess, no rendering)
so they stay in the fast CI lane; the render-heavy raster suites live
in validation_test/radia_mcp/.
"""

from __future__ import annotations

import pytest

from radia_mcp.gmsh.raster import volume_raycast
from radia_mcp.gmsh.render import _build_clip


def _tiny_msh(path):
    """Smallest valid v4.1 file with one node view (never probed --
    the argument checks under test fire first)."""
    path.write_text(
        "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
        "$Entities\n0 0 0 1\n1 0 0 0 1 1 1 0 0\n$EndEntities\n"
        "$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n4\n"
        "0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n"
        "$Elements\n1 1 1 1\n3 1 4 1\n1 1 2 3 4\n$EndElements\n"
        "$NodeData\n1\n\"f\"\n1\n0\n3\n0\n1\n4\n"
        "1 1.0\n2 1.0\n3 1.0\n4 1.0\n$EndNodeData\n",
        encoding="utf-8")
    return path


def test_reversed_value_range_raises_before_the_probe(tmp_path):
    f = _tiny_msh(tmp_path / "t.msh")
    with pytest.raises(ValueError, match="value_range must"):
        volume_raycast(f, tmp_path / "x.png", view="f", grid=8,
                       image_size=64, value_range=[1.0, 0.0])
    with pytest.raises(ValueError, match="value_range must"):
        volume_raycast(f, tmp_path / "x.png", view="f", grid=8,
                       image_size=64, value_range=[0.5, 0.5])


def test_zero_clip_normal_raises_in_both_spec_forms():
    with pytest.raises(ValueError, match="normal must be nonzero"):
        _build_clip([{"normal": [0, 0, 0], "offset": 0.01}])
    with pytest.raises(ValueError, match="normal must be nonzero"):
        _build_clip([{"plane": [0.0, 0.0, 0.0, 0.01]}])
    # a genuine plane still passes and assembles [a, b, c, d]
    ok = _build_clip([{"normal": [0, 1, 0], "offset": -0.02}])
    assert ok[0]["plane"] == [0.0, 1.0, 0.0, -0.02]
