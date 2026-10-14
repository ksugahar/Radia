"""Golden for the ROBUST mesh-boolean printable-former STL output.

``_write_former_stl`` (trimesh + manifold3d) handles the 3D SELF-CROSSING
single-stroke wire that the OCCT STEP route (``_write_former_step``, netgen.occ
Pipe/Fuse AND build123d sweep -- one kernel) segfaults / silently under-cuts on
(channel ratio ~0.5% or a degenerate shell).  See
``memory/sf_printable_former_cad_status.md``.

This golden encodes exactly that OCCT-breaker geometry (a multi-turn wire whose
inter-turn connectors cross other turns -> a self-overlapping tube) and asserts
the mesh route (a) stays WATERTIGHT and (b) actually cuts the channel to ~the
swept-tube volume -- the two things OCCT fails.  It calls the helper directly
(no NGSolve), so it is a fast in-process unit golden; it SKIPS if trimesh /
manifold3d are not installed (they are an optional dep of --former-stl).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


REPO = Path(__file__).resolve().parents[2]
PANELS = REPO / "src" / "radia" / "panels"
sys.path.insert(0, str(PANELS))

import calc_streamfunction as csf  # noqa: E402


def _self_crossing_wire():
    """A multi-turn single-stroke wire (MRI-gradient scale) whose inter-turn
    connectors jump ACROSS the coil -> the swept tube self-overlaps.  This is
    the exact geometry class that makes the OCCT sweep/Fuse route degenerate."""
    R = 0.15
    turns, ppt = 6, 40
    z0, z1 = -0.02, 0.02
    chain, prev = [], None
    for j in range(turns):
        t = np.linspace(0, 2 * np.pi, ppt, endpoint=False)
        z = z0 + (z1 - z0) * (j + 0.5) / turns
        loop = np.column_stack([R * np.cos(t), R * np.sin(t),
                                z + 0.002 * np.sin(3 * t)])
        if prev is not None:                            # connector = a crossing
            chain.append(np.array([prev, loop[0]]))
        chain.append(loop)
        prev = loop[-1]
    return np.vstack(chain)


def test_write_former_stl_3d_self_crossing(tmp_path):
    pytest.importorskip("trimesh")
    pytest.importorskip("manifold3d")
    chain = _self_crossing_wire()
    out = tmp_path / "former.stl"

    info = csf._write_former_stl(
        chain, wire_diam=1.0e-3, clearance=3.0e-4, margin=5.0e-3, wall=3.0e-3,
        decimate=2.5e-3, filename=str(out))

    assert out.exists() and out.stat().st_size > 0
    # the WATERTIGHT result is the whole point (OCCT gives a non-manifold shell)
    assert info["watertight"] is True
    # the channel is actually cut to ~the swept-tube volume (OCCT gives ~0);
    # mesh faceting + joint double-count keep the ratio a little under 1.
    assert 0.75 < info["channel_ratio"] < 1.10, info["channel_ratio"]
    assert info["channel_segments"] > 0
    assert info["channel_volume_m3"] > 0.0
    assert info["former_volume_m3"] < info["plate_volume_m3"]


def test_write_former_stl_decimate_caps_segments(tmp_path):
    """The max_segments cap coarsens the resample step (reported, not silent)
    so a very long wire cannot explode the mesh-boolean cost."""
    pytest.importorskip("trimesh")
    pytest.importorskip("manifold3d")
    chain = _self_crossing_wire()
    out = tmp_path / "former.stl"

    info = csf._write_former_stl(
        chain, wire_diam=1.0e-3, clearance=3.0e-4, margin=5.0e-3, wall=3.0e-3,
        decimate=0.0, filename=str(out), max_segments=120)

    assert info["channel_segments"] <= 120
    assert info["decimate_step_m"] > 0.0            # cap raised the step
    assert info["watertight"] is True


def test_write_former_stl_missing_deps_message(monkeypatch, tmp_path):
    """No-Fallback: without trimesh/manifold3d the helper raises a CLEAR,
    actionable ImportError (not a cryptic downstream failure)."""
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name in ("trimesh", "manifold3d"):
            raise ImportError(f"no module named {name}")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="pip install trimesh manifold3d"):
        csf._write_former_stl(
            np.array([[0.0, 0.0, 0.0], [0.02, 0.0, 0.001], [0.03, 0.01, 0.0]]),
            wire_diam=1.0e-3, clearance=3.0e-4, margin=5.0e-3, wall=3.0e-3,
            decimate=0.0, filename=str(tmp_path / "x.stl"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
