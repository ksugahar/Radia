"""Surrogate tests for CoilBuilder closure helpers (close_direct / close_loop).

Guards: explicit direct-segment loop closure and mode validation without
requiring native _radia_pybind or CAD binaries.
Assumes: CoilBuilder pure-Python path geometry in src/radia/coil_builder.py.
"""

from __future__ import annotations

from pathlib import Path
import sys
import types

import numpy as np
import pytest


@pytest.fixture
def coil_modules(monkeypatch):
    radia_src = Path(__file__).resolve().parents[1] / "src" / "radia"
    pkg = types.ModuleType("radia")
    pkg.__path__ = [str(radia_src)]
    monkeypatch.setitem(sys.modules, "radia", pkg)
    names = ("radia.coil_builder", "radia.coil_profile", "radia.coil_from_step")
    for name in names:
        monkeypatch.setitem(sys.modules, name, None)
        del sys.modules[name]
    from radia import coil_builder, coil_profile, coil_from_step

    return coil_builder.CoilBuilder, coil_profile.RectProfile, coil_from_step


def test_close_direct_closes_open_rectangular_path(coil_modules):
    CoilBuilder, RectProfile, _ = coil_modules
    cb = CoilBuilder(current=100.0)
    cb.set_start([0.0, 0.0, 0.0])
    cb.set_profile(RectProfile(0.005, 0.005))
    cb.add_straight(0.10)
    cb.add_arc(radius=0.01, arc_angle=90.0)
    cb.add_straight(0.05)
    cb.add_arc(radius=0.01, arc_angle=90.0)
    cb.add_straight(0.10)
    assert cb.gap > 0.05
    assert not cb.is_closed

    cb.close_loop(mode="direct")
    assert cb.is_closed
    assert cb.gap < 1e-10
    assert len(cb.segments) == 6


def test_close_loop_rejects_unknown_mode(coil_modules):
    CoilBuilder, RectProfile, _ = coil_modules
    cb = CoilBuilder(current=10.0)
    cb.set_start([0.0, 0.0, 0.0])
    cb.set_profile(RectProfile(0.005, 0.005))
    cb.add_straight(0.05)
    with pytest.raises(ValueError):
        cb.close_loop(mode="magic_auto")


def test_step_entry_requires_explicit_current_and_forwards_it(coil_modules, monkeypatch, tmp_path):
    builder_class, profile_class, step = coil_modules
    path = tmp_path / "coil.step"
    extracted = object()
    original_segments = []
    seen = []

    def build(result, current):
        seen.append((result, current))
        builder = builder_class(current=current)
        builder.set_start([0.0, 0.0, 0.0])
        builder.set_profile(profile_class(0.005, 0.005))
        builder.add_straight(0.05)
        return builder, original_segments

    monkeypatch.setattr(step, "extract_centerline", lambda name: extracted if name == str(path) else None)
    monkeypatch.setattr(step, "to_coil_builder", build)
    with pytest.raises(TypeError, match="current_a"):
        step.coil_from_step_file(path)
    with pytest.raises(ValueError, match="finite"):
        step.coil_from_step_file(path, current_a=float("nan"))
    assert not seen
    builder, segments, summary = step.coil_from_step_file(path, current_a=-12, close_mode="direct")
    assert seen == [(extracted, -12.0)]
    assert segments is original_segments
    assert builder.is_closed and summary["segment_count"] == 2
    assert summary["current_a"] == -12.0


def test_missing_build123d_propagates_without_step_fallback(coil_modules, monkeypatch, tmp_path):
    _, _, step = coil_modules
    path = tmp_path / "coil.step"
    path.touch()
    monkeypatch.setitem(sys.modules, "build123d", None)
    with pytest.raises(ModuleNotFoundError, match="build123d"):
        step.load_step_solid(path)
    with pytest.raises(FileNotFoundError):
        step.load_step_solid(tmp_path / "absent.step")


def test_existing_step_loader_reads_fixture_via_build123d(coil_modules):
    pytest.importorskip("build123d")
    pytest.importorskip("netgen.occ")
    _, _, step = coil_modules
    fixture = (Path(__file__).resolve().parents[1] /
               "packages/radia-mcp/tests/mcp_server/fixtures/overlay_box.step")
    solid = step.load_step_solid(str(fixture))
    assert len(list(solid.solids)) == 1
    assert solid.mass > 0
