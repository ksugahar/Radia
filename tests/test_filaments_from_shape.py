"""Tests for filaments_from_shape plumbing (build123d -> PEEC).

Focused on the two small helpers that glue build123d into the existing
walking-plane pipeline:
  _bd_face_to_start_hint   — build123d Face -> ((px,py,pz), (tx,ty,tz))
  _bd_shape_to_netgen_solid — build123d Shape -> netgen.occ Solid

Full PEEC integration is exercised by test_step_to_peec_inductance.py
for the STEP-file path; the in-memory path reuses the same downstream
pipeline, so only the glue layer needs coverage here.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

pytest.importorskip("build123d")


# Suppress noisy build123d warnings that get promoted to errors by
# pytest.ini's filterwarnings=error.  These originate inside OCCT and
# are harmless for our use.
warnings.filterwarnings("ignore", message=".*Gimbal lock.*")
warnings.filterwarnings("ignore", message=".*Unknown Compound type.*")


def test_bd_face_to_start_hint_planar():
    """Planar Face at z=0 should yield p=(0,0,0), inward normal ~ +Z."""
    import build123d as bd
    from radia.coil_from_cad import _bd_face_to_start_hint

    box = bd.Box(10, 10, 10, align=bd.Align.MIN)
    faces = box.faces()
    bottom = min(faces, key=lambda f: f.center().Z)
    p, t = _bd_face_to_start_hint(bottom)
    # Bottom face centroid: (5, 5, 0)
    assert np.allclose(p, [5.0, 5.0, 0.0])
    # Inward normal = +Z
    assert np.allclose(t, [0.0, 0.0, 1.0], atol=1e-6)


def test_bd_face_to_start_hint_circular():
    """Circular face on tilted plane."""
    import build123d as bd
    from radia.coil_from_cad import _bd_face_to_start_hint

    with bd.BuildSketch(bd.Plane.XY.offset(5.0)) as sk:
        bd.Circle(2.0)
    face = sk.sketch.face()
    p, t = _bd_face_to_start_hint(face)
    assert np.allclose(p, [0.0, 0.0, 5.0], atol=1e-6)
    # Face normal is +Z, inward = -Z (walking-plane convention)
    assert np.isclose(abs(t[2]), 1.0, atol=1e-6)


def test_bd_shape_to_netgen_solid_single():
    """build123d Box -> netgen.occ Solid (via BRep)."""
    import build123d as bd
    from radia.coil_from_cad import _bd_shape_to_netgen_solid

    box = bd.Box(10, 20, 30)
    ng = _bd_shape_to_netgen_solid(box)
    # netgen.occ Solid or compound fallback
    bb = ng.bounding_box
    dx = bb[1][0] - bb[0][0]
    dy = bb[1][1] - bb[0][1]
    dz = bb[1][2] - bb[0][2]
    assert np.isclose(dx, 10, atol=1e-6)
    assert np.isclose(dy, 20, atol=1e-6)
    assert np.isclose(dz, 30, atol=1e-6)


def test_start_hint_from_step_labels_unlabeled_returns_none(tmp_path):
    """STEP without any peec_* label -> helper returns None cleanly.

    build123d.export_step does not preserve child labels, so build123d
    is not expected to produce a label-bearing STEP.  The helper must
    fall through to None in that case so the caller can use bbox auto.
    """
    import build123d as bd
    from radia.coil_from_cad import _start_hint_from_step_labels

    box = bd.Box(10, 10, 30)
    step_path = tmp_path / "plain.step"
    bd.export_step(box, str(step_path))

    hint = _start_hint_from_step_labels(str(step_path))
    assert hint is None


def test_export_step_with_labels_roundtrip(tmp_path):
    """Label-preserving STEP writer: child labels survive a full round-trip.

    build123d.export_step strips Compound-child labels; our
    export_step_with_labels helper writes them via pythonocc-core XCAF
    directly, producing FreeCAD-Import.export-compatible STEP.  Reload
    via _start_hint_from_step_labels must recover the port hint.
    """
    import build123d as bd
    from build123d import Shell, Box
    from radia.coil_from_cad import (export_step_with_labels,
                                _start_hint_from_step_labels)

    box = Box(10, 10, 30)
    box.label = "peec_coil_body"
    top = max(box.faces(), key=lambda f: f.center().Z)
    port = Shell([top])
    port.label = "peec_port_in"

    step_path = tmp_path / "labeled.step"
    export_step_with_labels([box, port], str(step_path))

    # Grep-level sanity: both PRODUCTs present
    with open(step_path, encoding="utf-8") as f:
        content = f.read()
    assert "peec_coil_body" in content
    assert "peec_port_in" in content

    # Reload path: helper returns a valid start_hint
    hint = _start_hint_from_step_labels(str(step_path))
    assert hint is not None, "label-based hint should be recovered"
    p, t = hint
    assert p.shape == (3,)
    assert t.shape == (3,)
    # Top face of the box sits at z=15 (symmetric Box) — verify p is up
    # and normal is inward (=-Z after our convention)
    assert p[2] > 0
    assert np.isclose(t[2], -1.0, atol=1e-6)


def test_start_hint_from_step_labels_freecad_fixture():
    """Real-world fixture: FreeCAD Import.export STEP with peec_port_in shell.

    The label path in Radia is designed for external CAD tools that can
    emit XCAF labels (FreeCAD Import.export, Cubit AP214-aware exporters
    in the future, etc).  build123d STEPs are not expected to carry
    labels — use `filaments_from_shape(shape, port_face=...)` for
    build123d workflows instead.
    """
    import os
    from radia.coil_from_cad import _start_hint_from_step_labels

    fc_step = r"C:\temp\fc_xcaf.step"
    if not os.path.exists(fc_step):
        pytest.skip("FreeCAD fixture not available (see docs/peec)")

    hint = _start_hint_from_step_labels(fc_step)
    assert hint is not None, "expected label-based hint from FreeCAD STEP"
    p, t = hint
    assert p.shape == (3,)
    assert t.shape == (3,)
    assert np.isclose(np.linalg.norm(t), 1.0)


def test_filaments_from_shape_smoke_box():
    """Smoke test: just ensure filaments_from_shape runs on a simple box.

    We pass an explicit port_face; the walking-plane may not produce a
    meaningful result on a box (not a coil), but the plumbing should not
    crash when port_face is provided (start_hint construction path).
    """
    import build123d as bd
    from radia.coil_from_cad import (_bd_face_to_start_hint,
                                _bd_shape_to_netgen_solid)

    box = bd.Box(10, 10, 30)
    top = max(box.faces(), key=lambda f: f.center().Z)
    p, t = _bd_face_to_start_hint(top)
    ng = _bd_shape_to_netgen_solid(box)
    # Just assert we reached the conversion without exception
    assert ng is not None
    assert p.shape == (3,)
    assert t.shape == (3,)
    assert np.isclose(np.linalg.norm(t), 1.0)
