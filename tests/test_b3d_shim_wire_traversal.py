"""Wire traversal contracts for the pythonocc shim.

``_sample_face_perimeter_in_pt_frame`` walks a face's outer wire edge by
edge and samples each one by arc length.  That only produces a simple
outline if the shim hands back the edges in wire order, oriented along
the traversal, and measures arc length in absolute CAD units.  The shim
did none of the three: a boolean section face of a revolved beak fin
came back with 5 of 7 junctions discontinuous, and any edge shorter than
one CAD unit was sampled at ``s * length`` instead of ``s`` because an
absolute distance inside [0, 1] was reinterpreted as a fraction.  The
beak tip arc is 0.685 mm, so the sampled section self-intersected and
the fin analyser rejected it.

build123d is the parity reference here: these same checks pass on its
faces, which is why only the shim-fed route was affected.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).parent / "coil_from_cad" / "fixtures"
CURVED = FIXTURES / "beak_fin_curved.step"


def _section_faces():
    from radia._b3d_shim import import_step
    from radia.coil_from_cad import _section_faces_from_solid

    solid = import_step(str(CURVED)).solids()[0]
    sections = _section_faces_from_solid(solid, 1000.0, 7)
    assert sections is not None, "curved beak fin failed to section"
    return sections[1]


def _traversal_gaps(face):
    from radia._b3d_shim import PositionMode

    edges = face.outer_wire().edges()
    ends = None
    gaps = []
    for edge in edges:
        start = edge.position_at(0.0, position_mode=PositionMode.LENGTH)
        stop = edge.position_at(float(edge.length),
                                position_mode=PositionMode.LENGTH)
        start = np.array([start.X, start.Y, start.Z])
        if ends is not None:
            gaps.append(float(np.linalg.norm(start - ends)))
        ends = np.array([stop.X, stop.Y, stop.Z])
    return edges, gaps


def test_outer_wire_edges_are_contiguous_in_traversal_order():
    for face in _section_faces():
        edges, gaps = _traversal_gaps(face)
        assert len(edges) >= 3
        assert max(gaps) < 1e-6, (
            f"outer wire is discontinuous by {max(gaps):.3g} CAD units; "
            f"edges are not in traversal order or not oriented")


def test_position_at_length_is_absolute_arc_length_on_a_short_edge():
    """The tip arc is 0.685 mm long, i.e. inside [0, 1] CAD units."""
    from radia._b3d_shim import PositionMode

    face = _section_faces()[3]
    edge = min(face.outer_wire().edges(), key=lambda e: float(e.length))
    length = float(edge.length)
    assert length < 1.0, "fixture no longer has a sub-unit tip arc"

    end = edge.position_at(length, position_mode=PositionMode.LENGTH)
    expected = edge.end_point()
    assert np.allclose([end.X, end.Y, end.Z],
                       [expected.X, expected.Y, expected.Z], atol=1e-9)

    # Arc length must advance uniformly, not quadratically.
    midpoint = edge.position_at(0.5 * length, position_mode=PositionMode.LENGTH)
    start = edge.start_point()
    first = np.linalg.norm([midpoint.X - start.X, midpoint.Y - start.Y,
                            midpoint.Z - start.Z])
    second = np.linalg.norm([end.X - midpoint.X, end.Y - midpoint.Y,
                             end.Z - midpoint.Z])
    assert first == pytest.approx(second, rel=1e-6)


def test_sampled_section_outline_is_simple_and_closes():
    """End to end: the perimeter sampler must return a simple polygon."""
    from radia.coil_from_cad import _sample_face_perimeter_in_pt_frame
    from radia.fin_section import analyze_section

    faces = _section_faces()
    for face in faces:
        centre = face.center()
        origin = np.array([centre.X, centre.Y, centre.Z])
        normal = face.normal_at(centre)
        n = np.array([normal.X, normal.Y, normal.Z], dtype=float)
        n /= np.linalg.norm(n)
        seed = np.array([0.0, 0.0, 1.0])
        if abs(seed @ n) > 0.9:
            seed = np.array([1.0, 0.0, 0.0])
        u = seed - (seed @ n) * n
        u /= np.linalg.norm(u)
        v = np.cross(n, u)

        uv = _sample_face_perimeter_in_pt_frame(face, origin, u, v, 512)
        step = np.linalg.norm(np.roll(uv, -1, axis=0) - uv, axis=1)
        # Equal-arc sampling: no sample may jump across the section.
        assert step.max() < 8.0 * step.mean(), (
            f"outline jumps {step.max() / step.mean():.1f}x the mean step")
        # The fin analyser is the real consumer and rejects a
        # self-intersecting outline outright.
        assert analyze_section(uv * 1e-3).primary is not None
