"""Rebuild-from-reference knowledge and its lint rules for radia_mcp.build123d.

Guards: the `rebuild_from_reference_step` and `assemblies_and_compounds`
claims about label export, sketch vertex frames and Face.radius still hold
for the installed build123d, and the two lint rules flag the patterns those
claims warn about without flagging the corrected forms.
Assumes: build123d and its OCP runtime are installed (skips otherwise).
"""

import pytest

build123d = pytest.importorskip("build123d")

from build123d import (  # noqa: E402
    Box, BuildSketch, Compound, Cylinder, Face, GeomType, Plane, Pos, Rectangle,
    export_step, import_step,
)

from radia_mcp.build123d.build123d_knowledge import get_build123d_documentation  # noqa: E402
from radia_mcp.build123d.rules import (  # noqa: E402
    check_positional_compound_labels, check_sketch_vertex_global_z,
)


def _child_labels(path):
    return [child.label for child in import_step(str(path)).children]


def _two_labelled_boxes():
    a = Box(1, 1, 1)
    a.label = "a"
    b = Pos(3, 0, 0) * Box(1, 1, 1)
    b.label = "b"
    return a, b


def test_children_compound_keeps_labels_and_positional_compound_drops_them(tmp_path):
    a, b = _two_labelled_boxes()
    export_step(Compound(children=[a, b]), tmp_path / "children.step")
    export_step(Compound([a, b]), tmp_path / "positional.step")

    assert _child_labels(tmp_path / "children.step") == ["a", "b"]
    assert _child_labels(tmp_path / "positional.step") == []
    solid_labels = [s.label for s in import_step(str(tmp_path / "children.step")).solids()]
    assert solid_labels == ["", ""], f"solids() labels: {solid_labels}"


def test_shared_instances_export_the_last_label(tmp_path):
    base = Box(1, 1, 1)
    first = Pos(0, 0, 0) * base
    first.label = "first"
    second = Pos(5, 0, 0) * base
    second.label = "second"
    export_step(Compound(children=[first, second]), tmp_path / "shared.step")
    assert _child_labels(tmp_path / "shared.step") == ["second", "second"]


def test_xz_sketch_vertices_are_local():
    with BuildSketch(Plane.XZ) as sk:
        Rectangle(2, 4)
    zs = {round(v.Z, 9) for v in sk.vertices()}
    ys = sorted({round(v.Y, 9) for v in sk.vertices()})
    assert zs == {0.0} and ys == [-2.0, 2.0], f"Z {zs}, Y {ys}"


def test_face_radius_is_none_on_trimmed_cylinder_surface():
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.Geom import Geom_CylindricalSurface, Geom_RectangularTrimmedSurface
    from OCP.gp import gp_Ax3

    primitive = next(f for f in Cylinder(2, 5).faces() if f.geom_type == GeomType.CYLINDER)
    assert primitive.radius == pytest.approx(2.0)

    trimmed = Geom_RectangularTrimmedSurface(
        Geom_CylindricalSurface(gp_Ax3(), 2.0), 0.0, 3.0, 0.0, 5.0)
    face = Face(BRepBuilderAPI_MakeFace(trimmed, 1e-7).Face())
    assert face.geom_type == GeomType.CYLINDER
    assert face.radius is None

    from OCP.BRepAdaptor import BRepAdaptor_Surface
    assert BRepAdaptor_Surface(face.wrapped).Cylinder().Radius() == pytest.approx(2.0)


def _lines(text):
    return [line + "\n" for line in text.strip("\n").splitlines()]


def test_sketch_vertex_global_z_rule():
    wrong = _lines("""
with BuildSketch(Plane.XZ) as sk:
    make_face()
    pick = lambda x, z: min(sk.vertices(), key=lambda p: (p.X - x) ** 2 + (p.Z - z) ** 2)
    low = sk.vertices().sort_by(Axis.Z)[0]
    tops = [v for v in sk.vertices() if v.Z > 1]
body = min(body.vertices(), key=lambda p: p.Z)
""")
    findings = check_sketch_vertex_global_z("x.py", wrong)
    assert [(f["line"], f["rule"]) for f in findings] == [
        (3, "sketch-vertex-global-z"), (4, "sketch-vertex-global-z"), (5, "sketch-vertex-global-z")]

    right = _lines("""
with BuildSketch(Plane.XZ) as sk:
    pick = lambda x, z: min(sk.vertices(), key=lambda p: (p.X - x) ** 2 + (p.Y - z) ** 2)
    top = max(sk.sketch.vertices(), key=lambda v: v.Z)
    ref = max(reference.vertices(), key=lambda v: v.Z)
    # min(sk.vertices(), key=lambda p: p.Z)
    note = "min(sk.vertices(), key=lambda p: p.Z)"
with BuildSketch(Plane.XY) as top_sk:
    low = min(top_sk.vertices(), key=lambda p: p.Z)
""")
    assert check_sketch_vertex_global_z("x.py", right) == []


def test_positional_compound_labels_rule():
    wrong = _lines("""
a.label = "a"
asm = Compound([a, b])
""")
    findings = check_positional_compound_labels("x.py", wrong)
    assert [(f["line"], f["rule"]) for f in findings] == [(2, "positional-compound-drops-labels")]

    for text in ("""
a.label = "a"
asm = Compound(children=[a, b])
""", "asm = Compound([a, b])", """
housing.label = "housing"
cutter = Compound([left, right])
# Compound([housing])
"""):
        assert check_positional_compound_labels("x.py", _lines(text)) == [], text


def test_rules_ignore_multiline_strings_and_unbound_names():
    docstring = _lines('''
def helper():
    """Example of the trap:

    with BuildSketch(Plane.XZ) as sk:
        pick = min(sk.vertices(), key=lambda p: p.Z)
    a.label = "a"
    asm = Compound([a, b])
    """
    return None
''')
    assert check_sketch_vertex_global_z("x.py", docstring) == []
    assert check_positional_compound_labels("x.py", docstring) == []

    lambda_then_other_z = _lines("""
with BuildSketch(Plane.XZ) as sk:
    pair = (min(sk.vertices(), key=lambda v: v.Y), fixture.v.Z)
    shifted = min(sk.vertices(), key=lambda v: v.Y + offset.Z)
""")
    assert check_sketch_vertex_global_z("x.py", lambda_then_other_z) == []

    name_inside_call = _lines("""
housing.label = "housing"
cutter = Compound([make_clearance_for(housing), right])
""")
    assert check_positional_compound_labels("x.py", name_inside_call) == []


def test_rules_follow_python_scopes():
    other_function = _lines("""
def labelled_part():
    part = Box(1, 1, 1)
    part.label = "outer"
    return part

def helper():
    part = Box(1, 1, 1)
    return Compound([part, Box(2, 2, 2)])
""")
    assert check_positional_compound_labels("x.py", other_function) == []

    same_function = _lines("""
def assembly():
    part = Box(1, 1, 1)
    part.label = "outer"
    return Compound([part, Box(2, 2, 2)])
""")
    assert [f["line"] for f in check_positional_compound_labels("x.py", same_function)] == [4]

    shadowed = _lines("""
with BuildSketch(Plane.XZ) as sk:
    def highest(sk):
        return max(sk.vertices(), key=lambda v: v.Z)
    pick = max(sk.vertices(), key=lambda v: v.Y + (lambda v: v.Z)(reference_vertex))
    rebound = [v.Z for v in sk.vertices() for v in other_vertices]
""")
    assert check_sketch_vertex_global_z("x.py", shadowed) == []

    in_function = _lines("""
def profile():
    with BuildSketch(Plane.XZ) as sk:
        return min(sk.vertices(), key=lambda v: v.Z)
""")
    assert [f["line"] for f in check_sketch_vertex_global_z("x.py", in_function)] == [3]


def test_documented_surface_dispatch_runs_on_each_surface_type():
    import re

    from build123d import Cone, Torus

    text = get_build123d_documentation("rebuild_from_reference_step")
    block = next(b for b in re.findall(r"```python\n(.*?)```", text, re.S)
                 if "BRepAdaptor_Surface" in b)
    cases = {GeomType.CYLINDER: (Cylinder(2, 5), "Radius", 2.0),
             GeomType.CONE: (Cone(2, 1, 3), "RefRadius", 2.0),
             GeomType.TORUS: (Torus(5, 1), "MinorRadius", 1.0)}
    for geom_type, (shape, getter, expected) in cases.items():
        face = next(f for f in shape.faces() if f.geom_type == geom_type)
        scope = {"face": face}
        exec(block, scope)
        assert getattr(scope["c"], getter)() == pytest.approx(expected), geom_type


def test_rebuild_topic_is_registered():
    text = get_build123d_documentation("rebuild_from_reference_step")
    assert text.startswith("# Rebuilding a reference STEP"), text[:80]
    assert "rebuild_from_reference_step" in get_build123d_documentation("unknown-topic")


_L_OUTLINE = [(0.0, 0.0), (5.0, 0.0), (5.0, 3.0), (4.0, 3.0), (4.0, 1.0), (0.0, 1.0)]  # CCW, re-entrant (4, 1)


def _sketch(points):
    from build123d import BuildLine, Polyline, make_face

    with BuildSketch() as sk:
        with BuildLine():
            Polyline(*points, close=True)
        make_face()
    return sk.sketch


def test_tapered_extrude_follows_winding_and_rounds_reentrant_corner():
    from build123d import extrude

    spans = {}
    for name, pts in (("ccw", _L_OUTLINE), ("cw", _L_OUTLINE[::-1])):
        solid = extrude(Pos(0, 0, 2) * _sketch(pts), 0.3, taper=30, dir=(0, 0, 1))
        bb = solid.bounding_box()
        spans[name] = (round(bb.min.Z, 6), round(bb.max.Z, 6))
        assert sum(f.geom_type == GeomType.CONE for f in solid.faces()) == 1, name
    assert spans == {"ccw": (2.0, 2.3), "cw": (1.7, 2.0)}, spans


def test_ruled_loft_of_polygons_has_bspline_sides():
    from build123d import Vector, Wire, loft

    def face(pts, z):
        return Face(Wire.make_polygon([Vector(x, y, z) for x, y in pts], close=True))

    inner = [(0.2, 0.2), (4.8, 0.2), (4.8, 2.8), (4.2, 2.8), (4.2, 0.8), (0.2, 0.8)]
    solid = loft([face(_L_OUTLINE, 0.0), face(inner, 0.3)], ruled=True)
    types = [f.geom_type for f in solid.faces()]
    assert types.count(GeomType.BSPLINE) == 6 and types.count(GeomType.PLANE) == 2, types


def test_documented_frustum_is_planar_and_has_prismatoid_volume():
    import re

    text = get_build123d_documentation("rebuild_from_reference_step")
    block = next(b for b in re.findall(r"```python\n(.*?)```", text, re.S) if "def frustum" in b)
    scope = {}
    exec(block, scope)

    def area(pts):
        return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]))

    d, h = 0.2, 0.5
    top = scope["mitre_inset"](_L_OUTLINE, d)
    expected_top = [(0.2, 0.2), (4.8, 0.2), (4.8, 2.8), (4.2, 2.8), (4.2, 0.8), (0.2, 0.8)]
    assert [c for p in top for c in p] == pytest.approx([c for p in expected_top for c in p]), top
    solid = scope["frustum"](_L_OUTLINE, 0.0, top, h)
    assert solid.is_valid
    assert all(f.geom_type == GeomType.PLANE for f in solid.faces())
    mid = scope["mitre_inset"](_L_OUTLINE, d / 2)
    expected = h / 6 * (area(_L_OUTLINE) + 4 * area(mid) + area(top))   # exact for planar ruled sides
    assert solid.volume == pytest.approx(expected, rel=1e-12), (solid.volume, expected)


def test_documented_recipe_rejects_inputs_outside_its_preconditions():
    import re

    text = get_build123d_documentation("rebuild_from_reference_step")
    block = next(b for b in re.findall(r"```python\n(.*?)```", text, re.S) if "def frustum" in b)
    scope = {}
    exec(block, scope)
    square = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
    assert len(scope["mitre_inset"](square, 1.9)) == 4
    with pytest.raises(ValueError, match="collapses"):
        scope["mitre_inset"](square, 3.0)          # collapsed at 2.0; the formula alone returns a 2x2 square
    with pytest.raises(ValueError, match="not parallel"):
        scope["frustum"](square, 0.0, [(0.0, 0.0), (4.0, 1.0), (4.0, 4.0), (0.0, 4.0)], 1.0)
