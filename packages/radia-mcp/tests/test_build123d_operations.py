# -*- coding: utf-8 -*-
r"""Tests for the generic solid-modelling operations promoted into radia_mcp.build123d.modeling
(sweep / revolve / loft / coil / strut / thicken / draft / shell / fillet / chamfer / grid & path array).
Each op gets a closed-form or topological invariant -- these are the reusable modeller verbs the marathon
corpus prototyped, now formal API."""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from build123d import Axis, Box, Circle, Rectangle, Pos, Line, Spline, GeomType
from radia_mcp.build123d.modeling import (swept, revolved, lofted, coil, strut, thicken, draft_extrude,
                                          shell, fillet_edges, chamfer_edges, grid_array, path_array)


def _close(a, b, rel=3e-3):
    assert abs(a - b) <= rel * abs(b) + 1e-12, f"{a} != {b}"


def test_swept_cylinder():
    s = swept(Circle(0.5), Line((0, 0, 0), (0, 0, 6)), label="pipe")
    _close(s.volume, math.pi * 0.5 ** 2 * 6)
    assert s.label == "pipe"


def test_revolved_torus():
    s = revolved(Pos(4, 0) * Circle(1.0), Axis.Y, 360, label="ring")
    _close(s.volume, 2 * math.pi ** 2 * 4 * 1.0 ** 2)
    assert s.label == "ring"


def test_lofted_frustum():
    s = lofted([Circle(2.0), Pos(0, 0, 4) * Circle(1.0)], label="cone")
    _close(s.volume, math.pi * 4 / 3 * (2 ** 2 + 2 * 1 + 1 ** 2), rel=5e-3)
    assert s.is_valid


def test_coil_helix():
    rw, R, pitch, h = 0.2, 2.0, 1.0, 5.0
    s = coil(Circle(rw), pitch, h, R, label="winding")
    hlen = math.sqrt(h ** 2 + (h / pitch * 2 * math.pi * R) ** 2)
    _close(s.volume, math.pi * rw ** 2 * hlen, rel=5e-3)


def test_strut_between_points():
    p0, p1, r = (0, 0, 0), (3, 4, 0), 0.25
    s = strut(p0, p1, r)
    _close(s.volume, math.pi * r ** 2 * 5.0)              # |p1-p0| = 5


def test_thicken_symmetric():
    s = thicken(Rectangle(4, 3), 0.5, label="plate")
    _close(s.volume, 4 * 3 * 0.5)
    bb = s.bounding_box()
    _close(bb.size.Z, 0.5)                                 # symmetric about z=0
    assert abs(s.center().Z) < 1e-9


def test_draft_extrude_frustum():
    r, h, taper = 2.0, 3.0, 10.0
    s = draft_extrude(Circle(r), h, taper)
    rt = r - h * math.tan(math.radians(taper))
    _close(s.volume, math.pi * h / 3 * (r ** 2 + r * rt + rt ** 2), rel=5e-3)


def test_shell_closed():
    w, t = 4.0, 0.3
    s = shell(Box(w, w, w), t)
    _close(s.volume, w ** 3 - (w - 2 * t) ** 3)


def test_shell_open_top():
    box = Box(4, 4, 4)
    top = box.faces().sort_by(Axis.Z)[-1]
    cup = shell(box, 0.3, openings=top)
    assert cup.is_valid and 0 < cup.volume < shell(Box(4, 4, 4), 0.3).volume + 1


def test_fillet_edges_all_and_subset():
    box = Box(3, 3, 3)
    f = fillet_edges(box, 0.4)
    assert f.is_valid and f.volume < box.volume
    fz = fillet_edges(box, 0.5, edge_filter=lambda e: e.filter_by(Axis.Z))
    assert fz.is_valid and fz.volume < box.volume


def test_chamfer_edges():
    box = Box(3, 3, 3)
    c = chamfer_edges(box, 0.4)
    assert c.is_valid and c.volume < box.volume


def test_grid_array():
    g = grid_array(Box(0.4, 0.4, 0.4), 3, 4, 1.0, 1.0, label="pin")
    assert len(g.solids()) == 12
    _close(g.volume, 12 * 0.4 ** 3)
    assert [c.label for c in g.children][0] == "pin_00_00"


def test_path_array():
    path = Spline((0, 0, 0), (3, 1, 0), (6, 0, 0))
    g = path_array(Box(0.3, 0.3, 0.3), path, 5, label="bead")
    assert len(g.solids()) == 5
    assert [c.label for c in g.children] == [f"bead_{k:02d}" for k in range(5)]


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok {name}")
    print("[operations] all promoted modelling verbs verified")


if __name__ == "__main__":
    main()
