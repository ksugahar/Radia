# -*- coding: utf-8 -*-
r"""build123d strengthening marathon -- 100 parametric solids exercising build123d's full algebra-mode
surface (2D sketching, extrude/revolve/loft/sweep, primitives, booleans, local ops, patterns, features,
curves/surfaces, mechanical + EM/CAE archetypes, robustness).  Each problem BUILDS a real solid and
asserts a quantitative invariant (closed-form volume, face/edge topology, bounding box, or Netgen
meshability) -- a regression + documentation corpus that also surfaces build123d gaps.

Run: ``pytest tests/test_build123d_marathon.py`` (100 parametrized cases), or ``python ...`` for a tally.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from build123d import (Align, Box, Cylinder, Sphere, Cone, Torus, Wedge, Circle, Rectangle, Ellipse,
                       GeomType, RegularPolygon, RectangleRounded, Triangle, Polygon, Text, SlotOverall,
                       Plane, Pos, Rot, Axis, Compound, extrude, revolve, loft, sweep, fillet, chamfer,
                       mirror, scale, make_face, offset, Helix, Spline, Polyline, Line, CenterArc)

PROBLEMS = []


def problem(pid, title):
    def deco(fn):
        PROBLEMS.append((pid, title, fn))
        return fn
    return deco


def close(a, b, rel=1e-3):
    assert abs(a - b) <= rel * abs(b) + 1e-12, f"{a} != {b} (rel {rel})"
    return True


# =====================================================================================================
# Batch 1 (P001-P010) -- 2D sketches -> thin extrudes
# =====================================================================================================
@problem("P001", "rectangular plate")
def _p():
    w, h, t = 4.0, 3.0, 0.5
    part = extrude(Rectangle(w, h), t)
    close(part.volume, w * h * t)
    return part


@problem("P002", "disk")
def _p():
    r, t = 2.0, 0.5
    part = extrude(Circle(r), t)
    close(part.volume, math.pi * r ** 2 * t)
    return part


@problem("P003", "elliptical plate")
def _p():
    a, b, t = 3.0, 2.0, 0.4
    part = extrude(Ellipse(a, b), t)
    close(part.volume, math.pi * a * b * t)
    return part


@problem("P004", "hexagonal prism")
def _p():
    r, t, n = 2.0, 1.0, 6
    part = extrude(RegularPolygon(r, n), t)
    area = 0.5 * n * r ** 2 * math.sin(2 * math.pi / n)
    close(part.volume, area * t)
    return part


@problem("P005", "rounded rectangle plate")
def _p():
    w, h, rf, t = 4.0, 3.0, 0.5, 1.0
    part = extrude(RectangleRounded(w, h, rf), t)
    close(part.volume, (w * h - (4 - math.pi) * rf ** 2) * t)
    return part


@problem("P006", "right-triangle prism")
def _p():
    a, b, t = 3.0, 4.0, 1.0
    part = extrude(Triangle(a=a, b=b, C=90), t)
    close(part.volume, 0.5 * a * b * t)
    return part


@problem("P007", "obround slot")
def _p():
    L, w, t = 6.0, 2.0, 1.0
    part = extrude(SlotOverall(L, w), t)
    area = (L - w) * w + math.pi * (w / 2) ** 2
    close(part.volume, area * t)
    return part


@problem("P008", "annular plate (washer)")
def _p():
    R, r, t = 3.0, 1.5, 0.5
    part = extrude(Circle(R) - Circle(r), t)
    close(part.volume, math.pi * (R ** 2 - r ** 2) * t)
    return part


@problem("P009", "L-shaped region")
def _p():
    t = 1.0
    face = make_face(Polyline((0, 0), (3, 0), (3, 1), (1, 1), (1, 2), (0, 2), (0, 0)))
    part = extrude(face, t)
    close(part.volume, (3 * 1 + 1 * 1) * t)         # 3x1 base + 1x1 upright = 4
    return part


@problem("P010", "2D union of overlapping disks")
def _p():
    r, d, t = 2.0, 2.0, 0.5
    part = extrude(Circle(r) + Pos(d, 0) * Circle(r), t)
    assert 0 < part.volume < 2 * math.pi * r ** 2 * t, "overlap removes some area"
    return part


# =====================================================================================================
# Batch 2 (P011-P020) -- extrude / revolve / loft / sweep family
# =====================================================================================================
@problem("P011", "two-sided extrude")
def _p():
    r, t = 1.5, 0.7
    part = extrude(Circle(r), t, both=True)
    close(part.volume, math.pi * r ** 2 * 2 * t)
    return part


@problem("P012", "tapered (draft) extrude -> frustum")
def _p():
    r, h, taper = 2.0, 3.0, 10.0
    part = extrude(Circle(r), h, taper=taper)
    r_top = r - h * math.tan(math.radians(taper))
    close(part.volume, math.pi * h / 3.0 * (r ** 2 + r * r_top + r_top ** 2), rel=5e-3)
    return part


@problem("P013", "full revolve -> torus")
def _p():
    R, r = 4.0, 1.0
    part = revolve(Pos(R, 0) * Circle(r), Axis.Y, 360)
    close(part.volume, 2 * math.pi ** 2 * R * r ** 2, rel=3e-3)
    return part


@problem("P014", "partial revolve (90 deg)")
def _p():
    R, r = 4.0, 1.0
    part = revolve(Pos(R, 0) * Circle(r), Axis.Y, 90)
    close(part.volume, 0.25 * 2 * math.pi ** 2 * R * r ** 2, rel=3e-3)
    return part


@problem("P015", "lofted frustum")
def _p():
    R, r, h = 2.0, 1.0, 4.0
    part = loft([Circle(R), Pos(0, 0, h) * Circle(r)])
    close(part.volume, math.pi * h / 3.0 * (R ** 2 + R * r + r ** 2), rel=5e-3)
    return part


@problem("P016", "square-to-circle transition loft")
def _p():
    a, r, h = 3.0, 1.0, 3.0
    part = loft([Rectangle(a, a), Pos(0, 0, h) * Circle(r)])
    assert part.is_valid and part.volume > 0
    return part


@problem("P017", "pipe bend (sweep circle along arc)")
def _p():
    r, R = 0.4, 5.0
    path = CenterArc((0, 0, 0), R, 0, 90)
    part = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(r), path=path)
    close(part.volume, math.pi * r ** 2 * (2 * math.pi * R / 4.0), rel=5e-3)
    return part


@problem("P018", "helical coil (sweep circle along helix)")
def _p():
    r, R, pitch, h = 0.3, 2.0, 4.0, 8.0
    hel = Helix(pitch=pitch, height=h, radius=R)
    part = sweep(Plane(origin=hel @ 0.0, z_dir=hel % 0.0) * Circle(r), path=hel)
    hlen = math.sqrt(h ** 2 + (h / pitch * 2 * math.pi * R) ** 2)
    close(part.volume, math.pi * r ** 2 * hlen, rel=5e-3)
    return part


@problem("P019", "sweep along a 3D spline")
def _p():
    r = 0.3
    path = Spline((0, 0, 0), (2, 1, 1), (4, -1, 2), (6, 0, 4))
    part = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(r), path=path)
    close(part.volume, math.pi * r ** 2 * path.length, rel=2e-2)
    return part


@problem("P020", "revolved profile (vase)")
def _p():
    prof = make_face(Polyline((1, 0), (1.5, 0), (1.2, 2), (1.6, 4), (0.0, 4), (0, 0), (1, 0)))
    part = revolve(prof, Axis.Y, 360)
    assert part.is_valid and part.volume > 0
    return part


# =====================================================================================================
# Batch 3 (P021-P030) -- 3D primitives & rigid transforms
# =====================================================================================================
@problem("P021", "box")
def _p():
    part = Box(2, 3, 4)
    close(part.volume, 24.0)
    return part


@problem("P022", "cylinder")
def _p():
    R, h = 1.5, 5.0
    part = Cylinder(R, h)
    close(part.volume, math.pi * R ** 2 * h)
    return part


@problem("P023", "sphere")
def _p():
    R = 2.0
    part = Sphere(R)
    close(part.volume, 4.0 / 3.0 * math.pi * R ** 3, rel=2e-3)
    return part


@problem("P024", "conical frustum")
def _p():
    rb, rt, h = 2.0, 0.5, 4.0
    part = Cone(rb, rt, h)
    close(part.volume, math.pi * h / 3.0 * (rb ** 2 + rb * rt + rt ** 2))
    return part


@problem("P025", "torus")
def _p():
    R, r = 10.0, 2.0
    part = Torus(R, r)
    close(part.volume, 2 * math.pi ** 2 * R * r ** 2, rel=3e-3)
    return part


@problem("P026", "wedge")
def _p():
    part = Wedge(4, 3, 2, 1, 1, 3, 2)
    assert part.is_valid and part.volume > 0
    return part


@problem("P027", "translate preserves volume, moves centroid")
def _p():
    base = Box(2, 2, 2)
    part = Pos(5, 0, 0) * base
    close(part.volume, base.volume)
    close(part.center().X, 5.0, rel=1e-6)
    return part


@problem("P028", "rotate preserves volume")
def _p():
    base = Box(2, 3, 4)
    part = Rot(0, 0, 37) * base
    close(part.volume, base.volume)
    return part


@problem("P029", "mirror preserves volume, flips centroid")
def _p():
    base = Pos(3, 0, 0) * Box(1, 1, 1)
    part = mirror(base, Plane.YZ)
    close(part.volume, base.volume)
    close(part.center().X, -3.0, rel=1e-6)
    return part


@problem("P030", "uniform scale cubes the volume")
def _p():
    base = Box(1, 2, 3)
    part = scale(base, 2.0)
    close(part.volume, base.volume * 8.0)
    return part


# =====================================================================================================
# Batch 4 (P031-P040) -- boolean & local ops
# =====================================================================================================
@problem("P031", "union of disjoint boxes")
def _p():
    a, b = Box(2, 2, 2), Pos(5, 0, 0) * Box(2, 2, 2)
    part = a + b
    close(part.volume, 16.0)
    return part


@problem("P032", "drill a through hole")
def _p():
    w, R = 4.0, 1.0
    part = Box(w, w, w) - Cylinder(R, w + 1)
    close(part.volume, w ** 3 - math.pi * R ** 2 * w)
    return part


@problem("P033", "intersect box and sphere")
def _p():
    box, sph = Box(3, 3, 3), Sphere(2)
    part = box & sph
    assert 0 < part.volume < min(box.volume, sph.volume)
    return part


@problem("P034", "fillet all edges of a cube")
def _p():
    box = Box(3, 3, 3)
    part = fillet(box.edges(), 0.5)
    assert part.is_valid and part.volume < box.volume
    return part


@problem("P035", "fillet only the vertical edges")
def _p():
    box = Box(3, 3, 3)
    part = fillet(box.edges().filter_by(Axis.Z), 0.6)
    assert part.is_valid and part.volume < box.volume
    return part


@problem("P036", "chamfer all edges of a cube")
def _p():
    box = Box(3, 3, 3)
    part = chamfer(box.edges(), 0.5)
    assert part.is_valid and part.volume < box.volume
    return part


@problem("P037", "closed hollow shell")
def _p():
    w, t = 4.0, 0.3
    part = Box(w, w, w) - offset(Box(w, w, w), amount=-t)
    close(part.volume, w ** 3 - (w - 2 * t) ** 3)
    return part


@problem("P038", "thin-walled tube")
def _p():
    from radia_mcp.build123d.modeling import tube
    R, t, h = 2.0, 0.3, 4.0
    part = tube(R - t, R, h)
    close(part.volume, math.pi * (R ** 2 - (R - t) ** 2) * h)
    return part


@problem("P039", "keep upper half by boolean cut")
def _p():
    part = Box(4, 4, 4) - Pos(0, 0, -2) * Box(8, 8, 4)
    close(part.volume, 32.0)
    assert part.center().Z > 0
    return part


@problem("P040", "tee junction of two cylinders")
def _p():
    a = Cylinder(0.6, 6)
    b = Rot(90, 0, 0) * Cylinder(0.6, 6)
    part = a + b
    assert part.is_valid and part.volume > a.volume
    return part


# =====================================================================================================
# Batch 5 (P041-P050) -- patterns & arrays
# =====================================================================================================
@problem("P041", "linear array of pins")
def _p():
    from radia_mcp.build123d.modeling import assembly
    n, dx = 6, 1.5
    unit = Cylinder(0.3, 2)
    part = assembly(*[Pos(i * dx, 0, 0) * unit for i in range(n)])
    assert len(part.solids()) == n
    close(part.volume, n * unit.volume)
    return part


@problem("P042", "2D grid array")
def _p():
    from radia_mcp.build123d.modeling import assembly
    nx, ny, p = 3, 4, 1.2
    unit = Box(0.4, 0.4, 0.4)
    part = assembly(*[Pos(i * p, j * p, 0) * unit for i in range(nx) for j in range(ny)])
    assert len(part.solids()) == nx * ny
    close(part.volume, nx * ny * unit.volume)
    return part


@problem("P043", "polar array (gear-blank spokes)")
def _p():
    from radia_mcp.build123d.modeling import polar_array
    n = 8
    spoke = Pos(2, 0, 0) * Box(2, 0.3, 0.5)
    part = polar_array(spoke, n, 360.0)
    assert len(part.solids()) == n
    close(part.volume, n * spoke.volume)
    return part


@problem("P044", "mirror pattern (L+R bracket)")
def _p():
    from radia_mcp.build123d.modeling import assembly
    left = Pos(2, 0, 0) * Box(1, 2, 1)
    part = assembly(left, mirror(left, Plane.YZ))
    assert len(part.solids()) == 2
    close(part.center().X, 0.0, rel=1e-6)
    return part


@problem("P045", "hexagonal pin grid")
def _p():
    from radia_mcp.build123d.modeling import assembly
    p, rows = 1.0, 3
    pts = []
    for j in range(rows):
        for i in range(rows):
            pts.append((i * p + (j % 2) * p / 2, j * p * math.sqrt(3) / 2))
    unit = Cylinder(0.2, 1)
    part = assembly(*[Pos(x, y, 0) * unit for x, y in pts])
    assert len(part.solids()) == rows * rows
    return part


@problem("P046", "cylindrical (radial x angular) array")
def _p():
    from radia_mcp.build123d.modeling import assembly
    na, nr = 6, 2
    unit = Box(0.3, 0.3, 0.3)
    parts = []
    for ir in range(nr):
        for ia in range(na):
            parts.append(Rot(0, 0, ia * 360 / na) * Pos(2 + ir, 0, 0) * unit)
    part = assembly(*parts)
    assert len(part.solids()) == na * nr
    return part


@problem("P047", "array along a spline path")
def _p():
    from radia_mcp.build123d.modeling import assembly
    path = Spline((0, 0, 0), (3, 1, 0), (6, -1, 0), (9, 0, 0))
    n = 7
    unit = Sphere(0.25)
    part = assembly(*[Pos(*tuple(path @ (k / (n - 1)))) * unit for k in range(n)])
    assert len(part.solids()) == n
    return part


@problem("P048", "ring of holes in a plate")
def _p():
    from radia_mcp.build123d.modeling import polar_array
    plate = extrude(Circle(4), 0.5)
    hole = Pos(2.5, 0, 0) * Cylinder(0.3, 2)
    holes = polar_array(hole, 8, 360.0)
    part = plate - holes
    assert part.is_valid and part.volume < plate.volume
    return part


@problem("P049", "offset brick rows (running bond)")
def _p():
    from radia_mcp.build123d.modeling import assembly
    brick = Box(2.0, 1.0, 0.8)
    bricks = []
    for row in range(3):
        off = (row % 2) * 1.0
        for col in range(3):
            bricks.append(Pos(col * 2.1 + off, 0, row * 0.9) * brick)
    part = assembly(*bricks)
    assert len(part.solids()) == 9
    return part


@problem("P050", "phyllotaxis (golden-angle) disk of pins")
def _p():
    from radia_mcp.build123d.modeling import assembly
    n = 24
    ga = math.radians(137.508)
    unit = Cylinder(0.12, 1)
    pts = [(0.5 * math.sqrt(k) * math.cos(k * ga), 0.5 * math.sqrt(k) * math.sin(k * ga)) for k in range(1, n + 1)]
    part = assembly(*[Pos(x, y, 0) * unit for x, y in pts])
    assert len(part.solids()) == n
    return part


# =====================================================================================================
# Batch 6 (P051-P060) -- features
# =====================================================================================================
@problem("P051", "counterbored hole")
def _p():
    w = 4.0
    part = Box(w, w, 2) - Cylinder(0.5, 3) - Pos(0, 0, 0.5) * Cylinder(1.0, 1.2)
    assert part.is_valid and part.volume < w * w * 2
    return part


@problem("P052", "countersunk hole")
def _p():
    w = 4.0
    part = Box(w, w, 2) - Cylinder(0.5, 3) - Pos(0, 0, 1.0) * Cone(1.0, 0.5, 0.6)
    assert part.is_valid
    return part


@problem("P053", "boss on a plate")
def _p():
    plate = extrude(Rectangle(5, 5), 0.5)
    boss = Pos(0, 0, 0.5) * Cylinder(1.0, 1.5, align=(Align.CENTER, Align.CENTER, Align.MIN))
    part = plate + boss
    close(part.volume, plate.volume + boss.volume, rel=5e-3)
    return part


@problem("P054", "stiffening rib")
def _p():
    base = extrude(Rectangle(6, 4), 0.4)
    rib = Pos(0, 0, 0.4) * Box(0.4, 4, 1.2, align=(Align.CENTER, Align.CENTER, Align.MIN))
    part = base + rib
    assert part.is_valid and part.volume > base.volume
    return part


@problem("P055", "revolved groove in a shaft")
def _p():
    shaft = Cylinder(1.0, 6)
    groove = Torus(1.0, 0.2)
    part = shaft - groove
    assert part.is_valid and part.volume < shaft.volume
    return part


@problem("P056", "keyway slot in a shaft")
def _p():
    shaft = Rot(0, 90, 0) * Cylinder(1.0, 6)
    key = Box(6, 0.4, 0.4)
    part = shaft - key
    assert part.is_valid and part.volume < shaft.volume
    return part


@problem("P057", "dovetail prism")
def _p():
    prof = Polygon((-1.5, 0), (1.5, 0), (1.0, 1.5), (-1.0, 1.5))
    part = extrude(prof, 3)
    assert part.is_valid and part.volume > 0
    return part


@problem("P058", "triangular gusset web")
def _p():
    web = extrude(Triangle(a=2, b=2, C=90), 0.3)
    assert web.is_valid and web.volume > 0
    return web


@problem("P059", "embossed text on a plate")
def _p():
    plate = extrude(Rectangle(8, 3), 0.5)
    letters = Pos(0, 0, 0.5) * extrude(Text("RADIA", 2), 0.3)
    part = plate + letters
    assert part.is_valid and part.volume > plate.volume
    return part


@problem("P060", "chamfered hole rim")
def _p():
    plate = Box(4, 4, 1)
    drilled = plate - Cylinder(0.8, 2)
    rim = drilled.edges().filter_by(GeomType.CIRCLE).sort_by(Axis.Z)[-1]
    part = chamfer(rim, 0.2)
    assert part.is_valid and part.volume < drilled.volume
    return part


# =====================================================================================================
# Batch 7 (P061-P070) -- curves & surfaces
# =====================================================================================================
@problem("P061", "teardrop (spline) region")
def _p():
    pts = [(0, 0)] + [(2 + 1.5 * math.cos(a), 1.5 * math.sin(a))
                      for a in [math.radians(d) for d in range(-150, 151, 20)]] + [(0, 0)]
    part = extrude(make_face(Polyline(*pts)), 0.5)
    assert part.is_valid and part.volume > 0
    return part


@problem("P062", "conical (tapered) spring")
def _p():
    N, h, r0, r1, a = 4, 8.0, 0.6, 2.0, 0.15
    m = 60 * N
    pts = [((r0 + (r1 - r0) * s / m) * math.cos(2 * math.pi * N * s / m),
            (r0 + (r1 - r0) * s / m) * math.sin(2 * math.pi * N * s / m), h * s / m) for s in range(m + 1)]
    path = Spline(*pts)
    part = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(a), path=path)
    assert part.is_valid and part.volume > 0
    return part


@problem("P063", "Archimedean flat spiral spring")
def _p():
    turns, b, a = 3, 0.25, 0.08
    m = 80 * turns
    pts = [(b * (2 * math.pi * turns * s / m) * math.cos(2 * math.pi * turns * s / m),
            b * (2 * math.pi * turns * s / m) * math.sin(2 * math.pi * turns * s / m), 0.0) for s in range(1, m + 1)]
    path = Spline(*pts)
    part = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(a), path=path)
    assert part.is_valid and part.volume > 0
    return part


@problem("P064", "involute gear tooth")
def _p():
    rb, top, alpha = 3.0, 1.0, math.radians(8)        # rotate flanks apart -> finite-width root
    tvals = [0.025 * k for k in range(45)]
    inv = [(rb * (math.cos(t) + t * math.sin(t)), rb * (math.sin(t) - t * math.cos(t))) for t in tvals]
    inv = [p for p in inv if math.hypot(*p) <= rb + top]
    def rot(p, ang):
        x, y = p
        return (x * math.cos(ang) - y * math.sin(ang), x * math.sin(ang) + y * math.cos(ang))
    left = [rot(p, alpha) for p in inv]                # root -> tip
    right = [rot((x, -y), -alpha) for x, y in inv]     # mirror, root -> tip
    pts = left + list(reversed(right)) + [left[0]]     # up left, down right, close along root base
    part = extrude(make_face(Polyline(*pts)), 1.0)
    assert part.is_valid and part.volume > 0
    return part


@problem("P065", "harmonic cam profile")
def _p():
    base, lift = 2.0, 0.8
    pts = [((base + lift * (1 - math.cos(math.radians(d))) / 2) * math.cos(math.radians(d)),
            (base + lift * (1 - math.cos(math.radians(d))) / 2) * math.sin(math.radians(d))) for d in range(0, 360, 5)]
    pts.append(pts[0])
    part = extrude(make_face(Polyline(*pts)), 0.6)
    assert part.is_valid and part.volume > 0
    return part


@problem("P066", "NACA 0012 airfoil section")
def _p():
    c, tk = 6.0, 0.12
    xs = [c * (0.5 - 0.5 * math.cos(math.pi * i / 40)) for i in range(41)]      # cosine spacing
    def yt(x):
        xn = x / c
        return 5 * tk * c * (0.2969 * math.sqrt(xn) - 0.1260 * xn - 0.3516 * xn ** 2
                             + 0.2843 * xn ** 3 - 0.1036 * xn ** 4)            # closed (sharp) TE coeff
    up = [(x, yt(x)) for x in xs]                       # LE(0,0) -> TE(c,0)
    lo = [(x, -yt(x)) for x in reversed(xs)]            # TE(c,0) -> LE(0,0)
    pts = up + lo[1:-1] + [up[0]]                       # drop the duplicate TE / LE points
    part = extrude(make_face(Polyline(*pts)), 1.5)
    assert part.is_valid and part.volume > 0
    return part


@problem("P067", "superellipse (squircle) prism")
def _p():
    a, b, n = 3.0, 2.0, 4.0
    def sgnpow(v, e):
        return (1 if v >= 0 else -1) * abs(v) ** e
    pts = [(a * sgnpow(math.cos(math.radians(d)), 2 / n), b * sgnpow(math.sin(math.radians(d)), 2 / n))
           for d in range(0, 360, 6)]
    pts.append(pts[0])
    part = extrude(make_face(Polyline(*pts)), 1.0)
    assert part.is_valid and part.volume > 0
    return part


@problem("P068", "twisted (lofted) prism")
def _p():
    part = loft([Rectangle(2, 2), Pos(0, 0, 4) * Rot(0, 0, 45) * Rectangle(2, 2)])
    assert part.is_valid and part.volume > 0
    return part


@problem("P069", "multi-section barrel loft")
def _p():
    part = loft([Circle(1.0), Pos(0, 0, 2) * Circle(1.6), Pos(0, 0, 4) * Circle(0.8)])
    assert part.is_valid and part.volume > 0
    return part


@problem("P070", "toroidal solenoid winding path")
def _p():
    R, a, rw, N, frac = 8.0, 1.2, 0.08, 8, 0.75       # open partial winding -> no closed-loop self-touch
    m = 40 * N
    pts = [(((R + a * math.cos(2 * math.pi * N * frac * s / m)) * math.cos(2 * math.pi * frac * s / m)),
            ((R + a * math.cos(2 * math.pi * N * frac * s / m)) * math.sin(2 * math.pi * frac * s / m)),
            a * math.sin(2 * math.pi * N * frac * s / m)) for s in range(m + 1)]
    path = Spline(*pts)
    part = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(rw), path=path)
    assert part.is_valid and part.volume > 0
    return part


# =====================================================================================================
# Batch 8 (P071-P080) -- mechanical archetypes
# =====================================================================================================
@problem("P071", "spur gear (toothed blank)")
def _p():
    from radia_mcp.build123d.modeling import polar_array
    nteeth = 12
    hub = Cylinder(3.0, 1.0)
    tooth = Pos(3.1, 0, 0) * Box(0.8, 0.5, 1.0)
    part = hub + polar_array(tooth, nteeth, 360.0)
    assert part.is_valid and part.volume > hub.volume
    return part


@problem("P072", "gear rack")
def _p():
    from radia_mcp.build123d.modeling import assembly
    bar = Box(12, 1.0, 1.0)
    teeth = [Pos(-5 + i * 1.2, 0.6, 0) * Rot(0, 0, 45) * Box(0.5, 0.5, 1.0) for i in range(9)]
    part = assembly(bar, *teeth)
    assert len(part.solids()) == 10
    return part


@problem("P073", "sprocket (toothed disk with gaps)")
def _p():
    from radia_mcp.build123d.modeling import polar_array
    disk = Cylinder(3.0, 0.8)
    gap = Pos(3.0, 0, 0) * Cylinder(0.4, 1.0)
    part = disk - polar_array(gap, 16, 360.0)
    assert part.is_valid and part.volume < disk.volume
    return part


@problem("P074", "V-belt pulley")
def _p():
    body = Cylinder(3.0, 2.0)
    groove = Torus(3.0, 0.6)                       # V approximated by a round groove ring
    part = body - groove - Cylinder(0.5, 3)        # plus a bore
    assert part.is_valid and part.volume < body.volume
    return part


@problem("P075", "timing pulley (axial teeth)")
def _p():
    from radia_mcp.build123d.modeling import polar_array
    body = Cylinder(2.5, 2.0)
    slot = Pos(2.5, 0, 0) * Box(0.4, 0.4, 2.2)
    part = body - polar_array(slot, 20, 360.0)
    assert part.is_valid and part.volume < body.volume
    return part


@problem("P076", "bearing race (grooved ring)")
def _p():
    from radia_mcp.build123d.modeling import tube
    ring = tube(2.0, 3.0, 1.5)
    groove = Torus(3.0, 0.35)
    part = ring - groove
    assert part.is_valid and part.volume < ring.volume
    return part


@problem("P077", "threaded rod (helical thread)")
def _p():
    core = Cylinder(1.0, 6, align=(Align.CENTER, Align.CENTER, Align.MIN))
    hel = Helix(pitch=1.0, height=6, radius=1.0)
    thread = sweep(Plane(origin=hel @ 0.0, z_dir=hel % 0.0) * Triangle(a=0.4, b=0.4, C=90), path=hel)
    part = core + thread
    assert part.is_valid and part.volume > core.volume
    return part


@problem("P078", "hex nut")
def _p():
    body = extrude(RegularPolygon(2.0, 6), 1.5)
    part = body - Cylinder(1.0, 3)
    part = chamfer(part.edges().filter_by(GeomType.CIRCLE), 0.2)
    assert part.is_valid and part.volume < body.volume
    return part


@problem("P079", "chamfered flat washer")
def _p():
    from radia_mcp.build123d.modeling import tube
    w = tube(1.2, 2.5, 0.4)
    part = chamfer(w.edges().filter_by(GeomType.CIRCLE).sort_by(Axis.Z)[-2:], 0.1)
    assert part.is_valid and part.volume < w.volume
    return part


@problem("P080", "coil spring")
def _p():
    rw, R, pitch, h = 0.25, 1.5, 1.2, 6.0
    hel = Helix(pitch=pitch, height=h, radius=R)
    part = sweep(Plane(origin=hel @ 0.0, z_dir=hel % 0.0) * Circle(rw), path=hel)
    hlen = math.sqrt(h ** 2 + (h / pitch * 2 * math.pi * R) ** 2)
    close(part.volume, math.pi * rw ** 2 * hlen, rel=5e-3)
    return part


# =====================================================================================================
# Batch 9 (P081-P090) -- EM / CAE archetypes
# =====================================================================================================
@problem("P081", "flanged coil bobbin (revolved spool)")
def _p():
    prof = make_face(Polyline((0.6, 0), (2, 0), (2, 0.4), (1.0, 0.4), (1.0, 3.6),
                              (2, 3.6), (2, 4), (0.6, 4), (0.6, 0)))
    part = revolve(prof, Axis.Y, 360)
    assert part.is_valid and part.volume > 0
    return part


@problem("P082", "rectangular-section toroid core")
def _p():
    R, w, h = 5.0, 1.0, 1.4
    part = revolve(Pos(R, 0) * Rectangle(w, h), Axis.Y, 360)
    close(part.volume, w * h * 2 * math.pi * R, rel=3e-3)        # Pappus
    return part


@problem("P083", "E-core")
def _p():
    yoke = Box(6, 1, 2)
    legs = [Pos(x, 1.0, 0) * Box(1, 1, 2, align=(Align.CENTER, Align.MIN, Align.CENTER)) for x in (-2.5, 0, 2.5)]
    part = yoke + legs[0] + legs[1] + legs[2]
    assert part.is_valid and part.volume > yoke.volume
    return part


@problem("P084", "U-core")
def _p():
    base = Box(5, 1, 2)
    legs = [Pos(x, 1.0, 0) * Box(1, 2.5, 2, align=(Align.CENTER, Align.MIN, Align.CENTER)) for x in (-2, 2)]
    part = base + legs[0] + legs[1]
    assert part.is_valid and part.volume > base.volume
    return part


@problem("P085", "tapered pole shoe")
def _p():
    prof = Polygon((-0.6, 0), (0.6, 0), (1.4, -1.0), (-1.4, -1.0))     # widens toward the gap
    part = extrude(prof, 2.0)
    assert part.is_valid and part.volume > 0
    return part


@problem("P086", "slot wedge (trapezoid)")
def _p():
    prof = Polygon((-1.0, 0), (1.0, 0), (0.6, 0.5), (-0.6, 0.5))
    part = extrude(prof, 3.0)
    assert part.is_valid and part.volume > 0
    return part


@problem("P087", "busbar with bends")
def _p():
    path = Spline((0, 0, 0), (4, 0, 0), (5, 2, 0), (5, 5, 1), (5, 5, 4))
    sec = Plane(origin=path @ 0.0, z_dir=path % 0.0) * Rectangle(0.8, 0.3)
    part = sweep(sec, path=path)
    assert part.is_valid and part.volume > 0
    return part


@problem("P088", "finned heat sink")
def _p():
    base = Box(8, 6, 0.5)
    fins = [Pos(-3 + i * 1.0, 0, 0.5) * Box(0.3, 6, 2.0, align=(Align.CENTER, Align.CENTER, Align.MIN))
            for i in range(7)]
    part = base
    for f in fins:
        part = part + f
    assert part.is_valid and part.volume > base.volume
    return part


@problem("P089", "rectangular waveguide section")
def _p():
    a, b, wall, L = 2.0, 1.0, 0.2, 5.0
    part = Box(a + 2 * wall, b + 2 * wall, L) - Box(a, b, L + 1)
    close(part.volume, ((a + 2 * wall) * (b + 2 * wall) - a * b) * L)
    return part


@problem("P090", "pyramidal horn")
def _p():
    outer = loft([Rectangle(2, 1), Pos(0, 0, 5) * Rectangle(6, 4)])
    inner = loft([Rectangle(1.6, 0.7), Pos(0, 0, 5) * Rectangle(5.6, 3.7)])
    part = outer - inner
    assert part.is_valid and 0 < part.volume < outer.volume
    return part


# =====================================================================================================
# Batch 10 (P091-P100) -- assemblies, sheet metal, robustness / CAE gates
# =====================================================================================================
@problem("P091", "labelled motor assembly")
def _p():
    from radia_mcp.build123d.modeling import assembly, tube
    shaft = Cylinder(0.5, 6); shaft.label = "shaft"
    rotor = Cylinder(2.0, 4); rotor.label = "rotor"
    stator = tube(3.0, 4.0, 4); stator.label = "stator"
    part = assembly(shaft, rotor, stator, label="motor")
    assert len(part.solids()) == 3 and {c.label for c in part.children} == {"shaft", "rotor", "stator"}
    return part


@problem("P092", "exploded assembly offsets")
def _p():
    from radia_mcp.build123d.modeling import assembly
    parts = [Pos(0, 0, 4 * i) * Box(2, 2, 1) for i in range(3)]
    part = assembly(*parts)
    zs = sorted(c.center().Z for c in part.children)
    assert len(part.solids()) == 3 and zs[-1] - zs[0] > 6
    return part


@problem("P093", "sheet-metal L-bracket")
def _p():
    t = 0.2
    web = Box(4, 3, t)
    flange = Pos(2, 0, 0) * Box(t, 3, 2, align=(Align.MIN, Align.CENTER, Align.MIN))
    part = web + flange
    assert part.is_valid and part.volume > 0
    return part


@problem("P094", "sheet-metal flanged tab")
def _p():
    t = 0.15
    base = Box(5, 3, t)
    tab = Pos(-2.5, 0, 0) * Box(t, 1.5, 1.5, align=(Align.MAX, Align.CENTER, Align.MIN))
    part = base + tab
    assert part.is_valid and part.volume > base.volume
    return part


@problem("P095", "STEP export round-trip")
def _p():
    import tempfile
    from build123d import export_step, import_step
    src = fillet(Box(3, 2, 1).edges(), 0.2)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "rt.step")
        export_step(src, f)
        back = import_step(f)
    close(back.volume, src.volume, rel=1e-6)
    return back


@problem("P096", "STL export")
def _p():
    import tempfile
    from build123d import export_stl
    part = Sphere(1.5)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "s.stl")
        ok = export_stl(part, f)
        assert os.path.isfile(f) and os.path.getsize(f) > 0
    return part


@problem("P097", "Netgen tet-mesh gate")
def _p():
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    part = fillet(Box(4, 3, 2).edges(), 0.4)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "m.step")
        export_step(part, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=1.0))
    assert mesh.ne > 50, f"should tet-mesh (got {mesh.ne})"
    return part


@problem("P098", "mass-properties consistency")
def _p():
    part = Box(2, 3, 4)
    close(part.volume, 24.0)
    bb = part.bounding_box()
    close(bb.size.X, 2.0); close(bb.size.Y, 3.0); close(bb.size.Z, 4.0)
    assert abs(part.center().X) < 1e-6 and abs(part.center().Y) < 1e-6
    return part


@problem("P099", "section (cut face area)")
def _p():
    w, h, d = 2.0, 3.0, 4.0
    part = Box(w, d, h) - Pos(0, 0, -h / 2) * Box(w + 1, d + 1, h)     # keep z>0 half
    bottom = part.faces().sort_by(Axis.Z)[0]
    close(bottom.area, w * d)                                          # the new cut face
    return part


@problem("P100", "interference check")
def _p():
    a = Box(2, 2, 2)
    overlap = (a & (Pos(1, 0, 0) * Box(2, 2, 2))).volume
    assert overlap > 0, "overlapping solids interfere"
    try:
        clear = (a & (Pos(5, 0, 0) * Box(2, 2, 2))).volume
    except Exception:
        clear = 0.0
    assert clear < 1e-9, "disjoint solids do not interfere"
    return a


# -----------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("pid,title,fn", PROBLEMS, ids=[p[0] for p in PROBLEMS])
def test_marathon(pid, title, fn):
    part = fn()
    assert part is not None and getattr(part, "is_valid", True), f"{pid} {title}: invalid"


def main():
    ok = 0
    for pid, title, fn in PROBLEMS:
        try:
            fn()
            ok += 1
        except Exception as e:
            print(f"  {pid} FAIL {title}: {type(e).__name__}: {e}")
    print(f"[marathon] {ok}/{len(PROBLEMS)} problems verified")


if __name__ == "__main__":
    main()
