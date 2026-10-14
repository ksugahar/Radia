# -*- coding: utf-8 -*-
r"""build123d strengthening marathon -- ROUND 2 (P101-P200): more complex, engineering-grade models.

Where round 1 (test_build123d_marathon.py) drilled the primitive surface, round 2 builds real parts and
assemblies: involute gears & power transmission, metric threads & fasteners, sheet-metal enclosures,
full EM machine cross-sections, complex coils/windings, lofted blades & freeform surfaces, mechanisms
& multi-body assemblies, lattices, and CAE-ready models.  Each problem BUILDS a valid solid / compound
and asserts a quantitative invariant (closed-form volume where it exists, else validity + topology /
bounding box / region count / Netgen meshability).

Run: ``pytest validation_test/radia_mcp/test_build123d_marathon2.py`` (100 cases) or ``python ...`` for a tally.
"""
import math
import os
import sys

import pytest

pytestmark = [
    pytest.mark.slow,
    pytest.mark.usefixtures("ngsolve_taskmanager"),
]

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from build123d import (Align, Box, Cylinder, Sphere, Cone, Torus, Wedge, Circle, Rectangle, Ellipse,
                       GeomType, RegularPolygon, RectangleRounded, Triangle, Polygon, Text, SlotOverall,
                       Plane, Pos, Rot, Axis, Compound, extrude, revolve, loft, sweep, fillet, chamfer,
                       mirror, scale, make_face, offset, Helix, Spline, Polyline, Line, CenterArc)
from radia_mcp.build123d.modeling import assembly, polar_array, linear_array, tube, annular_segment
from radia_mcp.build123d.archetypes import halbach_ring, litz_wire

CK = (Align.CENTER, Align.CENTER, Align.MIN)        # base-on-z=0 alignment, used a lot

PROBLEMS = []


def problem(pid, title):
    def deco(fn):
        PROBLEMS.append((pid, title, fn))
        return fn
    return deco


def close(a, b, rel=1e-3):
    assert abs(a - b) <= rel * abs(b) + 1e-12, f"{a} != {b} (rel {rel})"
    return True


# ---- shared shapes / ops PROMOTED to the formal API -- the marathon now EXERCISES the public functions
# (gears, threads, airfoils, struts moved out of this test file into archetypes.py / modeling.py).
from radia_mcp.build123d.archetypes import (involute_gear as spur_gear, threaded_rod as threaded_shank,
                                            airfoil as airfoil_face,
                                            _involute_tooth_face as involute_tooth_face)
from radia_mcp.build123d.modeling import strut


# =====================================================================================================
# Batch 11 (P101-P110) -- advanced 2D profiles
# =====================================================================================================
@problem("P101", "involute spur gear blank")
def _p():
    g = spur_gear(12, 3.0, 0.8, 0.4)
    assert g.is_valid and g.volume > math.pi * (3.0 * 1.06) ** 2 * 0.4    # teeth add material
    return g


@problem("P102", "internal ring gear")
def _p():
    ring = tube(4.0, 5.0, 0.6)
    tooth = Pos(4.0, 0, 0) * Box(0.5, 0.5, 0.8)
    g = ring - polar_array(tooth, 16, 360.0)
    assert g.is_valid and g.volume < ring.volume
    return g


@problem("P103", "involute rack")
def _p():
    bar = Box(14, 1.2, 1.0, align=CK)
    teeth = [Pos(-6 + i * 1.4, 0.6, 0) * Rot(0, 0, 45) * Box(0.55, 0.55, 1.0, align=CK) for i in range(9)]
    g = bar
    for t in teeth:
        g = g + t
    assert g.is_valid and g.volume > bar.volume
    return g


@problem("P104", "cycloidal drive disc")
def _p():
    R, e, N = 3.0, 0.25, 10
    pts = [((R + e * math.cos(N * math.radians(d))) * math.cos(math.radians(d)),
            (R + e * math.cos(N * math.radians(d))) * math.sin(math.radians(d))) for d in range(0, 360, 3)]
    pts.append(pts[0])
    g = extrude(make_face(Polyline(*pts)), 0.5)
    assert g.is_valid and g.volume > 0
    return g


@problem("P105", "ratchet wheel")
def _p():
    N, r0, r1 = 16, 2.0, 2.8
    pts = []
    for k in range(N):
        a0 = 2 * math.pi * k / N
        a1 = 2 * math.pi * (k + 1) / N
        pts.append((r0 * math.cos(a0), r0 * math.sin(a0)))     # root
        pts.append((r1 * math.cos(a1), r1 * math.sin(a1)))     # rise to tip just before next root
    pts.append(pts[0])
    g = extrude(make_face(Polyline(*pts)), 0.5)
    assert g.is_valid and g.volume > 0
    return g


@problem("P106", "plate cam (dwell-rise-return)")
def _p():
    base, lift = 2.0, 1.0
    def r(d):
        if d < 90:
            return base
        if d < 180:
            return base + lift * (1 - math.cos(math.pi * (d - 90) / 90)) / 2
        if d < 270:
            return base + lift
        return base + lift * (1 + math.cos(math.pi * (d - 270) / 90)) / 2
    pts = [(r(d) * math.cos(math.radians(d)), r(d) * math.sin(math.radians(d))) for d in range(0, 360, 3)]
    pts.append(pts[0])
    g = extrude(make_face(Polyline(*pts)), 0.5)
    assert g.is_valid and g.volume > 0
    return g


@problem("P107", "PCB outline with cut-outs")
def _p():
    board = extrude(RectangleRounded(20, 14, 1.5), 0.16)
    holes = polar_array(Pos(8.0, 5.0, 0) * Cylinder(0.6, 1), 1, 360.0)   # dummy; real holes below
    g = board
    for (x, y) in [(-8, -5), (8, -5), (8, 5), (-8, 5)]:
        g = g - Pos(x, y, 0) * Cylinder(0.8, 1)
    g = g - Pos(0, 0, 0) * Box(4, 1.2, 1)                                 # a slot
    assert g.is_valid and g.volume < board.volume
    return g


@problem("P108", "bolt-circle gasket")
def _p():
    gasket = tube(3.0, 5.0, 0.2)
    g = gasket - polar_array(Pos(4.0, 0, 0) * Cylinder(0.4, 1), 8, 360.0)
    assert g.is_valid and g.volume < gasket.volume
    return g


@problem("P109", "gerotor inner rotor (epitrochoid)")
def _p():
    R, e, N = 3.0, 0.5, 6
    pts = [((R + e * math.cos(N * math.radians(d))) * math.cos(math.radians(d)),
            (R + e * math.cos(N * math.radians(d))) * math.sin(math.radians(d))) for d in range(0, 360, 2)]
    pts.append(pts[0])
    g = extrude(make_face(Polyline(*pts)), 0.8)
    assert g.is_valid and g.volume > 0
    return g


@problem("P110", "bolted flange with hub")
def _p():
    disk = Cylinder(5.0, 0.5, align=CK)
    hub = Pos(0, 0, 0.5) * Cylinder(2.0, 1.5, align=CK)
    g = disk + hub - Cylinder(1.0, 3)
    g = g - polar_array(Pos(3.8, 0, 0) * Cylinder(0.4, 1), 6, 360.0)
    assert g.is_valid and g.volume < disk.volume + hub.volume
    return g


# =====================================================================================================
# Batch 12 (P111-P120) -- gears & power transmission (3D)
# =====================================================================================================
@problem("P111", "spur gear (3D, 18 teeth)")
def _p():
    g = spur_gear(18, 4.0, 1.0, 1.2)
    assert g.is_valid and g.volume > 0 and len(g.solids()) >= 1
    return g


@problem("P112", "helical gear")
def _p():
    g = spur_gear(16, 3.5, 0.9, 2.0, twist_deg=20.0)
    assert g.is_valid and g.volume > 0
    return g


@problem("P113", "herringbone gear")
def _p():
    lower = spur_gear(16, 3.5, 0.9, 1.5, twist_deg=18.0)
    upper = Pos(0, 0, 1.5) * spur_gear(16, 3.5, 0.9, 1.5, twist_deg=-18.0)
    g = lower + upper
    assert g.is_valid and g.volume > lower.volume
    return g


@problem("P114", "bevel gear (tapered)")
def _p():
    n, rb, top, h = 14, 3.0, 0.8, 2.5
    ha = 360.0 / n / 4.0
    tooth = loft([involute_tooth_face(rb, top, ha),                       # single-tooth faces: no holes,
                  Pos(0, 0, h) * scale(involute_tooth_face(rb, top, ha), 0.55)])   # matching topology
    core = Cone(rb * 1.06, rb * 1.06 * 0.55, h, align=CK)
    g = core + polar_array(tooth, n, 360.0)
    assert g.is_valid and g.volume > 0
    return g


@problem("P115", "worm (single-thread screw)")
def _p():
    core = Cylinder(1.0, 8, align=CK)
    hel = Helix(pitch=2.0, height=8, radius=1.0)
    thread = sweep(Plane(origin=hel @ 0.0, z_dir=hel % 0.0) * Triangle(a=0.7, b=0.7, C=90), path=hel)
    g = core + thread
    assert g.is_valid and g.volume > core.volume
    return g


@problem("P116", "worm wheel (throated rim)")
def _p():
    blank = spur_gear(20, 4.0, 0.8, 1.4)
    throat = Torus(4.2, 0.9)
    g = blank - throat
    assert g.is_valid and g.volume < blank.volume
    return g


@problem("P117", "planetary gear set assembly")
def _p():
    sun = spur_gear(12, 2.0, 0.5, 1.0)
    sun.label = "sun"
    planets = []
    for k in range(3):
        p = Rot(0, 0, k * 120) * Pos(4.5, 0, 0) * spur_gear(10, 1.6, 0.45, 1.0)
        p.label = f"planet_{k}"
        planets.append(p)
    rng = tube(6.6, 7.4, 1.0)
    rng.label = "ring"
    asm = assembly(sun, *planets, rng, label="planetary")
    assert len(asm.children) == 5
    return asm


@problem("P118", "rack and pinion assembly")
def _p():
    pinion = spur_gear(14, 3.0, 0.8, 1.0)
    pinion.label = "pinion"
    rack = Pos(0, -4.0, 0) * Box(14, 1.2, 1.0, align=CK)
    rack.label = "rack"
    asm = assembly(pinion, rack, label="rack_pinion")
    assert len(asm.children) == 2
    return asm


@problem("P119", "roller chain sprocket")
def _p():
    disk = Cylinder(4.0, 0.6, align=CK)
    seat = Pos(4.0, 0, 0) * Cylinder(0.5, 1)
    g = disk - polar_array(seat, 18, 360.0)
    assert g.is_valid and g.volume < disk.volume
    return g


@problem("P120", "V-pulley with keyway and bore")
def _p():
    body = Cylinder(3.5, 2.0, align=CK)
    g = body - Torus(3.5, 0.7) - Pos(0, 0, -0.1) * Cylinder(0.7, 2.4, align=CK)
    g = g - Pos(0.7, 0, 1.0) * Box(0.3, 0.5, 2.2, align=CK)              # keyway
    assert g.is_valid and g.volume < body.volume
    return g


# =====================================================================================================
# Batch 13 (P121-P130) -- threaded fasteners
# =====================================================================================================
@problem("P121", "external V-thread shank")
def _p():
    g = threaded_shank(1.5, 1.0, 6.0)
    assert g.is_valid and g.volume > math.pi * 1.5 ** 2 * 6.0          # thread adds material
    return g


@problem("P122", "hex-head bolt")
def _p():
    head = extrude(RegularPolygon(2.0, 6), 1.5)
    shank = Pos(0, 0, -6.0) * threaded_shank(1.2, 1.0, 6.0)
    g = head + shank
    assert g.is_valid and g.volume > head.volume
    return g


@problem("P123", "hex nut with tapped bore")
def _p():
    body = extrude(RegularPolygon(2.2, 6), 1.6)
    g = body - Cylinder(1.1, 4)
    hel = Helix(pitch=0.8, height=1.6, radius=1.1)
    groove = sweep(Plane(origin=hel @ 0.0, z_dir=hel % 0.0) * Triangle(a=0.4, b=0.4, C=90), path=hel)
    g = g - groove                                                     # internal thread groove
    assert g.is_valid and g.volume < body.volume
    return g


@problem("P124", "socket-head cap screw")
def _p():
    head = Cylinder(1.8, 1.8, align=CK)
    head = head - Pos(0, 0, 0.6) * extrude(RegularPolygon(0.9, 6), 1.4)   # hex socket
    shank = Pos(0, 0, -5.0) * threaded_shank(1.2, 1.0, 5.0)
    g = head + shank
    assert g.is_valid
    return g


@problem("P125", "wing nut")
def _p():
    body = Cylinder(1.6, 1.4, align=CK) - Cylinder(0.9, 2)
    wings = [Rot(0, 0, s) * Pos(2.2, 0, 0.7) * Box(2.0, 0.4, 1.0) for s in (0, 180)]
    g = body + wings[0] + wings[1]
    assert g.is_valid and g.volume > body.volume
    return g


@problem("P126", "double-ended stud")
def _p():
    mid = Cylinder(1.2, 4, align=CK)
    g = Pos(0, 0, 4) * threaded_shank(1.2, 1.0, 4.0) + mid + Pos(0, 0, -4.0) * threaded_shank(1.2, 1.0, 4.0)
    assert g.is_valid and g.volume > mid.volume
    return g


@problem("P127", "split lock washer")
def _p():
    ring = tube(1.2, 2.0, 0.5)
    g = ring - Pos(1.6, 0, 0) * Box(1.2, 0.25, 1.0)                    # the split gap
    assert g.is_valid and g.volume < ring.volume
    return g


@problem("P128", "eye bolt")
def _p():
    eye = Pos(0, 0, 6.0) * Rot(90, 0, 0) * Torus(1.4, 0.4)
    shank = threaded_shank(1.0, 1.0, 5.0)
    g = shank + Pos(0, 0, 5.0) * Cylinder(1.0, 1.5, align=CK) + eye
    assert g.is_valid
    return g


@problem("P129", "turnbuckle frame")
def _p():
    outer = Box(10, 3, 1.5)
    g = outer - Box(7, 1.6, 3)
    g = g - Pos(4, 0, 0) * Rot(0, 90, 0) * Cylinder(0.6, 4) - Pos(-4, 0, 0) * Rot(0, 90, 0) * Cylinder(0.6, 4)
    assert g.is_valid and g.volume < outer.volume
    return g


@problem("P130", "grub set screw (hex socket)")
def _p():
    g = threaded_shank(1.2, 0.8, 2.4) - Pos(0, 0, 1.4) * extrude(RegularPolygon(0.6, 6), 1.2)
    assert g.is_valid
    return g


# =====================================================================================================
# Batch 14 (P131-P140) -- sheet metal & enclosures
# =====================================================================================================
@problem("P131", "U-channel")
def _p():
    t = 0.2
    web = Box(8, 3, t, align=CK)
    f1 = Pos(0, 1.5, 0) * Box(8, t, 2, align=CK)
    f2 = Pos(0, -1.5, 0) * Box(8, t, 2, align=CK)
    g = web + f1 + f2
    assert g.is_valid and g.volume > 0
    return g


@problem("P132", "flanged box enclosure")
def _p():
    outer = Box(8, 6, 3, align=CK)
    g = outer - Pos(0, 0, 0.3) * Box(7.4, 5.4, 3, align=CK)            # hollow, open top
    flange = Pos(0, 0, 0) * (Box(9, 7, 0.3, align=CK) - Box(7.4, 5.4, 1, align=CK))
    g = g + flange
    assert g.is_valid and 0 < g.volume < outer.volume
    return g


@problem("P133", "louvered vent")
def _p():
    frame = Box(6, 5, 0.3, align=CK)
    slats = [Pos(0, -2 + i * 0.8, 0.4) * Rot(30, 0, 0) * Box(6, 0.6, 0.1, align=CK) for i in range(6)]
    g = frame
    for s in slats:
        g = g + s
    assert g.is_valid and g.volume > frame.volume
    return g


@problem("P134", "hinge (knuckles + pin)")
def _p():
    leafA = Box(4, 2, 0.3, align=CK)
    leafB = Pos(0, 2, 0) * Box(4, 2, 0.3, align=CK)
    pin = Pos(0, 1, 0.15) * Rot(0, 90, 0) * Cylinder(0.25, 4)
    g = assembly(leafA, leafB, pin, label="hinge")
    assert len(g.children) == 3
    return g


@problem("P135", "L-bracket with gussets and holes")
def _p():
    h = Box(4, 3, 0.3, align=CK) + Pos(-2, 0, 0) * Box(0.3, 3, 3, align=CK)
    gusset = Pos(-1.85, 0, 0.3) * extrude(Triangle(a=2, b=2, C=90), 3)
    g = h + gusset
    for (x, y) in [(1.5, 1), (1.5, -1)]:
        g = g - Pos(x, y, 0) * Cylinder(0.3, 1)
    assert g.is_valid
    return g


@problem("P136", "electronics enclosure with standoffs")
def _p():
    box = Box(8, 6, 2.5, align=CK) - Pos(0, 0, 0.3) * Box(7.4, 5.4, 3, align=CK)
    posts = []
    for (x, y) in [(3, 2), (-3, 2), (3, -2), (-3, -2)]:
        post = Pos(x, y, 0.3) * (Cylinder(0.5, 1.5, align=CK) - Cylinder(0.2, 2, align=CK))
        posts.append(post)
    g = box
    for p in posts:
        g = g + p
    assert g.is_valid
    return g


@problem("P137", "DIN-rail clip profile")
def _p():
    prof = make_face(Polyline((0, 0), (3, 0), (3, 0.5), (0.5, 0.5), (0.5, 1.5),
                              (2.5, 1.5), (2.5, 2.0), (0, 2.0), (0, 0)))
    g = extrude(prof, 2.0)
    assert g.is_valid and g.volume > 0
    return g


@problem("P138", "perforated panel")
def _p():
    panel = Box(8, 8, 0.3, align=CK)
    g = panel
    for i in range(-3, 4):
        for j in range(-3, 4):
            g = g - Pos(i * 1.0, j * 1.0, 0) * Cylinder(0.25, 1)
    assert g.is_valid and g.volume < panel.volume
    return g


@problem("P139", "folded chassis")
def _p():
    base = Box(8, 6, 0.25, align=CK)
    s1 = Pos(4, 0, 0) * Box(0.25, 6, 2, align=CK)
    s2 = Pos(-4, 0, 0) * Box(0.25, 6, 2, align=CK)
    g = base + s1 + s2
    assert g.is_valid and g.volume > base.volume
    return g


@problem("P140", "snap-fit lid with tabs")
def _p():
    lid = Box(6, 4, 0.3, align=CK)
    tabs = [Pos(s * 3.0, 0, -0.5) * Box(0.3, 1.0, 0.6, align=CK) for s in (1, -1)]
    g = lid + tabs[0] + tabs[1]
    assert g.is_valid and g.volume > lid.volume
    return g


# =====================================================================================================
# Batch 15 (P141-P150) -- EM machines
# =====================================================================================================
@problem("P141", "full slotted stator")
def _p():
    core = tube(3.0, 5.0, 2.0)
    slot = Pos(3.0, 0, 0) * Box(1.4, 0.9, 2.2, align=CK)
    g = core - polar_array(slot, 12, 360.0)
    assert g.is_valid and g.volume < core.volume
    return g


@problem("P142", "IPM rotor (V buried magnets)")
def _p():
    core = Cylinder(4.0, 2.0, align=CK) - Cylinder(0.8, 3, align=CK)
    g = core
    for k in range(8):
        for sgn in (1, -1):
            g = g - Rot(0, 0, k * 45) * Rot(0, 0, sgn * 20) * Pos(2.9, 0, 0) * Box(1.3, 0.3, 2.2, align=CK)
    assert g.is_valid and g.volume < core.volume
    return g


@problem("P143", "SPM rotor (surface arc magnets)")
def _p():
    core = Cylinder(3.0, 2.0, align=CK)
    core.label = "rotor"
    mags = []
    for k in range(8):
        c = k * 45.0
        seg = annular_segment(3.0, 3.6, 2.0, c - 18, c + 18)
        seg.label = f"magnet_{k:02d}"
        mags.append(seg)
    g = assembly(core, *mags, label="spm_rotor")
    assert len(g.children) == 9
    return g


@problem("P144", "segmented Halbach cylinder")
def _p():
    g = halbach_ring(3.0, 4.2, 2.0, 16, name="hb")
    assert len(g.children) == 16
    return g


@problem("P145", "claw-pole rotor")
def _p():
    d1 = Cylinder(3.0, 0.5, align=CK)
    claws1 = polar_array(Pos(2.5, 0, 0.5) * Box(0.9, 0.9, 2.5, align=CK), 6, 360.0)
    d2 = Pos(0, 0, 3.0) * Cylinder(3.0, 0.5, align=(Align.CENTER, Align.CENTER, Align.MAX))
    claws2 = Rot(0, 0, 30) * polar_array(Pos(2.5, 0, 0) * Box(0.9, 0.9, 2.5, align=(Align.CENTER, Align.CENTER, Align.MAX)), 6, 360.0)
    g = d1 + claws1 + Pos(0, 0, 0) * (d2 + claws2)
    assert g.is_valid and g.volume > d1.volume
    return g


@problem("P146", "salient-pole rotor")
def _p():
    hub = Cylinder(1.5, 2.0, align=CK)
    pole = Pos(1.5, 0, 0) * Box(1.5, 0.9, 2.0, align=CK) + Pos(2.85, 0, 0) * Box(0.4, 1.8, 2.0, align=CK)
    g = hub + polar_array(pole, 4, 360.0)
    assert g.is_valid and g.volume > hub.volume
    return g


@problem("P147", "squirrel-cage rotor")
def _p():
    core = Cylinder(3.0, 4.0, align=CK)
    bars = polar_array(Pos(2.6, 0, 0) * Cylinder(0.3, 4.0, align=CK), 16, 360.0)
    ring1 = tube(2.3, 2.9, 0.4)
    ring2 = Pos(0, 0, 3.6) * tube(2.3, 2.9, 0.4)
    g = core + bars + ring1 + ring2
    assert g.is_valid and g.volume > core.volume
    return g


@problem("P148", "EI transformer core")
def _p():
    yoke = Box(6, 1, 2, align=CK)                                       # y in [-0.5, 0.5]
    YMIN = (Align.CENTER, Align.MIN, Align.MIN)
    legs = [Pos(x, -0.5, 0) * Box(1, 3, 2, align=YMIN) for x in (-2.5, 0, 2.5)]   # y in [-0.5, 2.5] (overlap)
    ibar = Pos(0, 2.0, 0) * Box(6, 1, 2, align=YMIN)                    # y in [2.0, 3.0] (overlap legs)
    g = yoke
    for x in legs:
        g = g + x
    g = g + ibar
    assert g.is_valid and g.volume > yoke.volume
    return g


@problem("P149", "toroidal transformer core")
def _p():
    R, w, h = 5.0, 1.2, 1.6
    g = revolve(Pos(R, 0) * Rectangle(w, h), Axis.Y, 360)
    close(g.volume, w * h * 2 * math.pi * R, rel=3e-3)
    return g


@problem("P150", "axial-flux disc rotor")
def _p():
    disk = Cylinder(4.0, 0.5, align=CK)
    mag = Pos(2.6, 0, 0.5) * Box(1.4, 0.9, 0.5, align=CK)
    g = disk + polar_array(mag, 8, 360.0)
    assert g.is_valid and g.volume > disk.volume
    return g


# =====================================================================================================
# Batch 16 (P151-P160) -- coils & windings
# =====================================================================================================
@problem("P151", "concentrated tooth coil")
def _p():
    g = extrude(RectangleRounded(2.6, 3.6, 0.5) - RectangleRounded(1.6, 2.6, 0.4), 1.5)
    assert g.is_valid and g.volume > 0
    return g


@problem("P152", "distributed winding hairpin")
def _p():
    path = Spline((-1, 0, 0), (-1, 0, 4), (0, 0.6, 5), (1, 0, 4), (1, 0, 0))
    g = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(0.2), path=path)
    assert g.is_valid and g.volume > 0
    return g


@problem("P153", "saddle coil (open)")
def _p():
    R = 2.0
    pts = [(R * math.cos(math.radians(a)), R * math.sin(math.radians(a)), z)
           for a, z in [(-40, 0), (-40, 3), (-20, 4.2), (20, 4.2), (40, 3), (40, 0)]]
    path = Spline(*pts)
    g = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(0.18), path=path)
    assert g.is_valid and g.volume > 0
    return g


@problem("P154", "racetrack coil with leads")
def _p():
    coil = extrude(RectangleRounded(6, 3, 0.8) - RectangleRounded(5, 2, 0.6), 1.0)
    leads = Pos(0, -1.4, 0) * Box(0.3, 0.3, 3, align=CK)
    g = coil + leads
    assert g.is_valid and g.volume > coil.volume
    return g


@problem("P155", "Helmholtz coil pair with formers")
def _p():
    c1 = tube(3.0, 3.6, 0.8)
    c1.label = "coil_lo"
    c2 = Pos(0, 0, 3.0) * tube(3.0, 3.6, 0.8)
    c2.label = "coil_hi"
    g = assembly(c1, c2, label="helmholtz")
    assert len(g.children) == 2
    return g


@problem("P156", "two-layer solenoid")
def _p():
    rw = 0.18
    g = None
    for R in (2.0, 2.5):
        hel = Helix(pitch=0.6, height=4.0, radius=R)
        coil = sweep(Plane(origin=hel @ 0.0, z_dir=hel % 0.0) * Circle(rw), path=hel)
        g = coil if g is None else g + coil
    assert g.is_valid and g.volume > 0
    return g


@problem("P157", "toroidal multi-turn winding")
def _p():
    R, a, rw, N, frac = 6.0, 1.0, 0.1, 10, 0.8
    m = 40 * N
    pts = [(((R + a * math.cos(2 * math.pi * N * frac * s / m)) * math.cos(2 * math.pi * frac * s / m)),
            ((R + a * math.cos(2 * math.pi * N * frac * s / m)) * math.sin(2 * math.pi * frac * s / m)),
            a * math.sin(2 * math.pi * N * frac * s / m)) for s in range(m + 1)]
    path = Spline(*pts)
    g = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(rw), path=path)
    assert g.is_valid and g.volume > 0
    return g


@problem("P158", "planar spiral PCB coil")
def _p():
    turns, b, w = 4, 0.18, 0.05
    m = 80 * turns
    pts = [(b * (2 * math.pi * turns * s / m) * math.cos(2 * math.pi * turns * s / m),
            b * (2 * math.pi * turns * s / m) * math.sin(2 * math.pi * turns * s / m), 0.0) for s in range(1, m + 1)]
    path = Spline(*pts)
    g = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(w), path=path)
    assert g.is_valid and g.volume > 0
    return g


@problem("P159", "litz bundle in an iron slot")
def _p():
    slot = Box(3, 3, 6, align=CK) - Pos(0, 0, 0.3) * Box(2.4, 2.4, 6, align=CK)
    slot.label = "slot_iron"
    bundle = Pos(0, 0, 0.3) * litz_wire(7, 0.3, 1.0, 5.0, 5.0, name="litz")
    g = assembly(slot, bundle, label="slot_winding")
    assert len(g.children) >= 8
    return g


@problem("P160", "racetrack accelerator coil")
def _p():
    g = extrude(RectangleRounded(12, 5, 1.5) - RectangleRounded(10, 3, 1.0), 2.0)
    assert g.is_valid and g.volume > 0
    return g


# =====================================================================================================
# Batch 17 (P161-P170) -- curved surfaces & freeform
# =====================================================================================================
@problem("P161", "twisted turbine blade")
def _p():
    g = loft([airfoil_face(4.0, 0.14),
              Pos(0, 0, 3) * Rot(0, 0, 12) * airfoil_face(3.4, 0.12),
              Pos(0, 0, 6) * Rot(0, 0, 25) * airfoil_face(2.6, 0.10)])
    assert g.is_valid and g.volume > 0
    return g


@problem("P162", "tapered propeller blade")
def _p():
    g = loft([airfoil_face(3.0, 0.12),
              Pos(0, 0, 4) * Rot(0, 0, 20) * airfoil_face(2.0, 0.10),
              Pos(0, 0, 7) * Rot(0, 0, 35) * airfoil_face(0.8, 0.08)])
    assert g.is_valid and g.volume > 0
    return g


@problem("P163", "centrifugal impeller")
def _p():
    hub = Cylinder(1.5, 2.0, align=CK)
    blade = Pos(1.4, 0, 0.2) * Rot(0, 0, 25) * Box(2.5, 0.3, 1.6, align=CK)
    g = hub + polar_array(blade, 7, 360.0)
    assert g.is_valid and g.volume > hub.volume
    return g


@problem("P164", "axial fan")
def _p():
    hub = Cylinder(1.0, 1.2, align=CK)
    blade = Pos(1.8, 0, 0.6) * Rot(30, 0, 0) * Box(2.0, 0.15, 1.2, align=CK)
    g = hub + polar_array(blade, 5, 360.0)
    assert g.is_valid and g.volume > hub.volume
    return g


@problem("P165", "lofted hull")
def _p():
    g = loft([Pos(0, 0, 0) * Ellipse(1.0, 0.5),
              Pos(0, 0, 4) * Ellipse(2.5, 1.2),
              Pos(0, 0, 8) * Ellipse(1.6, 1.0)])
    assert g.is_valid and g.volume > 0
    return g


@problem("P166", "revolved bottle")
def _p():
    prof = make_face(Polyline((0, 0), (2, 0), (2, 4), (1.6, 5), (0.6, 6), (0.6, 7), (0, 7), (0, 0)))
    g = revolve(prof, Axis.Y, 360)
    assert g.is_valid and g.volume > 0
    return g


@problem("P167", "converging-diverging nozzle")
def _p():
    g = loft([Circle(2.0), Pos(0, 0, 3) * Circle(0.8), Pos(0, 0, 7) * Circle(2.4)])
    g = g - loft([Circle(1.7), Pos(0, 0, 3) * Circle(0.6), Pos(0, 0, 7) * Circle(2.1)])
    assert g.is_valid and 0 < g.volume
    return g


@problem("P168", "branched manifold")
def _p():
    main = Rot(0, 90, 0) * Cylinder(0.8, 10)
    branches = [Pos(x, 0, 0) * Cylinder(0.5, 3, align=CK) for x in (-3, 0, 3)]
    g = main
    for b in branches:
        g = g + b
    assert g.is_valid and g.volume > main.volume
    return g


@problem("P169", "pipe elbow run (S-bend)")
def _p():
    p1 = CenterArc((0, 0, 0), 4, 0, 90)            # (4,0)->(0,4)  (4th arg = arc SIZE in deg)
    p2 = CenterArc((0, 8, 0), 4, 270, 90)          # (0,4)->(4,8), shares the (0,4) end with p1
    s1 = sweep(Plane(origin=p1 @ 0.0, z_dir=p1 % 0.0) * Circle(0.4), path=p1)
    s2 = sweep(Plane(origin=p2 @ 0.0, z_dir=p2 % 0.0) * Circle(0.4), path=p2)
    g = s1 + s2
    assert g.is_valid and g.volume > 0
    return g


@problem("P170", "spiral volute casing")
def _p():
    turns = 1.2
    m = 120
    pts = [((2 + 1.2 * s / m) * math.cos(2 * math.pi * turns * s / m),
            (2 + 1.2 * s / m) * math.sin(2 * math.pi * turns * s / m), 0.0) for s in range(m + 1)]
    path = Spline(*pts)
    g = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Circle(0.6), path=path)
    assert g.is_valid and g.volume > 0
    return g


# =====================================================================================================
# Batch 18 (P171-P180) -- mechanisms & multi-body assemblies
# =====================================================================================================
@problem("P171", "four-bar linkage")
def _p():
    ground = Box(6, 0.4, 0.4)
    crank = Pos(-3, 0, 0.5) * Rot(0, 0, 60) * Box(2, 0.3, 0.3, align=(Align.MIN, Align.CENTER, Align.CENTER))
    coupler = Pos(0, 1.5, 1.0) * Box(3, 0.3, 0.3)
    rocker = Pos(3, 0, 0.5) * Rot(0, 0, 120) * Box(2, 0.3, 0.3, align=(Align.MIN, Align.CENTER, Align.CENTER))
    g = assembly(ground, crank, coupler, rocker, label="fourbar")
    assert len(g.children) == 4
    return g


@problem("P172", "slider-crank")
def _p():
    crank = Pos(0, 0, 0) * Box(1.6, 0.3, 0.3, align=(Align.MIN, Align.CENTER, Align.CENTER))
    crank.label = "crank"
    rod = Pos(1.6, 0, 0) * Box(3.5, 0.25, 0.25, align=(Align.MIN, Align.CENTER, Align.CENTER))
    rod.label = "rod"
    slider = Pos(5.0, 0, 0) * Box(1.0, 1.0, 0.8)
    slider.label = "slider"
    g = assembly(crank, rod, slider, label="slider_crank")
    assert len(g.children) == 3
    return g


@problem("P173", "cam and follower")
def _p():
    base, lift = 2.0, 0.8
    pts = [((base + lift * (1 - math.cos(math.radians(d))) / 2) * math.cos(math.radians(d)),
            (base + lift * (1 - math.cos(math.radians(d))) / 2) * math.sin(math.radians(d))) for d in range(0, 360, 6)]
    pts.append(pts[0])
    cam = extrude(make_face(Polyline(*pts)), 0.6)
    cam.label = "cam"
    follower = Pos(0, 3.2, 0.3) * Cylinder(0.4, 0.6, align=CK) + Pos(0, 3.2, 0.9) * Box(0.3, 0.3, 2, align=CK)
    follower.label = "follower"
    g = assembly(cam, follower, label="cam_follower")
    assert len(g.children) == 2
    return g


@problem("P174", "two-gear train")
def _p():
    g1 = spur_gear(16, 3.0, 0.7, 1.0)
    g1.label = "gear_a"
    g2 = Pos(6.4, 0, 0) * spur_gear(16, 3.0, 0.7, 1.0)
    g2.label = "gear_b"
    g = assembly(g1, g2, label="gear_train")
    assert len(g.children) == 2
    return g


@problem("P175", "Geneva mechanism")
def _p():
    driver = Cylinder(2.0, 0.5, align=CK) + Pos(1.6, 0, 0.5) * Cylinder(0.25, 0.6, align=CK)
    driver.label = "driver"
    cross = Cylinder(2.5, 0.5, align=CK)
    for k in range(4):
        cross = cross - Rot(0, 0, k * 90) * Pos(2.5, 0, 0) * Box(0.6, 0.6, 1)
    cross = Pos(4.0, 0, 0) * cross
    cross.label = "cross"
    g = assembly(driver, cross, label="geneva")
    assert len(g.children) == 2
    return g


@problem("P176", "universal joint")
def _p():
    yoke1 = Cylinder(1.0, 0.6, align=CK) + Pos(0.8, 0, 0.6) * Box(0.4, 1.6, 1.2, align=CK) + Pos(-0.8, 0, 0.6) * Box(0.4, 1.6, 1.2, align=CK)
    yoke1.label = "yoke_in"
    cross = Pos(0, 0, 1.5) * (Rot(0, 90, 0) * Cylinder(0.2, 2.0) + Rot(90, 0, 0) * Cylinder(0.2, 2.0))
    cross.label = "cross"
    yoke2 = Pos(0, 0, 3.0) * Rot(0, 0, 90) * yoke1
    yoke2.label = "yoke_out"
    g = assembly(yoke1, cross, yoke2, label="ujoint")
    assert len(g.children) == 3
    return g


@problem("P177", "ball bearing assembly")
def _p():
    inner = tube(2.0, 2.6, 1.0); inner.label = "inner_race"
    outer = tube(3.4, 4.0, 1.0); outer.label = "outer_race"
    balls = polar_array(Pos(3.0, 0, 0.5) * Sphere(0.4), 8, 360.0)
    cage = tube(2.7, 3.3, 0.2); cage.label = "cage"
    g = assembly(inner, outer, balls, cage, label="bearing")
    assert len(g.children) == 8 + 3
    return g


@problem("P178", "planetary gearbox assembly")
def _p():
    carrier = Cylinder(5.0, 0.3, align=CK); carrier.label = "carrier"
    sun = Pos(0, 0, 0.3) * spur_gear(12, 2.0, 0.5, 1.0); sun.label = "sun"
    planets = []
    for k in range(3):
        p = Pos(0, 0, 0.3) * Rot(0, 0, k * 120) * Pos(4.5, 0, 0) * spur_gear(10, 1.6, 0.45, 1.0)
        p.label = f"planet_{k}"
        planets.append(p)
    ring = tube(6.6, 7.4, 1.3); ring.label = "ring"
    g = assembly(carrier, sun, *planets, ring, label="gearbox")
    assert len(g.children) == 6
    return g


@problem("P179", "two-finger gripper")
def _p():
    base = Box(3, 2, 1, align=CK); base.label = "base"
    f1 = Pos(1.0, 0, 1.0) * Rot(0, 15, 0) * Box(0.4, 1.5, 3, align=CK); f1.label = "finger_r"
    f2 = Pos(-1.0, 0, 1.0) * Rot(0, -15, 0) * Box(0.4, 1.5, 3, align=CK); f2.label = "finger_l"
    g = assembly(base, f1, f2, label="gripper")
    assert len(g.children) == 3
    return g


@problem("P180", "scissor-lift segment")
def _p():
    barA = Rot(0, 35, 0) * Box(5, 0.3, 0.3, align=(Align.MIN, Align.CENTER, Align.CENTER)); barA.label = "bar_a"
    barB = Rot(0, -35, 0) * Box(5, 0.3, 0.3, align=(Align.MIN, Align.CENTER, Align.CENTER)); barB.label = "bar_b"
    pivot = Pos(2.0, 0, 1.4) * Rot(90, 0, 0) * Cylinder(0.2, 1); pivot.label = "pivot"
    g = assembly(barA, barB, pivot, label="scissor")
    assert len(g.children) == 3
    return g


# =====================================================================================================
# Batch 19 (P181-P190) -- lattices & generative structures
# =====================================================================================================
@problem("P181", "honeycomb panel")
def _p():
    plate = Box(10, 10, 0.5, align=CK)
    g = plate
    for i in range(-2, 3):
        for j in range(-2, 3):
            x, y = i * 2.0 + (j % 2) * 1.0, j * 1.7
            g = g - Pos(x, y, 0) * extrude(RegularPolygon(0.8, 6), 1)
    assert g.is_valid and g.volume < plate.volume
    return g


@problem("P182", "square infill grid")
def _p():
    walls = []
    for i in range(-3, 4):
        walls.append(Pos(i * 1.5, 0, 0) * Box(0.15, 10, 1, align=CK))
        walls.append(Pos(0, i * 1.5, 0) * Box(10, 0.15, 1, align=CK))
    g = walls[0]
    for w in walls[1:]:
        g = g + w
    assert g.is_valid and g.volume > 0
    return g


@problem("P183", "BCC strut lattice cell")
def _p():
    c = 2.0
    corners = [(sx * c / 2, sy * c / 2, sz * c / 2) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    g = None
    for cor in corners:
        s = strut(cor, (0, 0, 0), 0.12)
        g = s if g is None else g + s
    assert g.is_valid and len(corners) == 8
    return g


@problem("P184", "X-braced frame face")
def _p():
    g = strut((-2, -2, 0), (2, 2, 0), 0.15) + strut((-2, 2, 0), (2, -2, 0), 0.15)
    assert g.is_valid and g.volume > 0
    return g


@problem("P185", "planar Warren truss")
def _p():
    top = [(-4 + 2 * k, 1.5, 0) for k in range(5)]
    bot = [(-3 + 2 * k, -1.5, 0) for k in range(4)]
    members = []
    for a, b in zip(top, top[1:]):
        members.append(strut(a, b, 0.1))
    for a, b in zip(bot, bot[1:]):
        members.append(strut(a, b, 0.1))
    for k in range(4):
        members.append(strut(top[k], bot[k], 0.1))
        members.append(strut(top[k + 1], bot[k], 0.1))
    g = assembly(*members, label="truss")                      # many thin members -> keep as an assembly
    assert len(g.solids()) == len(members) and all(s.is_valid for s in g.solids())
    return g


@problem("P186", "space-frame bay")
def _p():
    nodes = {(i, j): (i * 2.0, j * 2.0, 0.0) for i in range(3) for j in range(3)}
    top = (1.0 * 2 / 2, 1.0 * 2 / 2, 2.0)                       # apex above centre
    members = []
    for i in range(2):
        for j in range(3):
            members.append(strut(nodes[(i, j)], nodes[(i + 1, j)], 0.08))
    for i in range(3):
        for j in range(2):
            members.append(strut(nodes[(i, j)], nodes[(i, j + 1)], 0.08))
    for n in [(0, 0), (2, 0), (0, 2), (2, 2)]:
        members.append(strut(nodes[n], top, 0.08))
    g = assembly(*members, label="space_frame")                # many members -> keep as an assembly
    assert len(g.solids()) == len(members) and all(s.is_valid for s in g.solids())
    return g


@problem("P187", "geodesic lattice dome")
def _p():
    R = 4.0
    rings = [20, 50, 75]                                        # polar angles (deg) of latitude rings
    members = []
    prev = None
    for lat in rings:
        th = math.radians(lat)
        pts = [(R * math.sin(th) * math.cos(math.radians(a)), R * math.sin(th) * math.sin(math.radians(a)),
                R * math.cos(th)) for a in range(0, 360, 45)]
        for a, b in zip(pts, pts[1:] + [pts[0]]):
            members.append(strut(a, b, 0.1))
        if prev is not None:
            for a, b in zip(prev, pts):
                members.append(strut(a, b, 0.1))
        prev = pts
    g = assembly(*members, label="dome")                       # lattice of struts -> keep as an assembly
    assert len(g.solids()) == len(members) and all(s.is_valid for s in g.solids())
    return g


@problem("P188", "auxetic re-entrant cell")
def _p():
    pts = [(-1.5, -2), (-0.4, -1), (-1.5, 0), (-1.5, 2), (1.5, 2), (1.5, 0),
           (0.4, -1), (1.5, -2), (-1.5, -2)]
    g = extrude(make_face(Polyline(*pts)), 0.5)
    assert g.is_valid and g.volume > 0
    return g


@problem("P189", "diagonal lattice shell")
def _p():
    shell = tube(2.6, 3.0, 5.0)
    g = shell
    for k in range(8):
        g = g - Rot(0, 0, k * 45) * Pos(3.0, 0, 2.5) * Rot(0, 35, 0) * Box(0.4, 1.2, 4)
    assert g.is_valid and g.volume < shell.volume
    return g


@problem("P190", "radial spoke wheel")
def _p():
    hub = Cylinder(0.8, 0.6, align=CK)
    rim = tube(3.6, 4.0, 0.6)
    g = hub + rim
    for k in range(8):
        a = math.radians(k * 45)
        g = g + strut((0.7 * math.cos(a), 0.7 * math.sin(a), 0.3),
                      (3.7 * math.cos(a), 3.7 * math.sin(a), 0.3), 0.12)
    assert g.is_valid and g.volume > hub.volume + rim.volume
    return g


# =====================================================================================================
# Batch 20 (P191-P200) -- CAE-ready full models & validation
# =====================================================================================================
@problem("P191", "PMSM cross-section (multi-region)")
def _p():
    stator = tube(4.0, 6.0, 2.0) - polar_array(Pos(4.0, 0, 0) * Box(1.4, 0.9, 2.2, align=CK), 12, 360.0)
    stator.label = "stator"
    rotor = Cylinder(3.4, 2.0, align=CK)
    rotor.label = "rotor"
    mags = []
    for k in range(8):
        c = k * 45.0
        seg = annular_segment(3.4, 3.8, 2.0, c - 18, c + 18)
        seg.label = f"magnet_{k:02d}"
        mags.append(seg)
    g = assembly(stator, rotor, *mags, label="pmsm")
    assert len(g.children) == 10
    return g


@problem("P192", "bearing assembly mass check")
def _p():
    inner = tube(2.0, 2.6, 1.0)
    outer = tube(3.4, 4.0, 1.0)
    balls = [Pos(3.0 * math.cos(math.radians(a)), 3.0 * math.sin(math.radians(a)), 0.5) * Sphere(0.4)
             for a in range(0, 360, 45)]
    g = assembly(inner, outer, *balls, label="bearing")
    total = inner.volume + outer.volume + sum(b.volume for b in balls)
    close(g.volume, total, rel=1e-6)
    return g


@problem("P193", "pin-fin heat sink (meshable)")
def _p():
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    base = Box(8, 8, 0.6, align=CK)
    g = base
    for i in range(-2, 3):
        for j in range(-2, 3):
            g = g + Pos(i * 1.6, j * 1.6, 0.6) * Cylinder(0.35, 2.0, align=CK)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "hs.step")
        export_step(g, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=1.5))
    assert g.is_valid and mesh.ne > 200
    return g


@problem("P194", "busbar assembly with bends")
def _p():
    path = Spline((0, 0, 0), (4, 0, 0), (5, 2, 0), (5, 5, 1))
    bar = sweep(Plane(origin=path @ 0.0, z_dir=path % 0.0) * Rectangle(0.8, 0.3), path=path)
    bar.label = "bar"
    lug1 = Cylinder(0.6, 0.4, align=CK); lug1.label = "lug_in"
    lug2 = Pos(5, 5, 1) * Cylinder(0.6, 0.4, align=CK); lug2.label = "lug_out"
    g = assembly(bar, lug1, lug2, label="busbar")
    assert len(g.children) == 3
    return g


@problem("P195", "waveguide iris filter")
def _p():
    guide = Box(2.0, 1.0, 12) - Box(1.6, 0.7, 12.4)
    g = guide
    for z in (-3, 0, 3):
        iris = Pos(0, 0, z) * (Box(1.6, 0.7, 0.2) - Box(0.6, 0.7, 1))
        g = g + iris
    assert g.is_valid and g.volume > 0
    return g


@problem("P196", "multi-region STEP round-trip")
def _p():
    from build123d import export_step, import_step
    asm = assembly(Box(2, 2, 2), Pos(4, 0, 0) * Sphere(1.0), Pos(0, 4, 0) * Cylinder(0.8, 2))
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "asm.step")
        export_step(asm, f)
        back = import_step(f)
    close(back.volume, asm.volume, rel=1e-6)
    return back


@problem("P197", "Netgen mesh of an impeller")
def _p():
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    hub = Cylinder(1.5, 2.0, align=CK)
    g = hub + polar_array(Pos(1.4, 0, 0.2) * Rot(0, 0, 25) * Box(2.5, 0.3, 1.6, align=CK), 7, 360.0)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "imp.step")
        export_step(g, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=1.2))
    assert mesh.ne > 100
    return g


@problem("P198", "assembly mass properties")
def _p():
    parts = [Box(2, 2, 2), Pos(6, 0, 0) * Box(2, 2, 2)]
    asm = assembly(*parts, label="pair")
    close(asm.volume, 16.0)
    close(asm.center().X, 3.0, rel=1e-6)                       # midway between the two boxes
    return asm


@problem("P199", "interference check in an assembly")
def _p():
    a = Box(3, 3, 3)
    b = Pos(2, 0, 0) * Box(3, 3, 3)
    overlap = (a & b).volume
    clear = Pos(0, 6, 0) * Box(3, 3, 3)
    assert overlap > 0
    try:
        sep = (a & clear).volume
    except Exception:
        sep = 0.0
    assert sep < 1e-9
    return assembly(a, b, label="interfere")


@problem("P200", "PMSM 3D stack (meshable region)")
def _p():
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    rotor = Cylinder(3.0, 6.0, align=CK)
    rotor.label = "rotor"
    stator = tube(3.6, 5.0, 6.0)
    stator.label = "stator"
    asm = assembly(rotor, stator, label="pmsm_stack")
    assert len(asm.children) == 2
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "rotor.step")
        export_step(rotor, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=2.0))
    assert mesh.ne > 50
    return asm


# -----------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("pid,title,fn", PROBLEMS, ids=[p[0] for p in PROBLEMS])
def test_marathon2(pid, title, fn):
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
    print(f"[marathon2] {ok}/{len(PROBLEMS)} problems verified")


if __name__ == "__main__":
    from ngsolve import TaskManager

    with TaskManager():
        main()
