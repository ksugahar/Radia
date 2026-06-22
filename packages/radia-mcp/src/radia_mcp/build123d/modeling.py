# -*- coding: utf-8 -*-
r"""Parametric solid-modeling operations for CAE/EM geometry in build123d.

A small, **tested** library of the parametric modeling operations a solid-modeller user reaches for
when authoring electromagnetic structures -- the building blocks that build123d has as raw kernel
calls but that are fiddly to get clean and meshable by hand: a robust **annular wedge/sector**, a
hollow **tube**, a **racetrack** (rounded-rectangle) coil, and the **array-with-copies** transforms
(polar / linear / mirror) plus a labelled multi-region **assembly**.  Every helper returns a clean,
``is_valid`` build123d solid (or a labelled :class:`~build123d.Compound`), CAE-safe (primitives +
boolean, no micro-edges) and ready for the Netgen -> Radia / NGSolve tet pipeline.

Design rules (see build123d_usage "cae_guidelines"):
  * primitives + boolean only -- no tiny fillets, no over-constrained sketches;
  * every solid carries a ``.label`` (maps to a material region / Gmsh physical group);
  * z-centred solids (span -h/2..h/2) so stacking and symmetry planes are simple;
  * pure build123d -- no radia_mcp dependency, so the helpers are equally usable inside an
    ``execute_build123d`` subprocess or copy-pasted into a user script.

These are intentionally generic CAD operations (sector, array, tube, racetrack); the EM *device*
archetypes that compose them (Halbach ring, dipole/quadrupole yokes, C-core, solenoid ...) live in
:mod:`radia_mcp.build123d.archetypes`.
"""
from __future__ import annotations

import math

from build123d import (Axis, BuildLine, BuildSketch, CenterArc, Circle, Compound, Cylinder, Helix, Keep,
                       Line, Mode, Plane, Pos, Rectangle, RectangleRounded, Spline, chamfer, extrude,
                       fillet, loft, make_face, offset, revolve, split, sweep)

__all__ = ["annular_segment", "tube", "racetrack_coil", "polar_array", "linear_array",
           "mirrored", "assembly",
           # generic solid-modelling operations (constructors / local mods / arrays)
           "swept", "revolved", "lofted", "coil", "strut", "thicken", "draft_extrude",
           "shell", "fillet_edges", "chamfer_edges", "grid_array", "path_array",
           # boolean / slice / sheet-metal
           "fuse", "cut", "common", "slice_solid", "bend_sheet"]


def annular_segment(r_in, r_out, h, start_angle=0.0, end_angle=90.0, label="segment"):
    r"""A radial **wedge of an annulus**: ``r_in < r < r_out``, angular span
    ``[start_angle, end_angle]`` degrees, height ``h`` (z-centred).

    Robust for any span ``< 360`` deg: the full annulus is INTERSECTED with a pie wedge swept from the
    centre, so there are no self-intersecting sketch boundaries (the failure mode of a hand-built
    sector polyline).  The canonical building block for segmented PM arrays (Halbach), pole sectors,
    and angular slices.  Volume ``= (span/360) * pi (r_out^2 - r_in^2) * h``.
    """
    if not (0.0 < r_in < r_out):
        raise ValueError("require 0 < r_in < r_out")
    span = end_angle - start_angle
    if not (0.0 < span < 360.0):
        raise ValueError("require 0 < (end_angle - start_angle) < 360")
    annulus = Cylinder(radius=r_out, height=h) - Cylinder(radius=r_in, height=h)
    big = r_out * 1.5
    a0, a1 = math.radians(start_angle), math.radians(end_angle)
    with BuildSketch(Plane.XY) as sk:
        with BuildLine():
            Line((0, 0), (big * math.cos(a0), big * math.sin(a0)))
            CenterArc((0, 0), big, start_angle, span)
            Line((big * math.cos(a1), big * math.sin(a1)), (0, 0))
        make_face()
    wedge = Pos(0, 0, -h / 2) * extrude(sk.sketch, amount=h)
    seg = (annulus & wedge).solid()        # single labelled Solid (label survives .solids() downstream)
    seg.label = label
    return seg


def tube(r_in, r_out, h, label="tube"):
    r"""Hollow cylinder (annulus extruded), z-centred.  Volume ``= pi (r_out^2 - r_in^2) h``."""
    if not (0.0 <= r_in < r_out):
        raise ValueError("require 0 <= r_in < r_out")
    t = (Cylinder(radius=r_out, height=h) - Cylinder(radius=r_in, height=h)).solid()
    t.label = label
    return t


def racetrack_coil(length, width, band, h, corner_radius, label="coil"):
    r"""A **racetrack** (rounded-rectangle) coil: an outer rounded rectangle minus an inner one,
    extruded by ``h`` (z-centred).  ``length`` x ``width`` is the outer envelope, ``band`` the
    conductor wall thickness, ``corner_radius`` the outer corner radius; the conductor cross-section is
    ``band`` (in-plane) x ``h`` (axial).  The standard flat racetrack / pancake coil outline used for
    accelerator and MRI windings.
    """
    if band <= 0 or band * 2 >= min(length, width):
        raise ValueError("require 0 < 2*band < min(length, width)")
    with BuildSketch(Plane.XY) as sk:
        RectangleRounded(length, width, corner_radius)
        RectangleRounded(length - 2 * band, width - 2 * band,
                         max(corner_radius - band, 0.1), mode=Mode.SUBTRACT)
    coil = (Pos(0, 0, -h / 2) * extrude(sk.sketch, amount=h)).solid()
    coil.label = label
    return coil


def polar_array(part, count, total_angle=360.0, axis=Axis.Z, label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` rotated copies of ``part`` about ``axis``, returned as a labelled
    :class:`~build123d.Compound`.  A full ``360`` deg ring spaces the copies by ``360/count``; a
    partial fan (``total_angle < 360``) spaces them by ``total_angle/(count-1)`` so the first and last
    copies sit at the fan ends.  The "rotate with copies" / segmented-array verb.
    """
    if count < 1:
        raise ValueError("count must be >= 1")
    base = label or (part.label or "part")
    full = abs(total_angle) >= 360.0 - 1e-9
    step = total_angle / count if full else (total_angle / (count - 1) if count > 1 else 0.0)
    children = []
    for k in range(count):
        c = part.rotate(axis, k * step)
        c.label = label_fmt.format(base=base, k=k)
        children.append(c)
    return Compound(children=children, label=base + "_array")


def linear_array(part, count, spacing, direction=(1, 0, 0), label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` translated copies of ``part`` along (normalised) ``direction`` at ``spacing``,
    returned as a labelled :class:`~build123d.Compound`.  The "translate with copies" verb.
    """
    if count < 1:
        raise ValueError("count must be >= 1")
    base = label or (part.label or "part")
    d = direction
    norm = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) or 1.0
    ux, uy, uz = d[0] / norm, d[1] / norm, d[2] / norm
    children = []
    for k in range(count):
        c = Pos(k * spacing * ux, k * spacing * uy, k * spacing * uz) * part
        c.label = label_fmt.format(base=base, k=k)
        children.append(c)
    return Compound(children=children, label=base + "_array")


def mirrored(part, about=Plane.XZ, keep_original=True, label=None):
    r"""The mirror image of ``part`` across plane ``about`` (default ``Plane.XZ``); with
    ``keep_original`` (default) returns a :class:`~build123d.Compound` of the original + its mirror --
    the symmetry-completion verb (build a quarter/half model, then mirror to whole).
    """
    from build123d import mirror as _mirror
    base = label or (part.label or "part")
    mir = _mirror(part, about=about)
    mir.label = base + "_mir"
    if not keep_original:
        return mir
    orig = part
    try:
        orig.label = base + "_orig"
    except Exception:
        pass
    return Compound(children=[orig, mir], label=base + "_mirrored")


def assembly(*parts, label="assembly"):
    r"""Group ``parts`` (labelled solids and/or Compounds) into a single multi-region
    :class:`~build123d.Compound`, KEEPING each part as a separate **labelled child** (do NOT fuse) so
    it maps to its own material region / Gmsh physical group downstream.  Flattens one level of nested
    Compounds via ``.children`` (NOT ``.solids()`` -- the latter drops the region labels).
    """
    children = []
    for p in parts:
        if isinstance(p, Compound) and p.children:
            children.extend(p.children)          # labelled sub-parts: keep them (solids() strips labels)
        elif isinstance(p, Compound):
            sols = p.solids()                    # raw Part (Box / extrude / ... has no children): its solids
            if p.label and len(sols) == 1:
                sols[0].label = p.label          # carry the Part's own label onto its (stripped) solid
            children.extend(sols)
        else:
            children.append(p)                   # a bare Solid
    return Compound(children=children, label=label)


# =====================================================================================================
# Generic solid-modelling operations -- the constructor / local-modification / array verbs a full 3D
# modeller exposes, wrapped as clean, labelled, CAE-safe build123d helpers.
# =====================================================================================================
def swept(profile, path, label="swept"):
    r"""**Sweep** a 2D ``profile`` (a face / sketch such as ``Circle(r)`` or ``Rectangle(w, h)``) along a
    ``path`` (a wire / edge -- a line, arc, spline or :class:`~build123d.Helix`), keeping the profile
    perpendicular to the path tangent at the start.  The "sweep curve" verb -- pipes, coil conductors,
    swept beams.  Returns a single labelled :class:`~build123d.Solid`."""
    sec = Plane(origin=path @ 0.0, z_dir=path % 0.0) * profile
    s = sweep(sec, path=path).solid()
    s.label = label
    return s


def revolved(profile, axis=Axis.Y, angle=360.0, label="revolved"):
    r"""**Revolve** (spin) a 2D ``profile`` about ``axis`` by ``angle`` degrees -- bodies of revolution
    (shafts, vases, toroidal cores, pulleys).  Returns a labelled :class:`~build123d.Solid`."""
    s = revolve(profile, axis, angle).solid()
    s.label = label
    return s


def lofted(sections, label="lofted"):
    r"""**Loft** (blend) a list of 2D ``sections`` (faces, positioned at their stations) into a solid that
    interpolates between them -- transitions, blades, ducts, frusta.  Returns a labelled
    :class:`~build123d.Solid`."""
    s = loft(sections).solid()
    s.label = label
    return s


def coil(profile, pitch, height, radius, label="coil"):
    r"""A **helical coil**: sweep a 2D ``profile`` (e.g. ``Circle(wire_radius)``) along a helix of the
    given ``pitch`` (axial advance per turn), ``height`` and ``radius``.  Solenoids, springs, helical
    conductors.  Returns a labelled :class:`~build123d.Solid` (volume ``= area(profile) * helix length``)."""
    return swept(profile, Helix(pitch=pitch, height=height, radius=radius), label=label)


def strut(p0, p1, radius, label="strut"):
    r"""A round **strut** (cylinder) of given ``radius`` between two arbitrary 3D points ``p0`` and
    ``p1`` -- the member primitive for trusses, space frames, spokes and lattices."""
    return swept(Circle(radius), Line(p0, p1), label=label)


def thicken(profile, thickness, label="sheet"):
    r"""**Thicken** a planar ``profile`` (face / sketch) into a solid of the given ``thickness``,
    symmetric about the profile plane -- the "thicken sheet" verb (plates, laminations, membranes)."""
    s = extrude(profile, amount=thickness / 2.0, both=True).solid()
    s.label = label
    return s


def draft_extrude(profile, height, taper_deg, label="draft"):
    r"""**Extrude with draft**: extrude a 2D ``profile`` by ``height`` while tapering its walls by
    ``taper_deg`` (a positive angle shrinks the top) -- moulded / cast features, frusta."""
    s = extrude(profile, amount=height, taper=taper_deg).solid()
    s.label = label
    return s


def shell(solid, thickness, openings=None, label="shell"):
    r"""**Shell** (hollow out) a ``solid`` to a wall of ``thickness``, removing material inward.  With
    ``openings`` (a face or list of faces of ``solid``) the shell is OPEN on those faces (a cup / box);
    without, it is a CLOSED hollow shell.  The "shell solid" verb."""
    if openings is None:
        sh = (solid - offset(solid, amount=-abs(thickness)))
    else:
        sh = offset(solid, amount=-abs(thickness), openings=openings)
    sh = sh.solid()
    sh.label = label
    return sh


def fillet_edges(solid, radius, edge_filter=None, label="filleted"):
    r"""**Blend (fillet)** the edges of ``solid`` with the given ``radius``.  ``edge_filter`` is an
    optional callable ``edges -> edges`` to pick a subset (e.g. ``lambda e: e.filter_by(Axis.Z)`` for the
    vertical edges only); without it, every edge is rounded."""
    edges = solid.edges() if edge_filter is None else edge_filter(solid.edges())
    f = fillet(edges, radius)
    try:
        f.label = label
    except Exception:
        pass
    return f


def chamfer_edges(solid, length, edge_filter=None, label="chamfered"):
    r"""**Chamfer** the edges of ``solid`` by the given ``length``.  ``edge_filter`` is an optional
    callable ``edges -> edges`` to pick a subset; without it, every edge is chamfered."""
    edges = solid.edges() if edge_filter is None else edge_filter(solid.edges())
    c = chamfer(edges, length)
    try:
        c.label = label
    except Exception:
        pass
    return c


def grid_array(part, nx, ny, dx, dy, label=None, label_fmt="{base}_{i:02d}_{j:02d}"):
    r"""An ``nx`` x ``ny`` **rectangular grid** of copies of ``part`` at pitches ``dx`` (x) and ``dy``
    (y), returned as a labelled :class:`~build123d.Compound` -- the 2D "translate with copies" verb
    (pin grids, hole patterns, fin / standoff arrays)."""
    if nx < 1 or ny < 1:
        raise ValueError("nx, ny must be >= 1")
    base = label or (part.label or "part")
    children = []
    for i in range(nx):
        for j in range(ny):
            c = Pos(i * dx, j * dy, 0) * part
            c.label = label_fmt.format(base=base, i=i, j=j)
            children.append(c)
    return Compound(children=children, label=base + "_grid")


def path_array(part, path, count, label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` copies of ``part`` placed at equal parameter steps along ``path`` (a wire / edge) --
    the "array along a curve" verb (beads on a wire, bolts along a slot, stations along a spline)."""
    if count < 1:
        raise ValueError("count must be >= 1")
    base = label or (part.label or "part")
    children = []
    for k in range(count):
        t = k / (count - 1) if count > 1 else 0.0
        c = Pos(*tuple(path @ t)) * part
        c.label = label_fmt.format(base=base, k=k)
        children.append(c)
    return Compound(children=children, label=base + "_patharray")


# ---- boolean / slice / sheet-metal verbs -----------------------------------------------------------
def fuse(*parts, label="fused"):
    r"""**Boolean union** (add): fuse ``parts`` into one body, relabelled.  Overlapping inputs merge into
    a single solid; disjoint inputs remain separate solids in the returned body."""
    if not parts:
        raise ValueError("fuse needs at least one part")
    r = parts[0]
    for p in parts[1:]:
        r = r + p
    try:
        r.label = label
    except Exception:
        pass
    return r


def cut(base, *tools, label="cut"):
    r"""**Boolean subtract**: ``base`` minus each of ``tools`` (drill holes / remove material),
    relabelled."""
    r = base
    for t in tools:
        r = r - t
    try:
        r.label = label
    except Exception:
        pass
    return r


def common(*parts, label="common"):
    r"""**Boolean intersect**: the common (overlapping) volume of ``parts``, relabelled."""
    if not parts:
        raise ValueError("common needs at least one part")
    r = parts[0]
    for p in parts[1:]:
        r = r & p
    try:
        r.label = label
    except Exception:
        pass
    return r


def slice_solid(solid, plane=Plane.XY, keep="top", label="slice"):
    r"""**Slice** a solid by a ``plane`` and keep one side -- ``keep`` is ``"top"`` (the +normal side),
    ``"bottom"`` (the -normal side) or ``"both"`` (returns both halves).  The plane-cut / split verb
    (half models on a symmetry plane, sectioning)."""
    keepmap = {"top": Keep.TOP, "bottom": Keep.BOTTOM, "both": Keep.BOTH}
    if keep not in keepmap:
        raise ValueError("keep must be 'top', 'bottom' or 'both'")
    r = split(solid, bisect_by=plane, keep=keepmap[keep])
    try:
        r.label = label
    except Exception:
        pass
    return r


def bend_sheet(leg1, leg2, width, thickness, angle_deg, radius, n_arc=24, label="bent"):
    r"""A **sheet-metal bend**: a strip of cross-section ``thickness`` (in the bend plane) x ``width``
    (out of plane) that runs a flat leg of length ``leg1``, a cylindrical bend of inner ``radius``
    through ``angle_deg``, then a flat leg of length ``leg2``.  Built by sweeping the cross-section along
    the neutral-axis path (straight -> arc -> straight), so it is one clean labelled solid -- the
    sheet-metal "bend" verb (brackets, flanges, folded chassis)."""
    a = math.radians(angle_deg)
    rn = radius + thickness / 2.0
    pts = [(-leg1, 0.0, 0.0), (-leg1 * 0.5, 0.0, 0.0), (0.0, 0.0, 0.0)]
    for k in range(1, n_arc + 1):
        t = a * k / n_arc
        pts.append((rn * math.sin(t), rn * (1.0 - math.cos(t)), 0.0))
    pa = pts[-1]
    pts.append((pa[0] + leg2 * math.cos(a), pa[1] + leg2 * math.sin(a), 0.0))
    return swept(Rectangle(thickness, width), Spline(*pts), label=label)
