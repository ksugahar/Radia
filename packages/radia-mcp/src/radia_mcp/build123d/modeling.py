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

from build123d import (Axis, BuildLine, BuildSketch, CenterArc, Compound, Cylinder, Line, Mode,
                       Plane, Pos, RectangleRounded, extrude, make_face)

__all__ = ["annular_segment", "tube", "racetrack_coil", "polar_array", "linear_array",
           "mirrored", "assembly"]


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
    copies sit at the fan ends.  The CST-modeller "rotate with copies" / Radia segmented-array verb.
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
